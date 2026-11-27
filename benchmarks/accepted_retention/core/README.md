# Retaining one skipped hunt minimum

The paired discovery experiment changed one native refinement input in 14 of 48 docking runs. None of the incoming candidates was close to the experimental pose. Top-1, top-5 and returned-pool success were unchanged at both computational budgets. A displaced native pose in one cross-dock improved from 16.08 to 4.59 Å during explicit optimization. Selecting a low-energy, unrepresented hunt minimum did not reliably identify the best use of an explicit-refinement slot. The candidate remains a research patch against `31d073fdfb4baf7859b90a926f5debc71de26528`.

## Core intervention

Native MC minimizes a proposal under the hunt cap, tests acceptance, then performs a second minimization under the authentic cap if the current energy is a historical improvement or the trajectory's saved pool has room. An accepted state bypassing this branch can influence later mutations while never entering refinement itself. Prior passive traces located an experimental-like pose on this path. The existing trajectory archive can append additional candidates. Here the native trajectory, second minimizations, merged pool and number of explicit inputs are fixed. See [the original trace](../../pair_gradient/stages/README.md) and [the frozen discovery protocol](PROTOCOL.md).

Each trajectory records at most one skipped accepted state. The running record must have finite hunt energy lower than its predecessor, and be distinguishable from every currently saved pose by the native fixed-index heavy-atom RMSD in the receptor frame. The default distance is 2 Å. Persistent storage of a proposal is one internal conformation and one heavy-atom coordinate array. Checking novelty costs at most O(KN) for K saved representatives and N heavy atoms and is avoided when the energy cannot improve the record. Each native output stores one scalar, the hunt energy before its authentic-cap optimization.

The final selection takes place immediately after ordinary redundancy filtering. With at least two native inputs, an eligible reserve has a nearest distance of at least the larger of the merge radius and the final redundancy radius. Exact coordinate duplicates are also rejected at radius zero. The lowest-energy eligible reserve competes with the input having the least favorable originating hunt energy. Both energies come from the same hunt objective. Equal energies keep the native pose, and ties between trajectories retain the earlier trajectory. The ablation removes only the final strict energy comparison. There is at most one substitution across all trajectories. Subsequent explicit minimization recomputes the energy, and native postprocessing sorts all refined poses before assigning the common intramolecular reference. The policy uses neither a crystal pose nor a symmetry mapping. Evaluation uses independent graph-symmetry-corrected heavy-atom RMSD.

## End-to-end discovery

Paired cells use human A2A (3EML, 3UZA), renin (3G6Z, 3SFC), and androgen receptor (1E3G, 2AM9), four cognate/cross combinations per family, two seeds and two limits. The results use the archived inputs, the fixed 24 Å box, one CPU and one trajectory, nine output modes, and the original scorer and refinement parameters. Case order and order of arms per case are deterministically shuffled. Costs include input parsing, maps, kernel startup, the entire trajectory and postprocessing. Every row was run on the same pinned logical CPU. Single-run timings include scheduling variation. Execution used GCC 14.2.0 and Boost 1.83.0.

| Nominal MC evaluations | Cells | Replaced | Native → candidate successes at top-1 / top-5 / oracle | Actual joint evaluations | Native → candidate seconds |
|---|---:|---:|---|---:|---:|
| 3,000 | 24 | 4 | 1 / 1 / 1 → 1 / 1 / 1 | 78,536 → 78,558 | 4.220 → 4.131 |
| 12,000 | 24 | 10 | 1 / 1 / 1 → 1 / 1 / 1 | 294,450 → 294,462 | 4.960 → 5.160 |

The native and candidate arms had identical hunt and second-minimization calls and evaluation counts in every paired cell. The explicit-input counts were 95 and 110 across the low- and high-budget groups, respectively, for every arm. Wall retries and BFGS line search are data-dependent after a substitution. The tables include these differences. The unconditional control returned byte-identical pose files to the candidate in all 48 cells. Its aggregate timings were 4.209 and 4.757 seconds. These measurements establish neither a speed benefit nor an effect of the final energy gate on this panel.

A2A used 8 substitutions and had zero top-1, top-5 or oracle successes in 16 cells. Renin used 6 substitutions and had zero successes in 16 cells. AR used none and supplied the two repeated successes. The native pools contained a successful pose in two of 48 cells, and both ranked first in all arms. Top-1 RMSD was identical in 47/48 cells and increased by 0.344 Å in renin 3SFC → 3G6Z, seed 98251 at 12,000. Best returned RMSD improved in six cells and worsened in four. No change crossed the 2 Å threshold. The top-ranked severe-overlap count, as defined by the paired scorer, totaled 28 in each arm. These repeated runs share ligands and receptors. Treat each target and direction separately before interpreting pooled counts. [Per-target statistics](results/summary.csv) · [Paired results](results/cells.csv) · [Every returned mode](results/modes.csv)

## Where the incoming pose was lost

Passive replays for all 14 substitutions checked pose-file byte identity and native counters against the timed runs. Every preselection pool was identical at full precision. For every unchanged slot the coordinates remained identical before and after explicit refinement. Across the 14 changed slots, the closest incoming candidate had RMSD 5.085 Å before and 5.147 Å after refinement. Incoming states were closer to the crystal than the displaced native state in 11/14 comparisons before refinement and in 8/14 afterwards. In 3EML → 3UZA (seed 98251, 12,000), the native state at 16.079 Å reached 4.585 Å. The incoming state at 8.621 Å reached 8.669 Å. This local branch changed strongly, even though no arm achieved a successful pose. [Per-slot before/after RMSDs](results/changed_slots.csv) · [Full-precision stage digests](results/stages.csv)

The lower budget produced 108 skipped accepted states; the higher budget produced 2,378. Their geometry was examined 98 and 931 times and the one-state reserve changed 22 and 172 times, respectively. A reserve can overwrite an earlier geometry or become redundant by the end. These data show limited coverage by the selected reserve. They do not attribute every failure to the absence of a correct pose anywhere in the trajectory.

## Reproduction and checks

Run from the repository root in Linux with a C++14 compiler, Boost thread/serialization/filesystem/program_options and the dependencies of the existing paired pose benchmark (Python 3.12 with NumPy, SciPy, Gemmi, Meeko and RDKit). Use fresh output paths. A non-system Boost installation can be selected with `--boost-prefix PATH` on both build helpers. The build helpers export the original source and apply only the named research patch inside the requested directory.

```bash
B=benchmarks/accepted_retention/core
W=$(mktemp -d)
for arm in parent candidate ungated; do
    python "$B/build.py" "$arm" "$W/$arm" --jobs 4
done
python "$B/replay.py" "$W/discovery" --parent "$W/parent/build/vina" \
    --candidate "$W/candidate/build/vina" --ungated "$W/ungated/build/vina"
python "$B/controls.py" "$W/discovery" "$W/parent/build/vina" \
    "$W/candidate/build/vina" "$W/discovery/controls.json"
for arm in parent candidate; do
    python "$B/build_stage.py" "$W/$arm" "$W/$arm/build" "$W/trace-$arm"
done
python "$B/stages.py" "$W/discovery" --parent "$W/trace-parent/vina" \
    --candidate "$W/trace-candidate/vina"
tar -cJf "$W/discovery/runs.tar.xz" -C "$W/discovery" runs stage_runs
python "$B/audit.py" "$W/discovery"
```

The offline audit checks every returned file digest against the archived raw PDBQT, all mode-based endpoints, invariant MC counts and slot capacity, strict acceptance criterion, and one-slot isolation in the passive traces. Run `python "$B/audit.py" "$B/results"` to reproduce the distributed summary from captured files. Additional boundary checks cover `K=1`, custom final radius 3.5 Å, zero merge/deduplication radius, and three trajectories dispatched to one or two CPU threads. The sole-slot pose stayed byte-identical, and the parallel dispatches returned byte-identical pools in each arm. [Boundary results](results/controls.json) · [Timed and stage run logs with coordinates](results/runs.tar.xz) · [Executable identities](results/builds.json)
