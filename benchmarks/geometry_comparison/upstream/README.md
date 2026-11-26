# Official Vina controls and search under a CPU budget

Compare official AutoDock Vina v1.2.7 (8eb40404), FastDock Native, Late with native refinement, Late/Armijo, Late/geometry with finite-difference lever arms (branch 142c676), and Late/geometry with analytic lever arms in this revision. All three executables used the same C++14 compiler, optimization flags, machine and Boost libraries. The CLI controls used one CPU, three independent trajectories, 0.375 Å spacing, a 24 Å cube for a known pocket, nine output modes and a 100 kcal/mol output window. Experimental ligand coordinates supplied evaluation references and pocket centers. Ligand preparation used the molecular graph for starting coordinates.

[Official tag](https://github.com/forlilab/AutoDock-Vina/tree/v1.2.7) · [Vina 2010](https://doi.org/10.1002/jcc.21334) · [Vina 1.2 paper](https://doi.org/10.1021/acs.jcim.1c00203) · [Frozen human input manifest](../pose/manifest.json)

## Evaluation caps held fixed

The first cohort contains 10 human cognate structures, eight distinct targets and two new random seeds (85231, 85237). The archived reciprocal panel contains five additional human targets (EPAS1, FKBP1A, RBP4, FABP4 and DYRK1A), two directions per target and the same two seeds. Each cohort contributes 20 paired cells. A hit has whole-ligand symmetry-corrected heavy-atom RMSD ≤2 Å in the receiver frame, after aligning the local protein Cα atoms for reciprocal docking. Top 5 and oracle take the best of the first five or all emitted poses. Multiple random seeds share the same experimental structure.

| Protocol | Search cap per trajectory | Arms | Top 1 / top 5 / oracle | Median top 1, official / native / Late / Armijo / old geometry / analytic geometry (Å) |
|---|---:|---|---|---|
| Cognate | 12,000 | All six | 1 / 2 / 3 of 20 in every arm | 7.027 / 7.027 / 7.027 / 7.034 / 7.063 / 7.063 |
| Reciprocal | 12,000 | All six | 7 / 9 / 9 of 20 in every arm | 2.567 / 2.567 / 2.567 / 2.566 / 2.580 / 2.580 |

On each cell the five fork arms used identical global search evaluation counts. Summed across the 20 cells the hunt and secondary minimization steps used 722,797 evaluations on cognate and 721,800 on reciprocal structures. Per-trajectory caps can be exceeded by the last local search call. Independent explicit refinement is counted separately in the logs. Official Vina exposes the requested cap at its CLI. At fixed settings the two geometry implementations retained the same endpoint counts, with a largest change of 0.057 Å in top-1 RMSD. Only two of 40 complete PDBQT pools were byte-identical between the two geometry implementations. The optimizer decisions and lower-ranked poses are sensitive to numerical differences.

Official and fork Native produced identical reported top-1 RMSD and score in 38 of 40 cells. On both BCL6 7QK0 seeds, running the identical commands again with `--no_refine` yielded byte-identical *complete* grid-scored pools, localizing the observed divergence to the path with explicit refinement. With explicit refinement the maximum top-score shift was 0.073 kcal/mol and the maximum change in RMSD to the experimental structure was 0.039 Å. The 2 Å classifications agree. [BCL6 output hashes and measurements](bcl6_no_refine.json)

The official release already supports flexible receptor side chains, AD4 maps and explicit-receptor refinement. The official develop snapshot bb64e5e differs from the fork's upstream starting revision 3c65c0b in the C++ source only in precomputed pair-table deduplication and cutoff sizing. FastDock already implements typed receptor tables, sparse pair tables and longer macrocycle glue-pair tables. The official change concerns memory use for existing interaction terms. [Upstream comparison](https://github.com/forlilab/AutoDock-Vina/compare/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645...bb64e5e26b520569a86360d69248069b5e665b3f)

## CPU-calibrated search depth

The separate calibration seed 74317 measured end-to-end child-process CPU time (user + system) on every query/receptor at official 12,000 and FastDock Native at 12,000 and 96,000 evaluations per trajectory. Per cohort, fit a straight line through the two FastDock measurements, solve for the official CPU time and round to the nearest thousand. Pose scores and reference RMSD never enter the fit. The frozen cap applies equally to FastDock Native and analytic geometry in the confirmation. The fixed-cap results above motivated this test. These first confirmation seeds repeat the fixed-cap panel. CPU includes startup, maps, search, refinement and writing. Jobs used one CPU per process and ran concurrently across cases. This measures process cost on this machine. Per-cell cost varies with the protein, ligand and random search path.

| Protocol | Official / FastDock cap | Official, CPU s and top 1 / top 5 / oracle | FastDock Native | Late / analytic geometry |
|---|---:|---|---|---|
| Cognate, 20 cells | 12,000 / 122,000 | 31.09 s, 1 / 2 / 3 | 34.61 s, 3 / 9 / 10 | 35.21 s, 3 / 9 / 10 |
| Reciprocal, 20 cells | 12,000 / 98,000 | 30.82 s, 7 / 9 / 9 | 27.74 s, 10 / 14 / 14 | 27.77 s, 10 / 14 / 14 |

The cognate Native top-1 gains are the two seeds on BPTF 7F5D. Native gained seven oracle hits on PHIP, BPTF, HSP90AA1 and BRDT. On reciprocal docking the three improved first poses belong to EPAS1, FKBP1A and FABP4, with no lost first-pose hits for this seed cohort. Per-cell Native CPU / official CPU ratios span 0.56–1.64 on cognate and 0.52–1.55 on reciprocal cases. Individual CPU times and actual gradient counts are in the CSV files. Three trajectories is a low sampling budget. These data examine search allocation for the prepared structures.

### Fresh seeds and per-complex budget calibration

Reuse the timing-only pilot, fit the line *separately* for each query and receptor, then run two untouched seeds 93529 and 93563. All slopes were positive. For several small molecules the predicted cap extrapolated beyond 96,000. Vina can reach its heuristic MC step limit before the cap, so a linear cost model is approximate. The per-complex fit was a follow-up to the cohort comparison. This replicates the sampling test across new seeds while reusing the chemical systems.

| Protocol | Official CPU s and top 1 / top 5 / oracle, out of 20 | FastDock Native | Late / analytic geometry |
|---|---|---|---|
| Cognate | 34.77 s, 2 / 3 / 3 | 37.86 s, 2 / 8 / 8 | 37.80 s, 2 / 8 / 8 |
| Reciprocal | 34.51 s, 5 / 7 / 7 | 27.88 s, 11 / 13 / 14 | 27.80 s, 11 / 14 / 15 |

On cognate receptors the deeper search gained a first-pose hit on BPTF and lost one on BAZ2A. The five Native oracle gains came from PHIP, BPTF and BRDT. On reciprocal receptors the six top-1 improvements came from EPAS1, FKBP1A, FABP4 and DYRK1A (no losses). The per-cell Native/official CPU ratio ranged from 0.65–1.84 (median 1.04) on cognate and 0.38–1.15 (median 0.85) on reciprocal directions. Timings and attained search counts are provided per case. The larger search budget repeats the gain in candidate coverage across both panels. The extra coverage led to a reproducible first-pose advantage only on the reciprocal targets tested here. The two optimized-engine controls are essentially tied for top 1 at the larger cap, so the evidence primarily supports allocating CPU to search depth.

## Reproduction

The cognate inputs include the complete preparation manifest and pre-generated graph-based ligands. The reciprocal inputs come from the previously shipped five-target reciprocal study. The complete raw PDBQT files and logs for all 480 reported searches are in the six raw.tar.xz files. Calibration CPU measurements and individual fitted caps are recorded separately. The matched runs reuse the seed set for the fixed-cap controls and share their experimental structures. `verify.py` reads the emitted coordinates, scores, and hashes, recomputes receptor-frame RMSD using a scalar distance reduction and a Kabsch protein alignment, and checks the counts. The two verification paths share the released atom mappings. NumPy and SciPy are sufficient for these prepared-fixture replay scripts.

From the repository root, substitute paths to the official v1.2.7 executable, the original 142c676 executable, and the analytic executable in this revision. Build the upstream release with equivalent optimization/Boost settings if comparing timing. `matched.py` recalibrates for local hardware.

```bash
D=benchmarks/geometry_comparison/upstream
mkdir cognate-fixtures reciprocal-fixtures
tar -xJf benchmarks/geometry_comparison/pose/inputs.tar.xz -C cognate-fixtures
tar -xJf benchmarks/cubic_grid/human_panel/prepared.tar.xz -C reciprocal-fixtures
for protocol in cognate reciprocal; do
  fixtures="${protocol}-fixtures"
  python "$D/replay.py" "$protocol" "$fixtures" "runs-$protocol" "$OFFICIAL" "$BRANCH" "$ANALYTIC"
  python "$D/matched.py" "$protocol" "$fixtures" "matched-$protocol" "$OFFICIAL" "$BRANCH" "$ANALYTIC"
  python "$D/matched.py" "$protocol" "$fixtures" "individual-$protocol" "$OFFICIAL" "$BRANCH" "$ANALYTIC" \
      --per-complex-calibration "matched-$protocol/calibration.csv"
  for stage in "$protocol" "$protocol/matched" "$protocol/per_complex"; do
    python "$D/verify.py" "$protocol" "$fixtures" "$D/$stage"
  done
done
```

For the 7QK0 diagnostic run the official and FastDock Native commands from `replay.py` using the cognate fixture, seed 85231 or 85237, plus `--no_refine`. Compare the full PDBQT hash against `bcl6_no_refine.json`. The experimental measurements here address structure recovery. A Ki/Kd ranking test needs comparable measured constants for several molecules on the same target.
