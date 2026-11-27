# Where the 3EML cross-dock succeeds or fails

This retrospective diagnostic replays 3EML ZMA in the 3UZA ADORA2A receptor using the two original seeds and the 12,000 nominal per-trajectory budget. One parent run failed where the linear-gradient candidate succeeded. The second seed provides a failure of both arms on the same direction. These four cells cannot estimate performance on a population of targets.

A diagnostic-only C++ patch copies the model for coordinate reconstruction and records the random start, every proposed conformation immediately before hunt minimization, every hunt endpoint and MC decision, every second-minimization endpoint and its immediate admission into the native trajectory pool, the trajectory exit, post-merge pool, input and endpoint of explicit refinement, and the independently serialized final output. It reads the engine's ligand atom serial map. No score or gradient call is added to these predictive trajectories. It checks pool admission using exact coordinate equality and energy. An admitted pose can still be evicted before trajectory exit. Internal BFGS line-search trial points and intermediate retry steps of explicit refinement were outside the trace.

RMSD uses the archived two chemical symmetry mappings for the 25 heavy atoms in the original receptor frame. The pocket transform matched 58 Cα atoms (0.5681 Å RMSD). The RMSDs of all final modes match the independently checked pose panel within 10⁻¹⁰ Å. Every traced PDBQT output and each search, acceptance, merge and explicit-refinement counter matched both the ordinary replay and the archived run.

| Arm and seed | Closest proposed | Closest accepted hunt | Closest rejected hunt | Best after merge | Best in output | Top-1 |
|---|---:|---:|---:|---:|---:|---:|
| Parent 75467 | 7.090 | 7.349 | 5.570 | 7.398 | 7.393 | 7.393 |
| Linear 75467 | 1.832 | 1.005 | 3.477 | 1.794 | 1.790 | 1.790 |
| Parent 98251 | 3.248 | 3.248 | 5.899 | 7.509 | 4.585 | 7.501 |
| Linear 98251 | 7.162 | 8.247 | 5.581 | 9.133 | 8.898 | 9.131 |

All distances are Å. Linear/75467 first reached the 2 Å criterion at accepted hunt **151**, from a 5.734 Å proposal to a 1.005 Å local minimum. The original save trigger did not request a second minimization at steps 151, 152 or 153. Accepted minima of 1.264 and 1.760 Å were observed at 152 and 153. At 154 the second minimization produced a 1.794 Å pose from a 1.774 Å hunt minimum. It survived native deduplication, merge, explicit refinement and ranking. Its final serialized RMSD was 1.790 Å. For parent/75467 the closest observed state was a rejected 5.570 Å hunt endpoint at step 144. No ≤2 Å pose occurred among its recorded stages. Both arms failed in the second seed. A 3.248 Å accepted parent hunt at step 154 never yielded a saved near-native pose. These findings distinguish a missing near-native among observed states from misranking an existing hit.

## Same-conformation control

The first successful linear hunt was subsequently replayed in isolation from its full internal state (translation, quaternion and five torsions, 17 significant digits). The script terminates the probe immediately after that one ordinary grid-hunt optimization. Both arms use the identical box, map objective, hunt cap, coordinate transform and parent BFGS settings. An extra initial-energy reading was made once per arm in this offline probe and is counted separately. Coordinates before minimization matched the captured proposal exactly in both arms. Every atom coordinate at the isolated linear endpoint matched its original predictive endpoint at 17-digit precision.

| Arm | Starting RMSD (Å) | Final RMSD (Å) | Starting objective | Final objective | Joint local calls |
|---|---:|---:|---:|---:|---:|
| Parent | 5.734172 | 1.351360 | 87.028693 | −7.107535 | 32 |
| Exact linear | 5.734172 | 1.005362 | 87.028693 | −7.462690 | 31 |

The derivative change improves this selected local endpoint by 0.346 Å with one fewer joint call. Both optimizers cross 2 Å from this start. Thus availability of this proposal along the actual trajectory is material to explaining the binary success difference. The unobserved alternatives to earlier steps are not identified by this control. The score reported here is the capped hunt objective, with no experimental affinity interpretation.

The full-trajectory diagnostic writes coordinates and copies model state. The four plain/traced time pairs, arm order and CPU pin are in `results/verification.json`. Each is a single timing. Three pairs had a positive trace overhead, and one plain run was substantially slower. No speed advantage is inferred. The trace writer and post-hoc fits are benchmark work outside production. Saved coordinates and configurations, per-event RMSDs, counter logs, optimizer source digests and each isolated endpoint are under `results/`.

## Reproduction

Start with the independent exports built by `../build_arms.py` and the Python dependencies of `../pose/paired.py`. Use the same compiler and Boost installation as those exports. The build helper leaves the uninstrumented executables in place. These scripts target one trajectory, one worker and a rigid single-ligand receptor.

```bash
B=benchmarks/pair_gradient
S="$B/stages"
W=$(mktemp -d)
tar xf "$B/pose/inputs.tar.xz" -C "$W"
tar xf "$B/pose/runs.tar.xz" -C "$W"
for arm in parent linear; do
  python "$S/build_trace.py" "$BUILDS" "$W/traced" --arm "$arm"
done
python "$S/replay.py" "$W" "$W/runs" "$BUILDS" "$W/traced" "$W/results"
python "$S/fixed_hunt.py" "$W" "$W/traced" "$W/results"
```
