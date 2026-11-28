# Cross-receptor transfer of all published HUP starts

An experimental pose generator for the two known human AChE pockets projects an entire graph-only search pool into the other receptor frame. Its prediction stage uses the independently crystallized GNT and E20 sites to superpose 48 pocket Cα positions. It uses no geometry or bound protein from the HUP crystal. All 21 PDBQT atoms and their original torsion tree are retained at export precision, including polar hydrogens; the 18 heavy atoms determine pose distances. There is one full stereochemical graph automorphism. The [protocol](../TRANSFER.md) was frozen after the individual 24k and 48k search results and before any pooled transfer or selection.

Every input is optimized in the receiving protein's dry native Vina potential and rescored at the reported coordinates using a common unbound energy of zero. Endpoints are ranked by the score printed by the receiver's native CLI (three decimal places). Ties retain the recipient first, then original pool order. Endpoint selection keeps up to nine, requires at least 1 Å heavy-atom separation and enforces a 100 kcal/mol window. Each cross-receptor arm pays for an ordinary 24k search on each protein. Its main comparator uses two distinct 24k search seeds on the same protein with identical local scoring and selection; the second seed is the next registered one cyclically. The shared 24k recipient-only arm and the unchanged single native 48k pool are additional controls. All search and optimization costs are recharged to every cell consuming the inputs, even where the actual archived calculation was reusable. Gradient calls use the observer only when it produces the exact native output. Local optimization can change the final order as well as the coordinates.

## Measured outcome

On 4EY7 at seed 1103503, the other-receptor search supplied a near-native candidate prospectively. Its direct receptor-frame error was 0.414 Å in the 4EY6 source, 0.418 Å after transfer and 0.550 Å after 4EY7 local optimization. Every endpoint descended from the original 4EY7 pool lay at least 3.317 Å from this candidate. It survived the score window, 1 Å duplicate filter and nine-mode cap at rank seven. Its recipient score was −7.355 nominal kcal/mol. The rank-one 5.744 Å endpoint scored −9.313, making the preference for the wrong pose 1.958 nominal kcal/mol. Across all six directional seed cases, this was the only distinct sub-2 Å endpoint produced by an alternate-receptor donor.

| Receiver and seed | Recipient-only Top-1 | Cross-receptor Top-1 | Two same-receptor Top-1 | Native 48k Top-1 | Cross-receptor Best@5 / Oracle |
|---|---:|---:|---:|---:|---:|
| 4EY6 1103501 | 3.623 | 3.623 | 0.414 | 4.562 | 2.799 / 2.598 |
| 4EY6 1103503 | 0.414 | 0.414 | 0.414 | 0.424 | 0.414 / 0.414 |
| 4EY6 1103507 | 3.627 | 4.578 | 3.623 | 0.408 | 2.301 / 2.301 |
| 4EY7 1103501 | 4.199 | 4.199 | 5.744 | 5.730 | 2.976 / 2.976 |
| 4EY7 1103503 | 5.744 | 5.744 | 3.063 | 3.117 | 2.981 / 0.550 |
| 4EY7 1103507 | 3.063 | 3.063 | 3.063 | 3.063 | 3.063 / 3.063 |

Switching the second search from the same receptor to the alternate one improved the winning RMSD in one cell, worsened it in three and left it unchanged in two. The paired mean increase was 1.388 Å on 4EY6 and 0.379 Å on 4EY7 (0.883 Å over the six cells). The one improvement at 4EY7 seed 1103501 came from preserving the recipient's original winner while the other same-receptor seed introduced a worse winner. In the opposite direction at seed 1103507, transferring 4EY7 modes into 4EY6 improved the available and published best RMSD from 3.627 to 2.301 Å, but a different donor endpoint won at 4.578 Å. These donor-specific observations distinguish access to a productive local basin from the ability to recognize it.

All 24 arm/seed cells published nine modes. Every one of the 108 unique local jobs completed, and all 108 observer output PDBQTs matched the production coordinates byte-for-byte. The optimizer made 2,535 observed gradient evaluations. The largest coordinate displacement from box center along any axis was 16.588 Å at input and 17.606 Å at output within the 18 Å half-width. Including two complete native searches, all local optimizations, explicit final scores, selection and output, median CPU per cross/same cell was 14.132/13.056 s on 4EY6 and 14.414/14.971 s on 4EY7. Corresponding native 48k medians were 0.390 and 0.356 s. This implementation pays the CLI setup cost for every local solve; throughput includes that expense. The deduplicated experimental predictor run, counting observations and packing raw evidence, took 169.28 wall and 168.75 CPU seconds.

## Files and commands

[predict.py](predict.py) consumes graph-only HUP and the two noncognate receptors, the blind search runtime/identity manifest, and the native pose pools. It never consumes a blind evaluation table or HUP crystal coordinates. [export.py](export.py) materializes selected models and origin tags into multi-model PDBQTs. [evaluate_results.py](evaluate_results.py) is an independent stage with access to the reference. It measures each source in its own receptor frame, then measures each transferred start and returned endpoint against the 4EY5 reference *directly* aligned to the recipient frame. Small changes on the protein alignment cycle can consequently occur before local optimization. The source and initial coordinates are the same for same-receptor controls. A donor endpoint is called distinct if every successfully optimized original recipient mode lies at least 1 Å away. Pose-related distances are in Å; native empirical scores are nominal kcal/mol. Temperature, chemical potentials, standard-state corrections and bound/unbound populations do not enter this search model.

- [Every local calculation](results/jobs.csv) includes donor, recipient, seed, mode, input/output SHA, optimization and final score, local displacement, gradient calls, box extents and measured runtime.
- [Raw models and native logs](results/raw.tar.xz) include start, returned and counting-observer coordinates, literal CLI arguments, score/log output and time by operation. In a replay working directory, place the unchanged receiver files under the raw root using 4EY6.pdbqt and 4EY7.pdbqt. The two receptor hashes and source PDBQT hashes are included in the [manifest](results/receipt.json).
- [Selection](results/selections.csv) retains every candidate and reason for keeping or removing it; [stage geometry](results/stages.csv) adds the structural measurements only after selection. [Published pools](results/published), [pool hashes](results/published.csv), [all pool metrics](results/cells.csv), [distributions](results/groups.csv), [paired differences](results/paired.csv) and [paired quantiles](results/paired_summary.csv) report both receiving proteins and all seeds.
- [Per-cell operation counts / costs](results/costs.csv) include native search, endpoint optimization, scoring, ranking and the separate validation runs. [Charged cost with publication](results/charged_cost.csv) includes writing the actual PDBQT. Separate timings for common inputs, the shared experimental execution, raw archival, evaluator and selection replay are in the receipts. Historical source receptor preparation is reused across all arms.

From the repository root use a clean output location. PRO points to the ordinary Pro executable; BUILD is its object directory and BOOST its Boost prefix. Use the same Python packages as for the [parent blind experiment](../blind/README.md). The raw archive includes timings specific to the local machine. The reconstruction check compares selections, models and geometric metrics instead.

```sh
D=benchmarks/ache_receptor_state/transfer
W=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
bash benchmarks/affinity_conditioning/build_refine_probes.sh "$BUILD" "$BOOST" "$W/observed"
python "$D/predict.py" "$W/work" "$W/result" "$PRO" "$W/observed/vina"
python "$D/export.py" "$W/result"
python "$D/evaluate_results.py" "$W/result"
for name in selections stages cells groups paired paired_summary; do
    cmp "$D/results/$name.csv" "$W/result/$name.csv"
done
diff -qr "$D/results/published" "$W/result/published"
```
