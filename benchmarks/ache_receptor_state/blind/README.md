# HUP graph-only cross-docking in two human AChE structures

The original 4EY6 result depended on the seed. With seed 1103503, the official 24k, Pro 24k and Pro 48k procedures all placed HUP within 0.43 Å at rank one in 4EY6. With seed 1103507 only Pro 48k did so. None of the 12 pools on 4EY7 published HUP within 2 Å. One target and ligand supply all 21 pools. The repeated seeds measure numerical search variability on those fixed inputs.

| Receiver, center | Search | Top-1 RMSD for seeds 1103501 / 1103503 / 1103507 (Å) | Mean | Median | Best pooled Oracle |
|---|---|---|---:|---:|---:|
| 4EY6, GNT | Official 24k | 4.563 / 0.425 / 3.623 | 2.870 | 3.623 | 0.425 |
| 4EY6, GNT | Pro 24k | 3.623 / 0.414 / 3.627 | 2.555 | 3.623 | 0.414 |
| 4EY6, GNT | Pro 48k | 4.562 / 0.424 / 0.408 | 1.798 | 0.424 | 0.408 |
| 4EY7, E20 | Official 24k | 4.211 / 3.119 / 4.304 | 3.878 | 4.211 | 2.994 |
| 4EY7, E20 | Pro 24k | 4.199 / 5.744 / 3.063 | 4.335 | 4.199 | 2.976 |
| 4EY7, E20 | Pro 48k | 5.730 / 3.117 / 3.063 | 3.970 | 3.117 | 2.976 |
| 4EY7, mapped GNT | Pro 24k | 3.163 / 3.155 / 5.753 | 4.024 | 3.163 | 3.035 |

On 4EY7, moving the center 4.812 Å from E20 to the mapped GNT center changed the three Pro 24k errors by −1.036, −2.588 and +2.689 Å. Moving just the center can change the error substantially in either direction. The two alternate centers both left the nearest published candidate at least 2.976 Å from the reference. In 4EY7 the ordinary and doubled Pro budgets shared exactly the same pool Oracle RMSD at *each* of the three seeds, while the Top-1 changed by +1.531, −2.627 and 0.000 Å. On 4EY6, seed 1103501 repeated the former rank-one coordinates exactly at rank two when the Pro budget doubled. The current observations establish losses at or before publication and from the ordering of published modes. The trajectory counters cannot identify the earliest visit or eviction of a particular structure.

All 21 searches succeeded and published 186 structures. Three of the 21 first poses exceed 5 Å and all three came from 4EY7. The raw mode count, receptor-frame symmetry-corrected error, all pools and all paired changes are available in [modes.csv](results/modes.csv), [cells.csv](results/cells.csv), [groups.csv](results/groups.csv), [paired.csv](results/paired.csv) and [paired_summary.csv](results/paired_summary.csv). The quantiles use linear interpolation within the three seeds of one receptor, arm and box. [Coordinate matches](results/coordinate_repeats.csv) compare the same heavy atoms without superposition. [Production logs, command lists and coordinates](results/raw) preserve the exact multiconformer PDBQT files.

## Inputs, resource use and verification

The search process regenerates a graph-only neutral HUP ligand from its archived CCD and requires exact agreement with the earlier PDBQT SHA-256. No HUP crystallographic position is loaded by this process. The receptor-specific ligand centroids use GNT or E20 and the alternate GNT center uses 48 common pocket Cα atoms of 4EY6 and 4EY7 (0.276 Å protein fit). The separate evaluator matches HUP to its 4EY5 structure via the full chemical graph and performs a direct source-to-receiver protein fit for each receiver. The [input manifest](results/inputs.json) and [evaluation receipt](results/evaluation_receipt.json) include the source digests, rigid transforms, atom equivalence count, and timing. The known-pocket library contains the experimental sites of other ligands.

The three 4EY6 seed-1103501 production PDBQTs match the predecessor exactly, including the two Pro pools produced with this environment's different Pro binary hash. The production runs consumed a combined 14.41 wall seconds and 14.38 CPU seconds across separate serial calls (14.93/14.90 for the complete search script including loading, ligand preparation and result writing). Fresh production replay repeated all 21 PDBQTs byte-for-byte. The official release does not report gradient counts. Each of its six matching upstream counting runs also produced an identical PDBQT. Sum of hunt, second and explicit gradient calls across the 15 production Pro cells was 1,515,687 evaluations. The corresponding count across the six official cells was 434,873. The six upstream counter runs added 8.03 wall seconds / 8.01 CPU seconds. [Costs and counters](results/runs.csv) · [Replay and counter data](results/replay/replays.csv). This accounting reuses previously prepared receptors; the preparation records reside in the parent archive. No new affinity measurement or binding free energy was calculated.

## Reproduce

From the repository root, define OFFICIAL as the Vina 1.2.7 Linux x86-64 release, PRO as the native build from the registered branch, and COUNTER as the upstream v1.2.7 source build with the repository's benchmarks/pair_gradient/upstream_counts.patch. The archive stores all build identities in the production and replay receipts. Each script expects an independent output location as shown. Python dependencies include NumPy, RDKit, Gemmi, Meeko and pdbfixer for the original ligand preparation module.

```sh
D=benchmarks/ache_receptor_state/blind
W=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/dock.py" "$W/results" "$OFFICIAL" "$PRO"
python "$D/evaluate.py" "$W/results"
python "$D/replay.py" "$W/results" "$W/replay" "$OFFICIAL" "$PRO" "$COUNTER"
for t in cells groups modes paired paired_summary coordinate_repeats; do
    cmp "$W/results/$t.csv" "$D/results/$t.csv"
done
```
