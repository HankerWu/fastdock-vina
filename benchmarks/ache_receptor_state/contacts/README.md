# Which contacts changed between the two dry AChE receptors?

The ordinary Pro explicit-receptor score was recomputed for 12 fixed HUP starts, six identities per admitted receiver. The 18,759 receptor–ligand atom pairs used actual native XS types and the two native lookup-table samples at each squared distance. All 12 full scores agree with their three-decimal CLI values to within 0.000373 nominal kcal/mol. The six registered contrasts for 4EY5 could not be evaluated because deposited DMS A605 fails the parent preparation criterion. This leaves 12/18 contrasts including the two identical-pose checks. [Registration](PROTOCOL.md) · [Every contrast](results/contrasts.csv) · [Reference scores and hashes](results/scores.csv)

| Fixed start minus original rank 1 | 4EY6 | 4EY7 | Change on 4EY7 |
|---|---:|---:|---:|
| Pro 24k Oracle (same pose) | 0.000 | 0.000 | 0.000 |
| Pro 48k Oracle | +1.445 | +0.651 | −0.794 |
| Official 24k Oracle | +3.387 | +4.164 | +0.777 |
| Projection minus Pro 24k winner | +1.323 | +2.726 | +1.402 |
| Projection minus Pro 48k winner | +2.768 | +3.377 | +0.608 |
| Projection minus official winner | +2.803 | +3.419 | +0.616 |

Positive means the test start scores worse under the same dry Pro objective. “Official” identifies the coordinate pool supplied to this diagnostic. The original search rank is held fixed on both sides. The projected start uses the direct query-to-receiver fit in each state, whereas original pool positions use the 4EY6-to-receiver fit. The projection contrast therefore includes the measured small transform-cycle difference described in the parent experiment.

For the official Oracle versus its rank 1, Tyr124 changed the contrast by +3.195 and Tyr337 by −2.608 nominal kcal/mol. Their repulsion components alone contributed +2.893 and −2.754. For Pro 48k the same two residues contributed +0.527 and −1.050, and Ser203 contributed −0.521 (−1.019 H-bond, +0.505 repulsion). For the crystal projection versus the original Pro 24k winner, Tyr337 contributed +2.409, while Gly120 contributed −0.669. The sign can depend on the pose pair even when the two protein structures are identical across comparisons. A matched residue's difference is an arithmetic partition of a fixed-coordinate empirical potential. These data do not identify what happens upon a mutation or measure any one residue's affinity contribution. [Signed residue components](results/state_changes.csv.gz) · [Within-state residue sums](results/residues.csv.gz) · [Individual pairs](results/pairs.csv.gz)

Residue numbers above are those of deposited chain A. The protein preparation renumbered residues, e.g. prepared Tyr121 = deposited Tyr124 and prepared Tyr328 = deposited Tyr337. Each state uses a one-to-one Cα position plus residue-name match. After the same pocket alignment used to transfer the ligand, Tyr124 OH and Ser203 OG move 1.390 and 2.035 Å; their Cα move only 0.224 and 0.198 Å. Distances in the motion table compare like PDB atom names, including convention-sensitive symmetric side-chain labels. [Number crosswalk](results/numbering.csv.gz) · [Matched heavy atom positions, displacements in Å](results/receptor_motion.csv.gz)

The Python calculation uses the native squared-distance bin scale of 32 Å⁻², radial cutoff 8 Å, two Gaussians, repulsion, hydrophobic and radial H-bond terms. After summing ligand-atom contacts, the authentic positive-energy cap is 1000 nominal kcal/mol per ligand atom. The nonlinear cap remainder is reported separately. The HUP model has zero active heavy torsions in the current Pro scoring divisor and zero scored mobile internal pairs in these starts. Within the diagnostic box the wall is zero. Printed CLI score differences have up to 0.001 nominal kcal/mol rounding uncertainty. Coordinates represent exported 0.001 Å snapshots and each static calculation is conditional on one preparation. The audit makes 12 independent score calculations, zero gradient evaluations and takes approximately 1 CPU second including parsing, typing, crosswalk and file output. [Receipt](results/receipt.json)

Run from the repository root after extracting the parent diagnostic raw archive. Pass the parent inputs.json as the fourth positional argument. BUILD and BOOST refer to the same build objects and Boost installation as used for Pro. Rebuilding and remeasuring the score tables is deterministic; incidental counters and hashes of local build paths in the receipt may differ.

```bash
D=benchmarks/ache_receptor_state/contacts
bash "$D/build.sh" "$BUILD" "$BOOST" "$WORK/typing"
OPENBLAS_NUM_THREADS=1 python "$D/decompose.py" "$RAW" "$WORK/typing" "$WORK/contacts" "$INPUTS"
for file in scores.csv contrasts.csv failures.csv pairs.csv.gz residues.csv.gz state_changes.csv.gz numbering.csv.gz receptor_motion.csv.gz; do
    cmp "$WORK/contacts/$file" "$D/results/$file"
done
```
