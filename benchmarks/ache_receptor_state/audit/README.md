# Independent coordinate and native-record checks

The independent evaluation reads the PDB and CCD deposits, constructs the HUP chemistry from the deposited atom/bond tables, checks stereochemistry from every exported geometry and reads atom identities from the PDBQT SMILES map. It uses SciPy's proper rotation fit on common pocket Cα sites. It does not import the new predictor, its ranking code, its coordinate reference, or its score table to construct answers. Table comparison occurs after the reconstruction. For the fixed-position energies it reads the energy decomposition printed by the native scorer, checks the common unbound reference, then calculates the originally registered differences. The independent atom-pair energy implementation is in the sibling contacts directory.

From the repository root, with NumPy, SciPy, pandas, Gemmi and RDKit installed:

```sh
D=benchmarks/ache_receptor_state
W=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/audit/audit.py" "$W"
python "$D/audit/previous_refine.py"
python "$D/audit/transfer.py" "$D/transfer/results" "$W"
python "$D/audit/trace.py" "$W/trace_receipt.json"
```

The three historical result commits are ancestors of the starting commit `27de934`. The earlier water envelope, geometry and AChE checks were actually rerun; their digest and time receipts are in [previous_replays.json](results/previous_replays.json). To reproduce the historical inputs and intermediate tables, use the commands in [the water envelope](../../water_rank_envelope/README.md), [the AChE transfer](../../water_rank_envelope/transfer/README.md) and [the two geometry panels](../../pose_error_geometry/README.md). The FKBP script independently reads the prior 59 native jobs and their gradient-observer logs, compares raw coordinate and graph hashes, and reconstructs both the net and substitution signs on the 26 registered cells from printed scores.

[Coordinate and cost checks](results/both_receipt.json) · [Frozen-pair contrasts](results/diagnostic_gaps.csv) · [Projection against all 26 pool endpoints](results/diagnostic_projection.csv)

The [transfer check](results/transfer_receipt.json) verifies every donor atom in recipient coordinates, exact native/observer output bytes, all 324 candidate decisions, all 216 exported models and the per-arm energy/gradient charges. The [trajectory check](results/trace_receipt.json) decodes both binary streams with the independently reconstructed frame, follows membership lists and checks the origins of all 18 released modes.
