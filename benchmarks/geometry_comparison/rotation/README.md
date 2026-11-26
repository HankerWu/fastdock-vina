# Body orientation as a numerical hypothesis

The three world-axis diagonal entries of the rotational displacement metric
generally change when the whole molecule is rigidly rotated. An experimental
control replaces those entries by their mean before inversion and clipping.
That trace is rotationally invariant. Torsional entries, search, energy,
candidate count, convergence tests and all tuning parameters are identical.
The relative anisotropy of the three axes is discarded. An improvement on
receptor-frame RMSD or same-assay order is an empirical question.

Frozen comparison: use every eligible primary human cognate case in
`../pose_cases.csv` (excluding the reciprocal validation cohort), both seeds
63419 and 63421, on the existing 24 Å, 3 × 12000 search protocol. Compare with
analytic coordinate-wise geometry on the same search pool. Report top-1 and
any-mode RMSD at 2 Å and per-cell results. Do not select cases or tune the
isotropic mean from outcomes. In a separate scratch checkout of the tested
source, apply `git apply benchmarks/geometry_comparison/rotation/isotropic.patch`,
build normally, and pass `--refinement_optimizer geometry` to this binary.

Both the diagonal and isotropic arms produced 0/16 top-1 and 0/16 emitted-pool
hits. This low-coverage test does not establish equivalence. The isotropic
control is available as a research patch. To replay, run
`python rotation/run.py PREPARED OUT ANALYTIC_VINA ISOTROPIC_VINA` from the
parent directory, using a fresh input directory populated by pose_prepare.py.
