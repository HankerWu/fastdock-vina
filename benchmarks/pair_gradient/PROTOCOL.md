# Frozen comparisons

Starting point ac0a2f6. The frozen hypothesis is that calculating a derivative
from the energy used by the line search improves crystal-pose recovery at the
same budget. The linear arm differentiates the current interpolated energy.
The monotone cubic arm also removes internal slope discontinuities. Both use
the original nodal samples, map spacing, final fast scorer, and table footprint.
There are no fitted coefficients. Use the full case lists for all arms. A
better gradient check alone does not justify a change to a docking default.

## Local diagnostic

Use all eight complete experimental cases in
`geometry_comparison/basins/starts_and_relaxations.tar.xz`, three starts
per case (one crystal-informed torsion projection and two historical search
outputs). Use each *initial* PDBQT identically in all arms, one strict Armijo
local optimization, a hard cap of 12 or 36 joint objective/gradient calls,
one CPU, 24 Å box, prepared receptor. The two arms and parent are compared at
both caps. The crystal-informed starts diagnose preservation of an accessible
basin. Count hits at 2 Å for each start class and report per-start before/after
RMSD over archived chemically allowed mappings in the receptor frame. Score
start and end coordinates via the same parent single-point scorer on the same
original ligand graph with unbound energy fixed to zero. These post-hoc calls
are evaluation of the predictions, outside the predictive budget. Record each
local objective/call count, latency, and raw output. The eight systems have
already been studied in this repository, so this panel is a diagnosis.

## Prediction on independent structures

The new PDB entries selected before any docking results are ADORA2A 3EML and
3UZA, REN 1RNE and 2V0Z. Attempt each redock and each within-target reciprocal
cross-dock. Log preparation failures. Seeds 75467 and 98251 with search caps
3000 and 12000, one trajectory/CPU, at most nine modes, 24 Å box. Compare
upstream, parent, exact-linear and monotone-cubic, including all failures.
See the pose subdirectory for the precise preparation and audit. Receptor
alignment and crystal symmetry labels may be used to evaluate or localize
the known site. No template-derived ranking term is supplied to the search.

## End-to-end affinity transfer

To test whether changing only the optimizer's radial table alters *predicted*
same-assay rankings through changes of pose, run all 13 compounds from the
published Gal-1 and Gal-3 Table 1 archive. The 4 receptor states, both seeds
from the structural panel, 12,000 nominal evaluations, 1 trajectory, 1 CPU,
and 9 output modes apply identically to the four arms. Use the native
selection policy throughout. No contact templates, rescoring, or additional
conformer generation. All 13 compounds appear in both assays. Aggregate
scores over seeds in each receptor and average the two receptor scores as a
descriptive consensus. Report per-receptor sensitivities, within-assay rank
correlations, all at-least-threefold pairs, the same-heavy-atom subset, and
the predefined 12-member substituted triazole series. Record joint evaluations
and wall time. Matched nominal caps can incur different actual work. This is
a new test of the gradient variants on an existing assay panel.

## Mechanism follow-up after the three-arm model-gradient audit

The largest legacy explicit-coordinate error in the numerical panel came from
three receptor contacts in the last two cells of the 8 Å pair table. Both
consistent interpolants still forced the last active sample toward the
truncated zero node. A separate fourth Pro build now shifts only the explicit
receptor optimization potential by its pair-specific left-limit energy at the
cutoff. Its fast scorer and grid searches are identical to the exact-linear
build. It is compared on the same local, structural and assay panels.
Those data are development/replay evidence for this follow-up. The separately
frozen androgen receptor panel in holdout/ tests transfer to a new family.
