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
