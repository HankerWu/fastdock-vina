# Experimental geometric preconditioning of final docking refinement

The default optimizer is unchanged. Select `--refinement_optimizer geometry` for
explicit rigid single-ligand Vina docking. Compare directly with
`--refinement_optimizer armijo`, which differs only in initialization of the
inverse Hessian. Both use the unchanged scoring function, fixed global search,
strict Armijo acceptance and ordinary BFGS update. Native is a separate
compatibility control.

With an explicit rigid receptor and one ligand, `--local_only` also accepts
either optimizer for a prepared pose. A positive `--max_evals` caps scoring
gradients across the entire local refinement, including box retries.
Read the actual count and the final optimization objective from
`Local refinement`. A third local-only arm, `rigid_metric`, measures translation
at the initial heavy-atom centroid and scales rotation by its radius of
gyration. It leaves torsion coordinates unchanged. Use the same receptor,
starting pose and gradient cap across arms.

For Cartesian heavy-atom positions r(q), define the diagonal metric
G_kk = (1/N) sum_i |d r_i / d q_k|^2. The three translational components are
one. Each angular entry is the sum of squared perpendicular distances
to the current axis over its rotating chemical heavy atoms, divided by the
heavy-atom count of the entire ligand. Subtree centroids and second central
moments are accumulated once in a postorder traversal; torsion coordinates
are recorded in the optimizer's preorder. Initializing this metric costs
O(N + number of torsions). Hydrogen and closure dummy atoms are omitted.

D_kk = 1/clip(G_kk, 0.25, 64). The first scalar secant rescale is
s^T y/(y^T D y). The metric is frozen for each refinement and does not
include off-diagonal correlations. Earlier revision 142c676 estimated
rotational entries using central coordinate probes with h = 1e-4 rad.
The numerical regression compares each analytic entry against that method
on transformed, branched and macrocycle ligand conformations. Floating-point
differences can change the local optimization trajectory.

Motivation: a positive geometry-related preconditioner can mitigate disparate
Cartesian displacement scales of rotations and torsions. Mones et al. used
force-field Hessian-derived positive preconditioners for molecular geometry
optimization (https://doi.org/10.1038/s41598-018-32105-x); the metric here is
a distinct, simpler hypothesis rather than a reproduction of their method.
A scaled-quadratic regression tests optimizer mechanics. It does not
establish improvement on docking or Ki/Kd. Experimental opt-in only until
paired results across withheld targets and comparable assays support it.
