# Experimental geometric preconditioning of final docking refinement

The default optimizer is unchanged. Select `--refinement_optimizer geometry` for
explicit rigid single-ligand Vina docking. Compare directly with
`--refinement_optimizer armijo`, which differs only in initialization of the
inverse Hessian. Both use the unchanged scoring function, fixed global search,
strict Armijo acceptance and ordinary BFGS update. Native is a separate
compatibility control.

For Cartesian heavy-atom positions r(q), define the diagonal metric
G_kk = (1/N) sum_i |d r_i / d q_k|^2. The three translational components are
one; angular components are estimated by central coordinate probes in the
same conformation chart as optimization, with h=1e-4 rad. Use
D_kk=1/clip(G_kk,0.25,64). The first scalar secant rescale is
s^T y/(y^T D y), preserving relative geometric scales. This diagonal
approximation ignores correlations among internal coordinates, is frozen
per refinement and has overhead O(N * number of DOFs). It is not a
Hessian of the scoring function. Coordinate probes do not count as score
evaluations, so both wall-clock time and actual evaluations must be compared.

Motivation: a positive geometry-related preconditioner can mitigate disparate
Cartesian displacement scales of rotations and torsions. Mones et al. used
force-field Hessian-derived positive preconditioners for molecular geometry
optimization (https://doi.org/10.1038/s41598-018-32105-x); the metric here is
a distinct, simpler hypothesis rather than a reproduction of their method.
A scaled-quadratic regression tests optimizer mechanics. It does not
establish improvement on docking or Ki/Kd. Experimental opt-in only until
paired results across withheld targets and comparable assays support it.
