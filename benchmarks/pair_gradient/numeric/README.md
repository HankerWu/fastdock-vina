# Radial table and conformation-gradient audit

The BFGS evaluator uses a trilinear receptor grid during normal Monte Carlo hunting and explicit receptor pairs during final refinement. Both stages use the by-atom tables for movable pairs. The explicit evaluator caps a positive receptor energy **after summing over the neighbors of one ligand atom**. The final reported Vina score and the map builder use the separate midpoint lookup. Numerical differentiation below always targets the objective returned in the same call as the gradient.

At scale `f` (32 Å⁻² by default), `q=f|x|²`, bin `i=floor(q)`, and within-bin coordinate `t=q-i`, the parent table returns `E=E_i+t(E_{i+1}-E_i)`. It independently interpolates precomputed radial slope divided by radius and multiplies this result by the interatomic displacement. The exact gradient inside a bin is `2f(E_{i+1}-E_i)x`. The exact-linear arm returns that gradient. At a shared bin edge it takes the right derivative. The nonzero slope generally jumps between bins. The cubic arm uses monotone Hermite interpolation in q and differentiates its own cubic. Its slope can vary rapidly within a narrow terminal bin. Evaluating three offsets in every bin at scales 16, 32 and 64 verified monotonicity, nodal interpolation and C¹ continuity up to rounding for all three audited chemistries.

## Experimental coordinates

Both receptor stages were evaluated at all 16 crystal-derived starts from the eight PDB entries in the torsion-proposals fixtures. They comprise cognate and protein-aligned reciprocal cross placements in four families. For each degree of freedom, the probe perturbs the original conformation through the production rigid-rotation and torsion parameterization. The reference differentiates the full energy returned by `model::eval_deriv` with paired central differences of ±10⁻⁵, ±10⁻⁶ and ±10⁻⁷ Å or radians. The CSV records the base energy, returned derivative and finite difference. The input coordinates are in the archived fixtures. The cap was 10 in the map hunt and 1000 for explicit refinement. These are correlated coordinate-level diagnostics on previously studied complexes.

Median absolute error against the same objective at a step of 10⁻⁷ is shown below. Units are kcal mol⁻¹ Å⁻¹ for translation and kcal mol⁻¹ rad⁻¹ for angles.

| Objective and variable | Components | Parent | Exact linear | Monotone cubic | Shifted explicit |
|---|---:|---:|---:|---:|---:|
| Hunt translation | 48 | 2.21e−8 | 2.21e−8 | 2.21e−8 | 2.21e−8 |
| Hunt global rotation | 48 | 3.68e−8 | 3.68e−8 | 3.71e−8 | 3.68e−8 |
| Hunt torsion | 92 | 5.48e−4 | 1.84e−8 | 1.82e−8 | 1.84e−8 |
| Explicit translation | 48 | 0.389 | 2.54e−8 | 2.10e−8 | 2.03e−8 |
| Explicit global rotation | 48 | 1.018 | 1.83e−7 | 2.13e−7 | 1.61e−7 |
| Explicit torsion | 92 | 0.133 | 1.60e−8 | 1.75e−8 | 1.29e−8 |

At the same step the largest absolute errors across global rotations were 9.80, 1.81e−4, 2.27e−6 and 1.81e−4, respectively. Increasing the central-difference step to 10⁻⁵ can straddle multiple radial bins in the sum over a receptor. The largest error for the exact-linear arm then reached 6.35 in global rotation. Nodal C¹ continuity in the cubic arm reduced this to 0.393. Neither gradient agreement nor smoothness establishes which pose or affinity ranking is better.

The pair-only test also uses derivatives of the *direct weighted empirical potential* as a second reference. For C_H–C_H at 2–6 Å and f=32, mean absolute radial error against this reference was 4.66e−5 for the parent, 1.05e−3 for exact linear, and 3.24e−5 for cubic (see unrounded values in the CSV). Thus exact linear can better differentiate its numerical objective while giving a poorer approximation to the continuous pair potential. For N_D–O_A the cubic mean direct error was lower than the parent's, although its maximum direct error was higher. The three within-bin offsets do not sample radial space uniformly. These errors describe the numerical probes rather than the distribution of atom pairs in proteins.

## The largest whole-model discrepancy

For the 2YWP cognate ligand, rigid y rotation leaves all internal pair distances invariant. Independent summation of all 1,661 explicit receptor pairs reproduces the reported y derivative. Three pairs in the last two radial bins contribute +9.8501 kcal mol⁻¹ rad⁻¹ to parent minus exact-linear torque. The remaining 1,658 contribute −0.05184. Parent torque is +3.72280 and exact-linear torque is −6.07547. The common energy gives −6.07546 by central difference. The cubic torque is −12.76768 for its modified energy. The shifted torque is +13.60965 with +13.60966 from finite differences of its own energy.

The largest contribution comes from ligand carbon 23 and Tyr A20 CZ (receptor serial 169), both XS C_P, separated by 7.999335847 Å. In the terminal bin the node energies of the original table are −0.003601434 and zero kcal/mol. The node at exactly 8 Å uses the truncated potential. The analytic potential approaches −0.003597219 kcal/mol from the left. Radial derivatives below are before atomwise capping, in kcal mol⁻¹ Å⁻¹.

| Direct Vina from inside the cutoff | Parent | Exact linear | Monotone cubic | Shifted explicit |
|---:|---:|---:|---:|---:|
| 0.002157996 | 0.922257 | 1.843781 | 2.481144 | 0.002157922 |

Subtracting the left limit separately for each receptor–ligand type pair allows the last bin to approach zero without spreading a hard potential jump over ~0.00195 Å. Inside the radius, the experimental explicit optimizer interpolates `U(r)-U(8−)`, with the last pre-truncation sample taken from the left. Outside, its contribution is zero. This gives a continuous energy at 8 Å and retains the small interior derivative at the left boundary. The outer derivative is zero. Pair-dependent offsets change the energy when the set of neighbors changes. The atomwise cap can also alter a summed gradient. This is a different refinement objective with numerical motivation. There is no direct measurement of a force or a binding free energy in this audit.

The shifted patch applies to explicit rigid-receptor Vina with its native hard 8 Å cutoff. All fast/midpoint samples, all grid evaluations and movable pair lookups use the exact-linear control. The test checks the midpoint values at every audited bin, the left reference for three XS pairs, the endpoint, and fallback behavior on a non-aligned 6.234 Å cutoff, smooth-cutoff Vina, and Vinardo. A shifted-potential cutoff is distinct from subtracting the boundary force. The molecular-dynamics literature studies both under different objectives. [Toxvaerd and Dyre (2011)](https://doi.org/10.1063/1.3558787)

The numeric baseline reproduced a near-cutoff residual previously seen in the repository's [cutoff-switch test](../../check_cutoff_switch.cpp). The new contribution is its localization on real generalized coordinates, its decomposition into physical atom pairs and comparison against an explicit-only correction. The earlier wider switch has separate [docking benchmarks](../../cutoff_switch).

## Replay

Use matching source trees and build object directories for the four arms described in the parent experiment. Set `--boost-prefix` for a custom installation. Use the compiler and ABI that built the docking objects. The script reads the published source and crystal-diagnostic fixtures in the repository and emits aggregated radial probes, every generalized derivative and the atom-pair decomposition. The verifier independently reconstructs the continuous C_P Gaussian tail and checks the control invariants, finite differences and summed torques.

```bash
D=benchmarks/pair_gradient/numeric
python "$D/replay.py" --native "$PARENT_SRC" "$PARENT_BUILD" \
  --linear "$LINEAR_SRC" "$LINEAR_BUILD" --hermite "$CUBIC_SRC" "$CUBIC_BUILD" \
  --shifted "$SHIFTED_SRC" "$SHIFTED_BUILD" --output "$OUT"
python "$D/verify.py" "$OUT"
```
