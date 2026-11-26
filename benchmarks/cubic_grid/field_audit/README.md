# Equal-node receptor-field audit

This diagnostic uses 198 crystallographic ligand heavy atoms from all eight experimental complexes in the archived torsion-proposal input fixture. BACE1, BRD9, CHEK1 and EPHX2 each contribute two distinct complexes. Every atom was evaluated at its experimental position and at three reproducible random directions with displacements 0.25, 0.75 and 1.5 Å. Sampling is paired across interpolation methods. Experimental coordinates localize physically relevant regions for a numerical diagnostic only. These 792 positions are correlated within ligands and proteins.

The intact prepared ligand supplies native Vina XS types, including polar and nonpolar carbon. The deposited coordinate map supplies the receptor-frame positions. Each receptor retains its original 24 Å docking box, grid origin and 0.375 Å spacing. Every sampled position admits a full interior stencil. At each of the 64 neighboring nodes, the probe sums the same cutoff-limited Vina midpoint lookup table over real receptor atoms in their original order. Both interpolants share the resulting nodes. The native grid generator includes distances at the cutoff, while the explicit evaluator excludes the endpoint. The audit preserves these conventions. The search cap is 10 and is applied after interpolating each probe atom's receptor energy. This isolates the receptor term. It does not measure ligand strain, torsional entropy, whole-pose recovery or affinity.

The primary energy reference directly sums native fast pair-table evaluations at the query coordinates. An additional reference uses the explicit optimizer's smooth radial tables; its gradient uses the interpolated radial derivative returned by those tables. This is the field seen by explicit refinement under the audit's common cap. It is not the finite difference of the piecewise-constant fast score. In production, the final refinement uses a larger cap than the hunt stage. Python independently evaluates both interpolants as tensor products, using the endpoint Hermite values and centered finite-difference slopes for the cubic case. For every point the script requires both energies and all three derivative components to agree with the corresponding compiled search-grid functions within 2 × 10⁻¹¹.

| Probe set | Number | Energy MAE, linear / cubic | Mean gradient vector error, linear / cubic |
|---|---:|---:|---:|
| All | 792 | 0.02926 / 0.00648 | 0.27305 / 0.16063 |
| A receptor heavy atom within 0.01 Å of the 8 Å cutoff shell | 375 | 0.02892 / 0.00685 | 0.34959 / 0.25242 |
| Remaining probes | 417 | 0.02958 / 0.00615 | 0.20422 / 0.07808 |

Energies and gradient errors use kcal/mol and kcal/mol/Å, respectively. All four targets had lower average energy and gradient errors with cubic interpolation. The 95th percentile of gradient-vector error over all probes decreased from 0.88064 to 0.83137. Native hard truncation leaves rapid variation near 8 Å that neither method can reconstruct accurately from this grid spacing. At 19 points the cubic value extended beyond the minimum or maximum of the eight cell corners; its greatest excursion was 0.03849 kcal/mol before capping. An extremum inside a cell is not evidence on its own that the interpolant is less accurate. The row-level energy and gradient comparisons permit direct checks.

## Official grid-spacing control

Changing the existing Vina grid spacing gives a separate control for discretization. At 0.25 and 0.1875 Å we repeated every coordinate, displacement, atom type, direct reference score and direct gradient from the 0.375 Å audit. The search box is exactly 24 Å wide at all three spacings. On this specific box a map has 65³, 97³ and 129³ nodes respectively. Native trilinear interpolation always reads eight corners, whereas the cubic stencil accesses 64 nodes.

| Spacing (Å) | Relative map nodes | Energy MAE, linear / cubic | Mean gradient vector error, linear / cubic |
|---:|---:|---:|---:|
| 0.375 | 1.00 | 0.02926 / 0.00648 | 0.27305 / 0.16063 |
| 0.25 | 3.32 | 0.01321 / 0.00402 | 0.22202 / 0.14589 |
| 0.1875 | 7.82 | 0.00812 / 0.00322 | 0.19440 / 0.14396 |

At these points coarse cubic interpolation gave lower average errors on both metrics than native trilinear interpolation with half the grid spacing. This control only measures the receptor-field approximation. Actual construction and evaluation time require end-to-end measurements.

## Reproduction

Build from repository root using the usual build script. The following assumes the default build directory, a system Boost installation, and Bash. If Boost lives elsewhere, pass its usual include, library and rpath flags when compiling. The Python script requires NumPy.

```bash
python tools/build.py --jobs 4
audit=benchmarks/cubic_grid/field_audit
out=$(mktemp -d)
tar -xJf benchmarks/torsion_proposals/fixtures/inputs.tar.xz -C "$out"
c++ -O3 -DNDEBUG -std=c++14 -Isrc/lib "$audit/field_probe.cpp" \
  build/obj/{parse_pdbqt,model,grid,coords,quaternion,random,utils,conf_independent}.o \
  -lboost_filesystem -lboost_system -o "$out/field_probe"
python "$audit/check_field.py" "$out/field_probe" "$out" "$out/results"
python "$audit/check_field.py" "$out/field_probe" "$out" "$out/results/spacing_025" --spacing 0.25
python "$audit/check_field.py" "$out/field_probe" "$out" "$out/results/spacing_01875" --spacing 0.1875
```

The local measurement used GCC 14.2.0, Boost 1.83 and NumPy 2.3.5. The raw output records receptor-frame coordinates, atom serials and XS types, fast/smooth explicit energies, predicted energies, gradient components, proximity to the hard radial cutoff, and corner overshoot. All observations, including the largest errors, are retained in `results/field_probes.csv.gz`. Summary statistics by receptor, target, displacement and proximity to the cutoff are in `results/field_summary.json`.

The two original Vina papers describe the scoring/search distinction and the 1.2 features. Their docking results provide no numerical validation of this interpolant. GridFF studied gradient interpolation and cost for rigid substrates using a different cubic spline and different physical systems. [Trott and Olson](https://doi.org/10.1002/jcc.21334) · [Eberhardt et al.](https://doi.org/10.1021/acs.jcim.1c00203) · [GridFF](https://doi.org/10.1021/acs.jctc.5c01223)
