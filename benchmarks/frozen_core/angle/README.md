# Frozen-pool donor-angle correction to the Vina H-bond score

The branch's existing polar_geometry uses an angle to heavy neighbors of the acceptor, a Gaussian radial weight, and the best partner per ligand site. This candidate replaces the directional fraction of Vina's existing intermolecular N/O hydrogen-bond term, once for each contacting heavy-atom pair in the current pose. Interactions among atoms in the same molecule and interactions with metal ions retain the original energies.

For each N/O donor–acceptor heavy pair, obtain its attractive native value e from ScoringFunction::vina_hbond. The underlying empirical radial ramp saturates at 0.7 Å of heavy-atom surface overlap, with default weight −0.587439. For each bonded polar H, measure angle θ=D–H–A. Define g=max(0,−cos θ)². Take the maximum allowed g across hydrogens and both possible donor directions. A collinear forward H retains the native attraction; angles ≤90° lose it. The added raw penalty is −e(1−g), summed over all intermolecular pairs and multiplied by Vina's own conformation-independent torsional reporting scale. The constant is inherited from Vina. No coefficients were fitted, and no change to search or local optimization was made. This particular quadratic gate is a test hypothesis. Kortemme et al. report evidence for directional dependence in protein hydrogen bonds and derive a different multidimensional potential. [Primary physical study](https://doi.org/10.1016/S0022-2836(03)00021-4) · [Vina scoring function](https://doi.org/10.1002/jcc.21334)

## Results

Eligibility comprised all Native models in the six archived fixed-pool strata from geometry_comparison/upstream. The first two cap protocols reuse seeds 85231 and 85237. The per-complex protocol uses two fresh seeds, 93529 and 93563. Cognate includes ten complexes from eight proteins; reciprocal includes ten ligand/receiver directions from five proteins. Every comparison holds the pool fixed and reports success at heavy-atom symmetry-aware RMSD ≤2 Å.

| Protocol | Budget | Pools | Native top 1 | Corrected top 1 | Native / corrected top 5 | Oracle |
|---|---|---:|---:|---:|---:|---:|
| Cognate | Fixed | 20 | 1 | 1 | 2 / 2 | 3 |
| Cognate | Cohort CPU-matched | 20 | 3 | 3 | 9 / 9 | 10 |
| Cognate | Per-complex, fresh seeds | 20 | 2 | 2 | 8 / 8 | 8 |
| Reciprocal | Fixed | 20 | 7 | 7 | 9 / 9 | 9 |
| Reciprocal | Cohort CPU-matched | 20 | 10 | 7 | 14 / 14 | 14 |
| Reciprocal | Per-complex, fresh seeds | 20 | 11 | 9 | 13 / 13 | 14 |

The correction switched rank one in 32 of 120 pools. At fresh seeds the reciprocal losses were both FABP4, ligand 3FR2 in 2NNQ. At seed 93529 the native rank-one pose is 0.990 Å with score −9.003. Four native H-bond contacts incur a combined normalized penalty of +0.658. The original third pose is 5.202 Å and scores −8.800 with a penalty of +0.232, enough to overtake it. In BRDT 7MRH at the same seed, the correction switches a 3.079 Å pose to 1.233 Å. The native advantage of 0.218 is outweighed by penalties of +0.291 and +0.036. Fixed polar-H coordinates, including rotatable OH and NH orientations, were generated in molecular preparation; they are generally not experimentally resolved. The exercise does not identify which hydrogen orientations occur in solution. This treatment is experimental and the cross-docking results do not support default use.

1,073 poses contributed 2,095 native N/O H-bond contacts, including 125 dual-role pairs. Hydrogen ownership in the C++ scorer comes exclusively from the original prepared Vina bond graph; coordinates are transferred from the archived pose without re-inferring bonds. The independent Python audit uses the original ligand's SMILES-index and H-parent annotations and checks the unique nearest parent within each prepared receptor residue. Contact counts agreed for every pose, and the maximum energy discrepancy from the independent radial and angular formulas was below 5e−12 kcal/mol. There are no new global searches.

## Reproduce

Compile the research scorer against the same objects and Boost installation as the branch executable. From the repository root with a built Vina and default system Boost:

```bash
B=benchmarks/frozen_core/angle
W=$(mktemp -d)
OBJ=build/obj
g++ -std=c++17 -O3 -DNDEBUG -Isrc/lib "$B/donor_angle.cpp" \
    $(find "$OBJ" -name '*.o' ! -name main.o) \
    -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem \
    -lboost_program_options -o "$W/donor_angle"
python "$B/evaluate.py" "$W/donor_angle" "$W/results"
python "$B/audit.py" "$W/results"
```

When using a custom Boost prefix, pass the same include directory, link directory and runtime library path used for the branch build. NumPy and SciPy are needed for RMSD evaluation and independent verification. To score another native Vina pool, invoke donor_angle with a prepared rigid receptor PDBQT, the original prepared ligand PDBQT and a multi-model Vina output. Its CSV supplies the original rank, native and corrected scores, raw radial and angular terms and torsional reporting scale. The executable assumes the branch's standard Vina weights and printed native Vina energies.
