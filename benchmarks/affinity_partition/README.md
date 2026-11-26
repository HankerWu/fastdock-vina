# Affinity ranks, the torsion divisor and a held-out Factor Xa series

Removing Vina's global divisor improves two small development panels. On the external Factor Xa series it loses 77 of 295 net pair credits, and on the same-size, same-charge subset it loses one. A local statistical-mechanical replacement also fails. Keep the native default. The tested CLI fix gives each scoring function its own torsion-weight parameter.

## Paired tests

Each within-assay pair has a measured Ki or Kd ratio of at least three. Predictions are compared on the same receptor, seed and native top pose. Take the median of the within-block score differences and give an exact tie half credit. The receptor-free size rule favors the ligand with more heavy atoms. The raw-energy arm cancels the global divisor and keeps the same intramolecular/unbound reference. Pair estimates share ligands. Receptors and seeds are matched repeats. Summing independent ligand counts or treating 295 pair signs as independent trials would give misleading uncertainty.

| Target and matched endpoint | Ligands | Unique pairs | Native | Raw energy | Size | Single-rotor integral |
|---|---:|---:|---:|---:|---:|---:|
| T4 lysozyme L99A, ITC Kd | 9 | 19 | 12 | 17 | 17 | 3 |
| BRD4 BD1, Kd | 4 | 5 | 1 | 3 | 5 | 2 |
| Human CA II, Ki | 5 | 7 | 4 | 4 | 2.5 | 2 |
| Human Factor Xa, Ki | 28 | 295 | 223 | 146 | 74 | 204 |

For equal heavy-atom count and docked formal charge, T4L gives 1/2 native and 2/2 raw. FXa gives 9/17 native, 8/17 raw, and 3/17 with the single-rotor integral. BRD4 and hCA II supply zero qualifying pairs for this comparison. FXa 2P16 alone gives 187/295 native and 131/295 raw. The second structure 1FJS gives 231.5/295 and 157/295. Across both structures, deleting one ligand at a time gives accuracy changes from −28.2 to −21.6 percentage points for removing the divisor. The effect reverses across datasets with different correlations between ligand size and potency. This controls one large confound. It does not identify the one microscopic cause of the observed ranks.

**External target selection.** We used one congeneric series from Zhao et al., one assay (CHEMBL661804, PubChem AID 51686), one protein, one endpoint and one unit. All 34 original rows are exact Ki in nM. Six have an unassigned stereocenter in the curated graph. We used all remaining 28, with Ki from 0.026 to 420 nM, without cutting on potency, rotors, molecular size or docking score. ChEMBL's description explicitly identifies human coagulation factor X and the paper identifies activated FXa. ChEMBL assigns confidence 8 (homologous target) and supplies no variant, pH, temperature or replicate errors. PubChem's deposited assay identifies human F10. The reported within-assay ratio requires no temperature conversion. No claim of experimental per-pair statistical significance is made. [Original paper](https://doi.org/10.1016/S0960-894X(00)00139-6) · [ChEMBL assay](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL661804.json) · [Deposited assay](https://pubchem.ncbi.nlm.nih.gov/bioassay/51686)

We probed an apixaban-bound human receptor (2P16), then repeated the entire test on the human receptor co-crystallized with another chemotype whose amidine makes the Asp189 salt bridge (1FJS). Each independent cocrystal defines only the center of its 26 Å common box. Preparation uses protein ATOM records of both chains, major alternatives, no waters or ions, and the same ligand conformer recipe throughout. Each ligand uses two protonated amidines and a carboxylate where present. These states are an explicit chemical assumption. Meeko preparation on the fixed states gives net +1 or +2. Receptor waters, polarization, alternate protonation and receptor solution populations remain unresolved. Only affinity is evaluated for FXa. The new data establish no FXa pose accuracy. [2P16](https://www.rcsb.org/structure/2P16) · [1FJS](https://www.rcsb.org/structure/1FJS)

## Mechanism and controls

Vina reports a raw score divided by 1 + 0.05846 Nrot. A common divisor per ligand preserves its own pose order, but its magnitude couples the across-ligand flexibility penalty to every interaction in the raw score. An alternative hypothesis predicts a better affinity rank when a rotor receives a penalty from the loss of its accessible angular states on binding. [Trott and Olson](https://doi.org/10.1002/jcc.21334)

For each declared PDBQT rotor j connecting two nontrivial heavy-atom groups, the probe rotates the smaller fragment around the bond. It samples both choices with equal weight when the heavy-atom counts tie. With all other coordinates fixed, it computes

    Wj = F(I + V) − F(I) − V(0)
    F(E) = −RT log mean(exp(−E(θ)/RT))

I denotes the native intraligand potential and V the receptor interaction. Both use the original atom mapping, the original native tables and the same angular grid. The top-pose candidate score is raw energy plus the sum of these conditional corrections. T4L uses its assay temperature, 302.15 K. The other panels use a declared 298.15 K reference. A separate additive control gives every rotor RT log(3). A secondary analysis can reselect a saved candidate with the new objective. Its hCA II rankings depend on grid resolution and no hCA II top pose reaches 2 Å. This endpoint is unsuitable for selecting a default. Geometry labels never select a held-out affinity winner.

The approximation omits correlated torsion changes, relative translation/rotation, shape relaxation, changes in water and receptor populations, and non-native intrinsic torsion barriers. A rigid receptor can overestimate the loss of angular width even at accurate numerical quadrature. The 732 deposited development modes' reconstructed inter and intra energies agree with the native terms to within 0.023 and 0.009 kcal/mol, respectively, including coordinate print precision. Moving the only receptor atom outside interaction range returns exactly zero binding correction on all nine n-butylbenzene poses. The base 72-point angular grid and 144/288-point controls yield identical top-pose unique-pair totals for all four targets. The 36-point FXa control gives 208/295 and 6/17 equal-size/equal-charge, versus 204/295 and 3/17 at each finer resolution. The largest absolute 144→288 change in correction is 0.003 (T4L), 0.017 (BRD4), 0.036 (hCA II) and 0.024 kcal/mol (FXa). See `results/rotor_summary.csv`, `results/size_controls.csv` and the underlying cell and edge files.

## Independent physical evidence

The nine neutral graphs and Morton ITC ΔG values match the T4L panel. The experiments used the C54T/C97A/L99A construct. Absolute free energies from Mobley et al. and benzene-relative results from Khavrutskii and Wallqvist can be compared after subtracting scores within a chemical pair. No study requires fitting to the present native scores. The published calculations are retrospective and use structural information, distinct starting receptor states and simulation protocols. The 300 K simulations differ slightly from the 302.15 K assay. The 2007 apo-state calculation was initialized from the stripped benzene complex, PDB 181L. The present native pair gaps use two experimental receptor states and two seeds. [ITC experiment](https://doi.org/10.1021/bi00027a006) · [Crystal structures](https://doi.org/10.1021/bi00027a007) · [2007 absolute calculations, Tables 1–3](https://doi.org/10.1016/j.jmb.2007.06.002) · [2011 relative calculations, Table 1](https://doi.org/10.1021/ct2003786)

| Published arm on 19 threefold edges | Direction credit | Gap MAE, kcal/mol |
|---|---:|---:|
| Native docking | 12 | 1.178 |
| 2007 single docked mode, AM1-CM2 | 13 | 2.228 |
| 2007 multiple modes | 15 | 1.640 |
| 2007 multiple modes and Val111 reorientation | 16 | 1.762 |
| 2007 multiple modes, Val111 and AM1-BCC | 18 | 2.012 |
| 2011 single-reference thermodynamic integration, AM1-BCC | 19 | 0.717 |
| 2011 Hamiltonian replica-exchange calculation, AM1-BCC | 19 | 0.749 |

The 2007 apo/holo choice shifted o-xylene by 5.37 kcal/mol before Val111 reorientation and by 0.64 afterward. The relative calculations establish correct directions for these 19 edges even under ordinary integration. Replica exchange alone cannot explain that result. Both publications use correlated simulation legs and the nine measurements recur in all comparisons. The source error columns and pairwise SD bounds [|σa−σb|, σa+σb] preserve unknown covariance. The bound is a sensitivity calculation, not a confidence interval. The 2011 authors also report a different charge model on three problem ligands. Those post-hoc mixed results are retained in the CSV and excluded from the headline comparison.

## Replay

Python dependencies are RDKit, Meeko and the dependencies used by the repository's native pool benchmark. The saved score tables suffice for the four-target charge/size comparison. Run from repository root.

```sh
probe_dir=benchmarks/affinity_partition
python "$probe_dir/reconcile.py"
python "$probe_dir/size_controls.py"
python "$probe_dir/check_rot_weights.py" build/vina
```

For the C++ probe, build the native binary with the repository's build tool. Link the probe to the same non-main objects and Boost version. A common Linux build with system Boost uses the following commands.

```sh
python tools/build.py --output-dir build/affinity
c++ -O2 -DNDEBUG -std=c++14 -Isrc/lib "$probe_dir/rotor_partition.cpp" \
  $(find build/affinity/obj -name '*.o' ! -name 'main.o') \
  -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options \
  -pthread -o build/affinity/rotor_partition
```

The optional `fetch` action re-reads the ChEMBL assay and recreates `fxa_series.csv`. The pinned table gives both entry and preparation SMILES, all 34 inclusion decisions, and database activity identifiers. Meeko 0.8.0 and RDKit 2026.09.2 produced the archived runs. The PDB digest is checked during preparation. Independent native Late runs used two fixed seeds, one CPU each, exhaustiveness 4 and max_evals 50000.

```sh
for rec in 2P16 1FJS; do
  python "$probe_dir/fxa_holdout.py" prepare "/tmp/fxa-$rec" --receptor "$rec"
  python "$probe_dir/fxa_holdout.py" dock "/tmp/fxa-$rec" --receptor "$rec" --binary build/affinity/vina
  python "$probe_dir/fxa_holdout.py" collect "/tmp/fxa-$rec" --receptor "$rec"
done
python "$probe_dir/fxa_rotors.py" build/affinity/rotor_partition /tmp/fxa-2P16 /tmp/fxa-1FJS
python "$probe_dir/fxa_rotors.py" build/affinity/rotor_partition /tmp/fxa-2P16 /tmp/fxa-1FJS --angles 288
python "$probe_dir/size_controls.py"
```

To reproduce all development scans from the original checked native pools, use the same evidence layout as `benchmarks/strain_transfer/evaluate.py` and run `run_rotors.py PROBE EVIDENCE OUT` or add `--angles 288`. The pool checksums, fixed ligand preparation, model identity and unbound reference are validated before accumulating results. The optional `--evidence EVIDENCE` argument to `size_controls.py` rebuilds the development summary directly from deposited energy records.
