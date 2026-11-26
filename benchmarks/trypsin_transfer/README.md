# Transfer of substituted amidiniums between bulk water and bound poses

A para-amino and a para-methoxy benzamidinium both reach Ser195 and Asp189 in guided bovine trypsin docking. The measured association constants are 127,000 and 26,500 M⁻¹. At 25.1 °C the methoxy ligand therefore binds less favorably by **+0.929 kcal/mol**. The median paired Vina difference from the eight released double-shell searches is **−0.0185 kcal/mol**. We tested whether allowing electronic relaxation in a chemically explicit piece of the pocket, and subtracting the same ligand's hydrated energy, recovers the difference. A full-protein fixed-charge continuum calculation acts as a separate comparator. On the unambiguous 4.79-fold edge the fine-angle local quantum estimate is essentially tied (−0.002 kcal/mol), while the whole-protein fixed-charge OBC2 endpoint gives +0.725 kcal/mol on these poses. The latter loses pairs on the independent FXa series. Only the first real pose of each stored search enters the comparison. No electronic energy or experimental label is used to choose a pose.

## Fixed inputs and observables

Five tabled para substituents (H, NH₂, CH₃, OCH₃ and CONH₂), two independently prepared receptors (1BTY, 3PTB), four matched seeds (271828, 314159, 161803, 141421) and two search arms (native, double-shell) give 80 first poses. All five ligands have charge +1. The same-buffer, 25.1 °C trypsin ITC series uses 50 mM Tris, 10 mM CaCl₂, pH 8.0. The original tables say CONH₂ for the fifth substituent, while their materials text says amidino. Four unambiguous compounds form six pairs in total and one pair above the threefold threshold (NH₂/OCH₃). Including CONH₂ yields three pairs above that threshold. Results from different seeds and receivers are matched repeats of the same chemical pair. All analyses are retrospective. [ITC study](https://doi.org/10.1046/j.1432-1327.2001.01991.x) · [Original series and search protocols](../directional_ion_pair/README.md)

The primary model uses two acetate fragments (Asp189, Asp102), a methanol (Ser195), and 4-methylimidazole (His57, neutral with Nδ1-H). The methyl positions are the experimental Cβ atoms; the carboxylate and heterocycle use the receptor side-chain heavy-atom positions. The cluster includes the ionic anchor and catalytic triad. Crystal coordinates remain frozen. Prepared ligand heavy atoms and donor hydrogens are transferred by the Meeko graph/parent map; RDKit fills omitted hydrogens and minimizes only those additional H positions with MMFF. The Ser195 OH angle around its Cβ–O bond is the sampled internal degree of freedom. The other protein atoms, structural waters and Ca²⁺ are absent in the quantum fragment. The full-protein comparator uses the archived dry receptor with PDBQT charges. Neither model simulates water exchange or receptor reorganization.

GFN2-xTB solves the semiempirical electronic structure at each geometry, with multipole electrostatics and density-dependent dispersion. ALPB represents continuum solvation. For the primary calculation the free ligand, free pocket and complex are each exposed to the water parameterization. A bound-octanol variant changes the continuum of the complex and free pocket only; this changes empirical cavity and dispersion-related contributions as well as dielectric response. It is a sensitivity test, with no assignment of octanol as the actual local solvent. Two smaller pockets (Asp189 only and Asp189+Ser195), the alternative neutral His tautomer, the protonated His, and the original single Ser orientation test which assumptions control the sign. [GFN2-xTB](https://doi.org/10.1021/acs.jctc.8b01176) · [ALPB parametrization](https://doi.org/10.1021/acs.jctc.1c00471)

For each ligand pose L we compute the fixed-geometry hydration leg, the fragment association leg relative to the same vacuum geometry, and their difference. Energies use kcal/mol. The angular log sum includes the same six evenly spaced Ser-OH rotations in each environment.

```math
F(X)=-RT\log\left(\frac{1}{6}\sum_{j=0}^{5}e^{-E(X,60j^\circ)/RT}\right),\qquad
W=E_L^{aq}-E_L^{vac},\qquad B=F(PL)-F(P)-E_L^{vac},\qquad T=B-W.
```

Asp189 alone has one state. We difference T between ligands inside each receptor/seed block, then take the median of these paired gaps. This removes a common vacuum energy for each composition. The free ligand is held in the same geometry as its own bound pose. The result omits free-solution conformational selection, translational and rotational standard-state terms, and most protein motion. No direct comparison of total single-point energies between ligands of different composition enters T. Absolute local transfer scores should not be interpreted as experimental standard binding free energies.

## Measured ranking and model dependence

| Double-shell search, one vote per chemical pair | NH₂ → OCH₃, kcal/mol | Correct on six unambiguous pairs | Mean absolute error, six pairs, kcal/mol |
|---|---:|---:|---:|
| Experiment | +0.929 | | |
| Native Vina | −0.0185 | 3 | 0.485 |
| Asp189 only, in water | +0.250 | 3 | 0.687 |
| Asp189 and Ser195, in water | +1.095 | 5 | 0.377 |
| Asp189, Ser195, His57(δH), Asp102, in water | +0.284 | 4 | 0.523 |
| Same, His57(εH), in water | +0.703 | 3 | 0.461 |
| Same, His57(δH), in octanol with aqueous ligand | −0.122 | 2 | 0.569 |
| Full receptor, fixed-charge Coulomb plus OBC2 | +0.725 | 5 | 0.302 |

The primary six-angle local quantum calculation gives the experimental sign for five of eight NH₂→OCH₃ repetitions. Per-receptor medians are +0.534 and −0.008 kcal/mol (1BTY and 3PTB). The mean gap of its bound leg is +0.031 and that of its water leg −0.124, giving +0.154 after subtraction. Its median of paired differences is +0.284; subtracting two independently computed medians gives a different statistic. For H→NH₂ the same calculation predicts +0.268 kcal/mol versus experimental −0.615. Jiao et al. published an independently simulated explicit-water polarizable complex/water difference of −0.28 for H→NH₂; their series does not report OCH₃. Both edges are needed to judge sign and transferability. [Jiao et al., Tables 1 and 3](https://doi.org/10.1002/jcc.21268)

**Angular convergence changes the mechanistic conclusion.** For all 16 NH₂/OCH₃ shell poses, an independent 48-point uniform grid from 0° through 352.5° supplies four nested resolutions. Recomputing the median matched primary score at 60°, 30°, 15° and 7.5° gives **+0.2844, −0.0258, −0.0015 and −0.0022 kcal/mol**. The finest grid gets 4/8 repeated directions correct, with per-receptor medians −0.183 and +0.284. The main 60° grid spuriously suggested a small positive effect. At the two finest resolutions the local model fails to resolve the observed +0.929 difference. The convergence study includes 768 complex evaluations; each new angle uses the corresponding apo reference.

The independent 295 K neutron/deposition 5MO0 (reported crystal pH 7.5) places H/D at His57 Nδ1 with occupancies 0.09/0.91 and an additional D at Nε2 with occupancy 0.35. The original prepared models place a proton at Nε2 alone. The deposited Ser195 Oγ has two locations at occupancies 0.47 and 0.53, with a resolved D at 0.53 in one of them. A repeated coarse-grid calculation on the doubly protonated His57 microstate gives +0.605 kcal/mol for NH₂→OCH₃. The full panel now contains 2,480 quantum cluster energies, 62 isolated pocket values, and 160 ligand values (bound-pose geometry in vacuum and water). The neutron crystal supports testing different microstates. X-ray/solvent H placement, H/D exchange, alternate heavy-atom locations and solution pH each affect interpretation: the 0.35 crystal occupancy is not a measured unbound-state weight for this assay. A state-specific energy shift alone cannot establish that His proton uptake causes the affinity difference. [Neutron paper](https://doi.org/10.1038/s41467-018-05769-2) · [Deposited sites](results/neutron_sites.csv)

The repository's previously implemented OBC2 endpoint uses the full protein, native PDBQT charges, solute dielectric 4, solvent dielectric 78.5, and zero surface-area coefficient. It gives positive gaps in **8/8** shell cells for NH₂→OCH₃, with per-receptor medians +0.743 and +0.725. Mean contributions are +0.483 from the direct Coulomb term and +0.276 from the net GB term. The ligand hydration term differs by −1.372 on average; the complex supplies a compensating term. Restricted OBC2 calculations select whole residues having any atom within 6, 10 or 16 Å of either Asp189 Cγ or Ser195 Oγ. These give +0.595, +0.684, +0.709, converging to the full +0.725. This comparison uses identical charge, radius and dielectric parameters for each truncation. Deleting residues also changes the dielectric envelope. Fixed-charge and quantum models differ in extent, charge model, solvation model and the interactions included, so the observed difference cannot be assigned uniquely to any of these. [Original OBC2 implementation and prior independent benchmark](../polar_endpoint/README.md)

Across all three ≥3-fold trypsin edges, double-shell Vina, OBC alone and their unweighted sum give 0, 2 and 1 correct. FXa controls reuse 28 exact-Ki ligands, four seeds, two receptors and both search arms (448 additional true first poses). For identical shell snapshots, Vina/OBC/sum give **234/181/230 of 295** unique within-assay Ki pairs; on 17 equal-heavy-count, equal-charge edges they give **13.5/6/11**. For unguided poses the three totals are **221/177/211 of 295**. The sign rescued in trypsin does not justify adopting a universal weight: the whole-protein polar component omits other contributions, while adding it to Vina can count part of a hydrogen bond twice. The control uses original uncorrected preparation for CHEMBL369480; readers can find its alternate three-dimensional input in the prior benchmark. The estimated energy gaps for FXa use nominal 298.15 K, as the assay temperature is unspecified. [FXa experimental series](https://doi.org/10.1016/S0960-894X(00)00139-6)

## Reproduction

The original pose trees are reproducible from [the directional ion-pair experiment](../directional_ion_pair/README.md). Set ION to its directory containing `trypsin/`, `fxa/`, and `raw/`; FJS and P16 refer to the source FXa preparation directories. Standard Python dependencies are NumPy, SciPy, RDKit (2026.09.2 in this run), OpenMM (8.6.1) and tblite (0.8.0). The independent calculation additionally uses the original Fortran xTB executable (6.7.1). Single precision of some CPU OpenMM internals is checked against the repository's independent double-precision implementation. Floating-point output and hydrogen placement depend on implementation versions.

```bash
D=benchmarks/trypsin_transfer
OUT=/tmp/trypsin-transfer
mkdir -p "$OUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python "$D/quantum.py" "$ION/trypsin" "$OUT/energies.csv" --jobs 4
python "$D/classical.py" "$ION/trypsin" "$OUT/classical.csv"
python "$D/evaluate.py" "$OUT"
python "$D/refine.py" "$ION/trypsin" "$OUT" --jobs 4
python "$D/transfer_control.py" "$ION" "$FJS" "$P16" "$OUT" --jobs 4
python "$D/verify.py" "$ION/trypsin" "$OUT" --xtb /path/to/xtb --neutron "$ION/raw/5MO0.pdb"
python "$D/verify_polar.py" "$ION/trypsin" "$OUT"
```

`results/energies.csv` stores both legs and original pose-file hashes. `angular_raw.csv` records the finer-angle calculations. `cells.csv`, `paired.csv`, and `summary.csv` preserve which finite sampling results enter each count. `regions.csv` and `region_summary.csv` contain all control-cell scores and the stratified FXa/trypsin comparisons. The verification script also re-reads original Vina output and independently sums all reported primary angular partition functions. Sixty single-point energies were recomputed in the original xTB executable from Cartesian XYZ, with a maximum discrepancy of 0.00552 kcal/mol. The double-precision check on all 16 decisive trypsin ligand/cell scores differed from OpenMM by at most 0.000046 kcal/mol.
