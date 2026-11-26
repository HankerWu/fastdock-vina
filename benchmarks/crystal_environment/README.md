# Crystal environment, alternative copies and transfer of experimental ligand shape

## Questions and decisions

Can contacts with other protein molecules in the crystal explain incorrect poses of flexible three-tailed hCA II sulfonamides? Can local diffraction density distinguish the two declared protonated-amine geometries? Can a protein-aligned different ligand provide a transferable shape term that improves on the production Gaussian overlap in `src/lib/pose_template.h`? Predictions were tested on the same saved candidates within each search arm. Query-bound maps and coordinates were also scored as explicit retrospective controls. The tables identify these as `own_*`. The candidate lists, docking scores and measured chemical pairs stayed fixed.

The data support local crystal-environment and donor effects, including strong selection within a narrowly related BRD4 pair. Normalizing volume or allowing anisotropic disorder did not improve pose selection and same-target chemical ordering consistently across three targets. We retain the benchmark and keep the existing production scoring.

## Real coordinate pools and endpoints

| Panel | Ligands | Receptors | Seeds | Compared poses | Main selection pool |
|---|---:|---:|---:|---:|---|
| hCA II | 5 | 2 | 2 | 359 | common and reflected prepared amine starts combined |
| BRD4 BD1 | 4 | 3 | 4 | 1,272 | each of native, one correspondence, multiple correspondences separately |
| T4 lysozyme L99A | 9 | 2 | 4 | 1,733 | each of the same three search arms separately |

The third panel tests shape alone. In hCA II the joint endpoint requires symmetry-corrected heavy-atom RMSD ≤2 Å and agreement with the deposited sign of the prepared three-carbon environment of the central nitrogen. The heavy-atom sign is a model-derived geometric label. The input microstate has a sulfonamide anion and a tertiary ammonium with three distinct attached carbons and one H. Its SMILES leaves the amine hand unspecified. The scaffold also contains a tertiary amide nitrogen. Misidentifying that nitrogen would test a different physical question. The two independently prepared starts do not establish a measured inversion barrier, pKa or equilibrium population. BRD4 and T4L use the established heavy-atom RMSD endpoint. Receptor and seed cells repeat molecules. Chemical ordering is counted only within a target/assay, with a threefold measured difference, one-half credit for ties. Unique-pair counts collapse median signed gaps across receptor/seed cells. All selections retain the solution-facing native affinity at the selected pose. Physical comparisons are conditional on the fixed microstates, frozen receptors and existing search arms.

### Crystal packing and reference copies

`packing.py` enumerates crystallographic space-group operations and all lattice translations whose transformed asymmetric units could come within 18 Å of the query box. This bound comes from the reciprocal basis of the unit cell. It includes copies from the identity operation with a nonzero translation and different chains of the reference asymmetric unit. Only the host protein in the identity position is excluded. All deposited nonzero-occupancy protein heavy-atom conformers contribute. The test adds the native Vina-style repulsion coefficient 0.840245 times the occupancy-weighted squared overlap below the sum of elemental XS radii. Water, ions and other small molecules are outside this particular packing test. The receptor lattice is held fixed for every ligand docked to that receptor. Crystal repulsion is used to select poses, never to estimate solution binding energy. Other possible packing mechanisms, including attractive electrostatics and side-chain reorganization, are outside this test.

The five active hCA II ligands (6WQ4/5/7/8/9, residue A303) all have unit modeled occupancy, a single set of active-site ligand coordinates and whole-ligand RSCC 0.956–0.974. The nearest other protein is always the screw-related mate with operation `−x,y+1/2,−z` and lattice translation (−1,0,0). Distances are 3.03, 2.65, 3.37, 3.22 and 3.12 Å. The neighboring residues are Glu236/238. 6WQ8 A304 is a second ligand at a separate location. Its nearest ligand-to-catalytic-Zn distance is 16.47 Å and RSCC 0.808, compared with 1.91 Å and 0.963 for A303. Only A303 evaluates an active-site pose. BRD4 has ten experimental copies across its four query structures. Each copy was transferred into the docking receptor frame using matched protein Cα atoms at positions 44–166. The input chemistry and symmetric heavy-atom permutations were preserved.

| Main pool, original reference | Native | Native + receptor packing | Mapping prior | Mapping prior + packing | Any good candidate |
|---|---:|---:|---:|---:|---:|
| hCA II, joint, 20 cells | 0 | 0 | 6 | 5 | 10 |
| BRD4, multi, 48 cells | 19 | 17 | 21 | 18 | 24 |
| BRD4, multi, any legitimate crystal copy | 19 | 18 | 23 | 21 | 24 |

The change in reference copy at most recovers a few model selections in BRD4; it does not recover the missing ligands in the multi pool. In the native BRD4 search arm, admitting the additional copies changes oracle coverage from 2 to 4 of 48. The hCA II native selector before/after packing preserves 13/28 comparable pair decisions. BRD4 multi preserves 25/60. Neither result establishes a solution affinity effect of crystal packing.

### Map support and donor selection

`density.py` uses the wwPDB validation map coefficients (`pdbx_FWT`, `pdbx_PHWT`), reconstructs 2mFo−DFc maps and normalizes each whole-cell grid to zero mean and unit variance. Local values are interpolated at heavy atoms (not atom integrals, independent electron counts or atom-specific confidence intervals). The maps have model-derived phases and are separate from the paper's omit maps. In 6WQ7 (compound 42) N15 samples −0.00025σ at B=48.47 Å². Its three covalently bonded carbons C13/C16/C17 sample 1.16/0.37/0.46σ. The whole ligand nevertheless has RSCC 0.959. The other four central nitrogens sample 3.22/1.39/1.74/0.93σ in 6WQ4/5/8/9 respectively. Whole-ligand scores alone give little information about the proposed local geometry in 42. The paper independently reports weak omit density at the ends of several tails. This does not determine the configuration, alternate occupancy, protonation state or free-energy cost. All scored heavy atoms remain in the RMSD.

The preset donor comes from `benchmarks/mapping_ensemble/results/*/graphs.json`. The original graph-only rule ranks non-query ligands by the number of common substructure atoms, then common fraction, then PDB ID. hCA II donor identity here inherits the historical graph constructed from the independently embedded starting ligand, which can carry an incidental amine hand. The selection never uses the deposited query ligand position. A separate control takes the minimum Morgan radius-2, 2048-bit Tanimoto among remaining ligands from the same series after removing both the query and preset donor, breaking ties by PDB ID. For 4O76 (TG101209), the preset/control Tc values are 0.681/0.057. For hCA II all control Tc values are at least 0.562 and several panels have an alternative donor with a higher Tc than the preset MCS donor. The column `distant` denotes this within-series alternative, without implying universally unrelated chemistry. A truly unliganded map (hCA II 3KS3, BRD4 2OSS) gives an independent negative control for density. The two hCA II docking receptors include cognate cases; the strictly cross-docked results below exclude these.

We transfer each pose into the donor frame by protein Cα superposition. hCA II uses all shared residue numbers, BRD4 uses 44–166, and T4L uses 80–160 excluding flexible 107–115. No ligand coordinates determine this registration. The mean density over predicted heavy-atom positions can include protein, ions or ordered water. The masked version keeps only positive grid density within 2 Å of a donor ligand atom and closer to that atom than to any donor ASU protein/water/metal atom. Both are schematic map scores. A coordinate-only control uses mean exp(−d²/2) with d the closest donor heavy atom distance in Å.

| Main pool, pure selection | Query's own bound map | Other ligand's map | Other ligand's masked map | Unliganded map | Other ligand's coordinates |
|---|---:|---:|---:|---:|---:|
| hCA II, 20 | 5 | 7 | 7 | 4 | 9 |
| BRD4, 48 | 20 | 21 | 22 | 11 | 23 |
| hCA II, strictly cross, 16 | — | 4 | 4 | 4 | 6 |
| BRD4, strictly cross, 40 | — | 19 | 19 | 10 | 20 |

The entire coordinate/map experiment uses the same 1,631 poses in this table. A sampled full-cell density at a few points does not distinguish ligand from protein scattering by itself. An intentionally leaked query-bound density score is a diagnostic control rather than a best-achievable upper bound; for an actual upper bound use the oracle.

### Shape and size normalization

The existing C++ score averages the two heavy-atom nearest-neighbor Gaussian directions at σ=1.5 Å using same-element matches, including zero contributions for missing elements. The ablation untypes the atoms, changes directionality, or changes width to 1 Å. A separate analytic approximation uses a sum of fixed-integral Gaussian atom densities. An element has sphere volume M=4πR³/3 from its Bondi radius and isotropic covariance Σ=M^(2/3)/(2π) I, giving unit peak height. For query and donor atoms the product integral has exponent −Δᵀ(Σq+Σd)⁻¹Δ/2 and coefficient proportional to Mq Md / √det(Σq+Σd). Summing atom pairs gives Iqd. We normalize by Iqq+Idd−Iqd. Additive atomic density and fixed protein-frame overlay differ from a Gaussian union or a freely optimized shape alignment. The displacement sensitivity test broadens the donor covariance by its deposited ANISOU tensor when available, or by B/(8π²) I, preserving atom-integrated volume. These displacement parameters describe the crystal model's thermal/static spread; they do not supply a statistically independent coordinate-error estimate or locate an absent conformer.

| Selector on main pools | hCA II /20 | BRD4 /48 | T4L /72 |
|---|---:|---:|---:|
| Native | 0 | 19 | 58 |
| Existing typed bidirectional overlap, alone | 9 | 22 | 49 |
| Untyped one-way σ=1 Å, alone | 9 | 23 | 47 |
| Normalized volume, alone | 8 | 23 | 47 |
| Normalized volume with deposited U, alone | 8 | 23 | 47 |
| Native − 6 × existing overlap | 2 | 22 | 57 |
| Native − 6 × one-way | 2 | 23 | 61 |
| Native − 6 × normalized volume | 2 | 24 | 58 |
| Native − 6 × normalized volume with U | 2 | 24 | 57 |
| Alternative donor, one-way alone | 5 | 14 | 46 |

Only these comparisons hold the search and printed candidates fixed. The energy column here refers to the raw ranking objective before each ligand's common torsion scaling. It tests the selection term, not the production select/refine dynamics. Under the stricter condition in which neither the donor nor docking receptor is the query crystal, the native scores give 0/16, 16/40 and 54/68. The preset/alternative one-way coordinate scores give 6/3, 20/10 and 43/43 hits in those same panels. Raw native − 6 × normalized volume gives 2/16, 20/40 and 54/68, and including the donor's U gives 2/16, 20/40 and 53/68.

BRD4's 23/48 geometric hits comprise 12/12 for 4O76 and 11/12 for 4OGJ. Neither 4O77 nor 4OGI has a good candidate among the 24 corresponding cells. In hCA II, 6WQ7 and 6WQ9 are the missing ligands; each has 0/4 oracle cells. A change to scoring can affect pose choice only where sampling already produced a good candidate. At matched conditions the preset one-way geometry-only selection gives 11/28, 16/60 and 113/152 signed pair credits, or 1/7, 1/5 and 14/19 after taking the median signed gap per distinct pair. These denominators distinguish search replicates from independent chemistry. The pose improvement in the two related BRD4 inhibitors fails the generalization test.

## Reproduction and independent checks

Requires Python, NumPy, SciPy, RDKit, Meeko, Gemmi, and a C++17 compiler with Boost headers. All scripts accept positional directories. Obtain the original prepared inputs and candidate outputs from the matching mapping ensemble evidence (`hcaii_declared` two-state run, BRD4, T4L). Here `CA`, `BRD`, `T4` each contain `raw/` and `prepared/`; the `*_RUN` directories contain the unchanged PDBQT models and benchmark tables. `OUT` is a new result directory. Validation XML, independent RCSB ligand-instance SDFs and mmCIF validation coefficients download on first use and are cached. Raw-source checksums use paths relative to the crystal evidence root. The saved CSVs contain individual candidates, selected model IDs, receptor/seed cells, donor assignments, chemical pair credits, per-atom densities and cognate/cross splits.

```bash
export OPENBLAS_NUM_THREADS=1
cd benchmarks/crystal_environment
python evaluate.py "$CA" "$BRD" "$CA_RUN" "$BRD_RUN" "$OUT"
python density.py "$CA" "$BRD" "$CA_RUN" "$BRD_RUN" "$OUT" "$MAPS"
python shape.py "$CA" "$BRD" "$T4" "$CA_RUN" "$BRD_RUN" "$T4_RUN" "$OUT"
python split.py "$OUT" "$T4_RUN"
python verify.py "$CA" "$BRD" "$CA_RUN" "$BRD_RUN" "$OUT" "$INSTANCES"
python verify_maps.py "$MAPS"
python verify_shape.py "$CA" "$BRD" "$T4" "$CA_RUN" "$BRD_RUN" "$T4_RUN" "$OUT"
```

The independent packing path uses deposited PDB SMTRY matrices, CRYST1 and explicit translation enumeration from −3 to 3 in each dimension. Independent BRD4 RMSD uses ten RCSB instance SDFs, independently computed receptor registrations and RDKit graph symmetry. The reciprocal-space validation applies space-group matrices and Friedel mates directly to measured map coefficients, constructs the reciprocal grid and takes the NumPy FFT. The C++ adapter runs the repository's actual overlap method over every donor/candidate pair. See `results/verification.txt`.

## Primary sources

* Bonardi et al., three-tailed sulfonamides, enzyme inhibition and deposited hCA structures. [Journal article](https://doi.org/10.1021/acs.jmedchem.0c00733)
* Ciceri et al., BRD4(1) biochemical profiling. [Journal article](https://doi.org/10.1038/nchembio.1471). 4O76 coordinates originate from a separate structural publication. [Original study](https://doi.org/10.1021/cb500072z)
* Morton et al., matched T4 lysozyme thermodynamics. [Journal article](https://doi.org/10.1021/bi00027a006). [Deposited L99A–benzene structure](https://www.rcsb.org/structure/181L)
* Grant, Gallardo and Pickup, analytic molecular Gaussian shape comparison. [Journal article](https://doi.org/10.1002/%28sici%291096-987x%2819961115%2917%3A14%3C1653%3A%3Aaid-jcc7%3E3.0.co%3B2-k)
* Liebschner et al., model bias and solvent effects in ligand OMIT maps. [Journal article](https://doi.org/10.1107/S2059798316018210)
