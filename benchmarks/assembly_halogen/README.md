# Protein assemblies and halogen contacts in experimental human complexes

Completing a deposited transthyretin (TTR) tetramer rescued a known binding pose of one ligand, VCP-6, on fresh seeds under two search caps. The same operation did not recover thyroxine in a separately determined TTR structure. A local, directional chlorine–carbonyl scoring correction failed to improve top-ranked pose recovery across nine chlorinated ligands. These tests use known, crystallographically localized pockets and fixed protein conformations. Scores and structural success here have no measured affinity endpoint.

## Assembly test

A crystal asymmetric unit may provide only part of a protein–ligand interface. The missing partner changes the sterically accessible volume and supplies interactions during both search and scoring. An electron-density model at a crystallographic symmetry axis can also describe multiple half-occupied orientations inside one physical pocket. Scoring only one orientation creates false structural misses. These two effects were tested separately.

The coordinate audit of all 12 original deposits located every nonidentity crystal-protein contact within 4 Å. Five complexes had such contacts. Two also had nonidentity transforms in their matching deposited biological assembly. For BCL6 7LZQ, 19 ligand atoms contacted crystal copies, and ten contacted the added biological dimer partner. For TTR 8DW5 all 11 contacted ligand atoms could be explained by the added biological chain. The remaining three are MDM2 7BJ6 (two atoms), USP5 7MS7 (nine), and galectin-1 7NML (ten). The assembly preparation uses only the specified biological assembly. Some packing contacts can still be missing from an assembly used for docking. Endpoints in 7LZQ and 8DW5 were frozen before their assembly runs. Six new seeds complemented the original three on these two complexes. A subsequent study specified three further seeds and a new ligand and deposition, TTR–thyroxine 2ROX, before its docking. At 2ROX the asymmetric unit already contains two chains. Building the tetramer supplies a contacting partner for 14 of 24 ligand atoms.

All runs in a comparison share prepared ligand, box (26 Å per side), CPU count (one), three trajectories, nine modes, native Late refinement, and random seed. Each cap is per trajectory. Conformations were generated from the chemical graph after removing the deposited coordinates. The binding-site center and evaluation use the crystallographic ligand. Deposited formal microstates were used throughout. The 8DW5 input contains a neutral carboxylic acid and the 2ROX input contains neutral amine and acid groups. The 8DW5 authors describe a carboxylate contact with Lys15. Neither this model nor the fixed receptor represents a titration or solvent equilibrium. [8DW5 structure and interactions](https://doi.org/10.3390/molecules27217206) · [2ROX source experiment](https://doi.org/10.1107/S0907444996003046) · [7LZQ structure](https://www.rcsb.org/structure/7LZQ)

| Native top-1 success at ≤2 Å | Seeds | Evaluation cap | Asymmetric unit | Biological assembly | Assembly oracle |
|---|---:|---:|---:|---:|---:|
| BCL6–OICR-4425, fixed pose | 9 | 60,000 | 0 | 0 | 1 |
| TTR–VCP-6, fixed pose | 9 | 60,000 | 0 | 1 | 2 |
| TTR–VCP-6, same pocket | 9 | 60,000 | 0 | 2 | 3 |
| BCL6–OICR-4425, fresh seeds | 3 | 60,000 / 240,000 | 0 / 0 | 0 / 0 | 0 / 0 |
| TTR–VCP-6, fresh seeds, fixed pose | 3 | 60,000 / 240,000 | 0 / 0 | 0 / 0 | 0 / 0 |
| TTR–VCP-6, fresh seeds, same pocket | 3 | 60,000 / 240,000 | 0 / 0 | 3 / 2 | 3 / 3 |
| TTR–thyroxine, fresh seeds, same pocket | 3 | 60,000 / 240,000 | 0 / 0 | 0 / 0 | 0 / 0 |

The two entries separated by a slash reuse the same three new seeds at the two different caps. The VCP-6 sites related by the twofold axis have overlapping ligand positions with occupancies 0.46 and centroid separation 1.327 Å. Alternative T4 poses in 2ROX have occupancy 0.50. The same-pocket outcome takes the minimum receptor-frame graph-symmetry RMSD over colocalized official deposited copies within 2 Å centroid separation. The orientation-specific outcome uses only the initially selected copy. For the fresh VCP-6 seeds, median top-1 same-pocket RMSD dropped from 12.46 to 1.17 Å at 60,000 evaluations and from 12.46 to 1.20 Å at 240,000. The median against the single orientation was 2.57 and 2.60 Å in the assembled receptor. BCL6 and T4 had zero oracle hits on every fresh arm. Extra budget did not establish a transferable rescue. The halogenated analogue and hormone are distinct physical complexes on the same protein. Seeds are repeated numerical runs.

Four identity-control and active assemblies were compared against independently downloaded wwPDB biological-unit coordinates, and the later T4 case was checked separately. The 7,322 prepared heavy atoms in the first comparison and another 3,622 for T4 coincide with deposited coordinates at 0.001 Å file precision. Carboxyl-terminal O/OXT names may exchange under preparation. Two monomeric controls (MDM2, BPTF) give six byte-identical pose-pool pairs. A separate implementation parsed PDB BIOMT row records, rebuilt ligand instances, read Meeko conformers, and enumerated graph permutations. All 756 mode RMSDs from the 84 assembly cells reproduced to within 0.000001 Å. Twelve of the 48 initial assembly cells reuse existing native outputs from the halogen test.

## Halogen scoring hypothesis

The directional σ-hole opposite a covalent C–Cl, C–Br or C–I bond can favor a protein carbonyl oxygen along its extension. Auffinger et al. measured element-specific van der Waals sums and characteristic angles, and Koebel et al. introduced element-specific weights and angular factors in VinaXB. We use their structural observations for a deliberately bounded diagnostic. The literature formula and scoring implementation differ from this one. Our test keeps native poses fixed and only reranks them. The original VinaXB study collected known interacting complexes and altered search as well. [Auffinger et al.](https://doi.org/10.1073/pnas.0407607101) · [Koebel et al.](https://doi.org/10.1186/s13321-016-0139-1)

For each covalent C–X donor, the correction subtracts its best matching oxygen's contribution from the printed native score. The radial factor is a Gaussian in the X···O distance with width 0.35 Å, centered at 3.27, 3.37 or 3.50 Å for Cl, Br or I. Those centers are the published vdW sums; using them as radial maxima is a modeling assumption. A cubic smoothstep raises the C–X···O angle factor from zero at 135° to one at 165°. Maximal corrections are 0.265, 0.32, and 0.40 score units. The distance-only control uses identical radii and amplitudes and sets the angular factor to one. Acceptors comprise backbone carbonyls (excluding carboxyl termini identified by OXT) and ASN/GLN side-chain amide oxygen. Saturation is per donor. Solvation, protonation, ligand-dependent polarizability, receptor lone-pair direction and acceptor saturation are unmodeled. No benchmark outcomes were used to tune the constants.

Nine chlorinated ligands, each on a different human target, yielded 243 poses (27 pools). Two fluorinated controls added 54 poses. All nine chlorinated targets produced nine candidates per seed. Native, distance-only, and directional rankings each produced 2/27 ≤2 Å top-1 hits against an 8/27 fixed-pool oracle. The directional term changed one winner in BCL6 (7.07 to 7.23 Å). The distance-only term changed four winners overall. Both controls retained identical scores and winners in the no-donor arms. The selected galectin-1 pocket omitted nearby DMSO and the AKR1B1 control omitted NADP within 3 Å. NOTUM was excluded before docking because preparation removed a proline atom 3.36 Å from the ligand. Excluding all original complexes with a symmetry-related protein contact within 4 Å leaves four chlorinated targets, BRD9, LCN2, SERPINA1 and BPTF. Every method has 2/12 top-1 hits on this subset, against an oracle of 6/12.

LCN2 6Z2C provides a direct counterexample to the sufficiency of a favorable short interaction. All three pools sample a near-native mode (0.50–0.88 Å); the selected modes are 2.88–2.90 Å. In two seeds the directed relative correction favors the near-native conformation by 0.158 and 0.163 units, below the native score gaps of 0.750 and 0.951. On the remaining seed it favors the decoy by 0.183. The deposited ligand has six chlorines, and both poses can retain at least one favorable contact. In the independent 3DY7 MEK1 positive control, C6–I1···O(Val127) is 177.810° and I···O is 3.3542 Å, reproducing Koebel's raw geometry; this case was a crystal geometry check. The retrospective fixed-pool challenge itself consists of chlorine and fluorine chemistry. Br and I did not enter this primary panel.

Two processing paths recovered mode RMSDs and corrections for all 297 poses to within 0.000001 Å and 0.000001 score units. Constraining graph symmetries to preserve stereochemistry changed some individual LCN2 values by up to 0.703 Å. No pose crossed the 2 Å threshold. A default binary replay on all 33 pools reproduced the original pose files byte for byte.

## Files and replay

The archive contains the prepared ligands, receptors, original PDBs and inference coordinates, official biological-assembly files, original mode pools, preparation and docking logs, and the 3DY7 control. Data were drawn from the 428-case version of the PoseBusters benchmark and from RCSB. This is a manually selected subset based on chemistry and human target identity, with one pocket per protein. It is not an evaluation of PoseBusters as a whole. The original PoseBusters inputs are attributable to Martin Buttenschoen, Garrett Morris and Charlotte Deane. [Dataset](https://zenodo.org/records/8278563) · [Paper](https://doi.org/10.1039/D3SC04185A)

Run from the repository root using Python 3.11, NumPy 2.4.6, SciPy 1.17.1, RDKit 2025.03.5, Meeko 0.6.1 and Gemmi 0.7.3. The assembly CLI requires a deposited PDB with REMARK 350 and an explicitly selected biological assembly. Input rigid receptor chains must be aligned to their original author-frame coordinates, and each required chain must be present. Prepared cofactors on included chains follow their chains. Crystallographic waters, omitted cofactors and unmodeled residues in the input cannot be recovered by assembly generation. The original-frame search-box definition remains usable without ligand coordinates. The halogen CLI requires complete Meeko SMILES index records and written MODEL blocks.

```bash
B=benchmarks/assembly_halogen
W=/tmp/fastdock-assembly-replay
mkdir -p "$W/results"
(cd "$B" && sha256sum -c fixtures.sha256)
tar -xJf "$B/fixtures.tar.xz" -C "$W"
python tools/build_receptor_assembly.py "$W/data/8DW5_FQ7/deposition.pdb" \
    "$W/data/8DW5_FQ7/receptor.pdbqt" "$W/tetramer.pdbqt" --assembly 1
python tools/rank_halogen_pool.py "$W/data/6Z2C_Q5E/receptor.pdbqt" \
    "$W/data/6Z2C_Q5E/poses-314159.pdbqt"
python "$B/analyze.py" "$W/data" "$W/results"
python "$B/verify.py" "$W/data" "$W/results" --binary build/vina --replayed "$W/replayed"
python "$B/environment.py" "$W/data" "$W/3DY7.pdb" "$W/results"
python "$B/analyze_assembly.py" "$W/data" "$W/assembly" "$W/results"
python "$B/analyze_assembly.py" "$W/data" "$W/assembly" "$W/results" --high
python "$B/verify_assembly.py" "$W/data" "$W/assembly" "$W/results"
```

For fresh assembly docking, point the two invocations of `run_assembly.py INPUT NEW_DEST BINARY` (one with `--high`) at the extracted `data` directory and a new destination. Copy the official `deposited-assembly.pdb.gz` files for the five targets from the archive to their matching new output subdirectories before the subsequent analysis. Rebuilding the original nine-ligand + two-control native pools uses `prepare_dock.py EXTRACTED_POSEBUSTERS_DIRECTORY WORK BINARY --redock`. Logs document the excluded twelfth choice. The previously archived default binary SHA-256 was `3b09a3b70c27e39b872f69479cddd14e979cc509cfbe09f83d52302cee7806f0`; the separate default replay fingerprint is recorded in `results/audit.json`.
