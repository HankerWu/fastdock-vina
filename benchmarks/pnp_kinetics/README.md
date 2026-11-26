# Initial and slow-onset human PNP inhibition constants

Seven immucillins from a single purified human purine nucleoside phosphorylase study have separate rate-based initial Ki and equilibrium Ki* at 25 and 37 °C. The authors fit initial rates and later steady-state rates to a competitive inhibition equation, with substrate concentration and depletion of free inhibitor included. Assays used 1 nM PNP, 50 mM phosphate pH 7.5, 1 or 10 mM inosine depending on inhibitor, and a coupled spectrophotometric readout. Inhibitors were studied in triplicate. The primary endpoint for this panel is initial Ki at 25 °C. Threefold comparisons for each other endpoint are reported separately. [Original experiment, Table 1 and Eqs. 1–2](https://pmc.ncbi.nlm.nih.gov/articles/PMC5612120/)

| Temperature and endpoint | Eligible unique pairs | Native | Balanced torsion | Prepared torsion starts |
| --- | ---: | ---: | ---: | ---: |
| 25 °C, initial Ki | 9 | 9 | 5 | 8 |
| 25 °C, slow-onset Ki* | 9 | 8 | 6 | 7 |
| 37 °C, initial Ki | 15 | 11 | 13 | 7 |
| 37 °C, slow-onset Ki* | 10 | 9 | 7 | 7 |

Each number uses the median of four matched top-pose scores per ligand. The primary two-seed halves give 9/4/8 and 8/8/8 out of nine for native/balanced/prepared. Excluding the cocrystallized ImmH leaves four primary pairs with 4/3/4 correct. Changing the endpoint or temperature can change both the experimental order and which pairs meet the threshold. The underlying docking scores are identical for these sensitivity analyses. Close ligand analogues and repeated measurements limit independence. [All pairs and gaps](pairs.csv) · [Cohorts and exclusion sensitivity](summary.csv)

The pocket is the ImmH site in the phosphate-containing 1RR6 human biological trimer. Three crystallographic copies are included with separate chain identifiers. The model retains all deposited phosphate heavy atoms as rigid P/acceptor-O sites, omits water, and uses the deposited blank or A alternate positions. PDBFixer required no missing residue-sidechain heavy atom additions to the observed protein. Receptor centers use the known ImmH coordinates. Every ligand starts from a sequence of graph-derived conformers prepared independently of bound coordinates with the branch ETKDGv3/MMFF/Meeko workflow. The four H-series chemical graphs come from RCSB Chemical Component Dictionary entries IMH, DIH, 229 and 22A. Three G-series counterparts carry the C2 amino substitution shown in the authors' Figure 1. Each sugar or aminoalcohol amine is modeled protonated, following the bound cationic transition-state analogue interpretation. A 22 Å cube, two trajectories, 25,000 evaluations per trajectory, nine modes, one CPU and four seeds are held constant. [Primary trimer](https://www.rcsb.org/structure/1RR6) · [Four generations of structures and catalytic phosphate interactions](https://pubmed.ncbi.nlm.nih.gov/20212140/)

This model tests native endpoint scoring in one rigid, inhibitor-stabilized state. The trimer exhibits binding-site cooperativity, its phosphate–cation interaction is important, and the native Vina objective uses coarse atom-type potentials. The phosphate treatment omits its protonation ensemble, donor geometry and explicit Coulomb energy. The structural template is closer to a relaxed inhibited complex than to an initial encounter ensemble. A fixed score cannot distinguish initial from equilibrated inhibition or predict its temperature dependence.

The 84 search pools contain 739 modes. Reading the independent stdout rank tables reproduces all 129 arm–pair decisions across the four endpoints. To rerun, extract [evidence.tar.xz](evidence.tar.xz) into a fresh directory, then run the command below with a compiled two-switch binary and the same Python environment as the Galectin panel. The archive includes all source inputs, prepared files and raw logs. Match original source bytes against [source_sha256.json](source_sha256.json). RCSB 1RR6 uses the decompressed biological assembly 1 PDB file; the four chemical records are the RCSB core chemcomp JSON for their component codes.

~~~bash
OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 python benchmarks/pnp_kinetics/run.py "$WORK" "$VINA" --jobs 3
~~~
