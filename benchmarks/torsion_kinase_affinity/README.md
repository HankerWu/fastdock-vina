# Matched experimental Kd on CHEK1 and Aurora A

The outcome is ordering exact dissociation constants from Davis et al., Supplementary Table 4, within one unmodified kinase assay. Table 1 confirms the identifiers and wild-type labels. Aurora A has 20 measured compounds, CHEK1 has 14, with 25 compounds in the union. Within each target every compound has a distinct typed Bemis–Murcko scaffold. The 217 target-specific pairs separated by at least threefold share molecules. All labels are measured Kd. [Original assay and supplements](https://doi.org/10.1038/nbt.1990) · [Experimental labels and chemical identities](panel.csv)

## Preparation and evaluation

Each ligand was prepared from its molecular graph with ETKDGv3, eight trials and a lowest-converged-energy MMFF input. The binding parent was retained for R406 and JNJ-28312141 salts. Neutral input graphs were used throughout. Protonation, tautomer populations, water and state populations lack experimental calibration here. The Aurora A states are 2X6E and 2J4Z; CHEK1 uses 2E9P and 2YWP. The prospective structure screen discarded 7BJJ, a ten-point mutant. AURKA 2X6E required six incompletely observed side chains to be modeled, with a closest repaired atom 4.05 Å from the box-defining ligand. Thirteen distal incomplete side chains were truncated to Ala, retaining their observed backbone and CB, where all observed atoms of the residue lay beyond 15 Å from that ligand. In 2J4Z the nearest repaired atom is 17.85 Å from the ligand. [2X6E](https://www.rcsb.org/structure/2X6E) · [2J4Z](https://www.rcsb.org/structure/2J4Z) · [2E9P](https://www.rcsb.org/structure/2E9P) · [2YWP](https://www.rcsb.org/structure/2YWP)

The three paired arms used native Late refinement, the same 26 Å box per receptor, four seeds, two trajectories, a 25,000-evaluation cap per trajectory, nine requested modes and one CPU. The 816 runs returned 7,237 modes. The primary consensus averages the two state-specific medians over seeds. Score ties receive half credit.

| Target | Unique ≥3× pairs | Native | Balanced proposal | Prepared starts | Heavy-atom count |
|---|---:|---:|---:|---:|---:|
| CHEK1 | 66 | 42 | 43 | 38 | 31 |
| AURKA | 151 | 110 | 106 | 103 | 96.5 |

Balanced gains two pairs and loses one on CHEK1. On AURKA it gains two and loses six. Whole-compound resampling gives a descriptive 95% interval for its change in correct fraction of −4.3 to +11.8 percentage points for CHEK1 and −9.9 to +3.7 for AURKA. Across all measured values the native / balanced / prepared Spearman coefficients are 0.288 / 0.304 / 0.145 on CHEK1 and 0.565 / 0.507 / 0.474 on AURKA. Under an equal-population exponential score average the counts become 46 / 42 / 43 and 113 / 110 / 105. Using a common ligand-specific free reference from the prepared conformer gives 37 / 36 / 32 and 109 / 99 / 97 with arithmetic averaging. [Cells](cells.csv) · [Individual states and sensitivity](summary.csv) · [Correlations and bootstrap](uncertainty.csv)

The first computational panel embedded one hand at the unspecified stereocenter of Ki20227. Its original authors explicitly distinguish racemic Ki20227 from its R and S enantiomers. We completed 48 additional searches with both hands explicitly specified. For a 1:1 racemic stock with one competitive site, low ligand depletion and fixed free fractions, the apparent association constant averages the two enantiomer association constants. Applying the corresponding log-sum-exp at 298.15 K within each receptor state yields 111 / 106 / 103 of 151 on AURKA, including 13 / 12 / 11 of 15 pairs involving the racemate. Omitting the racemate yields 98 / 94 / 92 of 136. The docking scores only approximate the energy quantities in this relationship. [Chemical definition](https://aacrjournals.org/mct/article-pdf/5/11/2634/1871659/2634.pdf) · [R/S results](mixtures.csv)

An independent evaluator read the author tables, molecular graphs and 816 raw mode pools. A separate C++ pass rescored every serialized top pose with its original prepared topology and unbound reference. The maximum difference from the printed score was 0.014 kcal/mol. Every arm's four-seed arithmetic-consensus pair signs reproduced. A second source matched 28/34 target–compound entries by CID and parent. BIBF-1120 (derivative) and HKI-272 were outside the selected compounds. Coverage is 20/20 exact AURKA measurements and 14/16 exact CHEK1 measurements from these rows. [Independent ranking](independent/independent_ranking.csv) · [Identity checks](independent/independent_chemistry.csv)

### Reusing the receptor coordinates

The original 2X6E coordinate intermediate omits residues 281–291 and 303–305. OpenMM reads its continuous chain as two additional C–N bonds of 20.723 and 9.342 Å. The current preparation reports the missing sequences from the deposition and separates the observed polymer segments in its coordinate export. The archived docking inputs and results used a rigid PDBQT. A flexible or solvated model needs an explicit treatment of the internal sequences and chemistry. [Topology, matched preparation runs and pose evidence](../receptor_linkage/aurora/internal_loops/README.md)

## Replay

The archive contains prepared PDBQT, receptor preparation records, all 864 raw runs including the enantiomers, and pair-level sensitivity. Run from the repository root with Python 3.12, RDKit 2026.09.2, Meeko 0.8.0 and a Boost build of Vina. Data extraction and the independent check also use xlrd.

```bash
B=benchmarks/torsion_kinase_affinity
W=$(mktemp -d)
(cd "$B" && sha256sum -c evidence.sha256)
tar -xJf "$B/evidence.tar.xz" -C "$W"
cp "$B"/{names.json,panel.csv} "$W/"
python "$B/run.py" "$W" build/vina --jobs 4
python "$B/analyze.py" "$W" build/vina --jobs 3
python "$B/uncertainty.py" "$W"
python "$B/racemate.py" "$PWD" "$W" build/vina
```

For source preparation use the previous benchmark's `benchmarks/kinase_state_affinity/source/fetch.py DAVIS` to obtain the source sheets, identity records and secondary data. Extract the structural benchmark's `fixtures/inputs.tar.xz` into POSE. The two raw Aurora A PDBs are in the evidence archive. Run the following with DAVIS and POSE set to those directories.

```bash
python "$B/prepare_receptors.py" "$W"
python "$B/data.py" "$W" "$DAVIS/raw" "$PWD" "$POSE"
python "$B/independent/verify.py" "$W" "$DAVIS/raw" "$W/independent-check"
```
