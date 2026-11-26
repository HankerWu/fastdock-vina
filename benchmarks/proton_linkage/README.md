# Binding-linked protonation in measured carbonic anhydrase series

A change in solution pH can change conditional affinity even when the bound pose stays fixed. Adding the full cost of selecting a ligand anion to an empirical docking score is a separate hypothesis. These tests use observed dissociation constants, independently measured ligand pKas, paired docking runs and buffer-dependent enthalpies. All comparisons are retrospective. No coefficient was trained on their affinity labels.

## A binding polynomial with a falsifiable bound

For a free monoprotic sulfonamide and one Zn-bound water site, let L be the acidic ligand and R the receptor. If the binding-competent species are the ligand anion and the receptor with water bound to Zn, their solution fractions are

```math
f_L=\frac{1}{1+10^{pK_L-pH}},\qquad f_R=\frac{1}{1+10^{pH-pK_R}}.
```

By mass action, the observed association constant is the intrinsic association constant multiplied by both fractions. In free energies this gives an additive cost of $-RT\ln(f_L f_R)$. With a measured reference affinity at the same temperature, the unmeasured intrinsic affinity cancels.

```math
K_d(h_2)=K_d(h_1)\frac{f_L(h_1) f_R(h_1)}{f_L(h_2) f_R(h_2)},\qquad
\frac{d\log_{10}K_d}{dh}=f_L-f_R.
```

The derivative is the net uptake of protons per association. Its magnitude is bounded by one. A shift of 0.5 pH unit can change the conditional affinity by at most $10^{0.5}=3.16$ fold under these assumptions. When ranking different ligands for the same receptor at a single pH, the receptor factor cancels. When comparing poses of one ionization state both factors cancel. Additional populated states require their own free and bound partition functions, with each physical bound state counted once. A common solution-state constant cannot distinguish two poses with reversed donor directions. These relationships assume equilibrated, independent free sites, the same dominant bound state, and unchanged intrinsic affinity across the conditions.

## Acidic transfer in human CA IX

Vaškevičius et al. measured the five ligand pKas by UV/vis titration at 37 °C. Figure S24 gives the numbers in the table below. Both Kd columns in supplementary Table S2 were measured by fluorescence-based thermal shift assay (FTSA), with binding constants reported at 37 °C. Only the measured pH 7 values enter the prediction. The paper's intrinsic Kd columns are back-calculated from the observed Kds and are excluded throughout. Protein pKa 7.01 is the catalytic-domain kinetic titration inflection from an independent study by Alterio et al. The pH 7 binding assay used phosphate and 100 mM NaCl. The pH 5 assay used a phosphate/acetate/borate mixture and 50 mM NaCl. The result includes any accompanying buffer and ionic-strength response. [FTSA and supporting information](https://doi.org/10.1021/acs.jmedchem.5c01142) · [Figure S24 and Table S2](https://doi.org/10.1021/acs.jmedchem.5c01142.s001) · [Independent protein titration](https://doi.org/10.1073/pnas.0908301106)

| Ligand | Ligand pKa | pH 7 reference Kd, nM | pH 5 measured Kd, nM | pH 5 predicted Kd, nM |
|---|---:|---:|---:|---:|
| 4c | 7.44 | 15 | 330 | 564 |
| 5c | 8.18 | 0.022 | 1.1 | 1.055 |
| 6 | 8.70 | 0.012 | 0.67 | 0.601 |
| 13 | 8.29 | 0.0046 | 0.17 | 0.224 |
| 14 | 8.44 | 0.0045 | 0.25 | 0.222 |

Mean absolute error across five ligands was **0.0938 log10 units**, against 1.6201 for carrying the pH 7 reference unchanged, 0.3387 for a ligand-only correction, and 0.1603 when predicting each shift from the median log shift of the other four ligands. Pairwise differences of log shifts cancel any shared receptor or common buffer contribution. Their mean absolute error was 0.1471 over the ten dependent pairs, versus 0.1972 for any single common shift. The pKas span only 1.26 units and these five molecules share a scaffold. For protein pKa 6.6 reported in earlier CA binding work, the absolute error increased to 0.1961. Using 6.49, measured for a construct including the proteoglycan domain, gave 0.2744. The kinetic inflection is a surrogate and its construct, temperature and assay matter. The series does not establish broad cross-target calibration or new-ligand prediction. [CA isoform pKas](https://doi.org/10.3390/ph15040477) · [Construct dependence](https://doi.org/10.1073/pnas.0908301106)

## Failure at alkaline pH

In the same 2025 paper, compound 13 has a pH 10 direct ITC Kd of **26.6 ± 7.50 nM**, read from Figure S25D. Transfer from the 0.0046 nM pH 7 FTSA measurement gives **0.113 nM**, approximately 235-fold too tight. Its experimental ratio of 5,783 exceeds the largest three-pH-unit change of 1,000 allowed by this model. Retaining the measured ligand pKa tightens the maximum ratio to 49.74 for **any** receptor pKa. Here the readout and buffer also change. [Figure 1 and methods](https://doi.org/10.1021/acs.jmedchem.5c01142) · [Figure S25](https://doi.org/10.1021/acs.jmedchem.5c01142.s001)

The two alkaline **same-ITC**, same-37 °C, same-Tris/NaCl fits provide a sharper test. Panel C gives 1.04 ± 0.684 nM. Its graphic labels pH 9.5, its caption says 9.6, and the methods say 9.7. Panel D labels pH 10.0 consistently. Carrying all three readings forward gives predicted final Kd 3.15, 2.54, or 2.03 nM. The observed ratio is 25.58. The most permissive of the three two-site bounds is 3.16. Moving both fitted values against the effect by one printed fit error still gives 11.08. These fit errors are not replicated confidence intervals. At 10 µM protein the tighter point has an ITC concentration-to-Kd ratio near 10,000, with correspondingly poor resolution. The experiment/model mismatch survives varying the single receptor pKa. Additional coupled sites or changes in intrinsic affinity and systematic assay error remain distinguishable explanations in principle. This test does not identify which occurred. [Supplementary calorimetry](https://doi.org/10.1021/acs.jmedchem.5c01142.s001)

## Seven ligand docks in bovine CA II

Krishnamurthy et al. measured all seven fluorobenzenesulfonamide free-ligand pKas by potentiometry (reported uncertainty ±0.2 pH units) and conditional Kd by ITC at pH 7.5, 298 K in 20 mM sodium phosphate with 0.6% D₂O. The Kd uncertainties in Table 2 are maximum deviations from the average of three measurements. Our energy conversion uses 298.15 K for 25 °C. All seven compounds are retained, including 4-F, discussed separately in the source's fit. Only directly observed Kds from Table 2 are labels. [Tables 1–2 and methods](https://doi.org/10.1002/asia.200600360)

Each ligand's monoanion was prepared from its graph with RDKit ETKDGv3 and MMFF, then Meeko. The 20 Å receptor boxes use the three Zn-coordinating histidines, independent of ligand coordinates. For 1V9E the three closest histidine donors in chain B are N at 2.13, 2.13 and 2.18 Å. Chain A has two N contacts and CD2 at 2.27 Å. Chain B labels are remapped to the existing preparation template. 6SKT chain A has three N at 2.01–2.06 Å. Crystallographic solvent and pre-existing ligands are removed, retaining Zn. These are distinct coordinate sets of the same protein. 6SKT used 1V9E as its molecular replacement model. [1V9E](https://www.rcsb.org/structure/1V9E) · [6SKT experimental details](https://www.rcsb.org/experimental/6SKT)

Two receptors × seven ligands × two seeds (314159, 271828) × native/point-guided search give 56 docks. Each uses four trajectories, max_evals 50,000 per trajectory, one CPU and up to nine output modes. The point-guided arm directs the sulfonamide donor towards the Zn tetrahedral vacancy at 2 Å during search. Both arms use the reported top pose's Vina affinity. Gaps are computed **before** taking the median across the four matched receptor/seed conditions. Each experimental pair with weaker/stronger Kd ≥3 counts once per arm. All 15 pairs share seven drug measurements.

| Scoring model | Native correct / 15 | Native gap MAE, kcal/mol | Point-guided correct / 15 | Point-guided gap MAE, kcal/mol |
|---|---:|---:|---:|---:|
| Vina | 14 | 0.760 | 14 | 0.700 |
| Vina plus ligand population cost | 13 | 0.931 | 13 | 0.872 |
| Ligand population cost alone | 13 | 0.609 | 13 | 0.609 |
| Both terms plus Zn surrogate slope | 13 | 0.723 | 13 | 0.663 |

The population-only metric can benefit from covariation of fluorination, pKa and hydrophobicity. At fixed fluorine count, each search arm drops from 2/3 to 1/3 after the population correction. For 3-F (75 nM) versus 2-F (230 nM), the native median gap favouring 3-F is +0.0395 kcal/mol, and the population contribution is −0.1355 kcal/mol. Their pKas differ by only 0.1, smaller than the measurement uncertainty. On the stricter subset in which the ratio stays ≥3 under both reported Kd deviations, all scoring models get all 10 pairs correct. Endpoint uncertainty and shared ligand identities limit any statistical claim.

A lower pKa increases the solution anion fraction, while a weaker Lewis base can bind Zn less strongly. As an exploratory stress test after inspecting the first results, we transferred the independently measured Zn-triamine slope of −0.40 kcal mol⁻¹ per pKa unit. This recovers some gap accuracy but still gets only 1/3 position-isomer pairs correct. The slope measured on the surrogate omits the protein H-bond network, and Vina already includes an empirical metal interaction. These data favour using pH transfer conditional on measured intrinsic/reference information and retaining conventional docking scores by default. [Zn surrogate experiment](https://doi.org/10.1021/ja00045a002) · [Conversion to the energy slope](https://doi.org/10.1002/asia.200600360)

## Independent measurement of proton uptake

For ethoxzolamide, Morkūnaitė et al. provide observed binding enthalpies in phosphate and Tris. At matched pH and 25 °C, the difference divided by the buffer deprotonation enthalpy difference (47.5 − 5.10 kJ/mol) estimates the average proton uptake, provided other buffer-specific effects cancel. We used every exact matched pH in Table 1 for this compound and protein, the free ligand pKa 8.0 and older, separately measured enzymatic pH inflections of 7.3 (human CA I) and 6.8 (human CA II). The model predicts the uptake sign in 5/6 conditions, with MAE **0.2175 protons** versus **0.6904** for zero uptake. For CA I at pH 7.9, the model predicts +0.242 but the buffers give −0.408. The study's own pH fits give different receptor pKas (8.4 for I and 7.1 for II). Re-expressing those same data using their fitted parameters would overstate independent validation. [ITC Table 1 and buffer enthalpies](https://doi.org/10.3109/14756366.2014.908291) · [Independent enzymatic titrations](https://doi.org/10.1016/S0021-9258(18)95800-X)

## Reproduction

```sh
python tools/ionization.py --ligand-pka 8.29 --receptor-pka 7.01 --reference-ph 7 --new-ph 5 --kelvin 310.15
python benchmarks/proton_linkage/evaluate.py
python benchmarks/proton_linkage/verify.py
python benchmarks/proton_linkage/dock.py /tmp/linkage build/vina --jobs 4
python benchmarks/proton_linkage/evaluate.py --runs /tmp/linkage/runs --out /tmp/linkage/results
python benchmarks/proton_linkage/verify.py --runs /tmp/linkage/runs
```

The first three commands use only the Python standard library. The last three require the branch's FastDock Vina binary with search-anchor support and the preparation dependencies used by the existing anchor benchmark, including RDKit (run version 2026.09.2), Meeko (0.8.0), Gemmi and NumPy. Source PDB SHA-256 digests are embedded in the driver. Rebuilding the inputs from the two original PDBs and rerunning all searches in a separate directory reproduced all 9 prepared PDBQTs and all 56 output PDBQTs byte for byte. The independently written checker evaluates the binding polynomial directly, differences its derivative, checks pH reversal and the alkaline bound, and recomputes the matched win counts from raw scores. The dock data do not contain per-ligand experimental pose references and contribute evidence on ordering only.
