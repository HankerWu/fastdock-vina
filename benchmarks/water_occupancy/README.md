# Recurrent and displaceable water on matched experimental complexes

A water oxygen can recur across many structures because the site has stable polar partners. Assigning a cost to vacating every such site predicts worse BRD4 poses. This experiment distinguishes an independently witnessed displacement from a guessed replacement hydrogen bond. It uses the same candidate coordinates and native scores in each comparison. All methods are opt-in rescorers.

## Mechanisms and controls

For each receptor, each external structure contributes at most one water per site. Protein-only sequence and local-backbone alignment, water occupancy of at least 0.5 and the 1.0 Å clustering convention come from `tools/extract_independent_receptor_evidence.py`. We require at least two source structures and one receptor donor/acceptor 2.4–3.5 Å from the site's center. We exclude sites closer than 2.3 Å to any receptor heavy atom or 3 Å to a metal. The BRD4 target receptor is chain B for 4OGJ and chain A elsewhere. All external sources use chain A. Source structures are disjoint from the evaluated ligand crystal entries. They can share protein constructs or ligand fragments. The masks are receptor-frame specific.

With a site observed in k of D external structures, the structural recurrence proxy is p=(k+0.5)/(D+1). Let A rise smoothly from zero to one between 2.0 and 2.6 Å from the closest ligand heavy atom. Let u be the single strongest geometrically compatible ligand contact, tapered over 2.4–3.5 Å with a receptor–water–ligand angle of 65–160°. The depth is one Vina hydrogen bond (−0.587439 kcal/mol). Each site contributes −RT log[(1−p)+p A exp(−u/RT)] relative to the unbound proxy, with T=298.15 K. The implementation sums contributions over sites. Exclusion-only sets u to zero. Contact-only sets A to one. We compare every variant to the uncorrected score in the same search pool.

A fixed follow-up hypothesis allows the same occluding ligand atom to replace a fraction f of the water's original receptor polar contacts. The added partition weight is p(1−A)f. We compute f from the radial hydrogen-bond window and the angle to the original water direction at the receptor atom (full credit through 20°, zero at 45°), then take the best single occluding atom. Separate controls accept either any atom or only a complementary donor/acceptor. The model treats explicit polar hydrogens and OA/NA PDBQT types as evidence of the atom's role. It omits explicit water orientation, coupled water–water motion, protein relaxation, solution-state populations and the chemical potential of displaced solvent. Structural recurrence and the chosen temperature are modeling inputs. They cannot establish a site's binding free energy or calibrate predictions in a different-temperature assay.

After the first comparison we tested two site masks. The diagnostic mask uses overlap within 2.3 Å of the withheld experimental ligand, available **only after scoring**. The predictive mask uses the same threshold on ligand heavy atoms from four *other* complexes chosen from a published structural study of 101 BRD4–ligand complexes. We froze these four sources before scoring their geometrical overlaps, after inspecting the initial failure. Divakaran et al. directly observed water-network reorganization in 6MH1. The 101-structure study identifies 6MH1 and 5I88 as displacers, and reports absence of a site-5 water in 4GPJ and 5DLZ. Our rule requires actual heavy-atom overlap in aligned coordinates. Absence of a water in a PDB file is insufficient. Test and witness compounds differ, but some share motifs, including a substituted imidazole. The study's published collection may contain members of the evaluation set. This follow-up is data-informed and requires fresh ligand families for prospective validation. [BRD4 structural survey](https://doi.org/10.1039/D1OB00658D) · [6MH1 experimental structure and thermodynamics](https://pubmed.ncbi.nlm.nih.gov/30253095/) · [Forli and Olson's explicitly displaceable water model](https://doi.org/10.1021/jm2005145)

## Experimental results

BRD4 uses same-protocol ITC Kd at 15 °C, pH 7.4 (37, 123, 164, 3,400 nM). CA II uses same-assay stopped-flow Ki at pH 7.5 (368.5, 57.0, 30.4, 2.4, 72.5 nM). CA I has a single sixfold pair, MZM 1.2 and AZM 7.2 µM by DSF at 57 °C and pH 6.1. Independent 25 °C ITC agrees on the sign of this CA I pair. No Ki/Kd values were pooled across temperatures, assay types or targets. In each receptor and seed the lower-scored molecule should be stronger. Threefold or greater experimental gaps are eligible, with half credit for predicted ties. Sequences of search seeds and receptor structures reuse the same measured ligand pairs. [BRD4 ITC](https://doi.org/10.1038/nchembio.1471) · [CA II stopped-flow](https://doi.org/10.1021/acs.jmedchem.0c00733) · [CA I DSF](https://doi.org/10.1038/s41598-018-37072-x) · [CA I preparation and independent measurements](../anchor_search/holdout.md)

Values below are first-pose successes / pair credits. The CA II dual-state criterion requires whole-ligand RMSD ≤2 Å and the bound ammonium orientation. For the other rows only whole-ligand RMSD ≤2 Å is required. “Native pool” describes how the coordinates were generated. “Dry” specifies the score used to select among them.

| Fixed candidate panel | Case / pair counts | Dry | Exclusion only | Contact only | Joint | Any replacement | Polar replacement |
|---|---:|---:|---:|---:|---:|---:|---:|
| BRD4 native pool | 24 / 30 | 1 / 10 | 1 / 14 | 1 / 12 | 1 / 14 | 1 / 10 | 1 / 10 |
| BRD4 transfer pool | 24 / 30 | 7 / 9 | 2 / 12 | 5 / 9 | 2 / 13 | 2 / 12 | 2 / 13 |
| CA II native pool | 20 / 28 | 0 / 14 | 0 / 13 | 0 / 13 | 0 / 14 | 0 / 18 | 0 / 18 |
| CA II scaffold pool | 20 / 28 | 5 / 8 | 6 / 17 | 5 / 9 | 5 / 17 | 5 / 18 | 5 / 18 |
| CA II dual-state scaffold pool | 20 / 28 | 1 / 9 | 4 / 16 | 1 / 11 | 3 / 17 | 3 / 18 | 3 / 18 |

Whole-pose ≤2 Å successes on the CA II dual-state pool were 4/20 with both dry and joint scores. Requiring the experimental ammonium orientation exposes the 1→3 difference. Eight pools contained a candidate meeting both criteria. For the same seven distinct CA II labels, median-gap pair credits rose from 3/7 dry to 5/7 joint. All ten corrected repeat comparisons involve the weakest compound 34, and both losses involve the strongest 46. Omitting any one water source gives 3–4/20 signed hits and 16–19/28 repeat credits. Improved counts do not establish broad ranking across the remaining chemistry.

The BRD4 transfer pool contains an acceptable candidate in 12/24 cases. Dry → joint changes include one crossing to success and six to failure. The aggregate +4 Kd credits comprise six gains involving the weakest compound SB-202190 and two losses involving the strongest BI-2536. The five distinct pairs yield 1/5 and 2/5 when ranked by the median gap across seeds and receptors. Leaving out one water donor gives 2–4 pose successes and 8–14 Kd credits.

Of 54 observations of unanimous BRD4 sites within 4 Å of the experimental ligand across 12 ligand–receptor combinations, 42 had a modeled test water within 1 Å and nine were occupied by an experimental ligand heavy atom within 2.3 Å. These are site–ligand–receptor observations, including repeated projections of each crystal. CA II gave 14 observed waters and ten ligand overlaps among 30 observations. In the 3MXF frame both TG compounds displace site 32 (TG-101348 OAM at 1.96 Å and TG-101209 C15 at 1.10 Å). SB-202190 O20 displaces site 15 at 0.36 Å. External 5DLZ witnesses site 32 and external 6MH1 witnesses site 15. Chemical-contact replacement recovered none of the BRD4 pose losses.

| BRD4 transfer selection | Top-1 ≤2 Å / 24 | Kd credits / 30 |
|---|---:|---:|
| Dry | 7 | 9 |
| Joint | 2 | 13 |
| External ligand displacement mask | 6 | 7 |
| Withheld crystal displacement mask (diagnostic) | 6 | 12 |

Restricting transfer to the 20 cross-docking cells gives 6, 1, 4 and 4 top-pose successes for dry, joint, external and diagnostic, respectively. The four cognate receptor cells give 1, 1, 2 and 2. Prepared receptor conformations can still contain induced-fit information even though their bound ligand and water coordinates were discarded before prediction. Relative to the joint model, both masks recovered the same four cells without losing another hit. Kd lost six repeat credits with the external mask and one with the diagnostic mask. Within the native BRD4 pool the corresponding four pose counts stayed at 1/24; Kd credits were 10, 14, 8 and 8. The gain from excluding known displaceable waters concerns pose selection and falls short of dry on the transfer pool. The binary external mask removes sites for every query, including cases where a water can stay bound. It supplies neither a ligand-specific displacement free energy nor a common solution for pose and affinity ranking.

| Receptor | Query | Seed | Dry RMSD (Å) | Joint | External mask | Diagnostic mask |
|---|---|---:|---:|---:|---:|---:|
| 3MXF | 4OGJ | 314159 | 2.855 | 3.054 | 1.552 | 1.552 |
| 4OGI | 4OGJ | 314159 | 1.705 | 6.607 | 1.705 | 1.705 |
| 4OGJ | 4O76 | 271828 | 0.938 | 2.214 | 0.938 | 0.938 |
| 4OGJ | 4OGJ | 314159 | 1.841 | 2.790 | 1.761 | 1.831 |

### Same-isoform holdout

We froze the same weights and radii for CA I and CA XII, then regenerated site maps from four entirely different external sources per isoform. Testing used 519 previously generated poses from native, Zn-directional point and radial-shell searches. Across both targets, first-pose counts (native / directional / shell) in the 36 Å box stayed at 4 / 9 / 5 of ten cases under the joint model. In the 20 Å box the counts were 6 / 5 / 5 dry and 6 / 6 / 5 joint. The net point gain hides a CA I loss (5→4 of eight) and two CA XII gains (0→2 of two). Across the three search arms, the CA I experimental pair gave 1 / 2 / 2 correct of four repeats in the wide box and 2 / 3 / 2 in the narrow box, with identical credit for dry, joint and both substitution variants. CA XII supplies one ligand and no rank evidence. Wide-box native pose selection fell from 4/10 to 3/10 for either substitution variant. These holdouts falsify a general same-assay rank improvement and give mixed pose results. No weight was refit on either holdout.

## Reproduction

Use the frozen PDBQT coordinate pools, native scores and experimental preparation from the linked repository benchmarks. NumPy, SciPy, Gemmi and Biopython are needed. Source files have SHA-256 digests. All shell paths below refer to the evidence bundle layout or the original working directory names.

```bash
B=benchmarks/water_occupancy
W=/tmp/fdv-water3
BRD=/tmp/fdv-brd4
CA=/tmp/fdv-second
BP=/tmp/fdv-search-brd4
CP=/tmp/fdv-anchor-scaffold
SP=/tmp/fdv-ammonium-scaffold
LG=/tmp/fdv-anchor-holdout
SM=/tmp/fdv-anchor-holdout20
python "$B/fetch_donors.py" "$W/donors"
python "$B/evaluate.py" --donors "$W/donors" --brd "$BRD" --ca "$CA" --brd-pools "$BP" --ca-pools "$CP" --state-pools "$SP" --out "$W/result"
python "$B/holdout.py" --donors "$W/donors" --large "$LG" --small "$SM" --out "$W/holdout"
python "$B/validate_sites.py" --brd "$BRD" --ca "$CA" --out "$W/result"
python "$B/oracle_ablation.py" --brd-pools "$BP" --ca-pools "$CP" --state-pools "$SP" --out "$W/result"
python "$B/displacement_witness.py" --donors "$W/donors" --brd "$BRD" --brd-pools "$BP" --out "$W/result"
python "$B/check.py" --out "$W/result" --holdout "$W/holdout" --brd-pools "$BP" --ca-pools "$CP" --state-pools "$SP" --large "$LG" --small "$SM"
python "$B/check_masks.py" --out "$W/result" --brd-pools "$BP"
```

The standalone research tool accepts any external source set passing the alignment criteria. `python tools/water_occupancy.py prepare --help` describes donor and optional displacement inputs. `python tools/water_occupancy.py score SITES.json POOL.pdbqt` emits energy corrections in kcal/mol for each stored MODEL. The four source IDs per target are listed in each preparation script. A caller must exclude its query crystal from the donor and displacement sets. BRD4 witnesses are external known binders, so the mask requires available structural information.

`results/poses.csv` and `holdout_results/candidates.csv` contain every scored candidate. Site JSON files retain aligned positions, individual donor IDs and residues, receptor partners, and fit diagnostics. `results/mask_comparison_cells.csv` and `results/mask_comparison_pairs.csv` contain every matched BRD4 mask transition. The separate PDBQT parser and scalar partition in `check.py` reproduce all five scores on 1,656 saved coordinates to 3.55×10⁻¹⁵ kcal/mol. `check_masks.py` sums retained sites directly, checks the chosen model and score in each of 48 cells, and recomputes all 60 repeated-pair outcomes for the four BRD4 methods. These checks validate arithmetic and bookkeeping; prospective solvent thermodynamics would need independent experimental evidence.
