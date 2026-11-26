# Torsion search on four additional human therapeutic targets

Eight X-ray complexes provide two reciprocal redockings and two cross-dockings per target. The eight ligands have seven distinct Bemis–Murcko scaffolds. The two CHEK1 ligands share a scaffold. The experiment measures known-site recovery with a rigid, mostly dry receptor. The original input and every native output are archived here.

| Target | X-ray structures and ligand | Reciprocal pocket Cα count and RMSD |
| --- | --- | --- |
| BACE1 | [2B8L / 5HA](https://www.rcsb.org/structure/2B8L), [2ZJM / F1M](https://www.rcsb.org/structure/2ZJM) | 47, 0.717 Å |
| EPHX2 | [3ANS / S38](https://www.rcsb.org/structure/3ANS), [1VJ5 / CIU](https://www.rcsb.org/structure/1VJ5) | 58, 0.495 Å |
| CHEK1 | [2E9P / 77A](https://www.rcsb.org/structure/2E9P), [2YWP / A42](https://www.rcsb.org/structure/2YWP) | 54, 0.584 Å |
| BRD9 | [6Y7L / QMG](https://www.rcsb.org/structure/6Y7L), [6YQW / 82I](https://www.rcsb.org/structure/6YQW) | 37, 0.392 Å |

Both BACE1 structures carry the same K136A/E138A mutations in canonical numbering. Neither deposited query has a covalent LINK to protein. CHEK1 uses wild-type kinase domains. Before docking we removed the proposed 7BJJ receiver because its ten mutations reproduce properties of the LRRK2 pocket. The proposed EPHX2 backup 3WKE was unnecessary because 1VJ5 passed the preparation checks. These substitutions were decided before viewing docking outcomes.

Query geometries were generated from the wwPDB CCD graph by ETKDGv3 followed by MMFF, using the same stereochemistry, atom types, charges and torsion tree for all search arms. The input contains no experimental ligand coordinates. Each receiver defines its own 24 Å box using its deposited ligand. Identity-matched pocket Cα atoms within 14 Å of the corresponding ligand centroids define the receptor frame for cross-docking. RMSD minimizes over exact chirality-preserving graph symmetries and uses no ligand superposition. Preparation uses chain A and keeps deposited Zn/Mg. Waters are omitted; other omitted non-water heavy atoms must be at least 8 Å away. Protonation follows the archived preparation. The crystal pose transferred to a different receptor is an aligned experimental reference and can clash or change its preferred mode in that receiver.

All arms use four trajectories, one CPU, Late refinement of up to 40 candidates, up to ten output modes and 100 kcal/mol energy range. Compare arms within the same protein, query, receiver, seed and evaluation cap. The threshold is 2 Å. The oracle is the best RMSD in the returned set after native filtering.

## Paired results

| Per-trajectory cap | Seed cohort | Comparison | Top-1 native → candidate | Oracle native → candidate | Paired top-1 gains / losses |
| ---: | --- | --- | --- | --- | --- |
| 3,000 | original | balanced torsion | 5 → 5 / 48 | 6 → 5 | 4 / 4 |
| 12,000 | original | balanced torsion | 8 → 5 / 48 | 9 → 7 | 2 / 5 |
| 3,000 | original | prepared starts | 5 → 9 / 48 | 6 → 9 | 5 / 1 |
| 12,000 | original | prepared starts | 8 → 13 / 48 | 9 → 14 | 5 / 0 |
| 24,000 | fresh | prepared starts | 8 → 8 / 48 | 11 → 9 | 2 / 2 |
| 12,000 | fresh, bridge | prepared starts | 8 → 7 / 48 | 11 → 8 | 2 / 3 |
| 24,000 | original, bridge | prepared starts | 9 → 13 / 48 | 10 → 14 | 4 / 0 |
| 250,000 | two additional seeds, cross only | prepared starts | 1 → 3 / 16 | 5 → 6 | 2 / 0 |

The prepared-start hypothesis was introduced after inspecting balanced-torsion results. The first fresh-seed cohort and its 24,000 cap were declared before measuring prepared starts. The two subsequent bridge comparisons keep each seed cohort fixed while changing the cap. Within this panel, the difference between cohorts persists at both caps. Shared seeds, repeated receptors and repeated molecules limit the number of independent examples. Outcome gains on one seed cohort cannot establish general transfer. For prepared starts, odd trajectories preserve all recorded native search counters. Measured search and explicit-refinement costs are in `counter_audit.csv`.

The higher-budget cross-docking control used fresh seeds 62347 and 62383. Median top-1 RMSD fell from 7.077 to 2.886 Å. With the full PoseBusters docking validity criteria, 0/16 native and 1/16 prepared poses combined validity with RMSD ≤2 Å. The additional physically valid hit was BRD9 82I in 6Y7L at 1.025 Å. Three otherwise near-native CHEK1 77A predictions triggered the intramolecular geometry check. Fifteen native and thirteen prepared top-1 poses were physically valid at any RMSD. Mean search evaluations were 1,000,112 and 1,000,119 respectively. Median per-invocation elapsed time was 4.27 and 4.47 seconds with six concurrent one-CPU invocations; machine load and target composition influence these timings. Eight receptor–query directions and two seeds limit precision.

## Mechanism diagnostics

The local experiments transfer the crystal torsions into the graph-generated query geometry, align the query for diagnostic initialization, then relax it with the native potential. They separately relax the crystallographic geometry starting in the same receptor frame. These use the answer to probe the mechanism and are excluded from prospective docking performance.

For CHEK1 A42 the cognate crystal geometry relaxes to 0.250 Å, while transferring its active torsions into the ordinary input yields 2.113 Å. In 2E9P these values are 0.376 and 2.467 Å. Geometry held fixed by the torsion tree can change which basin is accessible even with measured input torsions. For BACE1 5HA in its cognate receptor the transplanted geometry relaxes to 0.489 Å, with common-reference INTER + INTRA 2.16 kcal/mol lower than the best top pose among three native 12,000-cap seeds. For 5HA in the paired receiver the transplant still relaxes to 0.676 Å, but its score is 9.46 kcal/mol higher than the best sampled pose. The latter gap is conditional on this particular pose, preparation and potential. It does not bound all other near-native conformations. The local results support different explanations for the cognate and transferred cases. Tables retain all controls.

The independent audit reads deposited ligand coordinates and reconstructs their bond orders and stereochemistry using a separate PDB parser and the prepared SMILES, then independently refits the receiver frame with SciPy. Pocket correspondences are assigned by PDBe SIFTS. The two BRD9 deposits use different UniProt isoform numbering; all shared unique seven-residue sequence stretches imply the same conversion. Every mode is checked against the primary metric. PoseBusters 0.6.3 docking checks are evaluated separately on each top-1 pose using the actual prepared receiver atoms, without the co-crystallized ligand. See `results/joint.csv` for the combined threshold and physical-validity endpoint.

## Retrospective geometry filter on frozen pools

All 224 native and prepared-start pools from both 12,000-cap cohorts and the cross-only 250,000-cap cohort were rescored with the PoseBusters 0.6.3 intraligand distance-geometry check. Keep the highest native-ranked passing pose; count an empty passing set as an abstention. The cutoff is 70% of the RDKit geometry lower bound (default docking configuration), ignoring hydrogens. Recalculate native-frame RMSD without pose alignment. Per-mode details include the closest pair, bond-path length, separation, lower bound and verdict. Every pool had a passing pose. This exploratory choice followed inspection of the initial results.

| Cap | Cohort | Native search, original → filtered top-1 hits | Prepared start, original → filtered top-1 hits |
| ---: | --- | --- | --- |
| 12,000 | Original (48 per arm) | 8 → 8 | 13 → 12 |
| 12,000 | Fresh (48 per arm) | 8 → 8 | 7 → 7 |
| 250,000 | Deep, cross only (16 per arm) | 1 → 1 | 3 → 2 |

All eight graph-generated inputs and seven of eight deposited query geometries passed the isolated clash test. The CHEK1 A42 crystal itself has an intraligand N···N distance of 2.830 Å across four covalent bonds, below the automated cutoff of 2.904 Å. For the chemically related 77A crystal, the distance between a urea NH and a pyrimidine N is 2.959 Å and passes. Transplanting its measured torsions into the ordinary prepared geometry gives 2.735 Å initially and 2.853 Å after local relaxation; directly relaxing the crystal gives 3.003 Å. The 77A top poses flagged at the largest cap have 2.801–2.863 Å N···N distances. These are possible intramolecular donor–acceptor contacts. A binary distance-geometry verdict and this matched panel support reporting the contact as a diagnostic. They provide insufficient evidence for rejecting all such poses.

## Reproduce

Use the Python 3.12 preparation environment from [the preceding cross-docking benchmark](../paired_crossdocking/README.md), with RDKit 2026.09.2, Meeko 0.8.0, Gemmi 0.7.5, PDBFixer 1.12.0 and OpenMM 8.6.1. For PoseBusters use the [Python 3.11 physical-validation environment](../openff_search_replay/README.md). Random seeds 51913, 51941, 51973 define the original cohort; 52121, 52147, 52163 define the fresh cohort. From the repository root:

```sh
B=benchmarks/torsion_proposals
W=/tmp/fastdock-torsions
mkdir -p "$W"
(cd "$B" && sha256sum -c fixtures.sha256)
for a in "$B"/fixtures/*.tar.xz; do tar -xJf "$a" -C "$W"; done
cp "$B"/results/*.csv "$W"/
python "$B/audit.py" . "$W"
python "$B/summary.py" "$W"
python "$B/counters.py" "$W"
python "$B/finish.py" "$W"

# Optional fresh structure preparation and prospective reruns.
python "$B/prepare.py" . "$W"
for cap in 3000 12000; do
    for arm in native balanced prepared; do
        python "$B/paired.py" . "$W" build/vina "$arm" "$cap"
    done
done
for arm in native prepared; do
    python "$B/paired.py" . "$W" build/vina "$arm" 24000 --seeds 52121 52147 52163
    python "$B/paired.py" . "$W" build/vina "$arm" 12000 --seeds 52121 52147 52163 --label=-bridge
    python "$B/paired.py" . "$W" build/vina "$arm" 24000 --label=-bridge
    python "$B/paired.py" . "$W" build/vina "$arm" 250000 --seeds 62347 62383 --cross-only
done
python "$B/diagnostic.py" . "$W" build/vina
# Run in the PoseBusters environment.
python "$B/check_validity.py" "$W"
python "$B/check_validity.py" "$W" --bridge
python "$B/check_validity.py" "$W" --deep
python "$B/clash_control.py" "$W"
```
