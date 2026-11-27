# Why a unit-weight MMFF replacement transfers poorly between BACE1 and Factor Xa

On frozen native first poses, removing the **intraligand electrostatic part** of the MMFF94s bound–free difference leaves BACE1 at 48/72 pair credits, down from 54/72 with full MMFF and above 23/72 with native Vina. The same removal takes Factor Xa from 206/295 to 162/295, further below native's 223/295. A modification that deletes this term from the score has no basis here for deployment. The results identify a composition effect and a receptor dependence in the earlier full replacement.

## What was held fixed

We reconstructed 168 native first poses (42 ligands, two receptors and two seeds per ligand). BACE1 supplies two distinct Ki series, CHEMBL2015090 and CHEMBL2155063. FXa supplies the single series CHEMBL661804. Every heavy coordinate, formal charge and chemical graph comes from the archived native run. Each free reference is the lowest full-MMFF energy among eight ETKDGv3 conformers generated with seed 777. Each bound pose has only its hydrogens relaxed. Recalculated full scores, archived native scores and all free references match the original numeric ledgers exactly in the pinned runtime. We independently evaluated the electrostatic term as both a difference of two single-point energies and an explicit electrostatic-only single-point energy, requiring agreement within 1e-6 kcal/mol. [Pre-analysis protocol](PROTOCOL.md) · [BACE assay and receptor preparation](../../geometry_comparison/affinity/README.md) · [FXa paper, assay and microstates](../../affinity_partition/README.md)

Let B and F denote the selected bound and reference conformers. Write E(B) − E(F) = ΔE_elec + ΔE_other. The term ΔE_elec comes from the MMFF94s intramolecular charge model. ΔE_other contains bonded and van der Waals terms. The archived full replacement is S_inter + f(ΔE_elec + ΔE_other), with S_inter = S_native − f(I − U) and the ligand's original Vina flexibility multiplier f. The two ablations use S_inter + fΔE_other and S_inter + fΔE_elec. All four predictions use identical poses. Hydrogen coordinates and the free-conformer choice were held at their full-potential values across the ablations. Removing electrostatics before minimization would ask a different question. The protocol changes no search evaluations. The diagnostic needed 336 free conformer generation attempts and 168 bound hydrogen relaxations beyond native docking, plus component single points. On this Linux host a complete offline two-worker replay took 27.35 s, including archive reads. This timing is specific to the diagnostic; neither a speedup nor a same-cost docking result follows.

## Results by receptor and assay

Values are correct direction credits for exact same-assay Ki values separated by at least threefold. Ties receive half credit. The BACE rows have two-seed ligand medians. The historical FXa estimator uses the median of matched score differences across the stated receptor/seed blocks. There are only 36 distinct BACE pairs when duplicate receptors are collapsed. The distinct FXa set contains 295 pairs among 28 compounds. BACE and FXa pair counts therefore have different aggregation conventions.

| Assay | Receptor | Qualifying pairs | Native | Full MMFF | Without electrostatic term | Electrostatic term only |
|---|---|---:|---:|---:|---:|---:|
| BACE CHEMBL2015090 | 4DJW | 14 | 1 | 10 | 11 | 7 |
| BACE CHEMBL2015090 | 4GID | 14 | 8 | 7 | 4 | 10 |
| BACE CHEMBL2155063 | 4DJW | 22 | 4 | 15 | 15 | 8 |
| BACE CHEMBL2155063 | 4GID | 22 | 10 | 22 | 18 | 9 |
| FXa CHEMBL661804 | 1FJS | 295 | 231.5 | 217 | 161 | 204 |
| FXa CHEMBL661804 | 2P16 | 295 | 187 | 208 | 148 | 195 |
| FXa CHEMBL661804 | two receptors | 295 | 223 | 206 | 162 | 204 |

Deleting one FXa ligand and all its incident comparisons at a time leaves the full-minus-native pooled change between −12.09 and −1.10 percentage points. Removing electrostatics from full MMFF costs another 8.43–17.58 points under the same deletion. These sensitivity ranges describe dependence on single compounds. For the receptor-specific full replacement, Spearman/Kendall on 1FJS change from .631/.466 to .507/.383, while on 2P16 they change from .321/.181 to .485/.351. All six assay/receptor correlation sets, exact unrounded pair gaps, ≥10-fold results, within-scaffold/cross-scaffold subsets and exact equal-size/charge controls are available in [correlations.csv](correlations.csv), [pairs.csv.xz](pairs.csv.xz), [summary.csv](summary.csv) and [influence.csv](influence.csv). Scaffold assignment uses the RDKit Bemis–Murcko graph on the prepared ligand. These retrospective subsets do not constitute a scaffold holdout.

## Which molecules account for the discrepancy?

In the 90 FXa pairs comparing different formal charges, native/full/without-electrostatics give 69/72/20 credits. In the 205 same-charge pairs they give 154/134/142. The 198 pairs in which both ligands are +2 alone give 149/127/137. Thus the full replacement's net 17-credit loss consists of **20 credits lost within a charge class and 3 gained across classes**. In the 17 pairs of equal heavy-atom count and formal charge the same three arms give 9/7/11. Seven within-class +1 comparisons and 17 exact matched comparisons allow limited interpretation. No ligand was excluded by its prediction or score.

| Exact equal heavy-atom count and formal charge | Pairs | Native | Full MMFF | Without electrostatic term | Electrostatic term only |
|---|---:|---:|---:|---:|---:|
| BACE CHEMBL2015090, 4DJW | 4 | 0 | 2 | 3 | 1 |
| BACE CHEMBL2015090, 4GID | 4 | 3 | 2 | 1 | 2 |
| BACE CHEMBL2155063, either structure | 0 | — | — | — | — |
| FXa, 1FJS | 17 | 12 | 10 | 12 | 8 |
| FXa, 2P16 | 17 | 6 | 7 | 9 | 9 |
| FXa, two receptors | 17 | 9 | 7 | 11 | 8 |

Five FXa compounds were prepared with two cationic amidines and one carboxylate (net +1), and 23 with net +2. For each compound, take the median over its four bound poses, then the median across compounds in its class. The normalized electrostatic contribution fΔE_elec has medians **+76.29 and +10.03 model-score units** in the +1 and +2 groups, while median measured Ki is 370 and 2.8 nM respectively. The corresponding electrostatic medians in the two uniformly +1 BACE assays are +1.76 and +1.12. Deleting the FXa term removes a ligand-class signal that happens to correlate with measured potency. It does not repair within-class ranking.

A post-hoc coordinate check measured the shortest central-carbon distance from a carboxylate to either amidinium, avoiding resonance-dependent atom labels. Across those five zwitterions the individual full-MMFF free minima range from 3.75 to 4.35 Å, while the per-ligand medians of the four bound distances range from 6.45 to 10.95 Å. Every free minimum is more compact by this measure than its ligand's median bound pose. Compact opposite charges and a large differential vacuum electrostatic term can be directly observed in the inputs and energy calculation. The FXa records supply no assay pH, temperature or replicate errors. The charged microstates are preparation assignments, and neither the aqueous conformer populations nor the true bound poses of this series are known. A difference between group medians cannot identify the effect of formal charge on Ki separately from chemical substitutions. [Ligand-level energies and distances](scores.csv)

In BACE1, non-electrostatic terms alone retain 25 of the 31 net pair credits gained over native, and the two receptor structures give opposite signs in CHEMBL2015090. FXa can improve on 2P16 while deteriorating on 1FJS under the same scoring formula. These findings reject both a universal unit-weight replacement and the proposed explanation that simply deleting ligand Coulomb energy repairs transfer. They cannot distinguish each bonded or van der Waals contribution or show that solution-phase strain is irrelevant. An informative subsequent test would select microstate- and scaffold-matched series independently of this analysis and supply an independently validated solution reference. [Halgren on polar conformational energies](https://doi.org/10.1002/%28SICI%291096-987X%28199604%2917%3A5/6%3C587%3A%3AAID-JCC4%3E3.0.CO%3B2-Q) · [Gokcan and Isayev, solution-phase ligand strain](https://doi.org/10.1021/acs.jcim.6c02468)

## Reproduction

Use the pinned chemistry stack (Python 3.12, RDKit 2026.09.2, Meeko 0.8.0, gemmi 0.7.5, NumPy 2.3.5, SciPy 1.17.0). The original archived PDBQT files, assay labels, scales and baseline ledgers are already in the repository. From repository root:

```sh
D=benchmarks/frozen_core/strain_components
W=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 python "$D/decompose.py" "$W" --jobs 2
python "$D/summarize.py" "$W"
for table in scores.csv pairs.csv.xz correlations.csv influence.csv summary.csv; do
  cmp "$W/$table" "$D/$table"
done
```

The second command also runs independently on the committed score table, using only Python and SciPy. Its checks reconstruct each prior native/full pair count directly from the Ki values and the per-seed scores.
