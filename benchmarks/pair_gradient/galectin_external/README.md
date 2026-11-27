# Chemical transfer of the Gal3 ranking result

The exact pair derivative lost 10 pair-ranking credits against parent Pro across 403 reliable pairs in a different human Gal3 series. The equal-heavy-atom, equal-formal-charge subset gained six credits across 35 pairs. Receptor choice changed both the size and the direction of the effect. This series does not establish a transferable improvement in Gal3 affinity ranking.

## Independent experimental series

The [frozen protocol](PROTOCOL.md) and [47-compound manifest](manifest.csv) specify the chemical series, common input recipe, endpoints and budgets. Sörme et al. measured derivatives of LacNAc carrying aromatic amides in an intact human Gal3 competitive fluorescence polarization assay. Figure 4 specifies 8 µM inhibitor with individually marked exceptions. The experimental supplement, p. S6, specifies 1 µM galectin and 0.1 µM fluorescent probe with duplicate measurements. The FP temperature and individual measurement errors are unreported. The ITC measurements in the same paper use a separate protocol. Affinity-gap evidence is reported as the dimensionless quantity ΔΔG/RT = ln(Kd,weak / Kd,strong) for each pair. [Original article and Figure 4](https://doi.org/10.1021/ja043475p) · [Experimental supplement](https://doi.org/10.1021/ja043475p.s001)

The pinned public GFA source attributes 50 unique, exact Figure 4 structures to this DOI. Compound Nos. 4, 12 and 44 used other inhibitor concentrations and were excluded. Nos. 9, 15, 24, 41, 42 and 55 have no eligible GFA row attributed to this DOI. No. 55 is also censored. All 47 retained Kd values agree with an independent transcription of the primary figure. No retained isomeric molecular graph or Murcko scaffold occurs in the earlier 13 Peterson compounds. The new compounds have twelve Murcko scaffolds but share a closely related chemical series. GFA attributes a single formal charge state, neutral, to all 47 supplied input graphs. Four contain nominally neutral carboxylic acids. Their aqueous microstate populations were not measured in this benchmark. [GFA release](https://github.com/choutkaj/GFA-ligand-set/tree/04dd9b0024ba3a81730731d256b6e7b285977aec)

The two archived Gal3 receptor preparations and 24 Å boxes are identical across arms. Receptor 1KJL was reported by Sörme's group with LacNAc bound. The pocket box center is known to both algorithms. Query ligand input coordinates were regenerated from chemistry alone with the same eight-conformer ETKDG/MMFF protocol as the Peterson panel. Preparing each of the 47 ligands once took 105.50 s summed over processes, shared by both arms. This experiment tests chemical transfer on the same target and two receptor states. It does not provide an independently held-out receptor family.

## Measured ranking

Each arm uses the first pose's native score, seed 214039, one CPU, one trajectory, cap 12,000 and at most nine modes. The consensus averages the two first-pose scores for each ligand. The reliable subset requires at least threefold separation in published Kd. Exactly-threefold comparisons are calculated using decimal or integer experimental units. Tied predicted scores receive half credit. The simple size control ranks the ligand with more heavy atoms as stronger.

| Input receptor | All ≥3-fold pairs, parent → linear | Equal heavy atoms and charge, parent → linear | Spearman, parent → linear | Kendall τ-b, parent → linear |
|---|---:|---:|---:|---:|
| 1KJL | 242.5 → 197.5 / 403 | 21 → 18 / 35 | 0.121 → 0.011 | 0.083 → 0.020 |
| 3T1L | 174 → 204 / 403 | 7 → 19 / 35 | −0.134 → 0.044 | −0.087 → 0.035 |
| Two-state consensus | **198 → 188 / 403** | **15 → 21 / 35** | **−0.059 → −0.019** | **−0.020 → −0.011** |
| Size only | 259.5 / 403 | 17.5 / 35 | 0.191 | 0.134 |

The 35 size-controlled pairs involve 30 distinct compounds. Omitting any one compound leaves the linear-minus-parent change between +2 and +8 credits on this subset. On all 403 pairs, single-compound omissions span −25 to +16 credits. Compound pairs share ligands, and both receptor rows use exactly the same measured Kd values. These denominators do not supply independent binomial trials. The full [per-compound scores](scores.csv), [pair details](pairs.csv.xz), [rank statistics](summary.csv) and [leave-one-out results](influence.csv) preserve these dependencies.

A paired bootstrap resampled 47 ligand identities with replacement, preserving both receptor outputs and all incident pairs. Over 10,000 draws, the consensus accuracy change of −2.48 percentage points had a 95% percentile interval of −26.91 to +22.99 points. The size-controlled change of +17.14 points had an interval of −35.29 to +57.69 points (9,997 draws retained at least one qualifying pair). These intervals describe within-series sampling uncertainty. [Paired compound bootstrap](compound_bootstrap.csv)

Compounds 45 and 46 carry two LacNAc moieties. The intact-protein competition result and a single static carbohydrate-recognition site may represent multivalent ligands differently. As a labeled sensitivity, removing both leaves parent at 178/375, linear at 168/375 and size-only at 253.5/375. Docking both receptor conformations gives opposite effects even among the equal-size pairs. This comparison cannot resolve whether orientation, solvent, local chemistry or receptor geometry explains the difference.

## Measured cost and independent checks

| Jobs per arm | Joint energy/gradient evaluations | Local optimizations | Total elapsed |
|---|---:|---:|---:|
| Parent Pro, 94 | 1,168,128 | 18,234 | 15.371 s |
| Exact derivative, 94 | 1,171,993 | 18,280 | 15.414 s |
| Exact derivative at independently transported caps, 94 | 1,171,526 | 18,272 | 15.496 s |

Linear took 0.28% longer and used 0.33% more joint calls at the same cap. The 94 paired single-job wall ratios ranged from 0.744 to 1.343, with 71/94 within ±5%. The parent and linear arms alternated order across shuffled conditions on one pinned CPU. The preset timing control used caps of 11,993 (1KJL) and 12,000 (3T1L), derived from the earlier, independent chemical panel. The 3T1L jobs were reused. The transported caps left the consensus rankings and paired credits unchanged. The prespecified series-specific recalibration threshold was 5% aggregate cost difference and was not triggered.

An independent [raw-output audit](audit.py) reconstructed scores from integer millikcal/mol ticks, checked every atom serial/type and finite coordinate in 1,346 returned modes from 235 unique jobs, reparsed 1,410 counter fields, and reproduced all 36 affinity-summary rows. All 47 experimental labels were re-read from the primary-figure transcription in integer 0.01 µM units. Re-evaluating from the archived outputs under different Python hash seeds reproduced the result tables and leave-one-out rows exactly. [Audit](audit.json) · [Run and binary identity](run.json) · [Per-condition cost and scores](fixed.csv)

## Reproduction

Use Python 3.12 with the dependencies described in the [pose study](../pose/README.md). From the repository root, the first block verifies the frozen measurements without recomputing docking.

```bash
P=benchmarks/pair_gradient/galectin_external
W=/tmp/fastdock-gal3-external
mkdir -p "$W"
(cd "$P" && sha256sum -c SHA256SUMS)
tar -xJf "$P/raw.tar.xz" -C "$W"
cp "$P"/{manifest,fixed,transported}.csv "$W/"
python "$P/evaluate.py" "$W"
python "$P/audit.py" "$W"
```

To recompute from the release embedded in `benchmarks/galectin_affinity/evidence.tar.xz`, use a new empty work directory and set `B` to the pinned arms built by [build_arms.py](../build_arms.py). Supply a valid affinity index for your machine. The original measurement used CPU 8.

```bash
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
python "$P/panel.py" "$W" --prepare-jobs 5
python "$P/run.py" "$W" --parent "$B/parent/build/vina" --linear "$B/linear/build/vina" --cpu 8
python "$P/evaluate.py" "$W"
python "$P/audit.py" "$W"
```

Original numerical verification used the published Figure 4 page image (SHA-256 `439373874defc9a0937f94eaad356386dcc88fb46e9a715264427e63c237ff74`), the [author supplement](https://ndownloader.figshare.com/files/5138623) (`f81e641a75d2714c4daa3269a716656b687642a7d655a2c7c93f92282bfc959b`) and pinned GFA.csv (`2eb8a77960dac49e3f52b66bd0c9897ea8eae8d595be7e2cdb4cf20034953a71`). All checked-in results have [SHA-256 digests](SHA256SUMS).
