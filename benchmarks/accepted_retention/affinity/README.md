# Same-assay affinity under one-slot candidate replacement

Changing a single retained search pose affected the first pose of just one of 104 receptor–compound combinations. On human galectin-3 (Gal3), every change in the two-receptor consensus involves compound 32 on receptor 1KJL. Human carbonic anhydrase II (hCA II) showed no change in any compound's first score, and every emitted mode across all four arms missed its experimental Zn–N contact. These two series do not establish a transferable affinity gain.

## Experimental inputs and endpoints

The protocol was frozen at commit `3a14147`, before inspecting test outcomes. Each compound used the same receptor, ligand microstate, pocket, new seed 614039, one trajectory, one thread, nominal `max_evals=12000`, nine output modes and a 100 kcal/mol output window across arms. The patch keeps a single geometrically distinct skipped hunt conformation and can replace one native candidate after pool filtering. Both panels were previously used in the repository for different mechanisms. The mechanisms and seed of this comparison were new. These are two experimental series and one realization of coupled search randomness per cell. [Protocol](PROTOCOL.md) · [Candidate definition and base commit](../core/PROTOCOL.md)

- **Gal3.** All 47 unique graphs from the frozen 8 µM competitive fluorescence polarization series, with the exact prior PDBQT inputs for 1KJL and 3T1L. Sörme et al. used intact human galectin-3, 1 µM protein, 0.1 µM probe, and duplicate measurements. The figure's Kd values are estimates from the inhibition measurement at the stated inhibitor concentration. Individual standard errors and an FP temperature are unavailable. Earlier filtering excludes separately marked inhibitor concentrations. Compounds 45 and 46 are bivalent, and drawn neutral acid forms in this source carry microstate uncertainty. The compounds span 12 Murcko scaffolds, with 31 in one scaffold. [Original article, Fig. 4 and supplement p. S6](https://doi.org/10.1021/ja043475p) · [Prepared fixture, curation and exclusions](../../pair_gradient/galectin_external/README.md)
- **hCA II.** All five co-crystallized ligands from Bonardi et al. The same human stopped-flow CO₂ hydration assay at pH 7.5 reports Ki for compounds 34, 41, 42, 46 and 48 of 368.5, 57, 30.4, 2.4 and 72.5 nM. Means of at least three determinations have estimated 5–10% errors. The 6WQ4 and 6WQ8 receivers retain Zn. Every ligand uses the existing assigned microstate, including deprotonated coordinating sulfonamide. The independently generated 36 Å boxes use the Zn and coordinating histidines. This endpoint is separate from the Gal3 Kd series. No fallback panel was needed. [Original article, Tables 1–2 and methods](https://doi.org/10.1021/acs.jmedchem.0c00733) · [Source hashes, alignment and ligand preparation](../../cutoff_switch/hcaii/README.md)

Lower first-pose native Vina scores predict lower Ki or Kd. The primary consensus is the arithmetic mean of the two receptor scores for each ligand. Scores use integer milli-kcal/mol ticks for ranking, and the exact sum preserves the ordering of the mean. Decimal Ki/Kd ratios ≥3 (pK difference ≥log10(3)) determine eligible pairs. A tied predicted score receives half credit. Spearman ρ and Kendall τ-b use all compounds and native mean/marginal rankings. The simple size baseline predicts that the larger ligand binds more strongly. Full ranks include the ≥10-fold endpoint and size/charge-matched controls. There is no conversion to kcal/mol of observed binding free energies and no fitting of affinity labels. [Manifest fingerprints](inputs.csv) · [All endpoints](ranks.csv) · [Individual decisions](pairs.csv)

## Pair ranking and structural controls

| Assay and receiver | ≥3-fold pairs | Upstream | Parent | Replacement | Parent ρ / τ-b | Replacement ρ / τ-b |
|---|---:|---:|---:|---:|---:|---:|
| Gal3 1KJL | 403 | 217 | 222 | 228 | 0.021 / 0.027 | 0.041 / 0.040 |
| Gal3 3T1L | 403 | 202 | 196 | 196 | −0.061 / −0.065 | −0.061 / −0.065 |
| Gal3 mean | 403 | 222 | 221 | 225 | 0.020 / 0.026 | 0.038 / 0.039 |
| hCA II 6WQ4 | 7 | 4 | 4 | 4 | 0.10 / 0.00 | 0.10 / 0.00 |
| hCA II 6WQ8 | 7 | 3 | 3 | 3 | 0.20 / 0.20 | 0.20 / 0.20 |
| hCA II mean | 7 | 2 | 2 | 2 | −0.10 / 0.00 | −0.10 / 0.00 |

On Gal3, compound 32 (Kd 1.3 µM) changes from −3.913 to −4.154 kcal/mol on 1KJL. The mean score now orders it correctly against 8, 18, 20 and 47. The 403 comparisons reuse the same 47 ligands. A paired bootstrap drawing whole molecules with replacement, sharing their multiplicities across both receptors and arms, gives a change of +0.99 percentage points with an empirical percentile interval of [0.00, 4.06] across 10,000 draws (seed 614039). This describes resampling within the observed panel without incorporating between-assay, between-scaffold or seed uncertainty. All four new credits occur within the 31-member scaffold (87→91/160). All 243 between-scaffold comparisons together give 134 credits in both arms. The size-only comparator scores 259.5/403. The matched-heavy-atom/charge consensus remains 19/35; 1KJL goes from 15 to 16 and 3T1L remains 20. At ≥10-fold separation the consensus goes from 36 to 37/63; the two matched pairs both remain correct. [Molecule bootstrap](bootstrap.csv) · [Per-scaffold and between-scaffold tallies](galectin_scaffolds.csv) · [Per-compound changes](changes.csv)

On hCA II, gated and ungated replacement produce identical PDBQTs, including three changes in non-leading modes relative to the parent. Both 10-fold endpoints remain 1/5 on the mean score. The panel has zero ≥3-fold same-size/same-charge comparisons. All arms obtain 0/10 Top-1, Top-5 and oracle hits at symmetry-corrected heavy-atom RMSD ≤2 Å (minimum 7.113 Å); the conditional ranking denominator is zero. Every ligand's crystallographic coordinating nitrogen stays farther than 2.5 Å from the retained receptor Zn in every emitted mode (minimum 5.668 Å). The alignment uses chain-A receptor coordinates, fixed atom mappings, symmetry and conjugated-terminal-group equivalences. These results localize the geometric failure to the available output pools. They leave open which unsaved intermediates visited metal-binding geometry. [Every hCA II mode with RMSD and Zn–N distance](caii_geometry.csv)

## Actual computation and fixed-coordinate control

| Panel | Engine | Processes | Wall s | Joint evaluations | Explicit inputs | Explicit optimizer calls | Output modes |
|---|---|---:|---:|---:|---:|---:|---:|
| Gal3 | Upstream | 94 | 134.416 | 1,171,897 | 525 | 553 | 508 |
| Gal3 | Parent | 94 | 15.919 | 1,173,621 | 525 | 560 | 507 |
| Gal3 | Replacement | 94 | 16.551 | 1,171,445 | 525 | 548 | 513 |
| hCA II | Upstream | 10 | 24.211 | 123,728 | 65 | 65 | 65 |
| hCA II | Parent | 10 | 5.360 | 123,801 | 65 | 65 | 65 |
| hCA II | Replacement | 10 | 5.109 | 123,809 | 65 | 65 | 65 |
| hCA II | Ungated | 10 | 5.595 | 123,809 | 65 | 65 | 65 |

All 104 fork comparisons preserve each cell's hunt calls/evaluations, second calls/evaluations, and number of inputs to explicit refinement. The geometric and energy gate replaces one slot in 43/94 Gal3 and 3/10 hCA II cells. On Gal3, 12 fewer repeated explicit BFGS calls at the box boundary account for the difference between equal explicit-input counts and unequal call counts. Downstream filtering and geometry also change the number of emitted output modes. Joint counts are actual objective-with-gradient calls in hunt, second and explicit refinement. They exclude scalar-only final rescoring and distance checks, both included in wall time. The Gal3 candidate took 3.97% more measured wall time while using 0.19% fewer joint evaluations. hCA II took 4.68% less time with eight additional joint calls. Neither target exceeded the prespecified +5% cost escalation threshold. The nominal cap is per trajectory and particular calls can overshoot it. This run does not demonstrate a strictly equal-time increase in Gal3 ranking accuracy. [Cell timings, counters, hashes and replacement decisions](cells.csv) · [CPU-stage totals](cost.csv)

The timing covers one complete CLI invocation per cell through PDBQT output, including map construction. Parent and candidate rotate block order on CPU 7. Official upstream timing uses CPU 6; its passive cost-instrumented replay produced byte-identical PDBQT on all 104 cells. The experimental/source basis is unmodified official `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`, described by Git as `v1.2.7-27-g3c65c0b`; its `src/lib` and `src/main` trees equal those of the v1.2.7 release tag. The fork measurements use passive counters. CPU wall differences relative to the upstream executable reflect the whole program, including preexisting map work. The same 47 Gal3 prepared ligands were reused across arms (105.50 s summed ligand preparation in the source benchmark). A cached-source regeneration of all five hCA II structures on CPU 6 took 14.35 s, with identical ligand, receptor, box, graph and metadata bytes. Neither candidate requires additional chemical preprocessing.

Scoring each of the 507 Gal3 and 65 hCA II parent poses again with the original input topology, identical printed coordinates and common unbound energy 0 produced identical ten-line decompositions from parent and candidate in all 572 cases, including all 104 first poses. This checks fixed-input output to CLI precision, whereas the affinity table uses the normal emitted native score and ligand-specific reference. Each fixed-pose hash, score and decomposition digest is in `fixed_scores.csv` in the archive. [Production PDBQTs, logs, passive replays' counters, binary hashes and fixed-score audit](runs.tar.xz)

## Reproduce

Run from repository root. The reference analysis used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and RDKit 2025.09.3. Rebuilding prepared hCA II inputs additionally requires the tools documented in its preparation recipe. Gal3 inputs are extracted from the existing archive. For the fork binaries run `core/build.py parent`, `candidate` and `ungated` in separate build directories. For official and passive upstream, use the pinned upstream source and `pair_gradient/upstream_counts.patch` as in `pair_gradient/build_arms.py`. The examples use already built executables in environment variables.

```bash
D=benchmarks/accepted_retention/affinity
W=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python benchmarks/cutoff_switch/hcaii/prep.py --work "$W/ca"
for panel in CAII Gal3; do
  python "$D/run.py" upstream "$panel" "$W" --cpu 6 \
    --upstream "$UPSTREAM" --upstream-counts "$UPSTREAM_COUNTS"
  extra=()
  if [[ $panel == CAII ]]; then extra=(--ungated "$UNGATED"); fi
  python "$D/run.py" fixed "$panel" "$W" --cpu 7 \
    --parent "$PARENT" --candidate "$CANDIDATE" "${extra[@]}"
  python "$D/run.py" matched "$panel" "$W" --cpu 7 \
    --parent "$PARENT" --candidate "$CANDIDATE"
done
python "$D/check_scores.py" "$W" "$PARENT" "$CANDIDATE" --cpu 6
python "$D/evaluate.py" "$W" "$W/tables"
python "$D/geometry.py" "$W" "$W/tables"
```

For a replay of the published tables using saved poses, initialize the same work directory and hCA II fixtures, extract `runs.tar.xz` into it, then run only the last two analysis commands. Recomputing the fixed-coordinate scores additionally needs the two binaries. The archive contains 322 complete production jobs plus upstream count logs and single-point results. Common source input files are referenced in their original benchmark directories.
