# Fixed four-version affinity panel

The primary matrix reuses every compound and prepared input from [the existing measured affinity panel](../../geometry_comparison/affinity/README.md). The endpoint is the direction and ranking of exact Ki or Kd within one target, assay and receptor. These are five human targets, four kinase Kd assays plus two distinct human BACE1 Ki series, seven receptor states and nine assay–receptor strata. Forty-seven assay–compound rows give 70 assay–compound–receptor cells and 140 runs per version at two registered seeds. Some kinase inhibitors occur across targets. The panel has been used in earlier work. The [original study for kinase Kd](https://doi.org/10.1038/nbt.1990) and both BACE1 series ([Cumming et al.](https://doi.org/10.1016/j.bmcl.2012.02.013), [Stamford et al.](https://doi.org/10.1021/ml3001165)) provide different experimental contexts.

## Longitudinal readout

On the complete nominal matrix, all 140 official and genuine FastDock parent PDBQT files are byte-identical. All 140 first-ranked model blocks of the final Pro at 5847b1b are byte-identical to the starting Pro at f46aa37. Thus the official and parent columns for rank metrics coincide, as do the two Pro columns. Lower score predicts stronger binding. The table gives each version's Spearman ρ / Kendall τ-b, followed by correctly ranked compound pairs separated at least threefold. The fold threshold is a proxy from the historical protocol. Pair observations share ligand identities and cannot be treated as independent experimental replicates.

| Assay and receptor | Compounds | Official or parent ρ / τ-b | Starting or final Pro ρ / τ-b | Official or parent ≥3× | Starting or final Pro ≥3× |
|---|---:|---:|---:|---:|---:|
| ALK Kd, 7BTT | 8 | 0.810 / 0.571 | 0.357 / 0.214 | 17/19 | 12/19 |
| BACE1 Ki CHEMBL2015090, 4DJW | 6 | −0.829 / −0.600 | −0.886 / −0.733 | 2/14 | 1/14 |
| BACE1 Ki CHEMBL2015090, 4GID | 6 | −0.486 / −0.333 | 0.771 / 0.600 | 4/14 | 11/14 |
| BACE1 Ki CHEMBL2155063, 4DJW | 8 | −0.762 / −0.571 | −0.595 / −0.357 | 2/22 | 5/22 |
| BACE1 Ki CHEMBL2155063, 4GID | 8 | 0.286 / 0.143 | −0.214 / 0.000 | 14/22 | 10/22 |
| MAPK14 Kd, 1KV2 | 9 | 0.267 / 0.278 | 0.150 / 0.167 | 19/29 | 17/29 |
| MAPK14 Kd, 3GCP | 9 | −0.033 / 0.056 | −0.250 / −0.222 | 15/29 | 11/29 |
| RET Kd, 7DUA | 8 | −0.238 / −0.071 | 0.071 / 0.000 | 8/19 | 11/19 |
| TNK2 Kd, 7KP6 | 8 | 0.381 / 0.286 | 0.381 / 0.286 | 12/19 | 12/19 |

Relative to either fixed historical baseline, Pro gained 41 directions and lost 44 across 187 stratum-specific comparisons. Three strata improved in ≥3× accuracy, five regressed and one retained every direction. The second BACE1 series illustrates a size-dependent failure in one receptor: on 4DJW, measured pKi has ρ = −0.67 against heavy-atom count while Pro predicts +0.88, and the measured/predicted correlations for TORSDOF are −0.68/+0.68. On 4GID, the predicted heavy-atom correlation for exactly these molecules is +0.02. No per-receptor fitting or microstate selection was performed. The observations do not identify the relative contributions of conformational selection, search and the fixed dry potential. Formal charge is constant within these series, precluding a charge-effect estimate within them. All pair identities, comparisons and chemical descriptors are exported in [pairs.csv](results/pairs.csv), [ligands.csv](results/ligands.csv) and [assays.csv](results/assays.csv).

At the registered native CLI scientific defaults, all four versions got only 1/14 ≥3× pairs correct on the preselected six BACE1 ligands in 4DJW (one seed). Verified stock-source counter replays of all six **nominal-budget** first-seed ligands reproduced every stock PDBQT byte and yielded 147,230 summed hunt + second + explicit joint energy/gradient calls in each of official and parent. Their extra observer CPU times were 16.50 and 6.71 s respectively. Uninstrumented full-matrix counters for these versions remain unmeasured. The Pro counters apply to joint BFGS evaluations; energy-only final ranking and map precomputation are included in elapsed/CPU times. The clocks do not separately attribute those stages. There is no inference of a bound-versus-solution free energy.

## Cost and validation

Nominal cost over 140 fresh process launches, including parsing, grid generation, search, refinement and serialization, measured as sums of per-job CPU / wall seconds: official 360.83 / 465.84, parent 119.63 / 125.47, stable 57.01 / 58.88, final 74.97 / 96.93. Peak per-process RSS is in [costs.csv](results/costs.csv). The timing windows had different concurrent workloads and the final run used two concurrent single-CPU processes. Summed per-job wall times are not batch makespans. The deterministic fixed-budget top-score comparison is independent of these fluctuations.

The first CPU-control target came from the six single-seed stable nominal runs (1.962 s total). The four observed confirmation times were 12.79, 6.79, 4.03 and 2.91 s. None reached the ±15% window. In the registered follow-on experiment, the six-ligand target was 14.5 CPU seconds. All four final confirmations were inside the window. The relatively smaller cap for the official and larger caps for Pro absorb their different setup and search costs. These are six BACE1 analogues in one previously examined receptor and assay.

| Version | Calibrated cap per trajectory | Observed six-ligand CPU s | ρ | τ-b | Correct ≥3× |
|---|---:|---:|---:|---:|---:|
| Official | 57,600 | 13.975 | −0.829 | −0.600 | 2/14 |
| FastDock parent | 210,200 | 16.627 | −0.771 | −0.600 | 3/14 |
| Starting Pro | 299,500 | 16.489 | −0.771 | −0.600 | 3/14 |
| Final Pro | 384,000 | 14.077 | −0.771 | −0.600 | 3/14 |

The adaptive caps use CPU clocks only. In [common_cost.json](results/common_cost.json) every anchor and confirmation has its actual cost. Tables contain separate rows labeled common for the final confirmations; those rows link to the same saved runs and should be counted only once when totaling research costs. The receipt of the final assay audit counts 728 distinct raw output files plus 24 reused final-confirmation references, yielding 752 ledger rows across all stages, and validates 1,448 pair rows. The 12 separate stock counter replays have their own CPU and coordinates.

Input archive verification and extraction took 0.133 CPU / 0.133 wall seconds. A full replay of the existing preparation pipeline from its cached raw records took 43.55 process CPU / 79.18 wall seconds (peak RSS 193,180 KiB), including 37.86 CPU seconds for 32 graph-based ligand preparations and 3.88 CPU seconds in five receptor-preparation calls. Two MAPK14 receptors are taken from their earlier prepared archive: it would cost additional, unmeasured time to regenerate those two from deposition. All 39 delivered PDBQTs, seven box files, panel and exclusion tables matched byte for byte. Five auxiliary receptor metadata files now record additional gap information. No changed metadata supplied the docking coordinates. A preliminary timing run regenerated the preparations but failed only when comparing new debug-coordinate files to the old archive; it was repeated with proper handling of generated-only files.

The final-domain fix affected 21 whole nominal Pro output files at lower ranks. Of 1,205 shared output ranks, 28 model blocks differed and 12 printed scores changed; no rank-one coordinate block or score changed. The largest same-rank score difference was 26.671 kcal/mol. Rank-aligned blocks can represent different minima. [Pool comparison](results/pool_changes.csv) provides coverage and the raw archive allows independent inspection.

## Rerun

Use the [registered panel](PROTOCOL.md), [counter addendum](COUNTER_PROTOCOL.md), [CPU target addendum](COMMON_COST.md) and [four source identities and build options](../../geometry_comparison/upstream/four_versions/PROTOCOL.md). The costs above refer to cold Vina processes on cached inputs. Use Linux, Python 3.12, NumPy, SciPy and xlrd. The optional preparation replay additionally needs the dependency versions listed in the original panel. Set OFFICIAL, PARENT, STABLE and CANDIDATE to the four exact executables.

~~~sh
D=benchmarks/cumulative/affinity
W=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 python "$D/run.py" "$W" --prepare
for engine in official parent stable candidate; do
    case "$engine" in
        official) exe=$OFFICIAL ;;
        parent) exe=$PARENT ;;
        stable) exe=$STABLE ;;
        candidate) exe=$CANDIDATE ;;
    esac
    python "$D/run.py" "$W" --engine "$engine" --binary "$exe" --stage nominal low default
done
for engine in official parent stable candidate; do
    case "$engine" in
        official) exe=$OFFICIAL ;;
        parent) exe=$PARENT ;;
        stable) exe=$STABLE ;;
        candidate) exe=$CANDIDATE ;;
    esac
    python "$D/run.py" "$W" --engine "$engine" --binary "$exe" --stage matched
done
python "$D/common_cost.py" "$W" --official "$OFFICIAL" --parent "$PARENT" --stable "$STABLE" --candidate "$CANDIDATE"
python "$D/evaluate.py" "$W" "$W/results"
python "$D/verify.py" "$W" "$W/results"
python "$D/pool_changes.py" "$W" "$W/results"
~~~

Time-dependent caps will vary on another computer. Each receipt records the actual cap. For offline replay, extract results/raw.tar.xz into a freshly prepared W and run the three analysis commands. To reproduce the optional independent stock counter experiment, use a separate fresh source copy of each stock checkout, run instrument.py SOURCE, compile with the same compiler/flags/Boost, then run replay_counts.py W official OFFICIAL_COUNTER and replay_counts.py W parent PARENT_COUNTER. prepare_cost.py W NEW_SCRATCH RECEIPT_JSON measures preprocessing independently of docking. The prior score-domain-only candidate 673f226 is kept in results/intermediate.tar.xz with its original timing calibration. The release candidate was run separately on all fixed and default cells. Its calibrated CPU comparisons were remeasured against the same predeclared clocks with common_cost.py --replace-candidate, retaining the historical engines' actual records. Neither transaction imports an intermediate score as a final prediction.
