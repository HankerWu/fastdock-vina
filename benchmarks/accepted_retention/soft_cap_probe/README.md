# What misranks the skipped FKBP1A candidates?

Three of 66 frozen human-family runs had a recoverable skipped pose in the native MC trace. These cases were selected using the crystal reference to diagnose candidate retention. On each of their 409 skipped accepted states, we changed the objective at the exact same stored internal conformation. For the two FKBP1A runs, a harder positive-energy cap and then a smooth explicit atomic receptor both preserved the same wrong minimum. In both runs the receptor term favored the closest-to-crystal geometry over that minimum. The ligand-internal term outweighed this advantage. Discarding the internal term also displaced the successful pose in the FABP4 control. [Coverage and intervention](../coverage) · [All skipped conformations](../coverage/results/choices.csv)

## Soft and authentic grid energies

The native global entry point uses `hunt_cap=(10,10,10)` and a second-minimization cap of `(1000,1000,1000)`. For each positive contribution where Vina applies its curl, the energy becomes `v e / (v + e)` and the force is multiplied by `(v / (v + e))²`. Negative contributions pass through. For fixed positive e, raising the cap increases this energy. For the rigid receptor grid the curl applies to each ligand atom's interpolated grid value. Intraligand movable pair contributions are curled separately. Soft hunt energies are available from the original trajectory.

The passive C++ observer owns a private copy of the model and its gradient buffer. It resets that model from the accepted internal conformation, verifies zero coordinate difference against every corresponding live movable heavy atom, and calls `model::eval_deriv` on the same search grid using the authentic cap. It does not minimize any endpoint. Per step it verifies the original soft energy against the separately archived observer. A second independent replay measures receptor and ligand components below.

| Query → receptor (seed 75431) | Skipped | Novel at observation and exit | Lowest soft → lowest authentic grid, RMSD Å | Closest-to-crystal pose's soft → authentic energy gap to winner | Closest pose's rank |
|---|---:|---:|---|---:|---|
| FABP4 2NNQ → 3FR2 | 79 | 37 | #198, 1.686 → #198, 1.686 | +0.3024 → +0.2988 | 13 → 12 |
| FKBP1A 1FKG → 1FKG | 183 | 167 | #150, 4.976 → #150, 4.976 | +0.7709 → +0.7609 | 13 → 13 |
| FKBP1A 1FKH → 1FKG | 147 | 121 | #140, 6.746 → #140, 6.746 | +0.3233 → +0.3395 | 9 → 9 |

Steps are trajectory-local and one-based. For each stored conformation, the reference RMSD uses the archived complete symmetry mapping; eligibility uses the selector's fixed-index distance ≥2 Å. Its distance at observation was captured at full precision. A final distance computed from coordinates rounded individually to 0.001 Å can differ from the unrounded fixed-index RMSD by at most 0.001733 Å. No archived final distance fell within this margin of the threshold. Both online tournament winners survived at exit in all three cases. In the 1FKH cross-run, one of 122 initially novel states was redundant by exit. Counts in these tables are conditional on the retrospectively identified opportunities.

Raw pairwise ordering changed for 21/666, 92/13,861, and 27/7,260 eligible-state pairs after raising the cap. The maximum change to a single skipped state's energy was 3.813, 0.169, and 0.217 kcal/mol across the three runs. These are comparisons among correlated observations along single trajectories. No independent-sample p-value is appropriate. [Raw joint evaluations](results/endpoints.csv) · [Winners and costs](results/cases.csv) · [Pair reversals](results/ordering.csv)

## Grid versus explicit receptor, at the same geometry

We passed a read-only view of Vina's initialized native explicit receptor into an independent diagnostic replay. With the authentic cap, the grid joint value H and the native grid's receptor component G give the shared internal ligand contribution L = H − G. For a rigid receptor and one ligand, L contains all remaining intraligand terms used by this objective. The smooth atomic receptor derivative routine returned the intermolecular component A at identical coordinates. Its corresponding complete atomic objective is A + L. The first skipped state in each trajectory received an additional direct full-model evaluation with the explicit receptor. All three agreed with this reconstruction within 1e−8 kcal/mol. All 409 endpoints were inside the box. Component comparison uses the smooth pair-table energy evaluated by the optimizer throughout. Vina also uses midpoint tables for some reported energy values. Their different numerical convention is outside this comparison.

| Query → receptor | Lowest full explicit objective | Lowest explicit receptor component alone | Closest-to-crystal minus full explicit winner, receptor | Same difference, internal | Sum |
|---|---|---|---:|---:|---:|
| FABP4 2NNQ → 3FR2 | #199, 1.668 Å | #165, 4.574 Å | +0.051 | +0.235 | +0.285 |
| FKBP1A 1FKG → 1FKG | #150, 4.976 Å | #72, 1.921 Å | −0.278 | +0.988 | +0.710 |
| FKBP1A 1FKH → 1FKG | #140, 6.746 Å | #216, 2.193 Å | −0.527 | +0.931 | +0.404 |

Values in the last three columns are in kcal/mol. The closest-to-crystal steps were #193 (1.606 Å), #65 (1.821 Å) and #265 (1.592 Å). Several near-native poses occurred close together along a trajectory. At least one ≤2 Å state ranks first, seventh and eighth, respectively, using the full explicit objective. Favorable internal energy of a geometrically wrong state is therefore a specific scoring tradeoff to investigate for these two FKBP1A runs. This observation alone cannot distinguish physical strain in the crystal-near ligand, bias in Vina's internal movable-pair model, and missing receptor or solvent stabilization. The FABP4 result supplies a counterexample to uniform suppression of the internal term. [Complete component scores](results/atomic_endpoints.csv) · [Energy decompositions and ranks](results/atomic_cases.csv)

## Cost and reproduction

The grid-cap probe adds exactly 79, 183 and 147 full joint value/gradient evaluations, compared with 12,204, 12,311 and 12,322 in the ordinary candidate runs. The independent component probe adds 409 repeated grid joint evaluations, 409 grid receptor evaluations, 409 explicit receptor evaluations and three complete explicit joint value/gradient controls. Both probes include the overhead of a private model, input reconstruction and diagnostic output. Their pose files match every uninstrumented candidate's bytes, native counters reproduce, and the ordinary parent and candidate hashes match the previous frozen six-family replay. The diagnostic evaluations were charged only to these post hoc probes. None was used to score or select a reported predictive pose. [Grid probe logs and outputs](results/runs.tar.xz) · [Component probe logs](results/atomic_runs.tar.xz) · [Component replay timings](results/atomic_cost.csv)

Use [the core builder](../core/README.md) to produce fresh parent and candidate directories in an absolute work path `$W`. The Python dependencies are shared with that benchmark. If Boost uses a separate prefix, append `--boost-prefix PATH` to the two compilation commands.

```bash
B=benchmarks/accepted_retention
python "$B/coverage/build_choice.py" "$W/candidate" "$W/candidate/build/obj" \
    "$W/hard-observer" --patch "$B/soft_cap_probe/hard_grid.patch"
python "$B/soft_cap_probe/measure.py" "$W/hard-replay" "$W/results" \
    "$W/parent/build/vina" "$W/candidate/build/vina" "$W/hard-observer/observer"
python "$B/soft_cap_probe/build_atomic.py" "$W/candidate" "$W/candidate/build/obj" "$W/atomic-observer"
python "$B/soft_cap_probe/atomic.py" "$W/results" --work "$W/atomic-replay" \
    --candidate "$W/candidate/build/vina" --observer "$W/atomic-observer/observer"
python "$B/soft_cap_probe/audit.py" "$W/results"
python "$B/soft_cap_probe/atomic.py" "$W/results"
```

The last two commands can also be run directly against the archived `results/` directory to recheck the stored energies, ranks, pose digests and native counters from raw logs. The algorithm under study receives the receptor and query only. Ligand reference coordinates enter the post hoc evaluation and the choice of these three diagnostic cells.
