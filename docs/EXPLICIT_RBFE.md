# Explicit-solvent relative binding free energy

The Pro RBFE tool runs separately from Vina's search and scoring. It takes a JSON configuration and pins an edge list before comparing explicit-solvent simulation records with same-assay measurements. Use only comparable measurements and verify ligand atom graph, pose and net charge before reading archive estimates.

## Reanalyse archived simulation repeats

The tool consumes an OpenFE benchmark export with one result row per directed edge and per-repeat complex and water legs. It reads the published MBAR leg estimates and computes the difference between the two leg means. Reports preserve individual errors, minimum MBAR overlaps, replica mixing statistics, equilibration and effective production iterations. The result also retains the archive uncertainty and reports independent-leg quadrature estimates for replica SEM and within-replica MBAR SEM. Raw trajectories and reduced energies are needed to re-estimate MBAR itself.

```bash
tool=/path/to/engine/tools/explicit_rbfe.py
python "$tool" freeze --config "$config" --selection "$frozen"
python "$tool" compare --config "$config" --selection "$frozen" --output "$comparison"
```

The comparison config points to the OpenFE result export, the prepared OpenFE ligand SDF and receptor, the old Ki benchmark's frozen inputs and docking pools, and the old fixed-pose score table. Relative paths resolve against the configuration directory. The comparison stage refuses a changed configuration or archive. The chemistry comparison checks exact canonical isomeric SMILES and formal charge. It aligns the reference receptor Cα atoms by chain and residue name/number, then evaluates symmetry-aware pose RMSD without a second fit on the ligand. The eligibility output records excluded and selected edges.

Both computed and experimental differences are defined for the directed A→B edge. The binding difference is ΔG(complex) − ΔG(water), and the experimental reference is RT ln(Ki_B/Ki_A). Smaller values favor B. The latter uses the configured comparison temperature and a same-assay Ki as a relative binding proxy.

## Run new complex and water legs with OpenFE

These commands use the OpenFE Relative Hybrid Topology Protocol through its standard planner and runner. Configure OpenFE in its supported environment and choose the device/scheduler there. The supplied OpenFE 1.8.0 planner uses its own default force field and simulation settings. Its planning YAML accepts planner, mapper, solvent component and charge choices. Use explicit protocol settings and serialized Transformations from OpenFE's Python API for a match to a different archived force-field and solvation setup.

```bash
python "$tool" prepare-openfe --config "$config" --selection "$frozen" --output "$bundles"
pair="$bundles/60_to_65"
python "$tool" plan-openfe --pair-dir "$pair" --openfe "$OPENFE" --repeats 3
```

The input preparation stage preserves original atom coordinates and preassigned charges in each two-ligand SDF. The planner must produce exactly one complex JSON and one solvent JSON. Inspect both state endpoints before proceeding. A planner's A/B direction controls the sign, so reverse the sign if it orders the endpoints oppositely to the comparison.

Run a real-input protocol smoke with its own reduced-budget Transformation files before production. Give each full leg its own fresh output directory, and run in separate tasks that can report their allocated device time and terminal state.

```bash
python "$tool" run-openfe --transformation "$complex_json" \
  --output "$pair/simulation/complex" --openfe "$OPENFE"
python "$tool" run-openfe --transformation "$solvent_json" \
  --output "$pair/simulation/solvent" --openfe "$OPENFE"
complex_result="$pair/simulation/complex/openfe_result.json"
solvent_result="$pair/simulation/solvent/openfe_result.json"
"$OPENFE" gather "$complex_result" "$solvent_result" --report raw -o "$pair/raw_legs.tsv"
"$OPENFE" gather "$complex_result" "$solvent_result" --report ddg -o "$pair/binding_differences.tsv"
```

Create the simulation parent before launching the separate legs. The per-leg output contains the quickrun command, execution status and actual wall time along with OpenFE's result and working files. Retain each leg's trajectories, checkpoints and analysis artifacts. If a stage fails, inspect that leg and follow OpenFE's documented resume flow with its existing execution files. A missing or unfinished leg leaves the edge incomplete. [OpenFE planning](https://docs.openfree.energy/en/stable/tutorials/rbfe_cli_tutorial.html) · [Quickrun outputs](https://docs.openfree.energy/en/latest/how_to_guides/run_quickrun.html)

For a directed smoke or an explicitly shortened profile, use the OpenFE
environment's Python so the native OpenFE API can serialize the two directed
Transformations. The profile sets salt concentration, repeats, iteration,
equilibration, production and minimization budgets. Checkpoint and position
frequencies follow the iteration. The protocol otherwise uses OpenFE's
native hybrid-topology defaults.

```bash
"$OPENFE_PYTHON" "$tool" plan-openfe-hybrid --pair-dir "$pair" \
  --config "$smoke_profile" --platform cpu --output "$pair/smoke_protocol"
python "$tool" run-openfe --transformation "$pair/smoke_protocol/complex.json" \
  --output "$pair/smoke_complex" --openfe "$OPENFE"
python "$tool" run-openfe --transformation "$pair/smoke_protocol/solvent.json" \
  --output "$pair/smoke_solvent" --openfe "$OPENFE"
```

A short run checks that the prepared input, settings, trajectory and analysis
work together. Interpret production binding differences using independent
replicas and adequate sampling of both legs.

If an executor exits during analysis after writing its simulation and
checkpoint, preserve that attempt. Reuse the saved reduced potentials and
coordinates with the same OpenFE installation. This invokes OpenFE's native
MBAR and structural analysis in a new directory and links the original inputs.

```bash
"$OPENFE_PYTHON" "$tool" analyse-openfe-leg --attempt "$attempt_dir" \
  --output "$fresh_analysis_dir"
```

Use the saved sampling files from a completed production attempt. The
`run-openfe` entry adds the directory containing the selected OpenFE
executable to its child PATH so OpenFE can find the corresponding
`openfe_analysis`.
