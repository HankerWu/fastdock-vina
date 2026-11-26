# Restricted rigid motion and PDBQT boundary replay

`vina --score_only --receptor REC --ligand PREP --rigid_window_pool POOL` evaluates independent rigid translations and rotations around every MODEL, with the chemical graph of PREP. An optional `--ligand_bonds FILE` supplies its calibrated complete graph. Each Cartesian translation is uniform in a 1 Å ball. Rotations are Haar-uniform inside 15°, around the chemical heavy-atom centroid. The two draws are independent under the target measure. The ligand-local frame makes finite-grid orientations follow nondegenerate input geometries. The centroid and frame for each pose are held fixed across probes. The linear- and single-atom cases use a fallback frame.

The evaluator uses the native explicit rigid receptor without a search-box wall. Pair energies, cutoff and the per-atom curl follow the selected native Vina parameters. Multiply the change in intermolecular energy by the native conformation-independent slope, then average its Boltzmann factor at 298.15 K. Intramolecular energy cancels during rigid motion. The output gives raw centre INTER and INTRA, slope, full correction, the nested half-grid correction, the outer-scale-only diagnostic and the effective number of weighted terms. Add the full correction to a compatible central score. Arbitrary zeroes across alternate ligand preparations can invalidate that addition. This integral omits ligand torsion, solvent motion, receptor internal states and the free-solution ensemble. The local windows can overlap; counting their volumes as independent binding modes would double-count some coordinates.

The numerical design uses three equally populated nested proposals with radius and angular limits scaled by 1, 1/2 and 1/4. For scale s the relative volume is s³(s a − sin(s a))/(a − sin(a)), a = 15° in radians. Weight each proposal point by the reciprocal of the mixture density relative to the outer window, then normalize by the summed weights. A positive correction means that the local Boltzmann-averaged score is less favorable than its central score. Negative corrections are allowed. The half-grid discrepancy and the effective sample count diagnose numerical stability. They are not confidence intervals. A fixed radius and temperature across all complexes test a single hypothesis, without molecule-specific calibration. Mechanical or standard-state binding free energies require further terms.

## Reproduction

From the repository root, unpack the accompanying experimental archive with its fdv-* directories under BASE. Set BIN to the absolute path of the compiled branch binary. Example commands (Bash):

```sh
python tools/build.py --output-dir build/physics --jobs 4
BIN="$PWD/build/physics/vina"
mkdir -p "$BASE/window-output"
for panel in hcaii brd4 t4l alkyl fxa fxa_pair bovine_ca trypsin; do
  python benchmarks/rigid_window/run.py "$panel" "$BASE" "$BIN" "$BASE/window-output/$panel.csv" --samples 2048 --jobs 4
done
python benchmarks/rigid_window/evaluate.py "$BASE/window-output"
python benchmarks/rigid_window/replay_searches.py "$BASE" "$BIN" "$BASE/replay" --jobs 4
```

A score file uses one row for each model. Results are paired within the same receiver, ligand preparation, random seed and candidate pool. Solvent flexibility and macrocycle closure are outside the scope of this diagnostic. The chemically distinct comparison metric retains pairs with at least threefold separation in the same assay, takes the median predicted within-block gap and gives half credit for an exact tie. Receptor structures and seeds are repeats. Only benzamidine itself has a whole-molecule experimental pose reference in the trypsin panel. FXa and bovine CA contribute no whole-ligand structural labels here. The hCA II joint endpoint requires a symmetry-aware heavy-atom RMSD ≤2 Å and the local bound orientation at the pyramidal centre used by the archived graph protocol. Its common and reflected pools share the unbound reference from their joint native minimum.

The C++ reference program builds a broad legacy non_cache grid and directly evaluates each rounded, independently transformed PDBQT, without the swept neighbour shortlist. Build it using the objects and Boost configuration used for BIN (omit main.o, include the same Boost linker flags as the engine). Then run:

```sh
python benchmarks/rigid_window/verify.py "$BASE" "$BIN" /path/to/reference "$BASE/window-output"
```

The checked tables and complete raw scores are archived separately. `summary.csv` and `numerics.csv` refer to 2048 points per scale, `sensitivity.csv` compares 128, 2048 and 8192 as available. The two diagnostic boundary re-reads are in `boundary.csv`. Replays compare all originally returned models using atomwise displacements, report additions and losses separately, and preserve both SHA-256 digests for every pool. Reproduction of optional probes requires Python, NumPy and SciPy. The independent reference uses a second implementation of the sampling distribution, rotation, full explicit potential evaluation and log-sum-exp.

The independent analytic.cpp checks geometric and numerical limits. Compile with C++14 and the same Boost headers. Run `python benchmarks/rigid_window/cli_checks.py BIN` for PDBQT reader coverage. Pass `--replayed "$BASE/replay"` to verify.py to independently re-score repaired poses through the normal entry point.
