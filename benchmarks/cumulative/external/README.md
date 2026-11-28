# Prospectively selected human structures

One of three prespecified structures was eligible for the historical dry, single-chain known-site protocol. EPHX2 and PARP1 were previously studied families. The three structure identifiers were new to the reports and tables inspected at registration. Ligand starting coordinates came from the CCD graph. The receiver ligand located the pocket; the RMSD reference stayed out of docking. These observations are now available to later developers.

| PDB | Query | Outcome |
|---|---|---|
| 4OCZ | EPHX2 ligand 2RU | Eligible. Closest omitted nonwater atom was a phosphate oxygen 48.719 Å away. |
| 6AUM | EPHX2 ligand BXV | Excluded. Ligand O25 was 4.129 Å from an omitted chloride. |
| 4R6E | PARP1 ligand 3JD (niraparib) | Excluded. Ligand CAG was 3.526 Å from an omitted glycerol C2. |

| 4OCZ seed | Official Top-1 / pool best (Å) | Parent | Initial Pro | Final Pro |
|---|---:|---:|---:|---:|
| 85231 | 9.088 / 7.561 | 9.088 / 7.561 | 1.762 / 1.762 | 1.762 / 1.762 |
| 85237 | 1.282 / 1.282 | 1.282 / 1.282 | 10.310 / 9.801 | 10.310 / 9.801 |

The final Pro improved one seed by 7.326 Å and worsened the other by 9.028 Å against both fixed baselines. Its change in mean Top-1 was +0.851 Å. Both final pools differ in bytes from initial Pro following the native wall fix. The 25-heavy-atom structural check enumerates 48 stereochemistry-compatible mappings and independently reads every rank directly from its published pool. See [audit.json](audit.json), [per-mode coordinates and scores](results/), [preparation and costs](preparation.json) and [the original registration](PROTOCOL.md). The four 4OCZ modes that motivated the wall repair became development observations at that point. Their post-repair numerical replay is a regression check. Search accuracy on a new family remains untested by this extension.

Sources are [4OCZ](https://www.rcsb.org/structure/4OCZ), [6AUM](https://www.rcsb.org/structure/6AUM), and [4R6E](https://www.rcsb.org/structure/4R6E). inputs.tar.xz contains the experimental PDB and component dictionaries, selected receptor, graph-generated ligand and excluded input records. The preparation clocks include downloading; extraction is an alternative to regeneration. All scoring, search and refinement used a fixed receptor and chemical state. No equilibrium populations are inferred.

## Reproduction

Build the four executables with the [four-version build script](../../geometry_comparison/upstream/four_versions/build.sh). From the repository root, with an empty work directory:

```sh
D=benchmarks/cumulative/external
for arm in official parent stable candidate; do
  python "$D/run.py" "$arm" "$builds/$arm/vina" "$work/runs/$arm"
done
python "$D/collect.py" "$work/runs" "$work/collected"
mkdir "$work/fixtures"
tar xf "$D/inputs.tar.xz" -C "$work/fixtures"
python "$D/audit.py" "$work/fixtures" "$work/runs"
```

Regenerate the chemistry with `OPENMM_CPU_THREADS=1 python "$D/prepare.py" "$work/new-inputs"`; the recorded chemistry environment was RDKit 2026.9.2, Meeko 0.8.0, Gemmi 0.7.5, OpenMM 8.6.1 and PDBFixer 1.12.0. The published input hashes are authoritative for the fixed benchmark. `intermediate/` retains the first candidate (673f226) separately.
