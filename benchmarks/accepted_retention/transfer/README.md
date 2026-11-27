# Independent one-slot retention test on NUDT1 and CHIT1

Gated substitution changed the final pool in 20 of 32 paired runs. All four engines returned zero poses within 2 Å, including the 24 runs with clash-free mapped references. The only change in the first pose occurred on a clashing NUDT1 cross-dock. These four compounds provide no independent evidence of pose-recovery gains. Conditional top-1 accuracy among successful pools is undefined for all arms.

## Deposited inputs

The [frozen protocol](PROTOCOL.md) preceded docking at commit `3a14147`. Each family contributes two ligand graphs, two cognate redocks and two reciprocal cross-docks; each direction uses both seeds and both budgets. The CHIT1 compounds share an aminotriazole-piperidine motif. Treat the two proteins as the generalization units. Multiple seeds, caps and receptors are repeated measurements.

| Human target | Structure / primary publication | Chain | Ligand | Resolution, Å | Ligand minimum occupancy | Fixed crystallization additives within 8 Å |
|---|---|---|---|---:|---:|---|
| NUDT1 | [4N1U](https://www.rcsb.org/structure/4N1U), [Gad et al.](https://doi.org/10.1038/nature13181) | A | 2GE, TH588 | 1.60 | 1.00 | Two sulfates at 3.73 and 5.38 Å |
| NUDT1 | [4C9X](https://www.rcsb.org/structure/4C9X), [Huber et al.](https://doi.org/10.1038/nature13194) | A | VHS, (S)-crizotinib | 1.20 | 0.75 | Chloride at 3.68 Å, occupancy 0.75 |
| CHIT1 | [5NRA](https://www.rcsb.org/structure/5NRA), [Mazur et al.](https://doi.org/10.1021/acs.jmedchem.7b01051) | A | 95K, compound 7g | 1.27 | 0.60 | Glycerol at 6.33 Å, occupancy 0.55 |
| CHIT1 | [6ZE8](https://www.rcsb.org/structure/6ZE8), [Koralewski et al.](https://doi.org/10.1021/acs.jmedchem.0c01179) | C | QGB, OATD-01 | 1.50 | 1.00 | None |

All arms receive the same prepared receiver and a single ETKDG/MMFF conformer from the neutral CCD chemical graph (seed 43103). The receiver crystal ligand defines the pocket center, shared in cognate and cross runs. Docking sees neither the bound query geometry nor a restraint to the measured ligand. Protein preparation uses PDBFixer and Meeko. All repaired heavy atoms are further than 9.6 Å from the ligand. Same-chain buffer positions are kept rigid, with CCD atom types and measured heavy coordinates. Pocket waters, symmetry mates, fractional-occupancy alternatives and microstate populations are omitted. Vina uses its native X-Score atom-type interactions for these additives; these fixed models cannot resolve ion screening or fractional occupancy. A packing neighbor approaches 5NRA to 4.12 Å. The deposited neutral ligand states, particularly the basic CHIT1 compounds, leave pH-dependent protonation as a possible limitation.

The protein alignment uses matching residue number, insertion code and residue identity in both 14 Å pockets. Reciprocal CHIT1 frames align 62 Cα atoms at 0.209 Å RMSD. Reciprocal NUDT1 frames align 52 at 0.697 Å. The TH588 reference in 4C9X has one severe overlap and a 1.993 Å minimum heavy distance. The mapped crizotinib reference in 4N1U has 19 severe overlaps and a 1.182 Å minimum. All four cognates and both CHIT1 transfers pass this test. An overlap is a pair closer than 75% of the summed tabulated radii. Each source's crystallographic and chemical evidence and raw atom coordinates are in [inputs.tar.xz](inputs.tar.xz). See [reference geometry](results/reference_geometry.csv).

## Blind pose results and compute

Each arm received one trajectory and one thread per run, two seeds, nominal search limits of 12,000 and 24,000 joint energy/gradient evaluations, and up to nine returned models. Top-1, top-5 and returned-pool oracle each scored 0/16 for NUDT1 and 0/16 for CHIT1 in every arm. The same zero result holds for each cap. Search and final explicit-refinement counts include work beyond the soft CLI cap. Joint evaluations count a function returning both objective and gradient once. The upstream counts-only binary reproduced each uninstrumented upstream output byte for byte. The official v1.2.7 release reproduced all 32 outputs of source-built upstream.

| Engine | Actual joint evaluations | Explicit inputs / calls | Returned poses | Primary end-to-end wall, s | Three-replay wall, s | Three-replay CPU, s |
|---|---:|---:|---:|---:|---:|---:|
| Source-built Vina | 586,223 | 155 / 163 | 151 | 40.332 | — | — |
| Production Pro parent | 586,657 | 155 / 171 | 150 | 7.789 | 7.145 | 7.110 |
| Gated substitution | 586,604 | 155 / 171 | 150 | 7.572 | 7.176 | 7.137 |
| Unconditional substitution | 586,604 | 155 / 171 | 150 | 8.001 | 7.319 | 7.288 |

Wall includes process startup, parsing, grids, search, final refinement and serialization. Preparation from local source files took 11.65 s once for all arms and reproduced all eight prepared PDBQTs byte for byte. Primary and Pro timing replays used logical CPU 8 of an AMD EPYC 9V74; official release validation used CPU 6. Independent tasks used other cores, so small timing differences have limited precision. Replay totals sum the three-run median within each cell. All 288 replay outputs were byte-identical to that arm's first output. Neither candidate exceeded the prespecified 5% aggregate trigger. Gated wall in the replays was 0.43% above parent. [Detailed timings](results/timing_replays.csv) and [preparation costs](results/preparation_cost.csv).

The gated and unconditional output files agree in every cell. Each candidate checked 2,208 opportunities across 5,995 skipped accepted hunts and replaced a terminal slot in 20/32 cells. The parent and candidates have identical hunt and second-stage counts cell by cell. Across 150 returned modes, 130 retain byte-identical heavy coordinates; 31 of 32 first coordinates also agree. The single changed top-1 is TH588 in 4C9X, second seed at 12,000, from 8.00 to 6.60 Å. Four oracle RMSDs improved and five worsened. The largest increase was compound 7g in its own 5NRA receptor, first seed at 12,000, from 5.21 to 6.09 Å. Every top-1 has zero severe heavy clashes by the specified radii. [Per-cell results](results/cells.csv), [continuous displacements](results/paired_changes.csv), [all modes](results/modes.csv) and [family, cap and direction summaries](results/summary.csv).

## Why the returned pool misses the labels

A post hoc diagnostic initialized parent local optimization at the deposited geometry or at crystal torsions projected onto the actual input's rigid fragments. This uses privileged structural labels and adds 16 local searches and 32 fixed-coordinate rescoring calls outside the blind comparison. Every projected start uses the input topology. The four cognate projected poses finished at 0.39, 0.55, 0.55 and 1.24 Å for 4N1U, 4C9X, 5NRA and 6ZE8, respectively. Seven of eight directions retain a pose within 2 Å. Projected crizotinib in 4N1U ends at 2.12 Å with common-graph inter+intra energy +17.20 kcal/mol, consistent with the severe clashes of the mapped label.

For TH588 redocking in 4N1U, the projected local pose has common-graph inter+intra energy −12.29; the four parent global top-1 energies span −8.12 to −6.94 kcal/mol. For OATD-01 in 6ZE8 the projected local pose yields −11.61 versus −9.76 to −5.65 kcal/mol. Both geometrically successful states are favorable on this model yet absent from its returned blind pools. For compound 7g in 5NRA the projected near-native state scores −10.76, while the four geometrically wrong parent first poses score −12.26 to −12.14. This example exposes an interaction between rigid-fragment geometry and scoring. Relaxation from the deposited 7g geometry retains 0.18 Å and scores −12.55, but that exact internal geometry was unavailable to the blind searches. These diagnostics establish accessible local minima under the stated initializations; intermediate trajectory poses would be required to distinguish failure to visit from failure to retain every feasible basin. [Local results](results/crystal_starts.csv) and [raw output coordinates and logs](runs.tar.xz).

## Reproduction

Use Python 3.12, RDKit 2025.9.3, Meeko 0.7.1, Gemmi 0.7.4, OpenMM 8.6.1, PDBFixer 1.12.0, NumPy and SciPy. The scripts reuse the repository's paired docking, chemistry and audit routines. From repository root, reconstruct the recorded evidence:

```bash
P=benchmarks/accepted_retention/transfer
W=/tmp/retention-transfer-audit
mkdir -p "$W"
(cd "$P" && sha256sum -c SHA256SUMS)
tar -xJf "$P/inputs.tar.xz" -C "$W"
tar -xJf "$P/runs.tar.xz" -C "$W"
cp "$P"/results/* "$W/"
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
python "$P/checks.py" "$W"
python "$P/summary.py" "$W"
cmp "$W/audit.json" "$P/results/audit.json"
cmp "$W/summary.csv" "$P/results/summary.csv"
```

For a fresh run use a different empty `W` and extract only inputs. Set `B` to a directory containing the three Pro builds produced by [the retention builder](../core/build.py). Set `U` to the upstream and upstream-counts build root from [the prior build helper](../../pair_gradient/build_arms.py). Select logical CPU indices available on your machine.

```bash
arms=(--parent "$B/parent/build/vina" --gated "$B/candidate/build/vina" --unconditional "$B/ungated/build/vina")
python "$P/experiment.py" "$W" --cpu 0 --upstream "$U/upstream/build/vina" --counts "$U/upstream-counts/build/vina" "${arms[@]}"
python "$P/experiment.py" "$W" --phase release --cpu 0 --release /path/to/vina-1.2.7
python "$P/experiment.py" "$W" --phase time --cpu 0 "${arms[@]}"
python "$P/timings.py" "$W" --cpu 0 "${arms[@]}"
python "$P/prepare_replay.py" "$W" --cpu 0
python "$P/local_reference.py" "$W" "$B/parent/build/vina"
python "$P/checks.py" "$W"
python "$P/summary.py" "$W"
```

The audit independently rebuilds the chemical graph and stereochemistry from PDB/CCD, places references using a second rigid transform implementation and verifies RMSDs and scores for 752 serialized output models (including release replicas). The maximum RMSD discrepancy is 2.84×10⁻¹⁴ Å. No inferred stereochemical tag changed and the maximum covalent bond-length change relative to the prepared graph is 0.00132 Å. [Audit outputs](results/audit.json).
