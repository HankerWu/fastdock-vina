# Can one receptor-directed hydrogen bond rescue the BRD4 scaffold-transfer failures?

Human BRD4(1) recognizes acetyllysine through the conserved Asn140 side chain. Crystal structures of kinase inhibitors in this pocket report interactions between several ligand chemotypes and Asn140. We tested a receptor-side alternative after graph-based transfer failed on SB-202190 (4O77) and BI-2536 (4OGI). The test and its 2.90 Å donor–acceptor separation were specified before generating this contact arm. These results are retrospective and target-specific. [Crystal structures and recognition chemistry](https://doi.org/10.1021/cb500072z) · [ITC affinities and crystal structures](https://doi.org/10.1038/nchembio.1471)

For each prepared receptor (3MXF, 4OGI, 4OGJ), we placed a point 2.90 Å from Asn140 ND2 along each prepared N–H unit vector (HD21 and HD22). For each query we paired every PDBQT NA, OA or SA atom with each point. Each pairing contributes one alternative with an equal prior. The resulting sets contain 10, 6, 10 and 12 assignments for 4O76, 4O77, 4OGI and 4OGJ respectively. These hydrogen positions are preparation-model coordinates. The guide reads the prepared receptor and query acceptor types, and uses no query crystal pose or other ligand template. The primary experiment does not assign a specific ligand acceptor, allow a bridging water, or reorient the rigid Asn side chain.

All three arms use the [matched prepared inputs](../cutoff_switch/brd4/README.md), four seeds (314159, 271828, 161803, 141421), exhaustiveness 4, a nominal 50,000 evaluations per trajectory, one CPU per run, late native refinement, nine requested modes and a 100 kcal/mol window. The OR assignments use the existing smooth search bias on all trajectories. Every final score comes from the native function with the search guide released. The full scaffold alternative and native controls come from the same four-seed mapping matrix and use the same search budget. There are twelve repeated receptor/seed conditions per ligand, including cognate receptor protein coordinates for two of the four compounds.

| Search guide | First pose ≤2 Å / 48 | Pool contains a pose ≤2 Å / 48 | Median first-pose RMSD in Å | Correct ≥3-fold Kd decisions / 60 |
|---|---:|---:|---:|---:|
| Native | 1 | 2 | 7.107 | 20 |
| All scaffold mappings | 19 | 24 | 6.352 | 25 |
| Asn140 point alternatives | 0 | 10 | 5.195 | 15 |

| Experimental complex | Native first/pool | Scaffold first/pool | Asn140 first/pool |
|---|---:|---:|---:|
| 4O76 (TG-101209) | 1 / 1 | 12 / 12 | 0 / 8 |
| 4O77 (SB-202190) | 0 / 1 | 0 / 0 | 0 / 0 |
| 4OGI (BI-2536) | 0 / 0 | 0 / 0 | 0 / 0 |
| 4OGJ (TG-101348) | 0 / 0 | 7 / 12 | 0 / 2 |

The two ligands that defeated the scaffold guide still gave zero successes among 24 first poses and 24 returned pools with the receptor-side guide. For 4O77 the nearest eligible crystal atom falls 2.74–3.46 Å away from the nearest constructed point across the three protein frames. For 4OGI the corresponding distances are 0.58–1.41 Å. Even a close ideal point in the cognate crystal failed to recover the global pose. Among the 427 contact-guided output poses, 342 have an eligible heavy atom within 1.0 Å of a constructed point, including 327 with whole-ligand RMSD >2 Å. Of the 48 first poses, 29 meet this distance diagnostic and all 48 exceed 2 Å global RMSD. One contact leaves ligand orientation and remaining torsions highly ambiguous. The 1.0 Å distance is a descriptive post-run diagnostic.

The four ITC Kd values are 37, 123, 164 and 3400 nM for 4OGI, 4O76, 4OGJ and 4O77. The measurements use BRD4(1), 15 °C, 50 mM HEPES pH 7.4 and 150 mM NaCl. Five distinct pairs differ by at least threefold. The same pairs recur in twelve conditions. The geometry and affinity data do not support enabling this prior by default.

## Reproduction

Prepare the BRD4 raw/ and prepared/ directories following the [existing panel](../cutoff_switch/brd4/README.md). Generate native and multi reference outputs following [the mapping benchmark](../mapping_ensemble/dock_matrix.py) using the four seeds listed above, then run the following from the repository root.

~~~bash
python benchmarks/receptor_contact/experiment.py /tmp/fdv-brd4/prepared /tmp/fdv-receptor-contact build/ensemble/vina --jobs 4
python benchmarks/receptor_contact/score.py /tmp/fdv-brd4 /tmp/fdv-receptor-contact /tmp/fdv-map-brd4 /tmp/fdv-receptor-contact/results
python benchmarks/receptor_contact/check_pairs.py /tmp/fdv-receptor-contact /tmp/fdv-map-brd4 /tmp/fdv-receptor-contact/results
~~~

The first script creates anchors and performs the 48 new searches. The second computes receptor-aligned symmetry-corrected RMSD from the deposited chemical graphs and audits every model against the independent RCSB ligand-instance SDF/Meeko/RDKit CalcRMS path. Across the 144 runs and 1,276 models, the maximum absolute RMSD discrepancy was 0.00000046 Å. The third script reproduces all 180 pair margins and credits directly from the released score lines and the frozen ITC labels.

[Summary](results/summary.csv) · [Each paired cell](results/cells.csv) · [Each pose and contact displacement](results/modes.csv) · [Crystal-to-point displacement](results/geometry.csv) · [Receptor-derived sites](results/sites.csv) · [Pair decisions](results/pairs.csv) · [Independent RMSD check](results/rmsd_check.txt) · [Independent ranking check](results/pair_check.json)
