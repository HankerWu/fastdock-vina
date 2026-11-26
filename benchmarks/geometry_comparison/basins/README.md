# Is a crystal-like pose reachable, and does the potential favor it?

All eight primary eligible human cognates from the rotation experiment are
included. For each ligand use exactly the serialized PDBQT tree and coordinates.
Move only its BRANCH torsions and overall rigid transform to minimize receptor-frame
RMSD to the first allowed reference mapping. Rigid alignment is solved by SVD
inside a least-squares fit from six deterministic torsion starts. This is a
crystal-informed diagnostic. The fitted distance is an achieved value; the
global optimum is unknown. The other two starts are the first emitted poses from
the dry analytic searches at seeds 63419 and 63421.

Give each start the same single native `--local_only` call. For every initial
and relaxed pose, transfer only its coordinates into the *original prepared
ligand* and score with `--score_only --pose_coordinates ... --unbound_energy 0`.
The energy in the table is intermolecular plus intramolecular Vina energy on that
one graph. Comparisons across starts have the same ligand and receiver. The
experimental pose never participates in predictive ranking.

| Receiver | Relaxed projected RMSD (Å) | Projected total | Wrong rank-one totals after local relaxation |
|---|---:|---:|---:|
| PHIP 5S8I | 0.339 | −5.749 | −6.601, −4.973 |
| BAZ2A 7BLA | 0.267 | −10.624 | −12.128, −11.666 |
| BPTF 7F5D | 0.243 | −10.248 | −8.612, −8.964 |
| BCL2L1 7JGW | 0.665 | −20.167 | −18.302, −16.488 |
| HSP90AA1 7LT0 | 0.698 | −11.435 | −11.760, −13.706 |
| BCL6 7LZQ | 0.396 | −5.533 | −7.310, −6.976 |
| BRDT 7MRH | 0.245 | −9.650 | −9.320, −8.989 |
| USP5 7MS7 | 1.264 | −8.198 | −7.385, −8.046 |

The projected basin wins both seeds in four cases, loses both in three, and
splits in PHIP. A separate scoring replay reissued all 96 CLI score calls
(original ligand graph vs reparsed local starting graph), recovered all 48
inter/intra rows at output precision and agreed on receptor-frame RMSD to
1.8e-15 Å. See the receptor-assembly audit before interpreting the BCL6 row.

Run `python basins/run.py PREPARED ROTATION_WORK VINA OUT` from the parent
benchmark directory after `rotation/run.py`. Outputs are energies.csv and
the projection/relaxation coordinates.

For independent verification from released artifacts, extract
`starts_and_relaxations.tar.xz` to POSES, then run
`python basins/verify.py PREPARED POSES VINA`.
