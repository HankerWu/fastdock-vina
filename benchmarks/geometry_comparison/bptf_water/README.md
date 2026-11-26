# Transferring ordered BPTF pocket waters across experimental receptors

The BPTF structural series identifies ordered pocket waters. We tested whether
immobilizing their positions in Vina helps recover two ligands in both cognate
and reciprocal structures, and whether it helps rank experimentally measured
binders. This is a rigid-receptor Vina protocol. The hydrated AD4 procedure
places displaceable water probes on the ligand and uses specialized maps, a
different hypothesis. [BPTF structural study](https://pubmed.ncbi.nlm.nih.gov/34515477/) ·
[Vina 1.2 methods](https://doi.org/10.1021/acs.jcim.1c00203)

## Evidence and model

Only water oxygens on chain A of the independent donor crystals 7K6R and
7M2E supplied the initial sites. A sequence-verified alignment of the local
protein backbone transfers water observations into 7F5D or 7RWO. All fits
gave 100% sequence identity over aligned residues; local backbone RMSDs were
0.49–0.59 Å. Require occupancy ≥0.8, support in both distinct donors within
the existing 1.0 Å clustering convention, distance ≤8 Å to the recognition
Asn side-chain ND2, and ≥2.5 Å to the nearest dry-receptor heavy atom. Each
receptor retained three sites. Donor ligand coordinates did not determine
site selection. Pocket centers and Cα reference alignments use the same
known-site protocol as the dry comparator.

Compare dry receptor, appended oxygen with no hydrogen-bond role, and
oxygen typed as donor and acceptor. The second arm retains the Vina
Gaussian and repulsion terms. The polar arm adds bonded HD records for XS
typing; Vina scores its oxygen isotropically. A compiled receptor check
confirmed all six oxygen types and bonds, with every original protein
atom's XS type and bond count preserved. Static occupancy, solvent
rearrangement and water binding free energy are separate questions.

In the pilot each arm used Late/analytic geometry, three trajectories,
seeds 84551/84559/84571 with cap 48,000 and 93701/93719/93731 with
cap 120,000. For 7F5D cognate, dry top-1 hits were 5/6; both fixed-water
arms gave 0/6. All arms gave 0/12 top-1 across the two reciprocal directions.

## Follow-up with two independent test ligands

At one 7F5D site the transferred mean oxygen lay 1.87 Å from the measured
ligand. The closest measured water was 0.91 Å away from the prediction
and 2.69 Å from the ligand. A second donor site had no measured water within
9 Å. All three 7RWO sites lay within 0.15 Å of a measured water. The post hoc
ligand distances diagnosed the pilot failure. Follow-up predictions used
only water observations, protein coordinates and known receiver boxes.

Select donor-supported sites with an observed receiver-chain water (occupancy
≥0.8, protein clearance ≥2.5 Å) within 1.25 Å. Compare *supported* (retain
the donor coordinates) and *registered* (replace with observed coordinates)
with *dry* and *all three polar donor means*. Match one-to-one by increasing
distance and choose ties deterministically. Both 7RWP and 7RWQ passed the
previous preparation rules as additional test ligands. All arms used new
seeds 95701/95713/95717 at cap 120,000. The two original cognates use water
observed with the evaluated ligand for the registered arm. All six cross
directions use a crystal of another ligand for the receptor water positions.

| Follow-up subset | Dry top1 / pool | All polar donor sites | Supported donor means | Registered coordinates |
|---|---:|---:|---:|---:|
| 2 cognates × 3 seeds | 3 / 4 of 6 | 0 / 3 | 2 / 6 | 2 / 3 |
| 6 cross directions × 3 seeds | 0 / 1 of 18 | 0 / 3 | 0 / 4 | 0 / 1 |
| Of those, two added ligands × two receivers | 0 / 1 of 12 | 0 / 2 | 0 / 2 | 0 / 0 |

Filtering the unsupported water already recovers two 7F5D cognate top hits.
Moving the sites to exact receiver coordinates gives no additional top hit.
This factorial comparison does not attribute the initial failure to the
short contact alone. See [same-assay Kd tests](affinity/README.md) and the
[BRDT BD1 replication](../water_transfer/README.md) for separate endpoints.

## Replay

From the repository root use the environment specified in ../pose/README.md
with Biopython installed. Extract ../pose/inputs.tar.xz to a fresh directory
PREPARED. Then run these commands (B is this benchmark directory).

```bash
python "$B/water.py" PREPARED WORK "$VINA"
python "$B/confirm.py" PREPARED WORK "$VINA" --jobs 2
python "$B/check_outputs.py" PREPARED WORK
```

To verify the released measurements without docking, extract raw.tar.xz to
WORK, copy the four measured CSVs and external_manifest.json there, and run
the last command. The independent scalar RMSD reconstruction checks all 628
pilot and 848 follow-up modes, pose hashes and serialized scores. The raw
archive includes the independent donor PDB files and additional test fixtures.
