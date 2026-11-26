# Does the shorter crystallographic Lys60 contact recover the measured order?

This post-hoc diagnostic tests the chain asymmetry of human BRDT BD1. In 7MRH the closest Lys60 NZ–sulfonyl O distance is 4.113 Å in chain A and 2.757 Å in chain B. The prepared monomer changes to the other complete chain, including every resolved residue. The experiment measures the effect of that receptor copy. With two ligands and repeated seeds, there is one independent experimental comparison. [Crystal structure](https://www.rcsb.org/structure/7MRH) · [Guan et al., Table 1](https://doi.org/10.1021/acs.jmedchem.2c00453)

All ligand-bearing copies in 7MRG and 7MRH are A or B. We reused the 12 dry A-chain dockings and ran six dry B-chain dockings (54 newly emitted poses). Both molecules used identical prepared ligand files across chains. The earlier dry Late/analytic geometry protocol used a 24 Å box, three trajectories capped at 120,000 evaluations, nine output modes and seeds 95701, 95713, 95717. 7MRH-B passed, with the nearest omitted nonwater atom at 10.06 Å and the nearest repaired heavy atom at 29.66 Å. 7MRG-B passed the 8 Å exclusion distance (nearest omitted ethylene glycol at 8.36 Å). Meeko could not assign templates to terminal Ala27 and Glu136, reporting excess inter-residue bonds. We excluded that receptor and retained its valid ligand, B-chain experimental pose and pocket backbone for cross-docking into 7MRH-B. The input archive includes the original dry-preparation diagnostic.

| Receptor copy | Sulfonamide, 77 nM, median score | Pyrrolidinone, 875 nM, median score | Predicted strong minus weak |
|---|---:|---:|---:|
| 7MRH-A | −7.766 | −9.421 | +1.655 |
| 7MRG-A | −8.210 | −9.767 | +1.557 |
| 7MRH-B | −7.682 | −9.684 | +2.002 |

Every eligible receptor favors the weaker ligand. In 7MRH-B, the sulfonamide rank-1 poses across the three seeds have whole-ligand RMSDs 6.59, 3.59, 3.00 Å and Lys60 distances 10.14, 3.91, 3.18 Å. The corresponding nearest-to-crystal poses rank fourth, fifth and second, with RMSDs 4.03, 3.21, 1.13 Å and Lys60 distances 3.14, 4.11, 3.01 Å. Thus one sampled near-crystal geometry reconstitutes the short contact in the more favorable receptor copy, while its real Kd order remains reversed under native Vina scoring. The measured affinity difference can include interaction energetics, state populations, solvation and effects elsewhere in the molecules. [All poses and contact minima](contacts.csv) · [Per-receiver pair estimates](orders.csv)

Reproduce after expanding the parent water-transfer archives and copying its cells.csv and modes.csv into the original directory. The copy-preparation wrapper interchanges A and B identifiers before running the unchanged preparation criteria. The generated files label the selected original B copy as A. The original depositions and their hashes establish the source identity.

```bash
D=benchmarks/geometry_comparison/water_transfer/chain_copies
mkdir copy-evidence
tar -xJf "$D/chain_B.tar.xz" -C copy-evidence
cp "$D"/{cells,modes,contacts,orders}.csv copy-evidence/
python "$D/verify.py" brdt copy-evidence
python "$D/copies.py" brdt fresh-copies /path/to/vina --jobs 4
```

The independent audit validates all 18 raw pools, 160 model scores, scalar whole-ligand RMSDs, top-1/closest-pose contact distances, and median score differences from the shipped fixtures.
