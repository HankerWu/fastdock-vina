# Independent transfer of BRDT BD1 pocket waters

Two related human BRDT BD1 ligands provide a second bromodomain test of recurrent water selection. The two ligand structures have reported affinities on the same isolated human protein in the same assay. Receivers were selected independently of hydration outcomes. BRDT 7MRH (ZMJ, NC-II-259-1) contributed one case in the initial human geometry panel. We added 7MRG (ZMM, NC-II-153-1), a monovalent analog from the same paper. Chain A passes the established dry receptor preparation in each structure. The closest omitted nonwater atom is 9.90 Å for 7MRH and 9.73 Å for 7MRG. The two structures are 1.98 and 1.99 Å resolution respectively. [7MRH](https://www.rcsb.org/structure/7MRH) · [7MRG](https://www.rcsb.org/structure/7MRG) · [Primary measurements](https://doi.org/10.1021/acs.jmedchem.2c00453)

Water donors were two separate human BRDT BD1 crystal structures, 4FLP with JQ1 and 9YCQ with a substituted diaminopyrimidine, at 2.23 and 1.40 Å resolution. An independent donor means a different PDB entry from either receiver or the other donor. Residue identity, sequence coverage and local backbone fit were checked before water transfer. The chemical anchor was pocket Asn109 ND2 on chain A. Selected oxygen sites had occupancy ≥0.8 in both donors, donor–donor separation ≤1.0 Å, averaged position within 8 Å of the anchor and distance ≥2.5 Å to any prepared protein heavy atom. The donor site search and fitting used protein coordinates. The supported arms also use the identified receiver's experimental water sites, including water coordinates measured in the query's bound state for cognate redocking. Each matched receiver oxygen needed occupancy ≥0.8, distance ≤1.25 Å from its donor mean and ≥2.5 Å from prepared protein heavy atoms. All occupied sites were represented as fixed donor/acceptor oxygens. The polar hydrogens set the Vina atom typing. [JQ1 crystal](https://www.rcsb.org/structure/4FLP) · [2026 donor crystal](https://www.rcsb.org/structure/9YCQ) · [2026 BRDT BD1 study](https://doi.org/10.1021/acs.jmedchem.6c00180)

## Paired pose results

One recurrent site on 7MRH and two on 7MRG qualified. All three also matched measured receiver waters. Their donor-mean versus receiver displacements were 0.558 Å (7MRH), 0.386 and 0.266 Å (7MRG). The all-donor and supported-donor arms generated identical receptor files and byte-identical output pools in every cell. This particular pair isolates positional relocation but provides no variation in filtering. Both cognate structures and both reciprocal directions were docked with seeds 95701, 95713 and 95717, three trajectories each capped at 120,000 evaluations, Late/analytic geometry refinement, 24 Å known-site box and nine requested modes. All arms used the same prepared ligands and search settings. Each cell's reference was its observed ligand geometry expressed in the receiver frame, with chemical symmetry and chirality preserved. The 48 runs emitted 428 poses.

| Water treatment | Cognate top 1 / top 5 / any within 2 Å, out of 6 | Reciprocal top 1 / top 5 / any within 2 Å, out of 6 |
|---|---:|---:|
| Dry | 0 / 2 / 2 | 1 / 3 / 4 |
| All donor means | 2 / 3 / 3 | 5 / 5 / 5 |
| Receiver-supported donor means | 2 / 3 / 3 | 5 / 5 / 5 |
| Receiver-supported oxygen at its measured position | 4 / 4 / 4 | 3 / 3 / 3 |

Independent recurrent waters improved top-1 selection at this pocket in both test types. Placing them on the measured receiver oxygens added two cognate successes but lost two reciprocal successes relative to donor averages. There are two ligands, two receivers and repeated seeds. Receiver-specific environments, small differences between structures and sampling all enter this comparison.

## Matched Kd and contact geometry

Guan et al. Table 1 reports qPCR-based DiscoverX BromoScan Kd from a single assay in duplicate (mean ± SD). The 7MRH sulfonamide measured 77 ± 23 nM, versus 875 ± 177 nM for the 7MRG pyrrolidinone, an 11.36-fold difference. In each receiver separately we took the median seed score for each ligand. All four protocols ranked this single chemical pair in the opposite order in both receivers. For example, the registered-water medians in 7MRH were −8.052 (sulfonamide) and −9.131 (pyrrolidinone), and in 7MRG were −8.063 and −9.643. These are repeat judgments of one pair. [Table 1 and structural discussion](https://doi.org/10.1021/acs.jmedchem.2c00453) · [Original NCBI full text XML](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=11727429)

The authors associate the sulfonamide with an additional interaction with Lys60. We measured the shorter distance from Lys60 NZ to ZMJ O02 or O03 in the deposited 7MRH structure as 4.113 Å in chain A and 2.757 Å in chain B. This experiment uses chain A. The three registered-water top poses in 7MRH had Lys NZ–sulfonyl O minima of 6.50, 5.96 and 6.99 Å even though two had whole-ligand RMSD below 2 Å. Oxygen serials 4 and 5 in the prepared molecule are the sulfonyl oxygens, both attached to sulfur serial 3. This contact diagnostic distinguishes functional-group placement from the aggregate RMSD. The two distances and the assay alone cannot assign an energetic contribution to the side chain.

## Offline audit and replay

The archives contain raw coordinate depositions and CCD records for both receivers, two donor structures, prepared inputs, all four modified receptor sets, complete emitted poses and logs. The offline audit checks each raw output hash against its cell record, parses scores and coordinates from PDBQT, recomputes protein-frame transformations using a quaternion fit, and independently recomputes every emitted RMSD and the chemical pair ordering. It uses the same stored chemical atom mappings. Auditing requires NumPy and SciPy. Regeneration also needs the repository's structural preparation dependencies and Biopython.

```bash
D=benchmarks/geometry_comparison/water_transfer
mkdir -p brdt
tar -xJf "$D/inputs.tar.xz" -C brdt
tar -xJf "$D/outputs.tar.xz" -C brdt
cp "$D/cells.csv" "$D/modes.csv" brdt/
python "$D/audit.py" brdt
python "$D/brdt.py" brdt /path/to/analytic/vina --jobs 4
```

An additional [crystal-chain diagnostic](chain_copies/README.md) tested the
shorter 7MRH-B Lys60 contact with the same ligand inputs. One rank-two pose
recovered a 3.01 Å contact at 1.13 Å ligand RMSD. Both ligand medians on
this receptor still predicted the wrong Kd order. The other B receptor failed
the same pre-existing preparation checks and its failure log is included.
