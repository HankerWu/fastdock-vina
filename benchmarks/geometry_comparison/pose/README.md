# Experimental structures, local refinement and receptor assembly

The human subset of the earlier frozen geometry panel and a supplementary human panel were prepared from experimental PDB depositions. Three original human entries failed the previously declared 8 Å dry, single-chain pocket rule. The closest omitted nonwater atoms were iodide at 6.172 Å in DDR1 5SB2, magnesium at 5.691 Å in PDE10A 5SIS and FMN at 4.181 Å in DHODH 6M2B. Entry 5SD5 is *M. tuberculosis* DHFR. The original multispecies panel was not completed here. The supplement and four potential reciprocal structures were specified before examining their docking predictions. Twenty human candidates yielded ten prepared structures across eight targets and ten distinct Bemis–Murcko scaffolds. All input IDs, preparation exclusions, deposited-source hashes, prepared-input hashes and ligand identities are in `../pose_cases.csv`, `manifest.json` and `ligands.csv`. The supplementary panel is a convenience sample.

## Protocol and observed performance

Each ligand started from the CCD chemical graph, ETKDG and MMFF, in a 24 Å cubic box centered on the receiving structure's measured ligand. The protocol assumes the binding site is known. PDB references were used only for box placement, structure eligibility, receptor alignment and evaluation. Known crystal ligand positions were never provided as a starting pose. The receptor was rigid, dry and single-chain for the primary comparison. The two reciprocal pairs are BPTF 7F5D/7RWO and BCL6 7LZQ/7QK0. Each pair contains different scaffolds. The local Cα alignments used 39 residues (0.627 Å) for BPTF and 27 residues (0.273 Å) for BCL6. The BPTF constructs have different author numbering. A unique whole-protein sequence register determined the ±2733 offset before the standard local Cα fit. No ligand superposition was used.

All methods used one CPU thread per run, identical seeds 63419 and 63421, three search trajectories, up to nine reported modes, a Vina energy window of 100 and Late explicit refinement. The comparator uses strict Armijo with ordinary scalar initialization. Both geometry variants use the same geometric metric with different derivative implementations. The lower cap was 12,000 evaluations per trajectory. After the low coverage was observed, the full comparison of the two production arms was repeated with a tenfold higher cap. Success means receptor-frame symmetry- and stereochemistry-aware heavy-atom RMSD ≤2 Å. The pool metric is the minimum over *emitted* modes. A seed is a repeated prediction on the same complex.

| Budget and optimizer | Cognate top / pool out of 20 | Cross top / pool out of 8 | Median top RMSD, 28 cells (Å) | Explicit evaluations, 28 cells |
|---|---:|---:|---:|---:|
| 12,000 Armijo | 0 / 1 | 0 / 0 | 7.275 | 6,542 |
| 12,000 finite-difference geometry | 0 / 1 | 0 / 0 | 7.158 | 8,419 |
| 12,000 analytic geometry | 0 / 1 | 0 / 0 | 7.169 | 8,288 |
| 120,000 Armijo | 1 / 4 | 0 / 0 | 7.128 | 6,425 |
| 120,000 analytic geometry | 1 / 4 | 0 / 0 | 7.097 | 8,400 |

Global search evaluation counts matched across arms within each cell (1,011,441 across 28 cells at 12,000 and 10,083,216 at 120,000). The explicit optimization paths diverged slightly between the two metric implementations. At 12,000, 4 of 28 output PDBQT files were byte-identical, the largest absolute difference in their top-pose RMSD was 0.0564 Å or less, and every top/pool threshold outcome agreed. Wall time in the CSVs includes receptor loading, grid construction, search, refinement and serialization. These jobs ran concurrently on a shared host, so their elapsed times cannot establish a speedup.

Increasing the cap recovered a BPTF 7F5D top pose at 0.510 Å with Armijo and 0.474 Å with analytic geometry; both recovered BRDT 7MRH near-native candidates twice but ranked them below first. For an independent check of cognate pose coordinates and atom correspondence, we reconstructed the molecules from output PDBQT with Meeko and used RDKit CalcRMS without superposition against PoseBusters SDFs. Graph identity agreed for all ten. Across 529 low-budget and 360 high-budget cognate modes the maximum absolute discrepancy was below 0.030 Å and all 889 success classifications agreed. The source archive MD5 was `f004ac7c4e68317b5348497d2bb6bee6`. [Authors' archive](https://zenodo.org/records/8278563) · [PoseBusters paper](https://doi.org/10.1039/D3SC04185A)

## Crystal packing and deposited assemblies

`packing.csv` enumerates every nonidentity-image protein heavy-atom site ≤4 Å from the selected experimental ligand, including unit-cell translations. Occupied alternative coordinates at one protein atom site count once per image. The Gemmi audit compared each space-group transform in Cartesian coordinates with deposited PDB biological assembly operators *relative to the receptor chain copy*. Exact translation matters. The selected ligand conformer uses all shared atoms and the first sorted nonblank label, including the D-only ligand in PHIP 5S8I.

| Structure | Closest mate, Å | Sites ≤4 Å | In deposited assembly among these sites |
|---|---:|---:|---:|
| PHIP 5S8I | 4.597 | 0 | 0 |
| BPTF 7F5D / 7RWO | 6.247 / 4.470 | 0 / 0 | 0 / 0 |
| HSP90AA1 7LT0 | 6.132 | 0 | 0 |
| BCL2L1 7JGW | 3.466 | 7 | 0 |
| BAZ2A 7BLA | 3.215 | 9 | 0 |
| BRDT 7MRH | 3.516 | 3 | 0 |
| BCL6 7LZQ | 2.939 | 17 | 10 |
| BCL6 7QK0 | 3.391 | 16 | 16 |
| USP5 7MS7 | 3.406 | 7 | 0 |

The six structures with direct contacts remained in the primary count. Restricting to directions in which neither the reference nor the receiving ligand contacted any packing mate within 4 Å leaves 12 paired cells. Armijo and analytic geometry both yielded 0 top and 1 pool hit at 12,000, or 1 top and 2 pool hits at 120,000.

The authors designated a dimer as biological assembly 1 in the REMARK 350 records of both BCL6 structures. In 7LZQ the closest image (2.939 Å) is a different lattice translation from its deposited dimer partner (3.50 Å). We rebuilt each dimer on the original coordinates with the repository's assembly tool and repeated both cognate and both reciprocal combinations under Armijo at 120,000. Four seeds included two new seeds 91571 and 91583. The monomer and the deposited dimer each produced 0/16 top and 0/16 pool hits. In 7LZQ cognate docking all four minimum pool RMSDs decreased from 6.42, 6.82, 6.55 and 6.40 Å to 4.99, 3.76, 3.84 and 5.89 Å. In 7QK0 cognate docking all four increased. Adding the separately identified packing image to the 7LZQ dimer produced one top/pool hit in eight trials with this receptor (1.82 Å for its own ligand, seed 91571). Both the image choice and pocket center use crystal information.

A separate four-seed repeat of the entire BCL6 matrix (127003, 127019, 127043, 127057) yielded no hits in eight dimer-plus-lattice trials. The author dimer yielded one near-native 7LZQ cognate candidate at 1.02 Å (seed 127043), while its top pose was 7.21 Å. Both four-seed series combined yielded 0/32 top and 0/32 pool hits for monomers, 0/32 top and 1/32 pool hits for author dimers, and 1/16 top and 1/16 pool hits for the 7LZQ-receiver dimer with its additional image. The dimer lowered the cognate 7LZQ oracle RMSD on all eight seeds, with no top-1 successes. The repeated reciprocal measurements all involve the same target and the same two ligands.

As a control, all four other contact-bearing structures with exclusively nonassembly mates were redocked with or without those copies, using the two new seeds. Adding the crystal copies changed top hits from 3/8 to 2/8 and pool hits from 4/8 to 2/8. In BCL2L1 the monomer produced 1.25 and 1.27 Å top poses versus 10.96 and 12.47 Å with three nearby crystal images. In BRDT top poses changed from 3.01 and 1.32 Å to 1.47 and 1.49 Å. BAZ2A and USP5 had no hits in either arm. These manipulations test the crystallographic condition of each individual structure; they do not estimate affinities or justify automatic inclusion of all crystal neighbors.

## Replay

From the repository root use the preparation environment documented for `benchmarks/paired_crossdocking` (Python 3.12, RDKit 2026.09.2, Meeko 0.8.0, PDBFixer 1.12.0, OpenMM 8.6.1, Gemmi 0.7.5, SciPy 1.17.0). Archives contain frozen fixtures and all reported output poses and logs. CSVs give the ranked modes, one row per completed cell, and the independent cognate audit. Use extracted inputs to avoid variability during receptor repairs. Or rerun `pose_prepare.py WORK --jobs 2` against RCSB to recreate the inputs. SHA-256 hashes are in `SHA256SUMS`.

```sh
P=benchmarks/geometry_comparison
mkdir -p work/inputs
tar -xJf "$P/pose/inputs.tar.xz" -C work/inputs
python "$P/pose_run.py" work/inputs /path/to/analytic/vina work/low --finite-binary /path/to/142c676/vina --jobs 2
python "$P/pose_run.py" work/inputs /path/to/analytic/vina work/high --max-evals 120000 --jobs 2
python "$P/pose_packing.py" work/inputs work/packing.csv
python "$P/pose_assembly.py" work/inputs /path/to/analytic/vina work/assembly --jobs 2
python "$P/pose_assembly.py" work/inputs /path/to/analytic/vina work/assembly-repeat --seeds 127003 127019 127043 127057 --jobs 2
python "$P/pose_lattice_control.py" work/inputs /path/to/analytic/vina work/lattice --jobs 2
python "$P/pose_audit.py" /path/to/posebusters_benchmark_set work/low work/sdf-audit.csv
(cd "$P/pose" && sha256sum -c SHA256SUMS)
```
