# Receptor hydroxyl orientations on four new human target structures

Frozen on `dec913f584a2bf188ce9704b0602f11c56069627` before any docking or outcome inspection. Use the chain-A cognate complexes below, with no substitution. These identifiers and gene names have no standalone matches in the parent checkout's tracked documentation, scripts, and result tables. This is an independent target panel for the receptor hydroxyl-rotation model. Holo receptor and the experimental ligand centroid define the 24 Å cubic pocket. Protein structure and ligand chemistry are known. Query torsions and starting coordinates come from canonical CCD SMILES, ETKDGv3 (43103), and MMFF as in `benchmarks/pair_gradient/pose/prepare.py`.

|Target|PDB|Ligand|
|---|---|---|
|TNKS2|4BFP|SWY|
|MALT1|4I1R|LZU|
|NR3C2|2AA2|AS4|
|MAP3K7|2YIY|YIY|

Reuse the preparation exclusions for incomplete/inconsistent ligand or stereochemistry, failed receptor preparation, and omitted nonwater or repaired heavy atoms within 8 Å. Report all exclusions. Inspect waters, other chains, metals and covalent links, structural resolution, and whether the symmetry/alternate-location mapping is complete. No new protonation adjustment. Use the unmodified parent Pro CLI, one CPU, exhaustiveness 4, 9 modes, 100 kcal/mol output range, seeds 43103 and 43117, budgets 3000 and 12000 evaluations per trajectory. Save the output before scoring. The same 3-decimal coordinates and native scores are used by every reranking arm. Native upstream at the same settings is a secondary comparator, using its own candidate pool.

Use the independently specified model in `../hydroxyl_partition` unmodified. Compare native printed order, fixed hydrogens, energy-minimized angular orientation and orientational partition correction. The primary comparison isolates the change from fixed orientation to the normalized orientational partition. Also report the net change from the actual native order. Preserve original rank on ties. Report symmetry-corrected heavy-atom RMSD in the protein frame for every output mode, Top-1 / Top-5 / pool-oracle success at RMSD ≤2 Å and Top-1 conditional on a successful pool. Report per target/budget/seed; repeats are not independent targets. Count rank changes even when no threshold crosses. No fitting, gating or register changes after pool unblinding. No affinity or unrestricted blind-docking inference.
