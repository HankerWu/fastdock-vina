# Experimental Kd source

Ryan et al., *J. Biol. Chem.* 2021, 296, 100251, Table 1, PARP1 kinetic Kd column. The ligand binding experiment used biotinylated human PARP1 residues 662–1011 on a streptavidin surface with a Biacore S200. The authors used referenced sensorgrams and a 1:1 Langmuir model. Table 1 reports 25 °C, 25 mM HEPES, 150 mM NaCl, 5% glycerol, 0.5 mM TCEP, 2% DMSO, 0.02% Tween 20, pH 7.4. The separate chip functionalization step describes 50 mM HEPES at pH 7.2 and 10 °C. Both receivers have the same human PARP1 catalytic-domain residue range in their PDB DBREF records. [Article](https://doi.org/10.1074/jbc.RA120.016573) · [Author manuscript and tables](https://pmc.ncbi.nlm.nih.gov/articles/PMC7948648/) · [Raw NCBI BioC XML](https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC7948648/unicode)

| Compound | PARP1 kinetic Kd, nM | Reported ±, nM | Stronger group for the ≥3-fold comparison |
|---|---:|---:|---|
| Talazoparib | 0.6 | 0.1 | Yes |
| Olaparib | 0.8 | 0.1 | Yes |
| Niraparib | 5.4 | 1.4 | No |
| Veliparib | 6.2 | 3.5 | No |

Each tabulated value averages two or more independent measurements. The table footnotes do not specify the error statistic. Eligibility uses the reported mean Kd, in keeping with the original ≥3-fold rule. The two within-group comparisons fail that threshold. The published article reports both PARP1 and TNKS1 measurements in the same table. The code selects its PARP1 column by position and checks the four values against an independent reading. The raw `table1.xml` was extracted from the NIH BioC passage (id `tbl1`, type `table`). Complete BioC response SHA-256 was `fb53e2aa9df61c1437f010c1beb550eea1f9a5e168dc0acd5ff794824cd11ec3`. An independently fetched copy of the journal PDF had SHA-256 `d622541f002154d6af8f67404165fce5801f3931dbac99a7d81f74d713134573`. [Journal PDF copy](https://pdfs.semanticscholar.org/6b1d/448f01506d1a0dcadf40c795e01a0c1fc936.pdf)
