# Same-assay transfer to human PARP1

Four inhibitors of the same human PARP1 catalytic construct gave four eligible ≥3-fold comparisons. Native Vina, raw Vina, the uniform-rotor control, the contact-rotor annotation and the unit-weight MMFF94s intraligand replacement all gave 4/4 on either receptor and on the average receptor predictions. A simple rule preferring the larger ligand also gave 4/4. This independently measured series reaches the native ceiling and provides no incremental ranking evidence for either correction.

## Cases and design

Ryan et al. measured kinetic Kd by SPR on the same human PARP1 662–1011 protein. We used every PARP1 ligand and the four exact means from Table 1. The two receivers prepared for the preceding pose transfer, 4UND and 7KK4, were fixed before scoring. Each of four compounds was docked to each receptor with seeds 98207 and 98213, cap 60,000, exhaustiveness 4, one CPU and nine saved modes. Searches used native mode, with the production `--contact_rotor_score` output annotation. All five comparisons use exactly the original rank-one pose for each receptor, ligand and seed. The molecular score per receptor is the median across seeds. The secondary receptor average counts each unique experimental pair once. [Assay table, conditions, uncertainties and provenance](SOURCE.md) · [Parent receptor and pose fixtures](../README.md)

| Ligand | Measured Kd, nM | Prepared charge | Heavy atoms | Native | Uniform | Contact | MMFF unit |
|---|---:|---:|---:|---:|---:|---:|---:|
| Talazoparib | 0.6 ± 0.1 | 0 | 28 | −11.427 | −11.460 | −11.623 | −7.936 |
| Olaparib | 0.8 ± 0.1 | 0 | 32 | −12.780 | −13.165 | −13.165 | −8.769 |
| Niraparib | 5.4 ± 1.4 | +1 | 24 | −9.703 | −9.451 | −9.451 | +10.140 |
| Veliparib | 6.2 ± 3.5 | +1 | 18 | −8.580 | −8.281 | −8.281 | +10.888 |

Table scores are means of the two receptor predictions. Each qualifying direction compares one of the two larger neutral ligands with one of the two smaller cationic ligands. Olaparib and talazoparib have measured Kd values less than 1.4-fold apart. Their nominal order is inverted by the receptor-averaged scorers. Niraparib and veliparib differ by less than 1.2-fold. The ≥3-fold analysis excludes both close pairs. Receptors, seeds and overlapping chemical pairs provide matched repeats. They do not increase the number of independently measured compounds. The large positive MMFF scores of both cationic ligands show why correct ordering here cannot establish an accurate free-energy scale.

The first two graph inputs and both receptor inputs are reused byte-for-byte from the pose transfer. Niraparib and veliparib use the stereo-specified CCD graphs of 7KK5/3JD and 7KK6/78P. For both secondary aliphatic ring amines we explicitly used the protonated form in the near-neutral assay. The niraparib label reports pKa 9.95. The 7KK6 legacy PDB HETNAM describes a pyrrolidinium but the current 78P CCD canonical SMILES and formal-charge field describe its neutral amine. The chosen +1 veliparib graph is an explicit chemical assumption. The benzimidazole tautomer follows the fixed CCD graph. Ligand starts for the two new compounds used ETKDGv3 (seed 43103) followed by MMFF, without crystal coordinates. The identity, prepared SMILES and checksums of all four are in the input archive. [Niraparib structure](https://www.rcsb.org/structure/7KK5) · [Veliparib structure](https://www.rcsb.org/structure/7KK6) · [Niraparib prescribing information, chemical description](https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/214876s003s004lbl.pdf)

Independent Python recomputation checked all 144 C++ rotor counts and scores on the printed coordinates. Maximum score difference from printed raw-energy reconstruction was 0.000889 kcal/mol. Fifteen of 16 original top poses classified every eligible rotor as constrained. The four query–receiver–seed runs that overlap the separate baseline pose experiment agree byte-for-byte after removing the new annotation lines. The specialist ran the production MMFF score API on the same pools with eight converged free reference conformers for each compound. Independently reconstructing its algebra on rounded Vina terms reproduced all 16 native-first results within 0.000853 score units. [Contact hypothesis and external controls](../../affinity/README.md) · [MMFF replacement protocol](../../strain/README.md)

## Replay

From the repository root, extract both sets of frozen fixtures. These commands reproduce scores and independent checks without network or new search. The optional MMFF command uses the current binary path supplied by the caller. The CSV supplied with the benchmark permits the remaining checks without compiling the native binary.

```sh
P=benchmarks/frozen_core/transfer
A=$P/parp_affinity
R=$(mktemp -d)
W=$(mktemp -d)
(cd "$P" && sha256sum -c SHA256SUMS)
(cd "$A" && sha256sum -c SHA256SUMS)
tar -xJf "$P/inputs.tar.xz" -C "$R"
tar -xJf "$P/outputs.tar.xz" -C "$R"
tar -xJf "$A/inputs.tar.xz" -C "$W"
tar -xJf "$A/outputs.tar.xz" -C "$W"
python "$A/summarize.py" "$W" "$R" --out "$W"
python "$A/audit.py" "$W" "$R" --tables "$W"
cmp "$W/summary.csv" "$A/summary.csv"
# Optional standalone MMFF recalculation
python "$A/mmff.py" "$W" "$R" /absolute/path/to/vina "$W/mmff_recheck.csv" --jobs 2
# Optional regeneration of all 16 full native pools
python "$A/run.py" "$W" "$R" --binary /absolute/path/to/vina
```
