# 独立单斜晶型中，VH032 的实测姿势改变了多少？

人 VHL 的 9GIO 晶体含有与 4W9H 完全相同的 VH032（CCD 3JF）。空间群由 P 41 2 2 换成 C 1 2 1，但四个 4W9H 配体拷贝对独立沉积的 9GIO 配体只有 0.153–0.195 Å 的全配体重原子 RMSD。在先前指定的 55 个口袋 Cα 配准下，历史盲搜里最好的 VH032 候选用新参照评价仍有 3.326 Å 误差。[事前方案](PROTOCOL.md) · [9GIO 实验结构](https://www.rcsb.org/structure/9GIO) · [Lucas 等原论文](https://doi.org/10.1021/acsmedchemlett.4c00582)

从原 mmCIF 的重原子和 CCD 化学键重建立体化学，InChIKey 为 GFVIEZBZIUKYOG-SVFBPWRDSA-N，33 个重原子均只有一个位形、占有率 1。全配体有 12 个保手性等价对应，在同一受体坐标系取最小 RMSD。9GIO 分辨率 1.486 Å，55 个口袋 Cα 全部存在且无备用位形，映射到原接收受体 5LLI 的 Cα RMSD 为 0.305 Å。分别将 9GIO 和 4W9H 四拷贝映射到 5LLI 后，4W9H 的 C/F/I/L 与 9GIO 相差 0.189/0.153/0.195/0.183 Å。直接把 9GIO 的蛋白配准到 4W9H 链 C 得到 0.189 Å。9GIO 还含与 Cys77 共价结合的 A1IMD，其最近重原子距 VH032 为 12.007 Å；这一观测同时改变晶格、共配体及实验条件。两个晶型仍不足以分离各因素贡献。[独立 mmCIF](9gio.cif.gz) · [化学与配准数据](results/comparison.json)

| 4W9H 查询 → 5LLI 盲搜的 VH032 坐标 | 对 4W9H 链 C RMSD，Å | 对 9GIO RMSD，Å | 对 9GIO 晶格蛋白的最近距离，Å |
|---|---:|---:|---:|
| 当前 Pro 首姿，种子 27183 | 10.203 | 10.234 | 3.432 |
| 当前 Pro 首姿，种子 91439 | 11.128 | 11.123 | 6.110 |
| 全部 43 个原池模式中最近的一个 | 3.319 | 3.326 | 3.482 |
| 旧参照引导的局优见证 | 0.761 | 0.779 | 4.659 |

八个原盲搜池合计 43 个模式，全部模式换参照的 RMSD 变化绝对值至多 0.067 Å。种子 91439 的父版与当前版首姿按新参照为 8.406 和 11.123 Å。完整键表 VH032 对照的 16 个模式另列，在该配体上与对应原结果重复。[所有模式、最近原子与晶格操作](results/forms.csv) · [逐池 Top-1 与 Oracle](results/pools.csv)

9GIO 实测配体到另一蛋白分子的最近重原子为 4.272 Å，来自相邻拷贝的 Val181 CG1。4W9H 的链 C 沉积配体最近蛋白为邻居 Arg182 NH1，距离 3.306 Å。43 个 VH032 盲搜模式中，16 个在 9GIO 晶格下出现 <2.3 Å 的接触。种子 91439 当前首姿在 5LLI/4W9H 两套四方晶格都产生明显穿插（0.318/0.235 Å），转换到 9GIO 后，最近晶格蛋白却在 6.110 Å 外，而相对于独立测定的配体姿势仍错 11.123 Å。因此，穿插是晶格特定的现象；即使固定同一配体，也无法用其解释多个晶型共同支持的全部姿势误差。[两份四方晶格的完整结果](../README.md)

RCSB 的 CCD 配体条目把 VH032（3JF）和 VH298（6Z3）分别列在三个和一个含非共价独立配体的 PDB 条目中；本轮读取的是额外的 VH032 单斜晶型。对 VH298 的 5LLI 结果尚无此类同配体异晶型验证。[3JF 条目](https://www.rcsb.org/ligand/3JF) · [6Z3 条目](https://www.rcsb.org/ligand/6Z3)

按上级说明解包工作目录后运行。Gemmi 0.7.5、NumPy 2.3.5、SciPy 1.17.0 和 RDKit 用于复算。在单线程下，从已解包输入到写表约 0.83 秒墙钟、0.83 秒 CPU，独立核验约 0.57 秒。原始 9GIO 压缩包来自 RCSB 并与 wwPDB 数据归档的内容相同，解压后 SHA-256 在 comparison.json 中。RMSD 数值由第二条路径根据 4W9H 原坐标恢复 CCD 原子名映射，以 Gemmi 重新配准，而晶格最近距离再用 Gemmi 周期邻居与直接距离计算。复核 61 行，最大 RMSD 差 7.11×10⁻¹⁵ Å，最大距离差 2.94×10⁻¹⁴ Å，冲突对数一致。[坐标下载](https://files.rcsb.org/download/9GIO.cif.gz) · [独立核对](results/verification.json)

```bash
D=benchmarks/accepted_bfgs/curvature/transfer/packing/independent_form
OPENBLAS_NUM_THREADS=1 python "$D/compare.py" "$W" "$W/9gio"
OPENBLAS_NUM_THREADS=1 python "$D/check.py" "$W" "$W/9gio"
cmp "$D/results/forms.csv" "$W/9gio/forms.csv"
cmp "$D/results/pools.csv" "$W/9gio/pools.csv"
```
