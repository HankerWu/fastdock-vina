# BIRB-796 在 p38α 两种受体状态下的结构核查

主实验中 BIRB-796 对一个种子的同靶排序影响很大，因此追加这项有明确选例偏倚的核查。其在 1KV2 中的晶体姿势转移到 1A9U 后，配体与蛋白发生严重重叠；在与晶体配体匹配的 1KV2 中，本轮两枚种子、三种回退方法也都未由搜索返回晶体姿势。

**结构证据。** 1KV2 是人 p38α 与 BIRB-796 的 2.8 Å 结构。1A9U 是与 SB203580 的 2.5 Å 结构。两者分别提供 DFG-out 和 DFG-in 状态。从两个口袋的 50 对对应 Cα 求变换（残差 RMSD 1.62 Å），将 BIRB-796 的晶体坐标转入 1A9U。39 个重原子都在原 26 Å 盒内，但有 33 对蛋白–配体重原子距离小于 0.75 倍范德华半径之和，其中 22 对涉及 Phe169，10 对涉及 Met109，另有 1 对涉及 Leu108。最近一对仅 0.39 Å。排除活化环后用其余 45 个口袋 Cα 对齐，仍有 34 对冲突，其中 Phe169 占 23 对；以 325 个全链公共 Cα 对齐，得到 50 对，其中 30 对涉及 Phe169。1KV2 中的原位晶体参考零冲突，最近 2.89 Å。以晶体坐标显式局部松弛一次，转入 1A9U 的构象移动 2.09 Å，原位 1KV2 仅 0.32 Å。两个受体分别使用已定义的 26 Å 盒，配体的原子参数和可旋转键相同。[Pargellis 等，2002](https://doi.org/10.1038/nsb770) · [1KV2](https://www.rcsb.org/structure/1KV2) · [1A9U](https://www.rcsb.org/structure/1A9U)

**真实对接。** 下表每格依次为首姿、返回池最佳的对称性校正重原子 RMSD（Å）。四个受体–种子组合共产生 12 个臂运行、69 个模式，其中六次 DFG-in 运行复用主实验的存档输出。12000 是命令行上限，实际调用数见逐运行表。2 Å 阈值下，所有单元的 Top-1、Top-5 及返回池 oracle 均为零。

| 受体 | 种子 | 当前 FastDock | 受保护二分 | 受保护二次 |
|---|---:|---:|---:|---:|
| 1A9U | 458921 | 14.85 / 14.26 | 11.38 / 10.69 | 11.51 / 9.81 |
| 1A9U | 703109 | 11.39 / 11.39 | 9.62 / 8.28 | 9.30 / 9.30 |
| 1KV2 | 458921 | 3.91 / 3.91 | 11.75 / 11.75 | 12.42 / 11.66 |
| 1KV2 | 703109 | 9.90 / 5.20 | 7.99 / 6.04 | 16.24 / 8.36 |

1A9U 上种子 458921 的原生回退与二次回退给出 −1.813 和 −9.424 的首姿分数，首姿仍分别偏离参考 14.85 和 11.51 Å。在 1KV2 的同一随机种子下，二次回退比当前 FastDock 的首姿偏差多 8.51 Å。实验向后判断的是已知晶体结合模式的回收；受体状态变更还伴随口袋几何、结晶条件、重建位置及盒中心变化，这些结果不能量化 DFG 状态对全部 Kd 排序的因果贡献，也不能从单个配体推断哪一套受体的分数整体更好。[逐运行结果](results/cells.csv) · [各模式、碰撞及最短距离](results/modes.csv)

**化学身份与结构局限。** BIRB-796 输入（CID 156422）与 B96 全部 39 个重原子的图和立体化学一致，共 24 个合法映射。SB-203580 的建模输入没有指定硫手性；实际对接 PDBQT 的三维构象与 1A9U 的 SB2 同为亚砜，硫手性相反。忽略手性有 8 个匹配，要求匹配则为零。两者的三维 InChIKey 分别为 CDMGBJANTYXAIV-HHHXNRCGSA-N 和 CDMGBJANTYXAIV-MHZLTWQESA-N。未据此推断实测材料是否为外消旋混合物或某一对映体，故省略了 SB-203580 晶体 RMSD。1KV2 的 115–122 和 170–184 区间无坐标；后一区间紧邻有晶体坐标的 Phe169。两种盒均省略水。1KV2 的 Tyr35 有 7 个重建侧链重原子距 BIRB-796 小于 8 Å，最近 5.29 Å，故 1KV2 作为探索性有缺失的模型单列，未并入预定结构面板。[身份检查](results/identity.json) · [几何检查](results/geometry.json)

**复现。** 从仓库根目录执行以下命令，使用与主实验相同的 Python 和 C++ 环境。P38 核查直接复用三个已编译实验臂。未生成酸碱状态或对映体备选，也未合并受体分数。已独立重建 SVD 对齐并以 RDKit CalcRMS 核对 69 个模式。

```sh
D=benchmarks/search_backtracking/p38_state
B=/tmp/fastdock-search
W=/tmp/mapk14-state-check
mkdir -p "$W"
tar -xJf "$D/evidence.tar.xz" -C "$W"
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
python "$D/experiment.py" "$W" "$B/native/vina" "$B/quadratic/vina" "$B/halving/vina" --cpu 6
python "$D/verify_modes.py" "$W"
```

从原始 PDB 重建新增受体时，可在仓库根目录执行 `PYTHONPATH=benchmarks/kinase_state_affinity/source python -c 'from pathlib import Path; from prepare import receptor; import sys; receptor(Path(sys.argv[1]), "1KV2", "B96")' "$W"`。已用独立制备逐字节复核 1KV2 PDBQT、盒与补建清单。原始 PDB、CCD、查询、两个受体、存档对照、预测和诊断坐标均在 [evidence.tar.xz](evidence.tar.xz)。[来源及 SHA-256](sources.json) · [独立复核](results/verification.json)
