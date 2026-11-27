# 受体末端酰胺方向的固定候选池迁移实验

使用分支提交 `fab52f715811133301dce9479c90890fa56965e0` 重新对接 46 次，产生 229 个输出姿势。开发实验所用的 HSP90α / ERRγ 均在本面板以外。历史迁移面板覆盖 DPP4、NR1H2、BRD9 和 THRB，另从 PDB 下载人 ACHE、BRD4、CREBBP、BAZ2B、WDR5 的 7 个结构。候选池已经固定，评分前登记的端点、两个种子和受体准入条件见 [PROTOCOL.md](PROTOCOL.md)。耦合项作为事后敏感性分析加入，未按迁移结果改动参数。共晶配体给出统一对接框位置，所有被测配体从化学图独立生成初始坐标。

| 在同一池中比较原生评分与双态积分 | 池数 / 姿势数 | 首姿平均 / 中位 Å | Q1 / Q3 Å | P90 / P95 Å | 前五最小 / 全部输出最小 / 排序损失的平均 Å |
|---|---:|---:|---:|---:|---:|
| 历史材料作为迁移测试 | 32 / 154 | 6.865 / 7.815 | 3.793 / 9.526 | 10.741 / 10.957 | 5.925 / 5.912 / 0.953 |
| 新下载结构，受体通过准备检查 | 14 / 75 | 7.427 / 6.348 | 4.566 / 10.468 | 12.738 / 15.367 | 6.335 / 6.335 / 1.091 |

双态积分在所有 46 个池中都保留了原生首姿，表中每个首姿配对差均为零。分位数采用经验分布及跳点均值，对应等权样本的第 2 型样本分位数。新结构面板存在不同数量的接收受体和查询，表中按运行等权。三个可评靶点各自的首姿均值为 ACHE 9.354 Å、BRD4 4.673 Å、CREBBP 5.225 Å，靶点等权均值 6.418 Å。每个结构对的两枚种子共享相同受体和配体，在该面板上对配对差做聚类重采样仍将得到全零，不能据此给出新家族总体误差界。

201/229 个姿势含至少一个 8 Å 内的可变末端酰胺。历史和新面板输出的归一化主模型修正分别在 −0.119 至 +0.015、−0.081 至 +0.173 个打分单位之间。新面板唯一的深层名次互换出现在 4NYX 第二枚种子的第 4 / 5 位。二者修正后只差约 0.0004，低于存档原生分数的 0.001 输出分辨率，无法确认该互换能在全精度引擎中复现。所有备选模型（全翻转、逐位点最小、均匀先验、显式耦合）的首姿也没有变化。历史 32 池中只有 6 池，新 14 池中只有 1 池含 RMSD ≤ 2 Å 的输出；输出最低 RMSD 的均值仍超过 5 Å。仅保留了导出候选，可追溯到局部优化、归档及裁剪之后，无法由此辨别搜索是否曾短暂访问过未导出的好构象。

**覆盖与结构解释。** 历史 8 个受体均通过准备，完成全部 16 个查询–受体组合。新面板 7 个受体中有 4 个达到原先规定的非水原子遗漏及近口袋修复标准，11 个预定组合中 7 个可对接，14/14 次返回合法坐标。3SMR WDR5 的未知位点距配体 3.939 Å，4RVR BAZ2B 的乙二醇距配体 7.878 Å，4J0S BRD4 的乙二醇距配体 4.081 Å，均触发既有的受体准备排除条件。4J0S 的配体仍用于 4QB3 的交叉对接。4QB3 和 4NYX 分别有距配体 3.441 和 4.419 Å 的晶体对称伴侣原子，未加入预测受体。4EY6 和 4EY7 的原始 PDB 均报告 pH 7.0；本轮采用 CCD 化学图对应的形式电荷 0，未枚举加兰他敏或多奈哌齐的胺质子化微态。两者较大的误差可能受到微态选择影响，本实验无法给出其贡献大小。原始坐标与逐结构审查保存在证据归档和 [receptor_coverage.csv](results/receptor_coverage.csv)。[RCSB 4EY6](https://www.rcsb.org/structure/4EY6) · [4EY7](https://www.rcsb.org/structure/4EY7) · [4NYX](https://www.rcsb.org/structure/4NYX) · [4QB3](https://www.rcsb.org/structure/4QB3) · [4J0S](https://www.rcsb.org/structure/4J0S) · [4RVR](https://www.rcsb.org/structure/4RVR) · [3SMR](https://www.rcsb.org/structure/3SMR)

**重原子评价与大环诊断。** 原有配对评价器读取 4EY6 加兰他敏时，无法识别 Meeko 标记的闭环碳。现将 CG0–CG3 作为真实碳，排除 G0–G3 构造点以及 H/HD；用于酰胺评分的供体氢仍由原始配体拓扑读取。RMSD 使用 CCD 化学图、立体化学和完整的实验配体坐标，保留全部有效同构映射，取重原子误差最小值。所有 7 个新结构均有完整配体，4EY7 有 4 个等价映射，逐结构映射数可在原始表中查看。交叉对接使用两侧口袋周围 14 Å 范围内共同的 Cα 确定旋转和平移。4EY6 的全部 11 个输出均保留 3 个指定手性中心。切断后重闭的 C–C 键在晶体中为 1.523 Å，预测为 1.544–1.570 Å。完整逐例结果见 [galanthamine_geometry.csv](results/galanthamine_geometry.csv)。

**复现与原始材料。** `transfer_evidence.tar.xz` 包含新结构 PDB 与配体 CCD、准备记录、每次对接的 PDBQT 候选、日志、模式表和时间。历史受体/配体由既有的三个归档自动恢复。源码位置为 `benchmarks/receptor_amide_states/transfer`。从仓库根目录执行下面的命令，其中 `VINA` 指向上述基线版本编译的二进制文件。编译使用 GCC 14.2.0 和 Boost 1.83.0，二进制 SHA-256 为 `82724f9518357fe6e7dc5d60445cb5e1dd1df976f476d3c0d2a75da0472f1176`。新材料使用 RDKit 2025.9.3、Meeko 0.7.1、OpenMM 8.6.1、PDBFixer 1.12.0 和 Gemmi 0.7.4。脚本需要 RDKit、Meeko、Gemmi、NumPy 与 SciPy。受体从头准备另需 OpenMM 和 PDBFixer。

~~~bash
D=benchmarks/receptor_amide_states/transfer
python "$D/run.py" restore /tmp/amide-transfer
OPENBLAS_NUM_THREADS=1 python "$D/run.py" score /tmp/amide-transfer --output /tmp/amide-transfer/scored

# 重跑全部 46 次对接（会更新该工作目录的候选及时间）
python "$D/run.py" dock /tmp/amide-transfer --binary "$VINA"
OPENBLAS_NUM_THREADS=1 python "$D/run.py" score /tmp/amide-transfer --output /tmp/amide-transfer/redocked
python benchmarks/pair_gradient/pose/test_heavy_atoms.py
~~~

若需要从 PDB / CCD 再生成输入，可在另一个工作目录运行 `OPENMM_CPU_THREADS=1 python "$D/run.py" prepare /tmp/amide-transfer-prep`。本轮已从归档恢复并逐池重新对接，46/46 个输出的 PDBQT SHA-256 重现，229 个模式的分数、RMSD 和 370 条位点贡献重现。计时单独存储于 [cost.csv](results/cost.csv)；历史面板首次对接最多并行执行 3 个单线程任务，相加的墙钟时间不能作为独占多核端到端耗时。逐池排序、完整 ECDF 和条件分布分别见 [cells.csv](results/cells.csv)、[ecdf.csv](results/ecdf.csv) 与 [summary.csv](results/summary.csv)。

针对 AChE 输入微态的追加反事实已完成。在两种化学状态下分别重跑对接，并在同一候选池上交叉评分。两个药物的结果方向和不同种子的效应大小差异较大，详见 [AChE 微态复测](ache/README.md)。
