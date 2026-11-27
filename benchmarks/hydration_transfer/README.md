# 两种人源受体上的有序水交叉对接验证

同一批 2 Å 候选在首位、前五和输出池中的出现次数，可以区分水重排改善了打分还是改善了采样。本实验在同一套原生 Pro 候选上分别使用原生、固定位点和变位水排斥重排。接收结构给出水的位置，查询结构只供评价。六个查询–接收配对、另两个自对接对照、预算和指标均在产生预测前公开冻结。[登记提交](https://github.com/HankerWu/fastdock-vina/commit/d8cd280618f7fc9ebbfa1724424a818352fe6c6b) · [物理模型与开发集](../mobile_water/README.md)

BRPF1 的晶体结构系列提供了保守水被配体置换的实验例子。验证使用 5C7N 链 A 的 1.75 Å 结构作为接收体。EPHX2 的接收体 4OD0 分辨率为 2.92 Å，解析出的水较少。每个接收体各有一个满占有率、无交替构象的供体配体定义盒中心。每个查询从 CCD 化学图准备。质子化态沿用该图，例如 3WKE 的 AUB 在这里以中性羧酸输入；实际 pH 下的羧酸盐比例需要另行确定。受体固定，输入保留 4OD0 的 Mg。跨体系的水参数全部相同。[BRPF1 原始结构研究](https://doi.org/10.1021/acs.jmedchem.6b00215) · [5C7N](https://www.rcsb.org/structure/5C7N) · [4OD0](https://www.rcsb.org/structure/4OD0) · [各查询、配置和信息权限](PROTOCOL.md)

## 前瞻冻结的交叉对接

下表各格为首位命中 / 前五命中 / 池内至少一项命中。每行三个查询、两个种子，共六次。官方 Vina、Pro、固定水及移动水四臂在这些阈值指标上得到同样的计数。原生有独立的搜索输出；两种水模型复用 Pro 的坐标。

| 接收体 | 每轨迹 3,000 次名义上限 | 每轨迹 12,000 次名义上限 |
|---|---:|---:|
| BRPF1 5C7N | 0 / 0 / 0 | 0 / 1 / 1 |
| EPHX2 4OD0 | 1 / 1 / 1 | 1 / 2 / 2 |

较高预算在 12 次交叉对接中得到 3 个有正确候选的池，1 个排首位。增加预算主要扩大池的覆盖。BRPF1 的低预算种子 43117 下，5FG5 的移动水首位 RMSD 从 4.463 增至 7.947 Å。5EPR 在高预算两种子下分别从 4.682 降至 3.931 Å、从 4.054 升至 4.662 Å。计数相同掩盖了错误姿势之间的换位。两个自对接对照单列在结果中。没有按靶点作统计推断，亦没有将同一配体的多个种子视作独立化合物。[逐姿势](results/poses.csv) · [逐池](results/cells.csv) · [汇总及条件首位命中](results/summary.csv)

## 晶体起点的反事实

独立运行局部优化。先把晶体扭转角转到准备好的查询分子，再整体对齐重原子。选择投影初始误差最小的合法化学图映射。比较局部优化结束后的同一目标 INTER + INTRA 与盲池最优值，使同一查询的参考零点一致。下表取四个预算/种子组合的范围。负值代表晶体辅助起点得到的局部姿势更有利。使用晶体信息的结果只作诊断。

| 查询 | 精修后投影姿势 RMSD，Å | 局部姿势减盲池最优总能的范围，kcal/mol |
|---|---:|---:|
| 5EPR | 0.818 | −0.11 至 +0.51 |
| 5EPS | 1.249 | −2.35 至 −0.10 |
| 5FG5 | 0.408 | −0.77 至 +1.44 |
| 3WKE | 0.622 | −5.30 至 +0.04 |
| 4X6Y | 1.115 | +3.77 至 +3.95 |
| 3ANT | 1.472 | +6.47 至 +7.37 |

5EPS 的四个干池都缺失几何正确的候选，而投影姿势在原目标中均优于池内最优。4X6Y 和 3ANT 的投影姿势同样位于 2 Å 内，原目标却偏爱池中其他姿势。3ANT 在直接使用晶体几何和在查询刚性片段上移植晶体扭转时，所得局部 INTER 分别为 −11.587 和 −4.313，INTRA 为 −0.121 和 −0.024 kcal/mol。主要差值在蛋白–配体相互作用项；两个局部终点和初始内部几何均有差异，这个对照无法单独确定某个键角或某一环构象的因果贡献。[各局部起点](results/crystal_starts.csv) · [逐次配对](results/oracle_starts.csv)

## 计算和复核

共运行 64 次单 CPU 对接，包括自对接，得到 560 个姿势，其中 Pro 279 个、官方 Vina 281 个。两档预算的原生 Pro 单查询中位墙钟为 0.181 和 0.325 秒，官方 Vina 为 1.059 和 1.203 秒。这些差异包含分支以往的工程优化。Pro 四条轨迹的实际搜索求值中位数为 12,191 和 48,099 次，显式精修中位数分别为 440 和 449.5 次。原生官方程序未输出等价计数。子进程峰值内存读自 Linux wait4，包含启动开销。[每个命令、CPU、墙钟与内存](results/cost.csv) · [二进制哈希](results/binaries.json)

5C7N 和 4OD0 受体准备分别耗时 1.90、6.37 秒，水筛选得到 12 和 3 个位点。11 阶积分网格预计算耗时 0.0079 和 0.0021 秒，另有解析与选水时间。279 个固定姿势的固定/11 阶移动排斥合计耗时 0.240 秒，17 阶检查耗时 1.036 秒。两个独立复核入口分别使用原始 PDB + CCD 和 RCSB 配体实例 SDF + PDBQT SMILES 映射。配体序号按准备后几何核实已声明立体中心。两条路径重算全部 560 个重原子 RMSD，最大差值 3.02×10⁻¹⁴ Å。[软件版本](results/environment.json) · [准备成本](results/preparation.csv) · [水成本](results/receivers.csv) · [坐标核验](results/verification.txt)

## 复现

在仓库根目录操作。归档含原始接收/查询 PDB、CCD、供独立验证的配体实例 SDF、准备好的输入、所有输出及局部诊断。环境版本见 environment.json。结果记录所用原生源码为 ccsb-scripps/AutoDock-Vina 的 3c65c0b，Pro 为本仓库 f0e6dbf。

```bash
B=benchmarks/hydration_transfer
W=$(mktemp -d)
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tar -xJf "$B/inputs.tar.xz" -C "$W"
tar -xJf "$B/outputs.tar.xz" -C "$W"
python "$B/evaluate.py" "$W" "$W/replay"
python "$B/verify.py" "$W" "$W/replay"
```

从头重算时在另一个空工作目录依次执行以下命令，NATIVE 和 PRO 指向上文两个版本。几何和分数可逐行对照，耗时取决于机器。清单中 3WKE 在旧研究曾被列为备选，这六个交叉配对仅对本水模型作留出；在据此设计后续搜索算法时它们归入开发集。

```bash
python "$B/panel.py" "$W" prepare
python "$B/panel.py" "$W" dock --native "$NATIVE" --pro "$PRO"
python "$B/panel.py" "$W" geometry
python "$B/oracle.py" "$W" "$PRO"
python "$B/evaluate.py" "$W" "$W/replay"
python "$B/verify.py" "$W" "$W/replay"
```
