# 相同受体上的官方 Vina、双预算 Pro 及独立坐标复核

用[已锁定的三个受体–配体对与种子](../mutation_scan/PROTOCOL.md)，补充官方 AutoDock Vina 1.2.7 与两倍每轨迹评价上限的默认 Pro。所有臂读取相同的受体、配体 PDBQT、盒中心及尺度。官方程序沿用其内部分型、拓扑推断、精修、去重与选姿；AT9283 的 Pro 臂另接收已核实的完整配体边表。两种程序的名义限额相等时，实际迭代数、求值次数与运行时间仍可不同。官方臂的真实目标函数求值计数未知。完整图 Pro 的两档上限（25,000/50,000，每次两条轨迹）才能在相同引擎和输入下单独考察预算。

## 结果

每项统计含三种配体各两枚新种子，均使用受体坐标系内、化学对称校正的重原子 RMSD，单位为 Å。历史 AT9283/61231 单列。空池均为零。两个 Aurora A 配体共享一个接收构建，BRD9 对应第三种配体、一种不同家族。组内复本相关，六个样本的 P90/P95 只描述实际尾部。Q25 到 P95 采用线性插值。[每配体及整体分布](results/distributions.csv) · [14 次原始坐标与日志](raw.tar.xz) · [逐运行与成本](results/cells.csv) · [每个模式的误差及键长](results/modes.csv)

| 方法 | 首姿均值 | Q25 | 中位数 | Q75 | P90 | P95 | Best@5 均值 | Oracle 均值 | Regret 均值 | 首姿 ≥8 Å | 每格平均墙秒 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 默认 Pro，完整图 | 8.601 | 8.118 | 8.501 | 8.983 | 10.447 | 11.158 | 7.668 | 7.321 | 1.280 | 5/6 | 0.262 |
| 无放回提议，完整图 | 8.596 | 8.170 | 8.487 | 9.967 | 11.739 | 12.408 | 7.327 | 7.327 | 1.269 | 5/6 | 0.299 |
| 双预算默认 Pro，完整图 | 6.994 | 6.403 | 7.664 | 8.668 | 8.867 | 8.873 | 6.775 | 6.616 | 0.378 | 3/6 | 0.386 |
| 官方 1.2.7，原生推断图 | 5.822 | 3.702 | 6.382 | 8.531 | 8.992 | 9.012 | 5.282 | 5.119 | 0.704 | 2/6 | 1.304 |

AT9283 历史种子的双预算 Pro 首姿由原来的 8.717 恶化至 13.196，Oracle 仍为 7.583。官方这格为 8.767/8.691。三种药物各自的表现和同一药物的种子差异均可直接检查。双预算 Pro 新种子中，danusertib 的 Oracle 分别为 2.602、7.135。官方在两枚 82I 种子上获得 2.993、1.082 的 Oracle，两个 AT9283 池最接近的也只有 6.990、8.057。三种经验搜索与有限模式池均未模拟溶液中的配分函数，因而这些数据只评价几何。

## 独立复核

逐字节验证三个体系的接收受体、配体 PDBQT 和显式图与主试验相同。根据原始晶体重新取蛋白口袋骨架，用 Gemmi 刚体配准，配体由 PDB CONECT 和元素恢复，RDKit 在固定受体框架计算 CalcRMS。BRD9 两份结构链 A 的共同结构域作者编号相差 116；同时要求每个 Cα 有相同残基名、距各自配体中心都小于 14 Å。37 个锚点配准 RMSD 为 0.392 Å。查询配体为 6YQW 中 1.50 Å 分辨率观察到的 82I；另一个配体 TP-472 所在的 6Y7L 提供 1.80 Å 分辨率接收受体。[6YQW 原始沉积](https://www.rcsb.org/structure/6YQW) · [6Y7L 原始沉积](https://www.rcsb.org/structure/6Y7L) · [82I 研究原文](https://doi.org/10.1021/acs.jmedchem.0c00075)

Aurora A 用 2X6E 的口袋 Cα 作为靶框架，2W1G 和 2J50 的 A 链作相应参考。2W1G 的 AT9283 配合物处于 DFG-in。6CPG 则含稳定抑制构象的 monobody，实验观察存在跨状态结合；本轮 2W1G 几何比较没有估计两构象的溶液占比。[2W1G 原始沉积](https://www.rcsb.org/structure/2W1G) · [两种构象与动力学实测](https://doi.org/10.7554/eLife.36656)

35 份完整池（官方/双预算 14 份、Pro/等价观察器/无放回 21 份）的 312 个发布模式全部复算。独立实现和主实现的 RMSD 以及衍生统计最大差 2.49×10⁻¹⁴ Å。在已检查的 SDF/SMILES 重键上，相对输入 PDBQT 的最大漂移为 0.00127 Å；90 个 danusertib 已指定中心的有向体积均保留原符号。等价观察器重复产物计入软件复核总数，没有把它加入独立统计样本。[独立复核脚本](verify.py) · [数值结果](results/verification.json)

另按原作者两份 XLS 的归档 SHA-256 核实了 Davis 面板六格。Table 1 中 AURKA、CHEK1、S6K1 分别列为 NP_003591.2、NP_001265.2、NP_003152.1，Mutant 均为 NO；按 Table 4 的同 accession 和同激酶行取值，六格依次为 3.9/1700、2.8/8700、48/69 nM。用 298.15 K 作为浓度比的换算约定，三个有向 RT ln 浓度比相加为 −8.580231575 kcal/mol。突变标记不确定检测构建体的截短范围、磷酸化或游离构象占比。重复执行幅度证书的独立验证脚本，12 份可行见证与 267 个覆盖区域全部通过。这部分核实了既有测量与约束，新搜索的交叉对接没有给这六格产生可比的自由能预测。[原始论文及附表](https://doi.org/10.1038/nbt.1990) · [六格复算](results/kinase_six.json) · [幅度证书](../../../../affinity_exchange/joint_reference/local_integral/magnitude/README.md)

## 已经访问的内部构象能否通过刚体运动达到晶体附近？

在阶段结论出来后追加这项诊断。直接读取已记录的全部真实提议和 hunt 局优返回，对每个重原子坐标集合在合法化学对应上取最小值。设原始受体框架误差为 $D$，两配体重原子质心距离为 $C$；只对齐质心后的误差为 $T$，允许任意平移和正行列式刚体旋转时的最小误差为 $S$。重心化的平方和恒等式给出 $D^2=C^2+T^2$。只允许绕预测质心改变取向的最小误差为 $\sqrt{C^2+S^2}$；报告的方向残差为 $\sqrt{T^2-S^2}$，单位均为 Å。$T$ 和 $S$ 可各自优化图的合法对称映射。这些叠合坐标读取了晶体标签，未经过受体排斥或弛豫检验，不能用于替代实际姿势指标或定义真实可结合构象。所有表同时保留受体框架原值。[按种子与阶段计数](results/placement.csv) · [每类极值的轨迹、步与误差分解](results/placement_witnesses.csv)

| 完整图，五种子的 hunt 返回 | 实际 ≤2 Å | 仅向晶体平移可达 ≤2 Å | 仅绕当前质心旋转可达 ≤2 Å | 同时转动及平移可达 ≤2 Å | 质心距 ≤2 Å |
|---|---:|---:|---:|---:|---:|
| AT9283，8,058 个返回 | 0 | 53 | 0 | 7,172 | 0 |
| danusertib，9,274 个返回 | 0 | 0 | 0 | 3,937 | 3 |

AT9283/61231 的 1,638 个完整图返回中，有 1,419 个内部构象经晶体引导的刚体放置可达到 2 Å。实际最近返回是轨迹 1 第 297 步，$D=7.238$、$C=6.285$、$T=3.591$、$S=1.390$ Å。只移动该姿势的质心仍超过 2 Å；另有四个返回只需正确平移即可达到这个标准。完整图五枚种子的返回中最近质心也相距 2.864 Å，但其中有一次原始提议（45104）质心距 1.631 Å，局优后没有保留如此接近的质心。对 danusertib 有三个返回质心够近，却没有一个返回姿势仅作平移或仅作转动即可达到 2 Å。将刚体放置作为下一次盲干预有明确的诊断动机，单独增加质心提议仍须跨分子检验，且应在相同评价次数内比较。这里的成千条记录发生在少数相关轨迹上，不能换算为跨分子成功率。

还核对了[交叉局优](../proposal_crossover/PROTOCOL.md)的全部 6,224 次起始与返回。相同来源提议在两张接收图上的 12 维广义坐标和三维原子坐标逐项相等。同图对角线与盲搜的对应原始轨迹逐步、逐原子完全一致。改用独立晶体配准复算，每个误差最大相差 3.29×10⁻¹⁴ Å。两张图均在历史漏图来源的轨迹 2 第 706–710 步返回 ≤2 Å，在完整图来源的 1,638 个提议上均返回零个。交叉实验没有把受体或原子坐标改变归到配体图上。[逐命中独立复核](results/crossover.json)

## 成本与复现

完整七条件上，官方七个成功进程累计 9.007 墙秒/8.988 CPU 秒，双预算七个成功进程累计 2.789/2.781 秒、701,911 次含末端精修的计数求值。成本包括每次启动与格点准备，独立读坐标、上一轮受体建模、构建与归档未计入这些进程时间。第一次官方 BRD9 运行在准备格点时被系统 SIGKILL；未生成坐标，未记录其独立成本和终止原因。相同输入在独立输出路径完整重跑；失败日志保留。[中断时的日志](results/interrupted_grid.log) · [程序与输入哈希](results/manifest.json) · [Vina 方法](https://doi.org/10.1002/jcc.21334) · [官方实现](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7)

以下在仓库根目录运行。WORK 是已按[上级说明](../REPRODUCE.md)解开的 Aurora 输入目录。T 为新的绝对路径。PRO、OFFICIAL 分别指默认 Pro 和发行版官方 1.2.7 二进制。SCAN 是本轮 mutation_scan/run.py 生成的目录，含 inputs.json、raw/、modes.csv、cells.csv。也可使用该试验归档的运行目录。独立复核可跳过两种程序的运行，使用已归档的 pools/、82I.bonds 和 results/ 内的 cells.csv、modes.csv 组装 control 目录即可。

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
D=benchmarks/receptor_linkage/aurora/internal_loops/mutation_controls
mkdir -p "$T/brd"
tar xf benchmarks/torsion_proposals/fixtures/inputs.tar.xz -C "$T/brd"
python "$D/run.py" "$WORK" "$T/brd" "$T/control" "$PRO" "$OFFICIAL"
python "$D/verify.py" "$WORK" "$T/brd" "$T/control" "$SCAN" "$T/checked"
cmp "$T/checked/verification.json" "$D/results/verification.json"
K=benchmarks/affinity_exchange/joint_reference/local_integral/magnitude
python "$K/verify.py" "$K/results"
# DAVIS 含原始 table1.xls 与 table4.xls，二者 URL 与哈希见下面的来源表。
python "$D/check_six.py" "$DAVIS" > "$T/kinase_six.json"
cmp "$T/kinase_six.json" "$D/results/kinase_six.json"
```

[测定原表下载地址与哈希](../../../../kinase_state_affinity/source/sources.json)


在 [basin_trace](../basin_trace) 复现命令生成未压缩的完整坐标流后，令 STREAMS 指向存有各个 *.bin 的目录，CROSSED 指向含四份 *_on_*.npz 的 proposal_crossover 结果：

```sh
B="$D/../basin_trace/results"
python "$D/placement.py" "$WORK" "$STREAMS" "$B" "$T/placed"
cmp "$T/placed/placement.csv" "$D/results/placement.csv"
python "$D/verify_crossover.py" "$WORK" "$CROSSED" "$STREAMS" > "$T/crossover.json"
cmp "$T/crossover.json" "$D/results/crossover.json"
```

位姿极值经 SVD 运算，不同 BLAS/平台可引起最后若干小数位差异。部署到不同数值环境时应逐行检查见证步和容差。
