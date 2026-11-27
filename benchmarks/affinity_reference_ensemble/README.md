# 受限刚体运动与单个配体扭转的耦合积分

固定已输出的首姿，通过测量势阱的宽度重新计算同一实验内配体的相对次序。现有 [rigid_window](../rigid_window/README.md) 已实现六维刚体积分。本轮新增一个从配体结构中选择的内部扭转，把它与平移和转动同时积分，并用相同角域中的游离配体内能作参照。所有节点均保留同一中心姿势、受体、化学图、质子化态和池内 unbound 参考。实验实现放在此目录。

## 物理模型

对由 PDBQT 声明的每个可旋转键，按拓扑切成两组重原子。各组绕该键轴的平方垂距之和为 $J_L,J_R$，单位为 Å²。选择 $J_LJ_R/(J_L+J_R)$ 最大的键，它在相对扭动时能移动两侧较多的原子。每次改变二面角 $\phi$，左、右部分分别绕同轴旋转 $-J_R\phi/(J_L+J_R)$、$J_L\phi/(J_L+J_R)$，然后恢复重原子质心。这种坐标选择最小化给定相对角变化时绕该键轴的一阶平方位移。其余内部角冻结。带大环断键标记的分子冻结全部内部几何，真实闭环碳参与受体相互作用，G 类辅助点从受体评分和质心中排除。

在重原子质心周围的 1 Å 平移球、15° 旋转帽和 ±30° 内部角域上取平均。平移服从均匀体积测度，外部旋转服从该帽内的 Haar 测度，转角累积分布为 $(\theta-\sin\theta)/(a-\sin a)$。内部二面角采用等宽中点。令 $E$ 为原生 Vina 显式受体相互作用，$I$ 为同一化学图的原生配体内部原子对能，$a_L$ 为原生配体柔性归一化斜率，$S_0$ 为中心姿势的原始输出。排序分数为

$$
S=S_0-RT\log
\frac{\langle \exp[-a_L(\Delta E(g,\phi)+\Delta I(\phi))/RT]\rangle_{g,\phi}}
     {\langle \exp[-a_L\Delta I(\phi)/RT]\rangle_{\phi}} .
$$

能量差均相对保存的中心姿势，单位 kcal/mol。计算温度取 298.15 K。分母使远离受体时配体形变贡献严格抵消。刚体与扭转可互相补偿，所以分别积分两个变量再相加，并不等价于联合积分。500 个构象中，这两种计算相差的绝对值中位数为 0.267 kcal/mol。既往研究也观察到整体运动与内部构象之间的耦合。[Chang、Chen 和 Gilson](https://pubmed.ncbi.nlm.nih.gov/17242351/) · [Chen 等人的多势阱方法](https://doi.org/10.1021/ct100245n) · [Vina 原始能量模型](https://doi.org/10.1002/jcc.21334)

这些积分仅覆盖一个已选盆地和一个局部内部角。分母仍以结合构象为中心，缺少全溶液自由配体、其他微态、水和受体重组。Vina 的经验坡度和受限能量不能校准为标准结合自由能。绝对域宽的选择会影响排名；本轮所有化合物使用相同设定。结论限于相对排序。

## 配对结果

Gal-1/3 使用 Peterson 等人的同一竞争荧光偏振实验。CHEK1、AURKA 和 S6K1 各自使用 Davis 等人同一靶点的直接报告 Kd。只对同靶点相差至少三倍的化合物对计分。每次先取同一受体的种子中位数，再对受体均分。Gal-1/3 各 13 个化合物，CHEK1 为 14 个，AURKA 去除外消旋 Ki-20227 后为 19 个。这四组历史池共 472 条运行，均有两个受体、四个种子。Gal-1/3 共享配体，几个激酶也共享部分配体。[Peterson 等，Table 1](https://doi.org/10.1021/acs.jmedchem.7b01626) · [Davis 等，Supplementary Tables 1、4](https://doi.org/10.1038/nbt.1990) · [历史激酶制备与测量](../torsion_kinase_affinity/README.md)

| 靶点 | 三倍差异对 | 原始 | 重原子数 | 原有刚体域¹ | 仅扭转 | 两项分别积分后相加 | 联合积分 | 同点集最低能量² | 高分辨率联合 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Gal-1 | 23 | 17 | 16.5 | 18 | 19 | 19 | 15 | 17 | 16 |
| Gal-3 | 35 | 28 | 34 | 32 | 27 | 33 | 34 | 28 | 34 |
| CHEK1 | 66 | 42 | 31 | 40 | 43 | 39 | 40 | 42 | 39 |
| AURKA | 136 | 98 | 87.5 | 91 | 97 | 95 | 92 | 98 | 91 |
| S6K1 | 58 | 31 | 14.5 | 37 | 35 | 36 | 37 | 31 | 36 |

¹ 为了比较相同提议，刚体消融用相同的 Sobol 外部坐标。S6K1 上额外用既有 CLI 的独立算法复核。² 取联合势能最低值，包括原始中心点，沿用同一原生参考。

Gal-3 增益可被分子大小基线追平，误差也未随排序改善，原始与联合积分的配对能隙 MAE 为 0.749 和 0.904 kcal/mol。AURKA 增益 4 对、退化 10 对。例如 PKC-412 对 ABT-869，原始预测的弱减强差为 +0.461，刚体修正后为 +0.098，联合修正后为 −0.027 kcal/mol。冻结坐标后跨过零点的最后一步是加入内部角积分。换用更密积分也未扭转激酶间互有得失的结果。当前数据支持保留独立诊断程序。

## 新靶点与算力对照

本轮增加了人 S6K1 的 3A60 链 A，采用不含点突变的未磷酸化激酶域，分辨率 2.8 Å。仅用共晶 staurosporine 的重原子质心定义 26 Å 盒子，并从评价中排除该分子。Davis 表中与此前 25 个已核验输入相交的 15 个精确 Kd 中，14 个进入对接。另有 10 个分子无精确 S6K1 数值。新的靶点没有用于选择扭转、设定系数或调整窗口，但分子在其他激酶上已使用。读数来自同一结合实验，磷酸化状态与测定构建体的实际构象仍不确定。结构中未观察到 73–84、130–140、243–255、372–399 段。最接近共晶配体的可见边界原子分别距 19.42、13.34、17.48、10.82 Å，实际缺失原子的位置未知。下载原始 PDB 和每个排除理由均可复查。[3A60](https://www.rcsb.org/structure/3A60) · [实验排除记录](results/s6k_panel.csv) · [坐标质量复核](results/s6k_structure_audit.json)

每个化合物使用相同的两个种子。上游 AutoDock Vina `3c65c0b` 和当前 Pro 均为两条轨迹、每条上限 25,000 次求值；较大预算 Pro 使用四条同样预算的轨迹。共进行了 84 次单 CPU 对接。上游、当前 Pro、较大预算 Pro 分别答对 33、31、31 对。同点集最佳能量仍为 31 对。联合积分在预设 1,024 点时增加 6 对，5,000 次按配体重采样的准确率差描述性区间为 0 至 +26.7 个百分点。原始、刚体和联合的 Spearman 分别为 0.068、0.286、0.290，能隙 MAE 为 2.216、1.874、1.926 kcal/mol。在 CEP-701（69 nM）对 PKC-412（1,300 nM）中，实验差对应 +1.740 kcal/mol，原始、刚体、联合依次给出 −0.025、+0.816、+0.978。该指标改善可由既有刚体域解释；加入内部角没有进一步改善全组计数。拿掉带断环标记的 CID 153999 后，S6K1 原始 29/48，联合与刚体均为 35/48；CHEK1 原始 38/58，联合与刚体均为 36/58。这是见到结果后增加的排除对照。[逐对结果](results/primary_pairs.csv.xz) · [大环复核](results/macrocycle_replay.csv)

三种搜索在本机三进程并发下每配体每种子墙钟中位数分别为 1.668 s、0.392 s 和 0.604 s。后处理只用 CPU，不产生新 docking pose。预设 1,024 个刚体/扭转节点的联合臂每个非环首姿新增 1,024 次受体相互作用能和 16 个内部角求值；同时运行刚体与单扭转控制还会分别增加相应求值。实现缓存各离散内部态。四个历史组含全部消融的最高 2,048 节点共耗 32.93 CPU 秒、11.35 墙钟秒（3 个进程），不含历史对接和准备。受体文件在同一靶点可复用，逐配体保留了准备图与对应原始池。本轮未重新计时历史配体准备或 S6K1 受体准备。另有一个预设化合物的各臂串行计时与峰值 RSS。[全量日志与制备](results/s6k_evidence.tar.xz) · [串行成本](results/s6k_serial_cost.csv)

## 数值核验与复现

两次积分分别采用 scramble 40721、16 个内部角、1,024 点和独立 scramble 66173、32 个内部角、8,192 点。500 个单细胞修正的绝对变化中位数 0.024，90 分位 0.081，最大 0.582 kcal/mol。后者足以改变接近打平的配体对，密集积分在 S6K1 使 37/58 变为 36/58。此检验同时改变角离散、节点数和 scramble，不能将差异唯一归因于其中某项。逐级节点扫描、相对排序、Kendall τ、带共享分子依赖的 bootstrap、逐分子结果和全对得分都在 `results/` 下。

零受体势时分子内部参考准确抵消；禁止刚体扰动时联合修正退化为单角积分；环关闭内自由度时联合修正退化为刚体积分。全部十个带标记的池首姿中，所有点对距离（包括 G 辅助点到对应闭环碳的距离）在刚体变换后的最大变动为 1.34×10⁻¹⁴ Å。用正式程序的 `--pose_coordinates` 和各池原始 unbound 复算这些中心总分，最大误差 0.007 kcal/mol。其余构象的显式 INTER 在全部历史池内与原始三位小数输出最大相差 0.013 kcal/mol。S6K1 中正式刚体积分支持的 26 个构象，原生与本实现的中心 INTER 最大相差 4.1×10⁻¹⁰；两种独立刚体积分相差最多 0.024 kcal/mol。原有刚体接口会拒绝另外两条带闭环标记的输出。

以仓库根目录为当前目录。下列编译使用仓库当前相同的 Boost 头文件与库；使用自定义前缀时，给两个编译步骤都提供相同前缀。从归档复核评分和结构质量无需联网。前三条实验使用的历史原始池已随仓库签入；新靶点的原始 PDB、输入、姿势和 stdout 位于 S6K1 压缩档。

```bash
B=benchmarks/affinity_reference_ensemble
W=$(mktemp -d)
python tools/build.py --output-dir "$W/build"
c++ -O3 -std=c++14 -DNDEBUG -Isrc/lib "$B/integrate.cpp" \
  $(find "$W/build/obj" -name '*.o' ! -name 'main.o') \
  -lboost_thread -lboost_filesystem -lboost_program_options -lboost_serialization \
  -pthread -o "$W/integrate"
python "$B/run.py" "$PWD" "$W/integrate" "$W/history" --workers 3
python "$B/checks.py" "$W/integrate" "$W/history"
python "$B/run.py" "$PWD" "$W/integrate" "$W/check" --workers 3 --nodes 8192 --bins 32 --scramble 66173
mkdir "$W/s6k"
tar -xf "$B/results/s6k_evidence.tar.xz" -C "$W/s6k"
cp "$B/results/s6k_panel.csv" "$W/s6k/panel.csv"
python "$B/s6k.py" dock "$PWD" "$W/s6k" --integrator "$W/integrate" --pro "$W/build/vina" --native /path/to/upstream/vina --workers 3
python "$B/s6k.py" sensitivity "$PWD" "$W/s6k" --integrator "$W/integrate" --workers 3
python "$B/audit.py" "$W/s6k" "$W/build/vina"
python "$B/macrocycle_audit.py" "$W/history" "$W/s6k" "$W/build/vina" "$W/macros"
python "$B/analyze.py" "$W/tables" "$W/history/cells.csv" "$W/s6k/cells.csv" --docks "$W/s6k/docks.csv"
python "$B/analyze.py" "$W/checktables" "$W/check/cells.csv" "$W/s6k/sensitivity.csv" --docks "$W/s6k/docks.csv" --detail-nodes 8192
```

新靶点制备可用 `s6k.py prepare REPO WORK` 重做，先把 [Davis 来源表](../kinase_state_affinity/source/sources.json) 中的 `table1.xls` 和 `table4.xls` 放入 `WORK/raw/`，与原始 3A60 PDB 同目录。程序核对三个文件的 SHA-256。制备用 Python 3.12、Meeko 0.7.1、OpenMM 8.6.1、PDBFixer 和 xlrd。历史配体直接复用已签入准备图，源档的散列也要通过。排名复算使用 NumPy 2.3.5、SciPy 1.17.0 与 pandas；本轮使用 x86-64 AMD EPYC 9V74，C++14，单线程受体能量求值。宏环自适应构象、绝对亲和力、实验构象配准及多受体热力学权重均需其他证据。
