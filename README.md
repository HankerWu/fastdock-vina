# FastDock Vina Pro

AutoDock Vina 的 CPU 优化和对接研发分支。

[性能实验记录](docs/BENCHMARKS.md) · [逐项数据](benchmarks/benchmark.csv) · [源码来源](docs/PROVENANCE.md) · [上游 AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)

FastDock Vina 基于 [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) 做 CPU 优化。Pro 研发分支继承这些优化，并提供可选的搜索、精修和姿势重排方法。

在线保存同一姿势的多个邻居时，默认实现已修复只比较最近代表造成的重复。在独立结构面板中，一个原已生成的 FABP4 正确候选由此进入了精修名额。两项改善都来自同一复合物，首姿命中没有提高。[候选保存、实测 Kd 与固定状态复算](benchmarks/archive_overlap/README.md)

此前在三个新的人源靶点检验了刚体主轴翻转，在三个独立大环复合物检验了协同扭转。ESR1 的交叉对接与 HSP90α–3VHA 有较低实耗下的局部收益；Galectin-3 的 47 配体 Kd 排序出现退化。生产版本继续沿用当前搜索。逐病例结果与复现入口见[刚体提议](benchmarks/shape_flip/README.md)、[独立结构](benchmarks/rigid_flip/README.md)、[大环协同转动](benchmarks/concerted_ring/README.md)和[游离态参照实验](benchmarks/frozen_core/solution_strain/README.md)。

上一轮的[大环提议与持续约束的等成本实验](benchmarks/macrocycle_proposal/README.md)比较了上游、当前 Pro 和三个 C++ 研究版本。持续约束消除了本面板中的严重断环，唯一新增的 2 Å 命中来自近晶体教程输入。新靶点与同系列实测 Kd 均不足以支持替换默认搜索。

[小环键距、大环闭合与同测定 Kd 的验证](benchmarks/bond_distance/README.md)复现了四元环原子编号改变能量和梯度的错误，默认实现已修正。同一组盲对接中出现了姿势退化。另一个研究变体在较少求值下消除了雷帕霉素的两次严重断环，尚未提高整分子姿势命中率。

[大环闭环编号修复](benchmarks/macrocycle_labels/README.md)消除了同一分子重编号导致的评分与对接差异，并修复了多配体之间的虚构闭环。[首步搜索与候选覆盖实验](benchmarks/hunt_overlap/README.md)测量了严重排斥区的网格梯度、盒外起点和同测定 Kd 排序。两项搜索方案均保留为研究记录。

[搜索阶段 BFGS 的固定预算验证](benchmarks/search_backtracking/README.md)复现了回溯耗尽时接受上坡试步的问题。修正后，三个靶点的独立结构面板仍为 3/40 个首姿命中；两项同测定激酶 Kd 的结果也不支持更换默认搜索。另一次[因子 Xa 的单配体几何修复与重对接](benchmarks/frozen_core/strain_transfer_audit/README.md)将原生与 MMFF 的排序差从 17 对缩小到 1 对。

[固定精修名额的候选保存研究](benchmarks/accepted_retention/README.md)追踪了已接受却未进入精修的姿势。在一条 FABP4 轨迹中，保留一个新颖接受态恢复了次名近原位姿势；八家族的新种子测试未提高首姿或候选池命中率。

上一轮[径向梯度与 8 Å 截断的配对实验](benchmarks/pair_gradient/README.md)定位了真实受体上的旋转导数符号错误，并检验了三种 C++ 修正。新结构及同测定 Kd 结果不足以替换默认算法。

![25 项测试的加速比及汇总指标](docs/figures/speedup-summary.png)

## 为什么用 FastDock Vina

在 CPU 集群上执行大批量虚拟筛选时，单次对接的耗时和内存需求会累积。常规 FastDock 的工程优化保留双精度和原搜索逻辑。

## 效果

下表记录发布版的历史结果。后续改动和复测见[工作线程复用](docs/BENCHMARKS.md#工作线程复用的本地复测)
和[受体查表与空间索引](docs/BENCHMARKS.md#受体查表与空间索引复用)。
近期对照见[刚性受体共享和建图复用](docs/BENCHMARKS.md#刚性受体共享和建图复用)、
[按使用需求生成网格](docs/BENCHMARKS.md#按使用需求生成网格)、
[按配体类型准备受体评分表](docs/BENCHMARKS.md#按配体类型准备受体评分表)、
[亲和力网格读取与姿势写出](docs/BENCHMARKS.md#亲和力网格读取与姿势写出)
和[局部搜索向量与网格索引](docs/BENCHMARKS.md#局部搜索向量与网格索引)。

测试平台为 Intel Xeon Platinum 8362，每个案例重复 3 次并取中位数：

| 指标 | 结果 |
| --- | --- |
| 25 项测试的等权几何平均加速比（19 项完整对接 + 6 项评分 / 局部优化） | **1.228×** |
| 19 项完整对接的合计耗时（原版 3242.525 秒 → 2989.231 秒） | 降低 **7.81%**（约 1.085×） |
| 峰值内存（RSS） | 降低 **28.71%–73.27%** |
| 对接输出 | 已测配置下与原版**逐字节一致** |

两个速度指标口径不同：1.228× 给每项测试相同权重；7.81% 是 19 项完整对接的总耗时变化。19 个完整对接案例都比原版快。所有结果只适用于已测的环境和输入，详细数据与验证范围见 [实验记录](docs/BENCHMARKS.md)。

![1iep 配体在受体中的位置，以及原版与 FastDock 构象的叠加](docs/figures/docking-overlay.png)

上图右侧放大了 MODEL 1：原版为蓝色，FastDock 为橙色，两者坐标完全重合（40 个原子的 RMSD 为 **0.000 Å**，未做结构对齐）。[对接输出](examples/basic/1iep_docked.pdbqt) · [图表数据与复现方法](docs/FIGURES.md)

## 做了哪些优化

- 调整 BFGS 矩阵布局，用 CPU SIMD 同时处理相互独立的行，并保留原有的累加顺序。
- 按实际进入受体空间索引和配体相互作用列表的类型对构建 Vina/Vinardo 查表。普通相互作用按实际截断距离采样，大环闭合项保留 20 Å 查表。
- 显式受体评分只为当前配体和柔性残基需要的类型建表。亲和力网格复用这些表，批处理遇到新类型时补表并继续复用已有表。
- 组装模型时确定原子对的查表索引。AD4 只采样参与评分的原子对，每对保留独立的电荷信息。局部搜索并行插值成对能量和径向导数，保留双精度和原子对累加顺序。
- 同一网格行的相邻点复用受体候选和两个距离分量，并跳过超出整段截断范围的受体原子。大环闭合碳复用普通碳的网格文件。
- 建图时，同一受体类型和距离对应的各配体类型能量连续存放。显式评分和分子内评分直接从相邻能量样本求均值。
- 建图按受体原子遍历网格行，仅处理截断距离内连续的一段格点。各格点仍按受体原子顺序累加。
- 局部搜索复用相邻原子对共享的首原子坐标，在该组结束时写回其力。每组评分表只保存一次距离采样倍率。
- 同一搜索轨迹或精修循环复用 BFGS 的构象、梯度和矩阵缓冲区。局部优化直接使用最后一次求值生成的坐标，回退时恢复起始坐标。
- 构象树直接写入 BFGS 的连续梯度数组，试探构象直接从起始构象和搜索方向生成。最后一步省去后续不会使用的逆 Hessian 更新。
- 网格插值只计算一次起始角点索引，沿现有行和层跨度读取七个相邻采样值。
- 搜索线程和配体模型共享刚性受体。柔性受体的固定锚点需要重编号时，才复制相应的受体拓扑。
- 每次建图为每个受体原子查找一次横向格点分界，然后投影到每段网格行。保留每个格点的受体原子累加顺序。

## 编译

需要 Python 3.8+、支持 C++14 的 GCC 或 Clang，以及 Boost（已验证 1.82）。

```bash
git clone --branch fastdock-vina-pro https://github.com/HankerWu/fastdock-vina.git
cd fastdock-vina
```

Ubuntu / Debian：

```bash
sudo apt-get install build-essential python3 libboost-all-dev
python3 tools/build.py --jobs 4
```

Boost 的头文件需包含 `boost/random/sobol.hpp`，例如 Ubuntu 22.04 的 Boost 1.74。

macOS（需要 Xcode 命令行工具和 Homebrew）：

```bash
brew install boost
python3 tools/build.py --boost-prefix "$(brew --prefix boost)" --jobs 4
```

## 运行

用仓库自带的 1iep 示例跑一次对接：

```bash
mkdir -p demo_out
./build/vina \
  --receptor examples/basic/1iep_receptor.pdbqt \
  --ligand examples/basic/1iep_ligand.pdbqt \
  --config examples/basic/1iep_receptor.box.txt \
  --scoring vina --cpu 1 --seed 20261022 \
  --exhaustiveness 8 --max_evals 0 --num_modes 9 \
  --out demo_out/1iep_out.pdbqt
```

结果保存在 `demo_out/1iep_out.pdbqt`，更多参数见 `./build/vina --help`。示例使用 `exhaustiveness=8`，上面的性能测试使用的是 `32`。

Vina 和 Vinardo 的普通相互作用截断距离默认是 8 Å。可用
`--interaction_cutoff 6.75` 覆盖，或在 `--config` 指向的文件中写
`interaction_cutoff = 6.75`。命令行显式值优先。取值大于 0 且最多 20 Å，
并统一用于普通评分势、受体及配体查表。自定义值用于直接读取受体并构建网格的运行。
预计算的 Vina/Vinardo 网格使用默认截断距离。AD4 和大环闭合项各使用其力场定义。

Vina 和 Vinardo 的 `--score_only` 可直接用刚性受体原子评分，省下亲和力网格的
计算和存储。`--ligand` 模式中省略六个盒参数时，固定姿势评分默认使用
`--autobox`。关闭显式评分（`--no_refine`）、默认参考能量涉及柔性残基、
以及搜索和导出时仍会构建所需网格。

从受体直接执行 Vina 或 Vinardo 的 `--batch` 时，各配体共享受体空间索引、
搜索盒和已建网格。首次用到某个配体或柔性残基的原子类型时才建图。批次带
`--write_maps` 时会一次导出完整的类型集合。`--maps` 加载预计算网格，
AD4 仍需对应的原子类型、电势和去溶剂化网格。导入会检查网格头部和采样值数量，
采样值须为有限数。姿势文件逐个写出模型，PDBQT 格式保持一致。

### 柔性受体的固定锚原子

为 Vina 或 Vinardo 提供同一受体的无重叠刚性和柔性 PDBQT 拆分。柔性文件里的
固定轴锚（通常是 Cα 和 Cβ）保留原位，和配体的相互作用进入显式评分与搜索梯度，
并计入 `Ligand - Flex side chains`。读取预计算地图时，使用同一刚性部分生成的地图。

用 `--score_only --unbound_energy 0` 复算固定坐标，读取 `Final Intermolecular Energy`
及其分项。全刚性评分和柔性显式配对的正能量软化聚合顺序不同。原生去重还包含可动侧链，
与实验配体做 RMSD 比较时只取每个模型中 `BEGIN_RES` 之前的配体原子。

### 扭转搜索、受体装配与同靶点实测排序

在 17 种人源蛋白上检验了搜索提议、受体装配和固定姿势排序。TTR 四聚体可恢复 8DW5 中实验支持的同口袋姿势，甲状腺素未复现。两种可选扭转提议跨靶点的首姿和 Ki/Kd 排序增益随受体、配体及实验条件变化。研究分支提供原始输入、逐种子配对、受体坐标独立复算及可执行的装配工具。[结果、公式和复现入口](benchmarks/torsion_search/README.md)

### 跨受体搜索与亲和力配对验证

默认显式精修已修正搜索盒外的梯度投影。可选的 `--cubic_search_grid` 在相同网格上降低受体势的插值误差；六个人源靶点的首姿命中为三线性 32/66、三次 30/66。另用两个同法 Kd 系列、一个同法表观 Ki 系列和官方 Vina、AD4、柔性侧链做了配对验证。[研究结果、原始姿势及复现入口](benchmarks/cubic_grid/README.md)

`--soft_receptor_search` 在偶数编号搜索轨迹的受体网格中将排斥权重减半，用同一套 Late 精修和评分选取姿势。四种子下，凝血酶五个有实验结构的配体跨晶体首位命中从 20 增至 34/52，同系列 Kd 三倍对均为 17/37。其他人源靶点的结果随受体、预算和实验系列改变。[21 个靶点的实验、逐例证据与复现入口](docs/PAIRED_POSE_AFFINITY_20261126.md)

### 从给定姿势局部精修

研发版使用 `--local_only --receptor receptor.pdbqt` 时，Vina 和 Vinardo 以显式受体原子
优化给定配体姿势，并用同一表示计算最终分数。加 `--no_refine` 或只用 `--maps` 时，
优化和打分都使用网格。AD4 沿用其专用地图。旧研发构建的
`--local_only --receptor receptor.pdbqt` 曾用网格优化，再用显式受体打分。

```bash
python3 benchmarks/check_local_only.py build/refine-pool/vina /path/to/pre-fix/vina
```

### 固定搜索预算下的核心评分

有独立配体共晶结构时，`--template_stage contact_rank` 按其受体接触模式重排现有候选。亲和力另有两项研究性读数，`--contact_rotor_score` 检验接触约束的扭转代价，`score_physical_pool.py --mmff-affinity` 在冻结姿势上替换配体内能差。五靶点 Ki/Kd 和因子 Xa 外部系列给出了不同的适用范围。[实测对照、新增人源靶点和复算说明](benchmarks/frozen_core/README.md)

### FastDock Vina Pro 姿势重排

在上述命令加上 `--rerank_gauss2_scale 0.5`，可用一半宽高斯（`gauss_2`）权重重排姿势。
取值在 0 到 1 之间，也可写入 `--config` 指向的文件，命令行值优先。
`0` 去掉这一项，`1` 和省略参数一样，使用原 Vina 路径及输出格式。

重排使用原 Vina 搜索、聚类和显式受体精修得到的整批候选，在应用 `--energy_range` 前执行。
`--num_modes` 仍约束原搜索保留的候选数。引擎保持每个候选坐标，用原有截断距离、盒边界惩罚、
原子对查表和扭转修正重新计算分子间和配体内部能量，仅将用户配置的 `weight_gauss2` 乘以该比例。
受体空间索引同时供精修和重排使用。重排沿用 Vina 的共同内能参照，所有候选减去按 Pro
分子间能量与内部能量之和最低的姿势的配体内能，再计算报告分。报告分相等时沿用原 Vina 池内顺序。

Pro 的每个 PDBQT 模型包含三条标记。

- `REMARK VINA PRO RESULT` 依次给出 Pro 分数、相对 Pro 首姿势的 RMSD 下界和上界。
- `REMARK VINA NATIVE` 依次给出原 Vina 分数和该姿势在原候选池中的名次，从 1 开始。
- `REMARK VINA PRO GAUSS2 SCALE` 给出实际采用的缩放比例。

`MODEL` 编号和 `--energy_range` 都按 Pro 顺序解释，`INTER`、`INTRA`、`INTER + INTRA` 和 `UNBOUND`
给出 Pro 能量分量。从已有对接模型拆出的配体可以重新对接。启用 Pro 重排、延后筛选或 Armijo 精修时，
PDBQT 只写本轮引擎评分，保留配体输入中的其他注释。

比例小于 1 时，命令行支持直接提供刚性受体的常规单配体 Vina 对接及逐配体 `--batch`。
命令行会拒绝载入或写出预计算网格、使用柔性受体、同时对接多个配体，以及 `--no_refine`、`--score_only`、
`--local_only`、`--randomize_only`、`--autobox` 和显式 `--unbound_energy`。
这是固定候选的打分消融，准确率需按具体终点验证。

可用仓库内置的 1IEP 输入，复查恒等输出、固定候选重排、能量窗口、姿势重新对接及参数约束。

```bash
python3 benchmarks/check_pro_rerank.py build/vina /path/to/pre-pro/vina
```

### 局部刚体积分、离子对几何与 pH 条件换算

通过 `--score_only --rigid_window_pool POOL` 可以在已生成的姿势上积分平移和转动。已知芳基脒与受体羧酸根配对时，`tools/make_ion_pair_anchors.py` 可生成只用于搜索的方向约束；有实测 pKa 和参考 pH 亲和力时，`tools/ionization.py` 可计算条件换算。适用范围、反例与同一候选池中的实测结果见[局部运动、质子化和晶体环境研究](docs/LOCAL_MOTION_AND_LINKAGE.md)。

在五个 hCA II 磺酰胺的固定候选池上，给可动非环羧酰胺的非平面扭转加上能量代价，将联合首姿命中从 0/20 提到 6/20。独立运行的等预算候选池从 2/20 提到 8/20。此项是可选的姿势重排，跨体系亲和力效果随候选池而变。[酰胺扭转与复现实验](benchmarks/amide_search/README.md) · [胰蛋白酶和因子 Xa 的水相转移检验](benchmarks/trypsin_transfer/README.md)

使用 `--rank_amide SITES` 可在原生 Late 精修后按酰胺氮侧各碳取代基的扭转角重排，输出同时记录原生亲和力和几何罚分。新种子上，五个 hCA II 配体在两个受体中的 20 个合并候选池，首位命中从 1 个增加到 9 个。将同一项用于搜索后，首位命中回落至 1 个；同靶 Ki 排序的变化也因系列而异。位点文件生成、独立起始状态的共同参考、直接命令行输出的检验和备选晶体配体的对照见[酰胺搜索与姿势排序](benchmarks/amide_search/guidance.md)。

### 显式精修的 Armijo 回退线搜索

几何预条件的解析实现、对照官方 v1.2.7 的 CPU 预算实验，以及新增人源靶点的姿势和实测 Ki/Kd 结果见[跨靶点配对验证](benchmarks/geometry_comparison/README.md)。

在常规刚性单配体 Vina 对接命令加 `--refinement_optimizer armijo`，可在最终显式受体精修中
启用严格 Armijo 下降条件。默认值 `native` 沿用原优化器。Armijo 可与 `--pose_selection late`
和 `--rerank_gauss2_scale` 分别组合。`--refinement_optimizer geometry` 使用同一
Armijo 线搜索，按重原子位移预条件化旋转和扭转。两种选项都沿用原评分函数。
[度量、开销与实验设计](benchmarks/geometry_metric/README.md)

对给定姿势，可在显式刚性单配体 Vina 的 `--local_only` 路径使用这两种优化器。
局部精修还支持 `--refinement_optimizer rigid_metric`。它在起始重原子质心换算平移和旋转步长，
沿用原配体扭转坐标。三种模式可作为同一显式评分目标的配对数值对照。
此时 `--max_evals` 为整次局部优化的梯度调用上限，包含初始梯度和盒恢复重试。
0 沿用原迭代和盒重试约束。`Local refinement` 日志给出实际调用数和最终优化目标。

```bash
./build/vina --receptor receptor.pdbqt --ligand pose.pdbqt --config box.txt \
  --local_only --refinement_optimizer geometry --max_evals 32 --cpu 1 --out refined.pdbqt
```

Armijo 只接收目标值和梯度均有限、满足充分下降且实际降低目标值的步长。最多尝试十个步长。
全部试步失败时，本次 BFGS 停在当前已接受点。收敛判断使用当前点梯度，只有有限的正曲率
才会用于缩放初始逆 Hessian。显式精修保留原有迭代预算和盒边界重试限制。

在引擎目录可重复数值回归和实际对接测试。若系统 Boost 位于非标准目录，
可在编译命令加上 `-I/path/to/boost/include`。

```bash
mkdir -p build/optimizer
c++ -O3 -DNDEBUG -std=c++14 -DBOOST_ERROR_CODE_HEADER_ONLY -Isrc/lib \
  benchmarks/check_armijo_numeric.cpp src/lib/quaternion.cpp src/lib/random.cpp \
  -o build/optimizer/check_armijo_numeric
build/optimizer/check_armijo_numeric
PYTHONDONTWRITEBYTECODE=1 python3 benchmarks/check_armijo_refinement.py build/optimizer/vina build/search/vina
```

### 延后筛选候选姿势

在常规刚性 Vina 对接命令加 `--pose_selection late`，可把跨轨迹筛选推迟到显式受体打分之后。
默认 `--pose_selection native` 沿用原有网格能量排名及跨轨迹 2 Å 合并。两种设置使用相同的
蒙特卡洛轨迹、每条轨迹的随机种子及局部优化预算，每条轨迹最多保存 `--num_modes` 个极小值。

延后筛选汇集每条轨迹已保存的极小值。它在原坐标上计算显式受体能与配体内能之和，用相同
索引和原评分权重按能量排序，然后依次保留满足 `--min_rmsd` 的姿势，默认至多送
`--num_modes` 个姿势进入原显式精修。精修后重新按能量和同一最小 RMSD 筛选。初筛分数相等时
沿用轨迹及轨迹内顺序，精修后分数相等时沿用送精修顺序。两轮距离按输入的化学重原子序号计算
坐标 RMSD，保留共同受体坐标系。氢和 `G0–G3`、`W` 等虚拟原子排除，`CG0` 等真实原子计入，
等效原子沿用输入序号对应。Monte Carlo 轨迹内、默认跨轨迹合并及输出表中 `l.b.` / `u.b.` 采用 Vina
原有非氢原子约定，含虚拟原子。
初筛不调用局部优化，精修沿用原有盒边界重试上限。普通详略级别会报告汇集、送精修、
实际精修调用及精修后保留数。

该选项支持直接提供刚性受体的单配体 Vina 对接及逐配体批处理。默认输出仍使用
`REMARK VINA RESULT` 和精修后首姿势的共同内能参照。和 `--rerank_gauss2_scale` 一起使用时，
后者在延后筛选得到的池上执行，模型的原版名次指该池的重排前名次。

可用 `--refine_candidates N` 单独提高末端显式精修的输入上限，`N` 至少为 `--num_modes`。
该参数用于原生精修的 Late Vina，可配合轨迹存档。模板、极性几何、宽高斯缩放和 Armijo 使用各自的
原有精修设置。每条轨迹仍各自保存最多 `--num_modes` 个姿势，精修后也全局保留最多这个数。
受体盒重试沿用每个候选最多五次，`Late pose selection` 的 `refine` 和
`refinement_calls` 分别报告实际输入数和含盒重试的调用总数。

```bash
python3 benchmarks/check_late_pose_selection.py build/vina /path/to/pre-pro/vina
```

构建 `vina_split` 后，可复验拆分配体评分头与带虚拟原子的 CLI 行为。第二个参数指修复前的对照构建。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 benchmarks/check_refinement_input_compatibility.py \
  build/verified/vina build/optimizer/vina build/verified/vina_split
```

[`check_late_pose_geometry.cpp`](benchmarks/check_late_pose_geometry.cpp) 验证 `G0–G3` / `W`、
`CG0` 与交换等效原子的固定索引处理。用引擎同一 Boost 配置将其链接到
`build/verified/obj/` 的 `coords`、`model`、`parse_pdbqt`、`quaternion`、`random` 和 `utils`。

已有直接锌配位证据时，可在对接后按三配位蛋白 Zn 的四面体空位与声明的配体原子筛选姿势。 [锌配位姿势选择](docs/ZINC_POSE_SELECTION.md)

有已对齐的外部共晶配体时，可在 Late 路径中分别检验模板参与排序、候选选择和显式精修的作用。 [模板引导对接](docs/POSE_TEMPLATE.md)

可在相同 Late 候选池上比较 N/O 极性接触的径向奖励、方向奖励和角度惩罚。 [极性几何重排](docs/POLAR_GEOMETRY.md)

极性埋藏评分在原生 Late 候选池上比较埋藏、方向缺陷和氢键补偿三种罚项，也可直接诊断输入姿势。 [极性埋藏重排](docs/POLAR_BURIAL.md)

可用同构象的溶剂可及面积差、原生氢键的方向替换和极性位点补偿重排同一 Late 候选池。补偿以直接蛋白伙伴和独立水桥的最佳质量为限。 [极性溶剂面积与氢键方向重排](docs/NET_POLAR.md) · [独立受体结构证据提取](docs/INDEPENDENT_RECEPTOR_EVIDENCE.md)

已保存的原生池还可比较配体应变和独立侧链构象。可导出 GNINA 的同坐标输入，随后核验输出并按 CNNscore 或它与原生分的合分排序。同一测法的化合物系列可用固定首姿计算快速分差，另以显式溶剂的复合物腿和水腿计算相对结合自由能。 [固定池重排与 GNINA 输入](docs/FIXED_PHYSICAL_POOL.md) · [OpenFE 双腿复算与运行](docs/EXPLICIT_RBFE.md)

可用轨迹存档比较 Monte Carlo 接受姿势的额外候选价值和跨轨迹支持排序。 [轨迹存档](docs/TRAJECTORY_ARCHIVE.md)

固定池可按 CCD 核验芳香环，以几何接触数重排原生近简并姿势。搜索阶段可试验保持质心的均匀刚体重定向，并报告实际梯度预算。 [芳香几何与固定池](docs/FIXED_PHYSICAL_POOL.md#芳香环几何与近简并排序) · [刚体重定向](docs/RIGID_REORIENTATION.md)

固定姿势重评分可保留准备态拓扑，也可从可信化学图指定单配体键表。 [固定姿势与配体拓扑](docs/FIXED_POSE_TOPOLOGY.md)

独立供体搜索盒可用于半数轨迹的配体质心初始放置。 [口袋质心的初始放置](docs/POCKET_START.md)

用 `--separate_hunt_chain` 可让搜索接受和保存触发使用同一软化网格目标。候选副本仍经原有的第二次网格精修。 [搜索与候选精修的接口](docs/HUNT_CHAIN.md)

原生候选筛选可用 `--use_min_rmsd_for_merge` 把命令给定的 `--min_rmsd` 用于跨轨迹合并。 [候选筛选和距离控制](docs/NATIVE_MERGE_RMSD.md)

可用 `--smooth_cutoff` 检验普通 Vina 势最后 1 Å 的截断，用 `--mixed_cutoff_search` 让半数轨迹在平滑势上搜索，再统一用原生势筛选和精修。四个靶点的配对实测尚未显示共同的姿势和排序收益。[算法、逐例结果与复核](docs/CUTOFF_SWITCH_RESEARCH.md)

`--search_anchor` 可将受体提供的配位空位或其他配体的共有骨架用于全局搜索。五个体系的配对试验检验了候选恢复、质子化叔胺的两种构型，以及相同实验条件下的 Ki/Kd 排序。[搜索引导与微态研究](docs/ANCHOR_AND_MICROSTATE_RESEARCH.md)

在结构先验可能只适用于部分配体时，可加 `--search_anchor_trajectories even` 或 `odd`，按从 1 开始的轨迹序号选择接受引导的那一半；其余轨迹在原生势上搜索，所有候选统一用 Late 原生精修与排序。默认 `all`。不同轨迹分配、正确与错误头部及实测 Ki/Kd 的结果见[等预算混合搜索实验](benchmarks/mixed_anchors/README.md)。

同一供体有多种合法原子对应时，搜索引导现支持在完整分配之间平滑切换。三体系的配对复算、独立烷基苯晶体系列及 28 个因子 Xa 配体给出了姿势与亲和力的适用范围。[多对应搜索、受体构象与亲和力实验](docs/CORRESPONDENCE_AND_AFFINITY.md)

## 引用

研究中使用时，请引用 Vina 的原始论文：

- Trott, O. & Olson, A. J. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* 31, 455–461 (2010). https://doi.org/10.1002/jcc.21334
- Eberhardt, J., Santos-Martins, D., Tillack, A. F. & Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* 61, 3891–3898 (2021). https://doi.org/10.1021/acs.jcim.1c00203

## Topic / Keywords

- **解决的问题**：AutoDock Vina 在 CPU 上做大规模分子对接和虚拟筛选时耗时长、内存占用高，而常见的加速方案依赖 GPU，或会改变计算结果。
- **实现方式**：在不改变双精度、搜索算法和累加顺序的前提下做 C++ 工程优化，包括共享评分查找表、缓存热循环指针、移动语义替代深拷贝、BFGS 矩阵 SIMD 化；用 25 项基准测试和逐字节输出比对来验证。
- **项目类别**：计算化学 · 分子对接 · 虚拟筛选 · 药物发现 · 高性能计算（HPC）· C++ 性能优化 · AutoDock Vina 分支
- **研发选项**：Pro 可选用不同的搜索、精修或重排方法。性能表记录常规版发布时的已验证条件。
- **中文关键词**：分子对接、AutoDock Vina 加速、虚拟筛选、药物发现、蛋白-配体对接、CPU 优化、SIMD、内存优化、可复现、结果一致
- **English keywords**: molecular docking, AutoDock Vina, faster Vina, Vina CPU optimization, virtual screening, drug discovery, protein-ligand docking, PDBQT, SIMD, BFGS, memory reduction, bit-identical results, reproducible benchmark, computational chemistry, HPC, C++

## TL;DR

> Pro 提供常规 FastDock 的 CPU 优化和可单独启用的对接实验。常规版发布时的 25 项对照几何平均加速 1.228×，该数据来自上表所述的历史环境。

## Link

- [FastDock Vina Pro](https://github.com/HankerWu/fastdock-vina/tree/fastdock-vina-pro)
- [常规 FastDock](https://github.com/holamochi/fastdock-vina)
- 上游 AutoDock Vina：https://github.com/ccsb-scripps/AutoDock-Vina
- 同类项目：
  - QuickVina 2（加速版 AutoDock Vina）：https://github.com/QVina/qvina
  - AutoDock-GPU（GPU 加速的 AutoDock）：https://github.com/forlilab/AutoDock-GPU
  - gnina（基于深度学习打分的分子对接）：https://github.com/gnina/gnina

## 开源许可

本项目遵循 [Apache License 2.0](LICENSE)，保留 AutoDock Vina 的版权与署名。源码来源见 [PROVENANCE](docs/PROVENANCE.md)。
