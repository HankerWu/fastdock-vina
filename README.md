# FastDock Vina

基于 [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) 的 CPU 优化版本，保留双精度计算和原有搜索设置，不使用 GPU。

## 相比原版的效果

下表记录发布版的历史结果。后续改动和复测见[工作线程复用](docs/BENCHMARKS.md#工作线程复用的本地复测)
和[受体查表与空间索引](docs/BENCHMARKS.md#受体查表与空间索引复用)。
近期对照见[刚性受体共享和建图复用](docs/BENCHMARKS.md#刚性受体共享和建图复用)、
[按使用需求生成网格](docs/BENCHMARKS.md#按使用需求生成网格)、
[按配体类型准备受体评分表](docs/BENCHMARKS.md#按配体类型准备受体评分表)、
[亲和力网格读取与姿势写出](docs/BENCHMARKS.md#亲和力网格读取与姿势写出)
和[局部搜索向量与网格索引](docs/BENCHMARKS.md#局部搜索向量与网格索引)。

在 Intel Xeon Platinum 8362 上测试 19 个完整对接案例，每个案例重复 3 次，取各案例耗时中位数后求和：

| 版本 | 总耗时 | 耗时降低 |
| --- | ---: | ---: |
| 原始 AutoDock Vina | 3242.525 秒 | — |
| FastDock Vina | 2989.231 秒 | **7.81%** |

全部 **25 项测试**（19 项完整对接 + 6 项评分／局部优化）的**等权几何平均加速比为 1.228×**。计算时，每项先取 3 次耗时中位数，以 `原版耗时 / FastDock 耗时` 得到加速比，再对 25 个比值取几何平均（[逐项数据](benchmarks/benchmark.csv)）。

两个指标的口径不同：**1.228×** 对每项测试赋予相同权重；**7.81%** 是 19 项完整对接的合计耗时降低，折合约 **1.085×**。几何平均值不代表每个案例的最低加速或单例最高加速。

19 个案例均快于原版，峰值内存（RSS）降低 **28.71%–73.27%**。已测配置的对接 PDBQT 输出与原版逐字节一致；结果仅适用于已测环境和输入。详细数据与验证范围见 [实验记录](docs/BENCHMARKS.md)。

![25 项测试的加速比及汇总指标](docs/figures/speedup-summary.png)

## 优化方法

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

需要 Python 3.8+、支持 C++14 的 GCC/Clang 和 Boost（已验证 1.82）。

```bash
git clone https://github.com/Drangon-Li/fastdock-vina.git
cd fastdock-vina
```

Ubuntu / Debian：

```bash
sudo apt-get install build-essential python3 libboost-all-dev
python3 tools/build.py --jobs 4
```

macOS（已安装 Xcode Command Line Tools 和 Homebrew）：

```bash
brew install boost
python3 tools/build.py --boost-prefix "$(brew --prefix boost)" --jobs 4
```

## 运行

使用仓库自带的 1iep 示例：

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

结果保存至 `demo_out/1iep_out.pdbqt`。更多参数见 `./build/vina --help`。示例使用 `exhaustiveness=8`，上述性能测试使用 `32`。

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

### 从给定姿势局部精修

研发版使用 `--local_only --receptor receptor.pdbqt` 时，Vina 和 Vinardo 以显式受体原子
优化给定配体姿势，并用同一表示计算最终分数。加 `--no_refine` 或只用 `--maps` 时，
优化和打分都使用网格。AD4 沿用其专用地图。旧研发构建的
`--local_only --receptor receptor.pdbqt` 曾用网格优化，再用显式受体打分。

```bash
python3 benchmarks/check_local_only.py build/refine-pool/vina /path/to/pre-fix/vina
```

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

### 显式精修的 Armijo 回退线搜索

在常规刚性单配体 Vina 对接命令加 `--refinement_optimizer armijo`，可在最终显式受体精修中
启用严格 Armijo 下降条件。默认值 `native` 沿用原优化器。Armijo 可与 `--pose_selection late`
和 `--rerank_gauss2_scale` 分别组合，Monte Carlo 搜索、单独局部优化和其他评分模式沿用原路由。

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

## 对接构象

![1iep 配体在受体中的位置及原版与 FastDock 的构象叠加](docs/figures/docking-overlay.png)

左图展示配体在受体中的位置，右图放大 MODEL 1：原版为蓝色，FastDock 为橙色，两者输出坐标完全重合（40 个原子的 RMSD 为 **0.000 Å**，未做结构对齐）。图中保留真实化学键，显示 37 个重原子。[对接输出](examples/basic/1iep_docked.pdbqt) · [图表数据与复现方法](docs/FIGURES.md)。

## 许可与致谢

遵循 [Apache 2.0](LICENSE)，保留 AutoDock Vina 的版权与署名。研究使用请引用 [Vina 原始论文](https://doi.org/10.1002/jcc.21334)和 [Vina 1.2.0 论文](https://doi.org/10.1021/acs.jcim.1c00203)。源码来源见 [PROVENANCE](docs/PROVENANCE.md)。
