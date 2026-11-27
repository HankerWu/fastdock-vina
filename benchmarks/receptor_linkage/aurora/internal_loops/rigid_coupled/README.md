# 固定一个配体重原子的刚体提议能否提高完整图盲搜覆盖？

三种配体、两枚新种子的六个配对条件中，将普通平移改为等原子位移的单枢轴转动，平均发布池 Oracle 从 8.162 升至 8.390 Å。AT9283 两次搜索略有改善，danusertib 的一个种子却在整个搜索中偏离实验位置。保留方案为独立研究补丁。下表给出在预先规定预算下的全部预测。[事前协议](PROTOCOL.md) · [完整模式、分数与重算](results/modes.csv) · [分项求值](results/counts.csv)

| 配体 | 种子 | 官方 25k | 官方 50k | Pro 25k | Pro 50k | 随机数匹配平移 | 单枢轴 |
|---|---:|---:|---:|---:|---:|---:|---:|
| AT9283 | 81371 | 8.981 / 7.168 | 9.110 / 7.168 | 8.998 / 7.693 | 8.998 / 7.693 | 13.197 / 7.504 | 9.217 / 7.377 |
| AT9283 | 81373 | 8.987 / 7.328 | 8.987 / 2.701 | 8.931 / 8.768 | 8.931 / 6.743 | 12.765 / 9.351 | 8.887 / 8.548 |
| danusertib | 81371 | 12.489 / 11.047 | 5.585 / 3.672 | 8.262 / 7.499 | 6.052 / 6.052 | 8.209 / 8.209 | 8.548 / 7.534 |
| danusertib | 81373 | 10.013 / 8.005 | 8.786 / 8.005 | 14.167 / 12.680 | 14.167 / 12.680 | 7.481 / 7.481 | 14.638 / 14.495 |
| BRD9 82I | 81371 | 8.097 / 6.581 | 6.110 / 6.110 | 8.103 / 5.390 | 8.103 / 5.390 | 8.676 / 7.354 | 8.110 / 6.188 |
| BRD9 82I | 81373 | 3.110 / 3.110 | 3.110 / 3.110 | 8.121 / 6.941 | 1.100 / 1.100 | 12.079 / 7.527 | 8.129 / 6.200 |

表格为受体坐标系的首位 RMSD / 实际 9 模式池的最小 RMSD（Å），配体化学等价映射取最小，不再作配体叠合。AT9283、danusertib 分别取 2W1G:A、2J50:A 配体，转入 2X6E 框架，在相同 1OL7 撤约束双环上对接；BRD9 查询为 6YQW:A 的 82I，投到 6Y7L 的 A 链口袋框架。各臂使用同一受体 PDBQT 与原始独立制备的配体 PDBQT。Pro 的四臂均用相同完整 SDF/CCD 边表和 late 发布策略。官方 1.2.7 没有显式键表选项，沿用自己的推断拓扑与发布策略，只作外部软件基线。冻结的配对算法归因在 Pro、随机数匹配平移和单枢轴间完成。对应的机理对照无额外受体、共晶模板、锚点或测试标签。六次预测来自三个分子和两个蛋白家族。小面板的尾部分位数受个别观测影响较大。[结构、制备与哈希](results/inputs.json) · [Aurora 原始资料与解包](../REPRODUCE.md) · [BRD9 来源](../../../../torsion_proposals/README.md)

| 方法（6 次） | Top-1 均值 | Q1 / 中位 / Q3 | P90 / P95 | 平均 Best@5 / Oracle / Regret | >5 / >10 Å | 空池 | 求值 | 墙钟 / CPU 秒 |
|---|---:|---|---|---|---|---:|---:|---|
| 官方 25k | 8.613 | 8.318 / 8.984 / 9.757 | 11.251 / 11.870 | 7.223 / 7.206 / 1.407 | 5 / 2 | 0 | 303262 | 8.978 / 8.939 |
| 官方 50k | 6.948 | 5.716 / 7.448 / 8.937 | 9.049 / 9.079 | 5.128 / 5.128 / 1.820 | 5 / 0 | 0 | 602727 | 8.089 / 8.033 |
| Pro 25k | 9.431 | 8.157 / 8.596 / 8.981 | 11.583 / 12.875 | 8.441 / 8.162 / 1.269 | 6 / 1 | 0 | 301535 | 1.677 / 1.661 |
| Pro 50k | 7.892 | 6.565 / 8.517 / 8.981 | 11.583 / 12.875 | 7.219 / 6.610 / 1.282 | 5 / 1 | 0 | 601428 | 2.786 / 2.746 |
| 随机数匹配平移 | 10.401 | 8.326 / 10.378 / 12.593 | 12.981 / 13.089 | 8.546 / 7.904 / 2.497 | 6 / 3 | 0 | 301380 | 1.587 / 1.581 |
| 单枢轴 | 9.588 | 8.234 / 8.718 / 9.134 | 11.927 / 13.282 | 8.410 / 8.390 / 1.198 | 6 / 1 | 0 | 301301 | 1.736 / 1.734 |

逐池 Regret 定义为 Top-1 − Oracle。名义预算是每轨迹的软上限，全部使用一个 CPU、两条轨迹。求值总数包含 hunt、二次局优与显式精修的实际目标和梯度联合调用；late 阶段的额外纯能量粗排另记在 late_prescore。程序时间包括启动、建图和输出。计时为各独立进程一次实测，不能从官方两档的非单调用时推断双预算更快。12 个官方池与上游源码计数版输出逐字节相同，计数版复跑另耗 15.526 秒墙钟与 15.496 秒 CPU。全部 324 个发布姿势重原子键长相对输入最大变化 0.001287 Å，108 个 danusertib 手性检查均保留，受体–配体最近重原子距离最小 2.405 Å。另一实现从 PDB 直接配准受体并用 RDKit 复核 RMSD，最大差 3.02×10⁻¹⁴ Å。[逐姿化学检查](results/qa.csv) · [官方计数验证](results/counter_replay.csv) · [逐池分位和成本](results/distributions.csv)

## 有限刚体移动和实际分叉

原生平移从局部球内生成扰动 $d u$，其中 $u$ 为单位方向，$d$ 以 Å 计。每次选择平移槽时，在当前配体的真实重原子中均匀抽取枢轴 $p$，对所有原子作 $x'_i=p+R_u(\theta)(x_i-p)$。若 $s^2=N^{-1}\sum_i\|u\times(x_i-p)\|^2$，本次转动的重原子 RMS 位移为 $2s\sin(\theta/2)$。以 $\theta=2\arcsin(\min(1,d/(2s)))$ 匹配原生位移，轴矩为零时退回平移，角单位弧度。该选择只使用当前构象与原随机流。保留原转动槽与全部扭转槽，后续局优仍可移动所有原子。单枢轴瞬间固定的是该原子的空间位置；与受体的单个距离不变，并不保证氢键方向或真实接触正确。随机数匹配平移组也抽取同一随机重原子，然后对所有原子施加原平移。两臂在首次不同提议之后可能在 MC 判定中消耗不同随机数，不能将后续同一步号解释为共同状态。[变异与几何实现](implementation.patch)

独立运行两个带只读观察器的二进制，从实际接受态开始逐事件检查。每个配体、种子和轨迹各有一对真实首次分叉（12 对）。首次不同提议之前，记录到的共同轨迹前缀逐字节相同。两臂的重原子 RMS 步长逐对相同，范围 1.321–1.996 Å、均值 1.788 Å。单枢轴质心移动 0.782–1.656 Å，均值 1.239 Å。实测最固定原子的最大偏移 2.87×10⁻¹⁴ Å，配体任意两原子间距离最大漂移 2.49×10⁻¹⁴ Å。该规则把一部分原有的质心位移消耗在角度变化上。这些共同起点离晶体 7.873–19.532 Å，不能据此估计近晶体盆地内的捕获概率。[每对共同前缀的分叉及接受信息](results/splits.csv) · [全部见证原始坐标](results/witnesses.npz)

| 配体 | 种子 | 匹配平移，提议 / hunt 返回 / 发布 | 单枢轴，提议 / hunt 返回 / 发布 |
|---|---:|---|---|
| AT9283 | 81371 | 6.809 / 7.261 / 7.504 | 5.996 / 6.335 / 7.377 |
| AT9283 | 81373 | 8.694 / 8.542 / 9.351 | 7.579 / 7.659 / 8.548 |
| danusertib | 81371 | 7.074 / 7.074 / 8.209 | 6.585 / 7.222 / 7.534 |
| danusertib | 81373 | 6.515 / 6.515 / 7.481 | 13.221 / 13.398 / 14.495 |
| 82I | 81371 | 5.867 / 7.290 / 7.354 | 5.470 / 6.141 / 6.188 |
| 82I | 81373 | 5.874 / 6.660 / 7.527 | 6.163 / 6.368 / 6.200 |

各列为相应阶段全部坐标的最小 RMSD，阶段极值可能来自不同步。13,698 个枢轴臂提议包含该臂的所有原生旋转、扭转与枢轴移动，12,559 个匹配平移臂提议也包含全部槽。没有任何原始提议、hunt 返回或发布姿势进入 3 Å；全体局优求值中，82I–81371 枢轴臂的两次回溯试探进入 3 Å，最近为轨迹 2 第 842 步第 14 次求值，2.914 Å，未获线搜索采纳。近晶体候选的缺失出现在最终排序之前。AT9283 的两个种子上，枢轴相对匹配平移均更接近晶体；danusertib–81373 的退化已发生在提议阶段。这个证据解释本轮观察到的路径，无法将质心位移缩短认作所有体系的唯一原因。[完整阶段计数及见证定位](results/stages.csv) · [646,824 条事件与误差](results/events) · [起点干预与捕获范围](../rigid_capture)

## 对既有机制、能量和成本的解释

既往刚体扫描中，5EPS 已提出近姿后被方向内粗排过滤，3WKE 则可在最终排序中失去首位。受体场一步梯度 lookahead 能改变进入显式精修的机会，也会优先稳定错误盆地。180° 主轴翻转保持质心，曾在部分 ESR1 例上改善取向。口袋起点实验直接规定部分轨迹的初始质心范围，供体共晶锚点实验则额外读取其他化合物的配体位置。当前实验针对另一环节，检验能否在真实 MC 中轻微转动已有内部构象并保留一个空间位置。DOCK 6 的 anchor-and-grow 先在受体球上定位一个刚性片段，再逐层生长，也需要解决最初锚点定位。这里均匀取一个尚未证实正确的原子，没有借用其他配体的位置或 DOCK 的片段搜索。[刚体扫描](../../../../rigid_translation_scan/README.md) · [受体场预搜索](../../../../rigid_lookahead/README.md) · [主轴翻转](../../../../rigid_flip/README.md) · [口袋起点核验](../../../../check_pocket_start.py) · [非查询模板转移](../../../../anchor_search/transfer_holdout.md) · [DOCK 6 原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4469538/) · [官方实现](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm)

受体、化学微态与水处理固定。显式自由度是配体刚体的平移、取向及 PDBQT 定义的扭角。受体侧链、环运动、溶剂/离子和质子交换都未采样；干态 Vina 分仅经验性地平均部分环境影响。MC 接受尺度 1.2 是经验分尺度。一次变异和 BFGS、截断保存与最终排序不能产生平衡占比。没有选用物理温度、溶液化学势、游离参考或实施 1 M 标准体积转换；本研究未测量绝对或相对结合自由能误差。[Vina 原始研究](https://doi.org/10.1002/jcc.21334) · [共同系综与标准态要求](../../../../../docs/BINDING_ENSEMBLE_AND_DOCKING.md)

36 个主进程总计 2,411,633 次联合求值、24.852 秒墙钟和 24.695 秒 CPU。12 个带观察器的重跑另耗 602,681 次求值、3.956 秒墙钟、3.918 秒 CPU；首轮观察器验证另用 3.870 秒墙钟、3.859 秒 CPU，全部 278 份数组和 203 个阶段的记录两轮一致。每轮的 12 个发布文件和各计数均与相应未插桩二进制相同。逐梯度双精度完整输出合计 401,166,392 字节，按 SHA 和下列命令可重建。归档另保存所有事件标签、步、能量、RMSD，203 份阶段最优坐标和 72 份首次分叉的原始父/子/返回坐标。独立 RDKit 查验这两类坐标，最大差分别为 3.02×10⁻¹⁴ 和 2.85×10⁻¹⁴ Å。编译、解压、历史制备与离线评价未计入对接时间；环境起初临时盘满，切换工作区临时目录后全量重编成功，墙钟 26.500 秒。派生两个变异二进制另用 4.044 秒墙钟。[构建记录](results/build.json) · [原始流的哈希](results/observation_cost.csv) · [公开坐标见证的独立核查](results/witness_check.json)

## 复现

在小端 Linux 上，以 WORK 指向按上级 REPRODUCE.md 解开且含原始 PDB、ligand PDBQT、键表与已补环受体的目录。T 取新的绝对路径；OFFICIAL 是 AutoDock Vina 1.2.7 官方 Linux 可执行文件；COUNTER 是按 benchmarks/pair_gradient/build_arms.py 的 upstream-counts 从上游 3c65c0b3 编译的计数版本。BOOST_PREFIX 为本机 Boost 前缀（系统安装可省略选项）。Python 包与版本见 build.json。

```sh
set -eu
R=$(git rev-parse --show-toplevel)
D="$R/benchmarks/receptor_linkage/aurora/internal_loops/rigid_coupled"
mkdir -p "$T/src" "$T/brd" "$T/tmp" "$T/build" "$T/check"
export TMPDIR="$T/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
git -C "$R" archive 8eafd0c src tools/build.py | tar xf - -C "$T/src"
tar xf "$R/benchmarks/torsion_proposals/fixtures/inputs.tar.xz" -C "$T/brd"
python "$T/src/tools/build.py" --boost-prefix "$BOOST_PREFIX" --output-dir "$T/build/base" >"$T/base.log"
python "$D/build_variants.py" "$T/base.log" "$T/build"
python "$D/study.py" run "$WORK" "$T/brd" "$T/run" --build "$T/build" --official "$OFFICIAL" --counter "$COUNTER"
cmp "$T/run/modes.csv" "$D/results/modes.csv"
cmp "$T/run/qa.csv" "$D/results/qa.csv"

# 仅复核已归档的发布坐标，无需二进制。
tar xf "$D/raw.tar.xz" -C "$T/check"
cp "$D/results/jobs.csv" "$D/results/counts.csv" "$T/check/"
python "$D/study.py" verify "$WORK" "$T/brd" "$T/check"
python "$D/verify_witnesses.py" "$WORK" "$T/brd" "$D/results"

# 重新生成所有逐事件二进制流、归档表与首次真实分叉。
cp -R "$T/src" "$T/observe-src"
cp "$D/../basin_trace/observer.h" "$T/observe-src/src/lib/basin_observer.h"
(cd "$T/observe-src" && GIT_CEILING_DIRECTORIES="$T" git apply "$D/../basin_trace/instrument.patch")
python "$T/observe-src/tools/build.py" --boost-prefix "$BOOST_PREFIX" --output-dir "$T/observe-base" >"$T/observe.log"
python "$D/build_variants.py" "$T/observe.log" "$T/observe-build"
python "$D/branches.py" "$WORK" "$T/brd" "$T/run" "$T/observe-build" "$T/trace" --keep-traces
python "$D/archive_events.py" "$WORK" "$T/brd" "$T/trace" "$T/events"
cmp "$T/trace/stages.csv" "$D/results/stages.csv"
cmp "$T/trace/splits.csv" "$D/results/splits.csv"
```
