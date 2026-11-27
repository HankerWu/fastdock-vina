# FKBP12 表观 Ki 与逐候选来源

**坐标复算勘误。** 首次发布时，本目录的独立 AD4 检查对 Pro 输出也逐通道限幅。当前复算程序按真实引擎选择限幅顺序。历史表保留原值，运行新版所得 Pro 坐标复算分数有细小差异。原生发布分和原始姿势可各自核对。[修复及全池复核](../../affinity_conditioning/README.md)

Holt 等在同一重组人 FKBP12、10 °C、pH 7.8 的 PPIase 测定中给出同一系列化合物 4–11 的表观抑制常数。重新核读第 9927 页 Table II 和同页正文，顺次为 186±24、110±20、12±5、250±30、10±1、7±2、300–600、300±40 nM。七个点估计用于秩相关，八个配体按点值和区间最近边界至少三倍的差距给出 15 对有向关系，误差栏随[原始转录](../../cubic_grid/human_panel/affinity/ligands.csv)保留。此定义没有把 ± 报告误差折入三倍边界。两个受体各预测八个化合物，两份受体的第一发布分等权平均。此权重没有物种丰度解释。同一系列、同一种子只计一次。 [原论文](https://doi.org/10.1021/ja00075a008) · [表 II 影印本](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1993%20%20%28vol%20115%29/22%20%20%2809879-10470%29/09925-09938.pdf#page=3)

在固定 524287 种子、两受体的第 1 条轨迹第 1 步，复用[上轮冻结库](../../cavity_triplet/fkbp/)的 1,024 个 matched 候选。每个候选沿用原来的重原子顺序、质心、取向和内部坐标。AD4 粗排所用 A、C、N、OA 类型逐个 PDBQT serial 映射到真实预制地图，只以每个重原子的类型通道求和。正值用 10 kcal/mol 限幅，越界夹到面上并加每轴 1e6 kcal/mol/Å 罚。两臂只交换类型和坐标的对应，真实/置换各查询 514,048 次，17,644 个原子位置越界并均在原分母内。按受体格统计，共 502 个原子条目中，306 个实际改变类别。重原子类型以外的电荷、显式氢及额外的电势与去溶剂化图未用于粗排；下游恢复完整原 AD4 输入。实际三条轨迹、单核、每条 24,000 软上限、9 模式。每格只改变一个刚体放置，平移、旋转及配体转动单键在后续搜索可以变化。受体、蛋白质子化状态和配体微态固定。坐标单位 Å，经验分单位 kcal/mol。[输入 SHA-256](results/source_hashes.json) · [每个图的制备与哈希](../../cubic_grid/fkbp_ad4/results/metadata.json) · [每格筛选成本](results/rank_cost.csv)

|方法|1FKG 正确/15|1FKH 正确/15|等权正确/15|等权 ρ / τ-b|8/10 弱减强|9/11 弱减强|
|---|---:|---:|---:|---:|---:|---:|
|官方 1.2.7|9|9|8|0.214 / 0.143|−0.208|+1.672|
|当前 Pro|6|5|7|0.107 / 0.048|−1.096|+0.761|
|原形状|9|4|7|0.143 / 0.048|−1.454|+0.761|
|真实类型|6|4|5|−0.179 / −0.238|−1.882|+0.761|
|置换类型|7|5|5|−0.107 / −0.143|−1.882|+1.088|
|Pro 双预算|10|4|9|0.143 / 0.048|+0.015|+0.234|
|官方双预算|9|4|7|−0.071 / −0.048|−0.394|+0.582|

16 格中的 13 格两种类型挑中完全相同的候选 ID。8 号配体的 1FKG 和 1FKH 冻结库分别有 12、5 个 ≤3 Å 的候选、均无 ≤2 Å 候选。原形状的首个 ≤3 Å 名次为 4、6，两种原子场均将其提到第 2，两种场在各受体都把同一个错误取向排第一，送入 hunt 时为 10.889、10.931 Å，返回 11.461、11.448 Å。9 号两库没有 ≤3 Å 候选；1FKH 的置换首名为 3.552 Å，二次局优到 3.137 Å 并入池，第 53 步被逐出。第 8 步已出现新的 ≤2 Å 入池提议，最后发布首名为 1.632 Å。三种方法的注入 ID 在 48 池中各自存活至本轨迹终点 1、4、3 次，实际全轨迹合池后均为 0；最终每个近姿须另查来源。[逐格来源和排名](results/stages.csv) · [逐发布坐标](results/modes.csv)

|实际晶体格|原形状 Top‑1 / Oracle (Å)|真实类型|置换类型|
|---|---:|---:|---:|
|1FKG / 8|5.109 / 5.109|5.938 / 5.938|5.938 / 5.938|
|1FKG / 9|4.590 / 2.398|4.590 / 2.398|4.590 / 2.398|
|1FKH / 8|10.877 / 3.821|9.346 / 8.494|9.346 / 8.494|
|1FKH / 9|6.973 / 6.097|6.973 / 6.097|1.632 / 1.632|
|等格平均|6.887 / 4.356|6.712 / 5.732|5.376 / 4.616|

同受体或局部 Cα 对齐受体框架中按保留化学键型与手性的等价重原子映射计算 RMSD。按四格，三方法 >10 Å 分别为 1、0、0；全部七臂合计 112 个非空池、1008 个模式，指定手性保持，最大重键长度漂移 0.001345 Å。新臂输出中最近配体–受体重原子距离均为 2.466 Å。由写出坐标重新插值后，所有七臂全部模式相对发布分的最大偏差 0.00343 kcal/mol，七臂 × 三类受体汇总 × 15 对的方向与输出分完全相同。真实类型相对 Pro 的八个分子分差与重原子数及 TORSDOF 的 ρ 为 0.094、−0.629，相对置换为 0.378、−0.293。八分子的形式电荷均为零。分数只用于同测定表观强弱排序；本面板的实验和取样无法估计结合自由能误差。[按分子的分布和严重错姿](results/per_molecule.csv) · [全部分位](results/pose_distributions.csv) · [逐受体逐分子的分数](results/scores.csv) · [每对分数及实测倍数](results/pairs.csv) · [偏差](results/bias_contrasts.csv) · [最近接触](results/min_contacts.csv)

## 成本与复核

16 个真实类型池合计 1,153,048 次 hunt + 第二次局优目标/梯度求值、108,544,170 次下游 AD4 通道插值。置换为 1,153,034 次和 108,542,739 次。两个臂分别再加 514,048 次粗场查询。以冻结库为输入，从 Python 新进程启动、加载地图、恢复候选、排名、启动各自 16 个引擎进程到写出最后一个池，端到端真实/置换分别是 20.891/18.392 墙秒、20.842/18.368 CPU 秒。计时重放的 32 份完整池、坐标档、轨迹历史、干预日志皆与原运行逐字节一致。单独原跑的各 16 次对接进程合计 19.136/17.199 墙秒，含建图读取和观察器，其余 Python 工作单列。历史 Pro/形状/双预算 Pro 的目标评估分别为 1,152,989/1,152,990/2,304,840；下游通道插值 108,538,317/108,538,878/216,957,906；对接子进程为 18.374/16.516/29.835 墙秒及 18.282/16.462/29.703 CPU 秒。官方/官方双预算对应 1,153,375/2,305,623 次求值、108,575,646/217,031,541 次插值、30.235/50.673 墙秒、30.085/50.454 CPU 秒。全部 CPU 计进程用户态和内核态。时钟受本机并行负载影响。[本次端到端](results/end_to_end.csv) · [下游原子/通道数](results/downstream_field.csv) · [所有池耗时及求值](results/runs.csv)

主配对实验重用归档源构象。另以准备好的配体、受体 PDBQT、AD4 原生地图和配置从干净目录完整重跑两种类型臂，未给进程历史池、源构象或候选。每臂先做 16 次捕获，再建腔、拟合全部刚体变换、保留生成器的形状粗分、算类型分、进相同下游搜索并发布。真实类型实测 39.885 墙秒 / 39.809 CPU 秒，置换 40.003 / 39.946。每个引擎单线程，两种全流程同时启动。真实/置换合计局优 2,306,037 / 2,306,023 次，下游通道插值加粗筛插值合计 217,596,535 / 217,595,104，几何最近邻查询另计。验证了两份受体、16 源构象、全部有序候选、全部干预发布坐标与归档完全相同。单核双预算 Pro 的原始 29.835 墙秒对照更接近该原型的实际求值预算。这里的普通输入从已制备的 AD4 地图开始，七臂共用地图；早期化学制备和 AutoGrid 建图不计入任一臂。 [真实类型完整成本](results/ordinary-typed.json) · [置换完整成本](results/ordinary-perm.json)

首次构造归档库的分项成本如下。16 次源轨迹重放 17.508 墙秒 / 17.427 CPU 秒，两份口袋建点 0.051 / 0.051 秒，16 格三点拟合 0.097 / 0.097 秒；其后可供三臂共享。原几何粗排另耗 0.221 / 0.220 秒及 514,048 次几何最近邻查询。没有把已有库视为普通对接免费的先验。合池身份观察额外重放 48 池，合计 55.004 墙秒 / 54.766 CPU 秒。双预算官方上游计数副本 16 池和发行版输出逐字节一致，额外 49.917 墙秒 / 49.749 CPU 秒。[历史源轨迹](../../cavity_triplet/fkbp/results/capture.csv) · [各阶段重放](results/merge_cost.csv) · [双预算官方对照](results/official_high_replay.csv)

首次读取官方双预算计数时，进度条与计数行共用一行，行首带星导致日志解析失败。保留了该次第一格的发行版输出和两个日志（归档中 failed_parse/）。改为定位明确的计数标识后全组重跑。合池重放脚本初次启动时与历史 run 模块重名，在产生对接输出前报错；随后将本轮模块更名为 frozen_assay.py。冻结机制、阈值和分母均沿用事前方案。

## 复现

从仓库根目录运行，WORK、VERIFY、BUILD、BOOST_PREFIX 均为绝对路径。只复核归档时使用第一段。原始包含两场臂的 16,384 组分数（每组各两分）、按序 16,384 候选逐例评价、实际原生输出池、注入、局优、逐出、合池坐标及全部阶段审核结果。官方普通、Pro、形状与双预算 Pro 的完整池通过 prepare 从已登记来源解压。已在空目录解包、逐文件核对 SHA 并重新生成分子标签、分数、几何和分位表。

```bash
D=$PWD/benchmarks/cavity_field_rank/fkbp
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/frozen_assay.py" prepare "$WORK"
(cd "$D" && sha256sum -c ARCHIVE.sha256)
tar xf "$D/field_forks.tar.xz" -C "$WORK"
python "$D/check.py" "$WORK" "$VERIFY"
cmp "$VERIFY/assay.csv" "$D/results/assay.csv"
cmp "$VERIFY/stages.csv" "$D/results/stages.csv"
```

从头重算排名和执行搜索时，先执行 prepare，然后运行以下命令。TRIPLET 为上轮三点基准目录，UP_SOURCE 为官方 v1.2.7 与对照一致的[带计数源码及对象](../../pair_gradient/build_arms.py)，UP_OBJECTS 指向该构建的对象目录。OFFICIAL 为 Vina 1.2.7 官方发行版，SHA 在程序内核对。

```bash
TRIPLET=$PWD/benchmarks/cavity_triplet
sh "$TRIPLET/build_observer.sh" "$BUILD" --boost-prefix "$BOOST_PREFIX"
bash "$TRIPLET/build_field_counter.sh" "$BUILD/src" "$BUILD/bin/obj" "$BUILD/field" "$BOOST_PREFIX"
bash "$TRIPLET/build_field_counter.sh" "$UP_SOURCE" "$UP_OBJECTS" "$BUILD/upcount" "$BOOST_PREFIX" "$TRIPLET/upstream_field.patch"
bash "$D/build_merge.sh" "$BUILD/src" "$BUILD/bin/obj" "$BOOST_PREFIX" "$BUILD/merge"
python "$D/frozen_assay.py" rank "$WORK"
python "$D/frozen_assay.py" run "$WORK" --exe "$BUILD/field/vina"
python "$D/frozen_assay.py" official-high "$WORK" --exe "$OFFICIAL" --counter "$BUILD/upcount/vina"
python "$D/merge_replay.py" "$WORK" "$BUILD/merge/vina"
python "$D/timed_replay.py" "$WORK" "$BUILD/timed" "$BUILD/field/vina"
python "$D/ordinary_inputs.py" measure "$BUILD/ordinary-typed" typed "$BUILD/field/vina" --reference "$WORK"
python "$D/ordinary_inputs.py" measure "$BUILD/ordinary-perm" perm "$BUILD/field/vina" --reference "$WORK"
python "$D/check.py" "$WORK" "$VERIFY"
```
