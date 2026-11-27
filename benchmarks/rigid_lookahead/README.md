# 在粗网格筛选前沿刚体受力方向走一步，能否保住正确结合姿势？

人 BRPF1 的 5EPS 配体提供了明确的开发反例。原扫描已采到近晶体位置，但按未松弛的受体相互作用能筛选时将其淘汰。我们把每个初始姿势在同一受体场上允许的一次刚体松弛放在候选竞争之前。5EPS 两个种子的真实结合盆地因此得到精修机会，最终分别排第七和第五。开发面板整体没有首位增益，新的人 RBP4 面板也没有 2 Å 命中。普通 Pro 继续作为默认对接方法。实验实现可单独运行。[固定协议](PROTOCOL.md) · [初始问题与之前的实验](../rigid_translation_scan/README.md) · [完整原子场算法](../../tools/rigid_lookahead.py)

## 能量与更新

对一个独立制备的配体构象，用重原子几何中心定义质心坐标，所有原子单位权重。令初始朝向下第 $i$ 个重原子的偏移为 $r_i$，其和为零。给定平移 $p$，粗分为 $E(p)=\sum_i V_{a_i}(p+r_i)$。$V$ 为 0.75 Å 原生 XS 类型图的三线性插值，单位 kcal/mol，数据及数字精度继承上一轮。计算笛卡尔梯度 $g_i$ 后定义

$$
I=\sum_i(\|r_i\|^2\mathbf{1}-r_i r_i^\mathsf{T}),\qquad
v=-\frac{1}{N}\sum_i g_i,\qquad
\omega=-I^{+}\sum_i r_i\times g_i .
$$

上式给出最小二乘意义下对负笛卡尔梯度的刚体投影。沿该方向，原子的一阶位移为 $v+\omega\times r_i$。剔除不可观察的转轴，用同一正因子使均方位移的平方根等于一个网格间距。将该因子乘入 $v,\omega$ 后，实际试探位置为 $p+v+\exp([\omega]_\times)r_i$。只接受留在图内且同一势能严格降低的候选。这个长度定义在切空间，有限旋转后的真实笛卡尔位移可能略小。一步粗筛随后接原分支显式受体、可动配体扭转的局部精修，同一配体按共同参考的 INTER+INTRA 排名。此计算操作只针对搜索，不估计自由能。[Vina 评分与精修定义](https://doi.org/10.1002/jcc.21334) · [刚体流形优化的既有研究](https://doi.org/10.1021/ct300272j)

raw 沿用相同的原始格点。sham 在相同可观察刚体切空间中抽取各向同性方向，先按重原子均方度量白化，再取同样的一阶步长。它与梯度臂在完全相同的父姿势上都计算基线梯度和一次试探能量，也都能保留原位。两者的随机源与父姿势表分开。每个受体的已知配体中心定义 24 Å 盒，查询配体的实验坐标只在预测结束后用于评估或单列的有标签干预。受体静止，粗搜索只用一个查询内部构象。

## 从刚体候选到输出

下表每格按“Top-1 / Top-5 / 全池”列出对称性校正重原子 RMSD ≤2 Å 的运行个数。六个开发方向和四个二级验证方向都是交叉对接；新结构 RBP4 的两个不同配体各做原位和互换受体预测。两颗种子均逐独立运行报告，同一结构的种子和受体间互换不作为独立蛋白样本。SHBG 按预先登记的缺失邻近原子规则排除。[逐单元 RMSD](results/cells.csv) · [基线](results/baseline_cells.csv) · [RBP4 实验及结构来源](transfer/README.md)

| 方法 | BRPF1、EPHX2 开发，12 次 | RBP4 新结构，8 次 | PGR、PPARα 二级验证，8 次 |
|---|---:|---:|---:|
| 原版 Vina，每轨迹 3,000 | 1 / 1 / 1 | 0 / 1 / 1 | 5 / 5 / 5 |
| 原版 Vina，每轨迹 12,000 | 1 / 3 / 3 | 1 / 2 / 2 | 5 / 6 / 6 |
| 当前实验起点 Pro，每轨迹 3,000 | 1 / 1 / 1 | 0 / 1 / 1 | 5 / 5 / 5 |
| 当前实验起点 Pro，每轨迹 12,000 | 1 / 3 / 3 | 2 / 3 / 3 | 4 / 5 / 5 |
| 原始格点，9 次精修 | 0 / 1 / 1 | 0 / 0 / 0 | 2 / 3 / 3 |
| 随机一步，9 次精修 | 2 / 3 / 3 | 0 / 0 / 0 | 2 / 3 / 4 |
| 梯度一步，9 次精修 | 0 / 2 / 3 | 0 / 0 / 0 | 3 / 5 / 5 |

开发数据上，梯度与随机方向各有三个池命中，其条件首位命中分别为 0/3 和 2/3。二级验证上，梯度三次首位命中均出现在 PGR；PPARα 虽有一次精修到 1.77 Å，首位仍约 12 Å。排除有修补原子距共晶配体 7.690 Å 的接收体 1SQN 后，九槽的原始、随机和梯度在剩余六次中的首位分别为 1、1 和 2 次。留出覆盖有限。[二级验证事前约定](BUDGET.md)

**相反效应来自同一个筛选步骤。** 5EPS 种子 43103 的同一父姿势由方向第 723 升到第一，随后精修至 1.641 Å。将没有经过试步的同一坐标人为放入显式精修，结果为 1.236 Å，表明该算法在此提供了进入精修的机会。3WKE 种子 43117 的原始近姿势本来是粗能量全局第一；旋转和平移后它的分数也降低，但另一远姿势从 37.605 降至 7.728 kcal/mol，把近姿势压到方向第七。新入选者最后为 6.571 Å，而被舍弃的父、子坐标分别能精修至 1.802、1.924 Å。另一个 5EPS 种子以 2.152 Å 的起点精修至 1.523 Å，说明一个粗位置略超 2 Å 也可能处于正确盆地。[逐步重放、逐原子变化及强制保留对照](diagnostics/README.md)

**改变计算预算。** 查看首轮结果后固定了两档追加实验。24 槽精修各方向的一名胜者。48 槽在每一方向保留原始胜者以及本臂排序最高的另一组不同坐标；原始组也得到 48 次精修。开发集梯度臂三档 Top-1 均为 0/12，Oracle 为 3、3、4；随机臂首位均为 2/12，Oracle 同为 3、3、4。48 槽确实将 3WKE 的 1.802 Å 解找回，最终列第六。二级验证扩大到 48 槽后，三臂全部得到 3 / 5 / 5。更宽的池可追平九槽梯度的数量指标，逐例入选的姿势和成功种子仍有不同。二级验证曾由其他历史机制使用，本轮在预算策略确定后才运行。[预算与实耗](results/summary.csv)

**单一配体构象的边界。** RBP4–A1120 独立构象与实验结构在所有刚体配准中的最小 RMSD 为 2.225 Å。即使利用晶体把它直接放到最佳刚体位置，在两个接收体上也仅精修到 2.539 与 2.741 Å。对视黄醇，起始刚体下界为 1.654 Å。在其原位受体中，同类标签辅助对照可精修至 0.989 Å，显式总能量也优于盲选首位；换成 A1120 受体则精修到 2.431 Å。该下界只限制精修前的单构象刚体搜索。A1120 的中性羧酸微态以及局部蛋白环的作用仍有待区分。[配准与精修记录](transfer/results/fit.csv)

## 成本和复核

在一个 CPU 线程上将每个方法、预算和面板单独启动，以空输出目录重建该面板的受体地图。时间包含 Python 导入、建图、类型提取、试探、显式精修、读写和 RMSD 评估。受体及配体使用已准备的公用输入。开发集中，9 槽原始、随机、梯度各自耗时 3.48、5.95、5.61 秒（各 12 次对接），三者精修梯度求值为 4,399、4,875、4,860。梯度与随机各有 437,376 个父姿势、16,320,960 次原子场插值和 108 次显式精修。梯度的 48 槽运行耗时 6.65 秒、26,014 次精修梯度求值。18 个独立进程的峰值 RSS 约 143–146 MiB。Pro 的开发集两档 12 个 CLI 进程耗时分别为 2.13、3.08 秒，共 152,263、582,748 次搜索加显式梯度求值。两类求值单次成本不同，Pro 每次命令重复初始化，两类计时也包括不同的文件及标签评估开销。BRPF1、EPHX2 两个既有受体准备分别耗时 1.90、6.37 秒，六个开发查询的配体准备为每个 0.010–0.201 秒。新 RBP4 两个结构连同配体的共同准备分别耗时 10.23 和 9.19 秒，另行计费。[既有准备时间](../hydration_transfer/results/preparation.csv) · [新结构准备时间](transfer/sources.json)。全部运行使用 CPU，未调用 MD、外部学习势或结合自由能方法。[分段时间与查值次数](results/cost.csv) · [CLI 两档成本](results/baseline_cost.csv)

若以首位命中数和整面板测得的墙钟作两个目标，在开发集上 3,000 档 Pro（2.13 秒，1/12）和随机九槽（5.95 秒，2/12）体现了所测成本与命中数之间的取舍。RBP4 的 12,000 档 Pro 用 2.49 秒得到 2/8 首位，九槽梯度用 4.53 秒得到 0/8。PGR、PPARα 的 3,000 档 Pro 用 1.63 秒得到 5/8，九槽梯度用 4.88 秒得到 3/8。每个面板都存在首位更多且用时更少的既有算法配置。若将开发集终点改为全池覆盖，原始 Pro 增至 12,000 档也能得到九槽梯度的 3/12。以上是当前实现、硬件与有限测点下的比较，单次计时差异不宜作微小速度优势的依据。

独立实现采用插值单元的面差分和完整原子雅可比矩阵的 SVD，在 158,400 个实际父姿势上核对受力投影。两份代码接受后能量的最大差小于 1.0×10⁻¹² kcal/mol。闭区间图边界、单点、共线配体及一般配体另行测试。用实际解析器复核了 5,076 个精修输入的 PDBQT 类型、键、转轴及带查表索引的相互作用配对，全部与相应原始输入相同。RBP4 的 717 个 RMSD 由另一实现从坐标重算。将终版二进制的普通 --local_only 入口用于三个结构面板全部 252 个九槽梯度起点，输出 PDBQT 与基线观察器的精修输出逐字节相同。[发布程序复跑结果](results/release_parity.json)[数值检查](diagnostics/results/numeric.json) · [PDBQT 拓扑复核](topology/README.md) · [新结构独立审计](transfer/results/audit.json)

## 复算

在仓库根目录执行。使用 Python、NumPy、SciPy 以及仓库原有评估器的 RDKit、Gemmi、Meeko 依赖，版本可见两个输入归档。需有 C++14 编译器与 Boost。若 Boost 在非标准位置，为两个构建命令均添加同一 --boost-prefix 参数。实验固定父版本，不将同轮另外修订的亲和力分母带入假设比较。UPSTREAM 指向仓库记录的原版 3c65c0b 二进制。所有受体及配体从已发布的输入归档重建。[输入来源](../hydration_transfer/README.md) · [新结构来源和版本](transfer/sources.json)

~~~bash
R=$PWD
D=$R/benchmarks/rigid_lookahead
W=$(mktemp -d)
mkdir -p "$W/dev" "$W/secondary"
git worktree add --detach "$W/base" 2e610c46
python "$W/base/tools/build.py" --output-dir "$W/build"
python "$D/../rigid_translation_scan/build_probe.py" "$W/build/obj" "$W/build/probe" --source-root "$W/base"
tar -xJf "$R/benchmarks/hydration_transfer/inputs.tar.xz" -C "$W/dev"
tar -xJf "$R/benchmarks/search_backtracking/pose/inputs.tar.xz" -C "$W/secondary"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/replay_budget.py" --dev "$W/dev" --secondary "$W/secondary" --work "$W/sweep" --probe "$W/build/probe" --vina "$W/build/vina"
python "$R/benchmarks/hydration_transfer/panel.py" "$W/dev" dock --native "$UPSTREAM" --pro "$W/build/vina"
python "$R/benchmarks/hydration_transfer/panel.py" "$W/dev" geometry
python "$R/benchmarks/rigid_translation_scan/baseline.py" --dev "$W/dev" --hold "$W/secondary" --work "$W/baselines" --pro "$W/build/vina" --upstream "$UPSTREAM"
python "$D/collect.py" "$W/sweep" "$W/results" --baselines "$W/dev" "$W/baselines"
python "$D/diagnostics/check.py" "$W/numeric.json"
(cd "$D/results" && sha256sum -c SHA256SUMS)
~~~

results 下六个按面板和预算命名的压缩包包含逐姿势输入、精修输出、完整日志、各臂逐候选及逐池排序。baselines.tar.xz 收录 80 次同条件上游和 Pro 的原始预测和日志。结果表可直接查看。各项不重复抽样构建置信区间，两个开发靶点不能支撑广泛外推。[新结构复现](transfer/README.md) · [固定起点的因果复现](diagnostics/README.md)

复核发布程序：完成三个面板并依照 transfer/README.md 准备好 INPUT_RBP4 和 RUN_RBP4 后，VINA_FIXED 指向包含本轮转子修复的可执行文件。

~~~bash
python "$D/verify_release.py" "$VINA_FIXED" "$W/dev" "$W/secondary" "$W/sweep" "$INPUT_RBP4" "$RUN_RBP4" "$W/release_parity.json"
~~~
