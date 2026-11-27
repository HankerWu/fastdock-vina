# 整配体刚体覆盖在六个新姿势实例中的迁移

DYRK1A 的两种晶体作双向交叉对接，种子为 96107 和 96131。这两种结构曾用于仓库其他实验，未用于本轮筛选器设计。RBP4 的 A1120 羧酸根和中性酸分别对接到 1RBP，种子均为 96107。六格包含三个分子、两个蛋白家族。盒为接收晶体配体已知口袋的 24 Å 立方体，预测程序只使用接收者、盒、制备配体、边表、随机种子。目标晶体仅在离线评价中打开。输入来源、散列、评分器、预算见[事前协议](../PROTOCOL.md)、[结构与 CCD 来源](../../cubic_grid/human_panel/sources.csv)、[准备输入](results/inputs.json)及[二进制](results/binary_context.json)。

## 发布质量

逐个输出模式从 PDBQT 读取重原子，用口袋 14 Å 内同编号、同插入码、同残基名的配对 Cα 将查询晶体变换到接收者坐标系。晶体与预测配体按 CCD/SMILES 全化学图及立体化学枚举等价映射，RBP4 阴离子允许羧酸根等价氧交换。整配体在此坐标系中直接计算 RMSD。独立复算另用 SciPy 旋转拟合受体，并从沉积坐标和原发布文件重建每个模式误差。Best@5 取已发布前五，Oracle 取实际发布池的最小值，Regret 为首位减去池 Oracle。所有池均有输出。

| 方法 | 均值 | Q1 | 中位 | Q3 | P90 | P95 | Best@5均值 | Oracle均值 | Regret均值 | 首位 >5 / >10 Å |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 原边长门 | 1.208 | 0.549 | 0.767 | 1.983 | 2.387 | 2.387 | 1.208 | 1.208 | 0.000 | 0 / 0 |
| 盒内边长门 | 2.249 | 0.549 | 0.767 | 1.975 | 5.511 | 7.078 | 2.246 | 2.239 | 0.010 | 1 / 0 |
| 最远点覆盖 | 1.201 | 0.540 | 0.733 | 1.977 | 2.390 | 2.391 | 1.201 | 1.201 | 0.000 | 0 / 0 |
| 均匀抽样、等工均匀抽样 | 1.212 | 0.571 | 0.769 | 1.974 | 2.381 | 2.384 | 1.212 | 1.212 | 0.000 | 0 / 0 |
| Pro 25k | 4.356 | 0.548 | 0.754 | 8.622 | 11.887 | 12.210 | 4.169 | 4.169 | 0.188 | 2 / 2 |
| Pro 50k | 4.350 | 0.548 | 0.754 | 8.622 | 11.816 | 12.104 | 4.186 | 4.156 | 0.194 | 2 / 2 |
| 官方 Vina 1.2.7 25k | 5.281 | 0.715 | 3.076 | 7.840 | 12.177 | 13.942 | 4.047 | 3.884 | 1.398 | 3 / 1 |
| 官方 50k | 4.454 | 0.545 | 0.721 | 6.643 | 12.159 | 13.933 | 3.689 | 3.575 | 0.880 | 2 / 1 |
| 预登记双轨迹、最近 CPU 的 Pro | 3.183 | 0.545 | 0.738 | 1.974 | 8.336 | 11.312 | 2.627 | 2.627 | 0.556 | 1 / 1 |
| 追加实际成本对照 | 3.180 | 0.534 | 0.738 | 1.974 | 8.336 | 11.312 | 2.624 | 2.624 | 0.556 | 1 / 1 |

表中单位 Å。每种分子恰占两格，逐分子等权的均值与实例等权相同；分位数另按三个分子计算。DYRK1A 四格的普通 Pro、原边长门、覆盖、均匀抽样首位均值 0.591、0.619、0.606、0.627 Å。两个 A1120 微态分别平均为 11.887、2.387、2.390、2.381 Å；中性酸在盒内边长门一格为 8.646 Å。主协议两轨迹成本最近的 Pro 对中性酸给出 14.287 Å、Oracle 10.952 Å。该分子贡献了不同方法均值间的大部分差异。均匀抽样与等工均匀抽样六格的源态、分数、所有步骤、首位和发布池字节相同。[逐格](results/cells.csv) · [相对于原边长门的配对差](results/paired.csv) · [分子及实例分位数](results/distributions.csv) · [全部发布模式](results/modes.csv)

## 从候选到最终模式

最远点覆盖在六格都缩小了盒内全集到所选库的**最大**最近刚体距离。六格最大距离的平均值为原边长 9.409、盒内边长 9.473、覆盖 3.147、均匀 6.025 Å；中位最近距离的格均值为覆盖 2.130、均匀 1.947 Å。最坏点获得更多预算，典型点密度可下降。覆盖库的 6,144 个候选中，6,080 个至少有一对配体与受体重原子相距 <2 Å，原边长有 5,267 个，均匀抽样有 6,016 个。这里检查的是初始放置，局优可消除碰撞。闭盒布尔量对落在盒面上的坐标末位敏感；1 格的两条独立 SVD 算法路径在 1 个 ID 的判定上存在差异，覆盖统计使用实际预测的 SVD 可行域。另见[事前登记的 Horn 灵敏度检验](../SENSITIVITY.md)。[所有格的门前、门后、覆盖距离和近接](results/gates.csv)

覆盖的库内最近晶体误差在六格**全部升高**，均值由原边长的 1.835 升至 2.461 Å；均匀为 1.844 Å。六格原边长、盒内边长、覆盖、均匀库分别有 4、4、1、4 格含 ≤2 Å；近姿总数分别 21、23、1、10。五类库的原生场首名在每格唯一。覆盖首名平均 5.470 Å，经 hunt 后为 4.969 Å；原边长 4.612→4.279 Å，均匀 3.653→2.959 Å。全部 30 个注入 ID 经过二次局优并入池。最远点覆盖的 6 个 ID 均在第 2–6 步逐出。仅 3ANR–96131 的原边长与盒内边长两臂能将注入 ID 追到最终第七位。生成的近姿计数、粗排、hunt、二次局优和原 ID 的结局均可逐个复查。[每个库的每个候选](results/candidates.csv.gz) · [注入事件](results/rank_stages.csv) · [阶段坐标误差](results/stages.csv) · [逐步容器](results/residency.csv)

主比较共读出 528 个模式（包含均匀两种重复发布及成本对照），已指定手性中心反转为零，最大重键长度漂移 0.001373 Å，所有首位最近受体距离 ≥2.546 Å。宽能窗低位存在穿插，阴离子的官方第六位重原子接触为 0.930 Å。每个输入的官方和 Pro 解析器所见的有序邻接、类型及转动标记逐行一致。官方计数程序的全部 12 个池与原版发行程序逐字节相同。[化学检查](results/chemistry.csv) · [解析器证据](results/graphs.json)

## 计算费用和运行限制

每格 1 核单进程运行。计时从制备完毕的输入开始，包含 Python 启动、原生建场、在线源态、配准、筛选、形状查询、原生类型粗排、两条搜索轨迹、精修及输出。状态准备输入、二进制编译、离线晶体评价单列。关闭逐阶段快照的复测逐字节复现发布池及原生选择。原边长、盒内边长、覆盖、均匀、等工均匀与普通 Pro 的六格墙钟合计为 4.755、5.956、10.075、5.231、10.091、1.969 秒，整棵进程树 CPU 为 4.744、5.944、10.057、5.209、10.076、1.946 秒。覆盖较原边长用 2.12 倍 CPU，较普通 Pro 为 5.17 倍；与等工均匀相差 0.2%。覆盖与等工均匀各做 178,404 次三点配准、171,547,648 次完整配体均方位移比较、129,024 次形状查询、129,024 次额外类型场插值。六格覆盖与普通 Pro 的含梯度场查询为 6,304,667 与 6,304,220。原生 hunt、二次局优、显式精修求值分别为覆盖 275156、25071、1529，普通 Pro 280728、19478、1701。受体和化学微态均固定，全分子三点变换不改扭角，随后原生优化允许改变扭角。能量读数采用 Vina 的 kcal/mol 约定。体积和热力学参考态没有进入这个筛选器。[全部操作计数](results/cost_cells.csv) · [总耗时](results/cost_total.csv)

普通 Pro 的 max_evals 是每条轨迹提前终止条件。3ANR 的第一次对照放宽到 12,800,000，搜索仍在每条轨迹 14,385 步停止。首次运行把“未达到覆盖 CPU”当作异常，保留原始日志、八个实际池和回溯于归档的 failure/。该次未保存完整逐次时钟。修订驱动检测两条轨迹均在 cap 以下跑完，继续倍增已无法改变实际搜索。在六格，双轨迹 CPU 最近对照合计 8.735 秒，为覆盖的 86.9%；四个 DYRK1A 格都未夹住覆盖。按独立冻结的追加方案增加轨迹，四格均夹住目标，实际选出的轨迹数为 4、3、2、2（分别按两个 3ANR、两个 7AKL 格）。保留两酸原生双轨迹后，六格总 10.178 CPU 秒、10.246 墙秒；相对覆盖分别 +1.2%、+1.7%，逐格 CPU 比为 0.843–1.335。追加比较的收益主要来自 A⁻，中性 HA 仍是严重错姿。[追加决策与预测](results/extra_trajectory_cells.csv) · [追加原生计数](results/extra_trajectory_cost.csv)

归档保留快照预测 30 次、对照和官方计数重放分别 24/12 次、Pro 发布程序控制重放 12 次、交错实耗重放 36 次、双轨迹扩大预算 22 次与计数核验 6 次、追加轨迹 6 次与计数核验 2 次。各组累计墙钟/CPU 秒分别为 37.98/35.43、17.34/17.24、12.25/12.20、7.62/5.07、38.08/37.98、21.14/21.11、8.73/8.73、10.99/10.87、3.32/3.26。不同任务可能并行，总墙钟是逐进程之和。首次正式姿势评价耗 2.564 墙秒。evaluation-launch 和 evaluation-import 保存了环境中外部 time 工具缺失及早期评价程序模块重名的失败，两次均发生在计算评价表之前，最终评价使用 shell 计时与无冲突导入。

## 重放与重新预测

从仓库根目录运行，W 和 OUT 取新的绝对路径。输入的 SHA 以及归档内每个文件均可核验，官方和 Pro 化学图位于归档的 official_graph/。预测完成后，评价程序使用仓库已有的 PDB 与 CCD；第一次下载和生成 PDBQT 的成本未纳入本轮时间。

~~~bash
D=$PWD/benchmarks/rigid_pose_coverage/pose
mkdir -p "$W" "$OUT" "$W/tmp"
export TMPDIR="$W/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
(cd "$D" && sha256sum --check --status ARCHIVE.sha256)
tar xf "$D/raw.tar.xz" -C "$W"
(cd "$W" && sha256sum --check --status "$D/results/raw_files.sha256")
python "$D/evaluate.py" "$W" "$OUT"
python "$D/cost.py" "$W" "$OUT"
python "$D/extend_cost.py" evaluate "$W" --out "$OUT"
python "$D/recheck.py" "$W" "$OUT"
for f in cells modes chemistry rank_stages stages residency distributions paired cost_cells cost_total source_identity gates extra_trajectory_modes extra_trajectory_cost extra_trajectory_cells; do
  cmp "$OUT/$f.csv" "$D/results/$f.csv"
done
cmp <(gzip -dc "$OUT/candidates.csv.gz") <(gzip -dc "$D/results/candidates.csv.gz")
~~~

从普通输入重新预测，另建空 W。OBSERVER 使用上一轮在线 FIFO 观测器，PRO 为普通 Pro，PRO_COUNTER 为带额外场计数的未干预 Pro；OFFICIAL 和 OFFICIAL_COUNTER 是官方 1.2.7 及对应计数重放程序。环境与散列见 binary_context.json，构建见[父实验](../../fragment_anchor/pose/README.md)。按机器允许的 CPU 设置 CORE。

~~~bash
mkdir -p "$W/tmp" "$OUT"
export TMPDIR="$W/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$D/run_panel.py" prepare "$W" --core "$CORE"
python "$D/run_panel.py" controls "$W" --core "$CORE" --production "$PRO" --pro_counter "$PRO_COUNTER" --official "$OFFICIAL" --counter "$OFFICIAL_COUNTER"
python "$D/run_panel.py" dock "$W" --core "$CORE" --observer "$OBSERVER"
python "$D/run_panel.py" timing "$W" --core "$CORE" --observer "$OBSERVER" --production "$PRO"
python "$D/run_panel.py" ladder "$W" --core "$CORE" --pro_counter "$PRO_COUNTER" --production "$PRO"
python "$D/extend_cost.py" predict "$W" --core "$CORE" --pro_counter "$PRO_COUNTER" --production "$PRO"
python "$D/evaluate.py" "$W" "$OUT"
python "$D/cost.py" "$W" "$OUT"
python "$D/extend_cost.py" evaluate "$W" --out "$OUT"
python "$D/recheck.py" "$W" "$OUT"
~~~

时钟决定扩大预算的比较对象，重新计时可使被选预算发生变化。固定原始输出的重放使用归档的选定项。
