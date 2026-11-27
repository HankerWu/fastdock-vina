# 同坐标三种排序的十格姿势实验

轨迹 1 第 1 次变异后，每格三种排序看到逐位相同的 1,024 组质心、取向、扭角和编号，各让一个候选进入同一原生 hunt。shape 是原形状排序；true 和 perm 均调用同受体、同盒的 Vina XS 图。perm 将 PDBQT 重原子序号排序后的类型循环移动半个分子，BRD9 实际换了 7/11 个类型，A1120 两种微态均为 18/28，两个 EPAS1 配体为 20/21、22/26。晶体查询姿势在评价进程中使用。[事前协议](../PROTOCOL.md) · [逐原子图](results/graphs.csv)

## BRD9 排进了第几名

| 种子 | 排法 | 首个 ≤2 Å 名次 | 获选 ID（从 0 起） | 获选 → hunt（Å） | 原 ID 逐出步 | 发布首位 / Oracle（Å） |
|---|---|---:|---:|---:|---:|---:|
| 91813 | shape | 8 | 45 | 11.933 → 11.933 | 4 | 6.129 / 6.129 |
| 91813 | true | 2 | 359 | 3.332 → 3.201 | 2 | 3.080 / 3.080 |
| 91813 | perm | 3 | 55 | 4.004 → 3.912 | 2 | 3.019 / 2.581 |
| 91823 | shape | 34 | 426 | 12.296 → 12.296 | 2 | 6.130 / 5.229 |
| 91823 | true | 12 | 399 | 3.156 → 3.070 | 4 | 1.140 / 1.140 |
| 91823 | perm | 15 | 569 | 2.658 → 2.766 | 2 | 2.755 / 2.743 |

两场臂都避开了形状首名的约 12 Å 错姿。真实类型相对置换把首个 ≤2 Å 候选再提前 1、3 名，均未取得唯一的局优名额。91813 的 ID 324 是上轮有标签分叉证实可捕获的 1.224 Å 候选，本轮 true 排第二。91823 的首个近姿在形状排序为 ID 115，在两场排序为 ID 51。91823 的 true 注入 ID 在第 4 步逐出，第 40 步才出现新的 1.141 Å 入池提议，最终首位来自后续采样。[逐候选分数和排名](results/candidates.csv.gz) · [各阶段与 ID](results/rank_stages.csv)

## 外推与微态

EPAS1 的两个配体晶体曾用于另一项插值研究，本轮在冻结排序器之后测试。两个接收受体均属 R247E 构建。4GHI→5TBM 的 93241 库中，两场均选中库最近项 1.179 Å，hunt 后达 0.928 Å。原 ID 于第 2 步逐出，最终首位同为 1.012 Å，形状排序在这格也选中 ≤2 Å 提议。93257 的库最近仅 2.261 Å，两场首位变为 13.492 / 12.093 Å。逆向两库均没有 ≤2 Å 候选，两场各选同一项，最终首位均超过 15 Å。5TBM→4GHI–93241 的真实类型把池 Oracle 从形状的 13.534 降至 7.317 Å，首位却由 13.567 升至 16.222 Å。四格平均 Top‑1 为形状 10.689、真实 11.635、置换 11.285 Å。[逐实例](results/cells.csv)

A1120 两种微态都进入 1RBP。HA 两个种子两场均选中同一项。A⁻–91813 的形状与真实都选 ID 591，最终 Top‑1 / Oracle 为 8.653 / 8.645 Å，置换选 ID 64，最终为 2.371 / 2.371 Å。新种子 93241 的阴离子三种排序都选 ID 696，得到 2.371 / 2.315 Å。四份库均无 ≤2 Å 的初始姿势。A⁻ 评价保留羧酸根氧的共振等价，HA 保持羰氧与羟氧身份；同种子的初始内部构象可随微态改变。这项对照没有测量微态在溶液中的占比。[逐分子分布、尾部和空池](results/distributions.csv)

三种干预的 10 个初始 ID 均接受并入池，最终直接发布数为形状 1、真实 1、置换 2。Pro 解析器的邻接顺序、XS 类型和转动标志在五份配体的推断图与显式图间均逐项一致。官方 1.2.7 使用自身的解析与发布，计数构建的 20 个池与官方发行版逐字节一致。594 个发布模式无已定义中心反转，最大重键长度漂移为 0.001216 Å。宽能窗下 A⁻–93241 官方第七位有 0.918 Å 受体重原子接触和 +42.264 kcal/mol 高经验分。[化学几何检查](results/chemistry.csv) · [每次入池和逐出](results/residency.csv) · [逐模式分数](results/modes.csv)

## 成本

输入为刚性受体 PDBQT、含可动扭角的配体、边表和已知口袋的 24 Å 立方盒。配体整体平移用 Å，取向为单位四元数，原生扭角用 rad。形状粗分无量纲，Vina 场与经验读分采用 kcal/mol 约定。库内只有刚体运动，hunt 和二次局优允许扭角变化。单 CPU、两轨迹、每轨迹 25,000 或 50,000 名义求值，显式精修另计。没有计算受体、水或微态的平衡概率，没有自由能误差结果。

[两遍实现冷进程计时](results/cold_timing.csv)包含 Python 启动、完整普通搜索提源、建场、构造和写出候选，以及完整干预搜索。每臂承担独立采源成本。计时关闭阶段与分数快照，30 个池及 30 份 bank 文件逐字节复现。十格每种排序各做 233,472 次形状近邻查询，另有 49,130 次建腔最近邻和 8,398 次半径查询。true、perm 各做 233,472 次额外标量 XS 插值，后续搜索各做 11,409,715 与 11,409,864 次含梯度插值，独立采源还各需 11,407,365 次。[逐格搜索和局优计数](results/dock_cost.csv) · [官方及 Pro 对照](results/control_cost.csv) · [干净生产进程计时](results/production_timing.csv)

初次将 controls 和 dock 重叠启动，dock 核验尚未落盘的 5TBM–93241 Pro 对照时提前退出。在完整对照落盘后重新顺序运行全部 30 个干预。汇总使用完整轮，提前退出的部分未单独计时。一次独立重放显示批量四元数运算的末位舍入差可能改变后续 MC 轨迹。最终实现保留原标量运算顺序，逐库校验全精度文本哈希。

## 复算

在仓库根目录执行。W 使用全新绝对路径。压缩档含实际 prepared 输入、全精度 bank、变换、源构象、原始输出池、逐候选分数和容器快照。四个继承库逐数组等同上轮归档。新库按同一规则冻结。[输入与盒](results/inputs.json) · [冻结核验](results/freeze.json) · [逐文件 SHA-256](results/raw_files.sha256)

~~~bash
D=$PWD/benchmarks/cavity_field_rank/pose
T=$PWD/benchmarks/cavity_triplet
mkdir -p "$W/tmp" "$W/output"
export TMPDIR="$W/tmp" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
(cd "$D" && sha256sum --check ARCHIVE.sha256)
tar xf "$D/raw.tar.xz" -C "$W"
(cd "$W" && sha256sum --check --status "$D/results/raw_files.sha256")
cp "$D/results/"{inputs.json,control_cost.csv,dock_cost.csv} "$W/"
python "$D/evaluate.py" "$W" "$W/output"
for file in cells modes chemistry stages residency rank_stages distributions; do
  cmp "$D/results/$file.csv" "$W/output/$file.csv"
done
~~~

从普通输入重新预测另选空目录 W。OFFICIAL 为 control_binaries.json 中的官方 1.2.7 发行版。OBSERVER 支持逐候选评分、插值计数及轨迹记录，未开启新排序时可作 Pro 对照。图探针和输入源码见相邻的 rigid_graph_control 实验，BOOST_PREFIX 为 Boost 路径。

~~~bash
mkdir -p "$W/tmp"
export TMPDIR="$W/tmp" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
bash benchmarks/cavity_field_rank/engine/build.sh "$W/build" --boost-prefix "$BOOST_PREFIX" --jobs 1
python tools/build.py --output-dir "$W/pro" --boost-prefix "$BOOST_PREFIX" --jobs 1
python benchmarks/pair_gradient/build_arms.py "$W/upstream" --arms upstream-counts --boost-prefix "$BOOST_PREFIX" --jobs 1
UP="$W/upstream/upstream-counts"
bash "$T/build_field_counter.sh" "$UP" "$UP/build/obj" "$W/official-counter" "$BOOST_PREFIX" "$T/upstream_field.patch"
OBSERVER="$W/build/bin/vina"
python "$D/run.py" prepare "$W"
python "$D/run.py" freeze "$W" --observer "$OBSERVER"
python "$D/run.py" controls "$W" --pro "$OBSERVER" --official "$OFFICIAL" --counter "$W/official-counter/vina"
python "$D/run.py" dock "$W" --observer "$OBSERVER"
python "$D/evaluate.py" "$W" "$W/output"
python "$D/timing.py" "$W" "$OBSERVER" --production "$W/pro/vina"
~~~

图探针用 [read_graph.cpp](../../receptor_linkage/aurora/internal_loops/rigid_graph_control/read_graph.cpp)；按同目录说明与 Pro 对象文件链接后可执行下面的命令。

~~~bash
python "$D/check_graphs.py" "$W" "$PROBE" "$W/output"
~~~
