# PDBQT 声明的转轴必须进入共价键图

对人凝血因子 Xa 抑制剂 CHEMBL369480 的原始准备态，记录 BRANCH 23 24 指定了一根真实化学单键。运动树绕其旋转，历史几何成键程序却从共价邻接中删掉了它。这会让四对本应排除的近键重原子进入配体内部非键势。上一轮已校正该分子的转子归一化为 7；本轮修复的是内部作用对的成员资格。化学定义上的修正没有提高本系列合并两个受体后的 Ki 排序。[上一轮转轴计数](../terminal_rotor_score/README.md) · [源代码](../../src/lib/model.cpp)

## 机制和修复

parse_pdbqt.cpp 的 set_rotor() 将两个声明端点标为 DISTANCE_ROTOR；原 assign_bonds() 仍可能因距离筛选或第三原子遮挡而不创建对应的 bond。修复在几何推断后核对运动矩阵的声明边，仅补回缺边。几何推断成功的输入沿用原有邻接顺序。新边同时记录两个方向和输入键长，配体或可动侧链随后由同一个图生成类型与非键排除。该矩阵不包含刚性受体网格原子，因而不会把彼此靠近的配体、侧链、金属或其他分子误连起来。矩阵遍历对 n 个非网格原子耗时 O(n²)，无持久附加表。

常规 Vina 配体内项只保留相对距离可动且最短键距大于三的普通原子对。固定坐标下，移出的项及其对应导数一并消失；搜索中使用的软化与既有径向能量表均未修改。实测病例的配体内部非键表为 497→493 对，删除 1–24、11–25、11–26、12–24。原始准备态中该键长仍为 2.094 Å，其他 6 条未明确声明的键仍因严重畸变而漏掉。该组 27 个其他化合物均与其相应输入化学图一致，补边前后成键表和重原子非键表相同。[图的逐分子核对](results/graphs.csv) · [历史坐标、力场与测定来源](../frozen_core/strain_transfer_audit/README.md)

## 相同坐标与重新对接

在原始输入的四个历史首姿上读取同一组坐标。下表的自参考分数为 score_only 所用的当前构象内部参考，单位 kcal/mol；修改只改变这四姿的内部能量。

| 受体与种子 | 配体内能，父版 → 修复 | 配体–受体能，两臂相同 | 自参考分数，两臂相同 |
|---|---:|---:|---:|
| 1FJS，271828 | 2.092→2.032 | −13.587 | −9.642 |
| 1FJS，314159 | 2.423→2.487 | −13.570 | −9.629 |
| 2P16，271828 | 2.345→2.281 | −12.085 | −8.576 |
| 2P16，314159 | 1.954→1.990 | −12.244 | −8.689 |

随后对原始与已归档的正常准备态交叉运行父版、补边版、补边版 + 完整化学图，共 24 个搜索。完整图源自同分子正常几何的 Meeko 映射；先校验两准备态的 SMILES、序号、极性氢归属和 BRANCH 序列严格一致，再将图施加于原始坐标。全部四个正常几何条件下三臂返回的 PDBQT 逐字节一致。畸变输入中，补边和完整图使用的仍是同一组畸变坐标，结果不同说明其他六条缺边也确实影响能量。完整图无法纠正键长。[全部固定坐标反事实](results/fixed.csv) · [实际评分、计数与输出散列](results/cells.csv) · [原始搜索日志和姿势](raw.tar.xz)

参照同一研究中 28 个 FXa 精确 Ki，两个种子各自首姿先按受体取中位，再对受体等权平均。实验 Ki 相差至少三倍的 295 对中，有 21 对涉及 CHEMBL369480。其余 27 个化合物的归一化及势能图与上一轮相同，沿用该轮统一重算的首姿分数。因此只有病例相关的 21 对可能改变方向。

| 原始准备态 | 1FJS | 2P16 | 两受体平均 |
|---|---:|---:|---:|
| 父版 | 228/295 | 187/295 | 207/295 |
| 默认补回声明轴 | 225/295 | 187/295 | 207/295 |
| 同坐标 + 完整化学图 | 229/295 | 186/295 | 207/295 |

1FJS 因补边获救 1 对、损失 4 对；2P16 获救和损失各 2 对。跨受体平均没有任何一个三倍分子对改变方向，Spearman ρ=0.472、Kendall τ=0.340。两受体只是对接输入，实验未测定它们在溶液中的占比。单化合物的重复种子和涉及同一分子的分子对都不能当作独立靶点。该分子缺少此处可用的实验姿势，不能报告其 RMSD。[每对分子的方向与差值](results/incident_pairs.csv) · [全部分层相关与计数](results/assay.csv) · [原实验面板、微态与限制](../affinity_partition/README.md)

在原始输入的四个搜索格，父版与修正合计分别为 802353、802611 次评估；四次单线程墙钟合计约 5.04 和 5.04 秒。正常几何两臂均为 802363 次。每格只报一次墙钟。有意隔离的配体准备微基准每臂 10 轮、每轮 300 次内存解析，CPU 中位数在原始输入上由 106.33 变为 103.30 µs，在正常准备态上由 107.09 变为 112.50 µs。无法由此确认普遍的速度优势。[准备耗时](results/time.csv)

## 测试与复算

C++ 回归从正常配体、超出共价长度筛选范围的声明轴、超出原空间分桶半径的轴和带拉长轴的可动残基出发，分别用 XS、AD 类型验证双向邻接、保存的原始长度及无重复。只改变键表，输入坐标能否支持合理化学结构仍需另行核查。完整图和良好几何均可用时，完整图提供较强的输入依据。

从仓库根目录运行。需要 C++14、Boost、Python、RDKit、NumPy、SciPy。可选 BOOST_PREFIX 的用法见同目录 build_probe.sh。实验使用的上游及父版源代码和二进制散列在相邻研究的 binaries.json 中。

~~~bash
W=$(mktemp -d)
git worktree add --detach "$W/parent" dec913f
python "$W/parent/tools/build.py" --output-dir "$W/old"
python tools/build.py --output-dir "$W/new"
D=benchmarks/declared_rotor_bonds
G="$PWD/benchmarks/bond_distance/ligands/graph.cpp"
bash "$D/build_probe.sh" "$W/parent" "$W/old" "$G" "$W/old-graph"
bash "$D/build_probe.sh" . "$W/new" "$G" "$W/new-graph"
bash "$D/build_probe.sh" . "$W/new" "$PWD/$D/check.cpp" "$W/check"
"$W/check" "$D/flex.pdbqt"
python "$D/replay.py" "$W/run" "$W/old/vina" "$W/new/vina" "$W/old-graph" "$W/new-graph"
P="$PWD/benchmarks/bond_distance/preparation/probe.cpp"
bash "$D/build_probe.sh" "$W/parent" "$W/old" "$P" "$W/old-parser"
bash "$D/build_probe.sh" . "$W/new" "$P" "$W/new-parser"
python "$D/time_parsing.py" "$W/run" "$W/old-parser" "$W/new-parser"
~~~

C3 几何规则是单独的可证伪实验。四个具实验结构的环配体和一个开链反例的输入、上游对照与结果在[三元环比较](../cyclopropane_topology/README.md)。
