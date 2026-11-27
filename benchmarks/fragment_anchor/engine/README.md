# 在真实搜索源态上用片段三点放置全配体

已知量是普通配体 PDBQT、已制备化学边表、受体与盒。首次原生变异后，子进程将全配体内部自由度和重原子坐标送给父进程，等待本次提议。片段只决定配准所用的三个点。所有重原子、氢、键、类型与当前扭转角继续随原生配体共同运动。运行不需要供评价的参考配体。

设三点源坐标均值为 $a$，三点受体坐标均值为 $b$。原三点拟合求行列式为 +1 的旋转 $R$，使 $\sum_i \|R(x_i-a)-(y_i-b)\|^2$ 最小。对任一重原子 $z$，新坐标为 $R(z-a)+b$。发送的全分子质心为 $R(\bar z-a)+b$；原生根部发送 $R(t_{root}-\bar z)+t_{center}$，四元数在原根部姿态左侧复合。内部扭转未变。相比跨轴三角，局部三点的边长对片段外扭转严格不变；更短的边及更小的垂线高度也会放大受体腔点的离散定位误差。具体命中与损失取决于全分子位置及被保留的对应。

生成和注入沿用同一实现链。

1. [anchors.py](anchors.py) 从重原子图删除 PDBQT 的 BRANCH 轴，按[预登记规则](../PROTOCOL.md)选首个有非共线三点的片段。
2. [Pocket.place](../../../tools/pocket_triangle.py) 使用整分子、片段或打乱的片段变换。腔点及其有序三元组与上一轮相同。先按三边平方误差选 1,024 项，才做配准及整分子形状计算。候选变换、分量和对应编号均随库保存。在全部对应上调用同一几何变换仅用于事后计算误差。
3. [live.py](live.py) 通过 [FIFO 控制器](../../cavity_field_rank/engine/live.py) 将 1,024 份全配体刚体变换按原顺序发回。Vina 当前类型场为每个可评分重原子进行一次有界插值。正值按当前 hunt cap 处理，越盒加墙势。首个严格最小者替换源态的刚体分量。整分子的形状分和类型分均有存档，本轮由后者决定赢家。[原生筛选](../../cavity_field_rank/engine/native_rank.h)
4. 原生链继续运行第一次 hunt、接受、第二次局优、轨迹容器、合池、显式精修及发布。观察器从生产源码构建并沿用父实验的补丁；生产源树为 `6c6f65b0b7445cb364415f2079906e0d67e015bb`。

`bank.npz` 的 `matched_rotation` 是本臂实际使用的旋转，含 sham 的移位；`center` 与 `site_triangles` 始终按保留该质心的对应排列。sham 的方向来自索引 `(i - shift) % 1024`。锚点、片段、是否回退分别见 `anchor_serials`、`component_serials`、`anchor_fallback`；`triangle_id` 为门前合法有序三元组表下标。库下标从 0 开始。几何误差与 xyz 以 Å 计、三边残差为 Å²，四元数和形状粗分无量纲。hunt 的经验分采用 kcal/mol 刻度。

没有合格的刚性分量时，片段及其打乱臂改用整分子三点并记录回退；打乱臂仍照常循环移动方向。全分子也无非共线三点、权威图与源重原子序号不符，或可用不共线腔点对应总数小于 1,024 时，此实验退出并保留已有输出。当前姿势面板与测定面板均有足够对应。

## 化学图与类型

姿势实验另有[官方解析器导出和对照](../pose/check_official_graph.sh)。以官方 1.2.7 源码的原解析函数生成已知输入的内存模型；重原子及氢的坐标与 PDBQT serial 唯一匹配。遍历邻接读出 XS 类型、邻居序号、邻居顺序和可转标志，对照 Pro 的普通推断和显式图。七个配体一致，其中 FABP4 的 2NNQ 与 3FR2 分别有 41 和 27 条边。[逐配体结果](../pose/results/graphs.json)

FKBP 的 AD4 搜索各自在原生解析器中读相同的普通 PDBQT，锚点只另外使用该配体的 SMILES 及其映射得到的化学边表。用上游 1.2.7 版本 `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` 和当前 Pro 的对象文件分别链接 AD4 读取器。两者直接读取各自解析出的 `atom.ad`、部分电荷及有序的原生邻接，并对应回 PDBQT serial。八份配体的原生表逐项相同，各自的重原子边也都与片段准备使用的边表一致。这项检查含两次探针编译共 3.494 墙秒、3.44 CPU 秒。[两实现的逐边实测值](ad4_parsed.csv.gz) · [差集及输入哈希](ad4_graphs.json)

## 独立复核与历史失败重放

48 个库，以不同路径从 PDBQT 嵌套刚体块确定片段，直接按点到直线的投影和距离确定配体三点，重算入选的 1,024 组有序对应。四元数 4×4 对称特征值问题独立重建刚体旋转。全部 49,152 份候选中与盲运行的最大逐坐标差为 2.17×10⁻¹³ Å，最大旋转矩阵分量误差为 5.51×10⁻¹⁴。48/48 次实际赢家均为存盘真实类型分的首个最小值。姿势面板没有触发片段回退。旧十格整分子原法的 8×10 个共享数组与父实验逐项严格相等。本次独立检查用时 3.403 墙秒、3.403 CPU 秒。[逐库复核](bank_reconstruction.csv)

重放了历史 BRD9–91813 在逐项和向量化的根部变换间出现的分叉。以同一普通输入、库和真实类型场开始，仅在候选 359 的根部 y 分量使用原存档的相邻数值，两者相差 1.11×10⁻¹⁶ Å。两遍均选择 359，但两条轨迹合计执行的 hunt 次数由 1,973 变为 1,823，重新得到的两个发布池分别与历史两份 PDBQT 逐字节相同。回放开销分别为 0.198/0.196 与 0.161/0.160 墙秒/CPU 秒。[输入输出哈希与计数](roundoff.json) · [本次实际日志与发布坐标](roundoff_raw.tar.xz)

## 复现

仓库根目录执行，以下路径变量均设为绝对路径。W 可按 pose 的复现步骤从本轮原始包解出。编译时把 Boost 安装路径传给构建器。实测 GCC 14.2.0、Boost 1.83.0、3 个编译任务，共 28.44 墙秒、59.29 CPU 秒，运行程序 SHA-256 为 `8ffe5fbcc37384c48f6d397adbc50db87558305c6dd00a488b5ef6587f48dfb9`。二进制哈希可受绝对构建路径影响。

~~~sh
D="$PWD/benchmarks/fragment_anchor/engine"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
sh benchmarks/cavity_field_rank/engine/build.sh "$B" --boost-prefix "$BOOST_PREFIX" --jobs 3
python "$D/verify_banks.py" "$W" "$W/reconstructed.csv"
cmp "$D/bank_reconstruction.csv" "$W/reconstructed.csv"
python "$D/replay_roundoff.py" "$W/replay" "$W" "$B/bin/vina"
(cd "$D" && sha256sum --check ARCHIVE.sha256)
# UP 为从官方 1.2.7 构建的源码与 build/obj，FW 是测定面板的工作目录
python "$D/check_assay_graph.py" "$UP/src/lib" "$UP/build/obj" \
  "$B/src/lib" "$B/bin/obj" "$BOOST_PREFIX" "$FW" "$W/ad4-probes" "$W/ad4-graphs"
~~~

从普通输入启动单实例，`--dump` 写入新目录。省略此选项即可使用同一预测路径仅发布结果。父进程仅在子进程结束成功后替换最终输出。

~~~sh
python "$D/live.py" --engine "$B/bin/vina" --arm fragment --dump "$W/one-trace" -- \
  --ligand "$W/inputs/82I-ligand.pdbqt" --ligand_bonds "$W/inputs/82I-ligand.bonds" \
  --receptor "$W/inputs/6Y7L.pdbqt" --center_x 128.0266 \
  --center_y -2.8662399999999995 --center_z 35.45272000000001 \
  --size_x 24 --size_y 24 --size_z 24 --seed 91813 --cpu 1 --exhaustiveness 2 \
  --max_evals 25000 --num_modes 9 --energy_range 100 --pose_selection late \
  --search_counters --out "$W/one.pdbqt"
~~~
