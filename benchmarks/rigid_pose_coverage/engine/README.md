# 以完整配体的刚体位移选取有序腔点对应

[coverage_gate.py](coverage_gate.py) 接收第一次真实变异后的重原子坐标、序号、受体与盒。三点来自源态中最远的一对重原子及到该连线最远的第三个重原子。沿用已有腔点、非共线有序对应、正旋转 SVD 配准、原生场排名和单个 hunt 名额。每个候选的原始 ID 是 `Pocket.triangles` 的行号，库内 index 对应实际送入场的顺序。盲预测独立于参考配体和测定数据。[预登记](../PROTOCOL.md)

## 距离和可保证的量

用列向量表示第 $i$ 个源重原子的去质心位置 $x_i$。变换 $T$ 的输出为 $Rx_i+t$，其中 $t$ 是变换后的全分子质心。令 $C=N^{-1}\sum_i x_i x_i^\mathsf T$。完整配体的固定原子身份均方位移为

$$
d^2(T,U)=\|t-u\|^2+\mathrm{tr}\big[(R-S)C(R-S)^\mathsf T\big].
$$

交叉项含 $\sum_i x_i=0$。对 $C=V\Lambda V^\mathsf T$ 取 $L=V\sqrt\Lambda$，令 $z(T)=(t,\mathrm{vec}(RL))$。两行 $z$ 的平方欧氏距离即上式，单位 Å²。算法始终使用受体坐标系，保留全部方向与平移差异。变换间不再最优叠合。各重原子等权，绕某轴的小角旋转贡献与垂直于该轴的均方距离成正比。一段远离转轴的末端能主导某些变换的差异。

全盒内的 $M$ 个有序对应构成有限域。若 $g$ 为保留集合，$r(g)=\max_{T}\min_{U\in g}d(T,U)$ 量化最差覆盖的对应。精确算术下的最远点遍历给出度量 $k$-center 问题的二倍半径上界：在所选 $k$ 点上追加离它们最远的一点，共 $k+1$ 点。其任意两点相距至少 $r(g)$。最优的 $k$ 球中至少一个覆盖这 $k+1$ 点中的两点，球直径至少为 $r(g)$。实现按预登记比较舍入至 $10^{-10}$ Å² 的最近平方距离以处理浮点并列，实际半径直接由距离测量。几何相同的不同对应保留不同 ID。[Gonzalez，1985](https://doi.org/10.1016/0304-3975(85)90224-5)

设 $q(T)$ 为对单个晶体的化学对称校正 RMSD。对每个允许的原子置换，由三角不等式得到 $|q(T)-q(U)|\le d(T,U)$，因此 $\min_{U\in g}q(U)\le\min_{T\text{ 可行}}q(T)+r(g)$。此界允许本来刚好小于 2 Å 的候选在离散化后丢失，也允许有多个等价原子时把物理相近的姿势算得很远。盒内碰撞、离体积质心很远的稀有方向和狭窄但化学不利的区域都参与覆盖。体积测度、盆地占有率、化学势、温度、游离参考等没有进入 $d$，热力学权重无法由此推断。

## 实际运行

`edge` 先筛三边残差再配准，只支付 1,024 个配准。其余四臂先配准所有原始有序腔点对应，以每个重原子在给定盒三个闭区间内为可行条件。每次按 512 个变换分块计算所有重原子的位置。按字典序生成对应，以其原始 ID 打破预筛并列。完整规则、种子和退出条件见预登记。保留 $M$ 行 12 维特征及 $M$ 行最近距离，每步对 $M$ 个点更新一次，1,024 步共 $1024M$ 次 12 维平方距。选点内存为 $O(M)$，无两两距离矩阵。用纠正浮点舍入后的半正定 $C$，允许平面、直线和重合几何的度量检查。实际预测要求非共线三点及足够对应。

五臂保留的完整配体接受相同形状查询（用于重放）以及相同真实类型粗排。Vina 粗排在子进程的搜索场中进行，AD4 面板的粗排沿用原面板逐重原子类型图的插值与正值 cap；电荷加权静电图、去溶剂化图和显式极性氢未参加这个粗排，完整 AD4 在下游搜索中计算。实际 hunt 接收全配体根部平移和旋转，保留原扭角、连通图、类型和电荷。观察器和逐项四元数合成沿用父实验。`uniform_work` 的 `radius_A2` 为其计算后舍弃的覆盖库半径。真实均匀门的几何性能另从实际选择 ID 测量。库内时间仅分阶段计，比较总成本应使用从制备完的输入启动到发布结束的进程树计时。

## 独立计算与失败重放

- 父面板 16 个源态，现有 `edge` 对原库的每一项均重算相等，包括 12 个历史数组或标量。合成几何包含三维、平面、线性、单点，每组 780 个变换对，和直接逐原子的最大平方距离差为 $2.51\times10^{-12}$ Å² 以内。另检查同分时原 ID 裁决、零距离补到预算和可行数不足时退出。[数值检查](results/metric_check.json)
- 六个新源态共 178,404 个门前刚体，按另一套 Horn 四元数配准直接得到所有重原子坐标。与正式配准全域最大坐标差 $3.49\times10^{-13}$ Å 以内。在正式运行使用的 SVD 可行域上，使用全部 $3N$ 坐标直接计算距离，重建五臂 30,720 个入选 ID 和次序，并核实 30/30 个原生场首位；首位均无精确并列。六格的两种 uniform 输出逐字节相同。连同下述替代支持域计算，独立复核实耗 15.599 墙秒、15.558 CPU 秒。[逐库复算](results/coordinate_check.csv) · [Horn，1987](https://opg.optica.org/josaa/abstract.cfm?uri=josaa-4-4-629)
- 初次独立复核直接用四元数配准重判严格闭盒。3ANR 96107 有八个对应在两种数值求解器间翻转，距盒面最多 $2.14\times10^{-14}$ Å；SVD 接受 28,735 个，Horn 接受 28,727 个。放回正式 SVD 支持域时直接几何复核通过。在替代支持域独立重选时，`feasible_edge`、`coverage`、`uniform` 分别与原库重合 1,023、335、419 项。前者主要改变截断边界；最远点有序依赖前序选择；均匀无放回在变长的可行 ID 数组中取样。其余五格两方法的可行域及入选 ID 全同。这一阶段只测替代数值解下的候选身份。[初次失败](results/initial_coordinate_check.log) · [每个翻转对应的盒面距离](results/boundary_flips.csv)
- 在结果解封前提交[数值敏感性追加方案](../SENSITIVITY.md)后，对同一格另做三次真实在线预测。逐候选保留生产 SVD 变换，只替换可行掩码；每个共同 ID 的质心、旋转、三边残差及形状分逐位相等。三臂仍分别选择原对应 15146、3452、3405。覆盖臂的库内编号从 109 变为 115，均匀臂从 113 变为 111。第一次 hunt 后分别为 2.415、3.027、3.013 Å；三个发布文件与各自主实验逐字节相同，首位误差 0.477、0.478、0.465 Å。附加实耗 3.709 墙秒、3.721 CPU 秒（包括额外四元数计算与存档）。这一格表明候选名单变化可在单一场赢家处截断。该结论仅限此源态与既定场。[源坐标到发布的复算](results/boundary_replay.csv) · [源态、两掩码差、刚体库、原生评分、轨迹和发布文件](results/boundary_replay.tar.xz)
- 再次重放 BRD9 91813 历史浮点分叉。两支均将 359 号刚体送入同一 hunt，根部一个坐标差约 $1.11\times10^{-16}$ Å。最终总 hunt 调用 1,973 和 1,823，两份发布 PDBQT 分别与最初历史原始件字节一致。合计 0.300 墙秒、0.298 CPU 秒，观察器 SHA-256 为 `8ffe5fbcc37384c48f6d397adbc50db87558305c6dd00a488b5ef6587f48dfb9`。[哈希与计数](results/roundoff.json) · [此次原始运行](results/roundoff_raw.tar.xz)

## 事后检查：原始候选的短接触

将保留库的每个全配体与固定接收体逐重原子查询，统计任一组跨分子重原子距离小于 2 Å 的候选。分母是变换数；相近对应可以重复计入，不能作独立分子解释。下表每项依次给出全部保留变换、以及其中晶体 RMSD ≤2 Å 变换的短接触计数。未凭此结果改门。[逐格计数](results/steric_counts.csv) · [全部近姿的最短间距](results/steric_near.csv)

|门|历史 16 格，全部|历史，近姿|新 6 格，全部|新，近姿|
|---|---:|---:|---:|---:|
|edge|14,763/16,384|11/42|5,267/6,144|6/21|
|feasible_edge|14,843/16,384|14/45|5,243/6,144|7/23|
|coverage|16,274/16,384|3/3|6,080/6,144|1/1|
|uniform|16,178/16,384|5/7|6,016/6,144|6/10|

A1120 A⁻ 93241 的恢复项 7672 在注入时有 1.236 Å 的最近重原子距离和 1.945 Å 晶体 RMSD。取既有知情分叉的实际坐标，第一次 hunt 将最近距离推至 2.497 Å，仍保有 1.941 Å RMSD；二次局优达到 2.726 Å 和 2.170 Å。对这批起始几何执行 2 Å 硬拒绝会删除这一个可修复的近姿，并删除覆盖法在两组中的其余三个近姿。[知情分叉逐态距离](results/steric_capture.csv)

## 复现

仓库根目录执行。W 为本轮姿势面板按其说明准备和生成的输出目录，P 为父轮 fragment_anchor/pose 归档的展开目录，TMP 为新的工作目录。C++ 观察程序直接沿用父轮构建入口。Python 3.12，NumPy 2.3.5，SciPy 1.17.0。只需运行数值检查时无需重新编译。

```sh
D="$PWD/benchmarks/rigid_pose_coverage/engine"
mkdir -p "$TMP/verify"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
sh benchmarks/cavity_field_rank/engine/build.sh "$B" --boost-prefix "$BOOST_PREFIX" --jobs 2
python "$D/check_metric.py" "$P" "$TMP/verify/metric_check.json"
python "$D/verify_coordinates.py" "$W" "$TMP/verify"
python "$D/steric_audit.py" "$P" "$W" "$TMP/verify"
cmp "$D/results/coordinate_check.csv" "$TMP/verify/coordinate_check.csv"
cmp "$D/results/boundary_flips.csv" "$TMP/verify/boundary_flips.csv"
python "$D/boundary_replay.py" predict "$W" "$TMP/sensitivity" "$B/bin/vina"
python "$D/boundary_replay.py" evaluate "$W" "$TMP/sensitivity" "$TMP/verify/sensitivity.csv"
python benchmarks/fragment_anchor/engine/replay_roundoff.py "$TMP/roundoff" "$P" "$B/bin/vina"
(cd "$D" && sha256sum --check ARCHIVE.sha256)
```

单个实时运行使用 `python "$D/live.py" --engine "$B/bin/vina" --arm coverage --dump "$TMP/new" -- ...`，双连字符后的参数与配对 `edge` 运行完全相同。完整命令由[姿势面板](../pose/run_panel.py)的 `live` 与父面板的 `native` 生成。
