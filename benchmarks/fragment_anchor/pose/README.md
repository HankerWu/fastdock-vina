# 同源态的三种定位方式与实际发布池

每个受体–配体–种子实例在第一条轨迹第一次变异后取得完全相同的内部构象。whole 用整分子三点，fragment 用一个刚性片段的三点，sham 将片段库的取向作非零循环移位。三臂用同一个真实受体类型场对完整配体评分，1,024 份尝试分配一个 hunt 名额。11–36 个重原子参与受体坐标系下的对称性校正 RMSD（Å）。位置允许全配体平移、正旋转及原生可动扭角；受体、化学微态、氢和水处理沿用准备好的输入。受体盒来自接收晶体已知口袋。查询晶体仅在预测完成后的独立评价进程使用。

| 主评价量，16 格 | 整分子 | 片段 | 片段打乱 | Pro 25k | Pro 50k | 官方 25k | 官方 50k |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top‑1 均值 | 5.653 | 6.759 | 5.851 | 9.017 | 8.305 | 5.889 | 5.046 |
| Q1 / 中位 / Q3 | 2.052 / 2.725 / 8.647 | 3.016 / 5.204 / 8.635 | 2.380 / 3.867 / 6.753 | 4.851 / 7.947 / 14.751 | 2.367 / 7.701 / 13.595 | 2.374 / 4.765 / 8.246 | 2.356 / 2.893 / 6.419 |
| P90 / P95 | 14.653 / 15.916 | 13.402 / 14.153 | 14.002 / 15.084 | 15.460 / 15.618 | 15.460 / 15.574 | 13.432 / 14.739 | 12.830 / 13.790 |
| Best@5 均值 | 4.830 | 5.919 | 5.107 | 8.400 | 7.652 | 4.472 | 3.951 |
| 实际发布池 Oracle 均值 | 4.800 | 5.891 | 4.783 | 8.177 | 7.531 | 3.830 | 3.785 |
| Regret 均值 | 0.853 | 0.868 | 1.068 | 0.840 | 0.774 | 2.060 | 1.261 |
| 首位 >5 / >10 Å，空池数 | 6 / 3，0 | 9 / 3，0 | 6 / 3，0 | 12 / 7，0 | 10 / 7，0 | 8 / 3，0 | 5 / 3，0 |

六个不同分子先各自平均种子与微态，再对分子等权，整分子、片段、打乱的 Top‑1 为 6.051、6.830、6.663 Å。六个新种子实例的三臂 Top‑1 为 3.990、4.499、3.819 Å，Oracle 为 3.415、3.864、2.945 Å。按未参与该机制选择的 FABP4 四格单列，三臂 Top‑1 为 3.235、3.997、4.540 Å，Oracle 为 2.471、3.130、3.310 Å。新种子的 A1120 两微态中，打乱臂平均首位 2.377 Å，整分子与片段为 5.501 和 5.503 Å。这些均值只有一个新家族和一组配对微态。[分子和实例分布](results/distributions.csv) · [逐格](results/cells.csv) · [逐模式](results/modes.csv)

## 近姿如何失去候选或排名

三个事前指明的失败格，保留库内全配体最小误差均下降。A1120 A⁻、HA（93241）分别 2.629→2.313 Å、2.633→2.372 Å，4GHI→5TBM（93257）2.261→2.165 Å。没有一格新进入 2 Å。4GHI→5TBM 的另一种子 93241 则从库内 6 个近姿、粗排第 1 位为 1.179 Å，变成库最近 2.954 Å、所选第 1 位 12.630 Å。两个 BRD9 的库内最近值还略有改善，但首个近姿的真实类型排名由 2、12 滑到 3、42。FABP4 的 3FR2→2NNQ（94727）片段库只保留 1 个 ≤2 Å 候选，排第 310，整分子库有 3 个且首个排第 16。各臂 16 格中覆盖 2 Å 的实例数为整分子 5、片段 3、打乱 4，共 42、25、21 个近姿。[49,152 个受保留候选](results/candidates.csv.gz) · [各阶段与选中项](results/rank_stages.csv)

全部 48 个盲首选经 hunt、二次局优后都曾入池。每臂仅 1 个注入 ID 可追到最终发布，整分子排第 2、片段第 9、打乱第 8。整分子 4GHI–93241 的近姿经 hunt 可至 0.928 Å，第 2 步逐出；后来该格的近姿属于独立提议。每一步的容器成员、原 ID 是否在场、合池和精修输入/输出均可由原始双精度坐标核对。[逐步容器](results/residency.csv) · [各轨迹及合池坐标误差](results/stages.csv)

两种新输入和五种既有输入中，显式 Pro 图和官方源码解析出的有序邻接表、XS 类型与转动标志逐项一致。额外探针使用官方计数构建的编译对象，其上游 SHA 为 3c65c0b，相关解析及原子类型源与官方 v1.2.7 相同。32 次官方计数重放与官方发行可执行文件的发布池字节相同。956 个发布模式均无已定义手性中心反转，重键长漂移最大 0.001297 Å。所有首位的受体最近重原子距离 ≥2.065 Å。宽窗的低排名仍有穿插，最短 0.918 Å。[化学检查](results/chemistry.csv) · [化学图](results/graphs.json)

## 代价

同一逻辑核、每格新进程、按种子轮换顺序交替计时。16 格从普通制备输入到发布，整分子、片段、打乱、普通 Pro 墙钟分别为 26.864、26.531、25.839、10.048 秒，整棵进程树 CPU 分别为 13.468、13.308、12.981、4.978 秒。片段相当于普通 Pro 的 2.64 倍墙钟和 2.67 倍 CPU。各臂均单次运行完成取源、建腔、生成、排序、局优和发布。计时关掉完整阶段与候选快照，逐格核对了输出字节、场首名与原生计数。主机在本次测量期间有调度竞争，墙钟与 CPU 差异较大，时间比仅适用于本次实现和环境。首次无轨迹复算与控制组在同一逻辑核重叠，48 格共 124.953 墙秒、47.019 CPU 秒，已存于 setup 归档的 contended_cost.csv；本段采用控制组完成后的 64 次交替计时。

每个提议臂在 16 格各作 411,648 次全分子形状查询、411,648 次额外类型场插值、78,608 次建腔最近邻与 12,832 次邻域基数查询。实际局优求值（hunt + 二次局优 + 显式精修）整分子 804,524、片段 804,706、打乱 804,643、普通 Pro 804,204。实际含梯度场查询为 20,114,758、20,113,961、20,113,455、20,112,055；形状/类型粗排计数以重原子位置为单位。前述几何门的数据遍历和 Kabsch 运算计入总时间，未折算成场查询。准备好的受体与化学图作为共同输入；历史数据制备和晶体下载没有重新计时。

归档另列正式带完整轨迹的 48 次运行（39.882 墙秒、39.704 CPU 秒）、64 次官方/Pro 及双预算运行（83.507、58.414 秒）、32 次官方计数对照（64.976、44.212 秒）、32 次普通 Pro 发布程序重放（42.944、14.083 秒）及最终 64 次实耗重放（89.282、44.735 秒）。墙钟均为逐进程时长之和，多组作业曾并行。正式坐标评价另耗 5.162 墙秒。输入、二进制、场计数及几何查询见[逐格实耗](results/cost_cells.csv)、[各臂总耗](results/cost_total.csv)、[输入散列](results/inputs.json)及[可执行文件与 CPU 条件](results/binary_context.json)。

## 运行和复算

从仓库根目录执行。仓库已有接收体、配体和参考结构的归档；原始 PDB/CCD 与准备来源在 [human_panel/sources.csv](../../cubic_grid/human_panel/sources.csv) 及 [BRD9 输入](../../torsion_proposals/fixtures/) 中。两种 A1120 输入来自 [acid_states](../../cubic_grid/acid_states/)。W 为空工作目录，OUT 为新的结果目录。源码、编译器及四组可执行文件按上轮的[编译步骤](../../cavity_field_rank/pose/README.md)定位，OBSERVER 按本轮 [engine](../engine/) 构建；PRO 为普通发布二进制，PRO_COUNTER 为带位置计数的同源 Pro。CORE 替换成本机允许的单一逻辑核索引。所有正式运行均为 CPU，局优计数包括显式精修，场查询计数按位置累计。

~~~bash
D=$PWD/benchmarks/fragment_anchor/pose
mkdir -p "$W/tmp" "$OUT"
export TMPDIR="$W/tmp" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/panel.py" prepare "$W"
python "$D/panel.py" controls "$W" --core "$CORE" --pro "$PRO" --observer "$PRO_COUNTER" --official "$OFFICIAL" --counter "$OFFICIAL_COUNTER"
python "$D/panel.py" dock "$W" --core "$CORE" --observer "$OBSERVER"
python "$D/panel.py" timing "$W" --core "$CORE" --observer "$OBSERVER" --pro "$PRO"
python "$D/evaluate.py" "$W" "$OUT"
python "$D/cost.py" "$W" "$OUT"
python benchmarks/cavity_field_rank/pose/check_graphs.py "$W" "$PRO_PROBE" "$OUT"
bash "$D/check_official_graph.sh" "$UPSTREAM_COUNTS" "$BOOST_PREFIX" "$W" "$PRO_PROBE"
~~~

仅重算结果时，用另一空目录 W，先解档后校验。raw 下的三种轨迹文件引用同一归档内 dump 中的原始快照；全部归档均须解开。

~~~bash
(cd "$D" && sha256sum --check --status ARCHIVE.sha256)
for archive in "$D"/archives/*.tar.xz; do tar xf "$archive" -C "$W"; done
(cd "$W" && sha256sum --check --status "$D/results/raw_files.sha256")
python "$D/evaluate.py" "$W" "$OUT"
python "$D/cost.py" "$W" "$OUT"
for f in cells modes rank_stages stages residency chemistry distributions cost_cells cost_total source_identity; do
    cmp "$D/results/$f.csv" "$OUT/$f.csv"
done
cmp <(gzip -dc "$D/results/candidates.csv.gz") <(gzip -dc "$OUT/candidates.csv.gz")
~~~
