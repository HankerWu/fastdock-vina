# 键输入与搜索器的配对对照

本项是在看过主实验后的解释性实验。将完整图 Pro 的原命令只去掉 `--ligand_bonds 文件`，分别用原平移和枢轴搜索，保留同一配体 PDBQT、受体、口袋盒、种子、每次 2×25,000 软限以及 late 发布。生产入口由 `Vina::set_ligand_from_file_with_bonds` 切换到常规 PDBQT 解析。两个图均把配体交给同一搜索器，从相同制备坐标、相同随机种子和预算设置开始；化学图可能改变随后的数值轨迹。全部 12 个新池均已发布。[追加范围](PROTOCOL.md) · [所有 CLI 调用](results/commands.json) · [主实验](../rigid_coupled/README.md)

生产解析器读出的 AT9283 显式图有 36 条键，推断图 35 条。推断漏掉制备序号 30 与 31 之间的环丙基 C–C 边。其余键、可转标志和原子类型一致。danusertib 的 41 条边、82I 的 12 条边逐项一致，包括每个原子的邻接顺序、XS 类型与可转标志。每个配体的显式边也从源 SDF 或 CCD 对应的化学结构核对。加入或省略键文件后，后两种配体的两种子、两搜索器共八个发布池字节完全一致。AT9283 四池全部变化。[实际双向邻接表](results/adjacency.csv) · [图差](results/graphs.json)

| AT9283 种子 | 搜索器 | Top‑1 完整→推断，Å | Oracle 完整→推断，Å | 官方推断图 Top‑1 / Oracle，Å |
|---|---|---:|---:|---:|
| 81371 | 当前 Pro | 8.998 → 8.983 | 7.693 → 7.366 | 8.981 / 7.168 |
| 81371 | 枢轴 | 9.217 → 8.926 | 7.377 → 8.060 | 8.981 / 7.168 |
| 81373 | 当前 Pro | 8.931 → 9.006 | 8.768 → 8.189 | 8.987 / 7.328 |
| 81373 | 枢轴 | 8.887 → 9.584 | 8.548 → 9.087 | 8.987 / 7.328 |

两枚新种子都没有重现历史种子 61231 的近晶体返回。AT9283 改用推断图后，原平移的发布池 Oracle 都略降，枢轴的 Oracle 都升高。两种算法对缺失同一条化学键的响应方向不同。在其他两种配体，官方与 Pro 的差异仍存在，例如 82I–81373 的官方首姿 3.110 Å，Pro 两种图都是 8.121 Å。这里官方采用原生输出选择，Pro 保留 late 选择，另外原生与 Pro 的搜索和数值实现也有差异。本对照可排除 82I 该格的化学图边集差异，无法把其剩余差异分配给某一个实现环节。[全部成对实例](results/paired.csv)

在推断接口上，六个结构–种子格的原平移 / 枢轴首姿均值分别为 9.441 / 9.656 Å，中位数 8.623 / 8.737，Q25–Q75 为 8.157–9.000 / 8.234–9.420，P90 为 11.587 / 12.111，P95 为 12.877 / 13.374。Best@5、Oracle、Regret 的均值分别为 8.275 / 8.709、8.011 / 8.594、1.430 / 1.062 Å。两臂均有 6/6 的首姿 ≥8 Å、1/6 >10 Å。累计 grad/energy 联合求值分别为 301,540 / 301,377，单线程墙钟 1.793 / 1.720 s，CPU 1.777 / 1.706 s。十二个池共 107 个模式，两条坐标计算路径的最大 RMSD 差低于 3×10⁻¹⁴ Å；用正确图检查所有重键，最大长度漂移 0.001240 Å，未发生已定义手性中心翻转，最短受体–配体重原子距离 2.405 Å。此处的几何和单次成本是已知三种配体的诊断结果。[分布](results/distributions.csv) · [逐池](results/cells.csv) · [逐姿势及化学检查](results/qa.csv)

## 复现

沿用[主实验的原始准备和构建说明](../rigid_coupled/README.md)。AURORA、BRD、BUILD、BASE_SRC、BOOST_PREFIX、WORK 为绝对路径，BUILD 含 base 与 pivot，BASE_OBJ 是构建 base 的 obj 目录。探针只调用生产的两个配体解析函数，通过制备态坐标唯一回连 PDBQT 序号，输出每个原子的有序邻接、类型和可转标志；始终独立于搜索命令运行。计时仅含真实对接的启动、读入、搜索与输出。

```sh
D=benchmarks/receptor_linkage/aurora/internal_loops/rigid_graph_control
mkdir -p "$WORK"
python "$D/graph_control.py" run "$AURORA" "$BRD" "$WORK" --build "$BUILD"
mapfile -t OBJECTS < <(find "$BASE_OBJ" -name '*.o' ! -name 'main.o')
c++ -O3 -DNDEBUG -std=c++14 -I"$BASE_SRC/src/lib" -I"$BOOST_PREFIX/include" \
  "$D/read_graph.cpp" "${OBJECTS[@]}" -L"$BOOST_PREFIX/lib" -Wl,-rpath,"$BOOST_PREFIX/lib" \
  -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options \
  -o "$WORK/read_graph"
python "$D/compare.py" "$AURORA" "$BRD" "$WORK" "$WORK/read_graph"
```

复核已归档坐标及化学图证据时解开本目录 raw.tar.xz 到新的 WORK，复制 results/jobs.csv、results/counts.csv 至 WORK，运行 `graph_control.py verify`，然后同样运行 compare.py。输入及两可执行文件、逐池哈希分别见 inputs.json、jobs.csv。主实验的完整图结果可由 compare.py 从其原始归档逐池核对。
