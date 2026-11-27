# Vina 刚体候选的原生场排序和单次在线搜索

在 MC 变异和 hunt 之间，1024 份刚体根部平移、四元数依次作用于同一个原态。每次保留全部扭角，在状态副本上用当前搜索的 `cache::eval` 计算与刚性受体的相互作用。每个可评分的重原子查询同一受体的一个 XS 通道，三线性插值后对正值 $u$ 使用 $10u/(10+u)$，再加各轴越界距离之和乘 $10^6$。Vina 场分的标度为 kcal/mol，坐标为 Å。当前面板中氢无 XS 通道。源状态的分子内项在 1024 个刚体变换中相同。通过稳定遍历选择最小的双精度分数，分数完全相等时取先出现的候选。

置换组在副本上按 PDBQT 重原子 serial 升序循环移动类型序列半圈。每次重排都记录修改前后的类型编号和实际进入三线性插值器的次数。第一个 hunt 使用所选的真实坐标、原有图和原生类型，后续 MC 正常继续。该研究编译沿用前轮观察器的轨迹、插入、逐出和合池跟踪，要求单进程单工作线程。未约束的游离配体、溶剂、受体构象与质子交换未在这个排序中采样。

## 构建和复现

仓库根目录下执行以下命令。B、S、P、W、T 均为绝对路径。W 是用本轮 pose/run.py 准备的工作目录，带有 inputs、banks 和 raw；P 是普通发布 Pro。B 和 S 分别为当前研究程序和原形状在线程序的构建目录。T 是独立核验目录。

```bash
D=$PWD/benchmarks/cavity_field_rank/engine
OLD=$PWD/benchmarks/cavity_triplet/online
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sh "$D/build.sh" "$B" --boost-prefix "$BOOST_PREFIX" --jobs 2
sh "$OLD/build.sh" "$S" --boost-prefix "$BOOST_PREFIX" --jobs 2
python "$D/replay_live.py" "$W" "$B/bin/vina" "$S/bin/vina" "$P" "$T/live"
python "$D/check_interpolation.py" "$W" "$P" "$T/maps" "$T/interpolation.json"
```

项目原始数据档中的 ligands、bonds、receptors 和盒即可启动单次在线搜索。`live.py` 使用独立的私有 FIFO；子进程保存指定步的原生状态，父进程读取重原子坐标并构造盲库，子进程筛选后继续原轨迹。一个候选进入 hunt。程序完成发布后才替换输出。以下是有实际数据的命令（W 指上一段的工作目录）。

```bash
python "$D/live.py" --engine "$B/bin/vina" --arm true -- \
    --ligand "$W/inputs/82I-ligand.pdbqt" --ligand_bonds "$W/inputs/82I-ligand.bonds" \
    --receptor "$W/inputs/6Y7L.pdbqt" \
    --center_x 128.0266 --center_y -2.8662399999999995 --center_z 35.45272000000001 \
    --size_x 24 --size_y 24 --size_z 24 --seed 91813 --cpu 1 --exhaustiveness 2 \
    --max_evals 25000 --num_modes 9 --energy_range 100 --pose_selection late \
    --search_counters --out "$W/82I-live.pdbqt"
```

`--arm perm` 使用预定类型置换；`--trajectory` 和 `--step` 从 1 计数，默认均为 1。对逐阶段观察和全库复算，pose/run.py 用 FD_FIELD_BANK 指向原始顺序的 1024×7 文本，FD_FIELD_ARM 选 true/perm，FD_FIELD_TRACE 导出以 0 起算的 index 和 score。FD_RIGID_INTERVENTION 规定在哪个事件执行。历史研究构建接受同样的源图/观察器配置。直播预测仅读取受体、盒和普通配体，使用暴露于生成器的重原子 serial 规定原子顺序。

## 重放、数值核验和实耗

10 个体系×种子实例在同一 Linux 逻辑核上交替运行普通 Pro、形状、真场、置换场四种流程，父控制器及其所有子程序均受同一 CPU affinity 约束。每格新启动一次，3 种提议均包含 Python 导入、取得真实源态、准备腔点和 1024 个形状查询候选、变换回传及发布。真场和置换还按冻结方案额外查询 1024×重原子数个原生场位置。成本不包含离线评价。形状使用前轮 socket 控制器，场排名使用 FIFO；两者准备和通信开销不同。同一环境的单次时间差不能解释为 C++ 场评分加速。

| 10 格 | 墙钟秒 | 整棵进程树 CPU 秒 |
|---|---:|---:|
| 普通 Pro | 2.768 | 2.737 |
| 形状 | 7.081 | 7.001 |
| 真实类型 | 6.992 | 6.955 |
| 置换类型 | 6.788 | 6.774 |

40 份发布池逐字节复现冻结实验，逐轨迹 hunt、二次局优、显式精修及接受计数相同。20 份类型筛选的赢家、分数、有效类型变更和实际查询也一致。[逐实例成本和二进制哈希](live_cost.csv) · [输出与日志](live_raw.tar.xz) · [原始文件散列](live_files.sha256)

用独立 Pro 的 `--write_maps` 导出原生场节点。配体文件用于列举通道，场值由相同受体和缓存参数产生；复核用正式注入日志里冻结图确定的 XS 序列。用存档的源坐标与三维矩阵直接恢复受体框架下的每个重原子，按实际原始盒中心计算坐标（地图头的 CENTER 文本精度较低）。SciPy 对各通道独立三线性插值，另算逐原子的正值 cap 和盒壁。对 20×1024 个分数，20/20 首名一致，最大绝对差为 0.003621 kcal/mol；原生文本导出将节点写成 4 位有效数。各库与运行时精确分数的 Spearman 相关最低为 0.999999899。诊断包括导图、读取和复算，共 3.693 墙秒、3.691 CPU 秒。[逐库核验与地图散列](interpolation.json)

首次在线适配尝试把逐项刚体运算改成批量运算。BRD9 91813 真类型仍选候选 359、场分相同，但该候选根部 y 平移差 $1.11\times10^{-16}$ Å，随后第一轨迹执行 836 步，对照 986 步，发布文件不同；第二轨迹计数保持一致。整个 1024×7 表中有 239 个分量变化，最大绝对差 $1.42\times10^{-14}$。恢复原逐项计算后，10 个刚体表的文本及上述 40 个池全部逐字节通过。该失败的日志、发布池和所选七个分量收录在归档 live_failure/。导图时曾组合 --ligand_bonds 和 --write_maps 触发 CLI 限制；又一次 --score_only 在成功写图后因探针原始坐标越盒退出。最终用普通 --write_maps 加 --randomize_only 从同一受体导图，随机探针坐标没有进入计算。两份失败日志均已收录。

在另一个干净构建目录全量重编，用两个编译工作进程共耗 35.637 墙秒、58.433 CPU 秒。该程序单次在线重放 82I 91813 的真场，发布池 SHA-256 与冻结预测相同。在两轨迹的同一实例请求第 3 条轨迹，控制器非零退出并保留原输出文件。[独立构建与退出状态](build_check.json)
