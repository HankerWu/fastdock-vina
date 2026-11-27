# 同一 MC 轨迹内在线放置

此 POSIX 单核研究入口从头启动一次普通 FastDock，在指定轨迹、指定变异之后取得真实重原子坐标与根部刚体位姿。父进程仅从此次受体、已有盒、当前分子形状和种子构造冻结的匹配库，对命令选定的一臂评分，通过私有 socketpair 回传一个刚体变换。子进程保留当时的扭转和随机数流，直接进入 hunt 并完成发布。源码改动仅在构建副本内；仓库的默认引擎及对接设置维持原状。

命令接受显式给出的受体、盒、种子、输出与 --cpu 1。--trajectory、--step 都从 1 计数，默认为 1。一次只有一条单配体刚性受体查询。在线入口先按 PDBQT 重原子 serial 排序，再选三点并破平，因而输入文件可使用不同的行顺序。对 24 个主姿势源、16 个 FKBP12 源及 6 个酸性微态源，从原始坐标头或 serial 数组核对，原顺序与升序全部一致。[原子顺序](atom_order.csv)进程异常、Vina 提前中止或未到选定步号均以失败退出；目标输出只在方案实际执行且引擎成功后原子替换。对 2 轨迹的 82I 要求第 20 轨迹已实测得到非零退出，已有目标文件仍完整。

~~~sh
# 仓库根目录；B 指向新的绝对构建目录，W 为主实验输入与已解开的原始姿势工作目录。
D=benchmarks/cavity_triplet/online
sh "$D/build.sh" "$B" --boost-prefix "$BOOST_PREFIX" --jobs 4
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/panel.py" "$W" "$B/bin/vina" "$W/online"
# 自有体系，以 BOX_X/Y/Z 等替换为实际坐标；其余普通 Vina 参数可以继续传入。
python "$D/dock.py" --engine "$B/bin/vina" --arm matched -- \
  --receptor "$RECEPTOR" --ligand "$LIGAND" --ligand_bonds "$BONDS" \
  --center_x "$BOX_X" --center_y "$BOX_Y" --center_z "$BOX_Z" \
  --size_x 26 --size_y 26 --size_z 26 --seed 91813 --cpu 1 --exhaustiveness 2 \
  --max_evals 25000 --num_modes 9 --energy_range 100 --pose_selection late --search_counters --out "$OUT"
~~~

若使用系统 Boost 可省略 --boost-prefix。在线版最后一次增量重编曾因 /tmp 空间耗尽失败，改用构建目录的 TMPDIR 后成功，最终二进制 SHA-256 见汇总。对照用同一新编译二进制，直接调用原生 CLI，没有开启在线放置。六格按系统和种子循环变更三臂运行顺序。所有测量均覆盖整个进程。两种在线臂包括启动 Python、导入科学库、启动引擎、实时取得坐标、建腔点、计算单臂粗分、继续 MC 至发布及退出；CPU 对 Python 和整个子进程树累加，每格都重新建图及构造腔点。[逐实例](jobs.csv) · [汇总及二进制哈希](summary.json) · [18 份原始池与日志](raw.tar.xz)

| 六格 | 墙钟秒 | CPU 秒 | 原生 hunt / second / explicit 求值 | 额外逐重原子形状查询 |
|---|---:|---:|---|---:|
| 普通 Pro | 1.7969 | 1.7648 | 291965 / 8323 / 1450 | 0 |
| 匹配 | 3.8825 | 3.8781 | 289452 / 10794 / 1324 | 151552 |
| 打乱 | 4.0012 | 3.9785 | 290852 / 9384 / 1371 | 151552 |

匹配、打乱各另用 29,478 次建场最近重原子查询及 6,158 次 8 Å 球计数。对照池与正式独立种子实验中对应的 18 套发布文件逐字节一致，hunt、二次局优、显式精修调用与求值、late 粗排计数全部相同。实测表明匹配流程约需普通 Pro 的 2.16 倍墙钟，打乱约 2.23 倍。这些单次测量描述当前实现成本，包含每格重新导入与准备。
