# RBP4 最后一个近姿在哪一步被删除

事前登记见 [PLAN.md](PLAN.md)。这是对已揭盲自由轨迹的因果诊断，研究臂使用上一轮未晋级的 `0bcae6e`。以正常 CLI 从头跑了 RBP4 两方向和 ESR1 两方向/两种子的 12 个格（各两臂的合计）。引擎在真实全局池每调用一次 `add_to_output_container` 后输出精确成员身份，保存全部到达时的二进制双精度坐标（用 17 位有效数字序列化）、构象及能量。12 个发布池的 PDBQT 与原预测逐字节相同，hunt、second、explicit 总调用和精修次数逐格复现。输入和输出哈希、实际命令见 [costs.json](raw/costs.json)。

**1RBP→3FMZ、种子 85237，错误发生于第 21 个末态到达时。** 27 个输入按三条轨迹顺序，各按轨迹内能量排序；表中的 ID 从 0 起算。候选 ID20 是第三条轨迹第 120 步保存、末态第 3 名的姿势。ID18 来自**同一条轨迹**第 383 步，排名该轨迹第一。ID20 到达时的九个居民依次为 18、0、1、5、9、12、7、15、16。

| 实际判断次序 | 被比较的值 | 结果 |
|---|---|---|
| 满池，比较最差居民 ID16 | ID20 的 −11.074503662 低于 ID16 的 +8.010726234 | 通过能量门槛 |
| 扫描优先居民 ID18 | ID18 的 −12.596351634 低于 ID20 的 −11.074503662 | 进入实际几何去重判断 |
| 原子索引配对的重原子间距 | ID20–ID18 距离为 1.892166366 Å，小于 2 Å | 拒绝 ID20，池中成员未变 |

同一对姿势可以同时通过轨迹内默认 1 Å 去重；全局合并使用 2 Å。ID18 保留到显式精修，成为发布首位。化学对称处理后的受体框架 RMSD，按旧记录的毫埃坐标复算为 1.862275 Å（ID20）、2.896610 Å（ID18）；按此次精确坐标复算为 1.862363 和 2.896519 Å。晶体标签只进入验后分析。将同一 27 项全局按能量稳定排序，再用同一容量和半径去重，九项身份和次序完全相同。由此否决了预先登记的“仅由容量能量比较引起”猜测，也未达到全局排序机制的注册晋级要求。

**12 格对照**共记录 278 次到达。Python 按原实现顺序重放 1,921 次条件判断，每次更新后的有序成员均与 C++ 一致，全部最终输入也逐个匹配旧精修记录。172 个消失的候选中有 127 个被更优或同分的重叠居民阻挡，3 个被更优新到达的重叠候选逐出，21 个在满池时能量不足，21 个在插入后的容量裁剪中逐出。正对照为各格在精修输入中晶体误差最低者的成员谱系。全量淘汰及全部远离晶体的候选均保留。按能量全局排序后，11/12 格有序身份相同；ESR1 的 1QKT→1X7R、85237 稳定臂多留一项（ID12），其余 8 项相对次序相同。所有 12 格的精修前 Oracle 均未改变。[逐格原生与全局排序结果](results/cells.csv) · [逐次有序池](results/arrivals.csv) · [逐次判断](results/comparisons.csv) · [每个被删候选的竞争者和去向](results/losses.csv)。

独立程序用能量优先级集合、欧氏范数和候选集合运算再次重建 278 次成员更新。6,046 个重原子坐标在此次精确引擎顺序与旧 PDBQT 记录顺序之间均唯一配对，所有晶体 RMSD 与旧坐标标签最大差为 0.000153 Å，跨过 2 Å 阈值的计数为零。[复核数值](verification.json)。表中的非距离比较未填距离；淘汰表另提供两者的诊断距离，包括容量短路时引擎无需计算的距离。

## 强制准入后，精修和终排会怎样处理它

[追加事前方案](REFINE_PLAN.md)在得知 ID18/20 的身份后登记。仅将原九个入池构象中的 ID18 换成 ID20，其余八个输入相同。正常流程重跑对照/替换两个进程，再用相反顺序独立重跑两个进程。四次均得到相同的 27 个源构象、能量、坐标和原生合并轨迹。两臂均保留原生显式精修、共享分子内参考的选取规则及输出舍入。为比较两个候选的分数，另对精修后的每个构象以**对照臂的同一分子内参考**做一次附加的只读评分，每次启动九项，这一观测不参与原生排序。

| 单槽干预 | 进入槽位 | 精修后发布 RMSD Å | 在各自九项中的发布名次 | 原生发布分 | 以对照参考重算的分数 | 该槽精修联合调用 |
|---|---|---:|---:|---:|---:|---:|
| 原九项 | ID18 | 2.896215 | 1 | −9.607 | −9.607070 | 43 |
| 只替换一项 | ID20 | 1.833424 | 1 | −8.541 | −8.421637 | 50 |

ID20 在原显式力场精修后仍小于 2 Å，并超过替换臂剩余八项（共同参考下第二名为 ID21 的 −8.372924）。在同一参考下，ID18 又比 ID20 低 1.185433 个 Vina 分数单位。即使两者都有机会经过显式精修，已实测的能量排序仍偏好较远的 ID18。两臂各自的原生分子内参考分别是 −0.678636183 和 −0.524434316；比较各臂的打印分会混入参考变化。此 Vina 分数用于比较特定构象，并未计算结合自由能。各臂分别重复时的完整原生合并轨迹、全精度精修轨迹和发布 PDBQT 逐字节相同。[九项的逐臂逐次能量分量与排名](refinement/scores.csv) · [四次运行收据](refinement/receipts.json)。

**成本。** 12 个归因作业含读入、建图、搜索、精修、发布和观测 I/O，合计 5.191136 CPU 秒，逐作业墙钟和 5.712780 秒，峰值 80,464 KiB。hunt/second/explicit 联合求值合计 331,325 / 101,319 / 2,373。四个干预及对照进程合计 2.204332 CPU 秒、2.278179 墙钟秒，峰值 60,972 KiB；每次 hunt/second 分别 31,273/4,750，explicit 从每对照的 198 增为每替换的 205。八个原有槽的逐项精修计数和共同参考分相同。输入复用原始制备文件；此处无新制备、打分参数搜索或生产重排。两种历史源码观察器的编译合计 184.977 CPU 秒、90.429 墙钟秒，另一个单槽观察器为 114.489 CPU 秒和 55.403 墙钟秒。[编译与安装失败记录](build_environment.json) · [原始编译日志](build.log) · [单槽编译](refine_build.log) · [独立读取与复算开销](analysis_costs.json)。

## 复现

仓库根目录开始。W 为空的临时工作目录。工具链和 Python 依赖沿用 stage_local_hunt；输入归档和参考完全复用。前五条验证命令不需要编译。每次真实运行的进程时间与构建绝对路径相应改变。

```sh
export OPENBLAS_NUM_THREADS=1
D=benchmarks/merge_replay
python "$D/replay.py" analyze "$D/raw" "$W/tables"
diff -qr "$D/results" "$W/tables"
python "$D/verify.py" "$W/verified.json"
cmp "$D/verification.json" "$W/verified.json"
python "$D/refinement.py" "$D/refinement"

sh "$D/build.sh" "$W/build-merge" "$BOOST"
python "$D/replay.py" run "$W/build-merge" "$W/merge-live"
python "$D/replay.py" analyze "$W/merge-live" "$W/fresh-tables"
diff -qr "$D/results" "$W/fresh-tables"

mkdir "$W/slot-source"
git archive 0bcae6ef09857f97656291bb4821d20b3b619bf5 src tools/build.py | tar -xC "$W/slot-source"
python "$D/observe.py" "$W/slot-source"
python "$D/refine_instrument.py" "$W/slot-source"
python "$W/slot-source/tools/build.py" --boost-prefix "$BOOST" --jobs 3 --output-dir "$W/slot-bin"
python "$D/refinement.py" "$W/slot-live" --binary "$W/slot-bin/vina"
python "$D/refinement.py" "$W/slot-live"
diff "$D/refinement/scores.csv" "$W/slot-live/scores.csv"
```
