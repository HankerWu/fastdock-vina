# 用连续访问的盆地选择一个候选

两条 FKBP1A 轨迹表明正确接受态已经存在，按最低 hunt 能量记录的却是另一个盆地。这里测试一个无需再打分、无需拟合的替代解释：连续几次局优都留在同一处，可能比孤立的一次低分更能识别值得精修的盆地。只考虑原生搜索已经接受且未送第二次优化的姿势，与已存姿势均相距至少原合并半径。新颖接受态落在同一段首个姿势的半径内就延续该段；选择最长段、以该段最低能量的姿势为代表。已提交的接受态或返回旧盆地时重置，拒绝不累加。定义、并列处理、固定 K 及留出种子见[事前协议](PROTOCOL.md)。Monte Carlo 的提议分布、接受温度、局部优化与单轨迹相关性都影响连续段长，因此它只是一种待检验的搜索信号。[盆地跳跃原始论文](https://doi.org/10.1021/jp970984n)

| 测试 | 选择规则 | Top-1 / Top-5 / 返回池命中 | 替换名额的运行 | 实际联合评估 | 端到端秒数之和 |
|---|---|---|---:|---:|---:|
| 已有 6 靶点，66 条轨迹 | 父版本 | 14 / 16 / 16 | 0 | 814,825 | 11.995 |
| 同上 | 最低能量 | 14 / 17 / 17 | 38 | 813,806 | 11.876 |
| 同上 | 连续访问 | 14 / 17 / 17 | 40 | 813,863 | 12.064 |
| 新种子，8 靶点、30 个方向 | 父版本 | 6 / 7 / 7 | 0 | 368,833 | 6.038 |
| 同上 | 最低能量 | 6 / 7 / 7 | 13 | 368,833 | 5.618 |
| 同上 | 连续访问 | 6 / 7 / 7 | 19 | 369,058 | 5.788 |

所有阈值变化都来自开发组的同一 FABP4 轨迹。新种子中三臂在同一批运行条件命中，含一个 NUDT1 命中；CHIT1 四个方向均未命中。增加了替换次数，未改善首姿或带来第二个新覆盖。该机制作为可构建的研究补丁保存。每个预测条件均为单线程、单轨迹，MC 提议、hunt 及第二次优化的计数均经日志核对一致，进入显式优化的名额也相同。一份候选的端到端时间只测了一次，短任务时间受调用顺序、机器负载等影响。[原始输出与逐靶点结果](results/summary.csv) · [实际运行及散列](results/cells.csv)

在仓库根目录，令 S 为此前六家族准备好的结构目录（含已解包的 cubic_grid/human_panel/prepared.tar.xz），N 为新家族的结构目录（transfer 的输入归档），W 为新的工作目录。在依赖环境中执行：

```bash
D=benchmarks/accepted_retention
python "$D/core/build.py" parent "$W/parent" --jobs 3
python "$D/core/build.py" candidate "$W/energy" --jobs 3
python "$D/core/build.py" recurrence "$W/recurrent" --jobs 3
OPENBLAS_NUM_THREADS=1 python "$D/persistence/paired.py" "$S" "$N" "$W/replay" \
  "$W/parent/build/vina" "$W/energy/build/vina" "$W/recurrent/build/vina" --cpu 3
python "$D/persistence/audit.py" "$S" "$N" "$D/persistence/results" "$W/summary.csv"
cmp "$W/summary.csv" "$D/persistence/results/summary.csv"
```
