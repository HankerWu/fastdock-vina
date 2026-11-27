# 同一搜索改动能否改善 FKBP12 表观 Ki 排序

在 Holt 等 Table II 的八个化合物上，用按块无放回变异替代每步独立抽取，三倍分离对从 5/15 增至 6/15。受体间差别很大。在 1FKG 上正确数从 3 增至 9，在 1FKH 上都是 5。每个受体先对三枚固定种子的 Top‑1 分数取中位数，再对两个受体等权。等权汇总新增 8–4、8–11 两对，丢失 8–7 一对。每个点都来自同一人 FKBP12、同一立体化学系列、同一 10 °C、pH 7.8 的酶抑制测定。化合物 10 的 300–600 nM 区间只参与有确定间隔的方向统计。其余七个点用于秩相关。[事前标准](PROTOCOL.md) · [原始 Table II](https://doi.org/10.1021/ja00075a008) · [逐对差值](results/pairs.csv)

| 等权汇总 | 三倍分离对 | Spearman ρ | Kendall τ-b | 8 强于 10 的预测差 | 9 强于 11 的预测差 |
|---|---:|---:|---:|---:|---:|
| 官方 Vina 1.2.7 的 AD4 模式（归档） | 6/15 | −0.179 | −0.143 | +0.421 | +0.1925 |
| 当前 Pro，AD4 | 5/15 | 0.071 | 0.048 | −0.894 | +0.154 |
| 无放回，AD4 | 6/15 | 0.321 | 0.238 | −1.0135 | +1.0325 |

正预测差表示较强配体的分数更低。实验 8 对 10 约强 30–60 倍，9 对 11 约强 43 倍。改动后，两对立体化学悬崖仍然只判对一对。八个化合物均为形式电荷 0。分数位移（新减旧）与重原子数、TORSDOF 的描述性 Spearman 分别为 −0.206、−0.753。相同 TORSDOF=11 的 8–11 号仍可向相反方向变化。本组只有一个骨架，无法估计净电荷影响、总体校准或未知分子的命中概率。[逐化合物描述符](results/descriptors.csv) · [逐受体排序](results/assay.csv)

## 同一批分子的结构检验

化合物 8/SB3 和 9/SBX 分别有 1FKG 和 1FKH 的共晶坐标。按含手性全图身份匹配查询，在原始受体框架重算这两个分子对两个受体、三个种子的所有发布姿势；跨晶体用事先确定的口袋 Cα 配准。含两个不同分子、四种配体–接收结构条件、每条件三次重复。其余六个化合物只参与亲和力测定。几何数据是本轮同测定实验的补充分析。[实验坐标 1FKG](https://www.rcsb.org/structure/1FKG) · [实验坐标 1FKH](https://www.rcsb.org/structure/1FKH) · [逐模式误差](results/geometry_modes.csv)

| 12 格，Å | Top‑1 均值 | Q25 | 中位数 | Q75 | P90 | P95 | Best@5 均值 | Oracle 均值 | Regret 均值 | ≥8 Å |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 当前 Pro | 5.073 | 3.679 | 5.443 | 6.074 | 7.279 | 8.668 | 4.352 | 4.352 | 0.721 | 1/12 |
| 无放回 | 5.250 | 2.348 | 5.528 | 6.882 | 8.719 | 10.062 | 3.806 | 3.717 | 1.533 | 2/12 |

两个分子的同晶体对照说明均值为何分开。8 在 1FKG 上的首姿均值从 5.728 降至 3.473 Å，9 在 1FKH 上从 2.673 升至 5.688 Å。后者发布池的 Oracle 也从 2.623 升至 3.554 Å。两分子的跨、同晶体合并后，池可提供更好坐标，首姿排序损失扩大。[逐格配对](results/geometry_cells.csv) · [完整分布](results/geometry_distributions.csv)

## 预算、化学与复核

沿用同一份 AutoGrid 4.2.6 地图、配体、电荷及 AD4 经验分、24 Å 盒与 0.375 Å 网格。详见[原输入及场的定义](../../cubic_grid/fkbp_ad4/README.md)。两臂各 48 次单线程盲搜，均有 432 个输出。每次 3 条轨迹、每轨迹名义 24,000 次求值，AD4 无末端显式精修。实际求值分别为 3,458,851 和 3,458,880 次，墙钟累计 55.506 和 54.368 s，CPU 累计 53.961 和 52.780 s，含读图和写出。沿用已生成的两套地图。旧 Pro 全部 48 个池与上次归档逐字节相同。读出 864 个模式，以 SMILES 原子序号恢复每个指定手性中心，身份改变为 0，重键长相对制备输入的最大漂移为 0.001398 Å。以上几何检查覆盖共价长度及手性。官方 48 个池由独立读数程序直接从原归档重算。[实际求值](results/counts.csv) · [每次耗时](results/runs.csv) · [独立重算排序及池坐标比较](results/independent.json) · [全部 96 个新池与日志](pools.tar.xz)

实验结果 CSV 中 parent 为当前生产版，fixed 为本次研究用无放回副本。按[研究补丁的构建指令](../../receptor_linkage/aurora/internal_loops/mutation_scan/README.md)分别准备 PRO、SCAN。两者均为绝对路径，WORK 选择新目录。运行前取消 FD_SCAN_IID 等研究环境开关。

```sh
D=benchmarks/proposal_sweep/fkbp
env -u FD_SCAN_IID -u FD_SCAN_TRACE python benchmarks/accepted_bfgs/fkbp/fkbp.py . "$WORK" "$PRO" "$SCAN" --search-counters
python benchmarks/accepted_bfgs/fkbp/check_pools.py "$WORK" benchmarks/cubic_grid "$WORK/independent.json"
python "$D/audit.py" "$WORK" "$WORK/audit"
cmp "$WORK/independent.json" "$D/results/independent.json"
```

只复算已发布数据时，将 pools.tar.xz 解至 WORK 并拷入 results/runs.csv 及 results/scores.csv，随后运行后两条 Python 命令。化学及晶体检查直接读取仓库原始制备归档。全程没有从启发式取样频数或经验分推算标准自由能。
