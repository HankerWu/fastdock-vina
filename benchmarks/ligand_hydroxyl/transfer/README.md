# 11 个新人源复合物的羟基方向验证

对接前冻结了人 ERRγ 的七个苯酚、双酚与萘酚配体，以及人 HSP90α 的四个含酚抑制剂。每个家族使用一个共同受体，另测各自的共晶受体。在 20 个不同查询–受体组合中，对每个组合运行两个种子、两档预算和两个实际 CLI，共 160 次。官方 AutoDock Vina 1.2.7 输出 706 个候选，分支起点 b991eb7 输出 701 个。随后仅改变冻结 Pro 候选池的排序函数。[预登记](PROTOCOL.md) · [结构、组分、分辨率、占据率、来源与 SHA-256](results/structures.csv) · [实际命令与成本](results/cost.csv)

每格依次列出首位 / 前五 / 池内命中次数，阈值为对称校正的重原子 RMSD ≤2 Å。家族之间分别汇总，同一配体的重复种子不代表两个独立样本。羟基模型为主协议登记的 MMFF94s 配体 C–OH 扭转积分加到原有报告分数。候选顺序变化之前，三维坐标和复合物的可搜索自由度均已固定。

| 受体与任务 | 名义每轨迹上限 | 池数 | 官方 Vina | 起点 Pro | Pro 加羟基积分 |
|---|---:|---:|---:|---:|---:|
| ERRγ 交叉对接 | 3,000 | 12 | 11/12/12 | 12/12/12 | 10/12/12 |
| ERRγ 交叉对接 | 12,000 | 12 | 12/12/12 | 12/12/12 | 9/12/12 |
| ERRγ 自对接 | 3,000 | 14 | 10/10/10 | 10/10/10 | 10/10/10 |
| ERRγ 自对接 | 12,000 | 14 | 12/13/13 | 12/13/13 | 12/13/13 |
| HSP90α 交叉对接 | 3,000 | 6 | 0/0/0 | 0/0/0 | 0/0/0 |
| HSP90α 交叉对接 | 12,000 | 6 | 1/1/1 | 1/1/1 | 1/1/1 |
| HSP90α 自对接 | 3,000 | 8 | 3/3/3 | 3/3/3 | 3/3/3 |
| HSP90α 自对接 | 12,000 | 8 | 4/4/4 | 4/4/4 | 4/4/4 |

正确候选在 ERRγ 交叉对接的全部池中都已存在。五个首位退化集中在 6I65 与 6I67 两个单酚。6I65 低预算种子 43103 的近晶体候选为 0.647 Å，原分 −6.549；另一候选为 4.755 Å、−6.437。加入羟基积分后，前者上升至 −6.272，后者仍为 −6.437。前者的酚氧距 Glu275 OE2 3.185 Å，距 Leu309 主链羰基氧 3.042 Å。只从新修正项除去单个受体接触，近晶体分数分别减少 0.105、0.186。两次干预分别进行，位点积分为非线性。五例的共同机制是近晶体姿势承担 0.261–0.354 的正修正，胜出的错误姿势没有可修正的受氢接触。原始晶体研究观察到 Glu275/Arg316 的极性锚定；本实验不能单独估计几何、背景介质与水的能量份额。HSP90α 的可排序正确候选本来就很少。[逐姿势](results/poses.csv) · [逐池](results/cells.csv) · [各对照及条件成功率](results/summary.csv) · [逐接触的冻结删除实验](results/loss_contacts.csv)

ERRγ 在同一 6I63 受体上与论文 Table 1 的直接 ITC Kd 对照。主分析按预案排除绝对立体化学不清楚的 4-sec-butylphenol 测定，共六个配体、七对至少三倍差异。低预算 Pro 与积分均 3/7，Spearman −0.029、Kendall 0.067；高预算均 4/7，Spearman 0.086、Kendall 0.200。纳入 6I66 后为 4/9 与 5/9。分子对共享配体。这一窄系列给出的排序证据不足以支持默认启用新项。HSP90α 仅用于姿势评价。[直接 ITC 数据与条件](https://pmc.ncbi.nlm.nih.gov/articles/PMC11105698/) · [逐对预测](results/affinity_pairs.csv) · [相关性](results/affinity.csv)

计算在 CPU 上进行，每个对接子进程使用一个线程。低/高预算 Pro 中位墙钟 0.170/0.216 s，搜索求值中位数 12,162/48,097，显式精修 395/394。本机官方 Vina 墙钟 0.847/0.941 s。资源记录包含进程启动，机器负载可影响时钟。十一种受体参数化耗时累计约 27.3 s，配体准备 0.58 s；额外的方向模型受体预处理 0.211 s，701 个 Pro 姿势计算各个对照共 0.554 s。上述事后计时包含全部比较臂和共享中间量。此处的附加计算没有带来准确率上限提升。[二进制哈希](results/binaries.json) · [软件版本](results/packages.json) · [准备时间](results/preparation.csv) · [额外打分时间](results/rescore_cost.csv)

3TUH 链 A 的 Thr176 在 PDB 缺少主链羰基 O，常规 PDBFixer 补出的 O 与 CA 碰撞，导致 Meeko 化合价检查中止。该输入使用相邻已解析 CA–C–N 肽平面，在 1.23 Å 键长的两个方向外角平分线上补 O。该 O 离配体最近重原子 19.90 Å。所有输入均重新准备一次，22 个配体/受体 PDBQT、11 个制备后蛋白 PDB 与参考坐标和映射的 metadata.json 逐字节一致。本次没有保留结晶水或有机添加剂。2XJJ 的甘油距配体 3.89 Å，HSP90α 的部分配体有可质子化胺；这些都是模型适用边界。所有配体交替构象和受体拟合均另用 PDB 距离成键、模板键级和 Meeko SMILES 序号计算 RMSD，对全部 1,407 个姿势的最大差 1.2×10⁻¹⁴ Å。

## 复现

在仓库根目录运行。版本见结果文件。下载并重做输入可用 prepare，重跑实际 CLI 需要将两个路径替换为表中版本的二进制，源数据和准备好的输入已一并收录。

```bash
B=benchmarks/ligand_hydroxyl/transfer
W=$(mktemp -d)
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tar -xJf "$B/inputs.tar.xz" -C "$W"
tar -xJf "$B/pools.tar.xz" -C "$W"
python "$B/run_transfer.py" geometry "$W"
python "$B/check_transfer.py" "$W"
python "$B/evaluate_transfer.py" "$W" "$W/report"
python "$B/diagnose.py" "$W" "$W/report"
SUM=$(realpath "$B/results/scientific.sha256")
(cd "$W/report" && sha256sum -c "$SUM")
```

从原始坐标开始可使用以下命令重建输入及候选池。费用、计时不参与科学结果的哈希。

```bash
python "$B/run_transfer.py" prepare "$W"
python "$B/run_transfer.py" dock "$W" --native "$NATIVE" --pro "$PRO"
```
