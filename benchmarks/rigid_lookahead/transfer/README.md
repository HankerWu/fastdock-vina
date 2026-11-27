# 人 RBP4 新结构上的独立验证

在预测前冻结两种维生素 A 转运蛋白 RBP4 配体与一个独立的人性激素结合球蛋白 SHBG 结构。RBP4 取 1RBP 的视黄醇（RTL，21 个重原子，2.0 Å）与 3FMZ 的非视黄醇类 A1120（2T1，28 个重原子，2.9 Å）。两者分别原位和相互交叉，两颗种子，共八个配对单元。两份沉积配体重原子均满占有率、无替代构象，口袋没有重建重原子。交叉对齐用 61 个共同口袋 Cα，残差 1.191 Å。3FMZ 的实验研究报道与视黄醇结合时不同的口袋环构象。SHBG 的二氢睾酮晶体结构 1D2S 中，处理流程会移除一个距配体 7.013 Å 的 Ca，按预先冻结的 8 Å 标准排除。新机制在这里取得的是一个家族、两个不同化学结构的结果。[预测前冻结的协议](https://github.com/HankerWu/fastdock-vina/blob/8f4ebac7363dc671066b15eb4c3533b12172ef03/benchmarks/rigid_lookahead/transfer/PROTOCOL.md) · [1RBP](https://www.rcsb.org/structure/1RBP) · [3FMZ 的原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC2658061/) · [1D2S](https://www.rcsb.org/structure/1D2S)

| 方法 | 2 Å 首位 / 前五 / 候选池命中（各 8 次） | 候选池命中时的首位成功 |
|---|---|---|
| 上游 Vina，3,000 / 12,000 | 0/1/1；1/2/2 | 0/1；1/2 |
| 分支基线 Pro，3,000 / 12,000 | 0/1/1；2/3/3 | 0/1；2/3 |
| 粗分原样 / 随机一步 / 梯度一步 | 三组均为 0/0/0 | 三组分母均为零 |

预算为单轨迹名义上限，Vina 和 Pro 每次运行四条轨迹。每条轨迹的实际停止点与显式精修消耗均计入逐格输出。三组刚体试验各评价 227,136 个父姿势，精修 72 个候选。梯度和随机对照均对相同父姿势计算梯度，再各执行一次试探评分。其精修分别消耗 3,598 和 3,991 次梯度求值，原样组为 4,067 次。独立进程分别重建两套受体图、全程执行单个刚体方案，原样、随机、梯度三组单核墙钟分别为 2.81、4.27、4.53 秒，峰值 RSS 为 143.0、143.3、143.4 MiB，各组全部 72 行最终排名与第一次运行一致。Pro 两档八次单核命令的墙钟和为 1.55、2.49 秒，实际梯度求值总数为 101,171、388,828；上游为 9.16、8.98 秒，实际 101,682、389,246 次。两类计时的复用范围不同。紧邻阈值的 1RBP 原位种子 43117 中，梯度解 RMSD 为 2.010 Å，仍计失败。3FMZ 配体进入 1RBP、种子 43103 的首位 RMSD 由原样 12.859 Å 降至 2.872 Å，随机对照得到 2.836 Å。[结构明细](results/cells.csv) · [分段成本](results/cost.csv) · [整进程成本](results/single_arm.csv) · [原生和 Pro 计时计数](results/baseline.csv)

## 固定配体内部坐标能达到多好的刚体配准

观察盲测结果后，加做了四个使用实验标签的诊断。保留原 CCD 化学图生成的独立输入构象，对晶体等价映射逐一求最佳刚体配准，再向同一个显式优化器交付一个最优起点。A1120 的初始形状与晶体之间，即使任意连续平移与旋转也至少相差 2.225 Å。改动可动扭转的局部优化在两个受体中分别得到 2.539 Å（1RBP）和 2.741 Å（3FMZ）。这为此单构象粗筛的表达能力提供定量上界。它仍允许其他初始构象、其他精修盆地或蛋白构象取得更低 RMSD。

视黄醇的刚体拟合下界为 1.654 Å。在原位 1RBP 中，同一输入投到最佳刚体位置后，显式精修达到 0.989 Å，分子间与内部能量之和为 −10.684 kcal/mol，优于梯度两种子的盲选首位 −4.271 与 −9.438。这里既存在正确局部盆地，也存在支持该盆地的显式排序。在 3FMZ 接收体中，视黄醇最佳刚体起点反而被同一精修器推至 2.431 Å。这四次计算使用实验坐标选择起点，只用于定位机制。视黄醇双键立体构型沿用 CCD，A1120 按 CCD 输入中性羧酸；本轮未测口袋微态。两个种子与互换的同系列受体共享信息，单个分子与受体之间的变动不能当作八项独立的生物证据。[四个诊断](results/fit.csv) · [逐姿势能量与 RMSD](results/ranked.csv)

## 复算

在仓库根目录使用 [输入登记与环境版本](sources.json)。复现现有预测可直接用归档；由 PDB、CCD 重新生成输入可执行 prepare_transfer.py。基线源为 2e610c4，上游为 3c65c0b。上游计数器使用 benchmarks/pair_gradient/upstream_counts.patch。每次上游和 Pro 的计数回放均要求输出 PDBQT 与未计数程序逐字节相同。probe 的构建参见 rigid_translation_scan/build_probe.py，链接同一份基线源码的目标文件。[二进制指纹](results/binaries.json)

~~~bash
D="$PWD/benchmarks/rigid_lookahead/transfer"
W=$(mktemp -d)
mkdir -p "$W/input"
tar -xJf "$D/inputs.tar.xz" -C "$W/input"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/replay.py" baseline "$W/input" "$W/base" --pro "$PRO" --upstream "$UPSTREAM" --counter "$COUNTER"
python "$D/replay.py" rigid "$W/input" "$W/rigid" --pro "$PRO" --probe "$PROBE"
python "$D/audit.py" "$W/input" "$W/base" "$W/rigid"
python "$D/best_rigid_fit.py" "$W/input" "$W/fit" "$PROBE"
~~~

outputs.tar.xz 保留 32 个基线文件及对应完整日志、计数、资源使用，216 个刚体初始与精修后的姿势，四个标签辅助诊断的输入、输出及日志。也可在不运行对接的情况下解包并复核 717 条坐标 RMSD。独立实现使用 SciPy 的向量配准，最大数值差为 7.11×10⁻¹⁵ Å。[审计结果](results/audit.json)

~~~bash
tar -xJf "$D/outputs.tar.xz" -C "$W"
cp "$D"/results/baseline{,_poses}.csv "$W/transfer-baselines/"
cp "$D"/results/{cells,poses,ranked}.csv "$W/transfer-rigid/"
python "$D/audit.py" "$W/input" "$W/transfer-baselines" "$W/transfer-rigid"
~~~
