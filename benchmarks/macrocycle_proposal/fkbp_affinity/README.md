# 人 FKBP12 的实测 Kd 与大环几何

FK506 与其近似物 APX879 在同一项人 FKBP12 等温滴定量热研究中的 Kd 分别为 2.7 ± 0.5 nM、119.4 ± 15.8 nM。两者在 25 °C、pH 7.0、50 mM 磷酸钠和 50 mM 氯化钠中测量，各有三次重复。由于亲和力与溶解度不同，研究者对两种药物使用了相反的滴定方向和不同注射浓度。44.22 倍来自同一化学系列的一对分子，换算的相对结合自由能差为 +2.245 kcal/mol。两个化合物均以电中性状态制备，重原子数为 57 与 61，本组没有严格等大小配对。两点秩相关仅编码这一个方向，记录逐条件方向和分差即可。[原文 Table 2 与 ITC 实验方法](https://doi.org/10.1128/mbio.03000-21) · [测量数值](measurements.csv)

两个查询分别来自 [1FKJ/FK5](https://www.rcsb.org/structure/1FKJ) 与 [6VCU/R27](https://www.rcsb.org/structure/6VCU) 的 CCD 化学图，以同一随机种子的一次 ETKDGv3 嵌入及 MMFF 松弛准备，初始配体未使用晶体坐标。6VCU 的 A 链配体与另一个配体的重原子最近 2.357 Å，依原规则未采用此受体。在任何对接运行前登记使用 [1FKF](https://www.rcsb.org/structure/1FKF) 作为预备受体。两个入选的人源受体都来自 FK506 晶体，对 FK506 是重对接，对 APX879 是交叉对接。APX879 晶体参考经 34 个共同口袋 Cα 配准到两个受体，Cα RMSD 为 0.307–0.348 Å，变换后最近蛋白重原子距离 2.539–2.577 Å。[初始协议](PROTOCOL.md) · [对照参考几何](results/reference_geometry.csv)

每个分子、受体、种子和预算上比较一条单线程轨迹，24 Å 盒子、至多九个输出、能量窗 100。3000 和 12000 为每轨迹名义上限。每个体系使用种子 740317 和 191009。父版、首步 CCD、首步恢复输入环、原生 Vina 1.2.7 在首次评估中运行。观察首步后又断环之后，同样的 16 个条件被用于复核持续固定环的实验臂。[追加实验时点](../LOCK_PROTOCOL.md)

| 实现 | 输出姿势数 | 首姿 / 全池中断环数 | 首姿 RMSD 中位，Å | 池内最佳 RMSD 中位，Å | 联合能量–梯度调用 |
|---|---:|---:|---:|---:|---:|
| 上游 1.2.7 | 53 | 8 / 21 | 9.55 | 7.96 | 149860 |
| 当前 Pro | 67 | 2 / 5 | 9.88 | 8.43 | 131929 |
| 首步 CCD | 53 | 2 / 8 | 9.44 | 8.96 | 131211 |
| 首步恢复输入环 | 55 | 0 / 0 | 9.92 | 8.54 | 130442 |
| 持续固定输入环 | 77 | 0 / 0 | 8.62 | 7.46 | 128242 |

严重断环按真实闭环碳距 >3 Å 计。每臂有 16 个池；五臂的 Top-1、Top-5、输出池 Oracle 均为 0/16，所以条件排序成功率无分母。首步 CCD 在 FK506 上使断环从 1/33 增至 8/25，在 APX879 上从 4/34 减至 0/28。持续固定环时 77 个输出姿势的 21 个环原子经独立刚体叠合后，与输入环最大差 0.000552 Å。输入环本身与三个实验参考的形状相差 1.493–1.583 Å。分子内这一诊断可单独叠合，同实验结合姿势的主要 RMSD 始终在受体坐标系计算。[逐模式](results/modes.csv) · [预算与几何](results/geometry_cost.csv) · [环形状](results/ring_shape.csv)

| 相同受体，两种子的各配体分数取中位 | 1FKJ，3000 | 1FKJ，12000 | 1FKF，3000 | 1FKF，12000 |
|---|---:|---:|---:|---:|
| 上游 1.2.7 | +54.364 | +0.507 | −3.428 | −41.762 |
| 当前 Pro | +0.376 | +0.494 | −0.331 | −1.259 |
| 首步 CCD | +0.167 | +0.731 | +0.157 | −1.239 |
| 首步恢复输入环 | −0.181 | +0.440 | −0.756 | −0.268 |
| 持续固定输入环 | −0.103 | −0.537 | +0.288 | +0.474 |

表内是 APX879 减 FK506 的模型分差，正号与实测 Kd 方向相符。CCD 仅在 1FKF 的低预算改变配对方向。持续固定环在 1FKJ 的两档从正确变为错误，在 1FKF 的两档从错误变为正确。受体、预算、种子反复使用同一化学对。该结果未验证亲和力改善。[所有种子及汇总](results/kd.csv)

固定父版的 67 个原坐标，用共同零内能参考独立复算十项打印分解；另对 16 个原生首姿使用自身内能参考。全部 83 条记录中父版、CCD、持续固定环逐项相同。按事先定义的规则，只把严重断环挪到既有池的后部，父版在 1FKF–APX879 同一个种子的两档预算能把断环首姿换成严重阈值内的姿势，首姿 RMSD 分别从 10.710 和 9.395 Å 变为 8.865 Å。CCD 在 1FKJ–FK506 同一种子两档首姿分别从 9.160 和 9.347 Å 变为 6.409 Å，且闭环碳距回落到 1.554 Å。上游换位五个池，另有三个池全部严重断环。两种恢复输入环的方案均无须换位。[固定坐标](results/fixed.csv) · [逐池固定顺序干预](results/fixed_closure_order.csv)

除首次计时以外，每个父版–CCD 和父版–持续固定环配对各重跑四轮，固定逻辑 CPU、交替执行次序。所有 256 份额外预测的整个输出文件及六个阶段调用计数逐个复现。各条件四次墙钟的中位数再求和，父版–CCD 为 3.490 / 3.331 秒，父版–持续固定环为 3.494 / 2.921 秒。进入实验前还要准备输入。利用已下载的 PDB/CCD 缓存，独立进程重做 1FKJ、6VCU、1FKF 准备分别耗时 3.016、1.482、2.947 秒；6VCU 的成本包括受体排除检查。重新生成的两份独立配体、两份实际受体与输入哈希相符。在 CCD 上核对了每个化合物 14 个指定立体中心及 2 或 3 个指定双键构型。[实际运行和二进制散列](results/runs.csv) · [CCD 逐次耗时](results/repeated_timing.csv) · [固定环逐次耗时](results/lock_timing.csv) · [输入复核](results/input_verification.csv)

## 复现

取仓库根为工作目录。使用 Python 3.12、RDKit 2025.09.3、Meeko 0.7.1、Gemmi 0.7.4、NumPy 2.3.5、SciPy 1.17.0、OpenMM 8.6.1、PDBFixer 1.12.0。脚本的第一个参数为仓库路径。以下命令使用与上一层相同的基线和核心实验二进制。计数版上游另行核对输出逐字节相同。可将 CPU 改为本机允许的逻辑核。

```bash
D=benchmarks/macrocycle_proposal/fkbp_affinity
W=$(mktemp -d)
tar xJf "$D/evidence.tar.xz" -C "$W"
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
python "$D/repeat_preparation.py" "$PWD" "$W/work" "$W/observed"
python "$D/run_assay.py" "$PWD" "$W/work" "$W/observed" "$UPSTREAM" "$UPSTREAM_COUNTS" "$PARENT" "$CCD" "$RESET" --cpu 4
python "$D/run_assay.py" "$PWD" "$W/work" "$W/observed" "$UPSTREAM" "$UPSTREAM_COUNTS" "$PARENT" "$CCD" "$RESET" --lock "$LOCK" --cpu 4
python "$D/analyze.py" "$PWD" "$W/work" "$W/observed"
python "$D/ring_shape.py" "$PWD" "$W/work" "$W/observed"
python "$D/fixed_assay.py" "$PWD" "$W/work" "$W/observed" "$PARENT" "$CCD" "$LOCK"
python "$D/repeat_cost.py" "$PWD" "$W/work" "$W/observed" "$PARENT" "$CCD" --cpu 4
python "$D/repeat_cost.py" "$PWD" "$W/work" "$W/observed" "$PARENT" "$CCD" --lock "$LOCK" --cpu 4
```

若仅复算几何和配对方向，可解开原始归档，复制 `results/runs.csv` 到 `$W/observed/`，执行 `analyze.py` 与 `ring_shape.py` 两步。归档内有实际输入、沉积 PDB/CCD、全部 80 个原始对接池、日志和 16 份上游计数复跑日志。
