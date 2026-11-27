# 修正原生局优线搜索后的 Aurora A 配对回放

冻结输入来自提交 `066b314`。保留两个历史受体、三个配体、四枚种子、两种 Pro 名义预算和原生 Late 候选选择，重新运行修正二进制的 48 次对接。官方 Vina 1.2.7 的 24 个结果复用已验证的字节级归档，进入结构对照和输出清单，其计时也沿用原数据。两种有实测坐标的配体每臂 16 格，首姿分别与 2J50 A 和 2W1G A 对齐；另列 2J50 B 与 6CPG A/D。Imatinib 只进入分数和计时。全部结构和状态假设取自[原协议](../PROTOCOL.md)。

| 方法 | 首姿均值，修复前 → 后 | 后中位数 | 后 Q25 / Q75 | 后 P90 / P95 | Best@5 均值，前 → 后 | Oracle 均值，前 → 后 | 平均 regret，前 → 后 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 官方 Vina，共用对照 | 9.199 → 9.199 | 9.661 | 8.038 / 10.818 | 12.572 / 12.725 | 7.038 → 7.038 | 6.748 → 6.748 | 2.451 → 2.451 |
| Pro，25k | 9.106 → 9.031 | 9.072 | 7.489 / 9.411 | 10.629 / 11.642 | 7.922 → 8.077 | 7.442 → 7.464 | 1.664 → 1.566 |
| Pro，50k | 8.896 → 9.117 | 9.081 | 8.907 / 10.134 | 10.629 / 10.898 | 8.118 → 8.401 | 6.590 → 7.141 | 2.305 → 1.976 |

单位 Å。25k 档的配对首姿误差 9 格下降、7 格上升，差值中位数 −0.598 Å；50k 档为 6 格下降、10 格上升，中位数 +0.012 Å。首姿超过 10 Å 的次数分别从 6/16 变为 4/16、从 4/16 变为 5/16。原生线搜索在蒙特卡洛和显式精修中都会运行；全部 48 个重新生成的池均发生变化，归档的逐轨迹步数或求值次数也均发生变化。只能把该实验的全流程差异归于修复后的算法，不能单独归于末端精修。独立 Gemmi 和 RDKit 从坐标与 PDB 图重新算出 1,076 个模式–参照误差，最大差 3.82×10⁻¹⁴ Å。[逐格新旧配对](results/solver_pairs.csv) · [逐模式](results/geometry_modes.csv) · [全部分布](results/geometry_distributions.csv) · [独立核验](results/geometry_independent.json)

AT9283 减 danusertib 的种子中位首姿分差，在两种受体之间仍取相反符号。25k 档 2J4Z 为 −0.270 → −1.700，2X6E 为 +1.005 → +0.611。50k 档分别为 −1.385 → −1.660、+0.796 → +0.514。两者同缓冲方法的实测亲和力相差不到三倍；排除 imatinib 不同 DMSO 条件后，可靠分离对的分母为零。没有可从这批数据推出的排序精度或标准态自由能增益。[分数](results/scores.csv) · [同方法的差值](results/differences.csv)

24 次 25k 搜索及显式精修的实际求值总数 1,214,088 → 1,206,536；50k 为 2,413,806 → 2,406,033。25k 全部进程的墙钟/CPU 总计从 10.118/9.975 秒到 8.419/8.394 秒，50k 从 14.132/14.015 秒到 13.145/13.116 秒。首次与修复后计时来自同一主机的不同时段，每个条件只计一次，不宜据此声明加速。两组均使用已准备受体和配体，时间不含输入制备及结构评估。[逐次时间、各阶段计数、原始输出哈希](results/solver_runs.csv) · [汇总](results/costs.csv)

## 复现

从仓库根目录运行，FIXED 指向包含原生线搜索修复的可执行文件，实际生成候选的二进制 SHA-256 为 `2b2f6afe30cd84cbff46187149939945f16c8ee75c42549a9d9001db4beda2a0`。重新运行只调用 Pro。24 个官方池连同来源哈希逐个校验。若只复核已有坐标，先将原实验 docking.tar.xz 解到 NEW，再将本目录 outputs.tar.xz 解入同一 NEW，无需执行 replay.py。

```sh
A=benchmarks/receptor_linkage/aurora
D="$A/solver_replay"
K=benchmarks/torsion_kinase_affinity/evidence.tar.xz
mkdir -p "$OLD" "$NEW" "$REF" "$OUT"
tar xf "$A/docking.tar.xz" -C "$OLD"
tar xf "$A/references.tar.xz" -C "$REF"
python "$D/replay.py" "$OLD" "$NEW" "$FIXED"
python "$A/geometry.py" "$NEW" "$REF" "$OUT" --kinase-archive "$K"
python "$A/frame_check.py" "$NEW" "$REF" "$OUT" --kinase-archive "$K"
python "$A/analyze.py" "$NEW" "$OUT"
python "$D/compare.py" "$A/results" "$OUT" "$OLD" "$NEW" "$OUT"
cmp "$OUT/solver_pairs.csv" "$D/results/solver_pairs.csv"
```
