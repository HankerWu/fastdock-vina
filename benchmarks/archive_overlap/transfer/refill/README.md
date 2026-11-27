# 归档回填的独立结构试验

在旧面板解释了第二次局部优化的预算分配后，单独登记了三个新结构。检查 PDB 的 SOURCE 记录后发现，登记时误将 1T7T 雄激素受体写为人源，实际为黑猩猩。报告保留预登记编号，分别列出人和黑猩猩数据。人 JNK3 1PMV 未通过 Meeko 制备，多处残基出现意外的残基间连接。保留输入及诊断日志；无替补。

| 来源、复合物 | 条件数 | 四种算法的 Top-1 / Top-5 / 输出池（2 Å） | 回填研究版跳过次数 |
|---|---:|---:|---:|
| 人 BCL-XL，2YXJ / ABT-737 | 2 | 0 / 0 / 0 | 0 |
| 黑猩猩雄激素受体，1T7T / 双氢睾酮 | 2 | 2 / 2 / 2 | 0 |
| 人 JNK3，1PMV / 537 | 未完成 | — | — |

每个结构使用两个种子，每种算法每次运行四条轨迹。16 条“归档一旦填满便不再切回填充状态”的研究轨迹均未跳过二次局优。三个 Pro 程序在每个合格条件得到逐字节一致的输出和相同的六项工作量；比较覆盖 4 个条件、2 种配体。三个 Pro 版本分别执行 113548 次 hunt、80398 次二次局优、3693 次显式精修的联合求值，总计均为 197639。四次 CLI 墙钟时间累计为父版 1.644 秒、修正版 1.602 秒、研究版 1.674 秒，上游为 7.962 秒。仅有一次交错计时，此处无法估计稳定的墙钟差。该面板检验了无触发时的行为一致性，触发后的收益需从其他数据判断。

人 BCL-XL 的 Pro 首姿为 8.384 和 12.690 Å，池中最近候选为 7.922 和 7.514 Å。粗粒度的输出端分析只能定位至“正确姿势未进入输出池”。黑猩猩雄激素受体的 Pro 首姿分别为 0.641 和 0.602 Å。用 RCSB 结合态 SDF、独立原子对应和坐标变换复算全部 112 个姿势，最大 RMSD 差异为 1.78×10⁻¹⁵ Å。各姿势的碰撞、RMSD、能量和计数见 modes.csv、cells.csv 及完整日志。

## 复现

在本目录执行下列命令。父目录 README 说明 Python 依赖、上游计数器和 CPU 绑定。

```bash
sha256sum -c SHA256SUMS
mkdir -p work
tar -xf inputs.tar.xz -C work
tar -xf outputs.tar.xz -C work
cp cells.csv modes.csv sources.json work/
python panel.py audit work
```

原始 RCSB 和 CCD 输入可用 `python panel.py prepare work` 重建，精确复算应使用归档冻结输入。下列命令同时运行四臂，并验证带分阶段及回填计数的回放与三个正式 Pro 程序及上游程序分别给出完全相同的输出。

```bash
OPENBLAS_NUM_THREADS=1 python panel.py dock work --cpu 4 --upstream "$UPSTREAM" \
  --counter "$UPSTREAM_COUNT" --parent "$PARENT" --candidate "$FIXED" \
  --hold "$HOLD" --hold-counter "$HOLD_COUNTER"
```
