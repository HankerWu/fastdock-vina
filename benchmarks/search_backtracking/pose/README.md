# 固定轨迹预算下的真实姿势回收

在三个靶点、十个配体–受体组合、两个种子、两档预算上，二次回退与现有 FastDock 各有 **3/40** 个条件的首姿达到 2 Å。相同的总数包含一次获得和一次丢失。不能据此更新默认搜索算法。

## 输入和对照

运行前协议见 [PROTOCOL.md](PROTOCOL.md)。六种从 CCD 化学图生成的查询配体共享给全部实验臂，不传递晶体姿势或受体接触模板。受体原配体提供口袋中心。1SQN 的制备过程补建了距离配体 7.69 Å 的重原子，因而根据预定 8 Å 规则排除该受体，保留其查询配体。其余五个受体均通过。源坐标、化学图、位点、配体骨架、晶体环境、来源 URL、序列化输入和摘要均收录在 [inputs.tar.xz](inputs.tar.xz) 和 [sources.json](results/sources.json)。

| 靶点 | 配体与 PDB | 查询重原子数 | PDBQT 扭转数 | 可评估组合 |
|---|---|---|---|---:|
| PPARα | GW409544（1K7L），GW735（2P54） | 38，33 | 13，9 | 4 |
| BRPF1 | 9F9（4UYE），OF-1（5FG4） | 29，26 | 4，4 | 4 |
| PGR | 孕酮（1A28），炔诺酮（1SQN） | 23，22 | 1，1 | 2 |

PPARα 配体均使用 CCD 描述的中性羧酸形式，实际结合时的电离状态未经测定。两个跨对接的参考位置在接收构象中分别有 1 和 2 对重原子接触低于 0.75 倍范德华半径和。BRPF1 的 5FG4 和 4UYE 配体距晶胞对称邻体最近分别为 2.98 和 3.36 Å，搜索时省略了该邻体和结晶水。这些条件限制了晶体姿势与孤立刚性受体之间的可比性。PGR 的两个配体具有相同的 Murcko 骨架。进一步泛化需要更多骨架与独立靶点。

## 结果

| 方法 | Top-1 | Top-5 | 返回池 oracle | oracle 命中时排第一 | 实际联合能量–梯度求值 | 局部优化调用 | 总墙钟秒 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 原生 Vina | 3/40 | 3/40 | 3/40 | 3/3 | 311085 | 6657 | 41.49 |
| 当前 FastDock | 3/40 | 3/40 | 3/40 | 3/3 | 310341 | 6648 | 6.87 |
| 失败时保留原点，二分回退 | 2/40 | 3/40 | 3/40 | 2/3 | 312355 | 14327 | 6.72 |
| 失败时保留原点，二次回退 | 3/40 | 3/40 | 3/40 | 3/3 | 311504 | 10232 | 6.46 |

CPU 绑定后，按病例改变臂的执行顺序，逐进程计时包含地图、搜索、精修与写出。表内总墙钟为一次交错实验；每次只有约 0.1–1 秒，不能把 0.41 秒总差视为稳定加速。六次预定复跑输出完全相同。Vina 参考是与 FastDock 使用相同构建流程的上游提交 3c65c0b，40 次均与官方 v1.2.7 Linux 二进制输出逐字节相同。其计数来自单独回放；40 次也均逐字节一致。FastDock 三臂共享同一研究二进制，基线的 40 个坐标文件及原有计数与未插桩的 bcb5d83 构建一致。报表不把诊断回放时间加入实际部署成本。共同配体和受体制备每套 2.01–3.61 秒，包括局部几何与来源处理，原始结构已下载。所有实验臂复用该开销。

| 条件 | 单元数 | 原生 Top-1 | 当前 Top-1 | 二分 Top-1 | 二次 Top-1 |
|---|---:|---:|---:|---:|---:|
| PPARα | 16 | 1 | 1 | 2 | 1 |
| BRPF1 | 16 | 0 | 0 | 0 | 0 |
| PGR | 8 | 2 | 2 | 0 | 2 |
| cap 3000（全部靶点） | 20 | 1 | 1 | 1 | 1 |
| cap 12000（全部靶点） | 20 | 2 | 2 | 1 | 2 |

二次回退把 hunt 和 save 的十步均失败次数从 7155 降至 2505，两阶段线搜索次数从 117929 升至 142101。实际联合求值增加 0.37%，局部优化调用增加 54%，反映相同预算内轨迹分配改变。PPARα 在种子 297149、预算 12000 时，1K7L 配体在自己受体上 RMSD 从 5.59 降至 0.58 Å；2P54 配体在同一受体上从 1.22 升至 2.56 Å。二分法在 PPARα 获益而 PGR 的两个首姿命中全部丢失。各臂首姿有受体冲突的单元分别为 4、5、4、5 个。706 个输出模式经独立化学图、受体对齐及坐标映射复算，最大 RMSD 差为 5.1×10⁻¹⁴ Å，未出现手性标记翻转，最大键长变化 0.00135 Å。这里的几何检查未评价全部构象应变。

## 正确姿势在哪一步丢失

120 次额外回放分别观察三个 FastDock 臂的 Monte Carlo 结束保存池、筛选后显式精修输入和精修输出。所有回放均复现主实验原始 PDBQT。每臂 40 个条件中只有 3 个在 MC 结束保存池含有 ≤2 Å 构象；这些命中均保留到最终输出。只有在二分组的 1K7L cognate（297149，12000）出现了近天然模式进入结果后被排在第二的情况（0.94 Å）。目前的阶段记录不能排除更早的 MC 轨迹曾短暂经过近天然构象；它们说明单改最终排序至多可在这组二分结果中增加一个首姿命中。逐候选记录见 [stages.csv](results/stages.csv)。

直接从晶体坐标启动显式局部优化，所有十个组合仍在 0.64 Å 内。起始配体的不可转动内部几何能否表达这个位置，取决于晶体扭转如何映射到等价原子。原有晶体诊断只使用第一个图同构。以 5FG4 为例，两个等价磺酰氧会改变用于复制二面角的参考向量。使用第一种映射，在相同查询内部几何上得到 2.27 Å 的投影，cognate 局部优化后为 9.29 Å；另一种等价映射得到 0.345 Å，优化后为 0.529 Å，在 4UYE 的交叉受体中为 0.620 Å。本轮诊断列出每个映射的**优化前**距离，按该距离选映射后再精修一次；全部十个组合的最终偏差都在 0.69 Å 内。旧诊断可能把参考坐标选择问题误判为配体内部几何不可达。这些晶体引导的运行独立于盲对接，首次实验证据不支持把它用作生产搜索初始化。[投影记录](results/projection_mappings.csv) · [局部优化结果](results/crystal_starts.csv)

## 复现

需要 Python 3.12、RDKit 2025.9.3、Meeko 0.7.1、Gemmi 0.7.4、PDBFixer 1.12.0、OpenMM 8.6.1、NumPy 2.3.5、SciPy 1.17.0 与 Biopython 1.88。以下命令从仓库根目录执行。把 CPU、Boost 前缀和官方 release 路径换成本机实际值。未启用 GPU。

```sh
D=benchmarks/search_backtracking/pose
B=/tmp/replay-backtracking
U=/tmp/replay-upstream
W=/tmp/replay-backtracking-poses
C=7
BOOST=/path/to/boost
mkdir -p "$W"
tar -xJf "$D/inputs.tar.xz" -C "$W"
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
python benchmarks/search_backtracking/core/build.py "$B" --boost-prefix "$BOOST" --jobs 2
python benchmarks/pair_gradient/build_arms.py "$U" --arms upstream upstream-counts --boost-prefix "$BOOST" --jobs 2
python "$D/experiment.py" "$W" --cpu "$C" --shared-core "$B/build/vina" \
  --upstream "$U/upstream/build/vina" --counts "$U/upstream-counts/build/vina" \
  --parent "$B/native/vina" --halving "$B/halving/vina" --quadratic "$B/quadratic/vina"
python "$D/summarize.py" "$W"
python "$D/verify.py" "$W" --parent "$B/native/vina"
python benchmarks/accepted_retention/core/build_stage.py "$B/tree" "$B/build" "$B/stages" --boost-prefix "$BOOST"
python "$D/search_stages.py" "$W" "$B/stages/vina" --cpu "$C"
```

使用 `python "$D/backtracking_panel.py" "$W"` 可从原始 PDB、CCD 重新制备。上游官方 release 与当前未修改生产二进制可通过 `parity.py "$W" --cpu "$C" --release /path/to/vina-1.2.7 --pristine /path/to/bcb5d83/vina` 独立比较。生产版本可在单独工作树检出 bcb5d83 后通过 `tools/build.py` 构建。已保存逐靶点、查询、受体、预算、姿势的全量结果及端到端/局部原始日志。[主实验坐标与日志](runs.tar.xz) · [分阶段坐标与日志](stage_runs.tar.xz) · [分组结果](results/summary.csv) · [逐个配对差值](results/contrasts.csv)
