# 搜索线搜索能否改善同一激酶内的 Kd 排序

在两个人源野生型激酶的同测定配对中，更稳定的局部线搜索给出了方向相反的结果。p38α 的第一枚种子在 12,000 搜索上限下从 16/37 提高到 24/37；用同样输入和上限，第二枚种子从 24/37 降到 9/37。EGFR 的相应变化是 8→2/24 和 3→3/24。第一枚种子上仅保护失败末态、继续用二分法，在 p38α 得到的 29/37 还高于二次插值。没有据此更换默认搜索或亲和力模型。[运行前协议](PROTOCOL.md) · [研究实现](../core/README.md)

## 数据与条件

从 Davis 等原论文补充表 1 确认非突变实验构建体，从补充表 4 提取以 nM 计的精确阳性 Kd。EGFR 在此前已核对身份的 13 个抑制剂中有 8 个实测值。p38α 在该组只有 sorafenib 与 dasatinib，因此在对接前登记扩充至该实验行的全部 10 个可用值。补充的 8 个药名按 PubChem 与公开 Davis 化学图逐一核对 InChIKey。两靶点有一个分子重合，合计 17 个不同分子。空白留空，每个靶点只比较同一行、相差至少三倍的配对。保留原两分子子集、去共晶子集、相同骨架、等大小等电荷和按大小排序对照。每个靶点都只有一对同时满足等重原子数、等形式电荷和三倍差；合格对中没有相同 Murcko 骨架的两分子。[Davis 原始论文及补充表](https://www.nature.com/articles/nbt.1990) · [精确来源与哈希](sources.json) · [逐分子结构、骨架与标签](panel.csv)

EGFR 用 1M17，p38α 用 1A9U。每个刚性受体保留蛋白 A 链的空白/A 替代位置，以各自沉积配体的重原子质心定义 26 Å 盒子。EGFR 的 Cys751 与 Asp831 在 8 Å 内还有未使用的 B 替代位置；两个口袋 8 Å 内均未补建重原子。省略水，只用一套固定化学微态；在多环系、碱性胺和受体诱导构象上均有物理局限。原有配体直接复用准备文件，新增者每个用相同八构象 ETKDGv3/MMFF 制备并选定一个收敛起点，全部实验臂共享，未调用查询晶体坐标。整个公共准备过程在源文件已下载后用时 13.93 秒。[EGFR 1M17](https://www.rcsb.org/structure/1M17) · [p38α 1A9U](https://www.rcsb.org/structure/1A9U) · [准备原始记录](preparation.json) · [固定输入](inputs.tar.xz)

各臂使用相同一条轨迹、一线程、种子和最多九个输出名额。顺序在输入区组内轮换，每个运行全新进程，计时含建图、显式精修和写出。表中的积分为强结合分子的原生首姿分数较低即得一分，同分得半分。Spearman 与 Kendall tau-b 对全部精确 Kd 及预测分数排序。差值误差使用 298.15 K 作为换算参考温度，按强弱对计算预测差与 RT ln(Kd弱/Kd强) 的绝对差。各对共享分子；附逐分子剔除结果，未将种子、受体或分子对当作独立靶点。[首轮全部配对](results/pairs.csv) · [逐分子敏感性](results/leave_one_out.csv) · [全部坐标、标准输出、计数与运行元数据](runs.tar.xz)

## 逐靶点结果

| 靶点 | 种子 | 名义预算 | 上游 Vina | 当前 Pro | 保护二分 | 保护二次 |
|---|---:|---:|---:|---:|---:|---:|
| EGFR，24 对 | 458921 | 3,000 | 5 | 6 | 7 | 9 |
| EGFR | 458921 | 12,000 | 8 | 8 | 9 | 2 |
| p38α，37 对 | 458921 | 3,000 | 26 | 14 | 24 | 26 |
| p38α | 458921 | 12,000 | 16 | 16 | 29 | 24 |
| EGFR | 703109 | 10,800 / 12,000 | 未重跑 | 3 / 3 | 7 / 7 | 3 / 3 |
| p38α | 703109 | 10,800 / 12,000 | 未重跑 | 24 / 24 | 24 / 24 | 9 / 9 |

首枚种子 12,000 档，EGFR 的 Pro / 二分 / 二次 Spearman 为 −0.333 / −0.357 / −0.881，Kendall 为 −0.286 / −0.214 / −0.714。p38α 相应为 −0.236 / 0.697 / 0.345 和 −0.200 / 0.556 / 0.244。仅按重原子数比较在两个靶点均得 12.5 分，两个分母分别为 24 与 37。去掉各自共晶药物后，EGFR 为 6 / 8 / 2（分母 18），p38α 为 12 / 22 / 17（分母 29）。p38α 的 BIRB-796 在本组实测 Kd 最低，为 0.45 nM，去掉它的九条相关比较后，该靶点为 16 / 20 / 15（分母 28）。二次法第一枚种子的净收益依赖这一个分子；第二枚种子改变了 19 个配对方向，其中 2 个转对、17 个转错。[首轮分层结果](results/summary.csv) · [首轮相关性与成本](results/cost_ranks.csv) · [第二种子的全部配对](cost_check/pairs.csv) · [第二种子的相关性与成本](cost_check/cost_ranks.csv)

## 预算、碰撞和固定坐标对照

同为 12,000 上限的首轮，p38α 的 Pro / 二分 / 二次实际联合能量–梯度求值为 123,504 / 123,231 / 124,560，十次端到端运行合计 2.513 / 2.590 / 2.559 秒。EGFR 相应为 99,381 / 99,491 / 99,868 与 1.695 / 1.803 / 1.767 秒。封顶检查发生在搜索外循环；末端显式精修也有实际花费。同一名义上限只限定预算范围。端到端时间另含非梯度打分与预处理，单进程结果含宿主机负载噪声。[逐次计数与时间](results/runs.csv)

第二枚种子的全因子验证同时运行 10,800、12,000 两档，避免把单纯缩短预算的影响算进算法。在 p38α，短预算二分为 112,015 次和 2.493 秒，低于全预算父版本的 123,894 次和 2.578 秒，都是 24/37；短预算二次为 111,885 次、2.632 秒、9/37。EGFR 的短预算二分为 89,737 次、1.793 秒、7/24，全预算父版本为 98,733 次、1.679 秒、3/24。短预算二次为 89,527 次、1.926 秒、3/24。该实验中没有一个改动在两靶点上同时满足准确率提高、求值数降低、实测时间降低。[全因子逐次记录](cost_check/runs.csv)

固定首轮 12,000 档父版本的整个返回池，18 个受体–配体组合共 102 个构象。原始 PDBQT 保持拓扑，读取每个候选的打印坐标，固定相同的游离参考能量零点。三个 Pro 臂在 102 个坐标的十项可见计分分解逐项相同（CLI 打印精度）。搜索所改变的首姿及其游离参考可改变分子间顺序。固定坐标单点评估不计入对接预算。[固定坐标与分量数值](results/fixed_scores.csv)

排序方向也会掩盖非物理首姿。二次法在 3,000 档给 p38α–nilotinib 的首姿报 +166.128 kcal/mol，18 个配体重原子距受体重原子不足 2 Å，最近 1.017 Å。该姿势的受体相互作用项为 +234.110 kcal/mol。在第二枚种子中，保护二分法的 EGFR–erlotinib 首姿也出现 1.340 Å 的最近重原子距离。这些原始值均保留在配对评价中。[全体首姿的最短距离](results/clashes.csv) · [从原始表、结构和输出的独立复算](independent_check.json)

## 复现

在仓库根目录运行。重算已归档的全部主指标只需要 Python、SciPy 与输入输出包。用存入 Git 的研究源码重新构建可执行文件并从源数据重建输入后，第二枚种子的 108 个完整返回池逐字节相同，搜索、精修计数和首姿分数也全部重现。[独立重放记录](replay.json)

重建制备用 Python 3.12、RDKit、Meeko、OpenMM、PDBFixer、xlrd，版本与结构面板一致。构建可指定机器上安装的 Boost，复算时替换可用 CPU 编号。

~~~bash
D=benchmarks/search_backtracking/assay
W=$(mktemp -d)
tar xJf "$D/inputs.tar.xz" -C "$W"
tar xJf "$D/runs.tar.xz" -C "$W"
python "$D/analyze.py" "$W" "$W/primary_summary"
python "$D/analyze.py" "$W" "$W/cost_summary" --cost-check
cmp "$W/primary_summary/pairs.csv" "$D/results/pairs.csv"
cmp "$W/cost_summary/pairs.csv" "$D/cost_check/pairs.csv"
~~~

从原始论文重建输入时，prepare_assay.py 会校验原始文件摘要。重新执行预测时用新目录存放输出；每个对照上游运行另用只加计数器的上游二进制回放，程序逐字节核对输出才接收计数。

~~~bash
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
python "$D/prepare_assay.py" "$W"
python "$D/independent_check.py" "$W" "$W/independent_check.json"
cmp "$W/independent_check.json" "$D/independent_check.json"
B="$W/research"
U="$W/upstream-builds"
python benchmarks/search_backtracking/core/build.py "$B" --jobs 3
python benchmarks/pair_gradient/build_arms.py "$U" --arms upstream upstream-counts --jobs 3
C=$(python -c 'import os; print(min(os.sched_getaffinity(0)))')
args=(--cpu "$C" --upstream "$U/upstream/build/vina" --upstream-counts "$U/upstream-counts/build/vina" --parent "$B/native/vina" --halving "$B/halving/vina" --quadratic "$B/quadratic/vina")
python "$D/assay_run.py" "$W" "${args[@]}"
python "$D/assay_run.py" "$W" "${args[@]}" --cost-check
python "$D/fixed_pool.py" "$W" "$W/primary_summary" "$B/native/vina" "$B/halving/vina" "$B/quadratic/vina"
~~~
