# 盲放置、逐候选复算与归因

受体与配体输入、盒及盲提议约定见[冻结协议](MATCH_PROTOCOL.md)，实际盲函数在 [pocket_triangle.py](../../tools/pocket_triangle.py)。输入与一对一的轨迹源事件由 [pose_runs.py](pose_runs.py) 保存。每个原生构象对应两个各 1,024 项的候选库。每项的平移由重原子质心记录，旋转矩阵作用于以该源态重原子质心为原点的列向量。候选库、固定的腔点几何和源坐标即可按双精度重建所有 49,152 项。源配体做整体刚体换系后，全部重建坐标最大分量差 8.80×10⁻¹⁴ Å，两臂各首名序号完全一致。各库 1,024 个质心均不同，取向有 1,016–1,024 个（量化到 1e−8），不区分顶点顺序时有 277–1,002 个不同三角。

- [全候选变换](candidate_transforms.tar.xz) 包含两个受体的 sites*.npy 和按 source 命名的 npz。center 是整个配体重原子质心，matched_rotation 和 sham_rotation 是左乘列向量的 SO(3) 矩阵，site_triangles 是三位腔点编号，ligand_triangle 是有序输入重原子数组索引，shift 是 sham 的非零循环移位量。matched_score 与 sham_score 是无量纲几何粗分。所有数组索引从 0 开始，表中 coarse_rank 从 1 开始。
- [逐候选 RMSD/粗分/排名](results/candidate_metrics.csv.gz)、[逐事件汇总](results/candidate_stages.csv)和[输入给生产引擎的 48 个首选变换](results/selections.csv)按同一 source、arm、index 相连。RMSD 在受体框架中最小化图等价原子映射；其中一条独立链从实验 PDB 读取重原子、重新配准受体并使用配体图枚举映射，另用 RDKit CalcRMS 查验各组首名和最近者。逐项 RMSD 两种算法最大差 3.56×10⁻¹⁴ Å。[独立误差](results/candidate_audit.json)
- [建场](results/site_cost.csv)和[候选成本](results/generation_cost.csv)分列。两个盒各建一次，两者共 9,826 次最近重原子查询、2,053 次半径 8 Å 的计数查询，耗时 0.0430 秒。24 个事件各检验 29,754 或 29,760 个非共线有序三元组；每次仅对 1,024 个做三点拟合，24 次共享配准与构造共 0.1463 秒。每臂对 667,648 个逐重原子坐标查询近邻形状场。匹配粗排 0.4097 秒墙钟、0.4098 秒 CPU，打乱 0.3928 秒墙钟、0.3928 秒 CPU。开发 18 个相关源态分别收费。分段计时不含诊断、导入、压缩、驱动程序启动或源态采样。另起一个全新 Python 进程，含导入、建场、两臂评分、所有候选的压缩与写盘，24 个源态合计 1.5475 秒墙钟、1.5203 秒 CPU，264 份候选数组和 48 个选择逐项一致。[完整进程复跑成本](results/generation_replay.json)
- [全对应诊断](results/posthoc_edge_gate.csv)是冻结盲候选库之后的补充。事后用全部腔点有序对应做同样三点拟合，按原边长残差排名记录完整误差数组（[NPZ](results/posthoc_rmsd_by_edge_rank.npz)，以 Å 计，按 float32 压缩存储，汇总由双精度计算）。在 AT9283 开发的 18 个源态中，所有对应仍均未进入 2 Å，但 11 个能进入 3 Å。在固定的前 1,024 个对应中 11 个全部消失，最低 RMSD 的事件中位数由 2.957 升至 3.547 Å。[事后计算](results/posthoc_cost.csv)另耗 5.155 秒墙钟、5.145 秒 CPU；未进入盲预测成本。
- BRD9 82I 两个种子的匹配候选分别有 18、13 个 ≤2 Å，最好粗排名为 8、34；打乱候选有 9、4 个，最好粗排名为 28、38。取这些已知近姿中粗排最高项的[四个晶体知情干预](results/posthoc_selections.csv)仅用于追踪后续捕获及保留。正式预测仍是各臂几何全局首名。

用仓库根目录执行下面的命令。W 是按主复现说明解开并生成 sources.csv 的绝对工作目录，包含 aurora、brd、sources 及原始输入。预测函数不读参考坐标。先运行 generate 并锁定全部选择，再评价；posthoc_selection 仅供知情分叉。

~~~sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
D=benchmarks/cavity_triplet
python "$D/generate.py" "$W/sources.csv" "$W/candidates"
python "$D/replay_generation.py" "$W" "$W/generation_replay.json"
python "$D/evaluate_candidates.py" "$W/sources.csv" "$W/candidates" "$W/aurora" "$W/brd" "$W/candidate_analysis"
python "$D/edge_gate.py" "$W" "$W/candidate_analysis"
python "$D/posthoc_selection.py" "$W"
~~~
