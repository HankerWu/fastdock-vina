# 三键邻域修复的真实结构对照

巴瑞替尼的四元氮杂环使深度优先遍历漏掉六个相隔三条键的原子对。逐层遍历把它们从非键能量中移除。初次盲对接中，最明显的变化是 JAK2 的一个近晶体姿势从保存的搜索池中消失。另一个预设种子下，父版和修复版都未找到 2 Å 以内的姿势。该实验支持图邻域修复的化学正确性，没有显示姿势成功率提高。

## 哪些结果发生了变化

首轮使用种子 817139，一个线程、一条轨迹，24 Å 立方盒，最多输出九个候选。所有输入都由同一套 CCD 化学图准备，搜索起点使用 ETKDGv3 加单构象 MMFF 优化。已知共晶配体只提供盒中心及评价所需的独立参考。将输出坐标保持在受体坐标系中，按化学图对称性和手性允许的映射取重原子 RMSD 的最小值。下表列出全部三个受影响组合，单位为 Å。

| 查询 → 受体 | 上限 | 首姿，父版 → 修复 | 输出池最佳，父版 → 修复 | 总评估，父版 → 修复 |
|---|---:|---:|---:|---:|
| JAK2，6VN8 → 6VN8 | 3,000 | 9.025 → 10.333 | 8.352 → 7.645 | 3,308 → 3,234 |
| JAK2，6VN8 → 6VN8 | 12,000 | 1.821 → 9.156 | 1.821 → 7.645 | 12,423 → 12,217 |
| JAK2，6VN8 → 6VGL | 3,000 | 10.519 → 9.665 | 10.519 → 9.574 | 3,183 → 3,284 |
| JAK2，6VN8 → 6VGL | 12,000 | 5.941 → 7.933 | 5.941 → 7.933 | 12,309 → 12,316 |
| BMP2K，4W9X → 5I3O | 3,000 | 8.346 → 8.813 | 8.157 → 8.813 | 3,133 → 3,177 |
| BMP2K，4W9X → 5I3O | 12,000 | 9.198 → 8.767 | 8.199 → 8.229 | 12,303 → 12,207 |

这三组共享巴瑞替尼。JAK2 交叉对接对齐 60 个口袋 Cα，残差 0.428 Å；BMP2K 为 54 个和 0.894 Å。BMP2K 4W9X 含两个工程突变，使用其配体参考及 5I3O 的受体。主实验单独替换广度优先遍历，配对时间只反映该项改动。另七组输入包括 CSF1R 的另一种氮杂环丁烷配体、BTK、CRABP2、FXR、另一种 BMP2K 配体及 JAK2–ruxolitinib 的原位与交叉对照。这七组在两个上限下的父版和修复版输出逐字节相同，每阶段计数相同。CSF1R 的 3,000 上限原位首姿在两臂均有 14 对重原子触及粗碰撞阈值。

对全部 20 个输入–预算单元，原生、父版和修复版的 Top-1、Top-5 及输出池 Oracle 命中数分别为 1、1、0。条件排名成功率在原生和父版为 1/1，修复版无可作分母的成功池。原生、父版和修复版联合能量–梯度调用分别合计 157,321、156,926 和 156,702 次，主运行实测端到端时间合计 20.962、3.635 和 3.598 秒。该计时包含准备网格、搜索和精修；之前两版实现之间已有的速度差异不能归因于这次修复。每次计时仅测一次，不据此声称加速。首轮有三项修复版的实际评估稍高；它们按 0.85 倍上限复测，实测次数低于对应父版，三项仍为零命中。输入、参数、完整阶段计数、运行时间和每个姿势见 [cells.csv](cells.csv)、[modes.csv](modes.csv) 和 [cost.csv](cost.csv)。

新种子 191009 对全部六个受影响单元进行父版、修复版配对，两臂均为 0/6，三个姿势指标一致。JAK2 原位 12,000 上限的首姿由 3.14 变为 2.77 Å，池最佳由 2.94 变为 2.77 Å。JAK2 交叉 12,000 上限的首姿由 3.56 变为 9.93 Å。BMP2K 同上限首姿由 10.56 变为 8.30 Å。两个预算、两个种子和两个受体中的相同配体不能算作独立药物样本。[确认运行](confirmation.csv)

## 近晶体姿势消失在哪一步

- 冻结首轮父版的六个完整候选池，共 29 个姿势。保持原子类型、可转动键、坐标和统一零游离态参考不变，父版和修复版均用显式受体重新打分。六个池的首选身份保持一致。6VN8 原位 12,000 上限的 1.821 Å 姿势在修复版中仍居首位；前两位相差仅 0.005 kcal/mol，使用命令行输出精度。[固定池结果](frozen_pool.csv)
- 在两套种子和三个受影响组合的 24 次运行上被动输出保存的蒙特卡罗候选、精修输入与精修结果。复跑的最终坐标文件和每阶段评估/调用数全部与正式运行相同。6VN8 原位首轮 12,000 上限下，父版保存池最佳为 1.819 Å，精修前为 1.819 Å，精修后为 1.821 Å。修复版对应为 7.647、7.647 和 7.645 Å。差异已出现在保存搜索池的时刻。日志没有记录所有被拒绝的瞬态提议，因此仍可能存在短暂访问但没有保存的情况。[逐阶段最佳姿势](stage_oracles.csv)
- 共晶坐标起步的局部精修另算。父版和修复版在 6VN8 原位分别收敛到 0.220 和 0.277 Å。保留盲输入的刚性几何，只向其投影实验扭转并按初始几何选择等价映射，两版均从 0.939 Å 收敛到 1.023 Å。晶体引导的三个组合都保留了小于 2 Å 的精修终点，这些额外信息仅用于判断局部盆地是否存在。[父版](guided_parent.csv) · [修复版](guided_fixed.csv)

## 输入、排除及复算

11 座输入结构有七个受体通过预先采用的质控。2XP2、3DCT、6WTO、4W9X 首次受体准备出现 Meeko 问题。沿用仓库的远端缺失侧链处理后，2XP2、6WTO 和 4W9X 仍准备失败；3DCT 会在距配体 4.088 Å 处重建原子。6WTO、4W9X 的第二种准备亦分别涉及距配体 6.753、6.003 Å 的修复。全部排除作为盲输入受体，4W9X 的实验配体参考可用于交叉评价。FXR 两结构先按 DBREF 校正为 UniProt Q96RI1 的共同编号，53 个口袋 Cα 配准残差 2.135 Å，超过既有协议 2 Å 上限，排除该交叉方向。[数据来源、分辨率、SMILES、骨架、最近遗漏组分、首轮准备耗时和输入 SHA-256](sources.json)

原始 PDB、CCD、输入 PDBQT、解析映射、准备失败日志、所有计数和姿势都在 observations.tar.xz。对上游源码和上游被动计数版本逐次核对坐标。六个受影响单元另用官方 1.2.7 Linux release 逐字节验证。审计脚本独立从模型坐标与原始日志复算主实验、低成本复测及第二种子合计 75 次运行、315 个姿势，最大 RMSD 差异小于 6×10⁻¹⁵ Å。[二进制哈希](binaries.json) · [官方版核对](replay.csv) · [审计结果](verification.json)

另外用提交 7b9bf89 的完整生产版本复跑了两个种子的全部 26 项修复臂。临近闭环编号的预索引与单独遍历修复输出逐字节相同，六种阶段计数逐项相等；生产版本的单次时间单列于 [production.csv](production.csv)。

从仓库根目录执行。归档输入可以直接使用；从头制备需 RDKit、Meeko、PDBFixer、OpenMM、Gemmi、NumPy 及 SciPy。主实验使用 GCC 14.2、C++14、-O3、Boost 1.83 和 AMD EPYC 9V74。编译时若使用自定义 Boost 安装，给各构建命令提供相同的 --boost-prefix。可通过系统工具 taskset 固定单个可用 CPU 核心。

```sh
P=benchmarks/bond_distance/pose
W=$(mktemp -d)
tar -xJf "$P/observations.tar.xz" -C "$W"
python "$P/audit.py" "$W"
git worktree add --detach "$W/parent-src" 36636171ae6aef42471b0de5d6b7a90b5f90732e
python "$W/parent-src/tools/build.py" --output-dir "$W/parent-build"
git worktree add --detach "$W/fixed-src" 36636171ae6aef42471b0de5d6b7a90b5f90732e
git -C "$W/fixed-src" apply "$PWD/benchmarks/bond_distance/preparation/bfs_only.patch"
python "$W/fixed-src/tools/build.py" --output-dir "$W/fixed-build"
git worktree add --detach "$W/production-src" 7b9bf89
python "$W/production-src/tools/build.py" --output-dir "$W/production-build"
python benchmarks/pair_gradient/build_arms.py "$W/up" --arms upstream upstream-counts
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
python "$P/panel.py" run "$W" --repo . --parent "$W/parent-build/vina" --fixed "$W/fixed-build/vina" \
  --upstream "$W/up/upstream/build/vina" --counts "$W/up/upstream-counts/build/vina"
python "$P/panel.py" confirm "$W" --repo . --parent "$W/parent-build/vina" --fixed "$W/fixed-build/vina"
python "$P/panel.py" production "$W" --repo . --fixed "$W/production-build/vina"
python "$P/panel.py" cost "$W" --repo . --fixed "$W/fixed-build/vina"
python "$P/frozen_pool.py" "$W" "$W/parent-build/vina" "$W/fixed-build/vina"
python "$P/guided.py" . "$W" "$W/parent-build/vina" "$W/fixed-build/vina"
python benchmarks/accepted_retention/core/build_stage.py "$W/parent-src" "$W/parent-build" "$W/parent-observer"
python benchmarks/accepted_retention/core/build_stage.py "$W/fixed-src" "$W/fixed-build" "$W/fixed-observer"
python "$P/stages.py" . "$W" "$W/parent-observer/vina" "$W/fixed-observer/vina"
python "$P/audit.py" "$W"
```

下载官方 1.2.7 可执行文件后，以 panel.py replay 单独核对六个主实验受影响输出，参数为工作目录、--repo . 和 --release 加可执行文件路径。从 RCSB 与 CCD 重建输入时，使用 panel.py prepare 加工作目录与 --repo .。若要复核四个失败受体的第二种准备，使用仓库 benchmarks/torsion_kinase_affinity/prepare_receptors.py 中的 prepare(work, pdb, component)，按 sources.json 的 PDB 和 CCD 逐一运行。
