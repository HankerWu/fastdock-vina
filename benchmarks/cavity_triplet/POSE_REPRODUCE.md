# 在真实 Monte Carlo 轨迹上重做放置、局优和发布

从仓库根目录执行。W 是新的空目录，使用绝对路径。原实验使用单线程 Linux、Python 3.12.14、NumPy、SciPy、RDKit、Gemmi 和 C++14。以源树哈希 `6c6f65b0b7445cb364415f2079906e0d67e015bb` 为生产对照。补丁均作用于独立研究构建，保留源提议和随机数序列，在指定轨迹的指定步仅改变六个刚体自由度。控制组、剩余轨迹与干预后的求解器均按各自原程序运行。未插桩生产 Pro 用当前源树的 tools/build.py 构建。OFFICIAL 是 Vina 1.2.7 Linux 发行版，COUNTER 是上游 3c65c0b3 的计数版本，构建方法见 [pair_gradient](../pair_gradient/build_arms.py)。实验用输入及四个二进制的散列见 [pose_inputs.json](results/pose_inputs.json)。

```bash
set -eu
D=$PWD/benchmarks/cavity_triplet
A=$PWD/benchmarks/receptor_linkage/aurora
O=$A/internal_loops
mkdir -p "$W/aurora" "$W/brd" "$W/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 TMPDIR="$W/tmp"
tar xf benchmarks/torsion_kinase_affinity/evidence.tar.xz -C "$W/aurora" raw/2X6E.pdb raw/2J4Z.pdb prepared/2X6E-protein.pdb prepared/2J4Z-protein.pdb
tar xf "$A/docking.tar.xz" -C "$W/aurora" prepared
tar xf "$A/references.tar.xz" -C "$W/aurora"
tar xf "$A/probe_raw.tar.xz" -C "$W/aurora" raw/at9283_chemical.bonds raw/danusertib_native.bonds
tar xf "$A/state_windows/coordinates.tar.xz" -C "$W/aurora"
tar xf "$O/evidence.tar.xz" -C "$W/aurora"
tar xf benchmarks/torsion_proposals/fixtures/inputs.tar.xz -C "$W/brd"
sh "$D/build_observer.sh" "$W/build" --boost-prefix "$BOOST_PREFIX"
python tools/build.py --output-dir "$W/production" --boost-prefix "$BOOST_PREFIX" --jobs 1
# COUNTER 可直接用固定上游源码编译。若本地已有相同构建，可略过该命令。
python benchmarks/pair_gradient/build_arms.py "$W/upstream_build" --arms upstream-counts --boost-prefix "$BOOST_PREFIX" --jobs 1
COUNTER=$W/upstream_build/upstream-counts/build/vina
BIN=$W/build/bin/vina
python "$D/pose_runs.py" sources "$W" --observer "$BIN"
python "$D/pose_runs.py" controls "$W" --observer "$BIN" --official "$OFFICIAL" --counter "$COUNTER" --pro "$W/production/vina"
python "$D/generate.py" "$W/sources.csv" "$W/candidates"
python "$D/pose_runs.py" branches "$W" --observer "$BIN"
python "$D/pose_runs.py" timing "$W" --observer "$BIN"
python "$D/pose_metrics.py" "$W"
python "$D/evaluate_candidates.py" "$W/sources.csv" "$W/candidates" "$W/aurora" "$W/brd" "$W/candidate_analysis"
python "$D/independent_pose.py" "$W" "$W/independent_pose.json"
cmp "$D/results/modes.csv" "$W/modes.csv"
cmp "$D/results/stages.csv" "$W/stages.csv"
cmp "$D/results/forks.csv" "$W/forks.csv"
cmp "$D/results/qa.csv" "$W/qa.csv"
cmp "$D/results/independent_pose.json" "$W/independent_pose.json"
```

`sources` 从既有完整图数据取 18 个按步号锁定的开发事件。六次独立种子的原生搜索各取第一条轨迹变异后的第一个构象，并继续到正常发布，以检查它们与干净生产二进制输出逐字节相同。这次源态提取共用了一份单独的完整搜索。公布的原型用两次进程和预先计算的放置实现一次干预，实际运行须另计采集与建场；这两个阶段没有调用晶体评价。`branches` 写出所有保存事件、投池坐标和末端快照。`timing` 对同一放置关闭这些快照并逐文件校验其发布输出。各进程计时从启动前开始，包含引擎内部建图与最后输出。姿势分布表的计时来自有阶段记录的运行。逐档成本见 [源态](results/source_cost.csv)、[候选](results/generation_cost.csv)、[观察运行](results/branch_cost.csv)、[关快照重放](results/prediction_cost.csv)、[外部基线](results/control_cost.csv)、[官方额外计数重放](results/counter_cost.csv)。控制组没有候选生成。正式实验结论用无晶体信息的粗排首名。

每一份 `.history` 逐行记录实际到达步、经验分、容量、排斥半径（Å）和操作后的全部有序成员 ID。每个 ID 是首次出现该精确构象的步号。首次淘汰以该 ID 离开容器为准。几何相似的后续构象可继续入池，分别报告直接身份和几何阈值。原 ID 的末端衔接使用双精度的已保存坐标和实际精修索引，最后与 PDBQT 按每个坐标 0.000501 Å 容差匹配。空直接排名表示未在该阶段追到同一构象。[逐次干预](results/forks.csv) · [每次入池后的列表](results/residency.csv) · [全部中间坐标与输出](pose_raw.tar.xz)

事后 BRD9 诊断在每个既有候选库中取粗排最前的 RMSD ≤2 Å 提议，独立记录两枚种子、两臂的四个干预，几何标签仅用于该诊断。按以下命令接续重跑。

```bash
python "$D/posthoc_selection.py" "$W"
python "$D/pose_runs.py" posthoc "$W" --observer "$BIN"
python "$D/posthoc_metrics.py" "$W"
cmp "$D/results/posthoc_branches.csv" "$W/posthoc.csv"
```

发布的归档包含源构象 NPZ、工作进程标准输出、原始 PDBQT、单次注入参数、全容器成员列表与各阶段双精度重原子坐标。每个原始文件有 [SHA-256 清单](results/pose_raw_sha256.json)。[候选变换](candidate_transforms.tar.xz) 保存选中的全部三点对应、刚体变换和粗分。使用生产输入重新生成路径索引，或按 [相对路径清单](results/pose_sources.csv) 重定向本地文件。官方计数重放在本面板 12 份输出与官方发行版全部相同。实验开始时系统临时盘已满，编译从开始就使用 W/tmp。本轮 48 个盲分叉与 4 个事后分叉均正常进程退出。首次生成便携输入索引时混用相对与绝对路径，改为相同绝对基准后重新导出，全部原生文件及已有散列未改变。

最后在空构建目录应用这些补丁，单线程全量编译耗时 56.203 墙秒、56.022 CPU 秒。用新程序重跑 48 个盲分叉及 4 个事后分叉，耗时 14.694 墙秒、14.591 CPU 秒、2,614,300 次主搜索目标与梯度评估。52 份发布池、52 份阶段坐标流、52 份历史、52 份注入日志及 52 份注入参数均与首次运行逐字节相同。逐项从压缩档读取并重算了 576 个原始文件和源构象文件的 SHA-256。[重编与重放](results/pose_clean_replay.json)
