# 氢专属扭转应当继续精修，还是重新分配？

在普通 Vina 势下，单独转动末端羟基氢不会改变重原子位置、原子类型或起点能量。搜索仍然选择这个扭转、从当前姿势重新启动 BFGS，再执行接受和保存。我们在真实 C++ Monte Carlo 循环中把这种提议的原有随机数转给其他扭转。无条件版本在 RXRα 和 PTP1B 的个别轨迹找到了此前未进入输出池的正确姿势。仅在接受态最后一次经过 hunt 势优化时才启用的版本，把这些新增命中全部漏掉。实验从 3d9dfb9 出发，默认搜索继续沿用父版。[预登记](PROTOCOL.md) · [C++ 实现及冻结观察器](proposal.patch) · [原生 Vina 对氢位置的说明](https://autodock-vina.readthedocs.io/en/latest/faq.html)

## 提议和因果对照

解析器把轴端的重原子存入父节点。对没有子节点且其自身覆盖原子都是氢的叶节点，Vina XS 的配体–受体项及活动内部原子对项均不依赖该角；轴端重原子不会随自己的轴转动。输入中的氢继续决定供受体分型。实验局限于普通 Vina 刚性受体单配体对接。AD4 对极性氢有真实地图项，不能套用这一推导。自动分类仅识别这个充分条件，未把绕轴对称或偶然共线的其他扭转判成零自由度。

原实体抽样和均匀角保留。如抽中氢专属转子，设已有均匀角对应的分位数为 u，其他扭转数为 K。每次替换臂用 floor(Ku) 选择编号，将 Ku 的小数部分线性映射回 [−π,π) 作为新角；K 为零时保留原提议。连续均匀变量的这两个输出分别在编号和角度上均匀，实际随机源有有限精度。条件臂在初始构象或刚完成保存用第二次局优时保留原提议，保留重新进入 hunt 势的机会；拒绝提议不改变已接受态的来源。各臂继续使用相同目标、梯度、接受温度、名额和最终排序。角度重分配只在原操作位置发生；轨迹分岔后接受分支可消耗不同随机数。带局优和历史门控的对接轨迹用于寻找极小值，不能解释成平衡系综采样。条件臂付出了拓扑分类和检测成本，门控零触发也可能改变耗时。

在原有 ESR1/PTP1B/PARP1 结构面板的 28 条高预算重放（两臂各 14）中，23 条出现首次替换。对每次替换，在同一传入构象及原始抽样下分别执行两种提议和全新的同配置 BFGS，接受检验使用同一生成器快照。原分支的重原子位移在 23 次提议中均为 0。即便如此，原分支的 hunt 终点与接受态的重原子 RMS 位移均值为 0.718 Å，最大 12.950 Å，包括从随机初态开始的事件。完整输出均逐字节复现。另在独立面板的七个触发事件中，原分支局优用了 26–46 次求值，移动 0.0013–0.3482 Å。因此不能把这些调用的全部成本都解释为零效用。[冻结状态和局优记录](replay/stages.csv) · [独立首分歧](transfer/first_hunt.csv)

## 已知口袋的结构实验

每格一条轨迹、一颗 CPU，输入、受体、最多九个模式、24 Å 盒以及 3000/12000 两档预算配对。首姿成功采用受体坐标系中对称性校正的重原子 RMSD ≤2 Å。两个面板内四个主臂的 Top-5 与输出池 Oracle 命中数各自均等于 Top-1，正确模式留在输出时均排首位。重复种子、预算和互换受体的条件共享分子。

| 面板 | 条件 | 官方原生 | Pro 父版 | 每次替换 | 条件替换 | 父版与每次替换事后择优 |
|---|---:|---:|---:|---:|---:|---:|
| 既有结构，5 个配体 | 28 | 1 | 1 | 3 | 1 | 3 |
| 冻结结构，4 个配体 | 24 | 5 | 4 | 5 | 4 | 5 |

按首姿改善或恶化至少 0.5 Å 计，开发面板每次替换的改善/恶化/其余为 10/5/13，条件替换为 2/2/24。原来 10 格改善中有 1 格在同格达到 0.5 Å 改善，原来 5 格恶化中有 3 格避开。冻结面板对应为 4/1/19 和 0/0/24。开发组的两个新增二值命中共享一个 PTP1B 配体、受体和种子，仅预算不同。独立组的一个新增二值命中来自 RXRα；原生在另一 RARγ 轨迹的正确姿势没有被找回。GSK3β 的 16 个阴性条件覆盖没有可重分配的重扭转、或没有氢专属扭转两种情况，两个候选均与父版逐字节相同。两种视黄酸使用共享的中性羧酸微态。该输入不能检验生理环境中离子态能否复现同一结果。[逐结构、预算、成本](structure_summary.csv) · [开发面板全部名次](primary/modes.csv) · [独立结构、分辨率和完整成本](transfer/README.md)

RXRα–1FBY 种子 51871 的 12000 档，父版首姿 5.414 Å，每次替换 0.700 Å，条件替换 5.414 Å。每次替换降至 10200 档仍为 0.700 Å，12,306→10,352 次实际联合求值，五次交错计时中位数 0.1722→0.1515 秒。第一次替换的姿势反而较远且被拒。到第 118 步，一次普通重扭转提议由 7.379 Å 经 hunt 到 1.576 Å，再经第二次优化以 1.496 Å 入池。相同种子的 3000 档是严格前缀，第 44 步的 2.569 Å hunt 被拒，最后首姿 5.206→7.934 Å。独立组八个可激活条件先逐格取五次中位再合计，父版/每次替换/条件替换为 1.282/1.288/1.328 秒。单格的等成本收益未迁移到全组。[五次计时](transfer/timings.csv)

PTP1B 的 1T49 配体转至 1T48、种子 202711，第 5 步普通重扭转的 hunt 终点及接受态仍为 2.762 Å，第二次局优才将它推进至 1.545 Å 并入池。高预算的首姿从 8.957 Å 降至约 1.414 Å。其后一个 0.850 Å 接受态在九个名额已满时没触发第二次局优；显式精修产生的 1.150 Å 姿势与最终首姿仅相距 0.940 Å，随后被 1 Å 输出间距过滤。由此可以分别看到正确盆地的进入、一个更佳候选的漏存，以及另一个更佳候选在精修后的去重。父版和条件版都没有在这些采样阶段进入 2 Å。低预算不等于新的独立配体。全过程双精度坐标与各阶段的判断见[12 次完整轨迹](path/README.md)。

## Galectin-3 的同测定 Kd

47 个化合物来自 Sörme 等人的同一荧光偏振实验，结构比较和自由能比较使用相同输入、评分与名额。本系列曾在早期 FastDock 研发中使用，本轮作为跨结构任务的同靶实测检验。每个受体取实际输出首姿，再对 1KJL/3T1L 的分数取均值。403 个至少三倍 Kd 差的配对按正确方向计 1、平局计 0.5，35 对还匹配重原子数和输入形式电荷。受体分别的结果、全部不等 Kd 配对、每对实验 ln(Kd弱/Kd强)、预测分差、删单分子影响和分子级配对自助法均归档。[原始论文](https://doi.org/10.1021/ja043475p) · [图 4 筛选和来源](../pair_gradient/galectin_external/README.md) · [追加预算和种子之前冻结的协议](AFFINITY_CONFIRM.md)

| 种子、上限 | 父版 / 每次替换 / 条件替换，403 对 | 同顺序，35 对 | 父版 / 每次替换 Spearman |
|---|---:|---:|---:|
| 214039，3000 | 116 / 202.5 / 116 | 11 / 19 / 11 | −0.359 / 0.008 |
| 751483，3000 | 127 / 209 / 127 | 8 / 15 / 8 | −0.303 / 0.023 |
| 214039，12000 | 193 / 220 / 159 | 13 / 18 / 9 | −0.077 / 0.061 |
| 751483，12000 | 223 / 216 / 187 | 20 / 27 / 13 | 0.057 / 0.020 |

首行的上游原生也是 116/403。单看 1KJL 与 3T1L，父版→每次替换分别为 146→213、119→174。把每次替换降额到 2550，均分结果仍为 182/403，联合求值 269,460（父版 310,279），wall 11.51 秒（父版 11.62），CPU 11.32 秒（父版 11.51）。这次计时是 94 个成功进程的时长之和，不能分辨百分之一量级的稳定加速。官方原生在首轮 94 个格子的进程峰值最大 1,865,412 KiB，资源竞争导致一个并发上游作业被系统终止，随后串行完成上游对照并核对计数重放的输出字节。失败尝试不计入单次成功预测。

原种子低预算时条件臂遇到 744 次氢专属提议，替换 0 次。新种子低预算为 546 次遇到、3 次替换。提高预算后分别替换 2194 和 1981 次，准确度分别比父版低 34 和 36 个配对分。这个门控既漏掉低预算中的机会，也没在高预算中保障改进。低预算的每次替换从负相关移至接近零；仅按重原子数预测大分子更强已达 259.5/403。两批低预算中删去任意单个化合物，每次替换相对父版仍至少净增 60.5 和 59 分；按整个分子配对自助，第二种子的分差区间覆盖零。相同分子和受体重复使用，不可把 403 对当作独立样本。原种子 94 个首姿中的严重蛋白–配体重叠由 0 个姿势增至 2 个，新种子的总碰撞对从 155 降到 82。碰撞减少不能统一解释两个低预算结果。

将 94 个父版 Gal3 首姿及全部 28 个开发结构父版输出池（114 个模式）逐一固定坐标，用相同游离态参照 0 重算。在共 208 个姿势上，父版、每次替换及条件替换所输出的总分和所有能量分项逐项完全一致，精度为 CLI 的 0.001 kcal/mol。研究实现仅替换搜索提议；亲和力排序的变化来自找到不同姿势及其常规定标，未增加可迁移的亲和力物理项。[固定坐标](fixed_pool.csv) · [低预算结果](affinity/primary/summary.csv) · [新种子高预算结果](affinity/confirm_751483-12000/summary.csv)

## 跳过无贡献的氢叶运动学

另一项实现保留每一次原提议和全部 BFGS 自由度，只在普通三线性/显式 XS 局优的内部试探步跳过仅氢叶的坐标变换，返回前恢复包括氢在内的完整构象。三个条件同时生效：单配体无柔性侧链、该力场下该叶没有对受体图或原子对的贡献、没有模板目标。AD4 继续更新氢。这个实现与改变抽样的两个实验独立。[实现](kinematics.patch) · [预先固定的计时](KINEMATICS.md)

原生轨迹对照中，122 个条件 × 两臂 × 各三次共 732 次运行，完整 PDBQT（包含氢）和搜索/精修求值数逐项相同。每格取三次中位数再合计，94 格 Gal3 的 wall 从 10.687 降至 10.622 秒，CPU 从 10.626 降至 10.539 秒；28 格结构的 wall 从 4.842 增至 4.980 秒，CPU 从 4.826 增至 4.955 秒。两组均未表现出足以推广的加速。44 组固定起点的验证包含严格 Armijo、显式/网格局优与对极性氢朝向敏感的 AD4，完整结果两版一致。该版本每个扭转树节点在本机多占 8 字节，保留为实验补丁。[逐次实耗](kinematics/timing.csv) · [独立局优回归](hydrogen_kinematics.csv)

## 复现

从仓库根目录出发，Python 环境使用既有结构面板依赖（RDKit、SciPy、Gemmi、NumPy，独立结构核查另需 Biopython 和 Meeko）。来自 Git 的父版、上游和研究补丁都可独立构建。若 Boost 未安装到系统路径，给构建命令附加同一份 `--boost-prefix`。本轮编译参数与父版一致。下文的 `$D` 指主工作树中包含本目录的绝对路径。

```bash
D="$PWD/benchmarks/silent_torsion"
W=$(mktemp -d)
git worktree add --detach "$W/base" 3d9dfb9
python "$W/base/tools/build.py" --output-dir "$W/parent" > "$W/parent.log"
python benchmarks/pair_gradient/build_arms.py "$W/official" --arms upstream upstream-counts
UP="$W/official/upstream/build/vina"
CT="$W/official/upstream-counts/build/vina"
P="$W/parent/vina"
for mode in always adaptive; do
  python "$D/build.py" "$W/base" "$W/parent.log" "$mode" "$W/$mode"
  python "$D/build.py" "$W/base" "$W/parent.log" "$mode" "$W/$mode-obs" --observer
done
for study in primary reduced; do
  mkdir "$W/$study"
  tar -xf benchmarks/rigid_flip/inputs.tar.xz -C "$W/$study"
  opts=(--always "$W/always" --adaptive "$W/adaptive")
  if [ "$study" = primary ]; then opts+=(--upstream "$UP" --counts "$CT" --parent "$P"); fi
  python benchmarks/rigid_flip/holdout.py dock "$W/$study" --study "$study" "${opts[@]}"
  python benchmarks/rigid_flip/holdout.py analyze "$W/$study"
done
python "$D/pose_summary.py" "$W"
python "$D/replay.py" "$W/primary" "$W/replay" --always "$W/always-obs" --adaptive "$W/adaptive-obs"
tar -xf benchmarks/pair_gradient/galectin_external/raw.tar.xz -C "$W" prepared
CPU=$(python -c 'import os; print(min(os.sched_getaffinity(0)))')
opts=(--inputs "$W/prepared" --parent "$P" --always "$W/always" --adaptive "$W/adaptive" --cpus "$CPU")
python "$D/affinity.py" run "$W/affinity" "${opts[@]}" --upstream "$UP" --counts "$CT"
for condition in 214039-12000 751483-3000 751483-12000; do
  python "$D/affinity.py" run "$W/confirm_$condition" "${opts[@]}" --confirm --seed "${condition%-*}" --cap "${condition#*-}"
done
for study in "$W/affinity" "$W"/confirm_*; do
  python "$D/uncertainty.py" "$study"
  python "$D/affinity_geometry.py" "$study" "$W/prepared"
done
python "$D/fixed_pool.py" "$W/primary" "$W/affinity" "$W/prepared" "$P" "$W/always" "$W/adaptive" "$W/fixed_pool.csv"
git worktree add --detach "$W/deferred" 3d9dfb9
git -C "$W/deferred" apply "$D/kinematics.patch"
python "$W/deferred/tools/build.py" --output-dir "$W/deferred-bin"
python "$D/kinematics.py" "$W/timing" "$W/primary" "$W/affinity" "$W/prepared" "$P" "$W/deferred-bin/vina" --cpus "$CPU"
```

运行 `python benchmarks/silent_torsion/verify_affinity.py` 可直接从四批归档复核全部 1410 个测定对接输出与实际求值、CPU/wall/RSS，以及按精确有理数和原始 Kd 计算的各项配对分。结构独立复核可将两个输入归档和各自的 raw.tar.xz 解包到独立目录，再执行 pose_summary.py 或 [transfer/verify.py](transfer/verify.py)。[完整轨迹观察器与数据重放](path/README.md)另附命令。提交中的低预算完整上游数值对照、其他预算的当前 Pro 数值对照使用不同分母；没有把不同随机重复或受体合并成额外实验分子。
