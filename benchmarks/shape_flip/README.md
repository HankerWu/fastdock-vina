# 在重复盆地中改变配体取向，能提高对接命中率吗？

在新的人源 ESR1 受体中，4-羟基他莫昔芬（OHT）一次交叉对接的首姿从 7.09 Å 改善到 1.54 Å，实际求值和耗时均下降。该方向在新种子中也出现改善。其他体系的收益随输入和种子改变；47 个 Galectin-3 配体的实测 Kd 排序退化。本分支保留可编译的研究版本，默认程序继续采用原搜索。

## 提议、触发条件和预算

沿用 Vina 原有的 MC 实体抽样、局部 BFGS、接受、二次局优和候选上限。对单配体、刚性受体的原生搜索，旋转被选中时，可把原本的小转动替换成绕化学重原子质心的一次 180° 转动。所有输入来自当前已接受构象。三个轴来自化学重原子等权坐标协方差矩阵。原有球内随机向量与哪个轴的绝对点积最大，就选哪个轴。三轴被选概率相等。简并或不可分辨的谱退回原提议。[研究源码](https://github.com/HankerWu/fastdock-vina/commit/b8f9aad6bd5d757c83d6bd6a4e0ad2ed4fa7bf95) · [预先固定的协议](PROTOCOL.md)

对单位特征向量 u，转动矩阵为 R=2uuᵀ−I，满足 R²=I、det(R)=1 以及 RΣRᵀ=Σ。这个操作保留质心、键长、手性和二阶形状，却能把不同化学基团放到口袋的另一侧。更高阶的形状、碰撞和供受体匹配仍会改变。标准主轴的符号简并已被用于分子形状筛选，本轮检验了将它用于可变构象下的实际搜索提议。[VAMS 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4140985/) · [Vina 搜索与评分](https://doi.org/10.1002/jcc.21334)

无条件版本替换每次符合条件的旋转。门控版本复用原生候选池已有的 RMSD 比较。已接受且经过第二次局优的姿势与保存池重复时才获得一次翻转机会。下一次完成的翻转或下一次接受的构象会清除旧机会；若这次接受又产生重复保存，可建立新机会。拒绝其他扰动时保留当前状态。每条轨迹只维护一个布尔值。比较沿用原有原子对应和 RMSD 阈值，与离线的化学对称 RMSD 分开。每次启用承担一次正向运动学、化学重原子坐标读出、O(N) 二阶矩累积和一个 3×3 对称特征问题。所有坐标均在局部优化前提出，BFGS 的目标、梯度和历史矩阵更新规则相同，每次 hunt 重新初始化 BFGS。单一构象上的翻转可逆；带历史、局优和筛选的 MC 是寻找极小值的启发式搜索，不能从轨迹频率推断平衡分布或结合自由能。

三个配体共 3,600 次显式配对抽样，323 次翻转，其余提议坐标对照逐点相同，全部提议的随机状态相同。重复翻转的最大重原子坐标误差为 3.98×10⁻¹² Å 以内。测试包含陈旧模型坐标、点状及线性输入。研究源码中的新运算没有可微能量项。[数值测试](check.cpp) · [逐配体误差](development/numeric.csv)

## 发生改善和退化的阶段

新面板的主试验、全体条件统一降额、两枚未用种子、原始输出、单核耗时和几何审计见[独立结构检验](../rigid_flip/README.md)。对 28 个主要条件，Pro、无条件和门控分别命中 1、4、2 个首姿。正确候选存在时都已排第一。以首姿变化至少 0.5 Å 描述幅度，门控将无条件的 12 次退化减至 4 次，同时把 10 次改善减至 3 次；逐一匹配，无条件的 10 次改善中仅有 2 次在同格保留下来。主面板的 9 个门控未启用条件，输出逐字节与 Pro 一致。

旧 5T9W 和 7MBO 作为开发反例。四臂各 12 个条件均无 2 Å 输出。在高档首次分岔上，离线复制整个接受状态、残留模型坐标、生成器并复用相同随机向量，再用同一局优起点生成规则和 hunt 目标分别最小化。复制后的生成器状态强制检查相等。观察器的最终输出逐一与未插桩版本比对。在发生翻转的 11 项追踪中，输出的 55 个源、提议和局优姿势最低 RMSD 仍为 9.83 Å。5T9W 种子 740317 的无条件首分岔中，原生和翻转 hunt 分别有 36.878 和 39.177 Å 的逐原子逐轴盒外距离之和；乘以本引擎 10⁶ 的盒壁斜率，已解释两个数千万量级的 hunt 目标。该状态的主导误差是盒外落位，重复保存不能诊断取向是否放反。[逐状态几何](development/stages.csv) · [能量、实耗与接受决定](development/decisions.csv)

另一项离线干预取上一轮持续固定环的三个高预算首姿，冻结每个原子的内部几何。以晶体答案分别纠正质心、绕质心的最佳取向以及二者，对照三个实际可提出的主轴半转。5T9W 即使同时纠正质心与取向，三个首姿仍为 2.79、3.34、3.24 Å。7MBO 仅纠正质心的范围为 2.36–8.66 Å，仅纠正取向为 2.69–6.96 Å；联合纠正为 1.11–1.45 Å，三者均无小于 2 Å 的受体重原子接触。联合放置的三个姿势却都比原错误首姿的固定参考分更高。种子 740317 从 −6.944 变为 −2.381 kcal/mol。这是同一内部几何、未重新松弛的分数，不能用来判断局优后的正确盆地稳定性。[全部 42 个干预](development/placement.csv)

## 亲和力

使用此前归档的 47 个中性输入微态、同一 Galectin-3 Kd 测定和两个受体，种子 391781、单轨迹、名义上限 3000、九个输出。两个受体首姿分取均值。403 个至少相差三倍的分子对，原生 Vina、Pro、无条件、门控分别为 220、201、176、171。等重原子数且等形式电荷的 35 对为 22、22、21、16。Pro、无条件、门控全体 Spearman 为 0.035、−0.061、−0.096；Kendall 为 0.020、−0.046、−0.078。依次去除任意一个配体，无条件的配对差始终是 −33 至 −12，门控是 −39 至 −12。[同测定数据来源](../pair_gradient/galectin_external/README.md) · [逐受体排序](gal3/ranks.csv) · [逐条件工作量](gal3/runs.csv)

在相同的 94 个 Pro 首姿上，三种 Pro 程序的十项打印评分完全一致。输出到小数三位的化合物 28 在 3T1L 上有一个原子越过原盒面约 0.0004 Å，三臂均拒绝原盒上的 score-only。仅对这个固定坐标诊断将三个方向盒宽同时增加 0.75 Å，完整搜索沿用原盒。开发集的全部 36 个父版候选也在同一分子内参考下复算一致。[固定首姿](gal3/fixed.csv) · [固定完整开发池](development/fixed_pools.csv)

## 复算

从包含研究提交的完整 Git 检出执行。首两条命令分别从归档原始输出核对亲和力、命中和结构面板，重新对接需要 Boost、C++14 以及独立验证目录 engines.json 记载的 Python 依赖。默认脚本另外构建两套离线观察器。上游计数版仅用于被动复算，实测成本来自无插桩的上游二进制。

```bash
D=$(pwd)/benchmarks/shape_flip
python "$D/audit.py"
python benchmarks/rigid_flip/verify.py
W=$(mktemp -d)
bash "$D/build.sh" "$W/core" --jobs 4
python benchmarks/pair_gradient/build_arms.py "$W/up" --arms upstream upstream-counts
mkdir "$W/inputs"
tar -xf benchmarks/macrocycle_proposal/external/inputs.tar.xz -C "$W/inputs"
U="$W/up/upstream/build/vina"
UC="$W/up/upstream-counts/build/vina"
engines=(--parent "$W/core/parent/vina" --always "$W/core/always/vina" --adaptive "$W/core/adaptive/vina")
python "$D/development.py" --inputs "$W/inputs" --output "$W/dev" "${engines[@]}" --upstream "$U" --upstream-counts "$UC"
python "$D/state_replay.py" --inputs "$W/inputs" --base "$W/dev" --output "$W/trace" --always "$W/core/observer-always/vina" --adaptive "$W/core/observer-adaptive/vina"
python "$D/placement.py" --inputs "$W/inputs" --parent "$W/core/parent/vina" --output "$W/placement.csv"
python "$D/fixed_pools.py" --inputs "$W/inputs" --base "$W/dev" --output "$W/fixed.csv" "${engines[@]}"
python "$D/assay.py" --work "$W/gal" "${engines[@]}" --upstream "$U" --upstream-counts "$UC"
```

本轮独立验证发现了等成本下可重复的局部收益，也确定了两种损失来源：重复保存无法指明需要改变哪个自由度；固定质心的翻转无法修复质心偏移。全体等成本、跨受体和同测定的结果尚不支持将两种策略纳入生产默认。
