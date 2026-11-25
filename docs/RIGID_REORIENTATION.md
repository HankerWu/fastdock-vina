# 保持配体质心的刚体重定向

用 `--rigid_reorientation mixed` 可在常规的单配体、刚性显式受体 Vina 搜索中改变旋转提议。默认 `native` 沿用原提议和随机流。逐配体批处理可以使用该选项。输入应提供受体和对接盒，并做全局搜索。两种提议均可和 Late 候选选择配合。

Monte Carlo 先沿用已有的平移、旋转和扭转自由度抽样。抽中旋转时，`mixed` 以 1/4 概率从均匀 SO(3) 抽取一个新的绝对方向，剩余概率使用原有局部旋转。均匀旋转由四个独立标准正态分量归一化成的单位四元数生成。新方向围绕当前已接受构象的化学重原子质心，保持扭转、重原子质心和内部几何。氢和虚拟点不参与质心，真实大环闭合碳计入。两个模式沿用既有网格优化、接受准则、轨迹数量、候选上限和显式末端精修规则。全局旋转有时同样会进入错误势阱，性能须用事先确定的新结构和相同预算实测。

```bash
build/rigid-reorientation/vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --config box.txt --scoring vina --cpu 1 --seed 42 --exhaustiveness 8 \
  --num_modes 9 --energy_range 100 --pose_selection late \
  --rigid_reorientation mixed --search_counters --out mixed.pdbqt
```

对照用同一二进制和参数，将模式改为 `native`，输出到另一文件。两臂都用 `--search_counters`。每条 `Search trajectory:` 记录实际 Monte Carlo 步数、每步首次网格最小化的调用及梯度计算数、接受后第二次网格最小化的调用及梯度计算数、以及全局旋转提议数（包含后来被拒绝的提议）。`Search counters:` 记录显式末端精修输入数、含盒重试的调用数、梯度计算数和 Late 筛选前显式分数的计算次数。梯度数包括每次 BFGS 的初始目标与梯度计算和线搜索的每个试值。后一条计数不含最终仅打分、输出去重或 RMSD。

`--max_evals` 的原始截止策略、局部步数和精修名额均保持。截止在完整的一轮优化后复查，实际梯度次数可因模式和收敛行为不同而变化。报告实计梯度数和墙钟，便于解释等参数比较。

[数值检查](../benchmarks/check_rigid_reorientation_numeric.cpp) 使用含可旋转键的 PDBQT，确认覆盖全局分支时重原子质心和内部距离保持，并故意给模型留下被拒构象以检查提议的坐标来源。

## 保持质心的局部转动

同一入口的 `--rigid_reorientation centered` 沿用原有局部角度提议。旋转中心由 PDBQT 根原子移至已接受姿势的化学重原子质心。只在抽中局部旋转时补偿根部平移，平移和扭转提议沿用原路径。若根原子位置为 `r`、重原子质心为 `c`、旋转矩阵为 `R`，按根原子旋转的质心位移为 `(R − I)(c − r)` Å。更新根部平移为其相反量，能保持质心在同一受体坐标内。该方法无全局方向重抽，`Search counters:` 会标出 `rigid_reorientation=centered`。

对照时用同一构建的 `--rigid_reorientation native`，匹配受体、配体、搜索盒、种子、梯度和显式精修名额。质心约束只作用于旋转提议，后续局优仍可移动质心。数值检查同时比较了两种模式在同一输入状态上的随机流、局部旋转四元数、扭转和平移。
