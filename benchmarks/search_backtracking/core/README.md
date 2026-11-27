# 网格搜索中的 BFGS 回溯

原生搜索会在连续十次试步都不满足下降条件时，把最后实算坐标及其较高能量交给 BFGS，供下一步使用。此时坐标沿方向移动了 \(2^{-9}p\)，返回给曲率更新的步长却是 \(2^{-10}\)。如果整次局部优化末态比初态更差，原流程仍会回到初态。在一次局部优化内部，耗尽会同时影响搜索路径和曲率估计。该边界行为此前已在最终精修的严格 Armijo 数值测试中记录。本实验量化它在 Monte Carlo 网格 hunt、save 两阶段的实际发生率和姿势后果。[原数值测试](../../check_armijo_numeric.cpp) · [原生局部优化器](../../../src/lib/bfgs.h)

## 实际拒绝的步长

在八个既有靶点、24 个固定起点的 native 实验中，共 34 条线搜索耗尽，其中 hunt 12 条，save 22 条。第十次试步比该条线搜索起点更高的有 **34/34**，差值为 0.00000318 至 5.44845 kcal/mol，中位 0.03995。测量的是对应阶段的网格目标函数。全部 34 条初始方向导数为负，末次实算方向导数为正。最明显的是 PHIP 5S8I，63419 起点，save 阶段第 9 条线搜索。起点目标值为 −6.20063 kcal/mol，末次候选提高 5.44845 kcal/mol。

沿这些已记录的方向，再单独探测 \(2^{-10},\ldots,2^{-19}\) 后，34 条线搜索均能找到满足原充分下降条件的步长。分别多缩小 1、2、3、4、5、6、9 次才首次成功的条数为 7、10、4、4、1、6、2。在最小步长处，方向有限差分与起点解析方向导数之比中位为 1.000025，范围 0.7807–1.0387，34 个符号全部一致。这说明十次尝试覆盖的步长区间对这些局部状态偏粗。整条线的光滑性、具体原子碰撞或导数跳变的因果贡献仍待判断。探测仅用于诊断，每条另外计算十次能量和梯度，再用一次计算恢复原坐标缓存，共 374 次。各阶段的 PDBQT、原计数及最终输出与不探测时逐一相同。[逐条原值与试探点](exhausted_lines.csv)

## 机制和对照

设已拒绝步长为 \(a\)，\(g_0=\nabla f(x)^Tp\)，\(\Delta=f(x+ap)-f(x)\)。单位步长仍首先尝试。二次模型 \(q(t)=f(x)+tg_0+bt^2\) 给出 \(b=(\Delta-ag_0)/a^2\)，下一步和上一步的比为

\[
\rho=-\frac{ag_0}{2(\Delta-ag_0)}.
\]

把这个无量纲比值限制在 [0.1, 0.5]。若斜率非负、残差非正或非有限，则使用 0.5。能量和原方向导数已由现有计算给出，每次拒绝只增加常数次标量运算。两条实验路径都在十次失败后结束本次局部优化并恢复最后接受点。Halving 仍逐次二分，quadratic 使用上述比值。两者保留原充分下降条件、其余 BFGS 更新、搜索预算和末端精修。该插值是已有数值优化方法，本次研究检验它与 Vina 搜索能量面的适配。[Dennis 与 Schnabel 的技术报告，算法 4.2.1](https://mate.uprh.edu/~pnegron/notas6682/TR87-23.pdf) · [四个核心文件的修改](https://github.com/HankerWu/fastdock-vina/commit/b26101b9028b24e488718024c8ba01e3c680d813)

## 从相同起点进入相同局部优化

每个靶点包含投影共晶和两个错误的既有原生首姿。相同坐标在各臂进入恰好一次 hunt 和一次 save。局部步数上限固定，达到停止条件的实际迭代数可不同。表中 FG 是一次联合能量与梯度评估，两次局部优化的计数相加。达到 2 Å 的判断使用受体坐标系内按已登记等价原子映射修正的重原子 RMSD。[预设方案](PROTOCOL.md) · [输入与映射](../../geometry_comparison/pose/README.md) · [逐靶点与起点](rows.csv) · [全部坐标及日志](raw.tar.xz)

| 回溯 | 24 个起点总 FG | 拒绝试步 | 耗尽线搜索 | 共晶起点终态中位 RMSD（Å） | 错误起点终态中位 RMSD（Å） |
|---|---:|---:|---:|---:|---:|
| Native | 2077 | 1229 | 34 | 0.575 | 7.212 |
| Halving | 1039 | 711 | 31 | 0.384 | 7.207 |
| Quadratic | 1639 | 796 | 4 | 0.452 | 7.202 |
| 固定缩至四分之一（探索性） | 1806 | 1021 | 7 | 0.480 | 7.217 |

原生投影起点中位 RMSD 为 0.201 Å，错误起点为 7.307 Å。四种路径都保住 8/8 个投影起点，都未救回 16 个错误起点。相对于 halving，quadratic 容许优化器继续做更多后续迭代，总 FG 增加 57.7%。USP5 投影姿势从 halving 的 1.451 Å 改善到 1.317 Å，而 BCL2L1 从 0.388 Å 退化到 0.678 Å。两例都保持在成功阈值内。

事后登记的固定四分之一对照在原二分序列中每次跳过一个步长。它将耗尽从 31 次降至 7 次，说明快速缩短步长已经解释了大部分耗尽减少。二次插值额外少耗尽 3 次，较该对照少用 167 次 FG，节省 9.2%，这组实验的姿势命中数未随之提高。该对照仅检验局部机制。[逐起点](quarter_rows.csv) · [原坐标和日志](quarter_raw.tar.xz)

在 PHIP、BRDT、USP5 的六个完整对接中，含计数器的 native 臂与同源码、同编译设置的原生产物给出了逐字节一致的姿势文件、搜索和显式精修计数。三个口袋均用原盒中心和 24 Å 边长，具体种子、预算、输出文件哈希和调用数可复核。[对照记录](native_parity.json)

## 复现

从仓库根目录运行。系统安装的 Boost 可省略编译命令中的 `--boost-prefix`。目录须为新路径。

```bash
python benchmarks/search_backtracking/core/build.py /tmp/backtrack-main --boost-prefix /path/to/boost --jobs 2
python benchmarks/search_backtracking/core/build.py /tmp/backtrack-local --first-hunt-replay --boost-prefix /path/to/boost --jobs 2
python benchmarks/search_backtracking/core/replay.py /tmp/backtrack-local/build/vina /tmp/local-three
python benchmarks/search_backtracking/core/replay.py /tmp/backtrack-local/build/vina /tmp/local-quarter --only quarter
python benchmarks/search_backtracking/core/build.py /tmp/backtrack-profile --profile-exhaustion --boost-prefix /path/to/boost --jobs 2
python benchmarks/search_backtracking/core/replay.py /tmp/backtrack-profile/build/vina /tmp/local-profiles --only native
```

第一条命令编译实际端到端研究用的同一可执行文件，对照由 `native/vina`、`halving/vina` 和 `quadratic/vina` 三个薄包装选择。构建脚本从固定父提交导出核心源码，并检查数值最小值、方向导数有限差分、收缩上下界、非有限试步、十次失败后的坐标恢复和调用数。局部回放的 `--max_evals 1` 仅用于强制原生 Monte Carlo 完成第一对 hunt/save 后退出；实耗是表中数值。局部回放和诊断探测二进制不用于端到端耗时对照。

复算六个 native 文件哈希时，从局部回放解包的各 PDB 目录取 `receptor.pdbqt`、`ligand.pdbqt` 和 `metadata.json` 的中心，显式指定 24 Å 各边、`--cpu 1 --exhaustiveness 1 --num_modes 9 --energy_range 100 --search_counters`，并按 `native_parity.json` 各行设置种子和预算。

逐轨迹的每阶段均满足 `FG = 局部优化次数 + 线搜索次数 + 被拒试步数 − 耗尽次数`。末端显式受体精修在三个全局臂中继续调用原优化器。单次完整计算的差异还包括由搜索结果改变而带来的精修输入与调用数。
