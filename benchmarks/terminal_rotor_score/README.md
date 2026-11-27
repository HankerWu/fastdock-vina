# 为什么同一构象会因为羟基氢能否转动而得到不同的 Vina 分数

普通 Vina 的显式氢决定供受体分型，受体–配体势只作用于重原子。对 ESR1–1SJ0 首姿的两个末端酚羟基，仅依次去掉控制氢转动的 PDBQT 分支标记。所有原子（包括氢）的坐标和类型保持相同。原实现给出 −4.947、−5.052、−5.161 kcal/mol，分子间项三次均为 −6.972，分子内项均为 −0.134。产生差异的第一步是转子描述符。从一个端点看需要另一个端点至少有两个非氢邻居，却没有对自身施加同一条件，然后两向求和除以二，每个末端 C–O 轴贡献 0.5。既有氢提议试验只改变 Monte Carlo 步骤，这里修复的是最终跨配体评分。[输入和反事实协议](PROTOCOL.md) · [官方原文，式 9 与表 1](https://doi.org/10.1002/jcc.21334) · [历史氢专属提议](../silent_torsion/README.md)

## 模型和反例

原生变换为 $g(U)=U/(1+wN)$，默认 Vina 的 $w=0.05846$。$U$ 为按本次运行的共有参考得到的分子间与内能总和减去参考内能。本文修订计数定义，保留原经验系数。对一条 PDBQT 声明的转轴，切开树得到两侧原子集合，各侧都包括自己的轴端。两个轴端均须为化学重原子，且每侧至少含有一个轴端以外的参考位点。氢与水代理 W 不提供这种参考。已断开的环键用 G 代理另一端真实碳的位置；它可以证明某个环轴两侧都有参考位点，但 G 自身不能成为化学轴端。解析配体时一次计算，多个配体的计数相加。Vinardo 使用相同的 XS 变换。AD4 按其自身系数直接乘输入 TORSDOF。

第一版要求距离推断的键图在两端均有重邻居。它通过酚羟基试验，却在因子 Xa 系列中找到反例。CHEMBL369480 的输入有一条 2.094 Å 的单键，坐标推断漏键。配体的声明运动树有七个双侧转轴，原算法计 5，第一版降至 4。最终实现直接从解析后的运动树计 7。对既有采用随机坐标起始重新制备的同一分子，旧版与最终版均计 7。对原始输入逐轴实际旋转也确认七次均改变了重原子内部距离。这个修复只解决转子定标，畸变输入的键长、局部原子分型和分子内相互作用仍需独立检验。[拒绝的实现与逐例结果](distance_gate) · [重新冻结的迁移协议](TRANSFER.md) · [两种准备态](results/distorted_input.csv) · [历史几何审计](../frozen_core/strain_transfer_audit/README.md)

修订后三种 ESR1 表示均给出 −5.161，Vinardo 均给出 −3.220。零场 AD4 对照的三个 0.2983 × TORSDOF 值在新旧二进制间逐项相同；真实 AD4 氢仍可参与势能。环十二烷原有的 11 个断环树转轴计数保留。独立脚本在 165 条输入上从 PDBQT 文本恢复转轴，并用 C++ 实际逐轴旋转核查。共线构型和仅移动断环代理的轴说明，双侧树计数衡量输入表示中的转轴，不能解释为独立力学自由度的维数或溶液扭转熵。原论文未对这个重新定义的描述符重新拟合 $w$。[固定坐标能量](results/recoding.csv) · [独立坐标核验](independent/README.md)

## 真实测定与完整搜索

三组分别使用各自同一实验的精确 Kd 或 Ki。每个配体先取同一受体内种子中位数，再对受体均值排序。以浓度相差至少三倍记方向命中，预测平局记半分；此阈值只表征浓度比。每一组既有父版首姿均在原始制备配体的拓扑上单独用新旧 CLI 复算，并独立用档案中三位小数的 INTER 和输入图计数重建全精度姿势上映射的近似得分。两种途径所得下表配对分逐项相同。输入、测定来源、筛选及测量条件见各原始面板。[Gal3 荧光偏振实验及入组](../pair_gradient/galectin_external/README.md) · [FXa Ki 来源与入组](../affinity_partition/README.md) · [S6K1 Kd、受体与入组](../affinity_reference_ensemble/README.md)

| 测定 | 分子 / 受体 / 首姿数 | 旧→新，三倍浓度差 | 改对 / 改错 | Spearman 旧→新 | Kendall τ-b 旧→新 |
|---|---:|---:|---:|---:|---:|
| Galectin-3 Kd | 47 / 2 / 94 | 116→113 / 403 | 1 / 4 | −0.359→−0.363 | −0.241→−0.239 |
| 因子 Xa Ki | 28 / 2 / 112 | 213→207 / 295 | 2 / 8 | 0.505→0.474 | 0.367→0.346 |
| S6K1 Kd，修订后冻结 | 14 / 1 / 28 | 31→32 / 58 | 1 / 0 | 0.068→0.103 | 0.044→0.066 |

Gal3 的 47 个糖配体都含半计数，每分子扣除 2.5 至 5 个旧计数，固定首姿变化最大 0.636 kcal/mol。1KJL 的三倍配对 146→147，3T1L 为 119→119，按分子平均两个受体后净减三对。FXa 的 27 个合规输入分数保持一致，全部变化发生在原始 CHEMBL369480 涉及的分子对。删除该分子，两个算法完全相同。此组的旧版 213 是按本协议先合成每个分子的点预测所得；既有文档中逐对取四个 receptor–seed 分差中位数得到的 223 属于另一个统计量。S6K1 的 14 种配体中 6 种计数改变。新校正将 CEP-701（69 nM）相对 PKC-412（1300 nM）的预测强弱差由 −0.024 变成 +0.608。三组逐分子删除后的净命中变化范围分别为 −4 至 −1、−7 至 0、0 至 +1。训练历史、化学系列、受体重复及共享分子使分子对不能作为独立样本。[234 行固定首姿](results/poses.csv) · [分受体、全部非平局对、数值精度、逐分子删除](results/summary.csv) · [逐对结果](results/pairs.csv.xz)

对 S6K1 又实际用新旧二进制完成 112 次完整对接，单 CPU，原始晚期选择，两颗种子，使用 2 或 4 条独立轨迹，每条最多 25,000 次名义求值。两档各有 28 个配对条件，新旧每个条件的完整输出坐标和搜索计数完全相同。2 条轨迹的 58 对为 31→32，4 条轨迹为 31→30。相同输入的已归档原版 Vina 在 2 条轨迹判对 33 对。前一档新旧各进行 1,418,605 次梯度求值，后一档各 2,820,650 次；新旧 CPU 秒从 11.29 到 11.39、17.08 到 17.24，单进程峰值均约 34 MiB。每个分子独立归一化会改变跨配体次序，提高搜索预算既可改对也可改错。此处的证据支持在默认路径修复输入表示依赖，未证明跨化学系列的亲和力预测改善。[逐查询与每次成本](results/s6k_search.csv) · [汇总](results/s6k_search_summary.csv) · [原始输出与日志](results/s6k_search.tar.xz)

结构检验另取 ESR1、PTP1B 和 PARP1 的五个共晶体系、两颗种子和 3,000/12,000 两个上限。原版 Vina、父版和修订版在这 20 格的 Top-1、Top-5 与已输出候选 Oracle 均为 2/20。使用 100 kcal/mol 输出窗，两个 Pro 版本的 80 个姿势坐标及 155,874 次总求值一致。五种输入原本已制备，复用成本归属上游基准；各次 CLI 包括读入、制图及搜索。两版在本机累计 wall 3.73/3.78 秒、CPU 3.71/3.77 秒、峰值约 25 MiB；原版累计 wall 26.50 秒、CPU 26.36 秒、峰值约 588 MiB。上述各臂各跑一次，时差不足以证明加速。默认 3 kcal/mol 输出窗使 1T49 种子 202727 高预算的第四名（RMSD 10.208 Å）越出窗口，输出总数 57→56，首姿保持不变。对同配体，这个正分母保持 native 排序，以绝对分差过滤的模式数仍可改变。额外的启发式重排器各有自己的分值定义。[逐结构结果](results/redocking.csv) · [完整 PDBQT 与运行日志](results/redocking.tar.xz)

## 重算

从仓库根目录使用 Python、NumPy、SciPy、C++14 和 Boost。所有实验读取本仓库先前冻结的结构与原始测定清单，无需重新下载。标准系统 Boost 路径的复现实例在下方；本轮通过同一 Boost 自定义前缀和等效链接选项构建。两份描述符探针各自链接同一版本的对象文件。上游二进制对应官方源码 3c65c0b，路径放入 VINA_NATIVE。需要独立复核时再运行 independent/README.md 的测试。

```bash
R=$PWD
D=$R/benchmarks/terminal_rotor_score
W=$(mktemp -d)
git worktree add --detach "$W/parent" 2e610c4
git worktree add --detach "$W/fixed" c0ce55f
for arm in parent fixed; do
  B="$W/$arm-build"
  python "$W/$arm/tools/build.py" --output-dir "$B"
  objects=()
  for obj in "$B"/obj/*.o; do
    case "$obj" in */main.o|*/split.o) continue ;; esac
    objects+=("$obj")
  done
  c++ -std=c++14 -pthread -O3 -DNDEBUG -I"$W/$arm/src/lib" \
    "$R/benchmarks/frozen_core/affinity/native_rotors.cpp" "${objects[@]}" \
    -lboost_thread -lboost_filesystem -lboost_serialization -lboost_program_options -o "$W/$arm-probe"
done
args=("$W/parent-build/vina" "$W/fixed-build/vina" "$W/parent-probe" "$W/fixed-probe")
python "$D/assay.py" "${args[@]}" "$W/results"
python "$D/check.py" "${args[@]}" "$VINA_NATIVE" "$W/results"
python "$D/s6k_search.py" "${args[0]}" "${args[1]}" "$W/results"
```
