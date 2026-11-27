# 主轴半转提议的结构留出验证

本目录对应首轮未参与机制设计的 ESR1、PTP1B 与 PARP1 实验结构。主搜索实现、分支准则与观察器见 [shape_flip](../shape_flip)。冻结和复验规定见 [HOLDOUT.md](HOLDOUT.md) 和 [CONFIRM.md](CONFIRM.md)。共用配体与受体准备耗时为每个结构 3.02–3.81 秒，具体看 inputs.json 中 seconds（含缓存原始文件读取）。每次运行都重新读输入、建网格。

| 研究 | 独立查询 / 靶点 | 每臂配对条件 | 原生 | Pro | 无条件 | 自适应 |
|---|---:|---:|---:|---:|---:|---:|
| 冻结面板 | 5 / 3 | 28 | 1 | 1 | 4 | 2 |
| 两个新种子，同一批结构，12000 | 5 / 3 | 14 | 1 | 1 | 3 | 2 |
| 原始种子、统一降额 15%，与首行 parent 比 | 5 / 3 | 28 | — | 1 | 4 | 2 |

表内为 Top-1 重原子 RMSD ≤2 Å 的配对条件数。三个阶段每臂的 Top-5 和返回池 oracle 都与各自 Top-1 相等，条件排序成功率均为 100%，分母极小。所有 RMSD 均作分子图对称校正，在蛋白受体坐标系计算。种子、预算和同一靶点内的交叉测量不能视作独立分子。输入均为 CCD 描述的中性微态，生理 pH 下的微态迁移未验证。PTP1B 两个受体附近有未纳入的晶体对称伙伴，距配体最近 6.20 和 5.40 Å；1T49 的蛋白准备保留远处 Mg。

* **ESR1。** OHT 从 3ERT 转移至 1SJ0，主面板 Pro 的 1/4 增至无条件 3/4、自适应 2/4。低、高预算的一部分新增命中来自同一个种子 202711。两个新增种子下两变体均 2/2，Pro 为 1/2。原始受体的口袋 Cα 对齐 RMSD 1.28 Å，转移后的晶体参考与 1SJ0 最近 1.17 Å、有两处重叠，局部优化使其偏离 1.17 Å。对 1SJ0 自身的 E4D，各臂均无命中。由输入构象重建共晶扭转、选择最优等价原子映射后，其局部优化结果为 2.46 Å；从沉积几何出发则为 0.52 Å。这个控制暴露了单一输入刚性片段几何对部分配体的限制。
* **PARP1。** Talazoparib 主面板种子 202727、高预算时，首姿从 Pro 的 7.25 Å 变为无条件 1.94 Å，自适应 4.69 Å。降额后无条件变体保住这一结果，求值数从 12154 降至 10443，三次交错重复的中位 wall 为 0.2228 s（parent 0.2139 s）。换两个种子均无命中。冻结起点的首次分岔在第 0 步，无条件半转端局优后 15.09 Å，原生分支 12.09 Å，均被接受。此例的终点收益涉及后续轨迹。
* **PTP1B。** 四种配体/受体组合在主面板四臂均无命中。新种子 91373 时，无条件提议把 1T48 的 BB3 送入 1T49 的参照位姿附近，1.85 Å；另一个新种子未复现。更大的 892 向 1T48 转移时，主面板 202711、高预算自适应使首姿从 8.96 恶化到 13.07 Å。首次同态重放的两个提议均被拒绝，局优分别耗费 107 和 203 次求值。最佳扭转投影起点在同一受体局优可稳定至 0.55 Å。

对主面板事后以 |ΔTop-1 RMSD|≥0.5 Å 描述端点变化，无条件提议改善/退化/近似不变为 10/12/6，自适应为 3/4/21。无条件与自适应分别在 27 和 19 个条件触发，累计执行 465 和 71 次半转。自适应未触发的九个条件逐字节复现 parent 输出。冻结面板的事后跨臂 oracle 最高 4/28；自适应保住全部原命中，并获得无条件相对 parent 新增三个条件中的一个。触发次数受先前决策影响，不能解释为直接测得的门控灵敏度。

ESR1 202711 高预算的无条件与门控两臂最终均达标。无条件首次切换（零起始步 2），原生/半转局优后分别 7.40/6.97 Å，均接受；门控首次切换（步 8）对应 7.44 Å 接受与 5.77 Å 拒绝。离线观察器重放共享随机数和实际首次操作前状态，每例核对完整最终输出字节；五个选定条件×两臂中一个门控未触发。已复算的首个分岔均没有一举进入 2 Å 内；其后的完整因果链未逐步冻结。仅从最终构象看，OHT 高预算 parent/自适应质心偏差 1.53/0.46 Å，刚体对齐后内部形状误差 2.02/0.94 Å，其余取向项 6.62/1.13 Å，三项平方和为重原子 RMSD 平方。取向项包含取向与内部形变的耦合。

三次交错重复比较全部 28 个原配置和降额变体，420 次复算输出全部逐字节一致。以下用每个条件三次中位 wall / CPU 后再求和。

| 臂 | 总梯度求值 | wall / CPU (s) |
|---|---:|---:|
| Pro | 218327 | 4.716 / 4.691 |
| 无条件 | 221491 | 4.837 / 4.823 |
| 自适应 | 219697 | 5.054 / 5.036 |
| 无条件，降额 | 189278 | 5.023 / 4.986 |
| 自适应，降额 | 188410 | 4.748 / 4.733 |

OHT 种子 202711 的高预算 Pro 在 0.1961 s 中使用 12167 次梯度求值，首姿 7.09 Å；自适应降额在 0.1714 s 中使用 10353 次，首姿 1.54 Å。此收益可同时满足本例的工作量、CPU 与 wall 控制。全体降额条件中，两个成本时间都不超过 parent 的有无条件 11/28、自适应 14/28，单个轨迹约为 0.1–0.3 s，不能外推为全批提速。原生 Vina 的 CPU 时间和 RSS 见 primary/cells.csv，28 次求值通过独立计数器实现且逐例与原生坐标字节一致。

## 复现

在仓库根目录安装 engines.json 记录的 Python 依赖后，执行独立评估。该命令对归档的 224 个主实验、降额与新种子结果，重新核算逐模型 RMSD、汇总、搜索计数、原始文件哈希。预期报告见 verification.json。

```sh
python benchmarks/rigid_flip/verify.py
```

完整重跑时，采用主实验所述构建方式得到三个 Pro 引擎。原生和计数器臂由已有脚本构建。下列示例假定 `PARENT`、`ALWAYS`、`ADAPTIVE` 分别指向三个 Pro 可执行文件。所有路径以新的工作目录为目标，防止已有缓存影响重算。

```sh
d=benchmarks/rigid_flip
w=$(mktemp -d)
python benchmarks/pair_gradient/build_arms.py "$w/up" --arms upstream upstream-counts
for study in primary reduced newseeds; do
  mkdir "$w/$study"
  tar -xJf "$d/inputs.tar.xz" -C "$w/$study"
done
for study in primary newseeds; do
  python "$d/holdout.py" dock "$w/$study" --study "$study" \
    --upstream "$w/up/upstream/build/vina" --counts "$w/up/upstream-counts/build/vina" \
    --parent "$PARENT" --always "$ALWAYS" --adaptive "$ADAPTIVE"
  python "$d/holdout.py" analyze "$w/$study"
done
python "$d/holdout.py" dock "$w/reduced" --study reduced --always "$ALWAYS" --adaptive "$ADAPTIVE"
python "$d/holdout.py" analyze "$w/reduced"
python "$d/geometry.py" "$w/primary" --crystal "$PARENT"
python "$d/cost_replay.py" "$w/primary" "$w/reduced" --parent "$PARENT" --always "$ALWAYS" --adaptive "$ADAPTIVE"
```

从原始结构重新准备可运行 `OPENMM_CPU_THREADS=1 python benchmarks/rigid_flip/holdout.py prepare NEW_DIR`。观察器通过 [observer.patch](../shape_flip/observer.patch) 构建后，用 `first_branch.py PRIMARY --always OBS_ALWAYS --adaptive OBS_ADAPTIVE` 重算轨迹首次分岔。单 CPU 计时在调用程序前以 taskset 固定至同一可用核心。
