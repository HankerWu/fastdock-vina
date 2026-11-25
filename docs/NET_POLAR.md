# 极性溶剂面积与氢键方向重排

`--net_polar native|direction|burial|combined|water|satisfied|satisfied_water` 适用于单配体、显式刚性受体的 Vina。常规对接保持原入口和默认评分。开启后，用 `--pose_selection late` 及相同搜索盒、预算和可选的 `--refine_candidates` 生成同一原生精修候选池。原生去重、截断和共同未结合能参照均先确定。最后按所选分数重排，同分沿用原生名次，`--energy_range` 仍基于原生分数过滤。与其它 Pro 重排、Armijo 和轨迹补存互斥。

## 分数

方向模式逐对处理原生 Vina XS 的配体–受体 N/O 氢键。它使用实际氢键权重、XS 距离势和原生快速表的两个采样端点，得到该氢键在原始总分中的距离贡献。对于有受体原子共价重原子方向的配对，按受体原子键与指向供体的方向之间的最小角给出乘数。≤80° 时为零，≥110° 时为一，中间对余弦线性插值。两种供受关系可用时取较大乘数。无法估计受体原子的共价方向时保留原生贡献。方向项从原生原子配对总和扣除原贡献并代入乘过方向系数的贡献，然后沿用原子级非线性能量压缩。原生配体内部氢键、金属相互作用和其它原生项保留。

埋藏模式对每个配体供体/受体 N/O 位点计算同一个配体构象在游离和复合物中的溶剂可及面积差。使用 Shrake–Rupley 球面点计数，256 个固定 Fibonacci 方向，1.4 Å 探针，重原子半径依 C、N、O、S、P、F、Cl、Br、I、Si、At、金属依次为 1.70、1.55、1.52、1.80、1.80、1.47、1.75、1.85、1.98、2.10、2.02、1.80 Å。氢原子和虚拟位点没有独立球面。单个采样点的面积为扩展球面面积除以 256。如果该点在游离配体中可及、且加入受体后被遮挡，则计入埋藏面积。固定方向使重复执行可复现，整体旋转可能改变离散采样结果。每埋藏 1 Å² 加罚 0.02 kcal/mol。

水桥臂使用 `--net_polar_waters FILE` 读取独立结构证据给出的水氧。输入为配准至受体 PDBQT 坐标系后的一行一组三个埃坐标，允许空行和 `#` 注释。水氧到配体和蛋白极性位点的两条边均在 2.4–3.5 Å，水处夹角在 75–135°，且水氧到每个重原子至少 2.2 Å。每条边距离质量为以 2.9 Å 为中心、宽度 0.35 Å 的高斯，纯受体位点沿用方向门控。供体端允许由水接受氢键。桥质量为两条边质量之积。每个配体极性位点最多取其最好的一个桥，以对应的质量比例抵消该位点的埋藏面积罚。水桥假定水氧已独立选定，可检验空间与极性匹配，不显式估计水占据率、取向或热力学。宜使用与外部供水证据分开的干受体输入，重合的内置水氧会由空间冲突拒绝。缺少水证据时无法选用水桥臂。

`satisfied` 以每个配体位点最好的直接蛋白伙伴质量 qᵢ 抵扣该位点面积罚。它复用极性几何评分中 4 Å 内互补 N/O 的距离高斯（中心 2.9 Å，宽度 0.45 Å）及上文的氢键受体角度门控。缺可用氢键几何时质量为零。若埋藏面积为 Aᵢ、归一化柔性因子为 F，增量等于方向替换增量加 `F × 0.02 × Σ Aᵢ(1−qᵢ)`。

`satisfied_water` 要求同样的独立水点文件。公式用 `max(qᵢ, wᵢ)` 代替 qᵢ，wᵢ 取上文最好的水桥质量。水增量抵扣面积为 `Σ Aᵢ max(0,wᵢ−qᵢ)`。同一位点的两条路径共享该位点埋藏代价上限。无水点时可用 `satisfied`。输入水点但未形成合格桥的姿势得到零水增量抵扣。

所有模式的基底都使用原生分数与同一个未结合能参照。方向增量替换原生氢键项。新增抵扣限于所加的面积罚。各增量沿用 Vina 的同一柔性缩放，权重是研发近似。

## 固定坐标诊断

```bash
vina --receptor receptor.pdbqt --ligand pose.pdbqt --score_only --net_polar combined --cpu 1 --verbosity 0
```

输入已有池的单个 PDBQT 构象时移除 `MODEL` 和 `ENDMDL`。此模式在输入坐标上计算，使用显式受体原生查表和零共同参照，不生成亲和力网格或执行优化。主输出每项均有十二位小数，计数为整数。

| 行前缀 | 字段顺序 |
| --- | --- |
| `Net polar terms:` | 原生及方向替换后的 N/O 氢键分量（原子级压缩前）、方向替换增量（原子级压缩后）、埋藏面积、获水桥补偿面积 |
| `Net polar deltas:` | 方向、埋藏、合用、合用带水的柔性归一化增量 |
| `Net polar native:` | 原生显式受体间能、原生配体内能、以零为共同参照的归一化两项和 |
| `Net polar normalization:` | 原生柔性乘数 |
| `Net polar water points:` | 输入水氧数，0 表示没有水证据 |
| `Net polar satisfaction:` | 直接抵扣面积、水在直接抵扣之后的新增面积、相应的两个归一化正抵扣量、`satisfied` 和 `satisfied_water` 的归一化分数增量。仅新模式输出 |

复排原生保存池时，保留文件既有原生分数，再加对应增量。以已保存坐标复算使用 PDBQT 精度，近似并列处可能和内存全精度排序有差异。无水证据时，旧 `deltas` 的末列等于 `combined`，新 `satisfaction` 的末列等于 `satisfied`。这两个末列在无水输入时只是诊断值。

全局重排输出保留标准 `REMARK VINA RESULT` 和能量分解。附加 `REMARK VINA NATIVE`（全精度原生分数和原生名次）及 `REMARK VINA NET POLAR RESULT`（重排分数和相对首位的 RMSD 上下界）、`TERMS`、`DELTAS`、`MODE`、`WATER POINTS`。两个新模式另外给出同字段顺序的 `REMARK VINA NET POLAR SATISFACTION`。将保存的 MODEL 拆出作为 Late 对接输入时，Late 输出会清理随输入携带的旧评分备注，重新写入当前评分备注。

在引擎根目录编译和验证。非系统 Boost 可使用构建器的 `--boost-prefix`。新测试的 Python 环境需要 [固定池工具](FIXED_PHYSICAL_POOL.md) 的依赖。

```bash
python3 tools/build.py --jobs 4 --output-dir build/satisfaction-final
python3 benchmarks/check_net_polar.py build/satisfaction-final/vina build/polar-burial/vina
python benchmarks/check_net_polar_satisfaction.py build/satisfaction-final/vina build/net-polar/vina
```
