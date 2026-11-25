# 三配位蛋白锌位点的姿势选择

当共晶、相近结构或同系列实验证据支持配体以单个已知原子直接配位 Zn 时，可用此工具从同一配体的对接候选中筛选姿势。工具读取刚性受体及候选 PDBQT 中已有的坐标和分数。输入应保留对接使用的坐标系，并用配体准备的原子映射确认配位原子的 PDBQT serial。官能团的成键、质子化和齿数由使用者根据化学证据判断。

```bash
python3 tools/select_zinc_pose.py \
  --receptor receptor.pdbqt --poses docked.pdbqt \
  --zinc-serial 101 --donor-serial 12 > zinc_pose_selection.json
```

命令中的序号是示意值。`--zinc-serial` 指受体 PDBQT 中 Zn 的原子序号，`--donor-serial` 指候选 PDBQT 各 MODEL 中同一配位原子的序号。提供一个配体、一个受体状态和一次 Vina 或 FastDock Vina Pro 运行的候选。生成候选时可增大输出的 `--energy_range`，以便检查更多已有姿势。

## 几何规则

AutoDock4Zn 用 N/O ≤2.50 Å、S ≤2.80 Å 分类配位壳。它使用 N、O、S 分别为 2.00、2.10、2.25 Å 的配位参考距离。原方法针对 `NA` 的方向项在四面体空位上放置距 Zn 2.00 Å 的 TZ 虚拟点。官方 `zinc_pseudo.py` 计算三个受体供体所张平面的法线，并选取 Zn 所在的一侧。 [AutoDock4Zn 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4144784/) · [AutoDock Vina 官方脚本](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/example/autodock_scripts/zinc_pseudo.py#L426-L447)

本工具按上述配位壳寻找恰好三个 N/O/S，并要求它们来自可识别的蛋白残基。每个 Zn–受体供体距离至少为对应元素的参考距离减 0.35 Å。Zn 与供体平面的距离须至少为 0.20 Å，已占方向和空位方向之间的六个夹角均须落在 109.47° ±25° 内。供体平面退化、保留水占位、四配位或受体羧酸有两个氧同时在壳内等情形会按这些条件判定适用性。工具直接使用原子的坐标。

选择器按大小写读取 AutoDock 类型。`N`、`NA`、`O`、`OA`、`S` 和 `SA` 分别按真实元素计入配位壳。`Na` 表示钠，`Zn` 表示锌。Vina 的评分参数将 `Se` 映射到 `S`，本工具的硫配位几何只适用于真实硫原子。 [Vina 原生类型与别名表](../src/lib/atom_constants.h)

使用者须确认受体位点只由三个蛋白 N/O/S 配位。受体 PDBQT 中的 `F`、`Cl`、`Br`、`I` 或 `At` 如距 Zn 不超过最大配位壳半径 2.80 Å，工具按可能有额外占位而弃选。其他非 N/O/S 配位物种须由位点证据排查。

默认的 `--mode hard` 要求每个姿势的声明供体落在 `[参考距离 − 0.35 Å, 配位壳上限]` 内，其 Zn–供体方向与空位方向的夹角须 ≤30°。有其他配体 N/O/S 同时进入配位壳时，该姿势未通过本工具的单齿条件。虚拟点放在对应元素参考距离与空位方向的交点上，虚拟点距离作为诊断输出。距离下限、平面偏移和角度容差是本工具预先固定的筛选阈值。

## 磺酰胺氮的连续重排

已知配体以磺酰胺氮直接单齿配位时，可在相同候选池上指定 `--mode continuous`。该选项由使用者根据共晶、相近结构或同系列实验确认适用性。工具通过 PDBQT 检查指定原子为 `N` 或 `NA`，对官能团的确认依赖输入证据和配体准备时的原子映射。

```bash
python3 tools/select_zinc_pose.py \
  --receptor receptor.pdbqt --poses docked.pdbqt \
  --zinc-serial 101 --donor-serial 12 --mode continuous > zinc_pose_continuous.json
```

连续模式沿用受体适用性和单齿条件。指定氮距锌至少 1.65 Å，并且没有第二个配体 N/O/S 进入配位壳，姿势即可参与排序。距 Zn 2.00 Å 的空位虚拟点给出偏差 Δ，也就是指定氮与虚拟点的距离。奖励为 `3 exp(−Δ² / (2 × 1.00²))`，其中 Δ 和宽度均以 Å 计。排序分为原 Vina 分数减去奖励，越低越优先。相同分数按原 Vina 名次打破并列。超过原距离上限或角度上限的姿势仍可参与，并按连续偏差获得奖励。使用 `--mode radial` 可做径向消融，它用 Zn–N 距离与 2.00 Å 的绝对差作为 Δ，沿用相同资格条件和排序公式。

## 读取结果

- `status` 为 `selected` 或 `abstain`。有合格姿势时，`selected_model` 给出原输入 MODEL 编号。硬规则比较合格姿势的 native 分数，按原 native 池内名次处理输出精度下的并列。无候选、位点不适用或候选均未通过时，`selected_model` 为 JSON `null`，并给出原因。
- `site` 记录受体配位原子序号和距离，以及可估计的平面偏移、六个四面体夹角和空位单位向量。`poses` 逐姿势记录 native 分数、native 名次、如有的 Pro 分数、Zn–供体距离、空位夹角、距虚拟点的偏差和未通过的原因。
- 连续和径向模式在顶层记录 `mode`，并逐姿势追加 `scoring_eligible`、`scoring_reason`、`scoring_delta_angstrom`、`geometry_reward` 和 `ranking_score`。无资格的姿势仍保留在输入顺序中，其 Δ、奖励与排序分为 JSON `null`。连续得分按浮点计算结果输出。原有 `passes_geometry` 和姿势内 `reason` 仍表示硬规则诊断，连续排序使用新增字段。
- Vina 与 Pro 均从各 MODEL 开头的引擎诊断读取分数。Pro 读取开头的 `VINA NATIVE` 作为选择依据。位置更靠后的 `VINA RESULT` 或 Pro 分数注释可能来自之前的对接，计入 `ignored_trailing_score_headers`。格式错误以非零退出码报告，几何 abstain 仍正常输出 JSON。

几何判定与连续重排均用于选择已生成的姿势。评价应分别报告位点适用率、候选池通过率，以及选中姿势相对实验结构的几何偏差。亲和力、跨靶选择性和其他锌配位模式需要各自的对照实验。

```bash
python3 benchmarks/check_zinc_pose_selection.py
```
