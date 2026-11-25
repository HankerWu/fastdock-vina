# 固定姿势的物理对照评分

`tools/score_physical_pool.py` 使用能报告 `--net_polar` 分项的研发构建。给它一个已保存完整 MODEL 的原生 Vina 池、池对应的 PDBQT 受体，以及可选的独立结构证据。它逐姿势执行显式受体 `--score_only`，保留已有原生分及顺序作为共同基底。CLI 不搜索，不改变候选坐标。Python 环境需要 RDKit、Meeko 和 NumPy。

| 臂 | 原生分上的处理 |
|---|---|
| native | 保存的原生 Vina 姿势分 |
| direction | 引擎替换原生氢键项后的方向差 |
| burial | 引擎未补偿极性面积罚项 |
| combined | 方向替换和埋藏修正 |
| water | 在 combined 上加入独立水补偿。表中同时保存 water 减 combined 的独立补偿差值 |
| strain | 用 MMFF94s 应变代理替换对应姿势的原生显式配体内能，用同配体的气相零点替换原生未结合能参照 |
| sidechain | 原生受体状态与证据允许的替代侧链之间按显式分子间能差和状态代价比较，选较小分 |
| joint | 在允许的受体状态间按可用极性、水、应变及侧链代价共同选较小分。缺项和弃权另有明细 |

配体图从 Meeko 写入的 SMILES 和原子映射复原。立体及微态保持独立准备时的定义。生成八个固定种子的 ETKDGv3 气相候选并做完整 MMFF94s 最小化，最低收敛能量作为该配体的共同零点。每个对接姿势锁定重原子，只松弛氢。应变替换在相同柔性因子下使用 `0.1 × (固定姿势的 MMFF94s 能量 − 气相零点) − (REMARK INTRA − REMARK UNBOUND)`。原生两项读取同精度保存的 PDBQT 字段，并验证所有姿势共享 UNBOUND，原生第一名内能相对于该值恰为零。正负值及每个构象的原始能量都保留。气相零点是小规模构象近似，权重为预定未拟合的模型设定。

替代受体须保留每个原子顺序、身份、电荷和 AutoDock 类型，仅修改证据支持的侧链坐标。使用 `−0.592 ln(p_alt / p_native)` 作为观察频率先验，并加入因位移改变的跨残基重原子排斥差。排斥取 `2 max(0, r_i+r_j−0.4Å−d_ij)^2`，同残基原子对排除。只计算至少一原子位移的配对。侧链分将显式分子间分差、频率先验和排斥差之和用真实 Vina 柔性因子归一化。联合臂在替代受体上重算极性分项。水臂或侧链臂遇到缺证据时弃权，联合臂由其余可用分项回退。

独立水点文件每行是受体坐标系下三个以空白分隔的氧坐标，单位 Å，允许 `#` 注释。供体实验结构、链映射和对齐应在评分之前另存。侧链文件格式如下，受体路径相对该文件目录，JSON 的来源字段由评测冻结保存。

```json
{
  "baseline_population": 0.6,
  "alternates": [
    {"id": "Ile61-alt1", "receptor": "independent_alternate.pdbqt", "population": 0.4}
  ],
  "provenance": {"source": "independent experimental structure and alignment record"}
}
```

```bash
python tools/score_physical_pool.py \
  --binary build/net-polar/vina --receptor receptor.pdbqt --pool native_poses.pdbqt \
  --waters aligned_waters.xyz --sidechains independent_alternates.json --output fixed_scores.json
```

用 `build/satisfaction-final/vina` 和 `--satisfaction` 可另外输出 `satisfied`、`satisfied_water` 的姿势分及排序。没有独立水文件时，`satisfied_water` 分数为 `null`、排序为空，`abstentions` 给出原因。每个姿势的 `engine` 项含两个抵扣面积（`direct_credit_area`、`additional_water_credit_area`）、两个正抵扣量（`direct_credit`、`additional_water_credit`）和两个总增量（`satisfied_delta`、`satisfied_water_delta`）。这两臂使用原生受体，现有联合臂继续沿用既有分项。省略开关时沿用原版 JSON 格式。 [极性评分公式](NET_POLAR.md)

```bash
python tools/score_physical_pool.py \
  --binary build/satisfaction-final/vina --receptor receptor.pdbqt --pool native_poses.pdbqt \
  --satisfaction --output satisfaction_scores.json
```

用于 GNINA 固定姿势对照时，`--export-gnina-sdf` 将输入保存为每个 MODEL 对应一条 SDF 记录。导出校验与下文 GNINA 导入共享重原子元素和坐标核对规则。大环闭合的 `CG0` 等是真实碳，计入核对。`G0` 等闭环辅助点及 `W` 水探针在核对时跳过。每条记录有 `native_model` 字段。全部记录通过校验后才创建目标文件，既有目标受到排他创建保护。

```bash
python tools/score_physical_pool.py --export-gnina-sdf \
  --pool native_poses.pdbqt --output new_native_poses.sdf
```

`--gnina-scores-sdf` 读取 GNINA 1.3 的同坐标 `--score_only --cnn_scoring rescore` 输出。每个原生 MODEL 需恰有一条记录，且重原子元素与数字坐标的多重集合要完全相同。氢和非原子辅助点不参与比对。有限且在 [0, 1] 内的原始 `CNNscore` 保存为 `gnina`，按分数降序排列。合分 `cnn_consensus` 用原生分减 `0.592 × ln(c/(1−c))`，c 为限制到 `[1e−6, 1−1e−6]` 的 CNNscore，按分数升序排列。两臂同分沿用原生名次。经验系数 0.592 用于排序对照。可指定引擎和受体，将这两臂与物理评分写入同一份 JSON。

```bash
python tools/score_physical_pool.py \
  --pool native_poses.pdbqt --gnina-scores-sdf gnina_scores.sdf --output cnn_scores.json
python benchmarks/check_gnina_consensus.py native_poses.pdbqt gnina_scores.sdf
```

## 芳香环几何与近简并排序

```bash
python tools/score_physical_pool.py --pool native_poses.pdbqt \
  --receptor receptor.pdbqt --aromatic-ccd component.cif --output aromatic.json
```

输入使用已准备、同坐标系的非共价对接受体和固定 MODEL 池，以及该配体独立的 CCD 单组分 CIF。此模式额外需要 Gemmi 来读取 CCD 的原子和键表。无需启动引擎。配体坐标取 Meeko 的 SMILES 原子映射，并对照当前 MODEL 的 PDBQT 重原子元素和坐标。先以 CCD 原子、电荷、氢计数及键级建立独立 RDKit 图，再检验完整重原子图及芳香环成员。保留所有图等价的 CCD 环名，仅在其给出同一组几何环时允许等价映射。分析不依赖晶体配体坐标。重原子匹配用于识别芳香环，立体构型仍应通过相应输入与验证协议检查。

受体使用完整的 Phe、Tyr、His 及常见命名变体、Trp 标准芳香环，Trp 的两个最小环分别计数。按标准模板的环顺序核验每条相邻边，包括闭环边。化学上合理的芳香 C–C 和 C–N 键长落在宽限范围 1.2–1.6 Å 内，最小二乘环面 RMS 需 ≤ 0.2 Å。原子残缺、重复、元素不一致、边长不合理或坐标退化的环不纳入几何计数，排除原因随结果保存。PDBQT 没有受体键表，边长检查只能验证几何连通性。若准备过程将一个完整芳香环整体赋予错误残基标签，仍需用同坐标系的原始受体核对所报告的残基身份。环对以 PLIP 的几何阈值诊断：环中心距离严格在 0.5 和 5.5 Å 之间，两种环面投影偏移的较小值严格小于 2.0 Å，折叠后的环面夹角小于 30° 为平行（包括完全平行和错位），大于 60° 为 T 型。保存两个偏移和拟合平面 RMS，包括在 5.5 Å 内未通过其他条件的配对。配体每条环键按 RDKit 元素共价半径和检验，接受区间为该和减 0.30 Å 至加 0.20 Å，以覆盖芳香缩短及 S 等杂原子，平面 RMS 使用同一 0.2 Å 上限。这些识别值是近邻几何代理，未包含短程排斥、取代基电子效应、质子化耦合、局部溶剂或阳离子–π。融合环按不同环对单独计数，一个残基与配体可以贡献多对。 [PLIP](https://doi.org/10.1093/nar/gkv315) · [几何实现](https://github.com/pharmai/plip/blob/master/plip/structure/detection.py)

只有按原生输出的十进制分值，距全池最低原生分不超过 **0.5 kcal/mol** 的模型参与重新排序，按识别出的环对数降序，同数沿用原生序，其余模型接在后面保持原生序。这个窗口和环对数排序是预设实验假设。模型的坐标、原生能量、未结合能参考和其他物理臂均沿用原输入。任何模型的配体映射失败、配体缺少芳香环或受体缺少可辨认的标准环时，`orders.aromatic` 返回全池原生顺序，`abstentions.aromatic` 记录原因。局部缺失受体环的计数仅代表可辨识部分。

`poses[].aromatic_contacts` 是几何环对数，`poses[].aromatic_ring_pairs` 保存逐环对测量。`aromatic_geometry` 保存窗口成员、已诊断模型、CCD 映射、受体环和排除原因。批量评测可调用 `score_aromatic_pool(pool, receptor, ccd)`，合并现有引擎物理打分时使用 `score_pool(..., aromatic_ccd=ccd)`。加入 `--gnina-scores-sdf` 可在同一份 JSON 输出中保留 GNINA 顺序。

## Vinardo 固定坐标排序

```bash
python tools/score_physical_pool.py --pool poses.pdbqt --receptor receptor.pdbqt \
  --vinardo-binary /path/to/vina_1.2.7 --vinardo-box box.txt --output vinardo.json
```

用同配体池、同一刚性受体和搜索盒。入口逐 MODEL 调用 AutoDock Vina 的 `--score_only --scoring vinardo --unbound_energy 0 --cpu 1`。它按显式分子间能量与内部能量之和升序排序，保留打印精度，完全同分时沿用原生次序。各模型的 `vinardo`、`vinardo_terms` 保存和分、分子间及内部能量、进程墙钟。模型没有移动，所用分量等价于给同一配体的姿势指定共同零参照，仅作为固定池选姿分。跨分子亲和力需要另行验证。[Vinardo 方法](https://doi.org/10.1371/journal.pone.0155183)

`orders.vinardo` 保存完整顺序。任一模型无法取到分量时整池采用原生排序，`abstentions.vinardo` 及 `vinardo_scoring.failures` 记录原因和模型；成功分量仍可查。 `vinardo_scoring` 也记录可执行文件、输入、总墙钟及覆盖。该选项可独立运行，或与既有物理、芳香和 GNINA 臂写入同一份 JSON。项目现有默认调用无需变动。

## 保留原准备配体的固定拓扑

默认继续从各姿势独立初始化历史固定评分。新构建支持将原准备配体作为唯一拓扑来源，用每个 MODEL 仅替换评分坐标。若 Meeko 输出的 SMILES、原子索引和显式氢归属完整，可进一步使用经验证的全配体邻接。导出要求一个无 MODEL 的准备配体，核查原子类别、唯一序号、显式氢合法归属、共价距离以及扭转树枢轴边。不支持带虚拟配体位点或拆环宏环类型的输入。邻接文件每行两个 PDBQT 原子序号，要求无向、完整且与导出图一致。各 MODEL 的原子记录、类型、电荷和树必须与同一个准备输入在坐标之外匹配。

```bash
python tools/score_physical_pool.py --export-ligand-bonds --prepared-ligand ligand.pdbqt --output bonds.txt
python tools/score_physical_pool.py --pool poses.pdbqt --receptor receptor.pdbqt \
  --binary /path/to/supporting/vina --prepared-ligand ligand.pdbqt \
  --ligand-bonds bonds.txt --satisfaction --output fixed-topology.json
```

省略 `--ligand-bonds` 时使用同一准备坐标推断的拓扑。提供邻接时，分析不会从每个候选的距离重新推断连接。物理评分以原始 native 排名第一的 MODEL 重算内能为同配体的共同参照，重新计算 `poses[].native` 和所有相依分数，包括 GNINA 共识。历史分和顺序保留在 `poses[].saved_native` 和 `fixed_topology.saved_native_order`。该选项可以和 Vinardo 同用。指定的 Vinardo 可执行文件也会收到固定配体、候选坐标及所给邻接；不支持这些选项的版本会按逐 MODEL 失败规则记录并整池原生回退。未启用固定拓扑时，公开 Vinardo 的原调用方式不变。芳香几何的原生分窗口沿用既有历史基准，不与改变该基准的物理固定拓扑重算混用。

[真实输入检查](../benchmarks/check_fixed_topology_pool.py) 读取 9FA3 与无推断键漂移的固定池，核对默认路径、稳定拓扑、同坐标 GNINA 共识、错误邻接与无支持 Vinardo 的回退。
