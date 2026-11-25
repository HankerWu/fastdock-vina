# 极性埋藏姿势重排

`--polar_burial burial|geometry|compensated` 用于单配体和显式刚性受体的 Late Vina 对接。搜索、原生显式精修、去重、`--num_modes` 截断和原生未结合能参照均由现有流程确定。可沿用 `--refine_candidates` 扩大精修输入。新增评分只对保留的构象计算，使用原生分数指定 `--energy_range`，同分时沿用原生先后顺序。与 Armijo、模板引导、高斯重标定、旧极性几何重排和轨迹补存互斥。

## 评分

对每个供体或受体配体 N/O 位点，统计 5.5 Å 内受体重原子。其埋藏值 `b` 为各重原子的 `clamp((5.5 - r)/2, 0, 1)` 之和除以 8，然后上限截为 1。沿用[极性几何](POLAR_GEOMETRY.md)的配对规则，分别取最佳径向质量 `d` 和最佳角向质量 `q`，无可评分伙伴时两者为 0。三个原始项为 `B = Σb`、`G = Σ(d-q)`、`U = Σb(1-q)`。因此氢键伙伴最多补偿该位点新增的全部埋藏罚，原生 Vina 氢键项保留。三臂分别对原生间能与内能之和扣除共同参照后的值加入 `0.6 B`、`0.6 G` 或 `0.6 U`，再做原生柔性归一化。该线性归一化等价于对原生报告分数加上同样归一化后的增量。罚项反映局部原子密度与可用的直接氢键几何。

PDBQT 保留三位小数的 `REMARK VINA RESULT` 和通常的原生能量分解。`REMARK VINA NATIVE` 记录全精度原生能量和原生名次。`REMARK VINA POLAR BURIAL RESULT` 记录选择的排序分数与相对输出首位的 RMSD 上下界。`REMARK VINA POLAR BURIAL TERMS` 按 B、G、U 顺序记录原始项，`REMARK VINA POLAR BURIAL DELTAS` 按同一顺序记录归一化增量。另有模式和原生名次备注。

## 固定坐标

```bash
vina --receptor receptor.pdbqt --ligand pose.pdbqt --score_only --polar_burial compensated --cpu 1 --verbosity 0
```

从多模型输出提取单个构象时去掉 `MODEL` 和 `ENDMDL`。固定坐标模式使用文件中的坐标，直接计算显式受体相互作用和配体内能，无搜索、优化、搜索盒罚或亲和力网格生成。任意一个极性埋藏模式值都输出同样的三行十二位小数诊断。

| 前缀 | 字段顺序 |
| --- | --- |
| `Polar burial terms:` | B、G、U |
| `Polar burial deltas:` | 埋藏、方向缺陷、补偿后增量 |
| `Polar burial native:` | 显式间能、配体内能、共同参照为 0 时归一化后的两项之和 |

对已写出的原生候选池，取文件中既有原生结果作为基底，加上固定坐标诊断增量再排序。原始 PDBQT 坐标和默认亲和力备注的截断精度可能影响近似并列的排名。

## 验证入口

在引擎根目录运行 `python3 benchmarks/check_polar_burial.py build/polar-burial/vina build/refine-pool/vina`。测试用可解析的几何构型核对三项罚分，并用 1IEP 验证原生输出、候选池和预算一致性、柔性归一化及固定坐标复算。
