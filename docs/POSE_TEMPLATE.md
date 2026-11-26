# 外部共晶配体引导的姿势选择

模板应来自不同的共晶配体，且预先对齐到本次刚性受体的坐标系。文本每行四列，用空格或制表符分隔，依次为大小写准确的元素符号和三个 Å 坐标。示例写两行，实际可使用参照配体的全部重原子。

```text
C 16.733 51.911 14.850
N 15.729 51.744 15.790
```

```bash
build/template/vina --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
  --pose_selection late --pose_template aligned_reference.txt --template_stage rank \
  --seed 42 --cpu 1 --exhaustiveness 8 --num_modes 9 --energy_range 100 --out template_rank.pdbqt
```

`--pose_template` 和 `--template_stage` 成对提供。下述三个元素重叠阶段支持一个刚性受体、一个配体、直接生成搜索网格、Vina 力场和原生显式受体精修。命令支持已有 Vina 权重参数。同一个对照使用相同权重和搜索预算。这三个阶段与 `--batch`、柔性受体、预计算网格、既有 Gauss2 重排、Armijo 精修、禁用精修、单点打分或局部优化模式互斥。

## 重叠项和三个作用阶段

每个重原子取另一组中最近的同元素重原子，记距离为 $d$。该原子的贡献为 $\exp[-d^2/(2\times1.5^2)]$。一个方向以出发组的全体重原子数为分母，无同元素匹配贡献零。两组相互计算后取两个方向的均值 $O$。有效模板和配体各至少要有一个真实重原子。PDBQT 中的氢、水伪原子 `W` 和大环闭合伪原子 `G0–G3` 排除，真实闭环碳 `CG0–CG3` 按 C。模板文本中的元素符号 `W` 表示钨元素。

| `--template_stage` | 送入精修前 | 显式受体精修 | 精修后去冗余优先级 |
| --- | --- | --- | --- |
| `rank` | 原 Vina 原子能量 | 原 Vina BFGS | 原 Vina 原子能量 |
| `select` | 原 Vina 原子能量减去 $6O$ | 原 Vina BFGS | 原 Vina 原子能量 |
| `refine` | 原 Vina 原子能量减去 $6O$ | 同一 $-6O$ 及其 Cartesian 梯度加入原 BFGS 目标 | 原 Vina 原子能量减去 $6O$ |

三个阶段均使用原 Monte Carlo 搜索和同一无模板的轨迹内局部优化。送入显式精修的名额、精修步数上限、盒外重试上限和两轮去冗余距离都使用原 Late 设置。精修后保留池的共同内能参照取原 Vina 原子能量最低的保留姿势。模板最后重排分为原 Vina 亲和力减去 $g(6O)-g(0)$，其中 $g$ 是对该配体应用的原 Vina 构象无关变换。重排分并列时采用原 Vina 配体池内名次。最近匹配在恰好同距时按首次遇到的项确定，精修使用该分段的梯度。

## 输出

- PDBQT 的 `REMARK VINA TEMPLATE RESULT` 依次给出模板排序分、相对本次输出首姿势的 RMSD 下界和上界。 `REMARK VINA NATIVE` 给出原 Vina 亲和力和保留池内原生名次，`REMARK VINA TEMPLATE OVERLAP` 给出 $O$，`REMARK VINA TEMPLATE STAGE` 给出实际阶段。计算模板排序分使用未舍入的原生亲和力和重叠值。模板排序分表示几何先验下的优先级，亲和力保留在原生分数字段中。
- `MODEL` 编号跟随模板排序，`--energy_range` 按保留池中最佳原 Vina 亲和力截断。这使 `rank` 及对应的无模板 Late 对照写出同一套候选。 `INTER`、`INTRA`、`INTER + INTRA` 和 `UNBOUND` 记录在实际输出坐标上重新计算的原 Vina 能量分量。
- 普通详略级别报告 Late 收集、送精修、显式精修调用和去冗余后的数量。 `Template pose scoring` 另报候选初筛、BFGS 目标评估及精修后用于排序的重叠计算次数。 `refine` 为每个已精修候选算一次重叠，再把已算值用于保留姿势的最终排序。

## 核验

```bash
python3 tools/build.py --output-dir build/template --jobs 4
PYTHONDONTWRITEBYTECODE=1 python3 benchmarks/check_pose_template.py \
  build/template/vina /path/to/verified-pre-template/vina
```

数值检查 `benchmarks/check_pose_template_numeric.cpp` 使用同一 Boost 配置链接
`build/template/obj/` 的 `cache`、`grid`、`coords`、`conf_independent`、`non_cache`、`szv_grid`、`model`、
`parse_pdbqt`、`quasi_newton`、`quaternion`、`random` 和 `utils`。运行该检查时传入
`examples/basic/1iep_ligand.pdbqt`，它检验双向归一、元素匹配、Cartesian 与姿势自由度数值梯度，以及 BFGS 对已知极小值的收敛。

### 按受体接触重排

为比较一种配体与另一种已知结合配体在同一受体上的接触，可用 `--template_stage contact_rank`。模板仍是受体坐标系中的元素和重原子坐标文件。默认 Native 和 Late 均支持该选项，也可沿用已有的显式 Armijo 或 geometry 精修。这个选项只改变精修后的排序。候选生成、局部优化、Vina 原生亲和力和原生能量窗使用原流程。

受体上每个重原子记录它到最近配体重原子的距离。4 Å 以内取 1，之后线性下降，6 Å 时归零。另给受体 N/O 记录到最近配体 N/O 的距离，在 3 至 4.5 Å 内线性下降。这是极性原子接近程度，未判断氢键方向。分别计算两类接触在查询与独立模板之间的连续 Tanimoto 相似度，取非空通道的平均。只有查询产生的接触也进入分母。保留原有模板排序的强度和扭转归一化，分数并列时按原生名次排序。不同结合模式、不同亚口袋或错误的模板对齐都可能削弱这个先验。

[同候选原生、元素模板、接触模板对照及复算方法](../benchmarks/frozen_core/contacts/README.md)
