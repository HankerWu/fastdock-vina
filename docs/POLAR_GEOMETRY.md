# Late 候选池的极性几何重排

用 `--polar_geometry distance`、`directional` 或 `penalty` 选择最终排序规则。三种模式共用原 Late 轨迹、候选选择、去重、显式受体精修和原 Vina 能量。省略该选项即可使用原有行为。

```bash
build/polar/vina --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
  --pose_selection late --polar_geometry directional \
  --seed 42 --cpu 1 --exhaustiveness 8 --num_modes 9 --energy_range 100 --out polar.pdbqt
```

入口支持一个刚性受体和一个配体的常规 Vina Late 对接。模板、Gauss2 重排、预生成网格、灵活受体、批处理和单点模式不属于此入口。设置 Armijo 时，两个对照应使用相同精修设置。

## 几何和排序

氢键两端限于引擎 XS 判为有供氢或受氢能力的真实 N/O。配体和蛋白都可扮演两种角色。受氢原子 Ac 必须有引擎共价拓扑中的同分子重邻 B。氢、伪原子和金属均不作 B。供氢原子 D 到 Ac 的距离 d 必须满足 0 < d ≤ 4.0 Å。径向权重为 `exp(-0.5 × ((d - 2.9) / 0.45)²)`。方向项从所有 B–Ac–D 角中取最小角 θ，门控为 `clamp((cos 80° − cos θ) / (cos 80° − cos 110°), 0, 1)`。

对每个配体极性原子，在所有匹配伙伴和两种角色上分别取最大径向权重与最大方向权重。方向权重等于径向权重乘门控值。对配体各原子的最大值分别求和，得到径向 R 和方向 A。

| 模式 | 最终排序分，越小越优先 |
| --- | --- |
| `distance` | 原 Vina 分 − [gV(R) − gV(0)] |
| `directional` | 原 Vina 分 − [gV(A) − gV(0)] |
| `penalty` | 原 Vina 分 + [gV(R − A) − gV(0)] |

gV 是同一配体的原 Vina 构象无关变换。系数 1.00 与几何参数固定。圆括号内的 R 和 A 是无量纲几何项。三种排序分用于在同一候选池选择姿势。并列使用原名次。没有合格配对时，几何项为零。

## 输出

| PDBQT 行 | 字段 |
| --- | --- |
| `REMARK VINA POLAR RESULT:` | 所选模式排序分、相对第一条输出的 RMSD 下界和上界 |
| `REMARK VINA NATIVE:` | 原 Vina 分、原生排序名次 |
| `REMARK VINA POLAR SCORES:` | `distance` 排序分、`directional` 排序分 |
| `REMARK VINA POLAR TERMS:` | 未经过 gV 的 R、A |
| `REMARK VINA POLAR MODE:` | `distance`、`directional` 或 `penalty` |

三种模式均保留相同的数值诊断。原 Vina、距离和方向分可在同一份输出上分别以原名次破平复算选择。Vina 在普通配体输入上的线性变换也允许用原 Vina 分加上方向分与距离分之差复算角度惩罚选择。数值字段输出双精度。原 Vina 能窗仍以池内最优原 Vina 分为基准，`MODEL` 则跟随当前模式的排序。输出接着列出原 Vina 的 `INTER`、`INTRA`、`INTER + INTRA` 与 `UNBOUND`。

## 核验

```bash
python3 tools/build.py --output-dir build/polar --jobs 4
PYTHONDONTWRITEBYTECODE=1 python3 benchmarks/check_polar_geometry.py \
  build/polar/vina /path/to/verified-pre-polar/vina
```

系统未提供已编译的 Boost 时，可用 `--boost-prefix` 指定现有的 Boost 安装目录。对照与复验应记录可执行文件的 SHA-256。

`benchmarks/check_polar_geometry_numeric.cpp` 通过真实 PDBQT 解析和拓扑构建测量距离边界、角度、双向饱和及原子资格。可用同一套编译选项和 Boost，将其链接到 `build/polar/obj` 中除 `main.o` 之外的对象。给可执行文件传入临时受体 PDBQT 路径。
