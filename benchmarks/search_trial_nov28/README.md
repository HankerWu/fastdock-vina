# 角度限幅搜索在两个新 BRD4 受体上的外检

两个合格的 BRD4 第一溴结构域受体来自 4BJX 和 5U2E。候选在两组外检中都降低了平均池内最优 RMSD，同时提高了实际首姿误差。八格中，候选相对轮前首姿两格改善、六格退化，八个首姿都超过 5 Å。三种结构预先登记，4WIV 在制备时因 Meeko 无法建立 GLN A:59 失败，已保留原始记录。两结构和骨架对本轮研发外留，溴结构域家族已被项目开发使用。[事前登记](PROTOCOL.md) · [三份原始沉积、配体字典与输入](external_inputs.tar.xz) · [制备状态与哈希](external_preparation.json)

| 外检条件，每臂 4 格 | 官方 / 父版 | 轮前 Pro | 候选 Pro | 候选减轮前 | 候选减官方 / 父版 |
|---|---:|---:|---:|---:|---:|
| 同晶平均首姿 (Å) | 4.325 | 6.373 | 7.449 | +1.076 | +3.124 |
| 同晶平均池最优 (Å) | 3.107 | 3.893 | 2.889 | −1.004 | −0.218 |
| 同晶平均排序损失 (Å) | 1.218 | 2.480 | 4.560 | +2.080 | +3.342 |
| 交叉平均首姿 (Å) | 6.851 | 6.579 | 8.150 | +1.571 | +1.299 |
| 交叉平均池最优 (Å) | 4.652 | 5.379 | 5.274 | −0.105 | +0.622 |
| 交叉平均排序损失 (Å) | 2.199 | 1.200 | 2.876 | +1.676 | +0.677 |

5U2E 同晶、种子 85237 最能区分两个终点。轮前首姿为 1.378 Å。候选生成了 1.211 Å 的模式，分数 −7.235，排名第二；首位为 7.190 Å，分数 −7.549 kcal/mol。4BJX→5U2E、种子 85237 则两端共同恶化，首姿从 5.584 升至 8.793 Å，池内最优从 4.965 升至 6.066 Å。搜索改变了供终排选择的实际池，终排函数在两个版本中相同。八格均值不能估计整个家族上的泛化误差。[逐格原始读数和计数](results/cells.csv) · [288 个模式的误差与分数](results/ranks.csv) · [分位数、严重误差及成本](results/summary.csv) · [逐格差值](results/paired.csv)

在实际 MC 中，未缩放的大旋转可能沿指数映射卷绕多圈；扭转坐标逐项对 2π 取模。候选从所有刚体旋转向量范数和所有扭转增量绝对值的最大值 A 计算第一试步。如果 A>π，令 α₀=π/A；否则令 α₀=1。α 同时乘于平移和全部角度，因而也可能缩短所需的平移。十次试步、接受准则、按实收步长更新的正曲率保护及失败恢复沿用既有代码。参数 π 是角坐标的几何界限。显式精修、打分势、软化、温度、种子和束缚体积均沿用基线。该局优算法作用于干受体经验能量及 MC 提议，不建立可用于平衡系综的采样权重。[候选代码 0bb827e](https://github.com/HankerWu/fastdock-vina/commit/0bb827ebfe69a32504b74488730f80a21495eccf)

两受体晶体号均映射至 UniProt O60885，实验残基编号与 UniProt 编号相符。跨受体参考经 36 个匹配口袋 Cα 配准，拟合误差为 0.673 Å。这些跨受体坐标是转移参考。评分时，两种配体分别枚举全部 4 个保留立体化学的图映射。独立核查从原始 PDB 和 CCD 重建分子图，按发布的 SMILES 索引解码姿势，采用 SVD 求旋转，对 288/288 个模式重算几何及分数；与第一次计算的最大 RMSD 差为 1.07×10⁻¹⁴ Å。全部 32 个命令独立再执行一次，32 个池逐字节一致，全部排名表逐字节一致。[独立几何核查](results/verification.json) · [原始发布和日志](results/raw.tar.xz) · [重放命令](results/replay_commands.json) · [重放成本和哈希](results/repeat.csv)

三份沉积在首轮制备共用 6.86 CPU 秒、34.07 墙钟秒（含下载）。使用固定沉积离线重做，两个合格受体和配体输入逐字节匹配，4WIV 在同一步失败，逐病例有 CPU、墙钟和峰值内存。外检的 8 次搜索每臂合计 CPU 秒为官方 17.829、父版 6.084、轮前 3.388、候选 3.388。两 Pro 的 hunt / 二次局优 / 显式精修联合求值总数分别为 269601 / 18788 / 1661 和 262181 / 26186 / 1821。官方和父版没有原生计数。以上含候选发布，使用已有的刚性受体输入，没有追加后处理重排。外检计时只支持本轮已测预算和机器配置。[离线制备复算](results/preparation_replay.json) · [四臂可执行文件哈希](results/execution.json)

## 复现

在仓库根目录运行，以下环境变量指向四个相应的真实源码构建。候选使用上述研究提交。普通 Python 环境需 NumPy、SciPy；原始化学输入制备和独立核验另需 RDKit、Meeko、Gemmi、OpenMM 和 PDBFixer。使用新的工作目录。

```bash
D=benchmarks/search_trial_nov28
mkdir -p "$WORK/inputs"
tar xf "$D/external_inputs.tar.xz" -C "$WORK/inputs"
OPENBLAS_NUM_THREADS=1 python "$D/external.py" "$WORK/inputs" "$WORK/run" \
  --official "$OFFICIAL" --parent "$PARENT" --stable "$STABLE" --candidate "$CANDIDATE"
OPENBLAS_NUM_THREADS=1 python "$D/verify_external.py" "$WORK/inputs" "$WORK/run"
OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 python "$D/replay_preparation.py" "$WORK/inputs" "$WORK/reprepared"
# 从公共原件完整重做制备（含下载）。
OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 python "$D/prepare_external.py" "$WORK/downloaded"
```

数值回归使用研究提交的头文件及同源码的 quaternion/random 对象；研究实现由 Git 祖先保留。多配体及柔性侧链的限幅、纯平移、成功或耗尽的实际返回步长、原有 Armijo 状态恢复和正曲率保护均已通过。

```bash
c++ -O3 -std=c++14 -I"$CANDIDATE_SRC/src/lib" -I"$BOOST/include" "$D/check_chart.cpp" \
  "$CANDIDATE_OBJ/quaternion.o" "$CANDIDATE_OBJ/random.o" -o "$WORK/check-chart"
"$WORK/check-chart"
c++ -O3 -std=c++14 -I"$CANDIDATE_SRC/src/lib" -I"$BOOST/include" benchmarks/check_armijo_numeric.cpp \
  "$CANDIDATE_OBJ/quaternion.o" "$CANDIDATE_OBJ/random.o" -o "$WORK/check-armijo"
"$WORK/check-armijo"
```
