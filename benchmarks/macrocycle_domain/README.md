# 网格边界应约束哪些原子

以前的边界检查在两条相反的路径上出错。Vina 与 Vinardo 的格点能量、显式受体能量及其梯度都跳过 G0–G3；这些点随配体移动，用于恢复开环键的几何关系。格点和显式精修的盒内判断却纳入了它们。另一个问题出现在 AD4，公开的固定打分与局部优化入口检查了未装载的 Vina 格点对象；三个方向的体素数都为零时，所有坐标都会通过。旧版在局部优化后又使用这张空图进行最终打分的域检查。前一个问题是假拒绝，后一个问题是假接受。闭合标记仍然参与内部势和梯度，未改变权重或搜索算法。[大环闭合方法原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10392634/)

模型现在共用一个盒内判据，保留原有每轴 0.0001 Å 容差。计入重原子、闭环碳 CG0–CG3 以及具有受体图的水探针 W，跳过氢及无受体图的 G0–G3。AD4 的输入和优化后检查使用实际装载的 AD4 格点。Late 为 PDBQT 三位小数边界作微平移时，也使用同一组受边界约束的原子。每次边界检查仍是 O(可动原子数)，遇到越界可提前退出。没有附加能量求值，BFGS 及 Metropolis 的目标与随机数均保持原实现。

## 冻结姿势的反事实

从归档初代 Pro 输出依文件名、排名取第一个闭环标记超出实体原子包围盒至少 0.05 Å 的姿势，直接使用原 PDBQT 的单姿坐标和原始配体拓扑。分别用 30 Å 盒包含全部原子、将一个盒面放在真实重原子极值和标记极值的中间、再把该面向内移到至少一个真实原子越界 0.25 Å。姿势、受体和物理相互作用完全相同。这是构造盒面的因果诊断。

| 人源结构 | G 超出实体包围盒 | 正常盒中的显式分数 | 仅 G 越界，父版 → 修正版 | 实体原子越界 |
|---|---:|---:|---|---|
| HSP90α，3VHA | 0.616 Å | −1.942 | 错拒 → −1.942 | 均拒绝 |
| 组织蛋白酶 D，6QBH | 0.818 Å | −5.136 | 错拒 → −5.136 | 均拒绝 |
| CDK2，8H6P | 无符合选取标准的归档姿势 | — | — | — |

在人 FKBP12 1FKH 实测系列的化合物 10 上，把已对接的原生 AD4 首姿沿 x 轴整体移出 100 Å。旧版固定打分返回成功和 3,332,132,362.570 kcal/mol，局部优化返回成功和 1,942.258 kcal/mol。修复后两个入口均报告越界。原图内，两版固定打分与局部优化均为 −7.462 kcal/mol。这里评估的是合法域，不能解释为改善了能量函数。[逐次运行](boundary.csv)

## 同预算对接与实测 Ki

三个大环体系、两个种子和两档预算，每种算法 12 个配置。每次一轨迹、一线程、五个名额、原始 26 Å 盒。所有 12 对父版和修正版输出逐字节一致。两版均输出 35 个姿势，实际 hunt、第二次最小化和显式精修的求值数分别为 66,703、24,477 和 3,935。选择原生 Vina 的 33 个输出作第三组，实际求值分别为 66,703、24,477 和 6,717，原生计数版的输出与上游可执行文件逐个比对。共同受体坐标系下，用沉积结构坐标和配体 SDF 分别作立体化学对称映射。三种算法各自的 Top-1、Top-5 和输出池 oracle 都是 0/12，无法估计找到正确候选时的排序成功率。重试、省时和精度普遍改善均未得到这组数据的支持。父版/修正版总墙钟 2.650/2.600 s，CPU 2.629/2.589 s，单次交错测量不能证明加速。[全部逐例求值和时间](structures.csv) · [每个姿势的独立 SDF 复算、碰撞和闭环残差](models.csv)

在两个固定姿势诊断盒中从头搜索，另外八对父版/新版本的 19 个输出也逐字节一致，两版显式精修各 27 次、2,978 次求值，Top-1/5 同为 0/8。改变搜索起点后，盒面测试中人工挑出的固定姿势不再保证会出现。这解释了“同一坐标的合法域确实改善”和“从头对接没有观察到收益”能够同时成立。[全部压力对接](stress_structures.csv)

固定 Holt 等人的八个中性 FKBP12 抑制剂在 1FKG、1FKH 中的 AD4 首姿，使用相同原格点复算三个种子的全部 48 个分数，两版分数完全相同。每个受体先取三种子的中位数，再对两个构象取等权平均。至少 3 倍且区间不重叠的唯一分子对为 15 对，父版和修正版各判对 6 对。七个点测定的 Spearman −0.179，Kendall −0.143。单个受体分别为 6/15 和 8/15，两个受体和三个种子不增加实验分子数。这只验证了有效域内同一 AD4 评分的可复现性。[实测系列、测定条件和数据来源](../cubic_grid/fkbp_ad4/README.md) · [固定首姿](fixed_affinity.csv) · [逐受体统计](fixed_affinity_stats.csv)

## 候选去重与闭环标记是另一个问题

在去除整文件重复后，扫描了归档父版输出中 39 个含至少两个姿势的文件，共 132 个文件内姿势对。1 Å 处有一对因纳入标记而跨过原子同序 RMSD 阈值，5K0X 的第一、二名从剔除 G 后的 0.9492 Å 变成 1.0147 Å。2 Å 处无跨越。这是已经释放的有限候选池，不能由此估计在线去重、保存或姿势排名的收益。没有把更改 native RMSD 定义与本次合法域修复绑定。[冻结距离](archived_distances.csv)

## 复现

从仓库根目录运行。父版为 3d9dfb9，原生对照为上游 3c65c0b。输入来源、制备和晶体坐标哈希见 [原结构记录](../concerted_ring/results/sources.json)。使用支持 NumPy、SciPy 与 RDKit 的 Python。将 CPU 编号设为当前可绑定的一核。如需指定 Boost 前缀，给构建命令追加 `--boost-prefix`。

```bash
D=benchmarks/macrocycle_domain
W=$(mktemp -d)
git worktree add --detach "$W/parent-src" 3d9dfb9
python "$W/parent-src/tools/build.py" --output-dir "$W/parent" --jobs 2
python tools/build.py --output-dir "$W/fixed" --jobs 2
python benchmarks/pair_gradient/build_arms.py "$W/up" --arms upstream upstream-counts
mkdir -p "$W/input"
tar -xf benchmarks/concerted_ring/new_inputs.tar.xz -C "$W/input"
p="$W/parent/vina"
f="$W/fixed/vina"
u="$W/up/upstream/build/vina"
uc="$W/up/upstream-counts/build/vina"
python "$D/compare.py" "$W/input" "$W/panel" "$u" "$uc" "$p" "$f"
python "$D/boundary.py" "$W/input" "$W/boundary" "$p" "$f"
python "$D/compare.py" "$W/input" "$W/stress" "$u" "$uc" "$p" "$f" --diagnostic-boxes "$W/boundary"
python "$D/fixed_affinity.py" "$W/affinity" "$p" "$f"
python benchmarks/check_local_only.py "$f" "$u"
```

[原始标准和压力对接坐标、日志、实际计数与运行配置](raw.tar.xz) · [运行所用二进制 SHA256](binaries.json)
