# 配体共价键、非键排除与大环闭合

默认解析器修正了四元环上的原子编号依赖。更准确地构建配对表，并未提高本轮盲对接的 2 Å 命中率。两个独立的大环实验保留为针对父版本的 C++ 补丁。一个恢复切断环键两侧完整的三键排除，另一个让相距 20 Å 以上的闭环端点仍受吸引。后者在较低预算下修复了雷帕霉素的两次严重断环，实验坐标恢复仍失败。

[初始协议](PROTOCOL.md) · [发现断环后的第二次实验](TAIL_PROTOCOL.md) · [父版 3663617](https://github.com/HankerWu/fastdock-vina/commit/36636171ae6aef42471b0de5d6b7a90b5f90732e) · [默认核心修复](https://github.com/HankerWu/fastdock-vina/commit/7b9bf89eec156bb452b91ce9b7f6ed62b0ae3b21)

## 默认路径的改动

Vina 搜索通过受体网格和显式配体内非键对算能量及梯度。MC 扰动后做局部 BFGS，接受的候选有条件地进入第二次局优和保存池，经筛选再作显式受体精修。对同一个刚性受体上的单配体，最终原生评分把共同参考姿势的分子内能从总能中减去，再用配体柔性因子归一化。修改配体非键表会改变局优的能量与力；在自参考首姿上，纯配体内能改变量可以抵消。默认使用的径向势、原子类型、网格、BFGS、搜索与输出预算均沿用父版。[原始评分和搜索方法](https://doi.org/10.1002/jcc.21334) · [配对构建与导数](../../src/lib/model.cpp) · [搜索](../../src/lib/monte_carlo.cpp) · [精修及最终计分](../../src/lib/vina.cpp)

给定已解析的成键图，满足“相对距离可动”与最短键距大于三的普通原子对才应进入非键表。旧遍历在第一次到达原子时标记它。若先沿环的长边到达，剩余步数可能不足；稍后从短边到达时又不再展开。新实现按距离层展开，到三层停止，使用原有距离、分型及参数。对闭环标记，在所属原子和直接邻居上预存四个编号位。测试独立验证 BFS 单独修复与加入缓存后的完整有序配对表在 XS/AD 及四种闭环重编号下均相同。缓存为每个有闭环标记的输入临时增加每原子一个字节，无闭环标记时不分配此表。一次性配体准备的实测 CPU 中位数：巴瑞替尼 81.29→58.12 µs，雷帕霉素 460.63→201.61 µs。[化学图与第二骨架](ligands/README.md) · [全部六个配体的时间与消融](preparation/README.md)

- 巴瑞替尼的 208 对中有 6 对实际相隔三根键。用 CCD 作独立核对，修复后 202 对全都合规。仅交换根部的四行原子记录，24 种编码中父版产生最高 0.768108 kcal/mol 的能量偏移和 2.809922 kcal/(mol·Å) 的单个梯度分量偏移；新结果对该置换在浮点精度内一致。考比替尼在独立骨架上复现同一缺陷，288→284 对。罗氟司特与仑伐替尼则各被几何成键推断漏掉一条真实三元环键，修正图遍历仍各有一对误计。[24 种原编码结果](root_parent.csv) · [修复结果](root_fixed.csv) · [逐原子 CCD 对照](ligands/results.csv)
- 初次盲评 10 个配体–受体组合、两档预算。原生/父版/修复的 Top-1、Top-5、输出池 Oracle 成功数分别均为 1/20、1/20、0/20。配体、种子和预算间有重复。JAK2 的唯一成功池在父版的已保存 MC 候选中就有 1.819 Å 姿势；修复版在同一阶段的最好姿势为 7.647 Å。冻结旧池，两种方法仍都把 1.821 Å 的最终姿势排第一。换种子六个受影响条件均为 0/6。考比替尼含实验辅因子/晶格的单列对照中，12,000 档的首姿由 6.343 退至 10.084 Å，两者都没有合格候选。[全部对接、补种子、控制预算与阶段追踪](pose/README.md) · [MEK1 晶体条件](ligands/README.md)
- EGFR 和 p38α 使用 Davis 同测定精确 Kd。在固定种子 740317 下，18 个受体–配体输入的父版与修复版输出和计数逐项相同。101 个固定父版坐标、18 个自参考首姿的分项评分相同。在相差至少三倍的实验分子对中，EGFR 三臂均为 3/24；p38α 上游为 30/37，父版和修复版均为 24/37。两个靶点父版的 Spearman 分别为 −0.905、0.273，Kendall 分别为 −0.786、0.289。各靶点的等大小等电荷子集仅一对，合格配对中无同 Murcko 骨架对。这两个面板主要检验回归，不足以判断热力学改进。18 个生产版重跑均复现研究版 PDBQT 及六项计数。[Davis 等的测定](https://doi.org/10.1038/nbt.1990) · [输入与测定构建体](../search_backtracking/assay/README.md) · [逐对 Kd、相对自由能与分层指标](kd/kd_summary.csv) · [每个固定姿势](kd/fixed.csv)

## 大环上两种机制的独立对照

闭环伪原子的用途是让断开的环在对接构象中重合。搜索可以用断开的扭转树表示环，但实际化学键仍决定两个原子是否属于 1–4 邻居。分别对 Vina 普通单环标签 G0 的 5QCT、1YET、1FKF 和 1FKB 检验。前两者是历史病例；1FKF/FK506 与 1FKB/rapamycin 是另登记的两个 FKBP12 配体。5QCT 教程的准备构象接近晶体，其余三个查询从 CCD 化学图生成。晶体口袋中心按所有方法共享。FK506 和雷帕霉素的配体与受体准备独立重跑分别约 2.17 和 2.56 秒，输入逐字节匹配。[大环方法](https://doi.org/10.1017/qrd.2022.18) · [实验输入来源](../macrocycle_labels/README.md) · [FK506 晶体](https://www.rcsb.org/structure/1FKF) · [雷帕霉素晶体](https://www.rcsb.org/structure/1FKB)

**恢复共价邻接。** 研究臂仅在计算最短键距时恢复同号两个真实闭环碳间的边。四个分子分别从非键表删除 3、4、3、3 对；各自两条闭环约束都保留。对 200 个重原子逐一映射后，从 CCD 或原教程 MOL2 独立证明这 13 条路径均经过被切开的真实键，且均为三根键。固定 16 个父版池的 44 个坐标，统一零参考下 42 个总分发生变化（−0.692 至 +0.029 kcal/mol，打印精度）。16 个池第一名均未改变，16 个自参考首姿分数均相同。三种方法在 16 个条件中均无 2 Å 候选；父版和研究臂分别用 127,671、127,934 次联合求值。1YET 的第二枚种子在高预算下从 5.763 退至 12.514 Å。[研究补丁](closed_graph.patch) · [逐对真实路径](chemical_pairs.csv) · [分体系、预算、RMSD 与开销](macro/structures.csv) · [固定池](macro/fixed.csv)

**延续长距离闭环力。** 父版线性闭环项在 20 Å 跌至零，在截断外做局优没有这项恢复力。1FKB 的 3000 档父版首姿两次给出 23.447 和 22.947 Å 的真实闭环碳距，均无配体重原子距受体小于 2 Å。研究臂仍用同一个闭环权重 w，在最后两个受截断影响的采样格及其外取 E = wr，沿距离的导数为 w，笛卡尔梯度系数为 w/r。原有正能量软化使用 C(E)=vE/(v+E)，力的乘数为 v²/(v+E)²。原点沿用原表，不引入除零。只对实际列入的 XS 闭环对作此处理，不给整套重原子打分添加远程项。有限差分覆盖两类软化、三个权重和接点两侧等共 42 次，最大径向误差 0.000182 kcal/(mol·Å)。普通 Vina 默认权重和 8 Å 普通相互作用截断是本实验的验证范围。[C++ 实验补丁](long_glue.patch) · [数值差分](tail/numeric.csv)

| 1FKB，同种子同输入 | 父版 3000 档 | 延续力 2000 档 |
|---|---:|---:|
| 740317，闭环碳距 / 整分子 RMSD，Å | 23.447 / 12.845 | 1.579 / 9.591 |
| 740317，实际总求值 / 秒 | 3436 / 0.179 | 2973 / 0.172 |
| 191009，闭环碳距 / 整分子 RMSD，Å | 22.947 / 11.836 | 1.525 / 13.655 |
| 191009，实际总求值 / 秒 | 3361 / 0.174 | 2199 / 0.136 |
| 458921，闭环碳距 / 整分子 RMSD，Å | 1.562 / 10.855 | 1.562 / 10.855 |

等名义预算的 24 个条件中，各臂共输出 67 个候选，真实闭环碳距大于 3 Å 的候选从 5 个降至 1 个，大于 20 Å 的从 3 个降至 0。剩余一个为 FK506 的 4.271 Å。过长闭环减少的直接实例来自雷帕霉素这一种配体。两臂的 2 Å 候选均为 0/24；父版/研究臂分别用 193,564/193,071 次联合求值。冻结父版 67 个坐标复算最终分数，闭环延续仅参与搜索/局优，候选最终打分维持原有定义。[全部首姿与成本](tail/structures.csv) · [下调预算的全矩阵](tail/cost/structures.csv) · [所有模式与闭环残差](tail/models.csv)

从相同断环坐标单独做严格 Armijo 局优，先核对重解析后的有序非键表与原输入逐字节相同。三条父版超过 20 Å 的断环均停在起点，延续力使其分别收至 13.349、1.566 和 4.414 Å，并未每次闭合。BACE1 共晶起点作为阴性对照，两臂均保持输入坐标。原有搜索中的正能量软化会削弱远距离回复力，有限扭转坐标和回溯也影响运动；这些数据不足以将剩余失败唯一归因于其中一项。下一项可检验的问题是：在同一扭转构象上构造满足闭环几何的提议，能否同时降低断环比例、保持相同实际开销，并提高有实验图证实的正确姿势覆盖。两项大环研究臂均未进入默认算法。[每个起点和求值](tail/local.csv)

亲和力另用 [Gal-1/Gal-3 的 MMFF 固定池迁移](../frozen_core/transfer/galectin_mmff/README.md) 验证。13 个化合物、208 个池的实验显示，真空 MMFF94s 的完整替换在 Gal-3 降低同测定配对准确率，进一步限制了从 BACE1 单一系列外推的范围。

## 复现

从本分支完整检出，在仓库根目录执行。重算归档结果需要 Python、NumPy、SciPy、RDKit 和 Gemmi。制备 FKBP 配体/受体另需 Meeko 0.7.1、PDBFixer、OpenMM，历史配体协议见上面的链接。生产代码使用 GCC 14.2、C++14、-O3 和 Boost 1.83。每个对接臂固定单核、单轨迹、相同搜索盒和种子，上游计数来自与无插桩输出逐字节相同的复跑。单次墙钟时间包括建图、运行和写出，具有宿主负载噪声。

```bash
D=benchmarks/bond_distance
bash "$D/replay.sh"
```

重新编译各研究臂并运行时使用以下命令。系统安装 Boost 的情形无需额外选项；自定义安装时给 build.py 添加 --boost-prefix。CPU 编号可用 run.py 的 --cpu 覆盖。输出写入新的临时目录。

```bash
D=$(pwd)/benchmarks/bond_distance
W=$(mktemp -d)
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
for arm in parent bfs closed tail; do
    git worktree add --detach "$W/$arm" 36636171ae6aef42471b0de5d6b7a90b5f90732e
done
git -C "$W/bfs" apply "$D/preparation/bfs_only.patch"
git -C "$W/closed" apply "$D/closed_graph.patch"
git -C "$W/tail" apply "$D/long_glue.patch"
for arm in parent bfs closed tail; do
    python "$W/$arm/tools/build.py" --output-dir "$W/bin-$arm"
done
python benchmarks/pair_gradient/build_arms.py "$W/up" --arms upstream upstream-counts
mkdir "$W/macroinput" "$W/kdinput"
tar xJf benchmarks/macrocycle_labels/observations.tar.xz -C "$W/macroinput"
tar xJf "$D/macro/extra_inputs.tar.xz" -C "$W/macroinput"
tar xJf benchmarks/search_backtracking/assay/inputs.tar.xz -C "$W/kdinput"
U=("$W/up/upstream/build/vina" "$W/up/upstream-counts/build/vina")
for arm in closed tail tail-cost; do
    bin=$arm
    study=$arm
    if [ "$arm" = closed ]; then study=macro; fi
    if [ "$arm" = tail-cost ]; then bin=tail; fi
    python "$D/run.py" "$study" "$W/macroinput" "$W/out-$arm" "${U[@]}" "$W/bin-parent/vina" "$W/bin-$bin/vina"
    python "$D/analyze.py" macro "$W/macroinput" "$W/out-$arm"
done
python "$D/fixed.py" macro "$W/macroinput" "$W/out-closed" "$W/bin-parent/vina" "$W/bin-closed/vina"
python "$D/fixed.py" macro "$W/macroinput" "$W/out-tail" "$W/bin-parent/vina" "$W/bin-tail/vina"
python "$D/run.py" kd "$W/kdinput" "$W/out-kd" "${U[@]}" "$W/bin-parent/vina" "$W/bin-bfs/vina"
python "$D/analyze.py" kd "$W/kdinput" "$W/out-kd"
python "$D/fixed.py" kd "$W/kdinput" "$W/out-kd" "$W/bin-parent/vina" "$W/bin-bfs/vina"
```

[独立小环姿势与阶段跟踪](pose/README.md) · [真实键图、MEK1 环境与单列病例](ligands/README.md) · [配体准备消融、探针构建方法](preparation/README.md)。局优重放用 local_tail.py，依次传入大环输入目录、out-tail、输出目录、父版程序、延续力程序以及按 preparation/README.md 编译的父版 probe。tail_numeric.cpp 链接延续力臂的对象文件，使用仓库环十二烷 fixture 检验差分。
