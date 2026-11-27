# 大环闭环标记如何改变姿势距离与候选保留

切开环键可以用树表示大环自由度。两端的真实碳 CG 与随另一端移动的 G 构造点之间有闭环势。G 记录断口应在的位置，实验分子没有与 G 对应的核。若用全部非氢点构造姿势 RMSD，就额外加权了断口附近，还把未闭合产生的 G 偏移计为分子的空间位移。原生轨迹存储、合并、精修后去重和打印的上下界都受其影响。[Holcomb 等的大环方法](https://doi.org/10.1017/qrd.2022.18) · [Vina 官方教程](https://autodock-vina.readthedocs.io/en/stable/docking_macrocycle.html)

设两个姿势中参与原生距离的真实可动重原子及水合位点为 S，闭环构造点为 G，坐标索引一致。旧平方距离是两个集合的均方位移按点数加权的平均。修复以 S 上的均方位移开方作为固定索引距离。零个符合条件的位点时距离为零。打印的下界继续用同种元素的最近点，取两个方向中较大值。精确化学对称映射另用于与实验结构比较。G0–G3 仍参与闭环势和输出；CG、可动受体重原子和原生水合 AD4 的 W 均保留。后期选择功能原有的纯化学几何会排除 W。这个定义区别使独立的水合位点保持可辨。

## 冻结构象的因果证据

b991eb7 和修正版共 18 个配对条件，三个结构、两个种子、三个预算。每个第二次局优事件都记录未量化至 PDBQT 精度的能量和全部非氢坐标（含 G）。两臂在这 376 次观察上逐项相等，预算间共享轨迹前缀。用记录的事件独立重建容量 5 的候选容器，每一步的保留能量列表都与相应原生程序一致。观察版生成的 36 份 PDBQT 及计数逐份匹配生产版。反事实比较给出了两次可重复定位的首个分歧（以下 MC 步从零计）。[原始轨迹](traces.tar.xz) · [标准预算](standard/trace.csv) · [高预算](higher/trace.csv)

- HSP90α 3VHA，种子 458921，第 11 与 12 步的候选在旧距离相差 1.086914 Å，在新距离相差 0.857520 Å。轨迹去重半径为 1 Å。第 12 步新候选替换第 11 步候选，释放一个名额保留较早候选。此分歧出现于三个预算；最终首姿和池内最佳 RMSD 相同。
- 新结构 HSP90α 4EGK，种子 191009、12000 档，第 49 和 90 步候选在原生合并时由 2.000188 Å 变为 1.990015 Å，跨过 2 Å 门槛。相同的 26 个输入事件，输出由 4 个变为 3 个，显式精修从 540 次求值降至 417 次。保留池最佳与首姿都为 12.671 Å。该差别依赖全精度坐标；打印为三位小数会掩盖旧距离距门槛仅 0.000188 Å 的事实。

## 真实分子、计算量和负结果

| 结构 | 来源和输入质量 | 48000 档两个种子各自最接近晶体的 RMSD，Å（修复版） |
|---|---|---|
| 3VHA 人 HSP90α | 历史准备，1.39 Å，27 个重原子 | 10.153、14.899 |
| 4EGK 人 HSP90α–radicicol | 本轮新增，1.69 Å，25 个重原子，一种化学映射；最近重建受体原子距配体 13.19 Å | 11.259、3.672 |
| 2J9M 人 CDK2–PY8 | 本轮新增，2.50 Å，22 个重原子，两种化学映射；无需重建受体重原子 | 7.280、11.568 |

查询配体从 CCD 经独立构象生成。采用已知共晶口袋中心，各方法共享同一准备和 26 Å 盒。两项新结构均通过预先写明的化学图和 8 Å 环境检查，准备计时约 7.28、6.14 秒，涵盖当次数据读取或下载。没有借用配体晶体坐标初始化。实验姿势同时用沉积 PDB 的原子映射和配体 SDF 的立体图映射独立复算。三种方法合计 177 个输出，两种参考得到的化学 RMSD 一致至线表精度。[3VHA](https://www.rcsb.org/structure/3VHA) · [4EGK](https://www.rcsb.org/structure/4EGK) · [2J9M](https://www.rcsb.org/structure/2J9M) · [输入、检查与结构来源](inputs.tar.xz)

官方 1.2.7、当轮父版和修复版在 3000/12000 档均为 Top-1、Top-5、oracle 0/12。事先登记的 48000 档追加对照均为 0/6。每个条件的父版与修复版首姿 RMSD、oracle RMSD 均相同。所有输出对首姿的 RMSD 上下界可以从 PDBQT 坐标独立复算，允许坐标和数值的打印舍入。父版上界中，两个预算阶段因计入 G 最大改变 0.395、0.438 Å。相对晶体的 2 Å 命中率没有改善，候选内条件排序成功率无可估计分母。没有测量亲和力。[标准档逐例](standard/structures.csv) · [高预算逐例](higher/structures.csv) · [各姿势和日志](predictions.tar.xz)

| 每臂总计 | 父版 | 修复版 |
|---|---:|---:|
| 3000/12000 档，hunt / 第二次局优 / 显式精修求值 | 68118 / 22836 / 4084 | 68118 / 22836 / 3961 |
| 3000/12000 档，壁钟秒 / CPU 秒 / 子进程最大 RSS KiB | 2.660 / 2.632 / 110140 | 2.605 / 2.578 / 110140 |
| 48000 档，hunt / 第二次局优 / 显式精修求值 | 263072 / 25182 / 2730 | 263072 / 25182 / 2730 |
| 48000 档，壁钟秒 / CPU 秒 / 子进程最大 RSS KiB | 1.951 / 1.904 / 109700 | 1.995 / 1.968 / 109700 |

每次单线程单轨迹，逐块轮换算法顺序。历史 Vina 在完整 MC 步后检查上限，实际求值另记。官方二进制仅有时间与内存记录，其编译环境不同。单次计时不能确定普遍加速。修复版在低档 39 个姿势中仍有一个闭环碳距达 11.002 Å；高档 21 个姿势均小于 1.67 Å，首姿仍远离晶体。对高预算每个条件的池内最佳构象，允许仅作刚体对齐，6 个条件有 5 个内部 RMSD 小于 2 Å，质心却分别偏离 5.17–14.24 Å。余下的 4EGK/458921 即使知道最佳刚体位姿，内部形状残差仍有 2.447 Å。刚体对齐使用答案，仅用于将当前误差拆解。数据支持继续区分整体放置与内部构象失败，尚不能区分正确刚体位置为何未被生成或保留。[每姿势形状分解](higher/shape.csv) · [闭环、碰撞和实验 RMSD](higher/models.csv)

在高预算结果上，另冻结这 6 个池内最佳构象的内部坐标，分别仅平移质心、或用允许的图映射作最小二乘刚体配准，再用原拓扑和共同零参考原位重算。5 个刚体配准后进入 2 Å 的构象，Vina 总分全部变差，恶化 4.848–11.743 kcal/mol；其中 3 个出现小于 2 Å 的蛋白重原子接触，另外两个最近接触是 2.068 和 2.165 Å。刚体不变的分子内项在 PDBQT 坐标舍入下最多偏移约 0.003 kcal/mol。剩下的 4EGK/458921 的最优配准仍为 2.447 Å，总分虽从 +23.891 到 +14.630，最终分子间项仍为 +18.672 kcal/mol。因此，给出正确的整体朝向仍不能保证这些内部构象在原有势上获得有利的局部接触。此处的固定坐标对照没有局部驰豫。下节进一步检查这些起点经显式受体局优后的相对分数。[配对坐标与原始评分日志](placement.tar.xz) · [逐构象能量和碰撞](higher/placement.csv)

### 从同一构象出发，近晶体位置能否经局优保留

对上面的六个配对各进行两次原生显式 `--local_only`。一路从已输出的真实候选开始，另一路从保持其内部构象、用参考答案作刚体配准的位置开始。12 个起点的 PDBQT 重解析均逐项复现原配体的 BRANCH 声明、分型、电荷、有效转轴、共价边、闭环对和有序非键评分表。动作仅重新定义当前坐标为局部搜索起点。每对使用同一受体、盒、原生 BFGS 及相同步数上限；配体 2J9M 上限 17 步，另两个为 18 步。回溯次数随起点而变，实际求值见表。局优后的总分都用原准备配体的拓扑、共同零参考和打印的最终坐标独立重算。RMSD 保留受体坐标系，并按配体化学图对称性取最小值。

| 结构与种子 | 原始起点局优后 RMSD / 总分 | 配准起点局优后 RMSD / 总分 | 局优求值次数（原始 / 配准） |
|---|---:|---:|---:|
| 2J9M，191009 | 7.280 / −4.340 | 1.127 / +5.228 | 144 / 87 |
| 2J9M，458921 | 11.568 / −3.290 | 0.887 / −2.461 | 145 / 66 |
| 3VHA，191009 | 10.153 / −2.700 | 1.432 / +2.557 | 181 / 85 |
| 3VHA，458921 | 14.899 / −2.551 | 1.210 / +1.047 | 156 / 82 |
| 4EGK，191009 | 11.236 / −2.479 | 1.872 / −1.002 | 161 / 120 |
| 4EGK，458921 | 3.012 / +5.302 | 2.909 / +6.907 | 42 / 37 |

RMSD 单位 Å，总分单位 kcal/mol（Vina 模型值）。5 个配准时位于 2 Å 内的构象经局优后仍位于 2 Å 内，分数相对各自同样局优的错误位置高 0.829–9.568 kcal/mol。三个近晶体结果仍有一个距离小于 2 Å 的配体重原子，另外两个无此类碰撞；最低能量差不能统一归因于单个短接触。第六个构象的原始起点经局优从 4 个重原子碰撞降至 0 个，总分从 +23.891 降至 +5.302；该降幅使配准前后的排名关系反转。近晶体放置可以形成在这次局优后仍靠近晶体的状态，候选生成和评分都有改进空间。这一次短程局优不能测定整个盆地最低能量或平衡占据。12 个计数观察版输出与生产版逐字节相同。[逐起点前后分项、碰撞、图摘要与实际求值](higher/relaxation.csv) · [坐标和原始日志](relaxation.tar.xz)

另外检查是否能把精确实验重原子坐标作为同一扭转树的一种构象。对每个允许的化学映射，至少有一对属于同一刚性片段的重原子，其间距与准备输入相差 0.060–0.087 Å，超过三位小数的舍入量。这些距离在树的刚体旋转、平移和扭转下保持不变。直接把全部实验坐标回填会改变刚性片段内部几何。刚体配准诊断始终保留查询构象自身的内部几何。[逐映射的不变量见证](higher/rigid_blocks.csv)

## 复现

使用 GCC 14.2、Boost 1.83、RDKit 2025.09.3、NumPy、SciPy。重新准备新增结构另需 Meeko 0.7.1、OpenMM 8.6.1、PDBFixer 和 Gemmi 0.7.4。以下从仓库根目录执行，BOOST 指向带有 include 和 lib 的 Boost 安装前缀。上游二进制应取官方 [1.2.7 release](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7)，SHA-256 为 f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644。

```bash
D=benchmarks/macrocycle_pose_distance
W=$(mktemp -d)
mkdir "$W/inputs"
tar -xf "$D/inputs.tar.xz" -C "$W/inputs"
git worktree add --detach "$W/parent" b991eb7
python "$W/parent/tools/build.py" --boost-prefix "$BOOST" --output-dir "$W/bin-parent"
python tools/build.py --boost-prefix "$BOOST" --output-dir "$W/bin-fixed"
python "$D/paired.py" "$W/inputs" "$W/standard" "$OFFICIAL_VINA" "$W/bin-parent/vina" "$W/bin-fixed/vina"
python "$D/paired.py" "$W/inputs" "$W/higher" "$OFFICIAL_VINA" "$W/bin-parent/vina" "$W/bin-fixed/vina" --higher-budget
python "$D/observe.py" "$W/parent" "$W/bin-parent" "$BOOST" "$W/trace-parent"
python "$D/observe.py" . "$W/bin-fixed" "$BOOST" "$W/trace-fixed"
for stage in standard higher; do
    python "$D/replay.py" "$W/inputs" "$W/$stage" "$W/events-$stage" "$W/trace-parent/vina" "$W/trace-fixed/vina"
    python "$D/reported.py" "$W/$stage"
    python "$D/shape.py" "$W/inputs" "$W/$stage"
done
python "$D/place.py" "$W/inputs" "$W/higher" "$W/bin-fixed/vina" "$W/placement"
c++ -O2 -std=c++14 -pthread -Isrc/lib -I"$BOOST/include" "$D/topology.cpp" \
    $(find "$W/bin-fixed/obj" -name '*.o' ! -name main.o) -L"$BOOST/lib" -Wl,-rpath,"$BOOST/lib" \
    -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -o "$W/topology"
python "$D/local_observer.py" . "$W/bin-fixed" "$BOOST" "$W/local-observer"
python "$D/relax.py" "$W/inputs" "$W/higher" "$W/placement" "$W/bin-fixed/vina" "$W/topology" "$W/local" --observer "$W/local-observer/vina"
python "$D/rigid_blocks.py" "$W/inputs" "$W/higher/rigid_blocks.csv"
python "$D/no_markers.py" "$W/bin-parent/vina" "$W/bin-fixed/vina" "$W/controls"
```

要重新准备两个新增结构，从空目录运行 `OPENMM_CPU_THREADS=1 python benchmarks/macrocycle_pose_distance/prepare.py /path/to/new-inputs`。既有 3VHA 可从 benchmarks/concerted_ring/new_inputs.tar.xz 提取。C++ 测试 benchmarks/check_late_pose_geometry.cpp 覆盖零个可观察位点、四种 G 标签、真实 CG 和水合 W，另包含 XS / AD 格点域检查。无闭环标签的 1IEP 和含可动侧链的固定测试在修改前后输出逐字节一致，六项计算计数逐项相同。[测试输出](parity.json)

本地用新构建的生产对象重放了 core-regression.yml 全部 11 个测试步骤，11/11 通过（另外完成构建，未在本地执行托管制品上传）。将相同的新增几何测试链接到 b991eb7 的生产对象，原生坐标数断言如期失败；链接修复对象后通过。检查记录列出各个测试与退出状态。[本地回归记录](local_checks.json)

修改原生姿势距离的远程提交已通过 [GitHub Actions 核心回归](https://github.com/HankerWu/fastdock-vina/actions/runs/45962205527)。
