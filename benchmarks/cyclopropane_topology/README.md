# 三元碳环的配体拓扑与真实姿势

在盲生成的四份配体输入上，C++ 成键推断会删除三元环最长边。用 CCD 化学图证实后补回这条边，使原来误列为普通非键相互作用的 9 个配体内原子对退出评分。相同名义预算下，2 Å 首姿、前五、输出池命中各由 1/16 升至 3/16。新增的两次命中都在此前已知漏键的 PDE4B–罗氟司特；预先指定的两个迁移靶点均未获得 2 Å 候选。可用于可靠输入的做法是给出完整化学图。只依赖距离的 C++ 实现留作可运行对照。

## 物理机制与因果对照

assign_bonds() 的几何遮挡判据移除端点之间存在第三个更近原子的候选键。一个三角环若有唯一最长边，第三个顶点与另两端的距离必定较短。删边改变了键距，真实 1–4 邻居会落入非键表。原生 Vina 对可变距离的常规配体内相互作用使用键距大于三的排除规则。误列的原子对同时进入内部势和导数，经局优、接受、保存和显式精修影响候选坐标。

实验例外只处理同一配体、三种相对距离均固定、三个 AD 类型全为普通脂肪碳 C 的三角形。三边均须落在 [d₀/1.1, 1.1d₀]，d₀ 为两个碳的共价半径和，1.1 沿用原代码的键长容差。输入坐标以 Å 表示。这个条件是一种几何推断，未拟合权重。另一路直接以完整化学图生成 --ligand_bonds，在同一父版程序上运行。四体系共 16 个分子–种子–预算格内，这条显式路径与实验 C++ 分支所得 PDBQT 逐字节一致，求值计数逐项一致。阴性图对照 6WTO、4LMN 均无边或非键表变化。[基线源码与归档补丁](cyclopropane.patch) · [逐配体图比对](results/graph.csv) · [Vina 原始评分和搜索](https://doi.org/10.1002/jcc.21334)

在 72 个父版输出姿势的相同坐标上，实验实现对每个姿势的配体–受体分项及自参考总分与父版相同。配体内能量差为 −0.721 至 +0.036 kcal/mol。用共同零参考的分子间、分子内能量和排序，16 个池有 1 个首位改变，此池最优 RMSD 仍大于 7 Å。两个新增正确姿势均出现在原父版池未命中时。这能定位候选生成或筛选环节的贡献，现有轨迹无法将贡献进一步归到某一条 MC 提议。[逐坐标分项](results/fixed.csv)

## 姿势和成本

每一行包含两个种子、两个预算。表中各 RMSD 为四格的中位数，分辨率与配体身份来自实验结构。每个配体的查询坐标由 CCD 化学图、固定种子 ETKDGv3 和 MMFF 优化生成；所有方法共享共晶配体定义的 24 Å 口袋中心。参考图允许的立体等价原子映射仅用于评估。

| 体系 | 分辨率 Å | 错误非键对 | 首姿成功格数，父版 → 实验 | 首姿 RMSD Å，父版 → 实验 | Oracle RMSD Å，父版 → 实验 |
|---|---:|---:|---:|---:|---:|
| PDE4B–罗氟司特，1XMU | 2.30 | 1 | 1→3 / 4 | 8.95→0.70 | 7.90→0.70 |
| VEGFR2–仑伐替尼，3WZD | 1.57 | 1 | 0→0 / 4 | 9.96→6.32 | 7.85→5.24 |
| PARP1–奥拉帕利，7AAD | 2.21 | 2 | 0→0 / 4 | 5.04→5.50 | 4.40→5.35 |
| OX2R–莱博雷生，7XRR | 2.89 | 5 | 0→0 / 4 | 6.86→5.99 | 6.00→5.99 |

原生 Vina 1.2.7 同协议也是 1/16。父版、C++ 实验和显式图三臂合计联合求值分别为 124493、123828、123828。16 个新进程墙钟之和分别为 3.46、3.81、3.45 秒，上游为 20.38 秒。每格只计时一次，速度差异可能有负载干扰。PDE4B 在种子 817139、3000 档的首姿由 8.616 变为 0.584 Å；种子 191009、12000 档由 9.284 变为 0.823 Å。PARP1 在 817139、12000 档则由 4.557 退至 8.451 Å。所有组的受体均为刚性。1XMU 保留 Zn 和 Mg；3WZD 未保留最近 2.600 Å 的甘油，7AAD 未保留最近 6.677 Å 的硫酸根。四份输入全新复算逐字节一致；本地缓存原始 PDB/CCD 的准备耗时约 3.25–4.20 秒/复合物。未开展分子动力学，未测试亲和力排序。[所有运行](results/runs.csv) · [全部姿势](results/modes.csv) · [受体制备限制](RECEPTORS.md) · [PARP1 结构](https://www.rcsb.org/structure/7AAD) · [OX2R 结构](https://www.rcsb.org/structure/7XRR)

## 几何反例与集成判断

用折叠成近等边三角形的开链丙烷检验歧义。三个碳均落入上述共价长度范围，输入声明的键只有 1–2、2–3。父版与化学图一致，实验实现额外构建 1–3。它是人为畸变输入，可直接复现歧义；它不用于估计真实配体的假阳性率。四个实测复合物支持可信配体图对部分姿势搜索有价值，无法据此自动推广几何例外。普通受体对接已有 --ligand_bonds 可接受完整准备图。另一项默认修复处理 PDBQT 已经给出的转轴边。[反例](ambiguous.pdbqt) · [声明轴的独立研究](../declared_rotor_bonds/README.md)

## 复现

从仓库根目录运行。Python 依赖为 RDKit 2025.09.3、Meeko 0.7.1、gemmi、NumPy、SciPy；重新制备需要 PDBFixer、OpenMM。实际二进制散列见 binaries.json。可令 BOOST_PREFIX 指向系统外的 Boost；构建 Vina 时相应传入 --boost-prefix。

~~~bash
D=benchmarks/cyclopropane_topology
W=$(mktemp -d)
git worktree add --detach "$W/parent" dec913f
git worktree add --detach "$W/triangle" dec913f
git -C "$W/triangle" apply "$PWD/$D/cyclopropane.patch"
for arm in parent triangle; do
  python "$W/$arm/tools/build.py" --output-dir "$W/b-$arm"
  bash benchmarks/declared_rotor_bonds/build_probe.sh "$W/$arm" "$W/b-$arm" \
    "$PWD/benchmarks/bond_distance/ligands/graph.cpp" "$W/graph-$arm"
done
mkdir -p "$W/run"
tar -xJf "$D/inputs.tar.xz" -C "$W/run"
python "$D/panel.py" "$W/run" "$UPSTREAM_VINA_1_2_7" "$W/b-parent/vina" \
  "$W/b-triangle/vina" "$W/graph-parent" "$W/graph-triangle"
"$W/graph-parent" "$D/ambiguous.pdbqt"
"$W/graph-triangle" "$D/ambiguous.pdbqt"
# 可选，从原始 PDB 和 CCD 重新制备（与归档四个受体、四个配体逐一比对）
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python "$D/prepare_inputs.py" "$W/fresh"
~~~
