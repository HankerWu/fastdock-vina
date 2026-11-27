# AD4 格点中的电荷、限幅与有效域

人 FKBP12（1FKG）的一个真实搜索起点完全在盒内。冻结全部坐标，旧 AD4 图能为 284.372 kcal/mol，只修几何壁仍为 284.372，按电荷组装原子势后为 280.930，最大笛卡尔梯度分量改变 0.376 kcal/(mol·Å)。同一配体另一个种子的首个构象部分越界，单改几何壁就减少 1,649,307 kcal/mol 的人工罚项。[八个共同起点的能量、梯度和逐原子同坐标核验](first_trials.csv)

## 算子与调用链

AutoGrid 类型图 A 已含短程相互作用和电荷无关的去溶剂化贡献。给定原子部分电荷 q（元电荷）、单位电荷静电图 E、单位绝对电荷去溶剂化图 D，原子受体势为 U=A+qE+|q|D，单位 kcal/mol。FastDock 延用搜索器的正能量限幅 C_v(u)=vu/(v+u) 当 u>0；u≤0 时保留 u，v=0 时正值归零。正支导数为 [v/(v+u)]²。hunt 使用 v=10，最终 AD4 打分使用 v=1000。限幅是一种优化近似。[AD4 原始力场](https://doi.org/10.1002/jcc.20634) · [Vina 1.2 的 AD4 方法](https://doi.org/10.1021/acs.jcim.1c00203)

旧 ad4cache 在三个通道上分别限幅，再施加电荷因子，计算的是 C(A)+qC(E)+|q|C(D)。它在知道配体电荷符号之前就抑制了正电势，同一单原子势在通道间重新分配后还会改变得分。三次格点求值各附带几何壁，旧能量在盒外还获得 s(1+q+|q|) 的壁斜率。上轮已改正固定打分的原子有效域，并用 BACE1 的 HD 边界坐标量出了这一系数，本轮改动发生在实际参与 BFGS、MC 及最终评分的能量与梯度。

令 y 是 x 在盒内的最近点。现在每原子贡献 C_v(U(y))+s‖x−y‖₁。三线性插值是线性的，先合并同一格元八个顶点的带权通道再插值一次，等价于先对各通道插值再线性组合。随后施加一次能量限幅与一次独立几何壁。图内梯度为 C′(U)∇U；越界坐标的法向导数为 ±s，切向导数保留。原子的化学类型、AD4 参数和搜索预算均保持原定义。配体、可动侧链与全部可动原子的分项及求导调用同一原子算子。回归涵盖正、负、零电荷，v=0、10、1000、无限幅，通道抵消，图内外三维中心差分，配体与柔性侧链分解。[数值及公共 API 回归](check_field.cpp)

## 真实姿势、反事实与计算成本

每格一条轨迹、一核、五个名额、100 kcal/mol 能窗，交错运行。FKBP12 开发集含两个受体与两个配体的四种 cognate/cross 组合，两种种子和两档预算。迁移集在运行前审查，人 CDK2–U32（2BTS，1.99 Å）与 PPARγ–BRL（2PRG，2.30 Å）通过。3MXF 的省略 DMS 距 JQ1 5.36 Å，4HG7 的硫酸根距 Nutlin 3.54 Å，未通过既有的 8 Å 准备门槛。CDK2 建图报告一次防止除零的警告，随后成功结束；各通道均为有限数，两次独立建图哈希相同。初始配体由化学图嵌入，晶体配体用于口袋中心与评价。cross-docking 用匹配的口袋 Cα 对齐受体。RMSD 为晶体配体等价重原子映射的最小值。[本地事先记录的方案与后续冻结说明](PROTOCOL.md) · [2BTS](https://www.rcsb.org/structure/2BTS) · [2PRG](https://www.rcsb.org/structure/2PRG)

表内为 Top-1 / Top-5 命中数，阈值为 2 Å。五个输出名额使 Top-5 等于输出池 oracle。种子、预算、受体重复均属于相关观测。

| 配置 | 上游 AD4 | 上游 Vina | 父版 Pro/AD4 | 仅修壁罚 | 合并原子势 |
|---|---:|---:|---:|---:|---:|
| FKBP12 开发，16 格 | 0/0 | 0/0 | 0/0 | 1/1 | 1/2 |
| CDK2、PPARγ 迁移，8 格 | 0/0 | 0/0 | 0/0 | 0/0 | 1/1 |
| 同两个迁移受体，新种子，12 格 | 2/2 | 2/2 | 2/2 | 1/3 | 2/2 |

CDK2 在 12,000 cap、种子 14591 下，合并原子势得到 0.521 Å 首姿。两个确认种子在相同预算下均未重现；增加至 48,000 后，父版与合并势都两次命中 CDK2。PPARγ 的十个配置中合并势均未产生 2 Å 候选。把父版、合并势各 36 个输出池交叉打分，72 个池中两种公式在相同坐标上选中的首名完全一致。观察到的首姿差异来自搜索期间势面改变。仅修壁罚在确认组的三个有正确候选的池中只排对一个。[逐预算统计及条件排名](summary.csv) · [固定池排名](fixed_pool_ranks.csv)

开发与迁移合计 24 格，父版 / 合并势的实际求值为 181,087 / 181,632，墙钟合计 11.453 / 10.852 秒，CPU 合计 11.346 / 10.769 秒。名义上限按整步检查，实际次数会略有出入。计时包含每次读图，未证明稳定加速比。两个新受体制备各约 3.47、3.14 秒，可复用的 AutoGrid 建图各约 6.86、7.11 秒。CSV 中 RSS 是子进程高水位，包含启动成本。较高预算对照尚未证明算法提高姿势预测上限。

从 411 个输出 PDBQT 直接以 SciPy 重建 AD4 能量。378 个误差低于 0.01 kcal/mol。其余 33 个用打印坐标 ±0.0005 Å 区间的格点顶点、单调限幅与壁罚界复算，报告分数均落入所得范围。壁斜率达 10⁶，靠近盒面的三位小数坐标会严重放大绝对能量的复核误差。[独立逐姿势审计](fixed_pool.csv) · [复算实现](audit.py)

## 输出有效域

411 个实验 AD4 姿势中，22 个坐标越界，含合并势的 5 个。固定打分接口会拒绝相同坐标，原先全局输出却会把超过 600 万 kcal/mol 的壁罚当作对接结果。现在每条轨迹搜索结束后、跨轨迹 RMSD 和容量合并之前，在各自尚存的池中剔除越界项，使剩余合法候选竞争名额。适用于 AD4、仅格点 Vina 与 Vinardo。通过与固定打分相同的、按力场选择原子的判据检查可动配体和侧链，原有每轴容差为 0.0001 Å。显式受体精修仍有后续机会将搜到的构象带回盒内，沿用其原有处理。

完整搜索后没有合法候选时，API 返回空结果，CLI 输出零模型，清除同路径的陈旧结果。尚未搜索、或刚更换输入的 API 实例对已有文件保留原有行为。重放全部 36 格，合并势输出 103 → 98，五格空结果计为失败；其余 98 个坐标及分数逐字节相同，36 格求值数相同。另两个三轨迹混合池均保留全部三个合法输出。该筛选保护现存候选，无法追回单条轨迹已丢弃的构象。[逐格重放](release.csv) · [多轨迹](multiple.csv) · [三种地图的可行和不可行输入](check_release_cli.py)

## 同一实验的 FKBP12 亲和力

Holt 等在 10 °C、pH 7.8 的同一偶联 PPIase 实验测得八个同系列酯类抑制剂的表观 Ki。化合物 10 保留 300–600 nM 区间。每个受体先取三个种子的中位首姿分数，再等权平均两个受体。15 个唯一分子对的区间边界仍相隔至少三倍，相关系数使用七个点测定化合物。冻结 48 个旧首姿，两种公式均判对 6/15。重新对接每格取 3 条轨迹、24,000 cap、9 个输出名额，两个种子组各 96 个搜索。[Holt 实验](https://doi.org/10.1021/ja00075a008) · [制备与实验值](../cubic_grid/fkbp_ad4/README.md)

| 受体 | 原种子组，父版 → 修正 | 新种子组，父版 → 修正 |
|---|---:|---:|
| 1FKG | 6 → 4 / 15 | 9 → 7 / 15 |
| 1FKH | 9 → 11 / 15 | 4 → 9 / 15 |
| 两受体平均 | 6 → 9 / 15 | 6 → 9 / 15 |
| 平均后 Spearman ρ | −0.036 → 0.286 | −0.357 → 0.071 |
| 平均后 Kendall τ | 0.048 → 0.143 | −0.333 → −0.048 |

两个种子组都改善三个合并分子对，且各自没有丢失此前判对的合并分子对。其中只有 9 强于 7 在两组中都获救。1FKG 的退化也在两组中重现。实际求值，原种子组为 3,460,051 → 3,460,095，新种子组为 3,459,930 → 3,460,384。可以确认的是本系列的多受体汇总结果，绝对排序仍弱。该实验没有测得结合自由能或受体平衡权重。下一步可用未参与本轮的新同靶系列检验构象选择和分子排名能否共同复现。[原种子逐对结果](affinity_primary/affinity_pairs.csv) · [新种子逐对结果](affinity_confirmation/affinity_pairs.csv)

## 复现

从本仓库根目录运行。需要 C++14、Boost、Python、NumPy、SciPy、RDKit、Meeko、gemmi。重新制备受体另需 OpenMM 与 PDBFixer。AutoGrid 为 Debian autogrid 4.2.6-10 (amd64)，软件包 SHA256 是 751e5256ce2d8452abdddff83abb2ce92d43e629ccd9cf84daa0e354f24149b7。包内二进制 SHA256 是 636df4a0075adeb1a04e9c26b0674ad8b77d49e0244f7a109aa302214a670846。独立重建时两个受体的每个地图都通过 SHA256 核验。上游 AutoDock Vina 基线采用提交 3c65c0b。[AutoGrid 包](https://deb.debian.org/debian/pool/main/a/autodocksuite/) · [上游源码](https://github.com/ccsb-scripps/AutoDock-Vina/tree/3c65c0b)

~~~bash
D=benchmarks/ad4_atom_field
W=$(mktemp -d)
mkdir -p "$W/input"
tar -xf benchmarks/cubic_grid/human_panel/prepared.tar.xz -C "$W/input"
tar -xf benchmarks/cubic_grid/human_panel/affinity-pools.tar.xz -C "$W/input" --wildcards --no-wildcards-match-slash 'affinity/*.pdbqt'
for r in 1FKG 1FKH; do tar -xf "benchmarks/cubic_grid/fkbp_ad4/maps-$r.tar.xz" -C "$W/input"; done
tar -xf "$D/prepared_transfer.tar.xz" -C "$W/input"
python "$D/prepare.py" "$W/input" "$(command -v autogrid4)" --maps-only

git worktree add --detach "$W/parent" f0e6dbf
git worktree add --detach "$W/wall" f0e6dbf
git -C "$W/wall" apply "$PWD/$D/wall_only.patch"
git worktree add --detach "$W/joint" HEAD
git -C "$W/joint" restore --source=f0e6dbf src/lib/parallel_mc.cpp src/lib/parallel_mc.h src/lib/vina.cpp src/lib/vina.h
for arm in parent wall joint; do python "$W/$arm/tools/build.py" --output-dir "$W/b-$arm"; done
python tools/build.py --output-dir "$W/b-release"
python benchmarks/pair_gradient/build_arms.py "$W/up" --arms upstream
p="$W/b-parent/vina"; w="$W/b-wall/vina"; j="$W/b-joint/vina"; r="$W/b-release/vina"; u="$W/up/upstream/build/vina"
for split in development transfer vina_reference confirmation; do
    python "$D/run_panel.py" "$W/input" "$W/panel" "$p" "$w" "$j" "$u" --split "$split"
done
for arm in parent wall joint; do
    python "$D/first_eval.py" "$W/$arm" "$W/b-$arm/obj" "$W/input" "$W/first" --arm "$arm"
done
python "$D/summarize.py" "$W/panel" "$W/first" "$W/panel"
python "$D/audit.py" "$W/input" "$W/panel" "$W/panel"
python "$D/release.py" "$W/input" "$W/panel" "$W/released" "$j" "$r"
python "$D/check_release_cli.py" "$r"
python benchmarks/macrocycle_domain/fixed_affinity.py "$W/fixed" "$p" "$r"
python "$D/affinity_search.py" "$W/input" "$W/assay-primary" "$p" "$r"
python "$D/affinity_search.py" "$W/input" "$W/assay-confirm" "$p" "$r" --confirm
python "$D/verify_final.py" "$W/input" "$r" "$W/verified.json"
~~~

binaries.json 保存结构实验的可执行文件指纹。structural_outputs.tar.xz 保存原始对接坐标和日志，first_trials.tar.xz 保存首次求值的全部可动原子坐标和梯度，released.tar.xz 单独保存有效域修复后的输出。两个 affinity 文件夹分别保存各自种子组的完整坐标、日志、实际求值及时间。重新计时会受宿主负载影响。
