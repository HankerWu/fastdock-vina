# AD4 未参数化原子与静电图使用了同一索引

父提交 `400d57d` 允许配体或柔性受体输入 Cu、Na、K、Hg、Co、U、Cd、Ni。Vina 和 Vinardo 为这些元素提供通用 XS 金属类型。AD4 的探针类型表没有它们的条目。缺失标志 `AD_TYPE_SIZE` 恰是通用静电图的下标。可动位点通过了图完整性检查，固定柔性锚点则可能被跳过显式配对。刚性蛋白预先制成的 AutoGrid 图无需从 PDBQT 重新解释其中的固定金属。

AD4 在图内按 $A_t(x)+q\Phi(x)+|q|D(x)$ 组合有类型的原子场、通用静电场和电荷相关的去溶剂化场，随后应用已有正能量软上限。旧程序把缺失的 $A_t$ 换成了 $\Phi$，使电荷设为零时仍有一个完整静电场。该非法替代缺少排斥、色散等有类型的贡献。仅添加一份自命名图也无法定义引擎所需的配体内和柔性受体成对参数。修正从解析时向配体、支链和柔性残基传递所选类型体系。AD4 遇到上述八种类型即返回带原始行的错误。Fe 有正式 AD4 参数，正常走 Fe 通道。在多配体 C++ 接口中，输入列表的所有配体成功解析后才替换上次完成的组合，避免失败后旧查表和新模型混用。[AD4 原始方法](https://doi.org/10.1002/jcc.21256) · [AutoDock Vina 的 AD4 建图与对接方法](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html) · [解析实现](../../../src/lib/parse_pdbqt.cpp)

## 两个真实位点的定点检验

采用人葡萄糖激酶的两份已归档晶体，按 PDB 沉积的元素和 A 链 600 号坐标构造单离子探针。两份记录均标注完整占据率，并记载该离子与 Met238、Val241、Val244、Gly246 的主链氧接触。干蛋白来自原有制备，未包含待测离子。使用 AutoGrid 4.2.6 在离子位置中心重新生成 16 × 16 × 16 体素图，间距 0.375 Å。原子坐标恰在中心格点，可以直接从文件读取 $\Phi$ 和 $D$。这里检验同一人源蛋白、同一对应位点的两种输入。沉积指认无法单独确定生理金属特异性或离子的结合自由能。[1V4S 原始坐标](https://www.rcsb.org/structure/1V4S) · [3IDH 原始坐标](https://www.rcsb.org/structure/3IDH) · [来源和受体制备](../../ligand_hydroxyl/glucose/README.md)

| 沉积位点 | 静电通道 $\Phi$ | 去溶剂化通道 $D$ | 旧 AD4 的 +1 探针分数 | 人工置零电荷的旧分数 | 修正后 AD4 |
|---|---:|---:|---:|---:|---|
| 1V4S Na，A600 | −2.306 | 1.131 | −3.481 | −2.306 | 拒绝该类型 |
| 3IDH K，A600 | −2.201 | 1.030 | −3.372 | −2.201 | 拒绝该类型 |

两处 +1 探针分数逐项复现 $2\Phi+D$。把同位置人工替换成有专用图的 C、Fe、Zn 后，六个父版/修正版分数逐一相同，且均可独立从格点值复算。原 Na/K 在 XS 的 Vina 与 Vinardo 上共四项分数相同。局优和从头对接的父版均成功返回，新版均在打分前拒绝非法 AD4 类型。该结果消除了无参数输入上的伪分数。没有为 Na/K 建立配位势，也没有产生可用于比较预测质量的新离子姿势。[全部定点评分](scores.csv) · [入口状态](entry_points.csv) · [原始图、建图配置、探针、日志和坐标](raw.tar.xz)

对已参数化有机配体另行运行两版生产程序，复算既有 FKBP12 同测定系列八个化合物、两个受体、三枚种子的全部 48 个冻结首姿。逐条分数和测定统计表与既有归档逐字节相同。按每个受体内种子中位数、再取受体均值，三倍且区间分离的分子对仍为 6/15，七个点值的 Spearman 为 −0.179，Kendall 为 −0.143。本修正发生在读取输入时，对有效配体没有新增对势求值；本项未作壁钟速度比较。[逐个首姿](../../mechanism_interactions/macrocycle_selection/fixed_affinity.csv) · [测定和分组](../../cubic_grid/fkbp_ad4/README.md)

## 复现

在仓库根目录运行。第一个参数是从 `400d57d` 构建的生产程序，第二个参数是修正版。原料从本分支既有存档读取。C++ 回归另覆盖八个缺参金属、Fe、柔性残基可动原子及固定锚点、多配体替换失败后的旧模型，父版在第一个 Cu 断言失败。

```bash
D=benchmarks/ad4_atom_field/metal_types
W=$(mktemp -d)
python "$D/replay.py" /path/to/parent/vina /path/to/fixed/vina /path/to/autogrid4 "$W/sites"
cmp "$W/sites/scores.csv" "$D/scores.csv"
cmp "$W/sites/entry_points.csv" "$D/entry_points.csv"
python benchmarks/macrocycle_domain/fixed_affinity.py "$W/fkbp" /path/to/parent/vina /path/to/fixed/vina
for t in fixed_affinity fixed_affinity_stats; do
  cmp "$W/fkbp/$t.csv" "benchmarks/mechanism_interactions/macrocycle_selection/$t.csv"
done
```
