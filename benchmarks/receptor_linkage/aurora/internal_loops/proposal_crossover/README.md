# 把每个真实提议交给两张配体图做同一轮 hunt 局优

AT9283 种子 61231 在漏边图上曾发布 1.892 Å 的首姿。对这次搜索的全部 1,474 个原始提议，完整 SDF 图能从同样的五个入口局优到 2 Å 内。反向实验中，完整图盲搜的 1,638 个原始提议交给任一图，局优均停在 7 Å 以外。此前的固定池重排检验了已经发布的产物，本实验向前追到每一步变异与首次格点局优之间。[判据与选择依据](PROTOCOL.md) · [双图逐提议坐标、求值次数和 RMSD](results/paired.csv)

| 提议来自哪条真实搜索 | 执行 hunt 的图 | 提议数 | 返回 RMSD 最小值，Å | ≤2 Å | ≤3 Å | hunt 实际求值 |
|---|---|---:|---:|---:|---:|---:|
| 历史推断图 | 推断图 | 1,474 | 1.752840 | 5 | 11 | 48,331 |
| 历史推断图 | 完整图 | 1,474 | 1.752838 | 5 | 11 | 48,595 |
| 完整图 | 推断图 | 1,638 | 7.248777 | 0 | 0 | 46,905 |
| 完整图 | 完整图 | 1,638 | 7.238488 | 0 | 0 | 47,652 |

五个近晶体返回均来自历史漏边图的第二条轨迹第 706–710 步，彼此连续、相关。第 706 步的原始提议为 2.494576 Å，两图均使用 34 次目标与梯度联合求值，首次将它局优到 2 Å 内。其余四个提议原本已在 2 Å 内。两图的十一条 ≤3 Å 返回也逐条对应。其他提议的局优会受化学图影响，例如同一提议的最终 RMSD 可相差超过 3 Å。证据将本例局限在进入这些提议前的搜索历史，尚未确定早先哪次接受、局优、位移或扭转对抵达这条轨迹充分或必要。[原始成功邻域和对照](results/inferred_on_complete.npz) · [完整图在自身历史上的响应](results/complete_on_complete.npz) · [独立重读四格](../mutation_controls/results/crossover.json)

## 源码、输入与计量

研究副本在 mutate_conf 后、quasi_newton 前按 17 位有效数字导出每个构象的三维平移、四元数和五个扭角。接收图使用原生产解析器、同一刚性受体、格点、hunt_cap、局优器和迭代数。交叉检验在原有的预计算对象和格点上，用独立模型副本执行；主搜索继续运行以复核仪器对照。用于交叉局优的输入只含真实提议与轨迹标识。用晶体计算对称性校正 RMSD 的代码在所有进程完成后运行。两条对角线的 3,112 个返回的每个笛卡尔坐标、原能量与真实 MC 轨迹逐比特相同。六次主搜索发布池逐字节重现。另一独立实现从实际 PDB 与不同图映射重算了 6,224 次交叉的起止误差，最大差约 3.3×10⁻¹⁴ Å。[调用点](calls.patch) · [导出与交叉局优实现](cross_hunt.h) · [运行与评价](run.py) · [独立复核代码](../mutation_controls/verify_crossover.py)

坐标均以受体为参考系，距离单位 Å，扭角为弧度，四元数无量纲。硬蛋白是 2X6E 加源自 1OL7 的撤约束静态环，281–291 与 303–305 两段在位。软件使用既有中性配体、原子分型及 Vina 的干态经验势。水、侧链、离子及化学态之间的交换未被采样。修改成键图会改变内部分数及相互作用排除，格点、姿势树、微态与身体自由度保持一致。这里计算的是局部优化算子对指定输入的响应。将这一串相邻输入做 Boltzmann 加权会把路径重访误作相空间体积。平衡盆地权重需要相同微观势和化学势下的结合/游离配分函数，参见[系综推导](../../../../../docs/BINDING_ENSEMBLE_AND_DOCKING.md)和 [Vina 原始方法](https://doi.org/10.1002/jcc.21334)。

六个进程累计墙钟 2.939 s，CPU 2.934 s，包含六次正常盲搜的 301,722 次求值，外加交叉 hunt 的 191,483 次。交叉局优从第一条轨迹开始前额外运行，只作诊断。二者相加为 493,205 次。单文件增量编译和链接另用墙钟 4.171 s、CPU 4.139 s。此前受体建模、生成完整轨迹、离线 RMSD、编解压的成本分开。每个 npz 保留双精度起止重原子坐标、PDBQT 序号、全体输入广义构象、原能量及求值数；数组行与 CSV 以来源、轨迹和步对应。[逐进程计时](results/costs.csv) · [常规搜索日志](logs.tar.xz) · [完整阶段回放](../basin_trace/README.md)

## 复现

先按上级 REPRODUCE.md 解开 WORK 内的原始输入；按 basin_trace 的说明在 TRACE 中生成全部回放（含 .bin、.pdbqt）。以下命令从仓库根目录执行，WORK、TRACE 与 BOOST_PREFIX 均为绝对路径。BUILD 可设为新的空目录。从科学结果提交或其后包含相同 src 的提交执行，研究副本使用父版生产源码。

```sh
R=$(pwd)
D="$R/benchmarks/receptor_linkage/aurora/internal_loops/proposal_crossover"
mkdir -p "$BUILD/source"
cp -R src tools "$BUILD/source/"
cp "$D/cross_hunt.h" "$BUILD/source/src/lib/"
(cd "$BUILD/source" && git apply "$D/calls.patch")
python "$BUILD/source/tools/build.py" --output-dir "$BUILD/bin" --boost-prefix "$BOOST_PREFIX" --jobs 2
python "$D/run.py" "$WORK" "$TRACE" "$BUILD/bin/vina" "$BUILD/replay"
cmp "$D/results/paired.csv" "$BUILD/replay/paired.csv"
```

全量干净重编耗时 35.240 s 墙钟、57.978 s CPU，重放并离线生成四份坐标和误差表耗时 4.390 s 墙钟、4.370 s CPU。四份坐标档案的所有数组均与首次计算逐元素精确相同，逐提议表字节相同。[重编复现计时](results/rebuild.json)

独立复核可直接使用已归档的四份 npz，也可指向重放目录：

```sh
python "$D/../mutation_controls/verify_crossover.py" "$WORK" "$D/results" "$TRACE" > "$BUILD/crossover.json"
cmp "$BUILD/crossover.json" "$D/../mutation_controls/results/crossover.json"
```
