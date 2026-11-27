# 同一大环姿势为何在搜索和最终打分时名次不同

环切开后，两个辅助位点 G 分别应与对侧闭环碳 CG 重合。闭环是已登记的共价拓扑约束。一个端点飞离 20 Å 不会使分子成为化学上合法的开链配体。分支父版 `3f21e02` 的线性约束在该距离归零，插值表末端还会混入零样本的导数。这个问题此前已在雷帕霉素上发现并以 XS 研究补丁处理。此次验证包含另一蛋白家族、另一配体和另一套评分函数，并把同一约束拓展到生产查表和解析尾部。[先前的独立研究与有限差分](../../bond_distance/README.md)

## 两种排序的反事实

对人 HSP90α（3VHA）、组织蛋白酶 D（6QBH）和 CDK2（8H6P）使用既有独立制备的配体和同一 26 Å 盒。主诊断为 Vina/Vinardo 各两枚种子、每种子 4 条轨迹，每轨迹名义上限 12,000，至多 20 个输出。栅格路径保留已释放的构象并比较两种列表顺序。一个使用引擎输出次序，另一个按三位小数报告的最终分数稳定排序。另运行相同搜索参数的显式精修、历史单轨迹小预算、无闭环配体 1IEP。面板包含已知案例，用于找出排序冲突的来源。小预算的 6QBH 一枚种子在两档均无合法盒内候选，表中保留空值。初步定位后才扩大到四轨迹、20 模式并加显式精修对照。这是机制诊断。

在父版四轨迹栅格搜索中，12 池共 188 个姿势。11 池出现显示分数逆序，235 对逆序中有 188 对至少含一个闭环碳距大于 3 Å 的构象（沿用历史严重断键定义）。从日志取保存目标，减去最终打印的未归一化 INTER + INTRA，然后独立按已四舍五入的坐标重建两个闭环项。扣除闭环项后，绝对残差中位数 0.068，最大 0.147；在 51 个严重断键上，中位数从约 611.35 变为 0.0086 个目标函数单位。残差包括分子内能插值/中点差和坐标、能量打印精度。12 个池的 60 个前五名中，直接按显示分数重排会使严重断键从 0 个增加到 6 个。逐池平均 Best@5 从 5.990 升至 6.093 Å，Top-1 均值仍为 6.531 Å。6QBH 的一枚种子中，重排会把真实共价键比晶体长 6.47 和 11.93 Å 的姿势排进前五。显式精修组按输出能量排序，0 对逆序，前五仍有两个严重断键。程序沿用原有经验分数、排序和未结合态约定。

## 限定修复和验证

令两个相应辅助位点的距离为 $r$，线性约束为 $U_c(r)=w r$。若原有成对上限 $M$ 有效，先取 $\min(M,w r)$。对正能量还沿用原局部搜索软化 $C_v(u)=vu/(v+u)$。在线性未触顶区，径向梯度为 $w v^2/(v+wr)^2$。内表保持现有采样和插值，表末用于插值/求导的额外节点也保留线性约束。超过最大表距时仅对已登记的 CG/G 对解析求值。AD4 普通静电和溶剂化作用仍在其原截止 20.48 Å 内按表计算。未将闭环罚加入最终亲和力。当前修复消除了远距离断崖；在近端查表、正能量软化以及树状扭转变量下，端点仍可能分离。

6QBH 父版 Vinardo 栅格池的一条姿势最大标记误差为 20.460 Å、真实闭环碳距为 19.061 Å，两个辅助项因越过 20 Å 都为零，记录的搜索能减最终未归一化能为 −0.003。相同实验系列的显式 Vina 另有 21.728 Å 的真实碳距。修复后所有输出的标记误差均低于 20 Å。以下每组均含 12 个运行条件。高预算组在两种力场、三个靶点、两枚种子间等权汇总；低预算为 Vina、两枚种子和 3,000/12,000 两档。种子、预算与算法路径相关。

| 父版 → 修复 | 栅格单轨迹 | 栅格四轨迹 | 显式四轨迹 |
|---|---:|---:|---:|
| 总输出 / 无输出条件 | 30/2 → 28/2 | 188/0 → 189/0 | 210/0 → 209/0 |
| 严重断键 / 总输出 | 8/30 → 8/28 | 51/188 → 45/189 | 10/210 → 8/209 |
| 最大辅助误差 ≥20 Å 的输出 | 0 → 0 | 2 → 0 | 2 → 0 |
| Top-1 均值，Å（有输出条件） | 12.207 → 11.686 | 6.531 → 6.531 | 6.532 → 6.532 |
| Best@5 均值，Å | 10.316 → 10.316 | 5.990 → 6.042 | 5.681 → 5.733 |
| 输出池 oracle 均值，Å | 10.316 → 10.316 | 5.534 → 5.534 | 5.024 → 5.024 |
| 搜索、二次局优和显式精修求值总计 | 91,180 → 91,417 | 581,094 → 581,063 | 604,534 → 603,385 |
| 进程 CPU 秒（单轮、含建图） | 2.79 → 2.54 | 4.78 → 4.66 | 5.11 → 5.27 |

低预算首姿的变化仅出现在同一 6QBH 种子的两个预算。四轨迹两个路径各有一项 Best@5 退化。每个完整输出均在共同受体框架内核对化学图对称重原子 RMSD，使用沉积坐标与独立结合态 SDF 两种映射复核。原始模型、极值、四分位数、P90/P95、各条件的 oracle 和 regret 全部保留。原有计算预算下能消除本组超远闭环，普遍姿势改善尚缺证据。[逐条件配对差](paired.csv) · [修复前的每个姿势](before/poses.csv) · [修复后的每个姿势](after/poses.csv) · [连续分布](after/distributions.csv)

1IEP 四个非大环对照共 51 个输出，两版逐字节一致。人 FKBP12 两个受体、八个化合物、三枚种子的 48 个冻结首姿，按 AD4 原流程比较同一 Holt Ki 系列，48 对分数逐一相同。先在各受体内取种子中位数，然后两个受体均分，三倍且区间分离的分子对为 6/15；七个点测量 Spearman −0.179、Kendall τ-b −0.143。此处固定坐标检验的是最终打分约定。没有估计绝对标准结合自由能。[逐条分数](fixed_affinity.csv) · [每个受体与汇总统计](fixed_affinity_stats.csv) · [测定和原始坐标](../../cubic_grid/fkbp_ad4/README.md)

## 复现

从仓库根目录执行。两个二进制分别由父提交 `3f21e02` 和包含本目录的修复版本用 `python tools/build.py` 构建。生成输入使用原有的冻结压缩包，无需下载 PDB。环境需 NumPy、SciPy、RDKit 和项目的 C++/Boost 工具链。

```bash
D=benchmarks/mechanism_interactions/macrocycle_selection
W=$(mktemp -d)
mkdir "$W/inputs"
tar xf benchmarks/concerted_ring/new_inputs.tar.xz -C "$W/inputs"
python "$D/run.py" "$W/inputs" "$W/before" --binary /path/to/parent/vina
python "$D/run.py" "$W/inputs" "$W/after" --binary /path/to/fixed/vina --continued
python benchmarks/macrocycle_domain/fixed_affinity.py "$W/affinity" /path/to/parent/vina /path/to/fixed/vina
```

仅复算发布的坐标时，将某一臂的 `runs.csv` 和 `raw.tar.xz` 复制到一个空输出目录，解包后去掉上面相应命令的 `--binary` 参数及路径。修复臂保留 `--continued`。原始日志保留实际搜索与精修求值和进程配置。宏环标签 C++ 回归另覆盖 Vina、Vinardo、AD4 三种打分、零权重、自定义权重、预计算上限、生效的软化上限、20/20.48 Å 两侧及更长距离，检查能量、径向导数、有限差分和成对等反力。该新增断言在父版失败，修复版通过；原有普通截断和光滑切换数值测试通过。

从已发布日志和坐标复现配对表时可以直接运行下列命令。这里复用发布的运行时间，重新计时的 CPU 秒、墙钟秒自然会不同。

```bash
D=benchmarks/mechanism_interactions/macrocycle_selection
W=$(mktemp -d)
mkdir "$W/inputs"
tar xf benchmarks/concerted_ring/new_inputs.tar.xz -C "$W/inputs"
for arm in before after; do
    mkdir "$W/$arm"
    cp "$D/$arm/runs.csv" "$W/$arm/"
    tar xf "$D/$arm/raw.tar.xz" -C "$W/$arm"
    flags=()
    if [ "$arm" = after ]; then flags+=(--continued); fi
    python "$D/run.py" "$W/inputs" "$W/$arm" "${flags[@]}"
    for table in poses cells distributions acyclic; do
        cmp "$W/$arm/$table.csv" "$D/$arm/$table.csv"
    done
done
python "$D/paired.py" "$W/before/cells.csv" "$W/after/cells.csv" "$W/paired.csv"
cmp "$W/paired.csv" "$D/paired.csv"
```

上面的 AD4 `fixed_affinity.py` 命令会在 `$W/affinity/` 同时输出 `fixed_affinity.csv` 和 `fixed_affinity_stats.csv`，分别与本目录的同名归档比较即可核对逐分数和逐受体相关性。
