# 仅在首轮局优中保留顺坡平移的四版本结构试验

[事前方案](../stage_local_attribution/PROTOCOL.md) 登记了两个新受体、固定面板、配时可达区间和晋级门槛。研究候选为 `0bcae6ef09857f97656291bb4821d20b3b619bf5`。固定五靶点交叉上，候选相对轮前的首姿均值、池 Oracle 和 P95 均降低。新 ESR1 受体上的同晶首姿均值下降，交叉首姿均值上升 2.266 Å，未通过外检门槛。真实官方和真实父版在固定交叉、外部交叉的平均首姿误差均小于候选。

| 条件 | 每臂格数 | 官方首姿均值 Å | 父版 | 轮前 | 候选 | 候选对轮前改善 / 持平 / 恶化 |
|---|---:|---:|---:|---:|---:|---:|
| 固定交叉，12k | 20 | 2.863 | 2.863 | 6.016 | 5.551 | 12 / 0 / 8 |
| 固定同晶，12k | 20 | 6.794 | 6.794 | 7.175 | 6.730 | 10 / 0 / 10 |
| 固定交叉，历史成本分配 | 20 | 4.709 | 2.536 | 4.746 | 4.741 | 10 / 0 / 10 |
| 原生默认 | 5 | 2.372 | 2.372 | 2.774 | 2.458 | 2 / 0 / 3 |
| ESR1 新受体同晶，12k | 4 | 3.936 | 3.936 | 3.721 | 2.144 | 2 / 0 / 2 |
| ESR1 新受体交叉，12k | 4 | 3.637 | 3.637 | 4.961 | 7.227 | 0 / 0 / 4 |

固定交叉候选首姿 Q1 / 中位 / Q3 / P90 / P95 分别为 2.746 / 5.371 / 6.540 / 11.150 / 13.820 Å。Best@5 4.871、池 Oracle 4.804、Regret 0.747 Å。轮前池 Oracle 4.925、Regret 1.091、P95 14.504 Å。首姿超过 5 Å 的次数在两个 Pro 中均为 11/20。相对官方和父版，候选首姿均高 2.688 Å。候选 65 个固定格全部输出，共 572 个模式。逐格原始首姿和各阶统计在 [fixed/cells.csv](fixed/cells.csv) 和 [fixed/summary.csv](fixed/summary.csv)，逐方向分布与三基线配对数在 [fixed/by_direction.csv](fixed/by_direction.csv)，各格每项配对差在 [fixed/paired.csv](fixed/paired.csv)。完整候选顺序、PDBQT 坐标、日志和据原 runner 恢复的命令在 [fixed/fixed](fixed/fixed)。三份对照沿用 [同一二进制和输入的源档](../geometry_comparison/upstream/four_versions)，包括官方原始超时和按原种子 900 秒补跑，二者均有收据。

## ESR1 两个外留受体

RCSB 沉积 [1QKT](https://www.rcsb.org/structure/1QKT) 为 17β-雌二醇（EST）复合物，原始链 A 的 SEQADV 记录 C381S、C417S、C530S。 [1X7R](https://www.rcsb.org/structure/1X7R) 为染料木素（GEN）复合物，链 A 使用人 ESR1 的 UniProt P03372 编号，链 B 共激活肽最近被省略的非水重原子距配体 15.855 Å。两份链 A 的 DBREF 均与 UniProt 编号一致。历史开发中已使用 ESR1 与 EST，两份受体坐标为本轮外留单位。按共有口袋残基编号和残基类型配对的 65 个 Cα，两个方向配准 RMSD 0.529 Å。交叉值衡量转移参考；两种配体在另一受体中均无本试验使用的直接实测坐标。

原始 PDB 和 CCD、去水制备、映射和搜索输入均随 [external/inputs.tar.xz](external/inputs.tar.xz) 保存。配体初始构象由 CCD 化学图独立生成，沉积配体给出口袋中心和验后参考位置。各臂使用同一 24 Å 立方盒、一个 CPU、三条轨迹、每条 12,000 名义评估、种子 85231/85237。两份输入皆通过，未替换结构。评价在接收受体坐标系中枚举化学等价且保留立体化学的配体映射；独立程序从原始沉积/CCD 和输出 SMILES 索引重建，分别得到 EST 的 1 个和 GEN 的 2 个合法映射。32 格 245 模式的最大 RMSD 复算差为 5.33×10⁻¹⁵ Å。逆序调用四个引擎又实际重复 32 次，所有池逐字节相同，rank、分数、计数全部复现。[全部结果](external/cells.csv) · [所有发布模式](external/ranks.csv) · [逐方向分布](external/by_direction.csv) · [原始日志和坐标](external/raw.tar.xz) · [二次运行](external/repeat)

| 交叉受体方向与种子 | 轮前首姿 → 候选首姿 Å | 轮前 Oracle → 候选 Oracle Å |
|---|---:|---:|
| 1QKT → 1X7R，85231 | 0.561 → 6.289 | 0.561 → 6.289 |
| 1QKT → 1X7R，85237 | 2.648 → 2.648 | 2.648 → 2.648 |
| 1X7R → 1QKT，85231 | 2.446 → 5.771 | 0.888 → 5.771 |
| 1X7R → 1QKT，85237 | 14.187 → 14.198 | 10.769 → 12.944 |

其中一格首姿仅差 0.000326 Å，一格差 0.010625 Å。主要退化发生在两个方向的 85231；即使使用候选池中最接近参考的一项，其 RMSD 也分别为 6.289 和 5.771 Å。四格的交叉 Oracle 从 3.716 升到 6.913 Å，Regret 从 1.244 降到 0.314 Å。一个家族和两个种子不足以估计总体尾部风险。

## 实际 CPU 和分母

每个方向各臂都独立测量 3k/12k（种子 74317）。正斜率的两点直线给出预算 [3000,240000] 的 CPU 可达范围，以四臂 12k 的最大 CPU 为目标。任一臂非正斜率或目标不可达即排除，保留计时与原因。80 个先导用 56.383 CPU 秒，五方向出现非正斜率，两方向超出稳定臂可达区间，仅三方向进入种子 93529 确认。12 个确认运行用 29.786 CPU 秒。候选/官方、父版、轮前的实际 CPU 比中位分别为 0.354、0.453、0.447，各 0/3 落入 [0.8,1.25]。没有前瞻近成本样本可估计误差差值。历史成本分配 20 格与上述标定分开，其实际匹配候选–官方、候选–父版、候选–轮前分母分别为 10/20、3/20、10/20，对应 Δ 首姿分别 −1.097、+2.292、−1.086 Å。[可达区间及排除](fixed/caps.json) · [实际配对](fixed/actual_cost_pairs.csv) · [近成本分母](fixed/actual_cost_summary.csv)

固定候选 65 格包含读入、建图、搜索、精修与发布，合计 175.291 CPU 秒、逐作业墙钟合计 213.443 秒，最大 RSS 43,756 KiB。五靶点固定交叉 20 格分别耗官方 41.964、父版 16.867、轮前 7.737、候选 9.019 CPU 秒。原生默认五格以各程序的搜索/返回默认值运行，保留共同口袋、种子、CPU 和 Pro 只读计数。归档的官方初次默认超时另耗 143.646 CPU 秒。每格 CPU、墙钟、RSS 以及两个 Pro 的 hunt / second / explicit 联合求值数在 cells.csv 中；官方/父版缺少原生计数，历史已独立核验的十个交叉格被动计数见四版本档案。每格完整阶段 CPU 不可分解。

新受体制备共 5.129 CPU 秒、18.931 墙钟秒（含下载），从固定原件离线重做用 6.993 CPU 秒，最高 RSS 184,940 KiB，两份受体/配体及坐标映射均逐字节一致。外部搜索每臂八格分别为官方 10.731、父版 3.667、轮前 2.268、候选 2.834 CPU 秒；候选 hunt / second / explicit 求值分别为 205121 / 83368 / 1801，轮前为 210704 / 77746 / 1614。无搜索后重排。外部搜索二次独立执行额外 25.053 CPU 秒。固定交叉输入的历史独立制备计时为 22.817 CPU 秒，同晶旧输入缺制备计时。[阶段原始时间](external/preparation_replay.json) · [独立复算成本](analysis_cost.json)

## 重建及复现

从仓库根目录运行。四个二进制、候选和轮前的源码树、version / help / 默认行为、编译器和 Boost 哈希见 [fixed/builds.json](fixed/builds.json)。四版本构建命令中的 candidate 对应 `5847b1b`，与本轮开始前 `666f49d` 的 src 树相同，作为此轮 stable。按原版构建脚本取得官方和父版，另用研究提交的代码单独编译候选。C++ 需要 GCC 14.2 及 Boost 1.83。Python 核查用 NumPy / SciPy，原始化学制备与独立图读数另需 RDKit、Gemmi、Meeko、OpenMM、PDBFixer。设置指向独立空目录的 WORK 和 Boost 前缀 BOOST。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
D=benchmarks/stage_local_structure
P=benchmarks/search_trial_panel/panel.py
UP=benchmarks/geometry_comparison/upstream/four_versions
REV=0bcae6ef09857f97656291bb4821d20b3b619bf5
mkdir -p "$WORK"/{inputs,source}
bash "$UP/build.sh" "$WORK/base" "$BOOST"
git archive "$REV" src tools/build.py | tar xf - -C "$WORK/source"
python "$WORK/source/tools/build.py" --boost-prefix "$BOOST" --output-dir "$WORK/new" --jobs 2
official="$WORK/base/official/vina"
parent="$WORK/base/parent/vina"
stable="$WORK/base/candidate/vina"
candidate="$WORK/new/vina"
tar xf "$D/external/inputs.tar.xz" -C "$WORK/inputs"
python "$D/reprepare.py" "$WORK/inputs" "$WORK/prep-repeat"
python "$D/external.py" check "$WORK/inputs" "$D/external"
python "$UP/offline.py"
python "$D/cost.py" check "$D/fixed"
python "$UP/run_panel.py" candidate "$candidate" "$WORK/runs/fixed" --jobs 3
for arm in official parent stable candidate; do
  python "$D/external.py" run "$WORK/inputs" "$WORK/external" "$arm" "${!arm}"
  python "$P" pilot "$WORK/runs" "$arm" "${!arm}"
done
python "$D/external.py" collect "$WORK/inputs" "$WORK/external" "$WORK/external-tables"
python "$D/cost.py" freeze "$WORK/runs"
# If at least one direction has caps, execute the matched stage for every arm.
if python -c 'import json,sys;sys.exit(not any("caps" in r for r in json.load(open(sys.argv[1]))["directions"].values()))' "$WORK/runs/caps.json"; then
  for arm in official parent stable candidate; do
    python "$P" matched "$WORK/runs" "$arm" "${!arm}"
  done
fi
python "$D/cost.py" collect "$WORK/runs" "$WORK/fixed-tables"
python "$D/directions.py" "$WORK/fixed-tables"
python "$D/directions.py" "$WORK/external-tables"
```

重复运行重新测得的 CPU 会改变可达范围与确认预算。原始配时结果的确定性复算使用归档收据和 caps；以上新执行单独记账。可以用 `python "$D/prepare.py" "$WORK/downloaded"` 从公共沉积重新生成输入。固定面板复算 65 个候选正式格、80 个先导格和 12 个确认格的全部 572、704、107 个模式，最大 RMSD 差 8.0×10⁻¹⁵ Å。外部二次运行可于另一个 work 目录按同一命令执行，然后调用 `repeat.py INPUTS SECOND_RUN_DIR PRIMARY_TABLES NEW_REPEAT_TABLES`。


## 生产源码恢复的独立全量重放

独立从公共提交 `1c3141238511ba98e13232ffff8f2ba51efef3d1` 导出 src 与构建脚本，在空目录编译。其 src 树与轮前相同（`f22a5de502749b4743ab0ce7fc2c9cbcbfdf7948`）。再执行完整 65 格，65/65 发布池、572 个模式和同种子 hunt/second/explicit 求值、精修调用及晚期预评分计数均与轮前完全一致。独立重放耗 149.865 CPU 秒，逐作业墙钟合计 155.202 秒，峰值 RSS 43,732 KiB。可执行文件含构建绝对路径，两次编译使用不同长度的目录。在二进制成品中替换源码路径和 ELF build ID 后哈希仍不同。随后从同一洁净导出重新编译，用 GCC 的 `-ffile-prefix-map` 将编译时源码路径映射为旧路径，两份 1,420,472 字节 ELF 完全一致，SHA-256 均为 `f6ff9e0d9d6f0bc6e77e65409f9c7f88fccd826eca064a49e9c91a9a73496f9f`。[路径复现编译收据](restored/path_reproducibility.json)[构建日志](restored/build.txt) · [源码/二进制/重放收据](restored/verification.json) · [原始坐标与日志](restored/raw.tar.xz)

```bash
mkdir "$WORK/release-source"
git archive 1c3141238511ba98e13232ffff8f2ba51efef3d1 src tools/build.py | tar xf - -C "$WORK/release-source"
python "$WORK/release-source/tools/build.py" --boost-prefix "$BOOST" --output-dir "$WORK/release" --jobs 4
python "$UP/run_panel.py" stable "$WORK/release/vina" "$WORK/release-runs" --jobs 4
python "$D/restore.py" "$stable" "$WORK/release/vina" "$WORK/release-runs" "$WORK/release-check"
```

另从两份真实上游提交独立导出并重新编译。复现历史路径字面量后，官方 859,976 字节和父版 860,096 字节的 ELF 均与参评文件逐字节相同，版本、普通及高级帮助、全部打印默认值一致。两次编译分别耗用 67.71 CPU 秒。[命令、源码树、编译器、默认值、费用及 ELF 哈希](restored/upstream_rebuild.json)
