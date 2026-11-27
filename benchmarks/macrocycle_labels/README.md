# 大环闭环编号的物理等价性

PDBQT 用同一分子中的 `G0…G3` 虚拟原子和 `CG0…CG3` 真实碳标记切开的共价键。每条断键有两条相向的闭环约束。只改变一组标签的数字，原子坐标、化学图和对接问题完全相同。在两个分别准备的分子中复用一个数字，这两条共价键依然各自闭合。

## 机制和取舍

父版本 `f04bd7f` 的 XS 枚举交替排列疏水碳、极性碳和虚拟原子；Vina 与 Vinardo 半径数组却先排八个碳再排四个虚拟原子，使 CG2 极性碳和两种 CG3 碳取得 0 Å 半径。缓存网格对闭环碳使用普通碳通道，显式精修和分子内对势则直接读取这张表。此外，`model::append` 根据新旧分子所用编号创建跨分子的虚构闭环约束，并省略两个真实分子之间本应计算的非键接触。

修复对齐了半径数组与枚举，并让独立输入之间只构建常规原子对。分子内部闭环对在各自解析时已经建立。12 元烃环回归测试要求两个分子只保留四条分子内闭环约束、全部 144 个跨分子碳对，四种重编号产生逐元素相同的显式能量和梯度。该检查覆盖 XS 和 AD 两种配对方式。势能/梯度实测采用 Vina，Vinardo 另检查半径映射。

另测了一个范围更大的化学改动（full）。父版把闭环标记碳和虚拟原子视作其相邻碳的杂原子邻居，疏水判据又未识别带标记的疏水碳。full 同时改变这两处，因此会改变最常用的 G0 能量面。这个变体在两个输入上都满足编号等价性，但两处 Top-1 均明显退化，补丁仅在实验归档中。默认仅采用保留普通单环 G0 结果的最小修复（fixed）。

机制依据见 [Holcomb 等的大环闭环方法](https://doi.org/10.1017/qrd.2022.18) 和 [Vina 1.2.0 的方法与 BACE 验证](https://doi.org/10.1021/acs.jcim.1c00203)。

## 真实结构和计算对照

所有方法使用同一输入、1 条轨迹、1 个线程、种子 458921、搜索上限 12000 和输出上限 5。旧版在整步完成后判断上限，实际次数另见 CSV。下面单位为 Å，每格为 **首姿 / 全部输出中最接近晶体的姿势**。两种编号是同一化学分子的编码重复，每个体系只算一个观测。

| 体系、编号 | 官方 1.2.7 | 父版 | 最小修复 | full 变体 |
|---|---:|---:|---:|---:|
| 5QCT BACE1，G0 | 4.723 / 4.723 | 4.723 / 4.723 | 4.723 / 4.723 | 10.232 / 10.232 |
| 同一输入，G3 | 9.311 / 8.017 | 9.335 / 7.992 | 4.723 / 4.723 | 10.232 / 10.232 |
| 1YET HSP90α，G0 | 4.566 / 4.566 | 4.567 / 4.567 | 4.567 / 4.567 | 9.388 / 9.388 |
| 同一输入，G3 | 9.941 / 6.969 | 9.941 / 3.662 | 4.567 / 4.567 | 9.388 / 9.388 |

- BACE1 使用官方教程的刚性受体、配体与 20 Å 对接盒。以 379 个共同 Cα 将 [5QCT](https://www.rcsb.org/structure/5QCT) 的 A 链配准到教程受体，残差 0.00049 Å。38 个重原子的教程输入与实验配体只差 0.00051 Å，故此例含有结合态起始几何信息。实验 E4Y 的带立体标记化学图只有一个自同构，姿势评价直接在受体坐标系中计算。实际输出的闭环碳键长和虚拟原子残差见 `geometry.csv`。
- 独立家族使用 [1YET](https://www.rcsb.org/structure/1YET) 的人 HSP90α–geldanamycin 复合物。配体从 CCD 化学图经 ETKDG、MMFF 和 Meeko 准备，初始几何独立于晶体。24 Å 盒中心取自已知口袋。使用既有数据准备协议，仅对参考答案匹配所有化学图允许的立体等价映射。准备状态及晶体重原子映射见归档 `hsp/1YET/preparation.json`。
- 对父版 G0 产生的全部 6 个候选冻结坐标、输入拓扑和统一的零游离态参考值，用 `--pose_coordinates` 逐一复算。只重命名为 G3，父版分数变化的绝对值达 1.955 kcal/mol。最小修复使全部 6 个候选在两种编号下得到相同分数，且保留各自父版 G0 值。相同坐标的 BACE 教程输入使用默认游离态时，父版和官方版给出 −7.628（G0）和 −7.459（G3），最小修复统一为 −7.628。此处是表示等价性试验，没有实测亲和力终点。
- 两个不同配体共同对接的技术对照使用教程 BACE 配体和仓库的环十二烷。父版在两个输入共用 G0 时生成 8 个虚构的跨分子闭环对，正常跨分子对由 456 个减少到 440 个；改用 G0/G1 后得到不同输出。修复后，两次输入的全部坐标与评估计数分别逐项相同，且复现父版独立编号的坐标。这个组合用于验证计算拓扑，无实验共结合标签。
- 对保留的 G0 基线，最小修复的单配体输出坐标及每阶段计数均与父版相同。BACE1 的搜索/精修评估为 12120/177，HSP90α 为 12139/533。父版的 G3 计数分别为 12131/312 和 12196/430。搜索和显式精修调用数与端到端单次壁钟时间见 `runs.csv`、`hsp90.csv`。官方二进制未输出相应计数。壁钟为单次计时，无法据此断言速度提升。
- 两个体系的所有已输出候选均超过 2 Å。HSP90α 原父版 G3 输出曾包含一个 3.662 Å 的候选，修复后池中最佳退为 4.567 Å。完整结果支持删除表示依赖，并保留已验证的 G0 行为；对大环的采样、2 Å 恢复率、分子排序或条件排名没有正收益证据。

## 复现

源码和数据下载地址及 SHA-256 在 `sources.json`。归档附完整原始输入、日志、PDBQT 和非默认变体补丁。官方二进制使用 [Vina 1.2.7 Linux x86_64](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7)，其 SHA-256 为 `f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644`。实验使用 GCC C++14、-O3、Boost 1.83，分析使用 RDKit 2025.9.3 与 NumPy 2.3.5。若重新准备 HSP90α，另需 Meeko 0.7.1、OpenMM 8.6.1、PDBFixer 1.12.0 和 Gemmi 0.7.4。

在仓库根目录运行以下命令。`OFFICIAL_VINA` 为官方二进制的绝对路径。系统没有标准 Boost 安装时，给三次构建增加 `--boost-prefix <路径>`。

```sh
D=benchmarks/macrocycle_labels
W=$(mktemp -d)
tar -xJf "$D/observations.tar.xz" -C "$W"
revision=$(git log -1 --format=%H --diff-filter=A HEAD -- "$D/observations.tar.xz")
git worktree add --detach "$W/parent" f04bd7f
git worktree add --detach "$W/fixed" "$revision"
git worktree add --detach "$W/full" "$revision"
git -C "$W/full" apply "$W/full-typing.patch"
for arm in parent fixed full; do
  python3 "$W/$arm/tools/build.py" --output-dir "$W/bin-$arm"
done
export OPENMM_CPU_THREADS=1
python3 "$D/run_labels.py" --upstream-data "$W/tutorial" --out "$W/bace" \
  --parent "$W/bin-parent/vina" --fixed "$W/bin-fixed/vina" --full "$W/bin-full/vina" --official "$OFFICIAL_VINA"
python3 "$D/geometry.py" --upstream-data "$W/tutorial" --structure "$W/bace/5QCT.pdb" --out "$W/bace"
python3 "$D/holdout.py" "$W/hsp" --parent "$W/bin-parent/vina" --fixed "$W/bin-fixed/vina" \
  --full "$W/bin-full/vina" --official "$OFFICIAL_VINA"
python3 "$D/fixed_pool.py" --bace-work "$W/bace" --hsp-work "$W/hsp" --upstream-data "$W/tutorial" \
  --parent "$W/bin-parent/vina" --fixed "$W/bin-fixed/vina" --full "$W/bin-full/vina"
```

若需从原始 1YET / GDM 数据重做输入，在运行 `holdout.py` 对接前执行 `python3 "$D/holdout.py" "$W/hsp" --prepare`。可单独重算归档的 `geometry.csv`；无需重新对接。标准 CI 从同一份源码构建 `check_macrocycle_labels.cpp` 并验证配对、半径与梯度。
