# 从原始输入、归档坐标及运行日志复算

以下命令均在仓库根目录执行。将 `WORK` 设为一个新的空目录的绝对路径。Python 3.12.14，OpenMM 8.6.1、PDBFixer 1.12.0、Meeko 0.7.1、RDKit 2025.09.3、NumPy 2.3.5、SciPy 1.17.0、Gemmi 0.7.4。使用单线程 CPU。其他平台和软件版本可能改变局优坐标。分阶段计时及各项预算见[计算成本](results/research_costs.json)。

## 解开原始输入与本轮证据

```sh
set -eu
: "${WORK:?Set WORK to a new empty directory}"
export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
D=benchmarks/receptor_linkage/aurora/internal_loops
A=benchmarks/receptor_linkage/aurora
mkdir -p "$WORK"
tar xf benchmarks/torsion_kinase_affinity/evidence.tar.xz -C "$WORK" \
  raw/2X6E.pdb raw/2J4Z.pdb prepared/2X6E-protein.pdb prepared/2J4Z-protein.pdb
tar xf "$A/docking.tar.xz" -C "$WORK" prepared
tar xf "$A/references.tar.xz" -C "$WORK"
tar xf "$A/probe_raw.tar.xz" -C "$WORK" raw/at9283_chemical.bonds raw/danusertib_native.bonds
tar xf "$A/state_windows/coordinates.tar.xz" -C "$WORK"
tar xf "$D/evidence.tar.xz" -C "$WORK"
```

`base/` 含 24 次缺段受体盲搜，`templates/` 含 24 次模板补环盲搜。初次的十次固定姿势干分在 `fixedbase/`。`finalfixed/` 保存统一比较的四个受体各十次定点评分。汇总时只使用后一份含干分的固定池目录，避免重复计入。

## 只用归档数据核查物理和评价

```sh
python "$D/verify.py" "$WORK" > "$WORK/ff.json"
cmp "$WORK/ff.json" "$D/results/ff14sb_independent.json"
python "$D/donor_check.py" "$WORK" "$WORK/donors_out" --models \
  loop43103/completed.pdb loop43104/completed.pdb loop43103/repaired.pdb loop43104/repaired.pdb
for f in donor_distances.csv donor_frames.csv donor_model_rmsd.csv donor_sources.json; do
  cmp "$WORK/donors_out/$f" "$D/results/$f"
done
python "$D/summarize.py" "$WORK/statistics" "$WORK/base" "$WORK/templates" "$WORK/finalfixed"
for f in cells modes distributions fixed; do
  cmp "$WORK/statistics/$f.csv" "$D/results/$f.csv"
done
python "$D/check_predictions.py" "$WORK" --results "$WORK/base" "$WORK/templates" "$WORK/finalfixed" > "$WORK/independent.json"
cmp "$WORK/independent.json" "$D/results/independent_predictions.json"
```

上述几步分别复核氨基酸序列、两根异常肽键的力和能量、七份实验供体的几何对应、所有发布姿势的 RMSD、定点评分差以及受体中新增原子的化学特征。读者还可单独检查[七供体实验](DONORS.md)。

## 重新生成原始受体的制备元数据

```sh
for w in "$WORK/producer_parent" "$WORK/producer"; do
  mkdir -p "$w/raw"
  cp "$WORK/raw/2X6E.pdb" "$WORK/raw/2J4Z.pdb" "$w/raw/"
done
git show 019c9ff:benchmarks/torsion_kinase_affinity/prepare_receptors.py > "$WORK/producer_parent/parent.py"
python "$WORK/producer_parent/parent.py" "$WORK/producer_parent"
python benchmarks/torsion_kinase_affinity/prepare_receptors.py "$WORK/producer"
python "$D/verify.py" "$WORK" --exports "$WORK/producer_parent" "$WORK/producer" > "$WORK/exports.json"
cmp "$WORK/exports.json" "$D/results/export.json"
```

## 从序列和独立实验模板重建受体

将重建放在新目录，保留原始归档供对照。PDBFixer 的两套提议在读取任何配体评分前未通过共价几何检查。对它们分别重做直接补环、保留原氢的撤约束，以及重新放置手性 CH 氢后撤约束。后两步的文件名分别为 `first_release.pdb` 和 `repaired.pdb`，详见[首次失败及修复对照](REPAIR_PROTOCOL.md)。

```sh
REBUILD="$WORK/rebuild"
mkdir -p "$REBUILD"
ln -s "$WORK/raw" "$REBUILD/raw"
ln -s "$WORK/prepared" "$REBUILD/prepared"
for seed in 43103 43104; do
  python "$D/build.py" "$REBUILD" "$seed"
  python "$D/repair.py" "$REBUILD" "$seed"
  python "$D/repair.py" "$REBUILD" "$seed" --reset-h
done
for donor in 1OL5 1OL7; do
  python "$D/template.py" "$REBUILD" "$donor" --initial
  python "$D/template.py" "$REBUILD" "$donor"
  python "$D/template.py" "$REBUILD" "$donor" --trans-model
done
for donor in 1OL5 1OL7; do
  target="template$donor/trans_model"
  cmp "$WORK/$target/closed.pdb" "$REBUILD/$target/closed.pdb"
  for stage in restrained released; do
    cmp "$WORK/$target/$stage/repaired.pdb" "$REBUILD/$target/$stage/repaired.pdb"
  done
done
for ref in 1OL5/trans_model/restrained 1OL7/trans_model/restrained 1OL7/trans_model/released; do
  cmp "$WORK/template$ref/receptor.pdbqt" "$REBUILD/template$ref/receptor.pdbqt"
done
test ! -e "$REBUILD/template1OL5/trans_model/released/receptor.pdbqt"
```

`--initial` 对应从供体出发的一次扭转闭合，默认模式进一步求三残基环的多个闭合解，`--trans-model` 固定先前的权重加入肽键和手性约束，随后撤去约束。脚本按阶段存储化学检查结果，只有通过验收且完成原子赋型的阶段生成受体 PDBQT。1OL5 撤约束构象中的 Arg304 翻转可从 `template1OL5/trans_model/released/template.json` 核查。[模型定义与验收条件](TEMPLATE_PROTOCOL.md)

## 重做固定池与盲搜

令 `PRO` 和 `OFFICIAL` 指向两个可执行文件的绝对路径。官方程序为 AutoDock Vina 1.2.7，SHA-256 为 `f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644`。本轮 Pro 实际二进制 SHA-256 为 `6a56ec3db8c23b2e8bea4f65eb94fe3a4ccc2b272b4ae27c0da8113c078b4304`，对应本轮基线的对接核心，可用 `tools/build.py` 构建。定点评分采用完整配体键表和共同未结合态零点。盲搜采用历史默认推断图、原盒及四枚种子。

```sh
: "${PRO:?Absolute path to the Pro binary}" "${OFFICIAL:?Absolute path to AutoDock Vina 1.2.7}"
OUT="$WORK/rerun"
mkdir -p "$OUT"
python "$D/replay.py" blind "$WORK" "$OUT/base" "$PRO" --official "$OFFICIAL" --arms dry official double
python "$D/replay.py" blind "$WORK" "$OUT/templates" "$PRO" --official "$OFFICIAL" \
  --arms 1OL5_restrained 1OL7_restrained 1OL7_released
python "$D/replay.py" fixed "$WORK" "$OUT/finalfixed" "$PRO" \
  --arms dry 1OL5_restrained 1OL7_restrained 1OL7_released
python "$D/summarize.py" "$OUT/summary" "$OUT/base" "$OUT/templates" "$OUT/finalfixed"
cmp "$OUT/summary/modes.csv" "$D/results/modes.csv"
python "$D/check_predictions.py" "$WORK" --results "$OUT/base" "$OUT/templates" "$OUT/finalfixed" > "$OUT/independent.json"
cmp "$OUT/independent.json" "$D/results/independent_predictions.json"
```

归档的逐次进程计时不可逐字节重现。最后两项比较使用逐姿表及独立几何/评分核查。[全池交叉评分、远移新增原子的数值对照及额外开销](pool_control/README.md)
