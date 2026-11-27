# 复算与重新对接

从完整 Git 检出执行。默认分支源码和各研究版本均能独立构建，不需要给生产 CLI 增加开关。研究目录位于仓库根目录下。下列第一条命令仅依赖 Python 标准库和 sha256sum，会核对完整 PDBQT 散列、日志中的六项工作计数、模式数、公开的结果归档，以及从逐姿势数据生成的统计和固定池干预。

```bash
D=$(pwd)/benchmarks/macrocycle_proposal
python "$D/verify.py"
W=$(mktemp -d)
bash "$D/build.sh" "$W" --jobs 4
```

实验使用 GCC 14.2、C++14、-O3、Boost 1.83。自定义 Boost 可额外传入 `--boost-prefix /path/to/boost`。上游固定于 `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`，用无插桩生产输出逐字节比对同源码计数版。离线编译上游可见 `benchmarks/pair_gradient/build_arms.py --help` 的本地 checkout 选项。结构指标依赖 NumPy、SciPy、RDKit。重新准备受体和配体另需 Meeko、Gemmi、PDBFixer、OpenMM，分试验的版本记录随各子目录提供。一次预测的原始时长包含 CLI 启动、输入解析、建图、MC、局优、精修和写出。输入制备成本单列，所有方法读取逐字节相同的预处理输入。

## 历史四个复合物

```bash
mkdir "$W/inputs"
tar -xJf benchmarks/macrocycle_labels/observations.tar.xz -C "$W/inputs"
tar -xJf benchmarks/bond_distance/macro/extra_inputs.tar.xz -C "$W/inputs"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
U="$W/up/upstream/build/vina"
UC="$W/up/upstream-counts/build/vina"
P="$W/bin-parent/vina"
C="$W/bin-ccd/vina"
R="$W/bin-reset/vina"
L="$W/bin-lock/vina"
python "$D/panel.py" "$W/inputs" "$W/first" "$U" "$UC" "$P" "$R" "$C"
python "$D/cost.py" "$W/inputs" "$W/first" "$W/first/timing.csv" "$P" "$R" "$C"
python "$D/fixed_scores.py" "$W/inputs" "$W/first" "$W/first/fixed_ccd.csv" "$P" ccd "$C"
python "$D/persistent.py" "$W/inputs" "$W/first" "$W/persistent" "$P" "$L"
```

数值探针用同一批生产对象文件链接。下例的两个额外输入来自外部面板，5K0X 来自下方双链归档。自定义 Boost 时为探针和观察器分别设置 `BOOST_PREFIX` 和 `--boost-prefix`。观察器运行中的额外坐标输出不计入生产速度。它必须重现未插桩的完整输出池和六项计数。阶段 `raw` 和 `proposal` 均在第一次普通突变后。对持续固定环臂，二者都已经位于受限空间。

```bash
mkdir "$W/ext" "$W/mer"
tar -xJf "$D/external/inputs.tar.xz" -C "$W/ext"
tar -xJf "$D/holdout/observations.tar.xz" -C "$W/mer"
for arm in parent ccd reset lock; do
    python "$D/observer.py" "$W/src-$arm" "$W/bin-$arm/obj" "$W/observe-$arm"
done
for arm in ccd reset lock; do
    items=("$W/ext/5T9W" "$W/ext/7MBO")
    if [ "$arm" = lock ]; then items+=("$W/mer/inputs/5K0X"); fi
    bash "$D/numeric/check_geometry.sh" "$arm" "$W/src-$arm" "$W/bin-$arm" "$W/inputs" "$W/$arm-check.csv" "${items[@]}"
done
python "$D/stages.py" "$W/inputs" "$W/first" "$W/stages" --parent "$W/observe-parent/vina" --ccd "$W/observe-ccd/vina" --reset "$W/observe-reset/vina"
python "$D/stages.py" "$W/inputs" "$W/persistent" "$W/lock-stages" --cap 12000 --lock "$W/observe-lock/vina"
```

直接解开 `historical/poses.tar.xz` 可检查所有历史输出池及阶段坐标。第一轮的 CSV 在 `historical/`，持续约束数据在 `historical/persistent/`，观察器输出在三个包含 stage 的子目录。raw 文件的散列与 runs.csv 一一对应。输出上限相同，实际产生的候选数可以不同。

## 新结构、同测定 Kd 与校准计时

MerTK 使用实验中提供的 A 链配体，用同一沉积的 B 链表达晶体接触。没有拿它代表单体生物学组装。下载的 PDB、CCD、结合态 SDF、两条受体链、排除记录和完整输出均在双链归档。以下重跑使用与归档相同的准备文件；若从源文件开始，先执行 `holdout/prepare.py WORK` 和 `holdout/assembly.py WORK`。

```bash
python "$D/holdout/run_assembly.py" "$W/mer/inputs" "$W/mer-run" "$U" "$UC" "$P" "$R" "$C" "$L"
python "$D/external/external_stages.py" "$W/mer/inputs" "$W/mer-run" "$W/mer-stages" --parent "$W/observe-parent/vina" --ccd "$W/observe-ccd/vina" --reset "$W/observe-reset/vina" --lock "$W/observe-lock/vina"
```

新 PPIA / FXIa 的两个独立 RMSD 算法、外部池、阶段验证与全部排除项见[外部复核](external/README.md)。纯计时选上限、三枚验证种子和每对的实际墙钟见[等时间实验](external/wall_match/README.md)。FKBP12 两配体两受体的重对接、交叉对接、固定姿势十项能量与 ITC 分差见[实测亲和力](fkbp_affinity/README.md)。保持环刚性而撤销提议重分配的对照见[空提议消融](null_control/README.md)。所有重复计时仅用于测量成本，未合并候选。

重新计算旧 MMFF 分量：`python benchmarks/frozen_core/strain_transfer_audit/replay.py "$W/mmff" --jobs 2`。本轮 42 个分子、168 个固定首姿得到的 components.csv、correlations.csv、panel.csv、summary.csv、pairs.csv.xz 均与原档案逐字节相同。
