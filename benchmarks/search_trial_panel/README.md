# 主角度约束的固定结构与实际成本对照

候选源 `0bb827ebfe69a32504b74488730f80a21495eccf` 将 MC 网格局优的首个 Armijo 步限制在旋转和扭转的主角度内。与当轮稳定源 `ed0d42d` 相比，五靶点双向交叉的 20 个固定预算格中，池内最低 RMSD 平均下降 1.064 Å，首姿平均下降 0.920 Å，排序损失平均增加 0.144 Å。首姿 P95 从 14.504 增至 14.860 Å，超过方案的晋级界限。相对真正的官方与父版，交叉首姿均值仍高 2.233 Å。完整预登记在 [PROTOCOL.md](PROTOCOL.md)，引擎 SHA、源码树、可执行文件 SHA-256、编译器、动态依赖与实际 CLI 默认值在 [builds.json](results/builds.json)。

| 设置（每臂格数） | 官方 Top-1 均值 Å | 父版 | 轮前稳定 Pro | 研究候选 | 候选相对稳定改善/持平/退化 |
|---|---:|---:|---:|---:|---:|
| 交叉，12k（20） | 2.863 | 2.863 | 6.016 | 5.096 | 13/0/7 |
| 同晶，12k（20） | 6.794 | 6.794 | 7.175 | 6.143 | 13/0/7 |
| 交叉，历史成本分配（20） | 4.709 | 2.536 | 4.746 | 2.791 | 14/0/6 |
| 原生默认（5，含官方原有补跑） | 2.372 | 2.372 | 2.774 | 2.364 | 2/0/3 |
| 交叉，前瞻成本分配（5 个可标定方向） | 5.786 | 3.344 | 3.469 | 2.827 | 3/0/2 |

前三臂复用上一轮相同源码、二进制、结构、种子与预算的实际发布和原始计时。稳定角色对应当轮开始前的分支，其 src 与 `5847b1b` 相同。默认组中官方 1FKG→1FKH 初次超时，配对表用原种子 900 秒补跑。初次 143.646 CPU 秒额外支出和未发布的事实留在[原始主面板](../geometry_comparison/upstream/four_versions/results/cells.csv)。候选独立执行固定设置合计 65 次，65 格均有输出，二次读取 572 个模式。原生默认在相同口袋、种子和一核下执行，附加 Pro 只读计数。

## 交叉误差的分布和来源

| 引擎 | 均值 | Q1 | 中位数 | Q3 | P90 | P95 | Best@5 | 池 Oracle | Regret | 首姿 >5 Å |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 官方/父版 | 2.863 | 0.779 | 2.567 | 5.065 | 5.649 | 5.867 | 2.502 | 2.482 | 0.380 | 5/20 |
| 轮前 | 6.016 | 2.564 | 5.400 | 8.242 | 12.602 | 14.504 | 5.061 | 4.925 | 1.091 | 11/20 |
| 候选 | 5.096 | 1.272 | 3.008 | 6.692 | 14.270 | 14.860 | 4.172 | 3.861 | 1.235 | 7/20 |

以上几何量单位为 Å；分位数沿用固定基准的 type-2 经验分位数。交叉方向以实验口袋 Cα 配准产生参考，在接收受体坐标系计算化学等价映射下的重原子 RMSD。同晶和交叉数据均已在历史开发中查看。

| 靶点（4 格） | 官方/父版首姿 | 轮前首姿 | 候选首姿 | 轮前→候选池 Oracle | 轮前→候选 Regret |
|---|---:|---:|---:|---:|---:|
| DYRK1A | 0.631 | 5.011 | 0.805 | 4.726→0.805 | 0.285→0.000 |
| EPAS1 | 2.043 | 5.338 | 8.562 | 3.180→4.465 | 2.158→4.097 |
| FABP4 | 3.434 | 5.544 | 3.380 | 4.697→2.584 | 0.848→0.795 |
| FKBP1A | 5.675 | 6.043 | 3.990 | 4.565→3.809 | 1.478→0.181 |
| RBP4 | 2.531 | 8.143 | 8.744 | 7.459→7.642 | 0.684→1.101 |

RBP4 的 3FMZ→1RBP 与反向在种子 85237 的候选池 Oracle 分别为 12.321 和 13.551 Å；这两个已发布池即便选择最接近参考的一项也有大的误差。EPAS1 的 5TBM→4GHI、种子 85231 则有 3.049 Å 的候选，发布首位为 14.297 Å。不同误差来源可在冻结池与固定源状态下分开测试。所有运行的配对、目标分位数、成本见 [cells.csv](results/cells.csv)、[summary.csv](results/summary.csv)、[paired.csv](results/paired.csv) 与 [paired_summary.csv](results/paired_summary.csv)。候选全部模式见 [fixed/ranks.csv](results/fixed/ranks.csv)，历史三臂的模式与坐标取自[四版本档案](../geometry_comparison/upstream/four_versions/)。

## 真实成本与标定覆盖

候选固定 65 格共用 159.623 CPU 秒，逐作业墙钟之和为 165.238 秒，最高 RSS 为 43,760 KiB。交叉固定 20 格官方、父版、轮前和候选分别为 41.964、16.867、7.737 和 8.569 CPU 秒。该组候选 hunt/second/explicit 联合能量梯度调用数是 619,526 / 101,445 / 4,627，轮前是 635,507 / 85,616 / 4,610。原始 Pro 每轨迹计数可从日志复核。官方/父版缺原生计数，历史已验证的 10 个交叉格计数见[被动计数重放](../geometry_comparison/upstream/four_versions/counter_results/)。以上 CPU 均计入读入、建图、搜索和发布。冻结交叉面板输入在上一轮独立制备共 22.817 CPU 秒，四臂共享；同晶历史制备时间未测。读数未包含编译、下载及事后 RMSD 计算。

成本标定用 80 个独立真实执行（64.468 CPU 秒，逐作业墙钟合计 65.996 秒）。固定种子 74317，在每条轨迹 3k/12k 两个预算测量 CPU，只用计时冻结单格预算。10 个方向中有 5 个方向至少一臂斜率非正；排除项与其涉及的引擎保存在 [caps.json](results/caps.json)。其余 5 个方向均由四臂按各自冻结预算独立在种子 93529 运行，共 20 次（50.578 CPU 秒，墙钟合计 53.938 秒）。运行含启动与建图，单次计时的分辨率与调度波动使此线性标定不稳定。实际候选/官方 CPU 比中位数 0.685，候选/父版 0.777，候选/轮前 0.809。按预登记 [0.8,1.25] 比率取到的配对分别为 2/5、1/5 和 2/5，候选减基线首姿均值分别为 −1.221、−2.149、+1.910 Å。当前近成本证据覆盖很少，且三次比较涉及的方向不完全相同。[逐对实际 CPU 比与入选](results/actual_cost_pairs.csv) · [分母、方向与幅度](results/actual_cost_summary.csv)

## 独立检查与复现

从仓库根目录运行。Python 依赖见历史结构面板，构建使用 GCC 14.2.0 和 Boost 1.83.0。

```bash
export OPENBLAS_NUM_THREADS=1
P=benchmarks/search_trial_panel
D=benchmarks/geometry_comparison/upstream/four_versions
python "$D/offline.py"
python "$P/panel.py" check "$P/results"
```

以下从四版源码重建并实际执行。预留 NEW 为尚不存在的临时目录，BOOST 为 Boost 的前缀。官方/父版/轮前的固定结构格来自上一轮原始真实执行。为复核其命令，可分别用 `"$D/run_panel.py"` 和相应二进制在新的工作目录独立运行。运行时间会影响重新标定的预算；原始成本实验的确定性核查应使用已归档的收据和 caps。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1
P=benchmarks/search_trial_panel
D=benchmarks/geometry_comparison/upstream/four_versions
TRIAL=0bb827ebfe69a32504b74488730f80a21495eccf
mkdir "$NEW"
bash "$D/build.sh" "$NEW/baselines" "$BOOST"
mkdir "$NEW/source"
git archive "$TRIAL" src tools/build.py | tar -x -C "$NEW/source"
python "$NEW/source/tools/build.py" --boost-prefix "$BOOST" --output-dir "$NEW/candidate" --jobs 2
W="$NEW/runs"
python "$D/run_panel.py" candidate "$NEW/candidate/vina" "$W/fixed" --jobs 3
for phase in pilot matched; do
  for arm in official parent stable candidate; do
    case "$arm" in
      stable) exe="$NEW/baselines/candidate/vina" ;;
      candidate) exe="$NEW/candidate/vina" ;;
      *) exe="$NEW/baselines/$arm/vina" ;;
    esac
    python "$P/panel.py" "$phase" "$W" "$arm" "$exe"
  done
  if test "$phase" = pilot; then python "$P/panel.py" freeze "$W"; fi
done
python "$P/panel.py" collect "$W" "$NEW/results"
python "$P/identity.py" "$NEW/baselines/official.git" "$TRIAL" "$NEW/results/builds.json" \
  --official "$NEW/baselines/official/vina" --parent "$NEW/baselines/parent/vina" \
  --stable "$NEW/baselines/candidate/vina" --candidate "$NEW/candidate/vina"
python "$P/panel.py" check "$NEW/results"
```

`check` 从原始坐标、命令与搜索日志复算固定研究候选、四臂成本试验及计数，验证原始得分、RMSD、预算、二进制身份和逐格配对统计，并只用标定 CPU 记录重建冻结预算。原有三个固定对照由 `offline.py` 独立读取。165 个新执行均成功，原始四组中最大 RMSD 数值复算差低于 1.1×10⁻¹⁴ Å。

## 逐方向分布及生产版本回放

[by_direction.csv](results/by_direction.csv) 给出每个有序配体–受体方向在相应种子上的首姿均值、中位数、Q1/Q3、P90/P95、Best@5、Oracle、Regret、空池数及相对三份基线的改善/持平/退化数。每个种子的精确读数在各实验的 cells.csv；单方向仅一至两枚种子时，分位数只描述这几次实际运行。表由 `python benchmarks/search_trial_panel/case_summary.py` 重建。

生产源码在 `d37e5e4743f2460c89726d18b6452487ae085096` 完整恢复至轮前。用相同工具链重新编译，再用该二进制执行含默认设置的全部 65 个固定格，65/65 个 PDBQT 与上一轮最终版原始输出逐字节一致，572 个模式、首轮局优、入池前二次局优及显式精修计数全部一致。两份二进制有 20 处构建目录字面量和一份 ELF build ID 不同；归一化后全部字节比较相等。这轮单独复现用时 123.274 CPU 秒，单进程墙钟合计 124.206 秒，最高 RSS 43724 KiB。[构建与重放收据](restored/verification.json) · [命令、坐标及日志](restored/raw.tar.xz)

```bash
python "$D/run_panel.py" stable "$RESTORED" "$NEW/restored-runs" --jobs 4
python "$P/restore.py" "$STABLE" "$RESTORED" "$NEW/restored-runs" "$NEW/restored-evidence"
```
