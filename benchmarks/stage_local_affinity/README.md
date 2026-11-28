# 固定 Ki/Kd 面板与生产源复算

[事前登记](../stage_local_attribution/PROTOCOL.md)固定五种人源靶标、九个测定–受体分层、七份受体和两枚随机种子（89513、89527）。发布的 Ki、Kd 必须在同一靶标、同一测定和同一受体内比较。每个配体的预测为两次首姿分数的中位数，分数越低预测结合越强。≥3 倍描述实测浓度分离；论文没有为全部分子对给出独立重复及置信区间。不同分子对可共享化合物。默认设置与相近 CPU 比较另用六个固定 BACE1 配体、4DJW、种子 89513。共同 CPU 预算为六个配体合计 14.5 秒，验收区间 12.325–16.675 秒，沿用原来的标定上限和最多两次确认。

## 只修改首轮局优后的盲预测

候选 0bcae6ef09857f97656291bb4821d20b3b619bf5 与轮前采用相同评分与端点。下表为 12k 次/轨迹、每格两条轨迹，两枚种子分别运行所得的结果。官方 v1.2.7 和真实父版的 140 份发布文件逐字节一致。187 个实测相差至少三倍的分子对按测定与受体计数。与轮前相比，42 对由错转对，36 对由对转错，109 对持平。五个分层增加正确数，四个减少。相对官方与父版都是 34 对改善、31 对恶化、122 对持平，四个分层改善、四个退化、一个持平。

| 靶点 / 测定 / 受体 | 官方、父版 ρ / τ-b / 正确数 | 轮前 ρ / τ-b / 正确数 | 候选 ρ / τ-b / 正确数 |
|---|---|---|---|
| ALK / Davis-ALK / 7BTT | 0.810 / 0.571 / 17/19 | 0.357 / 0.214 / 12/19 | 0.476 / 0.357 / 13/19 |
| BACE1 / CHEMBL2015090 / 4DJW | -0.829 / -0.600 / 2/14 | -0.886 / -0.733 / 1/14 | -0.029 / 0.067 / 7/14 |
| BACE1 / CHEMBL2015090 / 4GID | -0.486 / -0.333 / 4/14 | 0.771 / 0.600 / 11/14 | -0.257 / -0.200 / 5/14 |
| BACE1 / CHEMBL2155063 / 4DJW | -0.762 / -0.571 / 2/22 | -0.595 / -0.357 / 5/22 | -0.881 / -0.714 / 1/22 |
| BACE1 / CHEMBL2155063 / 4GID | 0.286 / 0.143 / 14/22 | -0.214 / 0.000 / 10/22 | 0.452 / 0.286 / 14/22 |
| MAPK14 / Davis-p38-alpha / 1KV2 | 0.267 / 0.278 / 19/29 | 0.150 / 0.167 / 17/29 | 0.617 / 0.500 / 22/29 |
| MAPK14 / Davis-p38-alpha / 3GCP | -0.033 / 0.056 / 15/29 | -0.250 / -0.222 / 11/29 | -0.050 / 0.000 / 14/29 |
| RET / Davis-RET / 7DUA | -0.238 / -0.071 / 8/19 | 0.071 / 0.000 / 11/19 | 0.190 / 0.143 / 10/19 |
| TNK2 / Davis-TNK2 / 7KP6 | 0.381 / 0.286 / 12/19 | 0.381 / 0.286 / 12/19 | -0.095 / -0.071 / 10/19 |

六个 BACE1 化合物（CHEMBL2015090）在 4DJW 的正确数从 1 升至 7，在 4GID 从 11 降至 5。先前轮前同一对实测在两种受体上方向翻转 12/14 次，候选翻转 6/14 次；候选有三对在两种受体都正确、五对都错。CHEMBL2011965 的 Ki 比 CHEMBL2011982 小 285.7 倍，候选却在 4DJW/4GID 分别给出弱者减强者 −0.3155/−1.816 kcal/mol。换受体后少翻转不能直接解释为更接近实测。另一 BACE1 八配体系列在 4DJW 上，候选的预测结合强度与重原子数、TORSDOF 的 ρ 分别为 +0.740、+0.678，实测 pKi 对应为 −0.667、−0.678；转到 4GID，预测相关性为 −0.145、−0.365。每个分层内部的形式电荷均为常量，无法估计电荷偏差。

原始[配体与测定分层](results/assays.csv)、[逐配体与化学骨架](results/ligands.csv)、[分骨架秩相关](results/scaffolds.csv)、[分子对方向](results/pairs.csv)、[相对三种对照的逐对变化](results/pair_changes.csv)、[分层配对计数](results/pair_summary.csv)和[BACE1 的逐对双受体预测](results/bace1_receptor_pairs.csv)可逐行核对。每臂 60 个测定–受体–骨架组合，其中十组有足够的非恒定数据计算秩相关。

## 默认、实际 CPU 和复算开销

六个 BACE1 配体在原生科学默认下四组均正确 1/14。事前冻结的共同 CPU 目标为 14.5 秒，六个配体在四臂的实际确认分别为 13.975、16.627、14.077、14.373 CPU 秒，4/4 组在区间内。正确方向分别是 2/14、3/14、3/14、1/14；候选 ρ 为 −0.943、τ-b 为 −0.867。候选三个标定点（3k、96k、384k）实测 2.666、5.458、15.243 秒，按插值得 cap 362100。早期以稳定版的低预算时长为目标的一次历史匹配原来就失败，这轮继续保留其实际开销及预测，不能纳入 14.5 秒确认。[逐次标定及确认](results/common_cost.json)

候选经正常 CLI 新执行 182 次，均返回有效发布文件。固定预算的 140 次耗费 68.492 CPU 秒、73.075 作业墙钟秒、峰值 34,272 KiB；首轮、二次局优、显式精修分别进行 3,025,026、339,812、27,211 次联合能量与梯度调用。候选全部运行合计 226.834 CPU / 258.670 作业墙钟秒，含低预算、旧匹配、默认、标定与最后确认，每份确认只计一次。抽取并验证共享制备输入耗用 0.140 CPU / 0.140 墙钟秒。原始输出、stderr/stdout 和逐进程成本见 [raw_candidate.tar.xz](results/raw_candidate.tar.xz)，四版本全成本及旧官方/父版被动计数对照见 [physical_costs.json](results/physical_costs.json)，分阶段统计见 [costs.csv](results/costs.csv)。原生搜索进程的读取、建图、重排和发布计入总 CPU；各个阶段的独立 CPU 未由这些计数器测量。

从候选原始产物和历史三个对照，在全新目录重新核验 47 个测定记录、728 份唯一 PDBQT、752 个运行标签、1,448 条分子对和 68 组排名摘要；四份新增分子对、受体切换和CPU表也逐字节复现。[独立验证](results/verification.json) · [验证的 CPU、墙钟和内存](independent_cost.json)

## 生产源的正交复算

恢复提交 d37e5e4 与轮前 666f49d、历史对照 5847b1b 的 src 树均为 f22a5de502749b4743ab0ce7fc2c9cbcbfdf7948。用恢复提交编译出的真实二进制重新启动 140 次。全部 1,205 个模式的 PDBQT 逐字节相同，逐格局优计数亦全部一致。此次首轮搜索、入池前二次局优、显式精修分别耗费 3,064,849、300,631、29,085 个联合能量/梯度调用。共耗用 57.130 CPU 秒，单作业墙钟求和 57.335 秒，逐进程峰值 RSS 34,316 KiB。参考与恢复二进制的 SHA256、原始命令、日志、所有坐标和收据分别见 [restored.json](restored/restored.json) 与 [restored.tar.xz](restored/restored.tar.xz)。两个编译路径的 ELF 身份对照见[历史规范化记录](../search_trial_panel/restored/verification.json)。

历史四臂从原始 ChEMBL Ki 记录和 Davis Kd 工作表，以及每次实际发布的 PDBQT 重新核对：10 个测定来源哈希、47 行测定–化合物记录、728 个独立首姿、752 个记录标签、68 个分层摘要及 1,448 条分子对。独立复算逐文件一致，耗时 2.414 CPU / 2.424 墙钟秒，峰值 RSS 92,524 KiB。历史主角度候选已退出生产。历史归档中带 candidate 字样的目录名由各轮次定义；本研究导入的 stable 只取 5847b1b 的完整参考池。所有输入来源与制备哈希见[历史面板](../cumulative/affinity/README.md)。该面板曾测得从缓存原件重做制备为 43.55 CPU / 79.18 墙钟秒；两份 MAPK14 用既有受体 PDBQT，其最初制备时间未知。

## 复算命令

构建身份见[四个真实引擎](../stage_local_structure/fixed/builds.json)。在仓库根目录运行。W 为新的临时目录。归档重建只需要 Python、NumPy、SciPy 和 xlrd，真实重复运行还需相应源码的 Vina 可执行文件。

```bash
S=benchmarks/search_trial_affinity/study.py
D=benchmarks/stage_local_affinity
W=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 python "$S" replay "$W/previous"
OPENBLAS_NUM_THREADS=1 python "$S" setup "$W/base"
mkdir "$W/restored"
tar -xJf "$D/restored/restored.tar.xz" -C "$W/restored"
python "$D/restored.py" "$W/restored" "$W/base" "$W/check"
cmp "$D/restored/restored.json" "$W/check/restored.json"
```

独立真实运行生产对照时，将上述解包替换为以下命令。

```bash
OPENBLAS_NUM_THREADS=1 python benchmarks/cumulative/affinity/run.py "$W/restored" \
    --prepare --engine stable --binary "$RESTORED" --stage nominal --jobs 2
```


对候选的实际四臂评价，以下命令沿用冻结的历史输入和三个基线运行（上例的 W/base）；CANDIDATE 指向 0bcae6e 的可执行文件。

```bash
OPENBLAS_NUM_THREADS=1 python "$S" run "$W/base" --binary "$CANDIDATE" --jobs 3
OPENBLAS_NUM_THREADS=1 python "$S" collect "$W/base" --artifacts "$W/new-results"
python "$D/paired.py" "$W/new-results"
python "$D/costs.py" "$W/base" "$W/new-results/physical_costs.json"
```

离线复核新研究的原始分值、所有阈值、相关性及分子对方向，以下命令每一步都从归档数据重建。

```bash
OPENBLAS_NUM_THREADS=1 python "$S" replay "$W/new-audit" --artifacts "$D/results"
python "$D/paired.py" "$W/new-audit/recomputed"
python "$D/costs.py" "$W/new-audit" "$W/new-audit/recomputed/physical_costs.json"
for file in pair_changes.csv pair_summary.csv bace1_receptor_pairs.csv physical_costs.json; do
    cmp "$D/results/$file" "$W/new-audit/recomputed/$file"
done
```

5U2E 已揭盲同图控制沿用[前轮定义的两个固定坐标和评分参考](../search_trial_nov28/POSTHOC_5U2E.md)。四次评分从独立进程取得，生产源和候选的十项打印能量逐项相同。1.211 Å 近姿相对 7.190 Å 远姿仍高 +0.314 kcal/mol。四次总计 0.096 CPU 秒；原始命令和输出在 [fixed_5u2e](fixed_5u2e/receipt.json)。按[原研究复算步骤](../search_trial_nov28/README.md)准备 FIXTURES 与 RUNS 两个目录后运行下列命令。

```bash
python benchmarks/search_trial_nov28/fixed_5u2e.py "$FIXTURES" "$RUNS" \
    "$RESTORED" "$CANDIDATE" "$W/5u2e"
```
