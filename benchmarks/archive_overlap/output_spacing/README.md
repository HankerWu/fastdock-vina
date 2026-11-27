# 显式精修后输出姿势的最小距离

在线候选去重修复后，两个原生 Gal3 文件仍然违反默认 `--min_rmsd 1` 的输出约束。3T1L 对化合物 10 的两次显式精修从相距 2.767203 Å 的第 4、9 个输入出发，落点只相距 0.509425 Å。化合物 26 的第 3、6 个输入从 2.131005 Å 缩至 0.706068 Å。两次观察器重放均复现修复前生产二进制的完整输出字节。这些近邻由精修形成。旧 Native 在网格阶段去重，历史 Late 才会在显式精修后再去重。CLI 将 `min_rmsd` 定义为输出姿势间的最小 RMSD。[全精度复算](full_precision.csv)

现在 Native 在显式精修及最终排名完成后，依照原有顺序保留距离合格的姿势。距离使用原生同序非氢原子、共同受体坐标系，严格小于阈值才抑制后一个姿势。普通原生能量排序在最终池中直接过滤一次。原生能量窗口与自定义排名可能有不同顺序，因而对这些模式先按调用者的能窗判断资格，再在合格姿势中去重，最后满足请求的输出数量。这样未通过能窗的高排名姿势不会抑制合格的后续姿势。RMSD 表、坐标 API、能量 API 和 PDBQT 共用排序筛选规则。控制台展示具有全部能量的排名，写文件继续遵守调用者提供的能窗。默认九个名额的下标存于栈内，既有 Late 重筛和无显式精修路径无需再比较坐标。

## 实际影响

离线检查本轮研究中来自原生产修复版的 440 个不同配置记录，包含结构开发和迁移、额外结构挑战、Gal3 单/四轨迹及降档预算，共 2,861 个已输出姿势。依据打印坐标预测会删除 4 个，3 个输出文件受影响。没有任何决定落在坐标保留三位小数造成的最坏误差范围（±√3×0.001 Å）内。3T1L/10 移除原第 6 名，3T1L/26 移除原第 4 名。黑猩猩雄激素受体 1T7T（438013）保留原第 1、2、4、6 名，从 6 个降到 4 个；移除的两个姿势 RMSD 均超过 6.6 Å，首姿仍为 0.641 Å。各条件不是统计独立样本。[全部条件](archived.csv)

对四轨迹 94 个 Gal3 条件重新运行新旧两个 CLI，逐一核对每个位姿主体的原始文本（包括所有原子、分数、RMSD、备注和扭转树），与根据原文件独立得到的顺序和子集完全一致。全部首姿、hunt/第二次最小化/显式精修及各自的求值数均相同。共输出 824→822 个模型。独立的 [full_assay.py](full_assay.py) 复算了全部 [282 个单轨迹和四轨迹 Gal3 条件](release.csv)。因此这些亲和力结果沿用原来的首姿及评分。另以不同配体的预测姿势构造技术测试用接触参考，实测自定义排名下 100、3、1、0.25 和 0 五档能窗，以及原生 0 和 2 Å 两档半径。所有输出均符合原排名约束。[边界检查](edge_cases.csv)

## 成本

两轮成本检查的条件不同，原始记录分别保存。最初带计数、打印模式表的 94 对运行，墙钟合计 16.513→17.503 秒。主计时用相同编译器、参数和从空目录构建的两版二进制，单独在 CPU 5 上对每个条件做三轮无计数、无模式表的整进程计时，轮内打乱条件并交换新旧次序。全部 564 个输出文件都核对了已验证的 SHA256。

| 94 个条件一轮 | 原版墙钟秒 | 修复版墙钟秒 | 原版 CPU 秒 | 修复版 CPU 秒 |
|---|---:|---:|---:|---:|
| 1 | 19.159 | 19.697 | 18.923 | 19.420 |
| 2 | 17.282 | 17.222 | 17.099 | 17.048 |
| 3 | 19.309 | 18.872 | 18.998 | 18.730 |
| 合计 | 55.750 | 55.791 | 55.020 | 55.198 |

合计墙钟 +0.074%，CPU +0.324%。三轮墙钟差异方向不一致。94 个真实池的内存中微基准测得 1 Å 几何筛选中位 2.38 μs/池，零半径对照 0.026 μs/池。微基准只覆盖已加载的坐标。整进程差值不能直接归因为几何比较耗时。此次改动恢复输出距离约束，尚无预测精度或吞吐量普遍提升的证据。[单次真实运行](four_trajectory.csv) · [三轮计时](timings.csv) · [筛选微基准](filter_cost.csv)

## 复现

从仓库根目录运行。`$OLD` 应为先前在线候选修复后的 40b25a5 版本，`$NEW` 为本次输出间距修复后的二进制。复测计时时，两者均使用同一版本的编译器、Boost 和 `tools/build.py` 在空目录构建。`$W` 为工作目录，`$CPU` 为可绑定的 CPU。两次精修观察使用仓库已有的 `benchmarks/accepted_retention/core/build_stage.py`，在 `$OLD` 的源树上构建。

```bash
D=benchmarks/archive_overlap/output_spacing
mkdir -p "$W/inputs" "$W/challenge" "$W/reference" "$W/evidence"
tar -xf benchmarks/pair_gradient/galectin_external/raw.tar.xz -C "$W/inputs" prepared
tar -xf benchmarks/archive_overlap/transfer/refill/inputs.tar.xz -C "$W/challenge"
tar -xf benchmarks/archive_overlap/affinity/four.tar.xz -C "$W/reference"
python "$D/compare.py" "$W/inputs/prepared" "$W/pair" "$OLD" "$NEW" --cpu "$CPU"
python "$D/timing.py" "$W/inputs/prepared" "$W/pair" "$W/time" "$OLD" "$NEW" --cpu "$CPU"
python "$D/edge_cases.py" "$W/inputs/prepared" "$W/reference/raw" "$W/challenge" "$W/edge" "$OLD" "$NEW"
tar -xf "$D/raw.tar.xz" -C "$W/evidence"
python benchmarks/archive_overlap/transfer/mechanism/check_spacing.py "$W/evidence/stage" > "$W/distances.csv"
diff "$D/full_precision.csv" "$W/distances.csv"
```

C++ 回归命令见仓库 `core-regression.yml` 的 “Test final pose spacing across energy windows and output APIs” 步骤。它包含自定义排名的能窗优先、等分、精确距离边界、零名额、超过九个名额，以及同一真实 3T1L 对接在 20 个数量和能窗组合下的三类输出一致性检查。`filter_cost.cpp` 可以和同版 `coords.o` 一起编译，输入为解包后的 94 个 `*_candidate_3000.pdbqt`。

BFGS 回溯修复改变了 Gal3 3T1L 与化合物 10 的搜索轨迹。同一种子和 3,000 次上限在父版返回八个姿势，最近间距 2.4508 Å；修复版返回九个，最近间距 2.6302 Å。两种输出都满足 1 Å 的最小间距要求。回归现从全精度坐标独立计算实际发布的每一对重原子间距，并在 25 组能窗与数量组合中核对三种接口的能量、次序、序列化坐标和候选上限。固定的小坐标池继续覆盖分数相同、恰在半径边界和能窗先于几何筛选等情况。[共享 BFGS 的真实回放](../../accepted_bfgs/README.md) · [接口回归](check.cpp)
