# 已接受却未保存的近原位构象

固定六个人源靶点、22 个方向、3 枚种子，全量追踪 66 条原生单轨迹。RMSD 按化学对称性校正，在受体坐标系中计算，以 ≤2 Å 为命中。输入来自 [六家族结构面板](../../cubic_grid/human_panel/README.md)，本轮的策略、轨迹数、名额和预算以[冻结协议](PROTOCOL.md)为准。原生常规、观察器、单名额实验各使用同一配体、盒、受体和随机流。原生与观察器在全部 66 个条件返回相同的 PDBQT 文件字节、整套搜索和精修计数。

| 靶点 | 条件 | 有 ≤2 Å 接受态 | 至少一次送第二次优化时已 ≤2 Å | 原生输出池命中 | 有未送第二次优化的正确姿势且原生池未命中 |
|---|---:|---:|---:|---:|---:|
| DYRK1A | 12 | 6 | 6 | 6 | 0 |
| EPAS1 | 12 | 5 | 5 | 5 | 0 |
| FABP4 | 12 | 4 | 3 | 3 | 1 |
| FKBP1A | 12 | 2 | 0 | 0 | 2 |
| JAK2 | 6 | 0 | 0 | 0 | 0 |
| RBP4 | 12 | 2 | 2 | 2 | 0 |
| 合计 | 66 | 19 | 16 | 16 | 3 |

44/66 在所有被记录的 hunt 局优终点（包括被拒绝者）及最终返回池中均无命中。各列描述阶段是否*曾经*命中；某个近原位输入也可能在局优中被带出 2 Å，远处输入也可能变好。表格没有观察 BFGS 线搜索每个中间点。所有条件、包括失败，见[逐运行记录](results/opportunities.csv)。

## 只替换原有末名的反事实

利用实验坐标挑出未送第二次优化的最近接受态，将其内部平移、四元数与全部扭转逐项恢复，替换原生最终去重后即将显式精修的末名。这额外用到了答案，只作机会诊断。搜索和原生候选合并的计数逐行相同，显式输入的 K 相同，实际精修评估分别多 4、15、3 次。

| 配体→受体，种子均为 75431 | 原生池最佳 Å | 漏存步 | 该步 RMSD→精修后 Å | 最终名次 |
|---|---:|---:|---:|---:|
| FABP4 2NNQ→3FR2 | 2.874 | 193 | 1.606→1.657 | 2 |
| FKBP1A 1FKG→1FKG | 2.850 | 65 | 1.821→1.803 | 4 |
| FKBP1A 1FKH→1FKG | 6.183 | 265 | 1.592→1.611 | 1 |

[原始 PDBQT、内部构象、日志与计数](oracle)保留了所有三项干预。在这三个明确的机会中，完美识别可以把返回池命中增加三次，但原生末端评分只让其中一个成为首姿。

## 不使用晶体的选择

最低软化能量选取的新颖接受态在 38/66 个运行中替换一个名额。它与移除最后一道能量门槛的消融在这 66 个条件产生完全相同的输出。父版本 Top-1/Top-5/Oracle 为 14/16/16，选择后为 14/17/17。联合评估 814,825→813,806，端到端耗时合计 12.113→12.051 秒。其余 28 个条件无替换时坐标逐字节不变。盲选救回的 FABP4 输出池最佳 RMSD 为 1.690 Å，首姿仍为 2.874 Å。

同一 FABP4 种子的额外成本对照覆盖所有改变成功状态的条件（本组只有这一项）。两臂均跑 11,900 和 12,000 两档，每档七次、交错顺序、同一个逻辑 CPU。父版本 12,000 上限实际 12,193 次联合评估，中位耗时 187.37 ms。新法 11,900 上限仍然命中，实际 12,089 次、188.96 ms；新法 12,000 上限 12,204 次、187.81 ms。坐标及联合计数在同档的各次重放完全相同。这是单个病例在略少评估及近似等时间下的候选覆盖收益。百分之一以内的耗时差不支持速度排名。[逐次成本与原始文件散列](results/cost_control.csv)

在三个反事实病例中，观察器读取了每个漏存接受态的真实 hunt 能量和与当时已存池的最近重原子距离。三例的最佳近原位接受态距当时已存姿势均超过 2 Å。FABP4 最终选择了同一近原位盆地的第 198 步。FKBP1A 的最低能量规则分别选中 RMSD 4.976、6.746 Å 的姿势；它们比上表近原位构象的软化能量仅低约 0.77、0.32 个评分单位。[409 个接受态的选择记录](results/choices.csv)

## 复现

Python 3.12、NumPy、RDKit、Gemmi 及依赖见原结构面板。以下在仓库根目录运行。B 指向 pair_gradient 的七臂构建输出；普通 parent 对应 31d073f 在默认构建下的原生调用链。父版本与这批追踪所用的源构建在这 66 项输出逐字节一致。共享配体及受体制备直接重用归档。程序按当前机器选择可用的 CPU；下面以 2 为例。

```bash
D=benchmarks/accepted_retention/coverage
W=$(mktemp -d)
python "$D/audit.py" "$W/audit" "$D/results" "$W/recomputed"
cmp "$D/results/coverage.csv" "$W/recomputed/coverage.csv"
cmp "$D/results/opportunities.csv" "$W/recomputed/opportunities.csv"

T=benchmarks/pair_gradient/stages
python "$T/build_trace.py" "$B" "$W/tracer" --arm parent
python "$T/build_trace.py" "$B" "$W/oracle-build" --arm parent --overlay "$D/oracle.patch"
python "$D/probe.py" "$W/raw" "$W/fresh" "$B/parent/build/vina" "$W/tracer/parent/vina" --cpu 2
python "$D/audit.py" "$W/replayed" "$W/fresh" "$W/fresh"
python "$D/intervene.py" "$W/replayed" "$W/fresh" "$W/oracle" "$W/oracle-build/parent/vina" --cpu 2
```

所有轨迹的坐标（整数千分之一 Å，PDBQT 精度）、阶段、全精度对称 RMSD、原生输出坐标与原生日志按靶点压缩归档。audit.py 逐项复算坐标、允许的舍入误差和成功判定。导出时取整的误差上界为 √3 × 0.0005 Å，全部归档分类复算一致。晚期精修中的名次按原生输出。各条件只对同一受体、配体和种子作配对，不把轨迹事件当成独立生物样本。

盲选程序 blind.py 以 work、results、parent、gated、ungated 和 --cpu 为参数，重放全部 198 项并验证 MC、最终名额及未触发时的字节等价。cost_control.py 以 work、results、parent、gated 和 --cpu 为参数。choice.py 以 work、results、observer、gated 为参数；build_choice.py 使用编译核心 gated 时的源码和对象目录构建诊断观察器。原始归档 blind.tar.xz 含全部盲选日志及所有不同于同条件原生输出的 PDBQT；与原生相同者由散列引用。
