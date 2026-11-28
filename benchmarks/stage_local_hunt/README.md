# 只在首轮网格局优保留下坡平移

研究源码为 `0bcae6ef09857f97656291bb4821d20b3b619bf5`，生产源码在 `1c3141238511ba98e13232ffff8f2ba51efef3d1` 恢复。原始假设及晋级规则见[事前方案](../stage_local_attribution/PROTOCOL.md)。源码只对 MC 首轮局优调整 BFGS 切向量。角位移限于主角度；平移方向的合计一阶功为非正时，首试保留完整平移，否则整条切向量共同缩放。线搜索检验新方向的下降和实际 Armijo 减少量，割线更新使用实际位移。二次局优、显式精修、排名所用能量和随机提议均沿用基线。几何分布、亲和力和外检见[结构四臂](../stage_local_structure/README.md)和[同测定面板](../stage_local_affinity/README.md)。

## 历史共同提议的真实求解

[补充方案](LOCKED_PROTOCOL.md)冻结后三个病例的六个历史驱动库共 7,130 个提议全部进入测试。代码从研究提交独立导出，在原 MC 提议所用网格、hunt_cap 和源码求解入口重放。全部提议的初始能量与归档相同，首个方向导数均为负。研究候选耗尽十次回溯 4,794 次，作 153,530 次联合能量/梯度调用。轮前分别为 5,079 和 144,229，上轮共同缩放角度与平移为 2,868 和 196,297。在 FABP4 的 2,691 个修复后驱动提议上，新策略有 2,690 次耗尽，虽然全库最佳 RMSD 从轮前 6.824 降至 5.559 Å，仍高于共同缩放的 3.917 Å。合计相对轮前，返回能量降低超过 1e−7 的有 1,971 个，升高的有 1,112 个；返回 RMSD 以 0.001 Å 计，1,320 个改善、1,687 个退化。其余在相应容差内。每条轨迹的提议共享历史，不能视为独立配体。[各库的调用、方向和完整返回](locked/locked.csv) · [原始命令、SHA 和成本](locked/costs.json)

独立观察器在同一 7,130 个提议上重新执行，六份原始返回逐字节相同。真正走过的 52,721 个迭代方向中有 7,462 个保留下坡平移，98 个因平移上坡同时缩短全部自由度，45,161 个无需限幅。变换前后均未发现非下降方向。这能核查实现和梯度的乘积，仍需沿真实坐标扰动才能检验梯度与局部势面是否相符。[逐方向数值稽核](../stage_local_attribution/candidate_direction_audit.json)

## 揭盲后对 RBP4 与 ESR1 的逐阶段核查

[揭盲后登记](POSTHOC_COVERAGE.md)选择两个固定 RBP4 格和全部外检交叉格。两个版本的被动观察器通过正常 CLI 重跑十二格，每格的 PDBQT 与未插桩的盲预测逐字节一致。下表为真实首轮局优终点达到 2 Å 的个数，以及各阶段连续最小误差。全部对应受体框架的对称校正重原子 RMSD，跨受体坐标为转移参考。

| 方向、种子 | 首轮近姿数，稳定→候选 | 首轮最小 Å | 发布池最小 Å |
|---|---:|---:|---:|
| 1RBP→3FMZ，85237 | 54→28 | 1.232→1.497 | 1.614→2.896 |
| 3FMZ→1RBP，85237 | 0→0 | 8.952→2.161 | 10.824→2.373 |
| 1QKT→1X7R，85231 | 314→0 | 0.280→6.268 | 0.561→6.289 |
| 1QKT→1X7R，85237 | 0→0 | 2.578→2.618 | 2.648→2.648 |
| 1X7R→1QKT，85231 | 379→0 | 0.733→5.268 | 0.888→5.771 |
| 1X7R→1QKT，85237 | 0→0 | 9.304→11.012 | 10.769→12.944 |

在 1RBP→3FMZ 的候选中，28 个近姿 hunt 有 24 个获接受，14 个近姿出现在二次局优终点，5 个近姿入池，4 个被逐出。三条轨迹的终池里仍有 1 个 1.862 Å 的近姿。合并后送显式精修的九个输入最低为 2.897 Å。这个最后存活的近姿在跨轨迹合并至显式精修输入的环节消失。其余九个输入的成对精修结果保存在 refinement.csv，单靠这些记录无法认定某一种聚类冲突或比较顺序。反向 3FMZ→1RBP 则在首轮最近距离就有 6.791 Å 的改善，两版在该种子的首轮均未达到 2 Å。对于 ESR1，85231 的两个方向合计 693 个近姿首轮终点在候选路径上降至零。计数描述重复访问；该观察不能反推出分叉前任一真实起点在另一求解器下的结果。

[每格分母、接受、入池、逐出、精修和发布](coverage/stages.csv) · [送精修的每个输入与其返回、发布排名](coverage/refinement.csv) · [观察器哈希、逐格成本和原始输出哈希](coverage/costs.json)

## 成本与复现

锁定重算六次真实启动共 1.943 CPU 秒、1.992 作业墙钟秒，最高 84,612 KiB；十二次自由轨迹诊断为 4.338 CPU 秒、4.912 作业墙钟秒，最高 85,536 KiB。诊断计时含建图和观察器 I/O。合格样本的制备成本与正常引擎成本记录在结构面板。离线核验及数值回归另计在 [checks.json](checks.json)。首次观察器导出命令从子目录使用 Git 相对路径而失败，修正为仓库根目录后完整导出并编译；未产生任何对接输出。

从仓库根目录开始，W 为一个新目录，BOOST 为可用 Boost 前缀。以下两项离线命令会重新读取压缩的原始坐标/状态与先前盲预测，并核对每一行汇总值。

```sh
export OPENBLAS_NUM_THREADS=1
D=benchmarks/stage_local_hunt
mkdir -p "$W"
python benchmarks/bfgs_source_replay/archive.py restore benchmarks/bfgs_source_replay/raw "$W/historical"
python "$D/locked.py" "$W/historical" "$D/locked"
python "$D/coverage.py" "$D/coverage"
```

重新构建观察器并真实执行，另用干净导出的研究源码运行数值回归。

```sh
python "$D/build_observer.py" "$W/observer" "$BOOST"
python "$D/locked.py" "$W/historical" "$W/new-locked" --binary "$W/observer/build/vina"
bash "$D/build_free.sh" "$W/free" "$BOOST"
python "$D/coverage.py" "$W/new-coverage" --builds "$W/free"
mkdir "$W/source"
git archive 0bcae6ef09857f97656291bb4821d20b3b619bf5 src tools/build.py | tar -xC "$W/source"
python "$W/source/tools/build.py" --boost-prefix "$BOOST" --output-dir "$W/build" --jobs 4
for test in "$D/check_direction.cpp" benchmarks/check_armijo_numeric.cpp; do
    c++ -O3 -std=c++14 -I"$W/source/src/lib" -I"$BOOST/include" "$test" \
        "$W/build/obj/quaternion.o" "$W/build/obj/random.o" -o "$W/test"
    "$W/test"
done
```
