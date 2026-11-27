# 两个新增命中的原因与代价

在 FABP4 与 EPAS1 的 16 个条件中，修复带来的两个阈值变化来自同一个 FABP4 共晶复合物 3P6E。它们分别发生在两枚种子的跨轨迹合并阶段。逐项比较四条轨迹给合并器的所有输入，两次运行的父版和修复版各有 36 个候选，每个候选的能量、顺序、每个原子的三个坐标均一致。FABP4 参考覆盖晶体两种占位和合法化学对称映射。两枚种子是同一套结构的重复测量。

- **684173**。第 28 个输入同时覆盖多个代表，修复版保留它后9 个名额中有 8 个被占用。第 34 个输入离实验坐标 1.829 Å，父版拒绝，修复版保留。它以第 9 名进入精修，精修后 1.800 Å，输出仍排第 9。池最佳由 3.889 降至 1.800 Å，但前五最佳由 3.889 升至 4.331 Å，首姿不变。这次改善体现了固定名额下的覆盖与排序之间的取舍。
- **961751**。第 31 个输入离实验坐标 4.212 Å，比两个重叠代表中的一个更差。原算法替换距离最近的较差邻居，修复版保留较优阻挡者并拒绝这次替换。第 32 个输入在两臂都留下来，精修前后 1.356→1.358 Å。按预测分数排在它前面的错误竞争者减少了一个，最终第 6 名升到第 5 名。两次运行都没有改变该正确姿势的局部优化或评分函数。

两例父版到修复版的 hunt 和第二次最小化计数逐项相同，显式精修分别少 8 次和 7 次求值。[每个输入的冻结重放](trace/merge_examples.csv)和[每个精修位点的几何、精确坐标哈希及输出名次](refinement/slots.csv)可直接复核这两个归因。名次按打印坐标逐个匹配，两臂共成功对应 282 个输出。未输出的精修终点名次留空。

| 每轨迹上限 | 在线轨迹保存 | 跨轨迹与精修前去重 | Top-1 / Top-5 / 输出池 | hunt + 第二次最小化 + 显式精修求值 |
|---:|---|---|---:|---:|
| 12000 | 父版 | 父版 | 8 / 9 / 13 | 774526 |
| 12000 | 修复 | 父版 | 8 / 9 / 13 | 774706 |
| 12000 | 父版 | 修复 | 8 / 10 / 14 | 774518 |
| 12000 | 修复 | 修复 | 8 / 10 / 14 | 774699 |
| 11400 | 父版 | 父版 | 8 / 9 / 13 | 736082 |
| 11400 | 修复 | 修复 | 8 / 10 / 14 | 736225 |

父版在 64 条轨迹中有 26 条曾违反 1 Å 保存半径，共 954 次更新后池不满足分离条件。在 16 次跨轨迹合并中，3 次出现低于 2 Å 的代表，共涉及 50 次更新后的状态。修复版全程无上述违反。两臂送入显式精修的池均满足 1 Å。EPAS1 的 3H82 共晶（684173）提供了另一类反例：父版的第 9、3 号输入原先相距 2.097514 Å，在显式精修后缩至 0.920646 Å。修复版同次运行的最小输出间距为 1.203452 Å。线上去重约束只作用于代表更新时刻，之后的优化可改变间距。[逐池距离](refinement/spacing.csv)

修复轨迹内保存后，默认逻辑会在空名额上继续第二次最小化。12000 档第二次最小化由 3239 次、127580 次求值升至 3469 次、136187 次求值，hunt 求值由 641670 降至 633127。11400 档修复版相对 12000 档父版，实际总求值少 4.95%，姿势阈值结果保留。在 CPU 5 上交替顺序、逐条件核对输出字节的三次无插桩计时中，墙钟合计 9.35075→9.54366 秒，CPU 合计 9.32775→9.50179 秒。三轮墙钟各自小幅增加。因而实际求值和总时间给出了不同的成本判断；此面板未证明完整修复兼有总墙钟优势。[全部配对计时](cost/paired_cost.csv)

## 复现

下列命令从仓库根目录运行。使用独立结构面板相同的 Python、编译器与依赖。`$B` 为用 `benchmarks/archive_overlap/build_research.sh` 构建的产物目录，`$W` 是可写的空工作目录，`$CPU` 是可绑定的 CPU 编号。

```bash
T=benchmarks/archive_overlap/transfer
M=$T/mechanism
mkdir -p "$W/panel"
tar -xf "$T/inputs.tar.xz" -C "$W/panel"
tar -xf "$T/outputs.tar.xz" -C "$W/panel"
python "$M/transfer_stages.py" "$W/panel" "$W/trace" trace --cpu "$CPU" --parent "$B/observed-parent/vina" --candidate "$B/observed-fixed/vina"
python "$M/transfer_stages.py" "$W/panel" "$W/reduced" reduced --cpu "$CPU" --parent "$B/parent/vina" --candidate "$B/fixed/vina"
python "$M/transfer_stages.py" "$W/panel" "$W/factorial" factorial --cpu "$CPU" --save "$B/save/vina" --merge "$B/merge/vina"
python "$M/budget_time.py" "$W/panel" "$W/reduced" "$W/cost" "$B/parent/vina" "$B/fixed/vina" --cpu "$CPU"
for arm in parent fixed; do
    python benchmarks/accepted_retention/core/build_stage.py "$B/$arm/src-tree" "$B/$arm" "$W/stage-$arm"
done
python "$M/transfer_stages.py" "$W/panel" "$W/refinement" refinement --cpu "$CPU" --parent "$W/stage-parent/vina" --candidate "$W/stage-fixed/vina"
```

使用非系统 Boost 时，可向两次 `build_stage.py` 构建及 `build_research.sh` 传 `--boost-prefix PATH`。原始观察器计时仅用于验证，表中时间均来自单独的未插桩构建。无需重新对接也可从原始坐标复核关键几何与事件。

```bash
tar -xf "$M/runs.tar.xz" -C "$W"
python "$M/check_spacing.py" "$W/refinement/runs" > "$W/spacing.csv"
python "$M/merge_examples.py" "$W/trace/runs" "$W/panel" > "$W/merge_examples.csv"
diff "$M/refinement/spacing.csv" "$W/spacing.csv"
diff "$M/trace/merge_examples.csv" "$W/merge_examples.csv"
```
