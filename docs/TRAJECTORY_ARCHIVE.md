# Monte Carlo 轨迹存档

在直接提供刚性受体和单个配体的普通 Vina 对接中，加 `--pose_selection late --trajectory_archive accepted`
可将 Monte Carlo 接受的部分姿势纳入原生显式受体精修。做支持度诊断时用 `support` 对同一最终候选池重新排序。
省略 `--trajectory_archive` 时采用原来的 Late 路径。存档模式使用原生优化器和评分权重，支持普通闭合环，
拒绝配体的虚拟原子、宏环闭合原子类型及与模板、极性几何、宽高斯缩放或 Armijo 组合。

```bash
./build/refine-pool/vina --receptor examples/basic/1iep_receptor.pdbqt \
  --ligand examples/basic/1iep_ligand.pdbqt --config examples/basic/1iep_receptor.box.txt \
  --pose_selection late --trajectory_archive accepted --refine_candidates 18 \
  --exhaustiveness 6 --max_evals 1100 --num_modes 6 --cpu 1 --seed 5 \
  --out 1iep_archive_accepted.pdbqt
```

每条独立随机轨迹额外保存最多 `--num_modes` 个已接受的姿势，取其初次局部优化后的能量、构象与化学重原子坐标。
采样节点位于原有第二次优化前，沿用同一随机序列和原搜索评分次数。存档用固定原子序号计算重原子 RMSD，
在低于 `--min_rmsd` 或坐标完全重合时保留较低能量。存档与历史轨迹的原有候选汇合后进入 Late 预评估，
完全相同的构象优先采用历史候选。末端精修上限可用 `--refine_candidates` 调整。 [Late 参数](../README.md#延后筛选候选姿势)

- `accepted` 按完整精度的最终 Vina 分数排序，同分时沿用 Late 精修筛选顺序。
- `support` 为每个姿势统计支持它的独立轨迹数。某轨迹可用它保存的任何一个历史或存档姿势投一票，
  位移按同一重原子序号计算，在共同受体坐标系中的 RMSD 必须小于或等于 2.0 Å。
  落在池内最佳原生分加 1.0 kcal/mol 以内的所有姿势按票数降序排列，同票保持原生名次。
  其余姿势随后按原生名次排列。

两种模式的 `--energy_range` 均基于池内最佳原生 Vina 分数。比较这两种排序时应固定末端精修预算，
可直接用 `support` 输出中的原生名次恢复 `accepted` 的排序。

每个 PDBQT 模型开头包含以下字段。浮点数使用可往返恢复的精度。

```text
REMARK VINA ARCHIVE RESULT: <native_score> <rmsd_lb> <rmsd_ub>
REMARK VINA NATIVE: <native_score> <native_rank>
REMARK VINA ARCHIVE SUPPORT: <votes> <eligible:0|1> <computed:0|1>
REMARK VINA ARCHIVE SOURCE: <historical|accepted>
REMARK VINA ARCHIVE MODE: <accepted|support>
```

原生名次从 1 开始。来源是该模型进入显式精修时的候选来源。
`eligible` 表示是否落在原生 +1.0 kcal/mol 窗口。仅 `support` 计算票数并将 `computed` 写为 1。
在 `accepted` 中，这两个计数字段为 0，`eligible` 仍按同一个窗口填写。

常规日志同时报告 `Late pose selection` 和 `Trajectory archive`。
`native` 与 `extra` 分别统计各轨迹原有与额外存档的候选。 `extra_prescore` 和 `extra_refine` 分别
统计实际参加显式预评估和精修的存档来源候选。 `native_search_evals` 和 `native_search_refinements` 统计
原有 Monte Carlo 路径的评分评估和接受后第二次优化。 `retained` 是输出能量窗前的姿势数量。

```bash
python3 benchmarks/check_trajectory_archive.py build/refine-pool/vina /path/to/pre-archive/vina
```

数值回归的入口在 [check_trajectory_archive_numeric.cpp](../benchmarks/check_trajectory_archive_numeric.cpp)。
