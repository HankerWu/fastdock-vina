# 控制原生跨轨迹合并距离

在单配体、刚性显式受体的 Vina 全局对接中，`--use_min_rmsd_for_merge` 使原生跨轨迹合并使用 `--min_rmsd`。默认路径在该阶段使用 2 Å。开关只用于原生候选筛选，也支持逐配体批处理。要求 `--min_rmsd` 为有限非负的 Å 值，命令未给该参数时取 1 Å。

```bash
build/native-merge-rmsd/vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --config box.txt --scoring vina --cpu 1 --seed 42 --exhaustiveness 8 \
  --max_evals 500000 --num_modes 9 --min_rmsd 1 --energy_range 100 \
  --use_min_rmsd_for_merge --search_counters --out requested.pdbqt
```

轨迹内保留、跨轨迹插入顺序、能量比较和精修名额沿用原有实现。同一编号原子的坐标 RMSD 小于阈值时参与合并。小于 2 Å 的阈值允许更多近邻竞争固定名额，也可能挤出较远的候选。大于 2 Å 的阈值可在精修前淘汰更远的近邻。该选项的距离作用于网格阶段的候选。精修后的分布仍取决于局部优化的落点。

对照省略 `--use_min_rmsd_for_merge`，其余参数相同。用 `--search_counters` 核对两臂 `Search trajectory:` 完全一致。 `Native merge:` 记录实际阈值和合并后候选数，`Search counters:` 记录实际进入显式精修的候选和梯度预算。指定 2 Å 时开关两臂在此合并阶段等价。
