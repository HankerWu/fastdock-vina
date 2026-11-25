# 分离搜索轨迹与候选网格精修

`--separate_hunt_chain` 是 FastDock Vina Pro 的可选搜索方式。它用于带显式刚性受体的单配体 Vina 全局对接，也支持逐配体批处理。可以与原生和 Late 候选选择组合。比较时用同一个构建，另一臂省略该选项。

```bash
build/separate-hunt-chain/vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --config box.txt --scoring vina --cpu 1 --seed 42 --exhaustiveness 8 \
  --max_evals 500000 --num_modes 9 --min_rmsd 1 --energy_range 100 \
  --separate_hunt_chain --search_counters --out separate.pdbqt
```

每步先对提议做网格局优。开启后，上一接受点与新局优点都用软化目标比较，候选库装满后的保存触发也比较该目标的历史最优值。普通小角旋转从上一接受构象计算配体的转动半径，首步使用该条轨迹的随机初始构象。候选的第二次网格局优作用于副本。最终显式受体精修使用原有势函数、盒边界和候选上限。

对正能量输入 `e`，网格中的正项采用 `C_v(e) = v·e/(v+e)`，其余项及盒外罚沿用原定义。第一次网格局优的 `v = 10`，进入候选库前的第二次为 `v = 1000`，参数和 `e` 同为评分能量单位。未经分离时，第二次局优会同时更新轨迹坐标和接受参考值。这样后续判断可用 `C_1000` 得到的历史值比较 `C_10` 得到的新值，分差还会受到两次局部优化落点不同的影响。分离后，仍可使用软化景观探索并对保存姿势按原能量精修。这个算子用于迭代局部搜索，接受率无热力学平衡概率含义。

`--search_counters` 输出逐轨迹的 `Hunt chain:`。其中 `accepted` 包含首次自动接受，`mixed_accept_checks` 统计使用第二次网格局优值作历史参考的 Metropolis 判断次数，包括连续拒绝时重复使用该值的判断。启用分离时该值为零。梯度计数、第二次网格局优和末端显式精修成本见原有 `Search trajectory:` 与 `Search counters:`。配对比较应汇报实计预算、逐查询姿势和候选覆盖。

[Vina 原论文](https://doi.org/10.1002/jcc.21334)
