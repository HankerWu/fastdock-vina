# 口袋质心的初始放置

`--pocket_start_radius` 在单配体、刚性受体的 Vina 或 Vinardo 全局对接中选择可选初始放置。半径单位为 Å。省略时沿用原初始放置与随机数流。启用时以实际搜索盒中心为球心，按 1 开始的轨迹编号，偶数轨迹将完成原始随机化的重原子质心平移到球内均匀采样点。每条轨迹的放置随机流独立于原搜索随机流。未锚定轨迹保留原初始坐标、取向、扭转和后续搜索。锚定轨迹保留原取向和扭转。整个搜索盒仍用于此后的随机移动、打分和精修。

半径须为有限正数且不超过实际盒最短边的一半。用实验供体或其他独立证据给出的受体和盒，可避免查询配体实验坐标参与预测。

```bash
build/pocket-start/vina --receptor receptor.pdbqt --config donor_box.txt \
  --ligand ligand.pdbqt --cpu 1 --seed 42 --exhaustiveness 8 \
  --pocket_start_radius 3 --search_counters --out anchored.pdbqt
```

`--search_counters` 对单配体刚性搜索另外输出 `Pocket start` 行，记录轨迹编号、是否锚定和随机化后的初始配体质心。原有搜索和精修计数另列。初始化对照时两臂均应保存计数。

```bash
python3 benchmarks/check_pocket_start.py build/pocket-start/vina build/trusted-ligand-bonds-serial/vina
```
