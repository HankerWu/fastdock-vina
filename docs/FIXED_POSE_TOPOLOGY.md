# 固定姿势与配体化学拓扑

Vina 从准备 PDBQT 的坐标和可动矩阵推断键，随后生成 XS 类型和内部评分配对。直接加载保存姿势时会重新执行推断。为了在已准备配体的拓扑上评价另一个姿势，可在同一配体上只替换实验坐标或保存坐标。

```bash
build/trusted-ligand-bonds-serial/vina --receptor receptor.pdbqt --config box.txt \
  --ligand prepared.pdbqt --pose_coordinates model.pdbqt --scoring vina \
  --score_only --unbound_energy 0
```

`--pose_coordinates` 接受单个裸配体 PDBQT 或只包裹一个完整 `MODEL … ENDMDL` 的 PDBQT。两份文件的非坐标 atom 识别字段、显式电荷和类型、分支结构、`TORSDOF` 及 Meeko 化学映射备注应一致。坐标须有限，原子 serial 须为唯一的正整数。旧评分备注可保留在姿势中。该选项限单配体、显式刚性受体的 Vina 或 Vinardo 固定评分。可使用相应构建本来支持的 Vina `--net_polar` 分项或普通评分。普通评分的 `--unbound_energy 0` 便于保存相同参照尺度。准备图推断有误时，坐标移植仍沿用该准备图。

## 使用已核实的配体键

```bash
build/trusted-ligand-bonds-serial/vina --receptor receptor.pdbqt --config box.txt \
  --ligand prepared.pdbqt --ligand_bonds ligand.bonds --scoring vina \
  --cpu 1 --seed 42 --out poses.pdbqt
```

`--ligand_bonds` 使用每行两个以空白分隔的正整数 PDBQT atom serial。每条化学键只列一次，端点顺序任意。允许空行和 `#` 注释。表须包含当前单个 PDBQT 配体所有实际原子之间的化学键，包括 PDBQT 显式氢原子的共价键。应从经核验的制备化学图和 atom serial 映射派生，省略没有写入 PDBQT 的隐式氢。引擎核对 serial、重复、连通性、端点坐标、显式氢的邻接数、原有 BRANCH 边及可动矩阵。表作为化学连接的来源。提供方负责保证同一刚性片段内的化学边完整，尤其是环闭合边。

该路径可用于 Vina 或 Vinardo 的单配体、刚性受体搜索、局部精修及固定评分，并可和 `--pose_coordinates` 配合。引擎沿用匹配的原 PDBQT 扭转树，用明确邻接生成 XS 类型和内部评分配对。宏环虚拟点和水探针、柔性受体以及预计算地图超出此路径。默认 CLI 沿用原来的几何键推断。受体沿用既有推断。

[9FA3 数值复验](../benchmarks/check_fixed_pose_topology.cpp) 链接相同引擎构建的 `coords`、`model`、`parse_pdbqt`、`quaternion`、`random` 和 `utils` 对象及 Boost，输入原准备 9FA3、保存的 9FA3/43 MODEL 2 裸坐标和正确的配体键表。它检查坐标移植保持原 41 对，显式图恢复全部三条环丙基键，且解析不同坐标和移植坐标两条路径均保留相同的 40 对。
