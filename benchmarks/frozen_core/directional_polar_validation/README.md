# 受体羟基方向积分的外部结构检验

第一组从四个人源靶点中按预定输入质量标准保留两个。它们的 8 个 Pro 输出池都缺少距晶体 ≤2 Å 的姿势。第二组专门选择 PPAR 激动剂口袋，提高羟基实际参与接触的概率。在 PPARG 的一个输出池，积分后首位姿势从 1.925 Å 退化至 9.093 Å；加入 MMFF94s 羟基扭转先验仍选中同一个错误姿势。两个协议分别冻结于 `0010f20` 和 `f17c265`，第二组选择时已知道第一组的结果。

## 样本、输入与观测独立性

|PDB / 靶点|输入审查|该组观测|
|---|---|---|
|4BFP / TNKS2|2.40 Å。另一条链的 Lys996 NZ 距配体 7.290 Å，排除。两枚沉积的 Zn 分别距查询配体 16.876 和 21.832 Å|未运行|
|2AA2 / NR3C2|1.95 Å。拟省略的甘油 O1 距配体 5.324 Å，排除|未运行|
|4I1R / MALT1|2.70 Å。25 个重原子、1 个图映射。噻利达嗪使用 CCD 所给中性结构|4 个 Pro 池共 36 个姿势，0 个接触可动受体羟基|
|2YIY / MAP3K7|2.49 Å。33 个重原子、12 个图映射。配体 4 Å 邻域的 7 个水氧在干燥受体中省略|4 个 Pro 池共 35 个姿势，3 个接触可动羟基|
|2PRG / PPARG|2.30 Å。25 个重原子、2 个映射。来自有历史预测的靶点和结构|36 个姿势，11 个接触可动羟基|
|3TKM / PPARD|1.95 Å。31 个重原子、6 个映射。旧研究涉及 PPARD 的其他结构；此 PDB 在原分支未见。GW0742 使用 CCD 所给中性羧酸|35 个姿势，12 个接触可动羟基|

晶体配体中心用于确定 24 Å 已知口袋，起始构象由分子图独立生成。四个可用复合物的查询配体各有一个完整沉积位置；同一蛋白坐标系内对所有化学图映射取最小重原子 RMSD。所有水分子省略；`sources.json` 列出各链各组分最近距离、保留状态、沉积 DOI、原始 PDB/CCD 地址和 SHA-256。没有实验氢坐标或溶液微态占比证据。第一组为新结构的跨靶点检验，第二组为按既有失效机制选择的同家族挑战。重复种子和预算不增加独立靶点数。

## 原始池与同坐标排序

原生 AutoDock Vina 1.2.7 和本轮开始时的 Pro 各自生成候选，两个预算（每条轨迹名义上限 3000/12000 次求值）、两个种子（43103/43117）、4 条轨迹、1 CPU、最多 9 个模式。所有方向性排序共用 Pro 的坐标、输出精度和 native 分数。按 `native_rotors.cpp` 从输入拓扑重算报告分母。各预算、种子和每一项能量都在分组结果 `modes.csv` 和 `cells.csv` 中。

|分组与引擎|Top-1|Top-5|池命中|存在正确候选时 Top-1|
|---|---:|---:|---:|---:|
|第一组 upstream / Pro / 固定氢 / 方向积分，各自|0/8|0/8|0/8|无可评价池|
|PPAR upstream / Pro / 固定氢，各自|1/8|1/8|1/8|1/1|
|PPAR 接触先验积分 / 加 MMFF 先验积分，各自|0/8|1/8|1/8|0/1|

2PRG、预算 12000、种子 43117 的原生第 1、2、3 名分别距晶体 1.925、2.740、9.093 Å。固定氢下正确姿势仅领先原生第 3 名 0.040 评分单位。逐个羟基回放将固定氢与方向积分的差异完全定位在 Ser289 OG，三个候选的原始惩罚增量分别为 0.424、0.097、0.072。在这个残基处，干燥受体的接触先验对旋转角度为常数。除以同一分母 1.40922 后，正确姿势落后错误姿势 0.209。MMFF 先验将正确姿势从第 3 名拉回第 2 名，仍落后 0.197。单独给 native 加羟基项、单独加其他供体方向项、固定氢以及各羟基取最低惩罚时，首位均仍为正确姿势。两个惩罚合用导致换首。单个成功池不足以估计该退化的发生率。原子级反事实见 `ppar/polar_switch.csv`。

### 区分搜索覆盖和末端排序

另从晶体坐标启动显式局优，并在原查询分子的固定内部几何上移植晶体扭转，按优化前距离选图映射后只精修一次。这一诊断使用答案，不计入预测准确率。用同一 Pro CLI、同一查询拓扑和 `--pose_coordinates` 分别重新评价精修后的投影和所有盲对接候选，比较导出的蛋白质–配体相互作用能和配体内能之和。在共用查询和同一 unbound 基线时，该和保留分子的候选排序。

|结构|直接精修晶体 RMSD|晶体扭转移植到查询后的 RMSD，精修前 → 后|投影精修结果与各盲池最佳分数之差（负值有利）|
|---|---:|---:|---:|
|4I1R|0.194 Å|0.947 → 2.302 Å|+1.711 至 +2.404|
|2YIY|0.275 Å|0.277 → 0.326 Å|−9.118 至 −4.211|
|2PRG|0.214 Å|0.380 → 0.491 Å|+0.216 至 +1.353|
|3TKM|0.531 Å|0.248 → 0.443 Å|−10.283 至 −0.872|

2YIY 和 3TKM 在相同查询几何内确有能量优于所有输出姿势的近晶体候选。此处只能定位到进入最终输出池前的覆盖缺口；没有轨迹级观测可进一步区分未访问和访问后丢弃。4I1R 的原始配体几何也影响了局优终点。2PRG 的投影是一个近晶体姿势，盲搜还可产生分数更低的其他近晶体姿势，因此不能把投影的分数当作整个晶体盆地的下界。

## 计算与复现

Pro 使用未改动的 `dec913f584a2bf188ce9704b0602f11c56069627` 构建，本机二进制 SHA-256 为 `8030f0b2e85c5a1931be37d3cc5e1c20f3618b1641ae318f11ccd07b405bb17a`。官方 Vina 1.2.7 Linux x86-64 二进制 SHA-256 为 `f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644`。Pro 每格墙钟 0.189–0.627 s，最大峰值 RSS 22.8 MiB；官方二进制 0.945–1.569 s，最大 529.0 MiB。二进制、编译和实现各异，不能用这里的差值归因羟基模型。父版本 4 条轨迹合计实际求值 12,426–48,725 次，分项见 `evaluation_counts.csv`。积分仅重排已输出姿势，受体预计算约 0.018–0.023 s，每池 0.0006–0.0034 s；PPAR 加入 MMFF 先验的额外受体成本 0.0044–0.0098 s，每池重新评分约 0.0014–0.0033 s。这些是单机观测。

每组 `raw.tar.xz` 包含下载的沉积数据、输入、准备日志、分子图映射、每个引擎的原始输出与运行日志，以及独立晶体诊断的输入输出。`binaries.json`、`sources.json`、`docking_costs.csv` 和 `scoring_costs.csv` 给出版本、采集与成本。准备脚本按历史标准同时记录排除项。可直接从归档重算；若要从头准备，运行 `prepare_validation.py WORK` 或 `prepare_validation.py WORK --ppar`。从归档另建目录完整重跑后，两个引擎的 32 次搜索、逐姿势分数及排名、晶体投影、候选插入比较和 Pro 求值计数均复现，两个阶段七张数值表逐字节一致，PPAR 单羟基反事实结果也完全一致。以下示例从仓库根目录运行，依赖版本见 `sources.json`，上游二进制路径在环境变量 `VINA_127` 中。

```bash
D="$PWD/benchmarks/frozen_core/directional_polar_validation"
B="$PWD/build-dec913"
W=$(mktemp -d)
git worktree add --detach "$W/base" dec913f584a2bf188ce9704b0602f11c56069627
python "$W/base/tools/build.py" --output-dir "$B"
objs=()
for obj in "$B"/obj/*.o; do
    case "$obj" in */main.o|*/split.o) continue ;; esac
    objs+=("$obj")
done
c++ -std=c++14 -O3 -DNDEBUG -pthread -I"$W/base/src/lib" \
    "$W/base/benchmarks/frozen_core/affinity/native_rotors.cpp" "${objs[@]}" \
    -lboost_thread -lboost_filesystem -lboost_serialization -lboost_program_options -o "$W/rotors"
for phase in first_stage ppar; do
    out="$W/$phase"
    mkdir -p "$out"
    tar xf "$D/$phase/raw.tar.xz" -C "$out"
    cp "$D/$phase/sources.json" "$out"
    extra=()
    if [ "$phase" = ppar ]; then extra=(--bonded); fi
    python "$D/run.py" "$out" --pro "$B/vina" --upstream "$VINA_127" \
        --rotor-probe "$W/rotors" "${extra[@]}"
    python "$D/reference_diagnostic.py" "$out" "$B/vina"
    python "$D/summarize.py" "$out"
    if [ "$phase" = ppar ]; then python "$D/explain_switch.py" "$out"; fi
    diff "$D/$phase/cells.csv" "$out/cells.csv"
done
```

评分模型和受体羟基假设分别见同级 `hydroxyl_partition`、`hydroxyl_evidence` 目录。这两个结构面板只提供姿势和成本证据，同 assay 亲和力检验另见 `hydroxyl_evidence`。
