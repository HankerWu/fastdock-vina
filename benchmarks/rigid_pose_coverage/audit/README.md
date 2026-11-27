# 从原始源态和实验记录复核筛选与捕获

父实验的八个原始姿势归档通过 SHA-256 检查，解档后的 1,006 个成员亦逐一匹配。用沉积结构、CCD 和原始搜索快照重算 16 格两臂、18 个 AT9283 源态两臂、32 行漏项和四次真正的原生知情分叉，指定的四张表与父提交全部逐字节相同。[核验收据](results/parent_replay.json)

## 几何与已知实验坐标

以下距离均用 Å。在同一可行域上计算每个有序对应到最近保留项的全配体有标签位移，分位数先在每格求，再对 16 格等权。原边长臂可有越盒保留项；其行也对同一盒内待覆盖域求最近邻。每格均有 1,024 个候选送入粗场。

|门函数|最近位移中位 / P90 / P95 / 最大|保留库晶体最小误差均值|有 ≤2 Å 的格|≤2 Å 的对应数|
|---|---|---:|---:|---:|
|原边长|3.287 / 4.736 / 5.131 / 8.009|2.500|5 / 16|42|
|盒内边长|3.256 / 4.735 / 5.155 / 8.075|2.468|5 / 16|45|
|盒内覆盖|2.197 / 2.794 / 2.923 / 3.257|2.690|3 / 16|3|
|盒内均匀|2.029 / 3.034 / 3.345 / 5.834|2.866|2 / 16|7|

门前共有 221 个近姿，闭盒保留全部 221 个。覆盖在 16/16 格缩小最大最近距离，恢复 A1120 A⁻，保住原有两格，丢掉 BRD9 的 91813 及 3FR2→2NNQ 两个种子。A1120 HA、5TBM→4GHI 没有恢复。三个漏例按照原始 ID 都可以在新门前产生；区别发生在把几万个全配体位置压缩为 1,024 个时。原门前有近姿的格仍为 8 格。AT9283 的 18 个旧源态在整分子有序对应穷尽时全部大于 2 Å；任何严格子集继承相同下界。[全部对应库的保留编号](results/selected_ids.npz) · [64 条逐格数据](results/historical.csv) · [逐近姿的最近保留项](results/near_neighbors.csv)

## 同一源态下的真实后续搜索

在历史库几何解盲后、运行下游前，[追加规则](CAPTURE.md)冻结了三个旧漏例及近姿入选规则。三个盲首选均由真实类型场选择。仅当保留库有近姿时，另行注入排名最高的近姿。

|病例与路径|库最低|入选 / hunt / 二次局优|逐出步|发布首位 / Oracle|原边长发布首位 / Oracle|
|---|---:|---|---:|---|---|
|A1120 A⁻，覆盖盲首选|1.945|13.210 / 13.899 / 13.899|4|15.479 / 14.807|2.371 / 2.315|
|A1120 A⁻，晶体知情第 11 名|1.945|1.945 / 1.941 / 2.170|2|2.356 / 2.356|同一边长对照|
|A1120 HA，覆盖盲首选|2.822|2.822 / 2.433 / 2.130|3|2.371 / 2.099|2.369 / 2.369|
|5TBM→4GHI，覆盖盲首选|2.947|12.535 / 12.535 / 12.535|6|13.958 / 12.789|16.222 / 7.317|

5TBM 的最终首位更近、池内最佳却更远。A⁻ 的知情分叉改变了整条后续轨迹，即使最终进入 2.4 Å，原始注入身份仍在第 2 步丢失。晶体知情结果单列。[逐阶段与身份](results/capture.csv) · [四次实耗](results/capture_cost.csv) · [原始输出](results/capture_raw.tar.xz)

## 同测定方向的独立复核

直接从每个发布文件第一份 REMARK VINA RESULT 取值，按 Holt Table II 重建七个点及一个区间，15 个可靠三倍方向在单受体与等权读法下分别复计。覆盖、两种边长、官方 24k 与实际成本 Pro 96k 在等权读法上恰好错同五对，均把实测 12 nM 的化合物 6 排在相对较弱的 4、5、7、10、11 之后。均匀法同为 10/15，却判对 6/10、判错 9/11。八个配体构成的 15 对共享配体。[逐对身份](results/pair_identity.csv) · [从所有发布坐标独立回算 AD4](../fkbp/results/scores.csv)

## 复现

从仓库根目录，PARENT、W 等为绝对路径；OBSERVER 按父实验编译，校验其 SHA。使用一线程 BLAS。归档打包时的绝对路径、日志计时和文件时间戳随机器变化。

```sh
D="$PWD/benchmarks/rigid_pose_coverage/audit"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$D/replay_parent.py" "$W/old" "$OBSERVER" "$W/parent_replay.json"
python "$D/frozen.py" "$W/old/parent" "$W/frames" "$W/recount"
python "$D/capture.py" "$W/old/parent" "$W/new-capture" "$OBSERVER" "$W/capture-results"
python "$D/pair_check.py" "$FKBP_WORK" "$D/../fkbp/results/assay.csv" "$W/pair_identity.csv"
cmp "$W/recount/historical.csv" "$D/results/historical.csv"
cmp "$W/capture-results/capture.csv" "$D/results/capture.csv"
cmp "$W/pair_identity.csv" "$D/results/pair_identity.csv"
```
