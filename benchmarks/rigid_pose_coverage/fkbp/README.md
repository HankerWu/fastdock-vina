# FKBP12 同测定的覆盖检验

八个分子使用 Holt 等论文第 9927 页 Table II 的表观 Ki。化合物 4–11 依次为 186±24、110±20、12±5、250±30、10±1、7±2、300–600、300±40 nM。作者使用重组人 FKBP12、α-胰凝乳蛋白酶耦联的 PPIase 测定、10 °C、35 mM HEPES、pH 7.8。七个点估计算秩相关，八个分子的最近区间端点相差至少三倍的 15 对算方向，平分各半。沿用前轮定义时 ± 栏不参与三倍筛选。两个受体 1FKG 和 1FKH 分开报告，再等权平均两张图的发布首分。对每个臂，只有一个种子 1047349，两个受体复用同一组标签。这是历史系列的新搜索源态。[原论文及表 II](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1993%20%20%28vol%20115%29/22%20%20%2809879-10470%29/09925-09938.pdf#page=3) · [事前登记](../PROTOCOL.md)

## 亲和力与姿势

分差是弱配体减强配体的经验分（kcal/mol），正号为方向正确。姿势误差由有独立晶体坐标的化合物 8、9 在两张受体的四格求均值，单位 Å。其余六个分子参加排序，保留实际发布坐标。

|方法|1FKG /15|1FKH /15|等权 /15|等权 ρ / τ-b|8→10 / 9→11 分差|Top-1 / Oracle / Regret|
|---|---:|---:|---:|---:|---:|---:|
|官方 24k|6|10|10|.357 / .143|+.641 / +.015|3.756 / 3.316 / .440|
|普通 Pro 24k|12|8|8|.214 / .143|−.079 / +.484|3.921 / 3.622 / .299|
|原边长|8|10|10|.286 / .048|+2.460 / +.304|3.758 / 3.038 / .720|
|盒内边长|8|10|10|.286 / .048|+2.819 / +.269|3.120 / 2.924 / .196|
|盒内覆盖|7|10|10|.357 / .143|+.652 / +.912|3.814 / 2.950 / .864|
|盒内均匀|10|10|10|.143 / .048|+2.437 / −.284|3.082 / 3.082 / .000|
|均匀加等量计算|10|10|10|.143 / .048|+2.437 / −.284|3.082 / 3.082 / .000|
|普通 Pro 48k|10|5|8|.214 / .143|−.414 / +.495|3.313 / 2.819 / .494|
|官方 48k|6|10|9|.143 / .048|+.790 / −.381|3.125 / 2.901 / .223|
|实测成本普通 Pro 96k|10|9|10|.357 / .143|+.500 / +.423|3.307 / 1.870 / 1.437|

覆盖法在 1FKG 相比普通 Pro 少对 5 对，在 1FKH 多对 2 对。等权判对数均为 10/15。与均匀法的逐对方向有两处差别，覆盖判对 9/11、判错 6/10，均匀法在这两对上相反。[逐对独立重算](../audit/results/pair_identity.csv)覆盖相对原边长的八个等权分差，与重原子数及 TORSDOF 的 Spearman ρ 为 −.165 和 +.179。相对普通 Pro 为 +.330 和 +.575。样本量为 8，全部形式电荷为 0。十臂的 160 个池全部发布 9 个模式。1,440 个坐标的图和指定立体化学均通过，最大重键漂移 0.00131 Å，最短受体–配体重原子距离 2.397 Å。独立按 AD4 原子类型图、电势、去溶剂化图及扭转项复计，最大分差 0.00455 kcal/mol。十臂 × 三受体汇总 × 15 对的 450 个方向与写出分一致。固定微态和受体，两张图的等权没有热力学丰度依据。数据支持同测定经验排序，缺少结合/游离两相参考、标准态变换和收敛采样以估计 ΔG 或 ΔΔG。[逐受体排序](results/assay.csv) · [全体分子对](results/pairs.csv) · [完整姿势分布及各分子分位数](results/distributions.csv) · [配对差](results/paired_pose.csv) · [化学与独立分数](results/modes.csv)

正式计算的中断曾留下旧命名管道。控制器现为每次调用建立私有管道目录。在刻意留下旧管道的干净目录，从普通输入重跑 1FKG 化合物 8 的五臂，40 个含源态、变换、候选、轨迹与发布坐标的文件逐字节一致。[中断重试实测](results/interrupted_retry.csv)

## 候选身份与完整成本

普通配体经原生搜索第一次变异后回传精确源坐标，所有臂复用同一轨迹的 RNG 分支及扭角。盒内资格要求所有重原子在原 AD4 盒内，未进行晶体匹配。三边、全配体形状场、真实类型 AD4 场和 FIFO 单次注入接续原实现，AD4 粗查的正值限幅为 10 kcal/mol，hunt/二次局优使用完整原生 AD4。每臂对 16×1,024 个候选共做 514,048 次重原子粗类型通道插值和相同数目的形状查询。原边长臂有 19,182 个越界原子位置并接受原壁势，新域的四臂均为零。盒内库每个源态有 23,662–27,825 个对应。覆盖逐次更新最短距离，共 419,308,544 个 12 维平方距离。等量计算的均匀臂做完全相同的覆盖选择及距离更新后丢弃。两均匀臂 16 格的候选 ID、变换、选中项及每个最终池逐字节一致。五臂首个粗分都没有并列。

四个有晶体格穷尽全部对应时，最小 RMSD 依次为 1FKG–8 3.057、1FKG–9 3.886、1FKH–8 3.189、1FKH–9 3.831 Å。盒内域保留这些最小值，任一筛法都无法从这些源态直接取得 ≤2 Å。五个方案的注入 ID 在第一轨迹结束时分别存活 2、3、3、3、3 次，均无直接发布。1FKH–8 的覆盖首选为 5.281 Å，经 hunt 变为 9.087 Å，二次局优到 8.786 Å，进入档案后在第 28 步逐出。第 225 步出现后来的 ≤2 Å 候选，身份与原注入项不同。[门前后及全域](results/gate_crystals.csv) · [候选与类型首位并列数](results/gates.csv) · [各阶段去向](results/stages.csv)

|16 格总计|墙秒 / CPU 秒|hunt + 二次局优求值|下游 AD4 通道插值|
|---|---:|---:|---:|
|官方 24k|29.93 / 29.83|1,153,306|108,569,106|
|普通 Pro 24k|25.59 / 22.70|1,152,971|108,537,069|
|原边长|50.27 / 47.39|1,153,075|108,546,147|
|盒内边长|40.49 / 40.42|1,153,060|108,544,728|
|盒内覆盖|67.04 / 62.72|1,152,986|108,539,262|
|盒内均匀|39.72 / 39.70|1,152,924|108,532,509|
|均匀加等量计算|53.63 / 53.61|1,152,924|108,532,509|
|普通 Pro 48k|39.57 / 38.20|2,304,969|216,969,453|
|官方 48k|59.50 / 58.45|2,305,485|217,015,959|
|实测成本普通 Pro 96k|58.31 / 58.17|4,608,963|433,832,349|

登记的线性估算为 3.583×24k，向上取整运行 96k。覆盖用时为实测成本普通 Pro 的 1.15 倍墙钟、1.08 倍 CPU，普通 24k 的 2.62/2.76 倍。任务调度和并行负载会影响计时。等量计算阴性对照实测用时低于覆盖，其中下游轨迹也不同。以上计时含读图、Python 启动（提议臂）、候选生成/粗排、原生下游与轨迹/坐标写出，普通 Pro 采用未加观察器正式程序。原地图及配体解包和八份键表生成另耗 0.856 墙秒、0.841 CPU 秒。离线建图费用无本次实测。每个普通 Pro 另跑计数构建，48 个池逐字节一致，32 个官方计数重放也与发行版一致。多臂同时启动时，三次在线控制器和两次官方计数重放收到 SIGKILL，一次旧 FIFO 造成空重试。原始部分输出、已完成冗余计时和失败日志在 replay 档中，进程被强杀的区间无完整总时，方法间成本比较只使用成功完成的整格。清洁复跑中，五臂对 1FKG–8 的 40 份源态、库、变换、原轨迹与发布文件逐字节一致。[正式及计数费用](results/cost.csv) · [未加观察器 Pro](results/unmodified_pro.csv) · [成本决策](results/cost_budget.json) · [清洁重放](results/fresh_replay.csv)

## 复现

下列命令从仓库根运行。WORK、BUILD、FRESH、RESTORE、VERIFY 为独立的绝对目录。OFFICIAL 为官方 Vina 1.2.7 Linux 发行版，PRO 为未加观察器的当前 Pro，UPCOUNT 为父实验的上游 1.2.7 计数构建。下游程序的构建、官方来源和 SHA 均可由以下文件追溯。[官方计数构建](../../cavity_field_rank/fkbp/README.md) · [本轮使用的原生观察器构建](../../fragment_anchor/fkbp/build.sh) · [每格命令的运行配置和文件哈希](raw/manifest.json) · [原始沉积、地图、配体和测定输入哈希](results/inputs.json)

```bash
D=$PWD/benchmarks/rigid_pose_coverage/fkbp
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
sh "$D/../../fragment_anchor/fkbp/build.sh" "$BUILD" --boost-prefix "$BOOST_PREFIX" --jobs 1
ENGINE=$BUILD/bin/vina
python "$D/trial.py" prepare "$WORK"
python "$D/trial.py" run "$WORK" "$OFFICIAL" "$ENGINE" "$UPCOUNT" --arms official base edge feasible_edge coverage uniform uniform_work high official_high
python "$D/trial.py" ordinary "$WORK" "$PRO"
python "$D/trial.py" budget "$WORK"
python "$D/trial.py" run "$WORK" "$OFFICIAL" "$ENGINE" "$UPCOUNT" --arms cost_pro
python "$D/trial.py" ordinary "$WORK" "$PRO" --cost
python "$D/trial.py" check "$WORK" "$VERIFY"
python "$D/gate_crystals.py" "$WORK" "$VERIFY/gate_crystals.csv"
python "$D/replay.py" "$WORK" "$FRESH" "$ENGINE" "$VERIFY/fresh_replay.csv"
```

重算已发布实验无需重新选择成本预算。下列流程重新从原档制备普通输入，解包全部已注册输出并校验逐文件 SHA-256，重新读取原子坐标、同测定标签及沉积配体所有拓扑/手性等价映射。受体以共同的 14 Å 邻域 Cα 配准。预测配体不参与拟合。

```bash
python "$D/trial.py" prepare "$RESTORE"
python "$D/trial.py" unpack "$RESTORE" "$D/raw"
python "$D/trial.py" check "$RESTORE" "$VERIFY"
python "$D/gate_crystals.py" "$RESTORE" "$VERIFY/gate_crystals.csv"
cmp "$VERIFY/assay.csv" "$D/results/assay.csv"
cmp "$VERIFY/gate_crystals.csv" "$D/results/gate_crystals.csv"
```
