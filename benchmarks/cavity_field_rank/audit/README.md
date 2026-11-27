# 从原始坐标、原测定和显式原子对复核

[direct_space.py](direct_space.py) 从 PDB 和 CCD 建立实验配体、等价映射及接收体局部 Cα 框架。BRD9 两条链的编号相差 116，按残基身份及双侧 14 Å 口袋筛选对应。酸根允许两个氧共振交换，HA 保留羰基和羟基身份。每个候选由源坐标与刚体变换重构。下表中的任意刚体最优配准只作标签辅助的表示能力下界；预测误差始终在固定受体框架内评估。独立路径共享存档的化学身份及常用 PDB/CCD 读取器，分别做重原子对应、四元数受体拟合、指标与阶段计算。

对固定重原子形状和保留化学身份的置换，平方 RMSD 分解为质心位移平方与居中坐标误差平方之和。任取一个已生成的质心和一个已生成的取向，二者的最小分项可以相加，得到全部 1,024² 交叉组合的严格最小平方误差。若再允许连续平移和旋转，每种等价映射各作一次无反射最优刚体拟合，最后对映射取最小。[每格三个边界](results/representation.csv)

10 格中 9 格的固定内部形状有能力到 2 Å，实际库有 3 格；保留的质心和取向全交叉后有 6 格。可以由这次交叉重新组合生成近姿的三个原空库是 A1120 两微态的新种子 93241 和 EPAS1 4GHI 的 93257。5TBM–93257 任意配准仍至少 2.431 Å。按 A⁻、HA、4GHI 的顺序，一格 1,048,576 种交叉中仅 10、1、7 对达到 2 Å。均匀随机选取一个保留全部质心和取向的双射时，命中数的数学期望分别为 10/1024、1/1024、7/1024，均小于 0.01。这些计数使用晶体，仅作可达性诊断。该分解只界定进入柔性 hunt 之前可用的几何。

## 移除空间格点的事后对照

对于相同候选逐个原子枚举距离小于 8 Å 的受体原子，用 Vina 五项原生函数及其权重，按 1/32 Å² 平方距离表的两个端点取均值。保持与网格臂相同的原子类型、每原子正能量限幅、盒壁及单位。直接在真实原子间分离处查能量，完全省去空间格点的三线性插值；保留原势的距离分箱。Python/SciPy 空间索引独立计算，与 [types.cpp](types.cpp) 调用原生距离表、遍历全部受体原子的 93 次抽查相差至多 1.9×10⁻⁹ kcal/mol。[逐候选五列数组](results/direct.npz) · [列名、时间、原生对照](results/time.json)

真实/置换各 10 格，共 20 个排名。移除空间格点后 19 个胜者相同；BRD9–91823 真实类型胜者由 399 变为 727，误差 3.156→3.560 Å。BRD9 两格真实类型的第一近姿名次由网格的 2/12 变为直接求和的 2/14。对于两个已知可捕获近姿，最近受体重原子分别为 2.949/2.662 Å，均无 <2 Å 接触。这里的最小距离是穿插诊断，不足以判断有利氢键、脱溶剂化、受体重排或某个盆地的实际权重。[完整对照](results/direct.csv)

## 带晶体标签的分叉

追加方案 bf6a052 在查到候选覆盖和盲结果后、执行干预前提交。三个有近姿的库有五个不同的“本排序内第一个 ≤2 Å”候选。分别从同一源态选出该候选，再完整搜索。五个注入 ID 都在第 2 步逐出。四个分叉的最终首位 ≤2 Å，另外一个的实际池最近为 2.570 Å。BRD9 第二枚种子的形状首个近姿与两场首个近姿属于不同的局优起点。后者（51）在真实和置换中排 12/15，经 hunt 达 1.852 Å；前者（115）在形状场排 34，经 hunt 达 2.343 Å。未将这些有标签的发布结果计入盲测。[方案](FORK_PROTOCOL.md) · [逐例与直接 ID](results/near_forks.csv) · [诊断成本](results/fork_cost.csv) · [原始分叉](fork_raw.tar.xz)

[verify_endpoints.py](verify_endpoints.py) 还从原始阶段坐标复算 7,129 项，从 70 个实际池复算全部 594 个模式和几何终点，从原始 Ki 表与每份首位 PDBQT 重新计算两受体及其等权读分的 21 组方向与相关。候选和发布几何差均小于 2×10⁻¹⁴ Å。[核验摘要](results/verification.json) · [独立同测定结果](results/independent_assay.csv)

## 复算命令

先按上一级 pose 和 fkbp 目录说明将各原始包恢复至 P 和 F，D 为仓库根目录下的基准路径，W 为新的绝对工作目录，ENGINE 为本轮构建的 XS 研究引擎。OBJ 指向同一生产源树的 tools/build.py 对象目录。如果 Boost 在系统标准路径，可省去 -I、-L 和 -Wl 参数。

~~~bash
D=$PWD/benchmarks/cavity_field_rank
mkdir -p "$W/tmp" "$W/results"
export TMPDIR="$W/tmp" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sh "$D/audit/build_types.sh" "$OBJ" "$W/types" -I"$BOOST_PREFIX/include" -L"$BOOST_PREFIX/lib" -Wl,--disable-new-dtags,-rpath,"$BOOST_PREFIX/lib"
python "$D/audit/direct_space.py" "$P" "$W/frames" "$W/types" "$W/results"
python "$D/audit/near_forks.py" "$P" "$D/pose/results" "$W/forks" "$ENGINE"
python "$D/audit/verify_endpoints.py" "$P" "$F" "$D/pose/results" "$D/fkbp/results" "$W/results"
~~~
