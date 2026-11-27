# 同测定半乳糖凝集素 Kd 排序

Gal-1 第一个随机种子的收益在第二个随机种子上反向。每个靶点对同一批 13 个中性化合物有 13 个 Kd，两个晶体构象各提供一次首姿预测。下表先对同一化合物的两个受体首姿打分取算术均值，再在各自测定内判断 Kd 至少相差三倍的分子对。两个靶点、多个种子及受体状态共享同一批药物，不能相加得到独立样本量。数据见 Peterson 等人的[表 1](https://doi.org/10.1021/acs.jmedchem.7b01626)，包括所有 5–17 号化合物；[归档](../../galectin_affinity/README.md)提供 PDB、受体准备和化学微态。

| 种子，单轨迹名义预算 | Gal-1 上游 / FastDock / 切线，23 对 | Gal-3 上游 / FastDock / 切线，35 对 | 等重原子数且等电荷，FastDock → 切线 |
|---|---|---|---|
| 912743，12000 | 16 / 13 / 19 | 17 / 17 / 17 | Gal-1 3→6/7，Gal-3 0→0/2 |
| 912749，12000 | 19 / 19 / 10 | 22 / 22 / 24 | Gal-1 5→2/7，Gal-3 2→1/2 |
| 912749，6000 | 17 / 17 / 8 | 17 / 16 / 14 | Gal-1 4→1/7，Gal-3 2→1/2 |

初次 12000 的切线对 Gal-1 新判对 7 对、判错 1 对，按化合物配对重采样的 95% 描述区间为 −10 到 +100 个百分点。Gal-3 有 6 对从错变对、另 6 对从对变错，区间 −48.3 到 +47.4 个百分点。第一个种子 Gal-1 的 Spearman / Kendall τ-b 从 −0.052 / 0.013 升到 0.426 / 0.323；第二个种子从 0.377 / 0.271 降到 −0.077 / −0.090。Gal-3 的 12000 结果分别从 −0.080 / −0.039 到 −0.096 / −0.039、从 0.256 / 0.219 到 0.267 / 0.194。按重原子数从大到小直接排序可判对 Gal-1 的 16.5/23、Gal-3 的 34/35，平局给半分。Gal-3 只有 2 对同时满足尺寸电荷匹配和三倍分离。[逐受体、核心系列、十倍差与实测误差分层](rankings.csv) · [按分子重采样](paired_bootstrap.csv)

固定 912743、12000 的 52 个父版本输出池，在正常拓扑和相同打印坐标上，用两臂各重新打分全部 295 个姿势。全部十项打印打分分解一致。此诊断共耗时 8.18 秒。每个单元还检查了全部姿势的糖环六原子到晶体供体同源糖环的无叠合 RMSD，以及最近受体重原子接触。供体配体与查询分子不同，因此这个指标只检验核心定位。首枚种子的 52 次运行中，糖环在首姿、前五、输出池任一姿势进入 2 Å 的次数分别由 3/3/3 变为 3/4/4；次枚 12000 由 5/5/5 降为 2/2/2。首姿均无小于 2 Å 的重原子间受体碰撞。实验亲和力变化与正确的整配体结合模式之间缺乏直接对应证据。

52 次 12000 搜索的总耗时，第一种子 7.548→7.251 秒，第二种子 7.334→7.714 秒；hunt、第二阶段、显式精修合计联合能量/梯度调用分别为 642815→641963、642410→642796 次；局优次数为 9063→9282、11928→11914。墙钟覆盖读入、建网格、对接、打分及写文件；化学输入完全复用同一压缩归档，本轮无新增构象或受体准备。静态打分和网格预处理耗时包含在墙钟内，计数只包含局优的联合调用。父版本与切线绑定同一逻辑核 7，逐配体轮换执行顺序。上游无修改程序与计数重放在另一逻辑核 8 上，跨程序的时间不能解释为切线加速。[逐受体分阶段计数和耗时](cost.csv)

第二种子 12000 的切线墙钟比对照高 5.18%，另以不参与评分的种子 912737 校准。12k 父版本耗时 6.821 秒，切线 6000/12000 耗时 6.146/7.050 秒。线性插值得到切线预算 10500。两枚评分种子实际成本分别为 7.128 / 6.842 秒，低于同种子父版本 12000 的 7.548 / 7.334 秒。较低成本切线的 Gal-1 判对 15、11 对（对照 13、19），Gal-3 判对 19、23 对（对照 17、22）。降低名义预算仍未使 Gal-1 跨种子表现稳定。原始打分、时间校准及其目标选择不使用条件式挑选分子。校准方案及追加分析时点见[协议](PROTOCOL.md)。

## 复算

从仓库根执行。运行所用 Python 3.12、NumPy 2.3.5、SciPy 1.17.0、RDKit 2026.09.2。解包后可独立从原始 PDBQT 和原论文实测表重建全部分子对，不依赖外部网络。压缩包中的 work 存放输出、日志、运行元数据、输入摘要及固定姿势打分；calibration 存放选择 10500 时的原始结果。共用输入位于原有 benchmarks/galectin_affinity/evidence.tar.xz，SHA-256 为 719501125503debd3557e9576d3683c5c4a9ac5ef2b49add181367936270efdb。

```sh
D=benchmarks/hunt_overlap/affinity
W=$(mktemp -d)
tar -xf "$D/runs.tar.xz" -C "$W"
OPENBLAS_NUM_THREADS=1 python "$D/evaluate.py" "$PWD" "$W/work"
python "$D/check.py" "$PWD" "$W/work"
cmp "$D/rankings.csv" "$W/work/rankings.csv"
cmp "$D/cost.csv" "$W/work/cost.csv"
python "$D/calibrate.py" "$W/calibration"
```

要重跑 C++，用 4ba83b5087452639de42502bc0194b6915363c26 的研究源码编译 FastDock 并生成两个命令包装。父程序清除 FASTDOCK_HUNT_TANGENT，实验程序设为 all，两者调用相同二进制。研究版的二进制 SHA-256 为 698dc4283df1cd64a574b5946f2351f7946b13cd7ad2929c483a2002be89ca92；上游 3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645 的发行版和被动计数版分别为 3b6492ecf51789c71195438b3bda5784455b5e7ab2bceb7a59a25c5c8d182ae8 和 62875fa8d1abc9e1af09f848f6d9cf052828359dd5246a6ef00d0300ca4672b2。该上游提交的 src/lib、src/main 与 v1.2.7 发布版相同；[计数补丁](../../pair_gradient/upstream_counts.patch)在全部 156 个交叉条件上输出逐字节相同的结果。工作流示例使用环境中已配置的四个命令绝对路径，核编号按主机调整。

```sh
for spec in '912743 12000' '912749 12000' '912749 6000'; do
  set -- $spec
  python "$D/assay.py" run "$PWD" "$W/new" --seed "$1" --cap "$2" --cpu 7 --parent "$PARENT" --candidate "$CANDIDATE"
  python "$D/assay.py" run "$PWD" "$W/new" --seed "$1" --cap "$2" --cpu 8 --upstream "$UPSTREAM" --counts "$UPSTREAM_COUNTS"
done
python "$D/assay.py" fixed "$PWD" "$W/new" --cpu 7 --parent "$PARENT" --candidate "$CANDIDATE"
for seed in 912743 912749; do
  python "$D/assay.py" run "$PWD" "$W/new" --seed "$seed" --cap 10500 --cpu 7 --candidate "$CANDIDATE"
done
OPENBLAS_NUM_THREADS=1 python "$D/evaluate.py" "$PWD" "$W/new"
python "$D/check.py" "$PWD" "$W/new"
```
