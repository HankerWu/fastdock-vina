# 同一批石杉碱甲坐标在两份干受体中的排序与局部保留

查询为形式中性的 (−)-huperzine A，固定水模型和受体权重均未参与。三批起点是种子 1103501 在 4EY6 的真实发布模式。另以准备好的配体刚性叠合至 4EY5 的 HUP，构造可由当前表示实现的近参照正对照。原配体有 18 个重原子和唯一的立体相容映射，唯一声明轴的下游重原子位于轴上。这组实验主要检验整体放置和环境的作用。[实验登记](../PROTOCOL.md) · [原子与几何变换](inputs.json)

4EY6、4EY7 源自分别与加兰他敏、多奈哌齐共晶的人 AChE。相同的 26 个起点用 48 个共同口袋 Cα 从 4EY6 移到 4EY7，拟合误差 0.276 Å。评价时各自从 4EY5 直接配准；用 HUP 18 个原子量得，经过 4EY6 的绕行配准与直达 4EY7 相差 0.023 Å。跨态固定姿势分差包含该拟合约定。受体单构象、固定坐标，所有水省略。4EY5 原位含 DMS A605，其 C1 距 HUP 最近重原子 4.694 Å，触发登记的 8 Å 排除规则。三份原计划受体中两份可运行，4EY5 的 27 个起点不计分。整个实验涉及一份配体晶体参照。[结构与配体沉积](https://www.rcsb.org/structure/4EY5)

## 可实现的近参照起点与真实终点

| 接收受体 | 近参照起点 → 原生局优终点 RMSD（Å） | 局优后导出坐标分数 | 该分减去各原池最佳局优分（Pro 24k / Pro 48k / 官方池） |
|---|---:|---:|---|
| 4EY6 | 0.1285 → 0.4292 | −9.577 | −0.165 / +1.280 / +1.315 |
| 4EY7 | 0.1285 → 0.4505 | −7.266 | +1.437 / +2.091 / +2.106 |

单位为 Vina 名义 kcal/mol。正差表示近参照端点排名靠后。三批姿势均在同一 Pro 算法与同一零点下重算。盲搜结果另行记录。对同一 4EY6 受体，投影局优点能否胜出，取决于候选池是否已包含更受该势偏好的竞争者。[逐格对照](contrasts.csv) · [所有端点、哈希及导数次数](endpoints.csv)

固定池原几何 Oracle 与原首位之间的分差，在 Pro 24k 始终为零（同一身份）；Pro 48k 从 4EY6 的 +1.445 变成 4EY7 的 +0.651，官方池从 +3.388 变成 +4.164。对两个端点分别局优后，4EY7 的对应分差为 +0.654、+4.546；4EY6 维持原数。没有一组不同身份在两份受体间翻转。仅迁移近参照几何，也无法保证当前目标函数在两个受体上都把它排到第一。[复核与逐残基势函数分解](../contacts/README.md)

Vina 的 score_only 在同一个 unbound=0 下用显式受体重算 INTER、INTRA 及构象独立修正。HUP 本例的 Pro INTRA 和有效重转子归一化均为零。local_only 经 Vina::optimize 和 quasi_newton，使用原生平滑查表目标；导出后坐标以 0.001 Å 表示，公开 score 再经快速表计算。保留这两类数值。例如 4EY7 投影的局优器返回 −7.266357，原生 local_only 的未序列化公开分为 −7.273，导出后重计为 −7.266。对照 26 个已发布姿势，以 Pro 重计的能量顺序三池全部复现。与原始日志 INTER 的绝对差有四行超过 0.010，最大 0.018；出现在 3 个不同的非首位坐标。官方发布还有不同的扭转归一化，跨引擎绝对打印分另列。[逐格原始 INTER、打印分及次序](same_state.csv) · [局优和计分源码](../../../src/lib/vina.cpp)

54 个合格起点均完成。局优目标没有升高，观察器记录 1,325 次能量/梯度调用，各终点均与原生二进制的 PDBQT 逐字节相符。初态评分、原生局优、导出后重评分分别耗用 2.483、33.007、2.446 墙秒和 2.471、32.923、2.436 CPU 秒。观察器重跑另用 33.198 墙秒、33.095 CPU 秒。36 Å 盒用于容纳原池远端错姿，全部起终点均通过显式范围检查。受体准备、文件读取和变换耗时另计在 inputs.json。无温度、平衡占比、标准态或约束释放的估计。[时间明细](costs.csv) · [原始结构、分数、日志及观察器输出](raw.tar.xz)

## 重放

在仓库根目录运行。Python 3.12、RDKit 2026.9.2、Meeko 0.8.0、Gemmi 0.7.5、NumPy/SciPy/pandas、OpenMM 8.6.1（CPU）和 PDBFixer 1.12.0；C++ 工具使用 GCC 14.2.0 与 Boost 1.83。下面用 BOOST 指向包含 include 与 lib 的 Boost 前缀。观察器只给现有计算加入记录。

```sh
D=benchmarks/ache_receptor_state/diagnostic
W=$(mktemp -d)
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python tools/build.py --boost-prefix "$BOOST" --output-dir "$W/pro" --jobs 4
bash benchmarks/affinity_conditioning/build_refine_probes.sh "$W/pro" "$BOOST" "$W/observer"
python "$D/measure.py" prepare "$W/data"
python "$D/measure.py" run "$W/data" --pro "$W/pro/vina" --observer "$W/observer/vina"
for f in endpoints.csv contrasts.csv same_state.csv; do cmp "$D/$f" "$W/data/$f"; done
```
