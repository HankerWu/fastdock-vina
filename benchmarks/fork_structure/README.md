# 第十次试探的四来源结构试验

候选源为 `1c672a8e2ce12c7a16e20172a3315ea5a19c23d8`，轮前为 `0449b1d0bc3adec768dd9a9f243f3ce5a29cef58`。固定交叉首姿平均降低 0.611 Å，MAPK14 外留受体交叉平均降低 1.250 Å。两组平均值均通过[结构事前判据](PLAN.md)。固定交叉仍高于官方与父版 2.542 Å，大部分外检改善来自同一个交叉方向。亲和力及综合晋级判断见[累计台账](../cumulative/README.md)。

| 面板与每臂格数 | 官方首姿均值 Å | 父版 | 轮前 | 候选 | 相对轮前改善／持平／恶化 |
|---|---:|---:|---:|---:|---:|
| 固定交叉 20 | 2.863 | 2.863 | 6.016 | 5.405 | 11／1／8 |
| 固定同晶 20 | 6.794 | 6.794 | 7.175 | 6.792 | 12／0／8 |
| 固定历史成本分配 20 | 4.709 | 2.536 | 4.746 | 3.891 | 13／0／7 |
| 原生默认 5 | 2.372 | 2.372 | 2.774 | 2.701 | 2／0／3 |
| MAPK14 新受体同晶 4 | 2.627 | 2.627 | 6.151 | 4.836 | 3／1／0 |
| MAPK14 新受体交叉 4 | 6.184 | 6.184 | 5.759 | 4.508 | 2／0／2 |

固定发布 65 格无空池，候选有 569 个模式。三份基线使用同一面板、种子、真实源码及同一 ELF 的历史原始执行，其中官方默认含既有超时和同种子 900 秒补跑。按每靶点与受体方向给出全部种子分布和配对，原始输入哈希、各臂冻结预算、候选的每次命令、PDBQT 和日志均可追溯。[固定逐靶分布](fixed/summary.csv) · [逐向分布](fixed/by_direction.csv) · [逐格配对](fixed/paired.csv) · [各格端点](fixed/cells.csv) · [全部候选发布坐标](fixed/fixed/raw.tar.xz) · [基线来源及输入](../geometry_comparison/upstream/four_versions/)

| 固定交叉分布 Å | Q1 | 中位数 | Q3 | P90 | P95 | Best@5 | Oracle | Regret | 首姿 >5 Å |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 官方／父版 | 0.779 | 2.567 | 5.065 | 5.649 | 5.867 | 2.502 | 2.482 | 0.380 | 5/20 |
| 轮前 | 2.564 | 5.400 | 8.242 | 12.602 | 14.504 | 5.061 | 4.925 | 1.091 | 11/20 |
| 候选 | 0.968 | 5.432 | 7.081 | 12.030 | 13.811 | 4.315 | 3.959 | 1.446 | 11/20 |

按靶点，固定交叉首姿均值从轮前到候选分别为 DYRK1A 5.011→1.838、EPAS1 5.338→4.124、FABP4 5.544→4.535、FKBP1A 6.043→6.228、RBP4 8.143→10.299 Å。池内近姿覆盖改善大于首姿改善，几何排序损失增加 0.355 Å。

## 冻结受体上的验证

[1BL7](https://www.rcsb.org/structure/1BL7) 的 SB4 与 [1DI9](https://www.rcsb.org/structure/1DI9) 的 MSQ 属于人 MAPK14。两份 A 链 DBREF 都直接使用 Q16539 残基编号。两套冻结输入通过原有排除规则，未替换。外检的独立单位是两份受体结构；MAPK14 家族有过历史开发。直接实验坐标可用于同晶参考，交叉值来自 50 个共同口袋 Cα 的最小二乘配准，拟合 RMSD 0.661 Å。搜索给定已知口袋中心，配体初始构象从化学图独立生成。各臂检查了命令、半径、CPU、反向配准和索引，独立从沉积 PDB/CCD 得到 SB4 的 4 个和 MSQ 的 1 个合法立体化学映射。所有 32 格共 288 个模式，最大复算差 1.87×10⁻¹⁴ Å。按逆臂顺序再次执行全部 32 格，发布池、分数、各阶几何及局优计数完全一致。[原始输入](external/inputs.tar.xz) · [原始坐标与日志](external/raw.tar.xz) · [逐格逐姿](external/ranks.csv) · [各向统计](external/by_direction.csv) · [独立复跑](external/repeat/)

| MAPK14 交叉与种子 | 轮前 → 候选首姿 Å | 轮前 → 候选 Oracle Å |
|---|---:|---:|
| 1BL7 → 1DI9，85231 | 6.576 → 6.579 | 3.616 → 4.024 |
| 1BL7 → 1DI9，85237 | 1.718 → 6.576 | 1.718 → 3.670 |
| 1DI9 → 1BL7，85231 | 7.464 → 2.622 | 3.521 → 2.622 |
| 1DI9 → 1BL7，85237 | 7.276 → 2.256 | 4.422 → 2.256 |

新受体交叉的候选均值优于三份基线、同晶均值仍落后两份上游。1BL7→1DI9 的两个种子都产生了更差的池 Oracle。该失败方向已经揭盲，不应再作为独立泛化证据。

## 可达成本、成本复现性与开销

240 个前瞻标定格用 406.270 CPU 秒，逐作业墙钟合计 514.251 秒，峰值 558,760 KiB。每个方向和每臂的 3k、12k、120k 都重复实测，120 对池全部逐字节一致。仅 51/120 对的 CPU 比落入注册的重复性区间，最大值／最小值的中位数 1.315、最大值 5.619。10/10 方向至少一项计时不稳定，7/10 至少一臂的三点中位数非严格递增，6/10 至少一臂的实测低高档未包围共同目标，原因可重叠。没有方向进入插值预算验证和独立种子确认。候选对各基线前瞻实际匹配均为 0 项，配对分母 0，不能估计等 CPU 改善。历史成本分配中的真实匹配分母分别为 12/20、2/20、9/20，平均首姿差分别为 −0.164、+3.169、−0.018 Å，三组涉及的方向和种子不同。运行平台和并行负载未实施隔离，这轮记录无法将重复性不足归因到某个硬件或运行阶段。[每一对重复计时](timing/repeats.csv) · [每一方向的排除](timing/caps.json) · [前瞻实际分母](timing/matched/actual_cost_summary.csv) · [历史实际分母](fixed/actual_cost_summary.csv)

| 本轮支出 | 运行数 | CPU 秒 | 逐作业墙钟秒 | 峰值 RSS KiB |
|---|---:|---:|---:|---:|
| 固定面板候选 | 65 | 245.367 | 368.569 | 43,748 |
| 新受体首次运行，官方 | 8 | 15.300 | 15.451 | 474,808 |
| 新受体首次运行，父版 | 8 | 6.898 | 7.043 | 61,836 |
| 新受体首次运行，轮前 | 8 | 9.040 | 10.614 | 59,424 |
| 新受体首次运行，候选 | 8 | 4.989 | 8.160 | 59,420 |
| 新受体逆序重跑四臂 | 32 | 41.194 | 50.411 | 474,828 |
| 新受体含下载制备 | 2 | 8.856 | 22.391 | 未测 |
| 新受体从沉积原件独立重制 | 2 | 12.274 | 12.957 | 196,128 |
| 因二进制身份不合格排除 | 68 | 58.249 | 59.758 | 63,072 |

候选在 65 格的 hunt、second、explicit 联合求值分别为 37,237,264、389,280、15,511。在新受体八格中分别为 266,474、22,128、2,209。读数覆盖输入、建图、搜索、原生精修与发布，未执行附加重排。逐阶段 CPU 无独立测量；官方和父版无原生 E/G 计数。固定输入的历史交叉制备另计 22.817 CPU 秒，同晶旧制备成本缺失。首次使用缓存路径时，标为稳定的 ELF 散列为 87b28d5d…，未通过身份检查。对应外部 8 格和先导 60 格已隔离并重做。第一次轮前洁净构建实际生成了正确 ELF，但检查脚本初次指向该错误缓存作为比较对象，当次 CPU 和 RSS 缺失。重新从轮前提交洁净构建并与已核验的旧制品逐字节比较，所得 SHA-256 为 f6ff9e0d9d6f0bc6e77e65409f9c7f88fccd826eca064a49e9c91a9a73496f9f。该次编译用 140.870 CPU 秒；研究候选的发布源码编译用 122.789 CPU 秒。上游两份 ELF 的 SHA 均与已有独立洁净重编相同。完整版本、二进制散列、帮助与默认值、构建命令和费用在[身份清单](builds.json)、[轮前重编](stable_build.json)、[候选编译](candidate_build.json)、[排除记录](rejected_identity/exclusion.json) 和[官方／父版重编](../stage_local_structure/restored/upstream_rebuild.json)。

## 复算与再执行

从仓库根目录运行。数值与制备依赖沿用前一轮，含 NumPy、SciPy、RDKit、Meeko、Gemmi、OpenMM 和 PDBFixer。工作目录必须全新，四个可执行文件使用 CLI 版本，官方与父版源码及构建命令沿用 four_versions/build.sh。轮前和候选分别导出 builds.json 的 Git revision 下的 src 和 tools/build.py，调用 python tools/build.py --boost-prefix BOOST --output-dir OUT --jobs 2 生成。路径字面量会影响 ELF 散列，轮前的无差异重编可用 rebuild_stable.py NEW_DIR HISTORICAL_STABLE BOOST RECEIPT.json。各版本的搜索和发布默认值详见 builds.json。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
D=benchmarks/fork_structure
U=benchmarks/geometry_comparison/upstream/four_versions
python "$D/verify.py"
python "$U/offline.py"
mkdir -p "$WORK/inputs"
tar xf "$D/external/inputs.tar.xz" -C "$WORK/inputs"
python "$D/reprepare.py" "$WORK/inputs" "$WORK/prep-repeat"
python "$U/run_panel.py" candidate "$candidate" "$WORK/fixed" --jobs 4
python "$D/evidence.py" fixed "$WORK/fixed" "$WORK/fixed-tables"
python benchmarks/stage_local_structure/directions.py "$WORK/fixed-tables"
for arm in official parent stable candidate; do
  python "$D/holdout.py" run "$WORK/inputs" "$WORK/external" "$arm" "${!arm}"
  python "$D/timing.py" probes "$WORK/timing" "$arm" "${!arm}" --jobs 3
done
python "$D/evidence.py" external "$WORK/inputs" "$WORK/external" "$WORK/external-tables"
python benchmarks/stage_local_structure/directions.py "$WORK/external-tables"
python "$D/timing.py" freeze "$WORK/timing"
for stage in validation matched; do
  for arm in official parent stable candidate; do
    python "$D/timing.py" "$stage" "$WORK/timing" "$arm" "${!arm}"
  done
  if test "$stage" = validation; then python "$D/timing.py" eligible "$WORK/timing"; fi
done
python "$D/evidence.py" timing "$WORK/timing" "$WORK/timing-tables"
python "$D/evidence.py" check "$WORK/timing-tables"
```

输入可另用 prepare.py EMPTY_DIR 从 RCSB 原始沉积下载。对外部第二轮按相反引擎顺序重新 run，再调用 benchmarks/stage_local_structure/repeat.py 输入目录 第二轮运行目录 第一轮结果表 第二轮输出目录。重新运行会得到新的实际成本与可达分母。确定性核查 verify.py 使用已归档计时和原始坐标恢复本轮数字。
