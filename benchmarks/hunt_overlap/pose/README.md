# 四组人源晶体姿势对照

原版 Vina、FastDock 父版本和两个搜索变体共运行 384 次盲对接。每组的一个单元对应同一查询、接收受体、种子与预算。搜索只读取受体、从 CCD 化学图生成的配体和受体已知口袋中心。共晶坐标用于定义口袋中心和事后评分，未用作初始构象。三个比较中的开发体系、独立靶点和同受体复跑分别统计。

| 实验与靶点 | 每臂单元 | 上游 Top-1 / 返回池 Oracle | FastDock 父版 | 主干预 | 消融 |
|---|---:|---:|---:|---:|---:|
| 切线延拓，THRB / BRD9 | 32 | 8 / 8 | 8 / 8 | 全部 hunt 4 / 4 | 首次 hunt 7 / 7 |
| 首次平移，THRB / BRD9 开发组 | 32 | 8 / 8 | 8 / 8 | 每次首次提议 5 / 5 | 仅修正越界首次提议 5 / 5 |
| 首次平移，NR1H2 独立靶点 | 16 | 1 / 2 | 0 / 0 | 2 / 2 | 2 / 2 |
| 新种子复核，NR1H2 / DPP4 | 16 | 0 / 2 | 0 / 0 | 0 / 0 | 0 / 0 |

以上各组 Top-5 均等于返回池 Oracle。每组的 summary.csv 列出按靶点、查询、预算和重对接/交叉对接拆分的结果，cells.csv 与 modes.csv 提供单元及逐姿势数据。成功线为受体坐标系内、考虑化学图对称性的重原子 RMSD ≤2 Å。每个查询使用两个接受受体和两档预算，同一种子在小预算和大预算下共享轨迹前缀。计数不能当作独立复合物样本量。每臂 1 线程、1 条轨迹、最多 9 姿势、24 Å 立方盒、3000 或 12000 搜索上限。实际联合能量–梯度求值为 hunt、第二次局优和显式精修的计数之和。

## 真实原因与反例

- **切线延拓。** 32 个单元内全程切线新增 1 次命中、丢失 5 次，均发生在 THRB 的 12000 上限，净减 4。此预算的实际求值为父版 196284 次、切线 196246 次，总用时分别为 2.558 和 2.574 秒。1NAX 重对接的一次首姿从 7.93 降至 0.99 Å，两次 3IMY → 1NAX 则从 0.50、0.77 增至 8.64、8.62 Å。仅首次 hunt 也净减一次。三个 FastDock 臂在 Monte Carlo 结束保存池、显式精修前后、返回池以及首位的命中次数逐单元一致；轨迹内更早出现过的构象未全部保存。BRD9 的 16 个父版单元无命中。将晶体扭转移植到原始查询的刚性几何后，四个 BRD9 方向的局优仍在 0.39 Å 以内，且每个单元的近天然构象比盲对接返回池最低的 inter+intra 有利 1.71–7.66 kcal/mol。对同一查询，共用分子内参考和相同扭转分母时，此量保留原生排序顺序。将投影近天然姿势单独放入各个固定池，BRD9 的 16 次排序均为第一。
- **首次平移。** 开发面板中，每次映射和只修正越界的两个臂在全部 32 个单元产生字节相同的输出。两臂在 BRD9 增加一次 5EU1 重对接命中（8.30 → 0.64 Å），在 THRB 从 8/16 降至 4/16。NR1H2 上的两次新命中来自同一 1P8D 配体、同一受体、同一种子 625237。换到 704521 后没有复现。在最初的种子，上游原生 Vina 两档预算都找到了正确模式，大预算时排在第二位。新种子上游仍有两次返回池命中，两个改动臂和 FastDock 父版均为零。
- **实际时间匹配。** 在 NR1H2 的种子 625231 上按查询–受体方向，以三个交错计时的中位数对候选预算作一次线性插值，目标为父版 12000，候选预算限制在 3000–12000。再对原实验另一种子 625237 测三次。四个方向各自时间中位数之和分别为父版 0.824 秒、全量首次映射 0.812 秒、仅修正越界 0.791 秒，实际求值分别为 48991、43803、40746，首姿命中分别为 0、1、1。具体到唯一命中方向，三者是 0.173、0.203、0.164 秒；只修正越界的方法在上限降至 8278 后保住命中。另一个方向实际超出目标的情况保留在 matched_trials.csv。这项追加检查用于排除一次命中完全来自多花时间的解释。独立种子仍然失败。
- **难例。** NR1H2 两个交叉方向的晶体配体映射到接收口袋后分别有 6、9 对严重碰撞。DPP4 两个交叉方向有 9、19 对，最近距离 1.06、0.99 Å。两组的全部重对接参考构象均无此类碰撞。DPP4 的两种输入配体比较刚性结构差异较大，两个原配体中心在对齐后相距 6.85 Å。3VJK → 4A5S 的参考姿势有一个重原子超出原盒 0.35 Å，因此保留其盲对接结果，省略无法在原盒运行的晶体起始局优。晶体诊断中，3VJK 重对接采用四种等价原子映射中初始投影误差最小的一种，精修后为 0.43 Å。使用第一种映射会得到 11.65 Å。相同查询在 12000 预算下，近天然姿势的 inter+intra 为 −10.792 kcal/mol，低于 2 Å 的姿势未进入返回池；返回池最好 RMSD 为 9.48 Å，但最优能量为 −10.949 kcal/mol。将同一个投影近天然姿势追加进固定父版返回池，按同一查询参考重新排序，3000 预算为 1/7，12000 预算为 2/5。两种名次在打印舍入允许的 ±0.002 kcal/mol 范围内均无并列。这是扩大搜索即使生成该近天然姿势也可能仍排错的具体反例。

立体几何审计独立重建 PDB/CCD、手性图匹配和受体局部坐标变换。各组分别核对 679、597、223、258 个输出模式，最大 RMSD 复算误差小于 2.5×10⁻¹⁴ Å，无已赋值手性标记变化，最大键长误差小于 0.0014 Å。严重碰撞指重原子距离小于两者范德华半径和的 0.75 倍。这些检查不包含结合态自由能、质子化平衡、所有键角或分子内应变。研究沿用 CCD 微态，THRB 的羧酸和 BRD9 的胺均以形式电荷为零的图制备。BRD9 存在距配体 5.39–5.94 Å 的晶格接触。每个结构的来源、原始文件哈希、骨架、分辨率、最近省略/补建原子和制备耗时见相应的 sources.json。首批共同制备 1.63–2.51 秒/结构，独立 NR1H2 2.76–2.98 秒/结构，新增 DPP4 7.65–7.83 秒/结构，运行时复用。盲对接的每次进程时间含地图、搜索、精修及写出。计数回放、晶体诊断、独立审计和观察性跟踪的开销单独记录。

## 复现

在仓库根目录执行。使用 Linux、C++14、Boost、Python 3.12、RDKit 2025.9.3、Meeko 0.7.1、Gemmi 0.7.4、PDBFixer 1.12.0、OpenMM 8.6.1、NumPy 2.3.5、SciPy 1.17.0、Biopython 1.88。计时需要可用的单个 CPU 核和同机串行运行各臂；替换下列 CPU 编号和 Boost 路径。上游可执行程序由 3c65c0b 构建。三个独立结构/种子面板的 64 个上游输出分别与计数回放、官方 v1.2.7 发行版逐字节相同，父版也与独立构建的 f04bd7f 对照逐字节相同。重复的开发集共用同一初始输入和参数。

```sh
P=benchmarks/hunt_overlap/pose
W=$(mktemp -d)
C=2
BOOST=/path/to/boost
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
python benchmarks/pair_gradient/build_arms.py "$W/up" --arms upstream upstream-counts --boost-prefix "$BOOST" --jobs 2
for pair in parent:f04bd7f6f555cf3b71d28e259552cf0740dc51a4 tangent:4ba83b5087452639de42502bc0194b6915363c26 placement:6e2513f213e23db4b0da41ac5aa89945f882ac6b; do
    name=${pair%:*}
    rev=${pair#*:}
    mkdir "$W/$name-src"
    git archive "$rev" src tools/build.py | tar -x -C "$W/$name-src"
    python "$W/$name-src/tools/build.py" --output-dir "$W/$name-build" --boost-prefix "$BOOST" --jobs 2
done
for study in tangent placement-dev placement-holdout placement-confirm; do
    dir="$W/$study"
    mkdir "$dir"
    case "$study" in
        tangent|placement-dev) input=tangent ;;
        *) input=placement-holdout ;;
    esac
    tar -xJf "$P/$input/inputs.tar.xz" -C "$dir"
    if [ "$study" = placement-confirm ]; then
        tar -xJf "$P/$study/dpp4_inputs.tar.xz" -C "$dir"
    fi
    variant=placement
    if [ "$study" = tangent ]; then variant=tangent; fi
    python "$P/panel.py" run "$dir" --study "$study" --cpu "$C" --core "$W/$variant-build/vina" --upstream "$W/up/upstream/build/vina" --counts "$W/up/upstream-counts/build/vina"
    python "$P/panel.py" verify "$dir" --study "$study"
    python "$P/panel.py" crystal "$dir" --study "$study" --pristine "$W/parent-build/vina"
    python "$P/summarize.py" "$dir"
done
python "$P/match_time.py" "$W/placement-holdout" "$W/placement-build/vina" --cpu "$C"
```

用同一目录运行 `panel.py prepare WORK --study STUDY` 可从原始输入重新制备，运行 `panel.py parity WORK --study STUDY --release OFFICIAL_VINA --pristine PARENT --cpu CPU` 可复核发行版坐标。使用仓库已有的被动跟踪补丁复现每一阶段。

```sh
python benchmarks/accepted_retention/core/build_stage.py "$W/placement-src" "$W/placement-build" "$W/observed" --boost-prefix "$BOOST"
python "$P/stages.py" "$W/placement-confirm" "$W/observed/vina" --study placement-confirm --cpu "$C"
```

切线实验使用相应的 tangent-src 与 tangent-build。每次被动回放均核对主实验坐标的 SHA-256 以及各环节实际计数。runs.tar.xz 包含盲对接坐标、日志和计时，trace.tar.xz 包含阶段回放。placement-confirm 的输入分别解压自 placement-holdout/inputs.tar.xz 与本目录 dpp4_inputs.tar.xz。仅验证已提交数据时，按上面的方法解压输入，再解压对应的 runs.tar.xz，复制 modes.csv，执行 `panel.py verify WORK --study STUDY` 即可。
