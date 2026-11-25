# FastDock Vina

基于 [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) 的 CPU 优化版本，保留双精度计算和原有搜索设置，不使用 GPU。

## 相比原版的效果

下表记录发布版的历史结果。后续改动和复测见[工作线程复用](docs/BENCHMARKS.md#工作线程复用的本地复测)
和[受体查表与空间索引](docs/BENCHMARKS.md#受体查表与空间索引复用)。
近期对照见[刚性受体共享和建图复用](docs/BENCHMARKS.md#刚性受体共享和建图复用)、
[按使用需求生成网格](docs/BENCHMARKS.md#按使用需求生成网格)、
[按配体类型准备受体评分表](docs/BENCHMARKS.md#按配体类型准备受体评分表)、
[亲和力网格读取与姿势写出](docs/BENCHMARKS.md#亲和力网格读取与姿势写出)
和[局部搜索向量与网格索引](docs/BENCHMARKS.md#局部搜索向量与网格索引)。

在 Intel Xeon Platinum 8362 上测试 19 个完整对接案例，每个案例重复 3 次，取各案例耗时中位数后求和：

| 版本 | 总耗时 | 耗时降低 |
| --- | ---: | ---: |
| 原始 AutoDock Vina | 3242.525 秒 | — |
| FastDock Vina | 2989.231 秒 | **7.81%** |

全部 **25 项测试**（19 项完整对接 + 6 项评分／局部优化）的**等权几何平均加速比为 1.228×**。计算时，每项先取 3 次耗时中位数，以 `原版耗时 / FastDock 耗时` 得到加速比，再对 25 个比值取几何平均（[逐项数据](benchmarks/benchmark.csv)）。

两个指标的口径不同：**1.228×** 对每项测试赋予相同权重；**7.81%** 是 19 项完整对接的合计耗时降低，折合约 **1.085×**。几何平均值不代表每个案例的最低加速或单例最高加速。

19 个案例均快于原版，峰值内存（RSS）降低 **28.71%–73.27%**。已测配置的对接 PDBQT 输出与原版逐字节一致；结果仅适用于已测环境和输入。详细数据与验证范围见 [实验记录](docs/BENCHMARKS.md)。

![25 项测试的加速比及汇总指标](docs/figures/speedup-summary.png)

## 优化方法

- 按实际进入受体空间索引和配体相互作用列表的类型对构建 Vina/Vinardo 查表。普通相互作用按实际截断距离采样，大环闭合项保留 20 Å 查表。
- 显式受体评分只为当前配体和柔性残基需要的类型建表。亲和力网格复用这些表，批处理遇到新类型时补表并继续复用已有表。
- 组装模型时确定原子对的查表索引。AD4 只采样参与评分的原子对，每对保留独立的电荷信息。局部搜索并行插值成对能量和径向导数，保留双精度和原子对累加顺序。
- 同一网格行的相邻点复用受体候选和两个距离分量，并跳过超出整段截断范围的受体原子。大环闭合碳复用普通碳的网格文件。
- 建图时，同一受体类型和距离对应的各配体类型能量连续存放。显式评分和分子内评分直接从相邻能量样本求均值。
- 建图按受体原子遍历网格行，仅处理截断距离内连续的一段格点。各格点仍按受体原子顺序累加。
- 局部搜索复用相邻原子对共享的首原子坐标，在该组结束时写回其力。每组评分表只保存一次距离采样倍率。
- 同一搜索轨迹或精修循环复用 BFGS 的构象、梯度和矩阵缓冲区。局部优化直接使用最后一次求值生成的坐标，回退时恢复起始坐标。
- 构象树直接写入 BFGS 的连续梯度数组，试探构象直接从起始构象和搜索方向生成。最后一步省去后续不会使用的逆 Hessian 更新。
- 网格插值只计算一次起始角点索引，沿现有行和层跨度读取七个相邻采样值。
- 搜索线程和配体模型共享刚性受体。柔性受体的固定锚点需要重编号时，才复制相应的受体拓扑。
- 每次建图为每个受体原子查找一次横向格点分界，然后投影到每段网格行。保留每个格点的受体原子累加顺序。

## 编译

需要 Python 3.8+、支持 C++14 的 GCC/Clang 和 Boost（已验证 1.82）。

```bash
git clone https://github.com/Drangon-Li/fastdock-vina.git
cd fastdock-vina
```

Ubuntu / Debian：

```bash
sudo apt-get install build-essential python3 libboost-all-dev
python3 tools/build.py --jobs 4
```

macOS（已安装 Xcode Command Line Tools 和 Homebrew）：

```bash
brew install boost
python3 tools/build.py --boost-prefix "$(brew --prefix boost)" --jobs 4
```

## 运行

使用仓库自带的 1iep 示例：

```bash
mkdir -p demo_out
./build/vina \
  --receptor examples/basic/1iep_receptor.pdbqt \
  --ligand examples/basic/1iep_ligand.pdbqt \
  --config examples/basic/1iep_receptor.box.txt \
  --scoring vina --cpu 1 --seed 20261022 \
  --exhaustiveness 8 --max_evals 0 --num_modes 9 \
  --out demo_out/1iep_out.pdbqt
```

结果保存至 `demo_out/1iep_out.pdbqt`。更多参数见 `./build/vina --help`。示例使用 `exhaustiveness=8`，上述性能测试使用 `32`。

Vina 和 Vinardo 的普通相互作用截断距离默认是 8 Å。可用
`--interaction_cutoff 6.75` 覆盖，或在 `--config` 指向的文件中写
`interaction_cutoff = 6.75`。命令行显式值优先。取值大于 0 且最多 20 Å，
并统一用于普通评分势、受体及配体查表。自定义值用于直接读取受体并构建网格的运行。
预计算的 Vina/Vinardo 网格使用默认截断距离。AD4 和大环闭合项各使用其力场定义。

Vina 和 Vinardo 的 `--score_only` 可直接用刚性受体原子评分，省下亲和力网格的
计算和存储。`--ligand` 模式中省略六个盒参数时，固定姿势评分默认使用
`--autobox`。关闭显式评分（`--no_refine`）、默认参考能量涉及柔性残基、
以及搜索和导出时仍会构建所需网格。

从受体直接执行 Vina 或 Vinardo 的 `--batch` 时，各配体共享受体空间索引、
搜索盒和已建网格。首次用到某个配体或柔性残基的原子类型时才建图。批次带
`--write_maps` 时会一次导出完整的类型集合。`--maps` 加载预计算网格，
AD4 仍需对应的原子类型、电势和去溶剂化网格。导入会检查网格头部和采样值数量，
采样值须为有限数。姿势文件逐个写出模型，PDBQT 格式保持一致。

## 对接构象

![1iep 配体在受体中的位置及原版与 FastDock 的构象叠加](docs/figures/docking-overlay.png)

左图展示配体在受体中的位置，右图放大 MODEL 1：原版为蓝色，FastDock 为橙色，两者输出坐标完全重合（40 个原子的 RMSD 为 **0.000 Å**，未做结构对齐）。图中保留真实化学键，显示 37 个重原子。[对接输出](examples/basic/1iep_docked.pdbqt) · [图表数据与复现方法](docs/FIGURES.md)。

## 许可与致谢

遵循 [Apache 2.0](LICENSE)，保留 AutoDock Vina 的版权与署名。研究使用请引用 [Vina 原始论文](https://doi.org/10.1002/jcc.21334)和 [Vina 1.2.0 论文](https://doi.org/10.1021/acs.jcim.1c00203)。源码来源见 [PROVENANCE](docs/PROVENANCE.md)。
