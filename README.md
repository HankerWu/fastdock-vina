# FastDock Vina

基于 [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) 的 CPU 优化版本，保留双精度计算和原有搜索设置，不使用 GPU。

## 相比原版的效果

下表记录发布版的历史结果。后续改动和复测见[工作线程复用](docs/BENCHMARKS.md#工作线程复用的本地复测)
和[受体查表与空间索引](docs/BENCHMARKS.md#受体查表与空间索引复用)。
近期的本地对照见[按需评分查表与网格构建](docs/BENCHMARKS.md#按需评分查表与网格构建)。

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
- 建图循环直接读各类型对的连续评分数组，沿用原子顺序与距离计算。
- 缓存热循环中的查表指针，通过移动对象避免深拷贝。
- 调整 BFGS 矩阵布局，使用 CPU SIMD 同时处理独立行，保留原有累加顺序。

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

## 对接构象

![1iep 配体在受体中的位置及原版与 FastDock 的构象叠加](docs/figures/docking-overlay.png)

左图展示配体在受体中的位置，右图放大 MODEL 1：原版为蓝色，FastDock 为橙色，两者输出坐标完全重合（40 个原子的 RMSD 为 **0.000 Å**，未做结构对齐）。图中保留真实化学键，显示 37 个重原子。[对接输出](examples/basic/1iep_docked.pdbqt) · [图表数据与复现方法](docs/FIGURES.md)。

## 许可与致谢

遵循 [Apache 2.0](LICENSE)，保留 AutoDock Vina 的版权与署名。研究使用请引用 [Vina 原始论文](https://doi.org/10.1002/jcc.21334)和 [Vina 1.2.0 论文](https://doi.org/10.1021/acs.jcim.1c00203)。源码来源见 [PROVENANCE](docs/PROVENANCE.md)。
