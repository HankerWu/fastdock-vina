# 配体准备耗时和配对表等价性

微基准覆盖 1IEP 伊马替尼、6WTO 巴瑞替尼和四个大环。全部文本在计时前读入内存。每个进程先预热八次，然后重复 192 次完整 PDBQT 解析、成键推断、原子分型和配对表构造。记录单进程 CPU 时间和壁钟时间。三个程序的六种运行顺序各用一次，各轮轮换配体次序。表中为六次的 CPU 中位数，单位为每次解析的微秒。一次生产对接一般只需准备一次配体，这些数值不能直接外推为整程加速倍数。

| 分子 | 父版 | 仅改为广度优先 | 广度优先与闭环位掩码 |
|---|---:|---:|---:|
| 1IEP | 148.65 | 150.89 | 85.25 |
| 巴瑞替尼 | 81.29 | 80.67 | 58.12 |
| 5QCT | 187.88 | 191.77 | 111.02 |
| 1YET | 196.75 | 207.12 | 109.87 |
| 1FKB | 460.63 | 451.88 | 201.61 |
| 1FKF | 347.05 | 356.70 | 169.39 |

闭环编号 k 在原子 i 的位掩码中出现，当且仅当 i 或其一个直接非网格邻居带有 CGk。因此两端掩码相交与原程序逐对检查两个邻域是否出现相同 CGk 等价。匹配的 Gk–CGk 闭环对继续在判断前得到豁免。无标签的整个配体可跳过闭环邻域查询。

在 XS、AD 两种表示下，分别对未标记的两个普通配体和四个大环的 G0–G3 编码比较全部有序输出，涵盖相互作用分类、两个端点和查表索引。36 组输入上，仅改遍历算法与加入掩码两臂的字节串完全相同。巴瑞替尼 XS 表由 208 个减为 202 个、AD 表由 224 个减为 218 个，差异均由广度优先修正产生。每一个表的 SHA-256 和行数在 tables.csv。timing.csv 含所有 108 个单独时间观测。不同运行顺序、同分子的四个等价标签均属于技术重复。

## 复现

先按上级目录的说明将四个大环输入准备到同一目录下，包含 bace/raw/single-0.pdbqt、hsp/1YET/ligand.pdbqt 和 extra 下的 1FKB、1FKF 两个子目录。以下命令从仓库根目录执行，使用系统 Boost。实测编译器为 GCC 14.2，Boost 为 1.83。自定义 Boost 安装另需向构建工具传入对应的 boost-prefix，编译探针时加上相应头文件、链接和运行时搜索路径。

```bash
REPO=$(pwd)
W=$(mktemp -d)
P="$REPO/benchmarks/bond_distance/preparation"
git worktree add --detach "$W/parent" 36636171ae6aef42471b0de5d6b7a90b5f90732e
git worktree add --detach "$W/bfs" 36636171ae6aef42471b0de5d6b7a90b5f90732e
git -C "$W/bfs" apply "$P/bfs_only.patch"
for arm in parent bfs release; do
  tree="$W/$arm"
  if [ "$arm" = release ]; then tree="$REPO"; fi
  python3 "$tree/tools/build.py" --output-dir "$W/bin-$arm"
  c++ -O3 -DNDEBUG -std=c++14 -pthread -I"$tree/src/lib" "$P/probe.cpp" \
    $(find "$W/bin-$arm/obj" -name '*.o' ! -name main.o) \
    -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem \
    -lboost_program_options -o "$W/probe-$arm"
done
python3 "$P/compare.py" --repo "$REPO" --macro-inputs "$MACRO_INPUTS" \
  --parent "$W/probe-parent" --bfs "$W/probe-bfs" --release "$W/probe-release" --out "$W/results"
```
