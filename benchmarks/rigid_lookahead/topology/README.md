# 精修起点的 PDBQT 拓扑复核

粗搜索会把起点写为 0.001 Å 精度的 PDBQT。显式精修从文件重新读入配体并推断键，因此某些接近成键阈值的输入可能因坐标舍入改变评分拓扑。用 **2e610c4 的实际解析器**逐一比较每份起点与其独立准备的配体。从开发面板取 2,916 份起点，从 PGR / PPARα 预算复核取 1,944 份，从 RBP4 新结构面板取 216 份，共 **5,076 份，拓扑差异为零**。前两项覆盖每臂的全部 9、24 和 48 份额及两颗种子。跨份额重复出现的起点也计入，其含义是检查所有精修输入。

按 PDBQT 原子序号比较每个原子的元素、AD 和 XS 类型、部分电荷、共价键的两个端点与每条可动键、BRANCH 树和 TORSDOF。对内部、闭环、分子间与其他相互作用配对分别比较端点和势能表索引，保留重复条目。共同共价边上，舍入带来的最大键长差为 0.001444 Å。附带的人工对照将 5EPS 的供氢原子单独拉离重原子，解析器按预期报告共价边及类型差异。

该复核检验准备图在这批重载起点中是否发生改变。它没有独立测定准备图本身的化学正确性，也不检验有限坐标精度对连续势能和梯度造成的微小变化。[逐配体、份额与方法核对](summary.csv) · [不匹配清单](mismatches.csv)

如已运行 replay_budget.py 与 transfer/replay.py，可直接使用它们的输出。也可仅解包已有结果进行复核。以下为 Bash 命令，BASE 与 BUILD 分别为 2e610c4 源码目录和对应 tools/build.py 生成的目录。若 Boost 位于非标准位置，相应加入头文件和库路径。

~~~bash
D="$PWD/benchmarks/rigid_lookahead"
W=$(mktemp -d)
mkdir -p "$W"/{dev,secondary,holdout,sweep}
tar -xJf benchmarks/hydration_transfer/inputs.tar.xz -C "$W/dev"
tar -xJf benchmarks/search_backtracking/pose/inputs.tar.xz -C "$W/secondary"
tar -xJf "$D/transfer/inputs.tar.xz" -C "$W/holdout"
tar -xJf "$D/transfer/outputs.tar.xz" -C "$W"
cp "$D/transfer/results/poses.csv" "$W/transfer-rigid/"
for stage in dev secondary; do
    for n in 9 24 48; do
        tar -xJf "$D/results/${stage}_${n}.tar.xz" -C "$W/sweep"
        for arm in raw sham gradient; do
            mv "$W/sweep/$arm" "$W/sweep/${stage}_${n}_${arm}"
        done
    done
done
c++ -DNDEBUG -O2 -pthread -std=c++14 -I"$BASE/src/lib" "$D/topology/probe.cpp" \
  "$BUILD"/obj/{coords,model,parse_pdbqt,quaternion,random,utils}.o \
  -lboost_thread -lboost_filesystem -lboost_serialization -lboost_program_options -lboost_system \
  -o "$W/probe"
python "$D/topology/check.py" "$W/probe" "$W/dev" "$W/secondary" "$W/sweep" \
    "$W/holdout" "$W/transfer-rigid" "$W/check"
diff "$W/check/summary.csv" "$D/topology/summary.csv"
diff "$W/check/mismatches.csv" "$D/topology/mismatches.csv"
~~~
