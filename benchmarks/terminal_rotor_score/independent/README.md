# BRANCH 计数与实际转动的独立核验

C++ 观察器使用与对接相同的解析器和 `model::set`。对每个树扭转单独增加 1.234 rad，其余全部坐标保持初值，记录所有**化学重原子内部距离**的最大变化（Å）。这个量与整分子刚体运动无关，不含闭环伪原子。另一个实现只读取原始 PDBQT 的 BRANCH/ENDBRANCH 文本区间和每条原子记录，各自计算规范规定的双侧计数。前者检验实际运动，后者核对实现的拓扑含义。

- Galectin-3 94、因子 Xa 56、S6K1 14 个配体–接收结构输入，加上一个大环，共 165 条输入记录。两种实现的计数在每条记录中相同。部分接收结构重复使用同一配体制备。实时改变构象的轴数和双侧计数在这 164 条测定输入中也相同（用 10⁻¹⁰ Å 区分浮点零）。CHEMBL369480 在两个接收体中都得到 7 个双侧轴，7 个实际转动均改变了化学重原子之间的相对距离。
- 对照包括五种几何及每种的两种根方向，共十项。重原子侧链和刻意伸长的枢轴键都计 1，只有氢或水代理位点的末端计 0。将三种配体按正序或逆序组装，允许它们重用相同的 PDBQT 原子编号，两个结果均为三配体合计 2 个轴。
- 一个严格共线的末端重原子构型，两种根方向都计 1，而任何化学重原子内部距离都没有改变。环十二烷用十一个 BRANCH 表示开环后的树，两个末端转动只移动闭环代理点，化学重原子内部距离保持不变。因此此描述符统计**两侧都有参考位点的声明轴**。它不能给出独立构象空间的维数或热力学构象熵。

[逐输入计数与全部扭转的有限移动](results/prepared.csv) · [定向和反向对照](results/controls.csv)

用仓库的 `tools/build.py` 编译正式代码。下面的 `OBJ` 为该构建的 obj 目录，Boost 安装在系统默认搜索路径；自定义安装需沿用构建时的 include、lib 和 rpath 参数。Python 依赖与同目录上一级 assay 复算相同。

```bash
D=benchmarks/terminal_rotor_score/independent
objects=()
for file in "$OBJ"/*.o; do
    case "$(basename "$file")" in main.o|split.o) continue ;; esac
    objects+=("$file")
done
c++ -std=c++14 -DNDEBUG -O2 -pthread -Isrc/lib "$D/axis_check.cpp" "${objects[@]}" \
    -lboost_filesystem -lboost_program_options -lboost_serialization -lboost_thread -o "$WORK/axis-check"
python "$D/verify.py" "$WORK/axis-check" "$WORK/results"
```
