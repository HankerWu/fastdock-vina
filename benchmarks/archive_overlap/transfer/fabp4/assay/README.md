# 人 FABP4 的五个同测定配体

在共同的两个受体上，上游 Vina、原 Pro 与在线档案修正版均正确区分 8/10 对配体，Spearman ρ 为 0.80，Kendall τ-b 为 0.60。八对实测解离常数相差至少三倍的配体全部排对。该系列中，单按重原子数排序也能排对全部八对，并能排对九个大小不同的配体对。仅有的同重原子数且同形式电荷对为布洛芬的 R/S 对映异构体，三个算法都将它们排反。

González 与 Fisher 的实验在 298 K、0.1 M 柠檬酸钠 pH 6.5、2% DMSO 下，对 1 µM 人 FABP4 进行 ANS 荧光竞争测定，以同一竞争模型拟合配体的解离常数 K_L。表 1 给出数值，表 2 给出晶体配体编号，图 1 及方法部分说明条件。每个化合物仅纳入一次，包括两种绝对构型明确的布洛芬。measured.csv 保留原文报道的误差，未将其解释为独立重复数或置信区间。[原论文，表 1、表 2 和图 1](https://doi.org/10.1107/S2053230X14027897)

| 化合物 | PDB / CCD | 解离常数 µM | 重原子数 | 父版首姿平均分 | 档案修正后 |
|---|---|---:|---:|---:|---:|
| L2 | 3P6D / ZGB | 0.737 ± 0.046 | 14 | −6.542 | −6.542 |
| L3 | 3P6E / ZGC | 73.0 ± 7.3 | 13 | −6.1555 | −6.1555 |
| L4 | 3P6F / FBZ | 148 ± 7.5 | 12 | −6.16325 | −6.16325 |
| L5，R-布洛芬 | 3P6G / IZP | 0.133 ± 0.0066 | 15 | −7.202 | −7.202 |
| L6，S-布洛芬 | 3P6H / IBP | 0.0949 ± 0.0059 | 15 | −6.847 | −6.80275 |

L2 与 L3 的解离常数相差 99.05 倍。按 RT ln(K_L3/K_L2) 计算，实验相对自由能差为 +2.721 kcal/mol，三个算法的分差均为 +0.3865 kcal/mol。两种布洛芬的实测顺序为 S 更强，差约 0.200 kcal/mol。父版和修正版分别预测 R 更强 0.355 和 0.399 kcal/mol。另一对误判来自 L3/L4（实测差约二倍，平均分差仅 0.00775 kcal/mol）。十对相对自由能差的平均绝对误差从父版 2.006 变为 2.024 kcal/mol，未作重新标定。十个分子对共用五个化合物，不能作为十项独立证据。

## 归因与成本

所有配体均使用同一制备流程、单一中性羧酸微态和两个原有受体 3P6D/E。每个分子共享两个种子和对应受体口袋中心，四条轨迹、单核、12000 名义求值/轨迹。L2/L3 复用了最先冻结的结构试验输出。L4–L6 从各自 CCD 化学图独立生成初始构象，核验沉积的手性映射。其共晶配体坐标仅用于离线识别与核验，无额外受体或配体初始姿势信息送入预测。主终点为两个受体、两枚种子上首姿分数的等权平均。两个受体分别得到 8/10 和 9/10，不能事后择优。该测定是对已知 FABP4 模型失误的回顾性验证，尚不足以代表跨家族亲和力泛化。

在原父版的 20 个首姿上，冻结所有原子坐标及该运行的游离态参考，分别运行三个程序的 score_only。20 组的总分、分子间项、内部项、扭转项和游离态参考逐项相同（CLI 输出精度 0.001 kcal/mol）。解锁搜索后，只有 L6 的三个条件改变了首姿分数，其中一项更低、两项更高。归档变化可以改变选中的姿势，在这个同测定系列中没有带来额外亲和力排序成功。

每臂 20 个预测的实际联合求值为上游 968086、父版 968082、修正版 968177，对应单核墙钟总计 15.544、3.581、3.501 秒。每项计数回放均核对完整 PDBQT 字节。该计时合并先前冻结条件与此次补充条件，不能从一次短程总时断言速度提升。sources.json 含三个新配体的图、原始来源、输入哈希、绝对构型和单次制备用时。其余两种化合物的输入与耗时见结构试验父目录。

## 复现

使用 transfer/verification.json 的依赖，在本目录执行。original 为先前结构面板；work 为本测定的三个新增查询，两者共享完全相同的两个受体。归档包含结构、输入、每次 CLI 的输出与计数、固定坐标及打分日志。下面的校验会从每个原始 PDBQT 读取首姿能量，重新汇总四条轨迹及精修计数，并检查三个程序的固定坐标结果。

```bash
sha256sum -c SHA256SUMS
mkdir -p original work
tar -xf ../../inputs.tar.xz -C original
tar -xf ../../outputs.tar.xz -C original
cp ../../cells.csv ../../modes.csv original/
tar -xf new_inputs.tar.xz -C work
tar -xf new_outputs.tar.xz -C work
cp -r original/3P6D original/3P6E work/
python verify.py original work .
```

重新运行三种算法的三种新增配体，再生成五配体的完整汇总。旧配体直接读取归档，避免重复求值。CPU 4 可换为系统允许绑定的核心。新配体的原始制备可用 `python assay.py prepare work --original original` 复算。

```bash
OPENBLAS_NUM_THREADS=1 python assay.py dock work --cpu 4 --upstream "$UPSTREAM" \
  --counter "$UPSTREAM_COUNT" --parent "$PARENT" --candidate "$ONLINE_FIXED"
OPENBLAS_NUM_THREADS=1 python summarize.py --work work --original original --out analysis \
  --fixed --cpu 4 --upstream "$UPSTREAM" --parent "$PARENT" --candidate "$ONLINE_FIXED"
```

## 最终发布程序的复验

在加入最终精修候选选择修正后，对全部 20 个配体–受体–种子条件独立重新执行。最终程序与在线档案修正版的 20 份输出 PDBQT 均逐字节一致，包括每一名姿势、首姿能量和所有坐标。六项调用/求值计数及显式精修输入数也全部一致。这一亲和力面板没有新增或删去低排名姿势。release_outputs.tar.xz 保留最终程序的实际资源使用和计数，release_cells.csv 给出逐条件核验结果。

```bash
OPENBLAS_NUM_THREADS=1 python verify_release.py --work final_run --original original \
  --new work --release "$FINAL_RELEASE" --cpu 4
```
