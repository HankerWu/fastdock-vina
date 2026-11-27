# 新大环复合物的验证与复核

最初登记的 5 个复合物都未通过输入准备。在没有读取对接结果的情况下，随后登记另一组 5 个。10 个样本及各自的失败原因、实验来源、坐标和 CCD 哈希、准备耗时列在 `sources.json`。SRC 备选 8K79 的实际物种为鸡，且受体制备失败。最终纳入两个不同人源靶点，各仅有一个配体。FXIa 采用 C500S 变体。三个种子衡量同一分子内的轨迹变化。

| 蛋白–配体 | PDB / CCD | 分辨率 | 闭环大小 | 最近晶格接触 |
|---|---|---:|---:|---:|
| Cyclophilin A 与简化 sanglifehrin 配体 | 5T9W / 78E | 2.00 Å | 22 | 3.32 Å |
| 凝血因子 XIa 与 milvexian | 7MBO / YXG | 0.92 Å | 12 | 3.31 Å |

采用 PDB 中 A 链，所有算法共享去水受体和以已知口袋为中心的 26 Å 正方体。制备保留了同一化学微态。配体从 CCD SMILES 单次 ETKDGv3 和 MMFF 优化生成，再由 Meeko 拆环。两例均存在靠近晶体配体的晶格邻居，独立采用单链受体时省略了这些接触。计算时间均为完整 CLI 启动到退出。预处理由各臂共享，`sources.json` 的耗时从配体/受体准备开始，原始 PDB 预下载耗时未计，下载 CCD 及并行准备间的竞争可能计入。数据版本在 `software.json` 中。

`paired` 是上游 1.2.7、父版、首步环扭转恢复及单次坐标下降的原始对照，`lower` 是后两臂以 2000/9000 上限复测的成本对照。12 个靶点–种子–预算组合中，每一臂的 Top‑1、Top‑5 和输出池 2 Å 成功数均为零。两个配体的完整协调验证在读取最初实验终点后参加了持续固定环的 `persistent` 复核。各臂的重复计时产生相同完整输出及六个计数器。低预算不一定等时，详见 `paired_cost.csv` 的逐单元实际比值。

对每个姿势，第一套 RMSD 使用冻结的 CCD 图和晶体 PDB 重原子对应。第二套使用 RCSB ModelServer 提供的结合态 SDF、自其中构建的带立体化学的配体图，以及 Meeko 写出的 SMILES IDX。两者均在受体坐标系中枚举对称映射。各个输出和中间阶段的两次计算在 CSV 的六位小数精度上完全相同。闭环两端在原始分子图中必须相邻。碰撞计数为距离任一受体重原子不足 2 Å 的配体重原子数。

`stages` 记录 3000 档的 18 条 Pro 轨迹。534 个阶段构象中没有 2 Å 成功；其中最接近实验位置的 hunt、第二次局优或保存构象为 4.973 Å。`persistent_stages` 记录六条持续固定环轨迹的 310 个阶段构象，最接近实验位置的构象为 4.817 Å。该臂的 raw/proposal 均为环已固定的第一步扰动状态。观察器输出不计入正式计时。每条轨迹都与无观察器二进制得到的完整最终输出和六个计数器逐一核对。阶段数据的 `rank` 表示记录顺序；中间阶段应读取整段最小 RMSD，而非把记录顺序看作预测排名。固定输出池按 3 Å 显著断环阈值排序的结果见 `paired/fixed_pools.csv`。48 个池中有 3 个换首姿，4 个池所有姿势均超阈值，没有新增晶体姿势成功。

## 复现

在本目录执行以下命令。依赖版本见 `software.json`。第一个归档包含下载得到的 PDB/CCD、独立结合态 SDF、入选样本的冻结 PDBQT/元数据及被筛除样本的记录。第二个归档包含每次运行的完整输出及日志，也包含阶段轨迹。解压不会运行新的对接。

```bash
sha256sum -c SHA256SUMS
mkdir -p replay
tar -xJf inputs.tar.xz -C replay
tar -xJf outputs.tar.xz
python external_score.py replay paired
python external_audit.py paired repeat
python external_cost.py lower paired repeat lower/paired_cost.csv
python external_score.py replay persistent
python external_audit.py persistent persistent_repeat
python external_cost.py persistent persistent persistent_repeat persistent/paired_cost.csv --low-repeat persistent_repeat
```

从原始化学图重新准备输入时，使用 `python prepare_panel.py REPOSITORY_PATH NEW_WORK` 和 `python manifest.py REPOSITORY_PATH NEW_WORK`。如果重新下载的 ModelServer SDF 有不同的生成时间元数据，可直接使用归档坐标复算。以下命令展示完整对接、降低上限和阶段回放；环境变量选择本实验相应的构建产物。输出写入新的空目录。`--cpu` 选择操作系统允许绑核的逻辑 CPU。原生计数二进制只读回放，需与未插桩上游产生相同的输出字节。

```bash
args=(--work replay --cpu 2 --upstream "$UPSTREAM" --upstream-counts "$UPSTREAM_COUNTS" --parent "$PARENT" --reset "$RESET" --ccd "$CCD")
python external_run.py "${args[@]}" --mode paired --output replay-paired
python external_run.py "${args[@]}" --mode repeat --output replay-repeat
python external_run.py "${args[@]}" --mode lower --output replay-lower
python external_run.py "${args[@]}" --lock "$LOCK" --mode persistent --output replay-persistent
python external_run.py "${args[@]}" --lock "$LOCK" --mode persistent --output replay-persistent-repeat
python external_stages.py replay replay-paired replay-stages --parent "$OBS_PARENT" --reset "$OBS_RESET" --ccd "$OBS_CCD"
python external_stages.py replay replay-persistent replay-lock-stages --lock "$OBS_LOCK"
```
