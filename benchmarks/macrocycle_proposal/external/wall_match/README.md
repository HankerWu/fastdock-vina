# 以独立计时种子校准固定环的预算

父版 max_evals = 12000。先用独立随机种子 281291，仅依据完整 CLI 的耗时与实际求值选择两蛋白各自的固定环上限，再在之前的三枚种子上做顺序平衡的四次运行。每次使用完全相同的受体、配体、搜索盒和 CPU。每个靶点只有一个配体，两个上限用于使实验计时可比。部署时不能假设相同的上限也适合别的配体、机器或输入规模。

`calibration/selection.json` 是读取验证终点前冻结的选择，`calibration/costs.csv` 保留所有试探预算，`validation/paired.csv` 给出六个验证单元的双方实际求值、时间和晶体姿势误差。`validation/all_runs.csv` 保留全四轮时间及计数，`validation/models.csv` 保留每个不同输出池的全部姿势。计时排除归档和解析阶段。协议与选择上限的依据见 `PROTOCOL.md`。这些验证种子曾用于正式实验，此处用于核验成本和已有方法的表现。校准进程和验证进程均为实验成本，使用者可以从各自的计时表重算。

按主目录说明准备相同的 PDBQT、晶体与 SDF 文件，生成父版及固定环二进制，设置 INPUTS、PARENT 和 LOCK 环境变量。下列命令以本目录为工作目录，`../external_score.py` 引用上级外部评测目录的评分程序。

```bash
sha256sum -c SHA256SUMS
tar -xJf raw.tar.xz
python ../external_score.py "$INPUTS" validation
python describe.py validation
python matched_time.py calibrate "$INPUTS" rerun-calibration "$PARENT" "$LOCK" --cpu 2
python matched_time.py validate "$INPUTS" rerun-validation "$PARENT" "$LOCK" --selection rerun-calibration --cpu 2
python ../external_score.py "$INPUTS" rerun-validation
python describe.py rerun-validation
```

使用本机可绑定的核心替换 2。新的校准可受机器负载影响，若要用此实验归档的上限复核，在二进制、输入、协议和核心编号一致的环境中，用 `--selection calibration` 指定原始选择记录。全部验证的重复只用于计时，每个算法、种子及预算仅对应一个确定性输出池。
