# 独立结构的水点和侧链构象

提取工具用同靶其它 PDB 的实验坐标，输出可传给配体姿势评分的受体证据。查询晶体只使用指定链的蛋白原子。提取区域来自预先声明的蛋白锚点或独立盒，需在用户提供的候选姿势和参考配体坐标揭晓前确定。供体的物种、靶标、构建、晶体条件和来源需由清单编制者核对。

在引擎根目录用已有包含 NumPy、Gemmi 和 Biopython 的 Python 环境运行：

```bash
python tools/extract_independent_receptor_evidence.py MANIFEST.json --output-dir NEW_DIRECTORY
```

清单中的路径相对清单所在目录。主要字段如下。

```json
{
  "targets": [
    {
      "id": "mdm2",
      "accession": "Q00987",
      "query": {
        "pdb_id": "4HG7",
        "structure": "source/4HG7.cif.gz",
        "chain": "A",
        "receptor": "prepared/mdm2.pdbqt",
        "protein_anchors": [{"resseq": 61, "atom": "CA"}],
        "box_size": [20, 20, 20]
      },
      "donors": [
        {
          "pdb_id": "1T4E",
          "structure": "source/1T4E.cif.gz",
          "chain": "A",
          "accession": "Q00987",
          "condition": "Independent MDM2 crystal with a different inhibitor"
        }
      ],
      "sidechain_sites": [{"resseq": 61, "icode": ""}]
    }
  ]
}
```

链和编号使用 PDB 作者编号，插入码默认为空。受体 PDBQT 使用其它链字符时可声明 `query.receptor_chain`。`box_size` 默认每轴 16 Å，中心为蛋白锚点的坐标均值。外部独立盒可用 `query.box` 的 `center` 和 `size`，顺序为中心后边长的六项数组，或含 Vina 中心和边长字段的文本文件。外部盒需同时声明 `query.box_basis`。由查询配体描出的对接盒不适用于证据提取。省略 `sidechain_sites` 时，选择盒中心 6.5 Å 内最近的至多三个支持且拓扑完整的侧链，距离按各侧链最近重原子计算。用空数组禁用侧链提取。

## 配准与证据规则

受体 PDBQT 的同残基 Cα 通过 0.75 Å 内的唯一坐标匹配关联查询晶体，允许受体重编号。若坐标尚未在同一参考系，使用名称和作者编号可唯一匹配的骨架原子做刚体配准。至少匹配 32 个骨架原子，配准后 RMSD 不高于 0.25 Å。原始查询坐标及蛋白锚点随变换进入实际受体参考系。共编号且按坐标无法区分的原子不参与映射。同名氢原子用物理行和最近的同组重原子区分。

供体必须与查询 PDB 不同，同一 PDB 最多接受一次，且清单 accession 一致。使用蛋白序列全局比对，末端缺口不罚。接受阈值为对齐位点序列一致性至少 90%、相对较短序列覆盖至少 75%、至少 25 个同残基映射。取盒中心 12 Å 内同残基的共同 N、CA、C、O，至少八个残基，以 Kabsch 方法做刚体配准，骨架 RMSD 不超过 1.25 Å。并列替代构象取按平均占据率排序的单个非空 altloc，保留共享原子。规则、实际配准和拒绝原因可在输出中检查。

供体水来自其首模型的实验水氧，占据率至少 0.5，配准后在盒各边外扩 1.5 Å 范围内。按 PDB、占据率和源水标识固定处理顺序，以离最近固定簇种子 1 Å 内为合并规则，同一 PDB 在同一地点只提供一个独立观测。簇坐标是独立 PDB 水氧坐标的均值，同 PDB 重复水和各原始贡献坐标在来源文件中保留。单个供体也可提供候选，来源计数体现证据深度。

侧链观察需要同一供体包含所有选择位点的完整 χ 读数。距原生任一 χ 超过 35° 才认作替代，替代簇代表与成员各 χ 距不超过 30°。Asp、Glu、Phe、Tyr 的末端对称 χ 按 180° 比较，其余按 360° 比较。将有差别的 χ 迁入原生 PDBQT 拓扑，绕原生共价键转动远端重原子和归属的显式氢。输出逐行原子身份、顺序、类型、电荷和未动原子坐标保持原样。每个受体状态人口用独立 PDB 观测数加一后归一化，原生也有一个先验计数。这个相对权重受提供的结构集合影响。没有可区分替代时，原生人口为 1。

## 输出

- `waters.xyz` 每行是对齐到受体的水氧 `x y z` 埃坐标。文件为空时跳过水桥评分。
- `sidechains.json` 包含 `baseline_population`、`alternates` 和侧链来源，每个 `receptor` 相对 JSON 文件定位。
- `provenance.json` 记录查询坐标映射、实际盒、供体接受或拒绝、局部拟合变换、各水点及独立 PDB 计数。
- `summary.json` 汇总每靶输出。缺少可接受证据时生成空水点和原生受体权重。

水点可提供给 `--net_polar_waters`。该文件保留几何候选，实际水桥仍由构象评分阶段检测。 [引擎评分入口](NET_POLAR.md)
