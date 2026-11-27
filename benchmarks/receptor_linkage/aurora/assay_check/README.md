# AURKA 三个预定配体的独立测定核对

Davis 等的原始 Supplementary Table 1 登记 AURKA (NP_003591.2)。Table 4 的 Imatinib 列在 AURKA 行为空，未给出精确 Kd。AT9283 和 danusertib 均未列入这张表，PHA-665752 是另一化合物。因此，三个固定配体在这套独立实验中可用的精确 Kd 和可计分三倍差异对均为零。空白未当作 10,000 nM，也未从另一个实验借用参照值。[原测定](https://doi.org/10.1038/nbt.1990)

来源核对另发现，BindingDB 的 AT-9283 人 AURKA 页面可见 3 和 26 nM 两条 IC50；该化合物其他条目中有 JAK2 的 Kd，均不满足本题终点和靶点。Klaeger 等使用细胞系裂解液 Kinobeads 剂量–反应法给出的量为表观 Kd，具有不同的实验定义。Lake 等在荧光标记的激酶及多个磷酸化／TPX2 条件中研究了 danusertib，本次未核实同状态下另一冻结配体的精确 Kd。检索范围内尚无可独立计分的追加配对。这不排除另有公开记录。[AT-9283 的靶点记录](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?enzyme=Aurora+kinase+A&fil=tg&monomerid=50243388&pcid=p1325&submit=summary&tag=rep) · [Klaeger 等](https://doi.org/10.1126/science.aan4368) · [Lake 等](https://doi.org/10.1073/pnas.1811158115)

仓库原有 benchmarks/kinase_state_affinity/source/sources.json 给出两份原始工作表的下载地址与 SHA-256。以 Python、pandas、xlrd 运行 `python davis.py WORKBOOK_DIR PATH_TO_sources.json davis.csv` 即可复算。[逐配体记录](davis.csv)
