"""Run and score the five FABP4 ligands from a single displacement assay."""
import argparse
import csv
import itertools
import math
from pathlib import Path
import statistics
import tarfile

from scipy.stats import kendalltau, spearmanr
from panel import HERE, GROUPS, SEEDS, dock_ad4, write_csv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("experiment", "flex", "old", "new"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    source = HERE.parent/"archive_overlap"/"transfer"
    labels = list(csv.DictReader((source/"fabp4"/"assay"/"measured.csv").open()))
    out = a.experiment.resolve()
    dest = out/"assay"
    dest.mkdir(exist_ok=True)
    prepared = {}
    for archive in (source/"inputs.tar.xz", source/"fabp4"/"assay"/"new_inputs.tar.xz"):
        with tarfile.open(archive) as t:
            for row in labels:
                name = row["pdb"]+"/ligand.pdbqt"
                if name in t.getnames():
                    path = dest/(row["id"]+".pdbqt")
                    path.write_bytes(t.extractfile(name).read())
                    prepared[row["id"]] = path
    cells = []
    for row,rx,seed in itertools.product(labels, GROUPS["FABP4"], SEEDS):
        ligand = row["id"]
        for arm,binary in (("old",a.old),("new",a.new)):
            path = dest/f"{ligand}-{rx}-{seed}-{arm}.pdbqt"
            cost = dock_ad4(binary, prepared[ligand], out/"maps"/rx/"receptor",
                            a.flex.resolve()/f"{rx}_flex.pdbqt", seed, path, True)
            if cost["status"]:
                raise RuntimeError(path.with_suffix(".log").read_text())
            modes = [float(line.split()[3]) for line in path.read_text().splitlines()
                     if line.startswith("REMARK VINA RESULT:")]
            cells.append(dict(id=ligand,receiver=rx,seed=seed,arm=arm,top1=modes[0] if modes else "",
                              reported=len(modes),**cost))
    write_csv(out/"assay_cells.csv",cells)
    predictions, pairs, metrics = [], [], []
    for rx,arm in itertools.product((*GROUPS["FABP4"],"mean"),("old","new")):
        scores = {}
        for label in labels:
            ligand=label["id"]
            vals=[r["top1"] for r in cells if r["arm"]==arm and r["id"]==ligand and
                  (rx=="mean" or rx==r["receiver"]) and r["top1"]!=""]
            value = statistics.mean(vals) if len(vals)==(6 if rx=="mean" else 3) else None
            predictions.append(dict(id=ligand,receiver=rx,arm=arm,top1_mean=value,
                                    pK=6-math.log10(float(label["Kd_uM"])),n=len(vals)))
            if value is not None:
                scores[ligand]=value
        concordance=[]
        for x,y in itertools.combinations(labels,2):
            kx,ky = [float(l["Kd_uM"]) for l in (x,y)]
            ex,ey = [float(l["reported_error_uM"]) for l in (x,y)]
            if abs(kx-ky)<=ex+ey or not {x["id"],y["id"]}<=scores.keys():
                continue
            delta=scores[y["id"]]-scores[x["id"]]
            agreement=.5 if delta==0 else int(delta*(ky-kx)>0)
            concordance.append(agreement)
            pairs.append(dict(receiver=rx,arm=arm,a=x["id"],b=y["id"],fold=max(kx/ky,ky/kx),
                              predicted_score_delta=delta,concordance=agreement))
        rows=[r for r in predictions if r["receiver"]==rx and r["arm"]==arm and r["id"] in scores]
        observed, predicted=[r["pK"] for r in rows],[-r["top1_mean"] for r in rows]
        metrics.append(dict(receiver=rx,arm=arm,ligands=len(rows),all_ligands=len(labels),
                            spearman=spearmanr(observed,predicted).statistic,kendall_b=kendalltau(observed,predicted).statistic,
                            correct=sum(concordance),eligible=len(concordance),total_pairs=math.comb(len(labels),2)))
    write_csv(out/"assay_predictions.csv",predictions)
    write_csv(out/"assay_pairs.csv",pairs)
    write_csv(out/"assay_metrics.csv",metrics)


if __name__ == "__main__":
    main()
