"""Rescore the archived first AD4 pose of all eight FKBP12 compounds on both maps."""
import argparse
import csv
import itertools
import statistics
import tarfile
from pathlib import Path

from scipy.stats import kendalltau, spearmanr

from boundary import BENCH, call, pose_blocks


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("output", "parent", "fixed"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    archive = BENCH / "cubic_grid/fkbp_ad4"
    for pdb in ("1FKG", "1FKH"):
        with tarfile.open(archive / f"maps-{pdb}.tar.xz") as tar:
            tar.extractall(a.output, filter="data")
    rows = []
    with tarfile.open(archive / "pools.tar.xz") as tar:
        for name in sorted(n for n in tar.getnames() if n.endswith("_ad4.pdbqt")):
            ligand, pdb, seed, arm = Path(name).stem.split("_")
            pose = a.output / "fixed.pdbqt"
            pose.write_text(pose_blocks(tar.extractfile(name).read().decode())[0])
            args = ["--scoring", "ad4", "--ligand", pose, "--maps", a.output/"maps"/pdb/"receptor", "--score_only"]
            for label in ("parent", "fixed"):
                rc, score = call(getattr(a, label), args)
                assert rc == 0
                rows.append(dict(pdb=pdb, compound=int(ligand), seed=seed, arm=label, score=score))
    write(a.output / "fixed_affinity.csv", rows)
    labels = {int(r["compound"]): r for r in csv.DictReader((BENCH/"cubic_grid/human_panel/affinity/ligands.csv").open())}
    unique = [n for n in labels if labels[n]["ki_low_nm"] == labels[n]["ki_high_nm"]]
    eligible = [(i,j) for i,j in itertools.permutations(labels, 2)
                if float(labels[j]["ki_low_nm"]) >= 3 * float(labels[i]["ki_high_nm"])]
    stats = []
    for arm in ("parent", "fixed"):
        medians = {(pdb,n): statistics.median(r["score"] for r in rows
                   if (r["pdb"],r["compound"],r["arm"]) == (pdb,n,arm))
                   for pdb in ("1FKG","1FKH") for n in labels}
        for pdb in ("1FKG","1FKH","receptor_mean"):
            score = {n: medians[pdb,n] if pdb != "receptor_mean" else (medians["1FKG",n]+medians["1FKH",n])/2
                     for n in labels}
            experimentally = [float(labels[n]["ki_low_nm"]) for n in unique]
            predicted = [score[n] for n in unique]
            stats.append(dict(pdb=pdb,arm=arm,eligible_pairs=len(eligible),
                              correct=sum((score[i] < score[j]) + 0.5*(score[i] == score[j]) for i,j in eligible),
                              points=len(unique),spearman=spearmanr(experimentally,predicted).statistic,
                              kendall=kendalltau(experimentally,predicted).statistic))
    write(a.output / "fixed_affinity_stats.csv", stats)


if __name__ == "__main__":
    main()
