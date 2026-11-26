"""Same-assay rank statistics within each shared receptor and search seed."""

import csv
import itertools
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent


def load_table(path):
    with (ROOT / path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def study(trial, target, rows, receptors, measure):
    groups = defaultdict(list)
    for row in rows:
        receptor = row.get("receptor", row.get("rec"))
        if receptor in receptors:
            groups[(receptor, int(row["seed"]), row["arm"])].append(row)

    for (rec, seed, arm), cohort in sorted(groups.items()):
        ordered = sorted((float(r[measure]), float(r["score"])) for r in cohort)
        unique = {r.get("ligand", r.get("lig")) for r in cohort}
        if len(unique) != len(cohort):
            raise ValueError("Duplicate ligand in receptor/seed block")
        pairs = [
            (kb / ka, float(a < b) + 0.5 * (a == b))
            for (ka, a), (kb, b) in itertools.combinations(ordered, 2)
            if ka < kb
        ]
        wide = [credit for ratio, credit in pairs if ratio >= 3]
        rho = spearmanr(*zip(*ordered)).statistic
        yield dict(
            experiment=trial,
            target=target,
            receptor=rec,
            seed=seed,
            arm=arm,
            ligands=len(cohort),
            all_pairs=len(pairs),
            all_credit=sum(c for _, c in pairs),
            threefold_pairs=len(wide),
            threefold_credit=sum(wide),
            spearman=round(rho, 6),
        )


def main():
    ca = [
        r
        for r in load_table("hcaii/rows.csv")
        if r["box_A"] == "36" and r["panel"] == "main"
    ]
    panels = [
        ("smooth", "T4L", load_table("t4l/cells.csv"), ("1L90", "186L"), "kd_uM"),
        ("smooth", "hCAII", ca, ("6WQ4", "6WQ8"), "Ki_nM"),
        ("mixed", "T4L", load_table("mixed_t4l.csv"), ("1L90", "186L"), "measurement"),
        (
            "mixed",
            "hCAII",
            load_table("mixed_hcaii.csv"),
            ("6WQ4", "6WQ8"),
            "measurement",
        ),
        (
            "mixed",
            "BRD4",
            load_table("brd4/cells.csv"),
            ("3MXF", "4OGI", "4OGJ"),
            "kd_uM",
        ),
        (
            "mixed_34A",
            "BRD4",
            load_table("brd4/cells_34A.csv"),
            ("3MXF", "4OGI", "4OGJ"),
            "kd_uM",
        ),
    ]
    rows = [row for panel in panels for row in study(*panel)]
    writer = csv.DictWriter(sys.stdout, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
