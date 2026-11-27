"""Compute unpooled structure-specific results and the predeclared frozen-pool intervention."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(out):
    out.mkdir(exist_ok=True, parents=True)
    groups = defaultdict(list)
    intervention = []
    panels = [
        ("historical", HERE / "historical", None),
        ("historical", HERE / "historical/persistent", {"lock"}),
        ("external", HERE / "external/paired", None),
        ("external", HERE / "external/persistent", {"lock"}),
        ("transfer", HERE / "holdout/results", None),
        ("ITC", HERE / "fkbp_affinity/results", None),
    ]
    for panel, path, arms in panels:
        is_itc = panel == "ITC"
        rr = read(path / ("runs.csv" if is_itc else "structures.csv"))
        bykey = defaultdict(list)
        for m in read(path / ("modes.csv" if is_itc else "models.csv")):
            bykey[m["key"]].append(m)
        for r in rr:
            if int(r["cap"]) not in (3000, 12000) or (arms and r["arm"] not in arms):
                continue
            modes = sorted(bykey[r["key"]], key=lambda x: int(x["rank"]))
            if len(modes) != int(r["modes"]):
                raise ValueError(r["key"])
            rms = [float(m.get("rmsd_A", m.get("rmsd"))) for m in modes]
            bonds = [
                float(
                    m.get("closure_bond_A", m.get("cg_distance_A", m.get("closure_A")))
                )
                for m in modes
            ]
            ok = [i for i, b in enumerate(bonds) if b <= 3]
            index = ok[0] if ok else 0
            name = r["receiver"] + "_" + r["query"] if is_itc else r["pdb"]
            intervention.append(
                dict(
                    panel=panel,
                    structure=name,
                    key=r["key"],
                    arm=r["arm"],
                    modes=len(modes),
                    all_broken=int(not ok),
                    selected_rank=index + 1,
                    first_rmsd_A=rms[0],
                    selected_rmsd_A=rms[index],
                    first_bond_A=bonds[0],
                    selected_bond_A=bonds[index],
                )
            )
            groups[panel, name, r["arm"]].append((rms, bonds, int(r["joint"])))
    summary = []
    for (panel, target, arm), obs in sorted(groups.items()):
        top = [r[0] for r, b, n in obs]
        oracle = [min(r) for r, b, n in obs]
        summary.append(
            dict(
                panel=panel,
                structure=target,
                arm=arm,
                cells=len(obs),
                top1=sum(x <= 2 for x in top),
                top5=sum(min(r[:5]) <= 2 for r, b, n in obs),
                oracle=sum(x <= 2 for x in oracle),
                median_top1_A=round(median(top), 6),
                median_oracle_A=round(median(oracle), 6),
                top1_broken=sum(b[0] > 3 for r, b, n in obs),
                all_broken_modes=sum(sum(v > 3 for v in b) for r, b, n in obs),
                modes=sum(len(r) for r, b, n in obs),
                joint_evals=sum(n for r, b, n in obs),
            )
        )
    write(out / "target_summary.csv", summary)
    write(out / "valid_first.csv", intervention)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    main(parser.parse_args().out)
