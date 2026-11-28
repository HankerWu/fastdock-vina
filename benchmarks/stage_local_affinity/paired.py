"""Compare frozen >=3-fold assay pairs across source versions and receptor structures."""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
from pathlib import Path


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def paired(out):
    # Every nominal (including missing) expected pair is emitted by evaluate.py.
    pairs = [
        row for row in read(out / "pairs.csv")
        if row["stage"] == "nominal" and Decimal(row["fold"]) >= 3
    ]
    groups = defaultdict(dict)
    for row in pairs:
        identity = tuple(row[k] for k in ("target", "assay", "receptor", "stronger", "weaker"))
        assert row["engine"] not in groups[identity]
        groups[identity][row["engine"]] = row

    direct, changes = [], defaultdict(list)
    for (target, assay, receptor, stronger, weaker), engines in sorted(groups.items()):
        assert set(engines) == {"official", "parent", "stable", "candidate"}
        candidate = engines["candidate"]
        for baseline in ("official", "parent", "stable"):
            other = engines[baseline]
            complete = candidate["credit"] != "" and other["credit"] != ""
            delta = Decimal(candidate["credit"]) - Decimal(other["credit"]) if complete else None
            changes[target, assay, receptor, baseline].append(delta)
            direct.append(dict(target=target, assay=assay, receptor=receptor, stronger=stronger,
                               weaker=weaker, fold=candidate["fold"], same_scaffold=candidate["same_scaffold"],
                               baseline=baseline, baseline_gap=other["weaker_minus_stronger"],
                               candidate_gap=candidate["weaker_minus_stronger"],
                               baseline_credit=other["credit"], candidate_credit=candidate["credit"],
                               delta_credit=delta))
    write(out / "pair_changes.csv", direct)
    write(out / "pair_summary.csv", [
        dict(target=t, assay=a, receptor=r, baseline=b, planned=len(vals),
             observed=sum(v is not None for v in vals),
             improved=sum(v > 0 for v in vals if v is not None),
             worsened=sum(v < 0 for v in vals if v is not None),
             tied=sum(v == 0 for v in vals if v is not None),
             credit_change=sum(v for v in vals if v is not None))
        for (t, a, r, b), vals in changes.items()
    ])

    # Compare the SAME measured molecular pair (one assay) in each prepared BACE1 receptor.
    receptor_changes = []
    for (target, assay, receptor, stronger, weaker), engines in sorted(groups.items()):
        if target != "BACE1" or receptor != "4DJW":
            continue
        second = groups[target, assay, "4GID", stronger, weaker]
        for engine, first in engines.items():
            other = second[engine]
            receptor_changes.append(dict(assay=assay, engine=engine, stronger=stronger,
                weaker=weaker, fold=first["fold"], gap_4DJW=first["weaker_minus_stronger"],
                gap_4GID=other["weaker_minus_stronger"], credit_4DJW=first["credit"],
                credit_4GID=other["credit"], same_scaffold=first["same_scaffold"]))
    write(out / "bace1_receptor_pairs.csv", receptor_changes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    paired(parser.parse_args().out)
