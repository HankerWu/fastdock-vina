"""Independently reconstruct first-pose ranks and matched pairs using integer units."""

import argparse
import csv
import hashlib
import itertools
import json
import re
from decimal import Decimal
from pathlib import Path

import numpy as np
from panel import FIG
from scipy.stats import kendalltau, rankdata


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit(work):
    manifest = {int(r["compound"]): r for r in rows(work / "manifest.csv")}
    kd = {i: int(Decimal(str(FIG[i - 4])) * 100) for i in manifest}
    assert all(
        Decimal(r["kd_uM"]) * 100 == n
        for (i, r), n in zip(manifest.items(), kd.values())
    )
    sources = {
        "parent": ("fixed.csv", "parent"),
        "linear": ("fixed.csv", "linear"),
        "transported": ("transported.csv", "linear"),
    }
    if (work / "matched.csv").is_file():
        sources["matched"] = ("matched.csv", "linear")
    expected = {}
    for r in rows(work / "summary.csv"):
        expected[r["method"], r["receiver"], r["subset"]] = r
    methods = {}
    keys = set()
    coordinates = hashlib.sha256()
    unique_modes = 0
    counters = 0
    for method, (file, arm) in sources.items():
        pred = {}
        for row in rows(work / file):
            if row["arm"] != arm:
                continue
            i = int(row["compound"])
            rec = row["receptor"]
            key = row["key"]
            path = work / "docking/raw" / (key + ".pdbqt")
            blocks = path.read_text().split("ENDMDL")[:-1]
            assert len(blocks) == int(row["modes"])
            ticks = int(
                Decimal(
                    next(
                        s.split()[3]
                        for s in blocks[0].splitlines()
                        if s.startswith("REMARK VINA RESULT:")
                    )
                )
                * 1000
            )
            assert Decimal(row["first_score"]) * 1000 == ticks
            pred[rec, i] = ticks
            if key in keys:
                continue
            keys.add(key)
            needed = {
                s[6:11]: s.split()[-1]
                for s in (work / "prepared" / (str(i) + ".pdbqt"))
                .read_text()
                .splitlines()
                if s.startswith(("ATOM  ", "HETATM"))
            }
            for block in blocks:
                ls = [
                    s for s in block.splitlines() if s.startswith(("ATOM  ", "HETATM"))
                ]
                assert {s[6:11]: s.split()[-1] for s in ls} == needed
                xyz = [[float(s[j : j + 8]) for j in (30, 38, 46)] for s in ls]
                assert np.isfinite(xyz).all()
                coordinates.update(
                    "".join(s[6:11] + s[30:54] for s in ls).encode("ascii")
                )
                unique_modes += 1
            log = (path.with_suffix(".log")).read_text()
            for name in (
                "hunt_calls",
                "hunt_evals",
                "second_calls",
                "second_evals",
                "explicit_calls",
                "explicit_evals",
            ):
                matches = re.findall(r"\b" + name + r"=(\d+)\b", log)
                assert len(matches) == 1 and int(matches[0]) == int(row[name])
                counters += 1
        for i in manifest:
            pred["consensus", i] = pred["1KJL", i] + pred["3T1L", i]
        methods[method] = pred
    methods["size"] = {
        (rec, i): -int(r["heavy"])
        for i, r in manifest.items()
        for rec in ("1KJL", "3T1L", "consensus")
    }
    for (m, rec, subset), row in expected.items():
        ids = sorted(
            i for i in manifest if subset != "without_bivalent" or i not in (45, 46)
        )
        pred = methods[m]
        credits = []
        for i, j in itertools.combinations(ids, 2):
            if max(kd[i], kd[j]) < 3 * min(kd[i], kd[j]):
                continue
            if subset == "equal_heavy_charge" and (
                manifest[i]["heavy"],
                manifest[i]["charge"],
            ) != (manifest[j]["heavy"], manifest[j]["charge"]):
                continue
            sign = (kd[i] - kd[j]) * (pred[rec, i] - pred[rec, j])
            credits.append(2 if sign > 0 else 0 if sign < 0 else 1)
        assert 2 * Decimal(row["credit"]) == sum(credits) and len(credits) == int(
            row["n_pairs"]
        )
        if row["spearman"]:
            x, y = np.array([kd[i] for i in ids]), np.array([pred[rec, i] for i in ids])
            rho = np.corrcoef(rankdata(x), rankdata(y))[0, 1]
            assert (
                abs(rho - float(row["spearman"])) < 1e-12
                and abs(kendalltau(x, y).statistic - float(row["kendall_tau_b"]))
                < 1e-12
            )
    return {
        "unique_runs": len(keys),
        "unique_modes_checked": unique_modes,
        "counter_fields": counters,
        "summary_rows_reproduced": len(expected),
        "all_pose_coordinates_sha256": coordinates.hexdigest(),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    work = p.parse_args().work
    result = audit(work)
    (work / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result)
