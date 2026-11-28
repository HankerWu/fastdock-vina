"""Consolidate replay results, retain every diagnostic log, and audit the winners."""

import argparse
import csv
import json
from pathlib import Path
import statistics
import tarfile


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


KEY = ("arm", "cohort", "query", "receiver")


def audit(pairs, scores):
    from collections import defaultdict

    by_pair = defaultdict(list)
    for row in scores:
        by_pair[tuple(row[k] for k in KEY)].append(row)
    checked = 0
    for p in pairs:
        modes = by_pair.pop(tuple(p[k] for k in KEY), [])
        if not modes:
            if p["status"] != "unavailable_search":
                raise ValueError("missing modes for completed search")
            continue
        if int(p["modes"]) != len(modes) or float(p["oracle"]) != min(
            float(m["rmsd"]) for m in modes
        ):
            raise ValueError("published pool size or oracle mismatch")
        for strategy in ("original", "archived", "batch", "fixed", "minus3"):
            if not p.get(strategy + "_winner"):
                if all(m.get(strategy) for m in modes):
                    raise ValueError("scorable pair missing from metrics")
                continue
            ordered = sorted(
                modes,
                key=lambda m: (float(m[strategy]), int(m["seed"]), int(m["rank"])),
            )
            best = ordered[0]
            if best["file"] != p[strategy + "_winner"] or float(best["rmsd"]) != float(
                p[strategy + "_A"]
            ):
                raise ValueError("raw score cannot reproduce selection")
            checked += 1
    if by_pair:
        raise ValueError("unmatched modes")
    return checked


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("work", type=Path, nargs="+")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    tables = {name: [] for name in ("pairs", "scores", "costs")}
    with tarfile.open(a.out / "logs.tar.xz", "w:xz") as tar:
        for folder in a.work:
            for key in tables:
                tables[key].extend(read(folder / (key + ".csv")))
            for log in sorted(folder.glob("*/*.log")):
                tar.add(log, arcname=log.relative_to(folder))
    if len({tuple(r[k] for k in KEY) for r in tables["pairs"]}) != len(tables["pairs"]):
        raise ValueError("duplicate pair")
    for name, rows in tables.items():
        write(a.out / (name + ".csv"), rows)
    report = []
    for key in sorted({(r["cohort"], r["arm"]) for r in tables["pairs"]}):
        group = [r for r in tables["pairs"] if (r["cohort"], r["arm"]) == key]
        valid = [r for r in group if r["status"] == "ok"]
        row = dict(cohort=key[0], arm=key[1], scheduled=len(group), complete=len(valid))
        if valid:
            for name in ("original", "archived", "batch", "fixed"):
                row["mean_" + name + "_A"] = statistics.mean(
                    float(r[name + "_A"]) for r in valid
                )
                row[name + "_change_vs_printed"] = sum(
                    r[name + "_winner"] != r["original_winner"] for r in valid
                )
            ds = [float(r["fixed_A"]) - float(r["original_A"]) for r in valid]
            row.update(
                mean_oracle_A=statistics.mean(float(r["oracle"]) for r in valid),
                improved=sum(v < -1e-8 for v in ds),
                worsened=sum(v > 1e-8 for v in ds),
                tied=sum(abs(v) <= 1e-8 for v in ds),
            )
        report.append(row)
    write(a.out / "summary.csv", report)
    receipt = {
        "pairs": len(tables["pairs"]),
        "coordinates": len(tables["scores"]),
        "verified_selections": audit(tables["pairs"], tables["scores"]),
        "probe_and_batch_binaries": json.loads(
            (a.work[0] / "binaries.json").read_text()
        ),
    }
    if any(
        json.loads((w / "binaries.json").read_text())
        != receipt["probe_and_batch_binaries"]
        for w in a.work
    ):
        raise ValueError("inconsistent rescoring executables")
    (a.out / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
