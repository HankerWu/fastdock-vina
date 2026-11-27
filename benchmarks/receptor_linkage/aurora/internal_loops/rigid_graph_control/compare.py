"""Read production parser graphs and compare individual arms with the primary output pools."""

import argparse
import csv
import io
import json
import subprocess
import tarfile
from pathlib import Path
from graph_control import scan


def csv_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def compare(a, b, out, probe):
    systems = scan.config(a, b, out)
    primary = Path(__file__).resolve().parent.parent / "rigid_coupled"
    all_graphs = []
    summaries = {}
    for drug, c in systems.items():
        graphs = []
        for label in ("inferred", "explicit"):
            args = [probe, c["ligand"]]
            if label == "explicit":
                args.append(c["bonds"])
            rows = list(
                csv.DictReader(io.StringIO(subprocess.check_output(args, text=True)))
            )
            all_graphs.extend(dict(drug=drug, graph=label, **r) for r in rows)
            bonds = {
                (
                    min(int(r["atom"]), int(r["neighbor"])),
                    max(int(r["atom"]), int(r["neighbor"])),
                )
                for r in rows
            }
            assert 2 * len(bonds) == len(rows)
            graphs.append(bonds)
        inf, exp = graphs
        expected = {
            tuple(sorted(map(int, line.split())))
            for line in c["bonds"].read_text().splitlines()
        }
        assert exp == expected
        summaries[drug] = dict(
            inferred=len(inf),
            explicit=len(exp),
            missing=sorted(exp - inf),
            extra=sorted(inf - exp),
        )
    (out / "graphs.json").write_text(json.dumps(summaries, indent=2) + "\n")
    scan.replay.geom.write_csv(out, "adjacency.csv", all_graphs)
    main_cells = {
        (r["drug"], r["seed"], r["arm"]): r
        for r in csv_rows(primary / "results/cells.csv")
    }
    pairs = []
    with tarfile.open(primary / "raw.tar.xz") as t:
        for new in csv_rows(out / "cells.csv"):
            drug, seed, arm = (new[k] for k in ("drug", "seed", "arm"))
            old = main_cells[drug, seed, arm]
            official = main_cells[drug, seed, "official"]
            fname = f"{drug}_{seed}_{arm}.pdbqt"
            member = next(m for m in t.getnames() if m.endswith("raw/" + fname))
            equal = (out / "raw" / fname).read_bytes() == t.extractfile(member).read()
            pairs.append(
                dict(
                    drug=drug,
                    seed=seed,
                    arm=arm,
                    byte_equal=int(equal),
                    explicit_top1=old["top1"],
                    inferred_top1=new["top1"],
                    official_inferred_top1=official["top1"],
                    explicit_oracle=old["oracle"],
                    inferred_oracle=new["oracle"],
                    official_inferred_oracle=official["oracle"],
                )
            )
    scan.replay.geom.write_csv(out, "paired.csv", pairs)
    print(json.dumps(summaries), sum(r["byte_equal"] for r in pairs))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("aurora", "brd9", "out", "probe"):
        p.add_argument(k, type=lambda x: Path(x).resolve())
    a = p.parse_args()
    compare(a.aurora, a.brd9, a.out, a.probe)
