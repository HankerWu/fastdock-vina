"""Archive the initial candidate or supplement predeclared official default timeouts."""

import argparse
import json
from pathlib import Path
import tarfile

from evaluate import read, verify
from run_panel import csvfile, digest, dock, panel


def archive(work, out, expected):
    out.mkdir(exist_ok=True, parents=True)
    with tarfile.open(out / "raw.tar.xz", "w:xz") as tar:
        for file in sorted((work / "raw").iterdir()):
            tar.add(file, arcname=file.name)
    for name in ("cells.csv", "ranks.csv"):
        (out / name).write_bytes((work / name).read_bytes())
    (out / "binary.sha256").write_bytes((work / "binary.sha256").read_bytes())
    audit = verify(out, expected)
    (out / "verification.json").write_text(json.dumps(audit, indent=2) + "\n")


def retry(query, binary, work):
    (work / "raw").mkdir(parents=True, exist_ok=False)
    cells = [
        c
        for c in panel(work / "fixtures")
        if c[0] == "default" and c[2]["pdb"] == query
    ]
    if len(cells) != 1:
        raise ValueError("query not in frozen default subset")
    row, ranks = dock(
        ("official", binary.resolve(), work / "raw", cells[0]), timeout_s=900
    )
    csvfile(work / "cells.csv", [row])
    csvfile(work / "ranks.csv", ranks)
    (work / "binary.sha256").write_text(digest(binary) + "\n")


def combine(primary, out, supplement):
    stock = [
        r
        for r in read(primary / "cells.csv")
        if r["cohort"] == "default" and r["arm"] == "official"
    ]
    replacement = {r["key"]: r for r in read(supplement / "cells.csv")}
    if set(replacement) != {r["key"] for r in stock if r["status"] == "timeout"}:
        raise ValueError(
            "supplement must correspond exactly to primary official default timeouts"
        )
    rows = []
    for r in stock:
        use = replacement.get(r["key"], r)
        rows.append(dict(**use, primary_status=r["status"], used_900s=use is not r))
    csvfile(out / "completed_defaults.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    stage = p.add_subparsers(dest="stage", required=True)
    t = stage.add_parser("retry")
    t.add_argument("query")
    t.add_argument("binary", type=Path)
    t.add_argument("work", type=Path)
    a = stage.add_parser("archive")
    a.add_argument("work", type=Path)
    a.add_argument("out", type=Path)
    a.add_argument("expected", type=int)
    c = stage.add_parser("combine")
    c.add_argument("primary", type=Path)
    c.add_argument("out", type=Path)
    c.add_argument("supplement", type=Path)
    args = vars(p.parse_args())
    cmd = args.pop("stage")
    {"retry": retry, "archive": archive, "combine": combine}[cmd](**args)
