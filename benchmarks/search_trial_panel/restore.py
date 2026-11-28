"""Compare the restored production executable and a full 65-cell rerun to the recorded release."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

HERE = Path(__file__).resolve().parent
UP = HERE.parent / "geometry_comparison/upstream/four_versions"


def digest(b):
    return hashlib.sha256(b).hexdigest()


def binary_record(path):
    raw = path.read_bytes()
    build = re.search(
        r"Build ID: ([a-f0-9]+)",
        subprocess.check_output(["readelf", "-n", path], text=True),
    )
    prefix = re.search(rb"/[^\x00\n]+?/src/lib/ad4cache\.cpp\x00", raw)
    source = prefix.group().rsplit(b"/src/", 1)[0]
    elf_id = bytes.fromhex(build.group(1))
    normalized = raw.replace(source, b"<source-root>").replace(
        elf_id, b"<link-build-id>"
    )
    return dict(
        sha256=digest(raw),
        compiled_path=source.decode(),
        path_literals=raw.count(source),
        elf_build_id=elf_id.hex(),
        normalized_sha256=digest(normalized),
    )


def read(p):
    with p.open() as f:
        return list(csv.DictReader(f))


def key(r):
    return tuple(r[k] for k in ("cohort", "query", "receiver", "seed"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "restored", "runs", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    versions = {n: binary_record(getattr(a, n)) for n in ("baseline", "restored")}
    baseline = {
        key(r): r for r in read(UP / "results/cells.csv") if r["arm"] == "candidate"
    }
    runs = read(a.runs / "cells.csv")
    assert len(runs) == len(baseline) == 65 and len(set(map(key, runs))) == 65
    for r in runs:
        assert (
            digest((a.runs / "raw" / (r["key"] + ".pdbqt")).read_bytes())
            == r["pool_sha256"]
        )
        assert all(
            r[k] == baseline[key(r)][k]
            for k in (
                "pool_sha256",
                "modes",
                "status",
                "hunt_evals",
                "second_evals",
                "explicit_evals",
                "explicit_calls",
                "late_prescore",
            )
        )
    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    receipt = dict(
        binaries=versions,
        source_commit="d37e5e4743f2460c89726d18b6452487ae085096",
        normalized_binary_identical=versions["baseline"]["normalized_sha256"]
        == versions["restored"]["normalized_sha256"],
        byte_identical_pools=len(runs),
        modes=sum(int(r["modes"]) for r in runs),
        cpu_s=sum(float(r["cpu_s"]) for r in runs),
        summed_wall_s=sum(float(r["wall_s"]) for r in runs),
        peak_rss_kib=max(int(r["peak_rss_kib"]) for r in runs),
    )
    for name in ("cells.csv", "ranks.csv"):
        shutil.copyfile(a.runs / name, out / name)
    with tarfile.open(out / "raw.tar.xz", "w:xz") as t:
        t.add(a.runs / "raw", arcname="raw")
    (out / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
