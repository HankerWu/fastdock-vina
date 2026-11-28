"""Verify every rerun of the production source against the historical nominal pools."""

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

COUNTS = ("hunt_evals", "second_evals", "explicit_evals", "late_prescore")
RESTORED_SHA256 = "181963d595b77eb7822e874a669882952765f685bef965f4090ef85951882aaa"
ARCHIVE_SHA256 = "f6ff9e0d9d6f0bc6e77e65409f9c7f88fccd826eca064a49e9c91a9a73496f9f"


def receipt(work, original, out):
    previous = original / "outputs/stable/nominal"
    fresh = work / "outputs/stable/nominal"
    new_files = sorted(fresh.glob("*.json"))
    assert len(new_files) == 140 and {p.name for p in new_files} == {p.name for p in previous.glob("*.json")}
    inputs = json.loads((work / "inputs.json").read_text())["files"]
    assert inputs == json.loads((original / "inputs.json").read_text())["files"]
    cpu_s = wall_s = 0.0
    peak = modes = 0
    totals = dict.fromkeys(COUNTS, 0)
    for path in new_files:
        now = json.loads(path.read_text())
        before = json.loads((previous / path.name).read_text())
        assert now["valid"] and before["valid"]
        assert now["binary_sha256"] == RESTORED_SHA256 and before["binary_sha256"] == ARCHIVE_SHA256
        for k in (*COUNTS, "seed", "assay", "compound", "cap", "top_score", "modes", "output_sha256"):
            assert now[k] == before[k], (path.name, k)
        body = path.with_suffix(".pdbqt").read_bytes()
        assert body == (previous / path.with_suffix(".pdbqt").name).read_bytes()
        assert hashlib.sha256(body).hexdigest() == now["output_sha256"]
        for k in COUNTS:
            totals[k] += now[k]
        cpu_s += now["user_s"] + now["sys_s"]
        wall_s += now["wall_s"]
        peak = max(peak, now["maxrss_kib"])
        modes += now["modes"]
    out.mkdir(parents=True, exist_ok=True)
    report = dict(binary_sha256=RESTORED_SHA256, reference_binary_sha256=ARCHIVE_SHA256,
                  identical_pools=len(new_files), total_modes=modes, counters=totals, cpu_s=cpu_s,
                  sum_job_wall_s=wall_s, max_rss_kib=peak)
    (out / "restored.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path, help="fresh nominal production runs")
    parser.add_argument("original", type=Path, help="study.py setup workspace containing archived stable controls")
    parser.add_argument("out", type=Path)
    parser.add_argument("--archive", action="store_true")
    args = parser.parse_args()
    receipt(args.work, args.original, args.out)
    if args.archive:
        with tarfile.open(args.out / "restored.tar.xz", "w:xz") as archive:
            for name in ("inputs.json", "outputs/stable/nominal"):
                archive.add(args.work / name, arcname=name)
