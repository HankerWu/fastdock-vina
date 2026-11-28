"""Run each executable independently on the frozen reciprocal and cognate panels."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent))
import replay

ARMS = ("official", "parent", "stable", "candidate")
PRO_COUNTERS = ("stable", "candidate", "e7031e9", "2326430")
LIMIT = 2 * 1024**3


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csvfile(path, rows):
    with path.open("w", newline="") as stream:
        columns = rows[0] if rows else ("key", "rank", "rmsd", "score")
        if not rows and path.name != "ranks.csv":
            raise ValueError(f"unexpected empty table: {path}")
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def panel(fixtures):
    manifest = json.loads((HERE / "inputs.json").read_text())
    limits = HERE.parent / "reciprocal/per_complex/frozen_cap.json"
    if digest(limits) != manifest["cost_cap_source_sha256"]:
        raise ValueError("changed timing allocation")
    caps = {
        (r["query"], r["receiver"]): r["cap"] for r in json.loads(limits.read_text())
    }
    for protocol in ("reciprocal", "cognate"):
        artifact = manifest if protocol == "reciprocal" else manifest["cognate"]
        archive = ROOT / artifact["archive"]
        if digest(archive) != artifact["sha256"]:
            raise ValueError(f"changed archive {archive}")
        inputs = fixtures / protocol
        inputs.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive) as tar:
            tar.extractall(inputs, filter="data")
        for pdb, fields in artifact["inputs"].items():
            for name, hash_ in fields.items():
                if digest(inputs / pdb / name) != hash_:
                    raise ValueError(f"changed input {pdb}/{name}")
        if (
            protocol == "cognate"
            and digest(inputs / "manifest.json") != artifact["manifest_sha256"]
        ):
            raise ValueError("changed cognate manifest")
        for q, r, ref in replay.experiments(protocol, inputs):
            for seed in replay.SEEDS:
                yield (protocol, inputs, q, r, ref, seed, 12000)
            if protocol == "reciprocal":
                for seed in (93529, 93563):
                    yield ("cost", inputs, q, r, ref, seed, caps[q["pdb"], r["pdb"]])
                if q["pdb"] < r["pdb"]:
                    yield ("default", inputs, q, r, ref, replay.SEEDS[0], None)


def bounded():
    resource.setrlimit(resource.RLIMIT_AS, (LIMIT, LIMIT))


def execute(cmd, log, timeout_s):
    start = time.monotonic()
    with log.open("w") as output:
        process = subprocess.Popen(
            cmd, stdout=output, stderr=subprocess.STDOUT, preexec_fn=bounded
        )
        timeout = False
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                break
            if time.monotonic() - start >= timeout_s:
                process.kill()
                _, status, usage = os.wait4(process.pid, 0)
                timeout = True
                break
            time.sleep(0.02)
        process.returncode = os.waitstatus_to_exitcode(status)
    return dict(
        status="timeout" if timeout else "failed" if process.returncode else "ok",
        exit=process.returncode,
        cpu_s=round(usage.ru_utime + usage.ru_stime, 6),
        user_s=round(usage.ru_utime, 6),
        system_s=round(usage.ru_stime, 6),
        wall_s=round(time.monotonic() - start, 6),
        peak_rss_kib=usage.ru_maxrss,
    )


def dock(job, timeout_s=180):
    arm, binary, out, (cohort, inputs, query, receiver, ref, seed, cap) = job
    q, r = query["pdb"], receiver["pdb"]
    key = f"{cohort}-{q}-{r}-{seed}-{arm}"
    pose, log = out / (key + ".pdbqt"), out / (key + ".log")
    cmd = [
        str(binary),
        "--ligand",
        str(inputs / q / "ligand.pdbqt"),
        "--receptor",
        str(inputs / r / "receptor.pdbqt"),
        "--out",
        str(pose),
        "--cpu",
        "1",
        "--seed",
        str(seed),
    ]
    for axis, value in zip("xyz", receiver["center"]):
        cmd.extend(["--size_" + axis, "24", "--center_" + axis, str(value)])
    if cap is not None:
        if arm == "official":
            cap = 12000
        cmd.extend(
            [
                "--exhaustiveness",
                "3",
                "--max_evals",
                str(cap),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--spacing",
                "0.375",
            ]
        )
    if arm in PRO_COUNTERS:
        cmd.append("--search_counters")
    result = execute(cmd, log, timeout_s)
    rmsds, scores = [], []
    published = pose.read_text() if pose.is_file() else ""
    complete = published.startswith("MODEL 1\n") and published.endswith("ENDMDL\n")
    if result["status"] == "ok" and complete:
        rmsds = replay.reciprocal.measure(pose, query["serial"], ref)
        scores = [
            float(l.split()[3])
            for l in published.splitlines()
            if l.startswith("REMARK VINA RESULT:")
        ]
        if len(rmsds) != len(scores):
            raise ValueError(f"pose/score mismatch: {key}")
    elif result["status"] == "ok":
        result["status"] = "malformed" if published.strip() else "empty"
    text = log.read_text()
    hunt = re.findall(
        r"Search trajectory: .*?hunt_evals=(\d+).*?second_evals=(\d+)", text
    )
    refine = re.search(
        r"Search counters: .*?explicit_inputs=(\d+) explicit_calls=(\d+) explicit_evals=(\d+) late_prescore=(\d+)",
        text,
    )
    if (
        rmsds
        and arm in PRO_COUNTERS
        and (len(hunt) != (8 if cap is None else 3) or not refine)
    ):
        raise ValueError(f"incomplete native counts: {key}")
    row = dict(
        cohort=cohort,
        target=query["target"],
        query=q,
        receiver=r,
        seed=seed,
        arm=arm,
        cap=cap,
        **result,
        modes=len(rmsds),
        first=rmsds[0] if rmsds else None,
        best5=min(rmsds[:5], default=None),
        best=min(rmsds, default=None),
        regret=rmsds[0] - min(rmsds) if rmsds else None,
        native_score=scores[0] if scores else None,
        hunt_evals=sum(int(h) for h, s in hunt) if hunt else None,
        second_evals=sum(int(s) for h, s in hunt) if hunt else None,
        explicit_evals=int(refine[3]) if refine else None,
        explicit_calls=int(refine[2]) if refine else None,
        late_prescore=int(refine[4]) if refine else None,
        pool_sha256=digest(pose) if pose.exists() else None,
        key=key,
    )
    return row, [
        dict(key=key, rank=i, rmsd=a, score=s)
        for i, (a, s) in enumerate(zip(rmsds, scores), 1)
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("arm", choices=ARMS)
    ap.add_argument("binary", type=Path)
    ap.add_argument("work", type=Path, help="fresh directory for this arm")
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()
    binary, out = args.binary.resolve(), args.work.resolve()
    (out / "raw").mkdir(parents=True, exist_ok=False)
    jobs = [(args.arm, binary, out / "raw", j) for j in panel(out / "fixtures")]
    if len(jobs) != 65:
        raise ValueError("frozen cell count mismatch")
    # Start the slower algorithm-default cohort first.
    jobs.sort(key=lambda j: j[-1][0] != "default")
    with ProcessPoolExecutor(args.jobs) as pool:
        results = list(pool.map(dock, jobs))
    csvfile(out / "cells.csv", [r for r, _ in results])
    csvfile(out / "ranks.csv", [rank for _, a in results for rank in a])
    (out / "binary.sha256").write_text(digest(binary) + "\n")


if __name__ == "__main__":
    main()
