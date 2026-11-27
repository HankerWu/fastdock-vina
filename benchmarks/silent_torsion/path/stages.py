"""Trace, or re-evaluate, four selected paired cases. The search has no reference coordinates."""

import argparse
import hashlib
import json
import lzma
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "transfer"))
from panel import holdout

# Retained success, its low-budget regression, and two independent chemical environments.
CELLS = (
    ("1FBY", "1FBY", 51871, 12000),
    ("1FBY", "1FBY", 51871, 3000),
    ("2LBD", "2LBD", 67213, 12000),
    ("1T49", "1T48", 202711, 12000),
)


def summarize(trace, query, recv, key):
    rot, trans, *_ = holdout.paired.base.alignment(query, recv)
    reference = np.asarray(query["mapped"]) @ rot + trans
    rows = []
    for line in trace:
        fields = line.strip().split(",")
        if fields[0] == "serial":
            serial = [int(v) for v in fields[1:]]
            order = [serial.index(i) for i in query["serial"]]
            if len(serial) != len(order):
                raise ValueError(key)
            continue
        step, stage, evals, energy, kind = fields[:5]
        # At proposal generation candidate.e is a reused buffer and has no physical meaning.
        coords = np.array(fields[5:], dtype=np.float64).reshape(-1, 3)[order]
        rmsd = np.sqrt(np.square(reference - coords).sum(axis=2).mean(axis=1)).min()
        rows.append(
            dict(
                key=key,
                step=int(step),
                stage=stage,
                evals=int(evals),
                energy=""
                if stage in ("proposal", "native_proposal")
                or (stage == "source" and step == "0")
                else float(energy),
                detail=kind,
                rmsd_A=rmsd,
                coordinates_sha256=hashlib.sha256(
                    coords.astype("<f8").tobytes()
                ).hexdigest(),
            )
        )
    return rows


def main(a):
    out = a.output.resolve()
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    all_events = []
    summary = []
    rows_by_key = {}
    for q, r, seed, cap in CELLS:
        data = (a.dev if q == "1T49" else a.transfer).resolve()
        m = {p: json.loads((data / p / "preparation.json").read_text()) for p in {q, r}}
        for arm in ("parent", "always", "adaptive"):
            key = f"{q}_{r}_{seed}_{cap}_{arm}"
            prefix = raw / key
            compressed = prefix.with_suffix(".trace.xz")
            log = prefix.with_suffix(".log")
            final = prefix.with_suffix(".pdbqt")
            plain = prefix.with_suffix(".trace")
            if a.action == "observe":
                if plain.exists():
                    plain.unlink()
                args = holdout.args_for(data, m, q, r, seed, cap, final) + [
                    "--search_counters"
                ]
                with log.open("w") as l:
                    subprocess.run(
                        [str((a.binaries / arm).resolve()), *args],
                        stdout=l,
                        stderr=subprocess.STDOUT,
                        env={**os.environ, "FASTDOCK_PATH_TRACE": str(plain)},
                        check=True,
                    )
                with (
                    plain.open("rb") as f,
                    lzma.open(compressed, "wb", preset=5) as dest,
                ):
                    shutil.copyfileobj(f, dest)
                plain.unlink()
            official = data / "runs" / final.name
            if final.read_bytes() != official.read_bytes():
                raise ValueError("public pose bytes changed: " + key)
            if holdout.paired.counts(log.read_text()) != holdout.paired.counts(
                (data / "runs" / log.name).read_text()
            ):
                raise ValueError("evaluation work changed: " + key)
            with lzma.open(compressed, "rt") as f:
                rows = summarize(f, m[q], m[r], key)
            all_events.extend(rows)
            rows_by_key[key] = rows
            by_stage = {
                x: [s for s in rows if s["stage"] == x]
                for x in (
                    "proposal",
                    "hunt",
                    "accepted",
                    "second",
                    "pool",
                    "explicit_in",
                    "explicit_out",
                )
            }
            cell = dict(
                key=key, output_sha256=hashlib.sha256(final.read_bytes()).hexdigest()
            )
            for stage, records in by_stage.items():
                hits = [rec for rec in records if rec["rmsd_A"] <= 2]
                cell[stage + "_best_A"] = min(
                    (rec["rmsd_A"] for rec in records), default=""
                )
                cell[stage + "_first"] = hits[0]["step"] if hits else ""
                cell[stage + "_hits"] = len({rec["step"] for rec in hits})
            official_stats = json.loads((data / "runs" / (key + ".json")).read_text())
            modes = holdout.paired.inspect(data / "runs", m[q], m[r], official_stats)
            cell["top1_A"] = modes[0]["rmsd_A"]
            cell["oracle_A"] = min(x["rmsd_A"] for x in modes)
            summary.append(cell)
    # The low-budget three-arm counterfactual is independently rerun with its own full output.
    # Confirm its MC records are the exact corresponding long run prefix.
    for arm in ("parent", "always", "adaptive"):

        def mc(cap):
            key = f"1FBY_1FBY_51871_{cap}_{arm}"
            return [
                {k: v for k, v in x.items() if k != "key"}
                for x in rows_by_key[key]
                if not x["stage"].startswith("explicit")
            ]

        short, long = mc(3000), mc(12000)
        if short != long[: len(short)]:
            raise ValueError("budget prefix mismatch: " + arm)
    holdout.paired.write(out / "events.csv", all_events)
    holdout.paired.write(out / "summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("observe", "analyze"))
    for opt in ("dev", "transfer", "output", "binaries"):
        p.add_argument(opt, type=Path)
    main(p.parse_args())
