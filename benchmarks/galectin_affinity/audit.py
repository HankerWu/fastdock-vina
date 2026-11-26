"""Independently rescore every serialized top pose and reconstruct paired Kd outcomes."""

import argparse
import csv
import hashlib
import itertools
import json
import re
import statistics
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def csv_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(text, label):
    return float(
        next(l for l in text.splitlines() if l.startswith(label))
        .split(":")[1]
        .split()[0]
    )


def audit(work, binary, out, jobs):
    states = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L")}
    kd = csv_rows(Path(__file__).with_name("table1.csv"))
    prep = work / "prepared"
    reference = {}
    for lig in range(5, 18):
        values = []
        for raw in (0, 100):
            output = subprocess.check_output(
                [
                    binary,
                    "--score_only",
                    "--cpu",
                    "1",
                    "--autobox",
                    "--unbound_energy",
                    str(raw),
                    "--receptor",
                    prep / "1GZW.pdbqt",
                    "--ligand",
                    prep / f"{lig}.pdbqt",
                ],
                text=True,
            )
            values.append(number(output, "Estimated Free Energy of Binding"))
        reference[lig] = (
            number(output, "(2) Final Total Internal Energy"),
            (values[0] - values[1]) / 100,
        )

    def one(path):
        state, c, se, arm = path.stem.split("_", 3)
        top = path.read_text().split("ENDMDL", 1)[0].split("\n", 1)[1]
        saved = number(top, "REMARK VINA RESULT:")
        pool = number(top, "REMARK UNBOUND:")
        common, scale = reference[int(c)]
        # Independent stdout rank and fixed-coordinate evaluation of the serialized geometry.
        log = path.with_suffix(".log").read_text()
        first = float(re.search(r"^\s*1\s+(-?\d+\.\d+)\s+", log, re.MULTILINE).group(1))
        if abs(first - saved) > 0.00501:
            raise ValueError(path)
        with tempfile.TemporaryDirectory() as temp:
            pose = Path(temp) / "pose.pdbqt"
            pose.write_text(top)
            res = subprocess.check_output(
                [
                    binary,
                    "--score_only",
                    "--cpu",
                    "1",
                    "--config",
                    prep / f"{state}.box",
                    "--ligand",
                    pose,
                    "--receptor",
                    prep / f"{state}.pdbqt",
                    "--unbound_energy",
                    str(common),
                ],
                text=True,
            )
        common_e = number(res, "Estimated Free Energy of Binding")
        replay = common_e + scale * (common - pool)
        return {
            "state": state,
            "compound": c,
            "seed": se,
            "arm": arm,
            "saved": saved,
            "log_rank1": first,
            "common_reference": common,
            "pool_reference": pool,
            "scaling": scale,
            "replayed": replay,
            "common_score": common_e,
        }

    paths = sorted((work / "pools").glob("*.pdbqt"))
    with ThreadPoolExecutor(max_workers=jobs) as e:
        entries = list(e.map(one, paths))
    if len(entries) != 1040:
        raise ValueError(f"Expected 1040 cells, found {len(entries)}")
    out.mkdir(exist_ok=True, parents=True)
    write(out / "readback.csv", entries)
    energies = defaultdict(list)
    for row in entries:
        for kind in ("saved", "replayed", "common_score"):
            energies[row["state"], row["compound"], row["arm"], kind].append(row[kind])
    summary = []
    independent = []
    for target, recs in states.items():
        K = {row["compound"]: float(row[f"{target}_Kd_uM"]) for row in kd}
        for arm in ("native", "balanced", "anchored", "balanced_anchored", "prepared"):
            for subset in ("all", "triazoles"):
                pairs = [
                    (a, b)
                    for a, b in itertools.permutations(K, 2)
                    if K[b] / K[a] >= 3 and (subset == "all" or "17" not in (a, b))
                ]
                decisions = {}
                for kind in ("saved", "replayed", "common_score"):
                    E = {
                        k: sum(
                            statistics.median(energies[r, k, arm, kind]) for r in recs
                        )
                        / 2
                        for k in K
                    }
                    decisions[kind] = [
                        1 if E[a] + 1e-9 < E[b] else 0 if E[b] + 1e-9 < E[a] else 0.5
                        for a, b in pairs
                    ]
                summary.append(
                    {
                        "target": target,
                        "arm": arm,
                        "subset": subset,
                        "n": len(pairs),
                        "saved": sum(decisions["saved"]),
                        "replayed": sum(decisions["replayed"]),
                        "changed": sum(
                            i != j
                            for i, j in zip(decisions["saved"], decisions["replayed"])
                        ),
                        "common": sum(decisions["common_score"]),
                    }
                )
                if subset == "all":
                    for (a, b), s in zip(pairs, decisions["saved"]):
                        independent.append((target, arm, a, b, s))
    write(out / "readback_pairs.csv", summary)
    computed = [
        (r["target"], r["arm"], r["strong"], r["weak"], float(r["correct"]))
        for r in csv_rows(out / "pairs.csv")
        if r["arm"] != "larger_heavy"
    ]
    if sorted(computed) != sorted(independent):
        raise ValueError("Primary paired decisions disagree")
    metrics = {
        "top_poses": len(entries),
        "unique_target_compound_pairs": len(
            {(target, a, b) for target, _, a, b, _ in independent}
        ),
        "arm_pair_decisions_checked": len(independent),
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "max_abs_score_readback_difference": max(
            abs(r["saved"] - r["replayed"]) for r in entries
        ),
        "max_abs_log_saved_difference": max(
            abs(r["saved"] - r["log_rank1"]) for r in entries
        ),
    }
    (out / "audit.json").write_text(json.dumps(metrics, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    audit(args.work, args.binary.resolve(), args.out, args.jobs)
