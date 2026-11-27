"""Independently rank identical hunt geometries with the smooth explicit receptor and shared ligand terms."""

import argparse
import hashlib
import io
import os
import tarfile
from pathlib import Path

from measure import COVERAGE, HERE, counters, load, run_binary, sha, table


def results(src, out, raw):
    native = {r["key"]: r for r in load(src / "cases.csv")}
    endpoints = load(src / "endpoints.csv")
    rows, summary = [], []
    with tarfile.open(raw) as tar:
        for key, case in native.items():
            log = tar.extractfile(key + ".stderr").read()
            plain = tar.extractfile(key + ".stdout").read()
            chk = tar.extractfile(key + ".reference").read()
            digest = hashlib.sha256(tar.extractfile(key + ".pdbqt").read()).hexdigest()
            if (
                digest != case["gated_sha256"]
                or counters(plain) != counters(chk)
                or len(counters(chk)) != 4
            ):
                raise ValueError("Atomic observer altered native work or pose file")
            parsed, soft = {}, {}
            for line in log.splitlines():
                if line.startswith(b"FDATOM "):
                    _, i, g, a, inside, c = line.split()
                    parsed[int(i)] = (float(g), float(a), int(inside), float(c))
                if line.startswith(b"FDCAPS "):
                    _, i, s, h, d = line.split()
                    if float(d) != 0:
                        raise ValueError("Private geometry changed")
                    soft[int(i)] = (float(s), float(h))
            domain = [r for r in endpoints if r["key"] == key]
            if len(parsed) != len(domain) or set(parsed) != set(soft):
                raise ValueError("Incomplete rescoring")
            pool = []
            for r in domain:
                i = int(r["step"])
                g, a, inside, c = parsed[i]
                s, h = soft[i]
                if s != float(r["soft"]) or h != float(r["authentic"]):
                    raise ValueError(
                        "Grid energies differ from independent prior replay"
                    )
                intra = h - g
                total = a + intra
                if i == min(parsed) and abs(total - c) > 1e-8:
                    raise ValueError(
                        "Complete atomic model objective does not match reconstructed energy"
                    )
                item = {
                    "key": key,
                    "step": i,
                    "rmsd_A": float(r["rmsd_A"]),
                    "novel": min(
                        float(r["nearest_saved_A"]), float(r["nearest_merged_A"])
                    )
                    >= 2,
                    "oracle": r["geometry_oracle"] == "True",
                    "grid_joint": h,
                    "grid_inter": g,
                    "atomic_inter": a,
                    "shared_ligand": intra,
                    "atomic_joint": total,
                    "inside": inside,
                    "atomic_joint_control": c if i == min(parsed) else "",
                }
                rows.append(item)
                if item["novel"]:
                    pool.append(item)
            oracle = next(r for r in pool if r["oracle"])
            best = {
                "key": key,
                "skipped": len(domain),
                "outside": sum(1 - x[2] for x in parsed.values()),
                "eligible": len(pool),
                "grid_joint_calls": len(domain),
                "extra_grid_receptor_calls": len(domain),
                "atomic_receptor_calls": len(domain),
                "extra_atomic_joint_controls": 1,
            }
            for metric in (
                "grid_joint",
                "atomic_joint",
                "grid_inter",
                "atomic_inter",
                "shared_ligand",
            ):
                order = sorted(pool, key=lambda r: (r[metric], r["step"]))
                winner = order[0]
                closest = next(r for r in order if r["rmsd_A"] <= 2)
                best.update(
                    {
                        metric + "_step": winner["step"],
                        metric + "_rmsd": winner["rmsd_A"],
                        metric + "_oracle_rank": order.index(oracle) + 1,
                        metric + "_best_hit_rank": order.index(closest) + 1,
                        metric + "_best_hit_gap": closest[metric] - winner[metric],
                    }
                )
            wrong = next(r for r in pool if r["step"] == best["atomic_joint_step"])
            for metric in (
                "grid_inter",
                "atomic_inter",
                "shared_ligand",
                "atomic_joint",
            ):
                best["oracle_minus_atomic_winner_" + metric] = (
                    oracle[metric] - wrong[metric]
                )
            summary.append(best)
    table(out / "atomic_endpoints.csv", rows)
    table(out / "atomic_cases.csv", summary)


def run(a):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0)) if a.cpu is None else a.cpu})
    a.work.mkdir(parents=True, exist_ok=True)
    with tarfile.open(
        HERE.parents[1] / "cubic_grid/human_panel/prepared.tar.xz"
    ) as tar:
        tar.extractall(a.work, filter="data")
    raw = a.results / "atomic_runs.tar.xz"
    cost = []
    with tarfile.open(raw, "w:xz") as archive:
        for r in load(COVERAGE / "results/opportunities.csv"):
            if not r["oracle_step"]:
                continue
            key = r["key"]
            pose = a.work / (key + ".pdbqt")
            ref, plain = run_binary(a.work, r, a.candidate, pose)
            p, obs = run_binary(a.work, r, a.observer, pose)
            for ext, data in (
                (".pdbqt", pose.read_bytes()),
                (".stdout", p.stdout),
                (".stderr", p.stderr),
                (".reference", ref.stdout),
            ):
                info = tarfile.TarInfo(key + ext)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
            cost.append(
                {
                    "key": key,
                    "candidate_s": plain,
                    "observer_s": obs,
                    "observer_binary_sha256": sha(a.observer),
                }
            )
    results(a.results, a.results, raw)
    table(a.results / "atomic_cost.csv", cost)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    p.add_argument("--work", type=Path)
    p.add_argument("--candidate", type=Path)
    p.add_argument("--observer", type=Path)
    p.add_argument("--cpu", type=int)
    a = p.parse_args()
    if a.work and a.observer and a.candidate:
        run(a)
    else:
        results(a.results, a.results, a.results / "atomic_runs.tar.xz")
