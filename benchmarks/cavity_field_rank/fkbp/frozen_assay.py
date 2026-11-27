"""Replay the frozen FKBP12 trial-1, step-1 candidate bank and its full AD4 searches."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import time

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TRIPLET = ROOT / "benchmarks/cavity_triplet/fkbp"
sys.path.insert(0, str(TRIPLET))
from place_assay import OLD, RECEIVERS, SEED, cpu_seconds, inputs, invoke, sha, table
sys.path.insert(0, str(OLD))
import evaluate as old

ARMS = ("typed", "perm")
CASES = [(r, n) for r in RECEIVERS for n in range(4, 12)]


def prepare(work):
    inputs(work)
    archive = TRIPLET / "blind_forks.tar.xz"
    assert sha(archive) == (TRIPLET / "ARCHIVE.sha256").read_text().split()[0]
    with tarfile.open(archive) as tar:
        tar.extractall(work, filter="data")
    paths = sorted([*work.glob("plans/*.npz"), *work.glob("maps/*/*.map"),
                    *work.glob("affinity/*.pdbqt"), *work.glob("*/receptor.pdbqt"),
                    *work.glob("capture/*"), *work.glob("pools/base/*.conf")])
    (work / "source_hashes.json").write_text(json.dumps({
        "upstream_candidate_archive_sha256": sha(archive),
        "files": {str(p.relative_to(work)): sha(p) for p in paths}
    }, indent=2) + "\n")


def atom_types(path, serials):
    keyed = {int(l[6:11]): l[77:].strip() for l in path.read_text().splitlines()
             if l.startswith(("ATOM", "HETATM")) and l[77:].strip() not in ("H", "HD", "HS")}
    assert set(keyed) == set(serials)
    serial_order = np.sort(serials)
    ordered = np.array([keyed[s] for s in serial_order])
    perm = np.roll(ordered, -(len(ordered) // 2))
    idx = np.searchsorted(serial_order, serials)
    return ordered[idx], perm[idx]


def field_score(xyz, types, fields):
    values = np.zeros(xyz.shape[:2])
    outside = np.zeros(xyz.shape[:2], dtype=bool)
    for kind in np.unique(types):
        indices = types == kind
        x = xyz[:, indices]
        g = fields[kind]
        limited = np.clip(x, [axis[0] for axis in g.grid], [axis[-1] for axis in g.grid])
        miss = np.abs(limited - x).sum(axis=-1)
        raw = g(limited)
        positive = raw > 0
        raw[positive] *= 10 / (10 + raw[positive])
        values[:, indices] = raw + 1e6 * miss
        outside[:, indices] = miss > 0
    return values.sum(axis=1), outside


def rank(work, arms=ARMS):
    rows = []
    directory = work / "ranks"
    directory.mkdir(exist_ok=True)
    for rec in RECEIVERS:
        start, cpu = time.perf_counter(), time.process_time()
        fields = {t: old.read_map(work / "maps" / rec / f"receptor.{t}.map")
                  for t in ("A", "C", "N", "OA")}
        rows.append(dict(receiver=rec, compound="", action="load_type_fields", nposes="", natom="",
                         changed_atoms="", out_of_bounds_atom_positions="", winner="", queries=0,
                         wall_s=time.perf_counter()-start, cpu_s=time.process_time()-cpu))
        for n in range(4, 12):
            stem = f"{rec}_{n}"
            with np.load(work / "plans" / (stem + ".npz")) as bank:
                start, cpu = time.perf_counter(), time.process_time()
                xyz, q = bank["xyz"], bank["q"]
                rot, center = bank["matched_rotation"], bank["center"]
                poses = (xyz - xyz.mean(0)) @ rot.transpose(0, 2, 1) + center[:, None]
                real, perm = atom_types(work / "affinity" / f"{n}.pdbqt", bank["serials"])
                assert len(poses) == 1024
                common = dict(receiver=rec, compound=n, nposes=len(poses), natom=len(real),
                              changed_atoms=int(np.sum(real != perm)))
                rows.append(dict(**common, action="restore_poses", winner="", queries=0,
                                 out_of_bounds_atom_positions="", wall_s=time.perf_counter()-start,
                                 cpu_s=time.process_time()-cpu))
                scores = {"true_types": real, "perm_types": perm}
                for arm in arms:
                    types = real if arm == "typed" else perm
                    start, cpu = time.perf_counter(), time.process_time()
                    scores[arm], outside = field_score(poses, types, fields)
                    i = int(np.argmin(scores[arm]))
                    rotation = Rotation.from_matrix(rot[i])
                    origin = rotation.apply(q[:3] - xyz.mean(0)) + center[i]
                    quaternion = rotation * Rotation.from_quat(q[[4, 5, 6, 3]])
                    pose = np.r_[origin, quaternion.as_quat()[[3, 0, 1, 2]]]
                    (directory / f"{stem}_{arm}.input").write_text(
                        "1 1 " + " ".join(format(v, ".17g") for v in pose) + "\n")
                    rows.append(dict(**common, action=arm, winner=i, queries=len(poses)*len(real),
                                     out_of_bounds_atom_positions=int(outside.sum()),
                                     wall_s=time.perf_counter()-start, cpu_s=time.process_time()-cpu))
                np.savez_compressed(directory / (stem + ".npz"), **scores)
    table(work / "rank_cost.csv", rows)


def run(work, exe, arms=ARMS):
    recorded = []
    for j, (rec, n) in enumerate(CASES):
        for arm in arms if j % 2 == 0 else reversed(arms):
            stem = f"{rec}_{n}"
            p = work / "pools" / arm / stem
            row = invoke(work, exe, rec, n, arm, {
                "FD_RIGID_INTERVENTION": work / "ranks" / f"{stem}_{arm}.input",
                "FD_RIGID_LOG": p.with_suffix(".stages"), "FD_RIGID_COORDS": p.with_suffix(".xyz"),
                "FD_RIGID_HISTORY": p.with_suffix(".history")})
            np.testing.assert_array_equal(
                np.fromstring(p.with_suffix(".stages").read_text().splitlines()[0][7:], sep=" "),
                np.loadtxt(work / "capture" / (stem + ".q"))[2:])
            recorded.append(row)
    for arm in arms:
        table(work / f"runs-{arm}.csv", [r for r in recorded if r["arm"] == arm])


def official_high(work, exe, counter):
    """Release executable followed by a byte-confirmed upstream observer replay per cell."""
    assert sha(exe) == "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"
    target = work / "pools/official_high"
    target.mkdir(exist_ok=True)
    rows, costs = [], []
    for rec, n in CASES:
        stem = f"{rec}_{n}"
        conf = target / (stem + ".conf")
        conf.write_text((work / f"pools/high/{stem}.conf").read_text().replace("pools/high/", "pools/official_high/"))
        costs_cell = {}
        p = conf.with_suffix(".pdbqt")
        for tag, program in (("release", exe), ("upstream_count", counter)):
            cmd = [str(program), "--config", str(conf.relative_to(work))]
            if tag == "upstream_count":
                cmd.extend(["--out", str(target / "replay.pdbqt")])
            start, cpu = time.perf_counter(), cpu_seconds()
            log = subprocess.run(cmd, cwd=work, env={k:v for k,v in os.environ.items() if not k.startswith("FD_")},
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout
            costs_cell[tag] = (time.perf_counter()-start, cpu_seconds()-cpu)
            (target / f"{stem}.{tag}.log").write_bytes(log)
        assert (target / "replay.pdbqt").read_bytes() == p.read_bytes()
        (target / "replay.pdbqt").unlink()
        replay = (target / f"{stem}.upstream_count.log").read_text()
        counts = re.findall(r"Search trajectory:[^\n]+", replay)
        channels = re.findall(r"GRID_CHANNEL derivative=(\d+) n=(\d+)", replay)
        assert len(counts) == 3 and len(channels) == 2
        # Legacy auditor expects counters in the same log. Preserve both original logs separately.
        conf.with_suffix(".log").write_bytes((target / f"{stem}.release.log").read_bytes() +
                                            ("\n# byte-identical upstream replay counters\n" + "\n".join(counts) + "\n").encode())
        counts_by = lambda label: sum(int(re.search(r"\b" + label + r"=(\d+)", s)[1]) for s in counts)
        rows.append(dict(arm="official_high", receiver=rec, compound=n, seed=SEED,
                         wall_s=costs_cell["release"][0], cpu_s=costs_cell["release"][1],
                         hunt_evals=counts_by("hunt_evals"), second_evals=counts_by("second_evals"),
                         steps=counts_by("hunt_calls"), modes=p.read_text().count("ENDMDL"),
                         sha256=sha(p), executable_sha256=sha(exe)))
        costs.append(dict(receiver=rec, compound=n, executable_sha256=sha(counter),
                          byte_equal=1, channel_interpolations=sum(int(v) for _,v in channels),
                          wall_s=costs_cell["upstream_count"][0], cpu_s=costs_cell["upstream_count"][1]))
    table(work / "runs-official_high.csv", rows)
    table(work / "official_high_replay.csv", costs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare","rank","run","official-high","pipeline"))
    p.add_argument("work", type=lambda s: Path(s).resolve())
    p.add_argument("--exe", type=lambda s: Path(s).resolve())
    p.add_argument("--counter", type=lambda s: Path(s).resolve())
    p.add_argument("--arm", choices=ARMS)
    a = p.parse_args()
    if a.action == "prepare": prepare(a.work)
    elif a.action == "rank": rank(a.work)
    elif a.action == "run": run(a.work,a.exe)
    elif a.action == "pipeline":
        if not a.arm: p.error("pipeline requires --arm")
        rank(a.work, (a.arm,))
        run(a.work, a.exe, (a.arm,))
    else: official_high(a.work,a.exe,a.counter)
