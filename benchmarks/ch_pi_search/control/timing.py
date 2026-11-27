"""Time a fresh single-core Python batch from prepared GAA inputs through persisted winners."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import shlex
import subprocess
import sys
from statistics import mean
from tempfile import TemporaryDirectory
from time import perf_counter

HERE = Path(__file__).resolve().parent
FIELD = "FASTDOCK_RESEARCH_CHPI"


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def manifest(gaa_results):
    selected = [r for r in rows(gaa_results / "costs.csv") if r["cap"] == "24000" and r["search"] == "native"]
    if len(selected) != 12:
        raise ValueError("Expected twelve matching GAA cases")
    return selected


def tag(record):
    return "_".join(record[key] for key in ("query", "receptor", "seed"))


def blocks(pool):
    return [s.strip() + "\nENDMDL\n" for s in pool.read_text().split("ENDMDL") if s.strip()]


def batch(arm, work, source, binary, out):
    if FIELD in os.environ:
        raise ValueError("Unexpected pre-existing research field")
    if arm == "joint":
        # Load the same typing, entropy divisor and radial interaction used in the four-arm experiment.
        sys.path.insert(0, str(HERE.parent))
        from study import CHPi, DEPTH, WEIGHT, CENTER, WIDTH, CUTOFF, count_torsions, np
    models = {}
    for r in manifest(source):
        args = shlex.split(r["args"])
        pool = out / (tag(r) + ".pdbqt")
        args[args.index("--out") + 1] = str(pool)
        env = os.environ.copy()
        if arm == "joint":
            key = (r["query"], r["receptor"])
            if key not in models:
                ligand = (work / args[args.index("--ligand") + 1]).read_text()
                receptor = (work / args[args.index("--receptor") + 1]).read_text()
                potential = CHPi(receptor, ligand)
                if potential.l_arom or not potential.l_ch:
                    raise ValueError("Search field requires ligand CH and no aromatic ligand atoms")
                scale = WEIGHT / (1 + .05846 * count_torsions(ligand.encode()))
                sidecar = out / ("_".join(key) + ".field")
                with sidecar.open("w") as stream:
                    print(DEPTH * WEIGHT, CENTER, WIDTH, CUTOFF, file=stream)
                    print(*potential.l_ch, file=stream)
                    np.savetxt(stream, potential.r_arom, fmt="%.3f")
                models[key] = potential, scale, sidecar
            potential, scale, field = models[key]
            env[FIELD] = str(field)
        else:
            args[args.index("--max_evals") + 1] = "48000"
        subprocess.run([str(binary), *args], cwd=work, env=env, stdout=subprocess.DEVNULL, check=True)
        poses = blocks(pool)
        if arm == "joint":
            def score(pose):
                native = next(float(line.split()[3]) for line in pose.splitlines() if line.startswith("REMARK VINA RESULT:"))
                return native + scale * sum(potential.components(pose))
            first = min(poses, key=score)
        else:
            first = poses[0]
        (out / (tag(r) + ".winner")).write_text(first)


def validate(arm, out, source, controls):
    if arm == "joint":
        original = [r for r in rows(source / "costs.csv") if r["search"] == "field" and r["cap"] == "24000"]
        ranks = {tag(r): int(r["winner"]) for r in rows(source / "cells.csv")
                 if r["search"] == "field" and r["rerank"] == "field" and r["cap"] == "24000"}
    else:
        original = [r for r in rows(controls / "costs.csv") if r["engine"] == "pro" and r["cap"] == "48000"]
        ranks = {tag(r): 1 for r in original}
    identities = {}
    for r in original:
        case = tag(r)
        pool = out / (case + ".pdbqt")
        winner = (out / (case + ".winner")).read_bytes()
        if hashlib.sha256(pool.read_bytes()).hexdigest() != r["sha256"] or winner != blocks(pool)[ranks[case]-1].encode():
            raise ValueError(f"Output differs from previous predictions at {arm}/{case}")
        identities[case] = hashlib.sha256(winner).hexdigest()
    if len(identities) != 12:
        raise ValueError("Incomplete validation")
    return identities


def cpu_children():
    stats = resource.getrusage(resource.RUSAGE_CHILDREN)
    return stats.ru_utime + stats.ru_stime


def measure(args):
    expected = json.loads((args.source / "binary.json").read_text())["sha256"]
    if hashlib.sha256(args.binary.read_bytes()).hexdigest() != expected:
        raise ValueError("Search executable differs from four-arm experiment")
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    times, hashes = [], {}
    for i, order in enumerate([("joint", "native"), ("native", "joint")], 1):
        for position, arm in enumerate(order, 1):
            with TemporaryDirectory() as temp:
                dest = Path(temp)
                command = [sys.executable, __file__, str(args.work), str(args.source), str(args.binary),
                           str(dest), "--worker", arm]
                c0 = cpu_children()
                t0 = perf_counter()
                subprocess.run(command, check=True, env=env)
                elapsed = perf_counter() - t0
                cpu_s = cpu_children() - c0
                identities = validate(arm, dest, args.source, args.controls)
                if hashes.setdefault(arm, identities) != identities:
                    raise ValueError("Different winners between batches")
                times.append(dict(round=i, position=position, arm=arm, wall_s=elapsed, cpu_s=cpu_s,
                                  identical_pools=12, identical_winners=12))
    aggregate = {}
    for arm in ("joint","native"):
        samples = [r for r in times if r["arm"] == arm]
        aggregate[arm] = {k:mean(r[k] for r in samples) for k in ("wall_s","cpu_s")}
    result = dict(jobs_per_batch=12, pinned_logical_cpu=cpu, search_binary_sha256=expected,
                  batches=times, means=aggregate,
                  wall_ratio_joint_over_native=aggregate["joint"]["wall_s"]/aggregate["native"]["wall_s"],
                  winner_sha256=hashes)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work","source","binary","output"):
        p.add_argument(name, type=lambda s:Path(s).resolve())
    p.add_argument("--controls", type=Path, default=HERE/"results")
    p.add_argument("--worker", choices=("joint","native"))
    a = p.parse_args()
    if a.worker:
        batch(a.worker, a.work, a.source, a.binary, a.output)
    else:
        measure(a)
