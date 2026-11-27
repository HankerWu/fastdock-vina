"""Run frozen five-arm FKBP12 panel. All paths in raw configs are relative to work."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import re
import resource
import subprocess
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
BENCH = ROOT / "benchmarks/cubic_grid"
SOURCE = BENCH / "fkbp_ad4"
ARMS = ("official", "base", "sham", "pivot", "high")
SEED = 524287
RECEIVERS = ("1FKG", "1FKH")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prepare(work):
    work.mkdir(parents=True, exist_ok=True)
    meta = json.loads((SOURCE / "results/metadata.json").read_text())
    source = BENCH / "human_panel/affinity-pools.tar.xz"
    with tarfile.open(source) as tar:
        for n, digest in meta["ligands"].items():
            name = f"affinity/{n}.pdbqt"
            tar.extract(name, work, filter="data")
            assert sha(work / name) == digest
    inputs = {str(source.relative_to(ROOT)): sha(source)}
    for rec, files in meta["maps"].items():
        archive = SOURCE / f"maps-{rec}.tar.xz"
        inputs[str(archive.relative_to(ROOT))] = sha(archive)
        with tarfile.open(archive) as tar:
            for name, digest in files.items():
                rel = f"maps/{rec}/{name}"
                tar.extract(rel, work, filter="data")
                assert sha(work / rel) == digest
    return inputs


def cpu_seconds():
    u = resource.getrusage(resource.RUSAGE_CHILDREN)
    return u.ru_utime + u.ru_stime


def run(work, binaries, only):
    inputs = prepare(work)
    known = {
        "official": "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"
    }
    paths = {k: p.resolve() for k, p in binaries.items() if k in only}
    digests = {k: sha(p) for k, p in paths.items()}
    if "official" in digests:
        assert digests["official"] == known["official"]
    measured = []
    for i, (r, n) in enumerate(itertools.product(RECEIVERS, range(4, 12))):
        # Cycle method order across cases.
        for arm in ARMS[i % len(ARMS) :] + ARMS[: i % len(ARMS)]:
            if arm not in only:
                continue
            pool = Path("pools") / arm
            (work / pool).mkdir(parents=True, exist_ok=True)
            base = pool / f"{r}_{n}"
            conf = base.with_suffix(".conf")
            config = {
                "ligand": f"affinity/{n}.pdbqt",
                "maps": f"maps/{r}/receptor",
                "scoring": "ad4",
                "num_modes": 9,
                "energy_range": 100,
                "cpu": 1,
                "exhaustiveness": 3,
                "max_evals": 48000 if arm == "high" else 24000,
                "seed": SEED,
                "out": base.with_suffix(".pdbqt").as_posix(),
            }
            (work / conf).write_text("".join(f"{k} = {v}\n" for k, v in config.items()))
            args = [str(paths[arm]), "--config", str(conf)]
            if arm != "official":
                args += ["--search_counters"]
            t, c = time.perf_counter(), cpu_seconds()
            logfile = work / base.with_suffix(".log")
            with logfile.open("w") as f:
                subprocess.run(
                    args,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    check=True,
                    cwd=work,
                    env={
                        k: v for k, v in os.environ.items() if not k.startswith("FD_")
                    },
                )
            elapsed, cpu = time.perf_counter() - t, cpu_seconds() - c
            out = work / config["out"]
            log = logfile.read_text()
            stages = [
                dict(re.findall(r"(\w+)=(\d+)", l))
                for l in log.splitlines()
                if l.startswith("Search trajectory:")
            ]
            if arm != "official":
                assert len(stages) == 3 and "explicit_evals=0" in log
            measured.append(
                dict(
                    arm=arm,
                    receiver=r,
                    compound=n,
                    seed=SEED,
                    wall_s=elapsed,
                    cpu_s=cpu,
                    hunt_evals=(
                        sum(int(s["hunt_evals"]) for s in stages) if stages else ""
                    ),
                    second_evals=(
                        sum(int(s["second_evals"]) for s in stages) if stages else ""
                    ),
                    steps=sum(int(s["steps"]) for s in stages) if stages else "",
                    modes=out.read_text().count("ENDMDL"),
                    sha256=sha(out),
                    executable_sha256=digests[arm],
                )
            )
    for arm in only:
        table(work / f"runs-{arm}.csv", [r for r in measured if r["arm"] == arm])
    (work / "inputs.json").write_text(json.dumps(inputs, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "official", "base", "sham", "pivot"):
        p.add_argument(k, type=Path)
    p.add_argument("--arms", choices=ARMS, nargs="+", default=ARMS)
    a = p.parse_args()
    run(
        a.work.resolve(),
        dict(official=a.official, base=a.base, sham=a.sham, pivot=a.pivot, high=a.base),
        a.arms,
    )
