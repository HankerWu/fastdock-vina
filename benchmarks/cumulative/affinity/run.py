"""Run the registered frozen panel from prepared inputs. One CPU per child."""

import argparse
import csv
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import tarfile
import time

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "geometry_comparison/affinity"
# This module is also named run.py. Load the historical protocol unambiguously.
import runpy

HISTORICAL = runpy.run_path(str(SOURCE / "run.py"))
STATES, SEEDS = HISTORICAL["STATES"], HISTORICAL["SEEDS"]
ASSET = "cf908040816ff85a7c5a7803e7f08fb3a8bafa8e2d40529b9d62d663b5d4f4c7"
LABELS = "bf2b5527aceb155f8ca3b69ba12ac0416d8fc23301114f6b535ca3ae41e41307"
ENGINES = ("official", "parent", "stable", "candidate")
STAGES = ("nominal", "low", "default", "matched")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_panel():
    path = SOURCE / "panel.csv"
    assert digest(path) == LABELS
    with path.open() as stream:
        return list(csv.DictReader(stream))


def prepare(work):
    archive = SOURCE / "inputs.tar.xz"
    pend = time.perf_counter()
    start = time.process_time()
    assert digest(archive) == ASSET
    with tarfile.open(archive) as dataset:
        dataset.extractall(work, filter="data")
    metadata = dict(
        archive_sha256=ASSET,
        panel_sha256=LABELS,
        wall_s=time.perf_counter() - pend,
        cpu_s=time.process_time() - start,
    )
    metadata["files"] = {p.name: digest(p) for p in sorted((work / "prepared").glob("*")) if p.is_file()}
    (work / "inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")


def cells(stage):
    for seed in SEEDS:
        for row in load_panel():
            for state in STATES[row["target"]]:
                if stage == "nominal" or (
                    row["assay"] == "CHEMBL2015090" and state == "4DJW" and seed == SEEDS[0]
                ):
                    yield state, row, seed


def run_cell(job, engine, binary, stage, cap, fingerprint):
    state, row, seed = job
    stem = Path("outputs") / engine / stage / f"{state}_{row['compound']}_{seed}"
    record = stem.with_suffix(".json")
    if record.exists():
        raise FileExistsError(record)
    output, log = stem.with_suffix(".pdbqt"), stem.with_suffix(".log")
    command = [
        str(binary),
        "--receptor",
        f"prepared/{state}.pdbqt",
        "--ligand",
        f"prepared/{row['compound']}.pdbqt",
        "--config",
        f"prepared/{state}.box",
        "--cpu",
        "1",
        "--seed",
        str(seed),
        "--out",
        str(output),
    ]
    if stage != "default":
        command.extend(
            [
                "--exhaustiveness",
                "2",
                "--max_evals",
                str(cap),
                "--energy_range",
                "100",
                "--num_modes",
                "9",
            ]
        )
    if engine in {"stable", "candidate"}:
        command.append("--search_counters")
    with log.open("wb") as trace:
        actions = [(os.POSIX_SPAWN_DUP2, trace.fileno(), fd) for fd in (1, 2)]
        start = time.perf_counter()
        pid = os.posix_spawn(command[0], command, os.environ, file_actions=actions)
        _, status, usage = os.wait4(pid, 0)
        elapsed = time.perf_counter() - start
    text = output.read_text() if output.exists() else ""
    scores = [float(l.split()[3]) for l in text.splitlines() if l.startswith("REMARK VINA RESULT:")]
    complete = text.startswith("MODEL 1\n") and text.endswith("ENDMDL\n") and bool(scores)
    returncode = os.waitstatus_to_exitcode(status)
    measurement = dict(
        engine=engine,
        stage=stage,
        receptor=state,
        target=row["target"],
        assay=row["assay"],
        compound=row["compound"],
        seed=seed,
        cap=cap if stage != "default" else None,
        command=command,
        binary_sha256=fingerprint,
        exit_status=returncode,
        valid=complete and returncode == 0,
        modes=len(scores),
        top_score=scores[0] if complete and returncode == 0 else None,
        output_sha256=digest(output) if output.exists() else None,
        wall_s=elapsed,
        user_s=usage.ru_utime,
        sys_s=usage.ru_stime,
        maxrss_kib=usage.ru_maxrss,
    )
    for counter in ("hunt_evals", "second_evals", "explicit_evals", "late_prescore"):
        matches = re.findall(r"\b" + counter + r"=(\d+)", log.read_text())
        measurement[counter] = sum(map(int, matches)) if matches else None
    record.write_text(json.dumps(measurement, indent=2) + "\n")


def total_cpu(engine, stage):
    measurements = []
    for state, row, seed in cells("low"):
        path = Path("outputs") / engine / stage / f"{state}_{row['compound']}_{seed}.json"
        r = json.loads(path.read_text())
        if not r["valid"]:
            raise ValueError(f"Incomplete CPU calibration: {path}")
        measurements.append(r["user_s"] + r["sys_s"])
    return sum(measurements)


def calibrate(engine):
    target = total_cpu("stable", "nominal")
    low, high = total_cpu(engine, "low"), total_cpu(engine, "nominal")
    if high > low:
        cap = round((3000 + 9000 * (target - low) / (high - low)) / 100) * 100
    else:
        cap = 3000 if abs(low - target) < abs(high - target) else 12000
    return min(48000, max(3000, cap))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--prepare", action="store_true")
    p.add_argument("--engine", choices=ENGINES)
    p.add_argument("--binary", type=Path)
    p.add_argument("--stage", choices=STAGES, nargs="+", default=["nominal", "low", "default"])
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    work = a.work.resolve()
    work.mkdir(exist_ok=True, parents=True)
    if a.prepare:
        prepare(work)
    if a.engine:
        binary = a.binary.resolve(strict=True)
        fingerprint = digest(binary)
        os.chdir(work)
        for stage in a.stage:
            cap = calibrate(a.engine) if stage == "matched" else 3000 if stage == "low" else 12000
            (Path("outputs") / a.engine / stage).mkdir(parents=True, exist_ok=True)
            with ThreadPoolExecutor(a.jobs) as ex:
                list(
                    ex.map(
                        lambda j: run_cell(j, a.engine, binary, stage, cap, fingerprint), cells(stage)
                    )
                )


if __name__ == "__main__":
    main()
