"""Run the prespecified transfer cells and the two executable reference arms."""

import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from prepare_transfer import BENCH, DIRECTIONS, SYSTEMS, prepare


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def references(args, meta, pairs):
    # Count upstream using its existing source-only observer, requiring byte parity.
    counters = load_module("holdout_counter", BENCH / "archive_overlap/transfer/run_panel.py")
    raw = args.work / "raw"
    raw.mkdir(exist_ok=True)
    rows, poses = [], []
    registry = {}
    for name in ("upstream", "counter", "pro"):
        binary = getattr(args, name)
        registry[name] = dict(sha256=sha(binary), version=subprocess.check_output(
            [binary, "--version"], text=True).strip())
    (args.work / "binaries.json").write_text(json.dumps(registry, indent=2) + "\n")
    for r, q in pairs:
        for cap, seed, engine in itertools.product((3000, 12000), (43103, 43117), ("upstream", "pro")):
            key = f"{r}_{q}_{cap}_{seed}_{engine}"
            common = [
                "--receptor", str(args.base / r / "receptor.pdbqt"),
                "--ligand", str(args.base / q / "ligand.pdbqt"),
                "--cpu", "1", "--exhaustiveness", "4", "--num_modes", "9",
                "--max_evals", str(cap), "--seed", str(seed), "--energy_range", "100",
            ]
            for axis, value in zip("xyz", meta[r]["center"]):
                common += ["--center_" + axis, str(value), "--size_" + axis, "24"]
            file = raw / (key + ".pdbqt")
            wall, cpu, rss = counters.cli(getattr(args, engine), common, file, raw / (key + ".log"))
            counted = raw / (key + ".counted.pdbqt")
            log = raw / (key + ".counts")
            counters.cli(args.counter if engine == "upstream" else args.pro, common, counted,
                         log, counters=engine == "pro")
            if file.read_bytes() != counted.read_bytes():
                raise ValueError("Counter changed output: " + key)
            counted.unlink()
            counts = counters.sum_counts(log.read_text())
            rot, shift, backbone, n = prepare.base.alignment(meta[q], meta[r])
            refs = np.asarray(meta[q]["mapped"]) @ rot + shift
            scores = []
            for model in file.read_text().split("ENDMDL"):
                if "MODEL" not in model:
                    continue
                atoms = {int(line[6:11]): prepare.base.coords(line)
                         for line in model.splitlines() if line.startswith(("ATOM", "HETATM"))}
                xyz = np.array([atoms[i] for i in meta[q]["serial"]])
                rmsd = float(np.sqrt(np.square(xyz - refs).sum(-1).mean(-1)).min())
                result = next(line for line in model.splitlines() if line.startswith("REMARK VINA RESULT:"))
                scores.append(rmsd)
                poses.append(dict(key=key, rank=len(scores), rmsd=rmsd, energy=float(result.split()[3])))
            rows.append(dict(key=key, receiver=r, query=q, cap=cap, seed=seed, engine=engine,
                             top=scores[0], top5=min(scores[:5]), oracle=min(scores), output_poses=len(scores),
                             receptor_ca_rmsd=backbone, receptor_ca_pairs=n, wall_s=wall,
                             cpu_s=cpu, peak_rss_KiB=rss, **counts, output_sha256=sha(file)))
    save(args.work / "baseline.csv", rows)
    save(args.work / "baseline_poses.csv", poses)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["baseline", "rigid"])
    p.add_argument("base", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--arms", nargs="+", choices=("raw", "gradient", "sham"),
                   default=["raw", "sham", "gradient"])
    for opt in ("pro", "upstream", "counter", "probe"):
        p.add_argument("--" + opt, type=Path, required=opt == "pro")
    args = p.parse_args()
    args.work.mkdir(exist_ok=True, parents=True)
    meta = {c: json.loads((args.base / c / "preparation.json").read_text()) for c in SYSTEMS}
    pairs = [(r, q) for r, q in DIRECTIONS if meta[r]["receptor_status"] == "ready"]
    if args.action == "baseline":
        if not args.upstream or not args.counter:
            p.error("baseline requires --upstream and --counter")
        sys.path.insert(0, str(BENCH / "archive_overlap/transfer"))
        references(args, meta, pairs)
    else:
        if not args.probe:
            p.error("rigid requires --probe")
        # Legacy prepare.py imports an unrelated module also named 'run'.
        sys.modules.pop("run", None)
        sys.path.insert(0, str(BENCH / "rigid_lookahead"))
        panel = load_module("holdout_panel", BENCH / "rigid_lookahead/panel.py")
        panel.run_panel(args.base, meta, pairs, [43103, 43117], args.arms,
                        args.work, args.probe, args.pro, "transfer")


if __name__ == "__main__":
    main()
