#!/usr/bin/env python3
"""Validate default parity and late Vina selection on the example docking fixture."""

import argparse
from itertools import combinations
import math
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COUNTS = re.compile(r"Late pose selection: saved=(\d+) refine=(\d+) refinement_calls=(\d+) retained=(\d+)")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dock(binary, directory, name, options=(), cpu=1, *, ligand=None):
    prefix = ROOT / "examples/basic"
    output = directory / f"{name}.pdbqt"
    command = [str(binary), "--receptor", str(prefix / "1iep_receptor.pdbqt"),
               "--ligand", str(ligand or prefix / "1iep_ligand.pdbqt"),
               "--config", str(prefix / "1iep_receptor.box.txt"), "--out", str(output),
               "--seed", "42", "--cpu", str(cpu), "--exhaustiveness", "3", "--max_evals", "1800",
               "--num_modes", "3", "--min_rmsd", "1", "--energy_range", "100"]
    completed = subprocess.run(command + list(options), text=True, capture_output=True, timeout=90)
    return completed, output


def parse(text, scale=None):
    poses = []
    for block in text.split("ENDMDL"):
        if not block.strip():
            continue
        lines = block.splitlines()
        atoms = tuple(line for line in lines if line.startswith(("ATOM  ", "HETATM")))
        key = "REMARK VINA RESULT:" if scale is None else "REMARK VINA PRO RESULT:"
        (result,) = [line for line in lines if line.startswith(key)]
        score = float(result.removeprefix(key).split()[0])
        native_rank = len(poses) + 1
        native_score = score
        if scale is not None:
            (native,) = [line for line in lines if line.startswith("REMARK VINA NATIVE:")]
            native_score, native_rank = native.split()[-2:]
            native_rank = int(native_rank)
            native_score = float(native_score)
            (rerank,) = [line for line in lines if line.startswith("REMARK VINA PRO GAUSS2 SCALE:")]
            require(float(rerank.split()[-1]) == scale, "Wrong rerank scale")
            require(not any(line.startswith("REMARK VINA RESULT:") for line in lines), "Found duplicate native result")
        require(atoms and math.isfinite(score) and math.isfinite(native_score), "Missing pose or score")
        poses.append((atoms, score, native_rank, native_score))
    require(poses and all(a[1] <= b[1] for a, b in zip(poses, poses[1:])), "Pose scores are not sorted")
    return poses


def assert_separation(poses):
    for (a, *_), (b, *_) in combinations(poses, 2):
        heavy_pairs = [(x, y) for x, y in zip(a, b) if x.split()[-1] not in {"H", "HD", "HS"}]
        squared_displacement = sum(
            (float(x[start:start + 8]) - float(y[start:start + 8])) ** 2
            for x, y in heavy_pairs for start in (30, 38, 46)
        )
        require(math.sqrt(squared_displacement / len(heavy_pairs)) >= 0.998,
                "Refined heavy-atom coordinates violate the 1 Å minimum")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    with tempfile.TemporaryDirectory() as workdir:
        workdir = Path(workdir)

        def run(name, options=(), cpu=1, engine=binary):
            proc, output = dock(engine, workdir, name, options, cpu=cpu)
            require(proc.returncode == 0, proc.stderr or proc.stdout)
            return proc.stdout, output.read_text()

        _, old = run("old", engine=baseline)
        _, default = run("default")
        _, explicit = run("native", ("--pose_selection", "native"))
        require(old == default == explicit, "Default/native output differs from the pre-change Vina output")

        options = ("--pose_selection", "late")
        log, late = run("late", options)
        (counts,) = COUNTS.findall(log)
        saved, refine, calls, retained = map(int, counts)
        require(saved > refine >= retained >= 1 and saved <= 9 and refine <= 3
                and refine <= calls <= 5 * refine, f"Invalid candidate bounds: {counts}")
        require(log.count("Computing Vina grid") == 1, "Expected one grid initialization")
        poses = parse(late)
        require(len(poses) == retained, "Retained pose count differs from the output")
        assert_separation(poses)

        thread_log, threaded = run("late-threaded", options, cpu=2)
        require(COUNTS.findall(thread_log) == [counts] and threaded == late, "Late output depends on worker count")

        both_log, both = run("late-and-rerank", options + ("--rerank_gauss2_scale", "0.5"))
        reranked = sorted(parse(both, scale=0.5), key=lambda pose: pose[2])
        require(COUNTS.findall(both_log) == [counts], "Gauss2 changed preceding selection counts")
        require([(x[0], x[2], x[3]) for x in reranked] == [(x[0], x[2], x[3]) for x in poses],
                "Gauss2 changed the late-selected native candidates")

        for name, extra, message in (
            ("invalid-selection", ("--pose_selection", "invalid"), "pose_selection must be native or late"),
            ("maps", options + ("--maps", "unused"), "late pose selection requires regular"),
            ("flex", options + ("--flex", "unused"), "late pose selection requires regular"),
            ("no-refine", options + ("--no_refine",), "late pose selection requires regular"),
            ("score-only", options + ("--score_only",), "late pose selection requires regular"),
            ("vinardo", options + ("--scoring", "vinardo"), "late pose selection requires regular"),
        ):
            result, output = dock(binary, workdir, name, extra)
            require(result.returncode != 0 and message in result.stderr and not output.exists(),
                    f"Unexpected handling for {name}: {result.stderr}")

    print(f"Default byte parity, deterministic late search, Gauss2 composition and invalid modes passed. "
          f"Saved {saved}, refined {refine} in {calls} calls, retained {retained}.")


if __name__ == "__main__":
    main()
