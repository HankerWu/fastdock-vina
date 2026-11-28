"""Check the fixed-score domain and released 4OCZ decimal coordinates alongside fresh blind pools."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DOMAIN = HERE.parent / "fixed_score_domain"
EXT = HERE.parent / "cumulative/external"
sys.path.insert(0, str(DOMAIN))
from publication import coords, csvsave, energy, run


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def score_modes(source, folder, scorer, binary, work, metadata):
    records = []
    inputs = work / "fixtures/4OCZ"
    center = np.asarray(metadata["center"])
    reference = np.asarray(metadata["mapped"])
    for seed in (85231, 85237):
        file = (
            folder / f"after-native-{seed}.pdbqt"
            if source == "archived"
            else folder / f"external-4OCZ-4OCZ-{seed}-candidate.pdbqt"
        )
        blocks = file.read_text().split("ENDMDL")
        assert not blocks.pop().strip()
        for rank, block in enumerate(blocks, 1):
            stem = f"{source}-{scorer}-{seed}-{rank}"
            pose = work / "scores" / (stem + ".pdbqt")
            lines = [s for s in block.strip().splitlines() if not s.startswith("MODEL")]
            pose.write_text("\n".join(lines) + "\n")
            log = pose.with_suffix(".log")
            box = [
                item
                for dim, v in zip("xyz", center)
                for item in ("--size_" + dim, "24", "--center_" + dim, str(v))
            ]
            receipt = run(
                [
                    binary,
                    "--ligand",
                    inputs / "ligand.pdbqt",
                    "--receptor",
                    inputs / "receptor.pdbqt",
                    "--cpu",
                    "1",
                    "--seed",
                    str(seed),
                    "--spacing",
                    "0.375",
                    *box,
                    "--score_only",
                    "--unbound_energy",
                    "0",
                    "--pose_coordinates",
                    pose,
                ],
                log,
            )
            assert receipt["exit"] == 0
            output = log.read_text()
            fixed = energy(output, "(1) Final Intermolecular") + energy(
                output, "(2) Final Total Internal"
            )
            atoms = {
                int(l[6:11]): coords(l)
                for l in lines
                if l.startswith(("ATOM  ", "HETATM"))
            }
            xyz = np.array([atoms[i] for i in metadata["serial"]])
            excess = np.maximum(0, np.abs(xyz - center) - 12)
            records.append(
                dict(
                    source=source,
                    scorer=scorer,
                    seed=seed,
                    rank=rank,
                    fixed_score=energy(output, "Estimated Free Energy"),
                    fixed_raw=fixed,
                    printed_score=float(
                        next(
                            l for l in lines if l.startswith("REMARK VINA RESULT")
                        ).split()[3]
                    ),
                    residual=fixed - energy(block, "REMARK INTER + INTRA:"),
                    max_excess_A=excess.max(),
                    rmsd=np.sqrt(
                        np.mean(np.sum((reference - xyz) ** 2, axis=2), axis=1)
                    ).min(),
                    coordinate_sha256=hashlib.sha256(pose.read_bytes()).hexdigest(),
                    **receipt,
                )
            )
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("prefix", "stable", "candidate", "work"):
        p.add_argument(key, type=Path)
    p.add_argument(
        "--existing",
        action="store_true",
        help="reuse completed domain and blind runs in work",
    )
    a = p.parse_args()
    work = a.work.resolve()
    work.mkdir(exist_ok=True)
    domain, fresh = work / "domain", work / "fresh"
    if not a.existing:
        subprocess.run(
            [sys.executable, DOMAIN / "run.py", a.prefix, a.candidate, domain],
            check=True,
        )
        subprocess.run(
            [sys.executable, EXT / "run.py", "candidate", a.candidate, fresh],
            check=True,
        )
    assert (domain / "summary.csv").read_bytes() == (
        DOMAIN / "results/summary.csv"
    ).read_bytes()
    with tarfile.open(DOMAIN / "publication/raw.tar.xz") as t:
        t.extractall(work / "archived", filter="data")
    (work / "scores").mkdir()
    (work / "fixtures").symlink_to(fresh / "fixtures", target_is_directory=True)
    meta = json.loads((work / "fixtures/4OCZ/metadata.json").read_text())
    scores = []
    for source, folder, scorer, binary in [
        ("archived", work / "archived/raw", "stable", a.stable),
        ("archived", work / "archived/raw", "candidate", a.candidate),
        ("fresh", fresh / "raw", "candidate", a.candidate),
    ]:
        scores.extend(score_modes(source, folder, scorer, binary.resolve(), work, meta))
    original = {(r["seed"], r["rank"]): r for r in scores if r["scorer"] == "stable"}
    for r in scores:
        if r["source"] == "archived":
            old = original[r["seed"], r["rank"]]
            assert all(
                r[k] == old[k]
                for k in (
                    "fixed_score",
                    "fixed_raw",
                    "coordinate_sha256",
                    "max_excess_A",
                )
            )
    assert all(r["max_excess_A"] < 1e-6 and abs(r["residual"]) < 0.5 for r in scores)
    out = work / "results"
    out.mkdir()
    csvsave(out / "roundtrip.csv", scores)
    receipts = json.loads((domain / "runs.json").read_text())
    meta = dict(
        domain_calls=len(receipts),
        domain_summary_equal=True,
        mode_score_calls=len(scores),
        max_abs_residual=max(abs(r["residual"]) for r in scores),
        max_excess_A=max(r["max_excess_A"] for r in scores),
    )
    for src in ("fresh", "archived"):
        pool = [r for r in scores if r["source"] == src and r["scorer"] == "candidate"]
        for k in ("fixed_score", "printed_score", "rmsd"):
            r = min(pool, key=lambda a: (a[k], a["seed"], a["rank"]))
            meta[src + "_" + k + "_winner"] = {
                key: r[key]
                for key in ("seed", "rank", "rmsd", "fixed_score", "printed_score")
            }
    meta["binary_sha256"] = {
        key: hashlib.sha256(getattr(a, key).read_bytes()).hexdigest()
        for key in ("prefix", "stable", "candidate")
    }
    (out / "verification.json").write_text(json.dumps(meta, indent=2) + "\n")
    for n in ("cells.csv", "ranks.csv"):
        shutil.copyfile(fresh / n, out / n)
    with tarfile.open(out / "control_raw.tar.xz", "w:xz", dereference=True) as t:
        for d in (domain, fresh / "raw", work / "scores"):
            t.add(d, arcname=d.relative_to(work))
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
