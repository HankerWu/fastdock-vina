"""Paired flexible Factor Xa docking, using archived same-assay compounds."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import itertools
import json
import lzma
import os
from pathlib import Path
import re
import subprocess
import tarfile
from time import monotonic

HERE = Path(__file__).resolve().parent
SERIES = HERE.parent / "affinity_partition/fxa_series.csv"
ARCHIVE = HERE.parent / "frozen_core/affinity"
RECS = ("1FJS", "2P16")
SEEDS = (271828, 314159)


def write(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=list(dict.fromkeys(k for r in rows for k in r)),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def prepare(out):
    series = list(csv.DictReader(SERIES.open()))
    data = json.loads(lzma.open(ARCHIVE / "fxa_activity.json.xz").read())["activities"]
    observed = {r["activity_id"]: r for r in data}
    assert len(series) == len(observed) == len(data) == 34
    assert len({r["ligand"] for r in series}) == len(series)
    for row in series:
        m = observed[int(row["activity"])]
        assert (
            m["standard_relation"],
            m["standard_type"],
            m["standard_units"],
            m["document_chembl_id"],
            m["assay_chembl_id"],
            m["target_chembl_id"],
            m["molecule_chembl_id"],
            float(m["standard_value"]),
        ) == (
            "=",
            "Ki",
            "nM",
            "CHEMBL1133358",
            "CHEMBL661804",
            "CHEMBL244",
            row["ligand"],
            float(row["ki_nm"]),
        )
    included = [r["ligand"] for r in series if r["included"] == "1"]
    assert len(included) == 28
    source = ARCHIVE / "fxa_frozen.tar.xz"
    assert (
        hashlib.sha256(source.read_bytes()).hexdigest()
        == "c7786f8f37562a339298f390574ff8aed15131f5aaabc17c11562e80989a5475"
    )
    with tarfile.open(source) as pack:
        for rec in RECS:
            dest = out / "inputs" / rec
            dest.mkdir(exist_ok=True, parents=True)
            for mol in included:
                (dest / (mol + ".pdbqt")).write_bytes(
                    pack.extractfile(f"{rec}/{mol}.pdbqt").read()
                )
            receptor = (
                pack.extractfile(f"{rec}/{rec}.pdbqt").read().decode().splitlines()
            )
            names = ("CA", "CB", "CG", "OD1", "OD2")
            residue = {
                line[12:16].strip(): line
                for line in receptor
                if line.startswith("ATOM")
                and line[17:20] == "ASP"
                and line[21:22] == "A"
                and line[22:26].strip() == "189"
                and line[12:16].strip() in names
            }
            assert set(residue) == set(names)
            ca, cb, cg, od1, od2 = (residue[k] for k in names)
            i, j, k = [int(line[6:11]) for line in (ca, cb, cg)]
            flex = [
                "BEGIN_RES ASP A 189",
                "ROOT",
                ca,
                "ENDROOT",
                f"BRANCH {i} {j}",
                cb,
                f"BRANCH {j} {k}",
                cg,
                od1,
                od2,
                f"ENDBRANCH {j} {k}",
                f"ENDBRANCH {i} {j}",
                "END_RES ASP A 189",
            ]
            rigid = [line for line in receptor if line not in residue.values()]
            assert sorted(receptor) == sorted(rigid + list(residue.values()))
            (dest / "rigid.pdbqt").write_text("\n".join(rigid) + "\n")
            (dest / "flex.pdbqt").write_text("\n".join(flex) + "\n")
            rows = (HERE / (rec + "_box_atoms.pdb")).read_text().splitlines()
            assert len(rows) == {"1FJS": 38, "2P16": 34}[rec]
            center = [
                sum(float(line[p : p + 8]) for line in rows) / len(rows)
                for p in (30, 38, 46)
            ]
            (dest / "box.txt").write_text(
                "".join(
                    f"center_{axis} = {x:.4f}\nsize_{axis} = 26\n"
                    for axis, x in zip("xyz", center)
                )
            )
    return included


def run_pair(job):
    index, rec, mol, seed, arms, out = job
    folder = out / "inputs" / rec
    rows = []
    for arm in list(arms) if index % 2 == 0 else list(reversed(arms)):
        stem = out / "raw" / f"{rec}-{mol}-{seed}-{arm}"
        pose, log = stem.with_suffix(".pdbqt"), stem.with_suffix(".log")
        cmd = [
            arms[arm],
            "--receptor",
            folder / "rigid.pdbqt",
            "--flex",
            folder / "flex.pdbqt",
            "--ligand",
            folder / (mol + ".pdbqt"),
            "--config",
            folder / "box.txt",
            "--cpu",
            "1",
            "--exhaustiveness",
            "1",
            "--seed",
            str(seed),
            "--max_evals",
            "12000",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--search_counters",
            "--out",
            pose,
        ]
        start = monotonic()
        with log.open("w") as stream:
            proc = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        row = dict(
            receptor=rec,
            compound=mol,
            seed=seed,
            arm=arm,
            status=proc.returncode,
            wall_s=monotonic() - start,
            cpu_s=usage.ru_utime + usage.ru_stime,
            peak_kib=usage.ru_maxrss,
        )
        logtext = log.read_text()
        if proc.returncode == 0 and pose.exists():
            scores = re.findall(
                r"^REMARK VINA RESULT:\s+([-+.0-9]+)", pose.read_text(), re.M
            )
            row.update(modes=len(scores), top_score=float(scores[0]) if scores else "")
            for col in ("hunt_evals", "second_evals", "explicit_evals"):
                row[col] = sum(
                    int(v) for v in re.findall(r"\b" + col + r"=(\d+)", logtext)
                )
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ("baseline", "fixed", "out"):
        ap.add_argument(key, type=Path)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    out = args.out.resolve()
    (out / "raw").mkdir(parents=True, exist_ok=True)
    start = monotonic()
    mols = prepare(out)
    arms = {k: getattr(args, k).resolve() for k in ("baseline", "fixed")}
    meta = {
        "engine_sha256": {
            k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in arms.items()
        },
        "jobs": args.jobs,
        "preparation_seconds": monotonic() - start,
    }
    tasks = [
        (i, *condition, arms, out)
        for i, condition in enumerate(itertools.product(RECS, mols, SEEDS))
    ]
    start = monotonic()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = [r for pair in pool.map(run_pair, tasks) for r in pair]
    meta["dock_wall_s"] = monotonic() - start
    write(out / "runs.csv", rows)
    (out / "run.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    main()
