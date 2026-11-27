"""Prepare the preregistered cross-docks, run actual CLI binaries and score geometry."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "paired_crossdocking"))
import run as cross

PANELS = {
    "BRPF1": {"5C7N": "BMF", "5EPR": "5QY", "5EPS": "5QX", "5FG5": "5XF"},
    "EPHX2": {"4OD0": "2RV", "3WKE": "AUB", "4X6Y": "S94", "3ANT": "S82"},
}


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(work):
    from importlib.metadata import version

    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    meta, records = {}, []
    for target, entries in PANELS.items():
        receiver = next(iter(entries))
        for entry, comp in entries.items():
            directory = work / entry
            directory.mkdir(exist_ok=True)
            pdb, ccd = raw / f"{entry}.pdb", raw / f"{comp}.cif"
            hashes = {}
            for filename, url in [
                (pdb, f"https://files.rcsb.org/download/{entry}.pdb"),
                (ccd, f"https://files.rcsb.org/ligands/download/{comp}.cif"),
            ]:
                hashes[url] = cross.fetch(filename, url)
            lines = pdb.read_text().splitlines()
            refs = cross.dep_records(lines, comp)
            instances = {
                l[22:27]
                for l in lines
                if l.startswith("HETATM") and l[21] == "A" and l[17:20] == comp
            }
            if len(instances) != 1:
                raise ValueError(f"Ambiguous chain A component in {entry}")
            center = np.mean(list(refs[0].values()), axis=0).tolist()
            start = time.perf_counter()
            ligand = cross.ligand(ccd, refs, directory / "ligand.pdbqt")
            ligand_time = time.perf_counter() - start
            start = time.perf_counter()
            receptor = (
                cross.receptor(lines, directory, entry, list(refs[0].values()))
                if entry == receiver
                else {}
            )
            rec_time = time.perf_counter() - start if entry == receiver else 0
            ca = {
                l[22:27].strip() + ":" + l[17:20]: cross.coords(l)
                for l in lines
                if l.startswith("ATOM  ")
                and l[21] == "A"
                and l[12:16].strip() == "CA"
                and l[16] in " A"
            }
            meta[entry] = dict(
                target=target,
                ligand=comp,
                center=center,
                ca=ca,
                sources=hashes,
                **ligand,
                **receptor,
            )
            records.append(
                {"entry": entry, "ligand_s": ligand_time, "receptor_s": rec_time}
            )
    (work / "metadata.json").write_text(json.dumps(meta, sort_keys=True))
    write_csv(work / "preparation.csv", records)
    (work / "environment.json").write_text(
        json.dumps(
            {
                "packages": {
                    p: version(p)
                    for p in [
                        "rdkit",
                        "meeko",
                        "gemmi",
                        "numpy",
                        "pdbfixer",
                        "openmm",
                        "scipy",
                    ]
                },
                "python": sys.version,
                "platform": platform.platform(),
                "threads": {
                    k: os.environ.get(k)
                    for k in [
                        "OPENMM_CPU_THREADS",
                        "OPENBLAS_NUM_THREADS",
                        "OMP_NUM_THREADS",
                    ]
                },
            },
            indent=2,
        )
        + "\n"
    )


def timed_process(cmd, work, log):
    start = time.perf_counter()
    with log.open("w") as stream:
        proc = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT, cwd=work)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    if proc.returncode:
        raise RuntimeError(f"{' '.join(cmd)} failed, see {log}")
    return {
        "wall_s": time.perf_counter() - start,
        "cpu_s": usage.ru_utime + usage.ru_stime,
        "peak_rss_kib": usage.ru_maxrss,
    }


def dock(work, native, pro):
    directory = work / "runs"
    directory.mkdir(exist_ok=True)
    metadata = json.loads((work / "metadata.json").read_text())
    binaries = {"upstream": native.resolve(), "pro": pro.resolve()}
    builds = {
        arm: {
            "sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
            "version": subprocess.check_output(
                [str(binary), "--version"], text=True
            ).strip(),
        }
        for arm, binary in binaries.items()
    }
    (work / "binaries.json").write_text(json.dumps(builds, indent=2) + "\n")
    records = []
    for target, entries in PANELS.items():
        receiver = next(iter(entries))
        for entry, cap, seed, (arm, binary) in itertools.product(
            entries, [3000, 12000], [43103, 43117], binaries.items()
        ):
            ident = f"{entry}_{receiver}_{cap}_{seed}_{arm}"
            cmd = [
                str(binary),
                "--receptor",
                f"{receiver}/receptor.pdbqt",
                "--ligand",
                f"{entry}/ligand.pdbqt",
                "--exhaustiveness",
                "4",
                "--max_evals",
                str(cap),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--out",
                f"runs/{ident}.pdbqt",
            ]
            for axis, center in zip("xyz", metadata[receiver]["center"]):
                cmd.extend(["--center_" + axis, str(center), "--size_" + axis, "24"])
            if arm == "pro":
                cmd.append("--search_counters")
            log = directory / (ident + ".log")
            cost = timed_process(cmd, work, log)
            output = log.read_text()
            search = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", output)
            final = re.findall(r"explicit_evals=(\d+)", output)
            records.append(
                dict(
                    target=target,
                    query=entry,
                    receiver=receiver,
                    cap=cap,
                    seed=seed,
                    engine=arm,
                    search_evals=sum(int(a) + int(b) for a, b in search)
                    if search
                    else "",
                    explicit_evals=sum(map(int, final)) if final else "",
                    **cost,
                    args=" ".join(cmd[1:]),
                )
            )
    write_csv(work / "cost.csv", records)


def geometry(work):
    meta = json.loads((work / "metadata.json").read_text())
    rows = []
    for path in sorted((work / "runs").glob("*.pdbqt")):
        query, receiver, cap, seed, engine = path.stem.split("_")
        rotation, translate, ca_rmsd, n_ca = cross.alignment(
            meta[query], meta[receiver]
        )
        refs = np.asarray(meta[query]["mapped"]) @ rotation + translate
        for block in path.read_text().split("ENDMDL"):
            ls = block.splitlines()
            atoms = {
                int(l[6:11]): cross.coords(l)
                for l in ls
                if l.startswith(("ATOM  ", "HETATM"))
            }
            if not atoms:
                continue
            score = next(
                float(l.split()[3]) for l in ls if l.startswith("REMARK VINA RESULT:")
            )
            model = next(int(l.split()[1]) for l in ls if l.startswith("MODEL"))
            xyz = [atoms[i] for i in meta[query]["serial"]]
            rms = np.sqrt(np.square(refs - xyz).sum(axis=2).mean(axis=1)).min()
            rows.append(
                {
                    "target": meta[query]["target"],
                    "query": query,
                    "receiver": receiver,
                    "cap": cap,
                    "seed": seed,
                    "engine": engine,
                    "model": model,
                    "dry": score,
                    "rmsd": rms,
                    "ca_rmsd": ca_rmsd,
                    "n_ca": n_ca,
                }
            )
    write_csv(work / "geometry.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("step", choices=["prepare", "dock", "geometry"])
    p.add_argument("--native", type=Path)
    p.add_argument("--pro", type=Path)
    args = p.parse_args()
    if args.step == "prepare":
        prepare(args.work.resolve())
    elif args.step == "dock":
        if args.native is None or args.pro is None:
            p.error("dock requires --native and --pro")
        dock(args.work.resolve(), args.native, args.pro)
    else:
        geometry(args.work.resolve())
