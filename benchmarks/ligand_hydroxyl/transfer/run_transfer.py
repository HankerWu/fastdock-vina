"""Reproduce the preregistered receptor preparation and paired CLI pose pools."""

import argparse
import hashlib
import itertools
import json
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "hydration_transfer"))
import panel as hydration

cross = hydration.cross
GROUPS = {
    "ERRg": {
        "6I63": "2OH",
        "6I61": "H3W",
        "6I62": "27N",
        "6I64": "H48",
        "6I65": "H3Z",
        "6I66": "H4Q",
        "6I67": "H42",
    },
    "HSP90a": {"2XJG": "XJG", "2XJJ": "L81", "2XJX": "XJX", "3TUH": "TUH"},
}
ASSAY = {
    "6I61": 569.0,
    "6I62": 157.5,
    "6I63": 98.6,
    "6I64": 48.7,
    "6I65": 390.4,
    "6I66": 182.5,
    "6I67": 201.4,
}


def missing_carbonyls(lines, query_xyz):
    """Complete an observed peptide plane using its two observed C bonds."""
    residues = {}
    for line in lines:
        if line.startswith("ATOM  ") and line[21] == "A" and line[16] in " A":
            residues.setdefault(line[22:27], {})[line[12:16].strip()] = line
    additions, details = [], []
    for (key, a), (_, b) in itertools.pairwise(residues.items()):
        if "O" in a or not {"C", "CA"}.issubset(a) or "N" not in b:
            continue
        c, ca, n = (np.array(cross.coords(l)) for l in (a["C"], a["CA"], b["N"]))
        d1, d2 = np.linalg.norm(c - ca), np.linalg.norm(c - n)
        if not (1.30 < d1 < 1.65 and 1.20 < d2 < 1.45):
            continue
        bisector = (c - ca) / d1 + (c - n) / d2
        o = c + 1.23 * bisector / np.linalg.norm(bisector)
        src = a["C"]
        additions.append(
            (
                src,
                src[:12]
                + " O  "
                + src[16:30]
                + "".join(f"{x:8.3f}" for x in o)
                + src[54:76]
                + " O  ",
            )
        )
        details.append(
            dict(
                residue=key.strip(),
                min_ligand_distance_A=float(
                    np.linalg.norm(query_xyz - o, axis=1).min()
                ),
            )
        )
    after = dict(additions)
    return [el for l in lines for el in ([l, after[l]] if l in after else [l])], details


def prepare(work):
    raw = work / "raw"
    raw.mkdir(exist_ok=True, parents=True)

    def one(task):
        target, entry, comp = task
        out = work / entry
        out.mkdir(exist_ok=True)
        pdb, ccd = raw / (entry + ".pdb"), raw / (comp + ".cif")
        sources = {}
        for path, url in [
            (pdb, f"https://files.rcsb.org/download/{entry}.pdb"),
            (ccd, f"https://files.rcsb.org/ligands/download/{comp}.cif"),
        ]:
            sources[url] = cross.fetch(path, url)
        lines = pdb.read_text().splitlines()
        hetero = [
            l
            for l in lines
            if l.startswith("HETATM")
            and l[21] == "A"
            and l[76:78].strip() not in ("H", "D")
        ]
        query = [l for l in hetero if l[17:20].strip() == comp]
        if len({l[22:27] for l in query}) != 1:
            raise ValueError(f"{entry}: missing or ambiguous ligand")
        refs = cross.dep_records(lines, comp)
        # Use a center averaged across deposited alternate conformers without
        # treating an alternate heavy-atom orientation as a second compound.
        center = np.mean([np.mean(list(r.values()), axis=0) for r in refs], axis=0)
        now = time.perf_counter()
        ligand = cross.ligand(ccd, refs, out / "ligand.pdbqt")
        ligand_s = time.perf_counter() - now
        now = time.perf_counter()
        receptor_lines, completed = missing_carbonyls(
            lines, np.array(list(refs[0].values()))
        )
        rec = cross.receptor(receptor_lines, out, entry, list(refs[0].values()))
        rec_s = time.perf_counter() - now
        cas = {
            l[22:27].strip() + ":" + l[17:20]: cross.coords(l)
            for l in lines
            if l.startswith("ATOM  ")
            and l[21] == "A"
            and l[12:16].strip() == "CA"
            and l[16] in " A"
        }
        others = {
            t: len({l[22:27] for l in hetero if l[17:20].strip() == t})
            for t in {l[17:20].strip() for l in hetero}
            if t != comp
        }
        nearest = {}
        for kind in others:
            pts = [cross.coords(l) for l in hetero if l[17:20].strip() == kind]
            nearest[kind] = float(cdist(pts, [cross.coords(l) for l in query]).min())
        occupancy = sorted({float(l[54:60]) for l in query})
        resolution = float(
            next(l.split()[3] for l in lines if l.startswith("REMARK   2 RESOLUTION."))
        )
        info = dict(
            target=target,
            ligand=comp,
            center=center.tolist(),
            ca=cas,
            completed_peptide_oxygens=completed,
            sources=sources,
            occupancy=occupancy,
            resolution=resolution,
            other_chain_A_residue_counts=others,
            nearest_hetero_A=nearest,
            **ligand,
            **rec,
        )
        return entry, info, dict(entry=entry, receptor_s=rec_s, ligand_s=ligand_s)

    tasks = [(t, e, c) for t, es in GROUPS.items() for e, c in es.items()]
    with ThreadPoolExecutor(2) as pool:
        rows = list(pool.map(one, tasks))
    (work / "metadata.json").write_text(
        json.dumps({e: m for e, m, _ in rows}, sort_keys=True)
    )
    hydration.write_csv(work / "preparation.csv", [r for _, _, r in rows])
    (work / "packages.json").write_text(
        json.dumps(
            {
                p: version(p)
                for p in [
                    "rdkit",
                    "meeko",
                    "gemmi",
                    "pdbfixer",
                    "openmm",
                    "numpy",
                    "scipy",
                ]
            },
            indent=2,
        )
        + "\n"
    )


def dock(work, binaries, jobs):
    (work / "runs").mkdir(exist_ok=True)
    meta = json.loads((work / "metadata.json").read_text())
    pairs = [
        (s, r)
        for entries in GROUPS.values()
        for s in entries
        for r in dict.fromkeys((next(iter(entries)), s))
    ]
    (work / "binaries.json").write_text(
        json.dumps(
            {
                k: {
                    "version": subprocess.check_output(
                        [str(v), "--version"], text=True
                    ).strip(),
                    "sha256": hashlib.sha256(v.read_bytes()).hexdigest(),
                }
                for k, v in binaries.items()
            },
            indent=2,
        )
        + "\n"
    )

    def one(task):
        (s, r), cap, seed, arm = task
        tag = f"{s}_{r}_{cap}_{seed}_{arm}"
        cmd = [
            str(binaries[arm]),
            "--receptor",
            f"{r}/receptor.pdbqt",
            "--ligand",
            f"{s}/ligand.pdbqt",
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
            f"runs/{tag}.pdbqt",
        ]
        for a, p in zip("xyz", meta[r]["center"]):
            cmd.extend(["--center_" + a, str(p), "--size_" + a, "24"])
        if arm == "pro":
            cmd.append("--search_counters")
        log = work / "runs" / (tag + ".log")
        cost = hydration.timed_process(cmd, work, log)
        content = log.read_text()
        ev = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", content)
        ef = re.findall(r"explicit_evals=(\d+)", content)
        return dict(
            query=s,
            receiver=r,
            cap=cap,
            seed=seed,
            engine=arm,
            **cost,
            search_evals=sum(int(x) + int(y) for x, y in ev) if ev else "",
            explicit_evals=sum(map(int, ef)) if ef else "",
            args=" ".join(cmd[1:]),
        )

    tasks = itertools.product(pairs, [3000, 12000], [43103, 43117], binaries)
    with ThreadPoolExecutor(jobs) as pool:
        rows = list(pool.map(one, tasks))
    hydration.write_csv(work / "cost.csv", rows)
    hydration.geometry(work)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("step", choices=["prepare", "dock", "geometry"])
    ap.add_argument("work", type=Path)
    ap.add_argument("--native", type=Path)
    ap.add_argument("--pro", type=Path)
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    work = args.work.resolve()
    if args.step == "prepare":
        prepare(work)
    elif args.step == "geometry":
        hydration.geometry(work)
    else:
        if args.native is None or args.pro is None:
            ap.error("dock needs both --native and --pro executables")
        dock(
            work,
            {"upstream": args.native.resolve(), "pro": args.pro.resolve()},
            args.jobs,
        )
