"""Run two search landscapes and replay two final rankers on each emitted pool."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(BENCH / "ch_pi_transfer/gba1"))
from evaluate import count_torsions, cross, write_csv
from experiment import timed_process
from model import CHPi, CUTOFF, CENTER, DEPTH, WEIGHT, WIDTH


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def jobs(gba, glucose):
    meta = json.loads((gba / "metadata.json").read_text())
    for r in read_csv(gba / "compute.csv"):
        if r["state"] != "protonated" or r["arm"] != "pro":
            continue
        src, dst = r["src"], r["dst"]
        rot, shift, _, _ = cross.alignment(meta[src], meta[dst])
        refs = np.array(meta[src]["mapped"]) @ rot + shift
        for cap in (12000, 24000):
            args = shlex.split(r["args"])
            args[args.index("--max_evals") + 1] = str(cap)
            arch = gba if cap == 12000 else gba / "budget_24000"
            yield (
                ("GBA1", src, dst, int(r["seed"]), cap),
                gba,
                args,
                arch / r["args"].split("--out ")[1].split()[0],
                meta[src]["protonated_serial"],
                refs,
            )
    meta = json.loads((glucose / "metadata.json").read_text())
    for r in read_csv(glucose / "compute.csv"):
        if r["arm"] != "pro":
            continue
        m = meta[r["case"]]
        args = shlex.split(r["args"])
        p = args[args.index("--out") + 1]
        panel = "GLUT3" if m["pdb"] == "4ZW9" else "GCK"
        yield (
            (panel, m["comp"], m["pdb"], int(r["seed"]), int(r["cap"])),
            glucose,
            args,
            glucose / p,
            m["serial"],
            np.array(m["mapped"]),
        )


def run_cases(cases, binary, out):
    if "FASTDOCK_RESEARCH_CHPI" in os.environ:
        raise ValueError("Clear FASTDOCK_RESEARCH_CHPI before the native controls")
    out.mkdir(parents=True, exist_ok=True)
    run = out / "pools"
    run.mkdir(exist_ok=True)
    sidecars = out / "typing"
    sidecars.mkdir(exist_ok=True)
    costs, modes, cells = [], [], []
    if len(cases) != len({c[0] for c in cases}):
        raise ValueError("Duplicate structural condition")
    for j, (key, work, args, archive, serial, refs) in enumerate(cases):
        panel, ligand, receptor, seed, cap = key
        case = dict(panel=panel, query=ligand, receptor=receptor, seed=seed, cap=cap)
        query = (work / args[args.index("--ligand") + 1]).read_text()
        protein = (work / args[args.index("--receptor") + 1]).read_text()
        potential = CHPi(protein, query)
        assert not potential.l_arom and potential.l_ch
        scale = 1 / (1 + 0.05846 * count_torsions(query.encode()))
        sidecar = sidecars / f"{panel}_{ligand}_{receptor}.txt"
        with sidecar.open("w") as s:
            print(DEPTH * WEIGHT, CENTER, WIDTH, CUTOFF, file=s)
            print(*potential.l_ch, file=s)
            np.savetxt(s, potential.r_arom, fmt="%.3f")
        for search in ("native", "field") if j % 2 == 0 else ("field", "native"):
            tag = "_".join(map(str, (*key, search)))
            pdbqt = run / (tag + ".pdbqt")
            log = run / (tag + ".log")
            options = args.copy()
            options[options.index("--out") + 1] = str(pdbqt)
            cmd = [str(binary), *options]
            if search == "field":
                cmd = ["env", "FASTDOCK_RESEARCH_CHPI=" + str(sidecar), *cmd]
            timing = timed_process(cmd, work, log)
            trace = log.read_text()
            counters = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", trace)
            explicit = re.findall(r"explicit_evals=(\d+)", trace)
            if len(counters) != 4:
                raise ValueError(f"Cannot verify the four search trajectories in {log}")
            identity = dict(**case, search=search)
            costs.append(
                dict(
                    **identity,
                    **timing,
                    search_evals=sum(int(a) + int(b) for a, b in counters),
                    explicit_evals=sum(map(int, explicit)),
                    sha256=digest(pdbqt),
                    archive_equal=(pdbqt.read_bytes() == archive.read_bytes())
                    if search == "native" and archive
                    else "",
                    args=shlex.join(args),
                )
            )
            pool = []
            start = time.perf_counter()
            for rank, block in enumerate(
                (b for b in pdbqt.read_text().split("ENDMDL") if b.strip()), 1
            ):
                lines = block.splitlines()
                xyz = {
                    int(l[6:11]): cross.coords(l)
                    for l in lines
                    if l.startswith(("ATOM  ", "HETATM"))
                }
                ligand_xyz = np.array([xyz[i] for i in serial])
                rmsd = float(
                    np.sqrt(((ligand_xyz - refs) ** 2).sum(axis=-1).mean(axis=-1)).min()
                )
                native = float(
                    next(
                        l.split()[3]
                        for l in lines
                        if l.startswith("REMARK VINA RESULT:")
                    )
                )
                fwd, rev = potential.components(block)
                assert rev == 0
                row = dict(
                    **identity,
                    rank=rank,
                    rmsd=rmsd,
                    native=native,
                    raw_chpi=fwd,
                    factor=scale,
                    corrected=native + WEIGHT * scale * fwd,
                )
                modes.append(row)
                pool.append(row)
            costs[-1]["analyze_pool_s"] = time.perf_counter() - start
            if not pool:
                raise ValueError(f"Missing released poses for {tag}")
            best = min(p["rmsd"] for p in pool)
            for rerank in ("native", "field"):
                ordered = (
                    pool
                    if rerank == "native"
                    else sorted(pool, key=lambda p: p["corrected"])
                )
                first = ordered[0]
                cells.append(
                    dict(
                        **identity,
                        rerank=rerank,
                        output_count=len(pool),
                        winner=first["rank"],
                        top1=first["rmsd"],
                        best5=min(p["rmsd"] for p in ordered[:5]),
                        pool_min=best,
                        regret=first["rmsd"] - best,
                    )
                )
    for name, rows in (("costs", costs), ("modes", modes), ("cells", cells)):
        write_csv(out / (name + ".csv"), rows)
    (out / "binary.json").write_text(
        json.dumps({"sha256": digest(binary)}, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gba", "glucose", "binary", "out"):
        parser.add_argument(name, type=lambda s: Path(s).resolve())
    a = parser.parse_args()
    cases = list(jobs(a.gba, a.glucose))
    if len(cases) != 70:
        raise ValueError("Incomplete historical manifest")
    run_cases(cases, a.binary, a.out)
