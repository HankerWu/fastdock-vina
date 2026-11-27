"""Compare all prespecified scorers on identical prepared coordinates and pose pools."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import resource
import runpy
import tarfile
import tempfile
import time

import numpy as np

from scorer import Receptor, atoms, score_pose

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
FROZEN = HERE.parent
ARMS = ("native", "fixed", "minimum", "uniform", "ensemble", "rotor_only")
UPSTREAM = BENCH / "geometry_comparison/upstream"


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def unpack(archive, work):
    with tarfile.open(archive) as f:
        f.extractall(work, filter="data")


def blocks(text):
    return [b for b in text.split("ENDMDL") if b.strip()]


def archive_cases(root):
    groups = defaultdict(list)
    for row in read(FROZEN / "angle/poses.csv"):
        groups[row["protocol"], row["design"], row["query"], row["receiver"], row["seed"]].append(row)
    for panel, inp in [("cognate", BENCH / "geometry_comparison/pose/inputs.tar.xz"),
                       ("reciprocal", BENCH / "cubic_grid/human_panel/prepared.tar.xz")]:
        sites = root / panel
        unpack(inp, sites)
        for design in ("fixed", "matched", "per_complex"):
            study = UPSTREAM / panel
            if design != "fixed":
                study /= design
            with tarfile.open(study / "raw.tar.xz") as t:
                for cell in read(study / "cells.csv"):
                    if cell["arm"] != "native":
                        continue
                    q, r, seed = (cell[k] for k in ("query","receiver","seed"))
                    poses = t.extractfile(f"{q}-{r}-{seed}-native.pdbqt").read()
                    assert hashlib.sha256(poses).hexdigest() == cell["pool_sha256"]
                    rows = groups[panel, design, q, r, seed]
                    base = dict(panel=panel, design=design, target=cell["target"], query=q, receiver=r,
                                seed=seed, cap=rows[0]["cap"])
                    yield base, sites, poses.decode(), rows
    sites = root / "transfer"
    unpack(FROZEN / "transfer/inputs.tar.xz", sites)
    unpack(FROZEN / "transfer/outputs.tar.xz", sites)
    for cell in read(FROZEN / "transfer/cells.csv"):
        q, r, seed = (cell[k] for k in ("query","receiver","seed"))
        poses = (sites / "pools" / f"{q}-{r}-{seed}.pdbqt").read_bytes()
        assert hashlib.sha256(poses).hexdigest() == cell["pool_sha256"]
        rows = [dict(**row, scale=1/float(row["divisor"])) for row in read(FROZEN/"transfer/modes.csv")
                if (q,r,seed) == (row["query"],row["receiver"],row["seed"])]
        base = dict(panel="transfer", design="cross", target=cell["target"], query=q, receiver=r, seed=seed, cap=60000)
        yield base, sites, poses.decode(), rows
    sites = root / "esr2"
    unpack(HERE / "esr2.tar.xz", sites)
    count = runpy.run_path(str(BENCH / "terminal_rotor_score/independent/verify.py"))["independently_count"]
    for path in sorted(sites.glob("*/*-*-*.pdbqt")):
        q, cap, seed = path.stem.split("-")
        meta = json.loads((sites/q/"preparation.json").read_text())
        data = path.read_bytes()
        assert hashlib.sha256(data).hexdigest() == json.loads(path.with_suffix(".json").read_text())["pool_sha256"]
        scale = 1/(1+.05846 * count((sites/q/"ligand.pdbqt").read_bytes()))
        rows = []
        for block in blocks(data.decode()):
            coords = atoms(block)
            xyz = np.array([coords[i].xyz for i in meta["serial"]])
            rmsd = float(np.sqrt(((xyz - np.array(meta["mapped"]))**2).sum(axis=-1).mean(axis=-1)).min())
            rows.append(dict(scale=scale,rmsd_A=rmsd))
        base = dict(panel="esr2", design=cap, target=meta["target"], query=q, receiver=q, seed=seed, cap=cap)
        yield base, sites, data.decode(), rows


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    receptors, modes, cells, groups, costs = {}, [], [], defaultdict(list), []
    start = time.perf_counter()
    max_fixed_diff = max_radial_diff = max_ensemble_diff = max_fixed_node_diff = max_phase_diff = 0.
    with tempfile.TemporaryDirectory() as tmp:
        for case, sites, poses, table in archive_cases(Path(tmp)):
            receptor_path = sites / case["receiver"] / "receptor.pdbqt"
            ligand = (sites / case["query"] / "ligand.pdbqt").read_text()
            if receptor_path not in receptors:
                text = receptor_path.read_text()
                a = time.perf_counter()
                rec = Receptor(text)
                setup_s = time.perf_counter() - a
                receptors[receptor_path] = rec, Receptor(text, n=144), Receptor(text, phase=.5)
                costs.append(dict(operation="prepare_receptor", panel=case["panel"], receiver=case["receiver"],
                                  query="",seed="",cap="",poses=0,wall_s=setup_s))
            rec, fine, shifted = receptors[receptor_path]
            pool = []
            duration = 0.
            assert len(blocks(poses)) == len(table)
            for rank, (block, row) in enumerate(zip(blocks(poses), table, strict=True),1):
                start_scoring = time.perf_counter()
                terms = score_pose(rec, ligand, block)
                duration += time.perf_counter() - start_scoring
                fine_terms = score_pose(fine, ligand, block)
                shifted_terms = score_pose(shifted, ligand, block)
                max_ensemble_diff = max(max_ensemble_diff, abs(terms["raw_ensemble"]-fine_terms["raw_ensemble"]))
                max_phase_diff = max(max_phase_diff, abs(terms["raw_ensemble"]-shifted_terms["raw_ensemble"]))
                max_fixed_node_diff = max(max_fixed_node_diff, abs(terms["raw_minimum"]-fine_terms["raw_minimum"]))
                if "raw_penalty" in row:
                    assert int(row["contacts"]) == terms["contacts"]
                    max_fixed_diff = max(max_fixed_diff, abs(terms["raw_fixed"]-float(row["raw_penalty"])))
                    max_radial_diff = max(max_radial_diff, abs(terms["raw_radial"]-float(row["raw_radial"])))
                scale = float(row["scale"])
                native = float(next(l.split()[3] for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")))
                if "native_score" in row:
                    assert native == float(row["native_score"])
                entry = dict(**case, rank=rank, scale=scale, rmsd_A=row["rmsd_A"], native=native, **terms)
                for arm in ARMS[1:]:
                    entry[arm] = native + terms["raw_"+arm]*scale
                assert terms["raw_minimum"] <= min(terms["raw_fixed"],terms["raw_ensemble"])+1e-10
                pool.append(entry)
            modes.extend(pool)
            costs.append(dict(operation="score_pool", panel=case["panel"], receiver=case["receiver"],
                              query=case["query"],seed=case["seed"],cap=case["cap"],poses=len(pool),wall_s=duration))
            cell = dict(**case, modes=len(pool), oracle=min(float(s["rmsd_A"]) for s in pool))
            for arm in ARMS:
                order = sorted(pool, key=lambda p:(p[arm],p["rank"]))
                cell[arm+"_rank"] = order[0]["rank"]
                cell[arm+"_top1"] = float(order[0]["rmsd_A"])
                cell[arm+"_top5"] = min(float(p["rmsd_A"]) for p in order[:5])
            cells.append(cell)
            groups[case["panel"],case["design"]].append(cell)
    summaries=[]
    for (panel, design), group in groups.items():
        for arm in ARMS:
            summaries.append(dict(panel=panel,design=design,arm=arm,cells=len(group),
                                  top1=sum(c[arm+"_top1"]<=2 for c in group),
                                  top5=sum(c[arm+"_top5"]<=2 for c in group), oracle=sum(c["oracle"]<=2 for c in group),
                                  switched=sum(c[arm+"_rank"]!=c["native_rank"] for c in group),
                                  rescued=sum(c[arm+"_top1"]<=2<c["native_top1"] for c in group),
                                  lost=sum(c["native_top1"]<=2<c[arm+"_top1"] for c in group)))
    assert max(max_fixed_diff,max_radial_diff)<1e-9
    for filename, rows in (("modes",modes),("cells",cells),("summary",summaries),("cost",costs)):
        write(out/(filename+".csv"),rows)
    audit = dict(poses=len(modes),pools=len(cells),max_fixed_diff=max_fixed_diff,max_radial_diff=max_radial_diff,
                 max_72_to_144_raw_ensemble_diff=max_ensemble_diff,max_half_step_raw_ensemble_diff=max_phase_diff,
                 max_72_to_144_raw_minimum_diff=max_fixed_node_diff,wall_including_diagnostics_s=time.perf_counter()-start,
                 peak_rss_linux_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(json.dumps(audit,indent=2))
    (out/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    main(parser.parse_args().out)
