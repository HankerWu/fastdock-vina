"""A paired test on archived real structures and three measured within-assay series."""
import argparse
from collections import defaultdict
import csv
from functools import lru_cache
import json
import lzma
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

from model import ARMS, BENCH, Receptor, score

from evaluate import archive_cases, blocks, write  # common historical pool reader

from importlib.util import module_from_spec, spec_from_file_location
def load(name, path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

old_affinity = load("hydroxyl_assays", BENCH/"frozen_core/hydroxyl_evidence/affinity.py")
from assay import WEIGHT, assay_inputs, read

ALL_ARMS = ("native", *ARMS)


def make_row(meta, native, scale, terms):
    return dict(**meta,native=native,scale=scale,
                **{arm:native+scale*terms[arm] for arm in ARMS},
                **{k:v for k,v in terms.items() if k not in ARMS})


def rank_cell(meta, rows):
    cell = dict(**meta,poses=len(rows),oracle=min(r["rmsd"] for r in rows))
    for arm in ALL_ARMS:
        order = sorted(rows, key=lambda r:(r[arm],r["rank"]))
        cell[arm+"_rank"] = order[0]["rank"]
        cell[arm+"_top1"] = order[0]["rmsd"]
        cell[arm+"_top5"] = min(r["rmsd"] for r in order[:5])
    return cell


def pose_test(work, out):
    costs, cells, poses = [], [], []
    max_fixed_error = max_grid_error = max_phase_error = 0.
    old_fixed = {(r["panel"],r["design"],r["query"],r["receiver"],r["seed"],int(r["rank"])):float(r["raw_fixed"])
                 for r in read(BENCH/"frozen_core/hydroxyl_partition/results/modes.csv")}
    recs = {}
    for base, sites, posefile, table in archive_cases(work):
        rpath = sites/base["receiver"]/"receptor.pdbqt"
        query = (sites/base["query"]/"ligand.pdbqt").read_text()
        if rpath not in recs:
            start = time.perf_counter()
            recs[rpath] = Receptor(rpath.read_text())
            costs.append(dict(**base,stage="receptor",count=1,seconds=time.perf_counter()-start))
        rec = recs[rpath]
        pool = []
        elapsed = 0.
        for rank, (block,a) in enumerate(zip(blocks(posefile),table,strict=True),1):
            native = float(next(l for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")).split()[3])
            start = time.perf_counter()
            terms = score(rec,query,block)
            elapsed += time.perf_counter()-start
            fine = score(rec,query,block,n=144)
            shifted = score(rec,query,block,phase=.5)
            max_grid_error = max(max_grid_error,abs(fine["bonded"]-terms["bonded"]))
            max_phase_error = max(max_phase_error,abs(shifted["bonded"]-terms["bonded"]))
            key = tuple(base[k] for k in ("panel","design","query","receiver","seed"))+(rank,)
            max_fixed_error = max(max_fixed_error,abs(terms["full_fixed"]-old_fixed.pop(key)))
            pool.append(make_row(dict(**base,rank=rank,rmsd=float(a["rmsd_A"])),native,float(a["scale"]),terms))
        costs.append(dict(**base,stage="score",count=len(pool),seconds=elapsed))
        cells.append(rank_cell(base,pool))
        poses.extend(pool)
    if old_fixed or max_fixed_error > 1.e-9:
        raise ValueError(f"Historical fixed-H mismatch: {max_fixed_error}, unmatched {len(old_fixed)}")
    write(out/"poses.csv",poses)
    write(out/"cells.csv",cells)
    write(out/"cost.csv",costs)
    groups = defaultdict(list)
    for c in cells:
        groups[c["panel"],c["design"]].append(c)
    sums = []
    for (panel,design), group in groups.items():
        for arm in ALL_ARMS:
            sums.append(dict(panel=panel,design=design,arm=arm,cells=len(group),
                             top1=sum(c[arm+"_top1"]<=2 for c in group),
                             top5=sum(c[arm+"_top5"]<=2 for c in group),
                             oracle=sum(c["oracle"]<=2 for c in group),
                             wins=sum(c[arm+"_top1"]<=2<c["native_top1"] for c in group),
                             losses=sum(c["native_top1"]<=2<c[arm+"_top1"] for c in group)))
    write(out/"summary.csv",sums)
    return dict(cells=len(cells),poses=len(poses),max_raw_full_fixed_error=max_fixed_error,
                max_raw_bonded_72_144=max_grid_error,max_raw_bonded_half_node=max_phase_error)


def assay_test(out):
    records = {(r["target"],r["ligand"],r["receptor"],int(r["seed"])):r for r in
               read(BENCH/"terminal_rotor_score/results/poses.csv")}
    rows = []
    calls = []
    @lru_cache(maxsize=None)
    def receptor(text):
        start = time.perf_counter()
        rec = Receptor(text.decode())
        calls.append(dict(stage="receptor",seconds=time.perf_counter()-start))
        return rec
    for target,ligand,meas,receiver,seed,protein,query,posefile in assay_inputs():
        source = records.pop((target,ligand,receiver,seed))
        native = float(source["fixed"])
        scale = 1/(1+WEIGHT*float(source["fixed_rotors"]))
        first = posefile.decode().split("ENDMDL",1)[0]
        if abs(native-scale*float(source["saved_inter"]))>1.e-10:
            raise ValueError("Inconsistent reconstructed native score")
        if float(first.split("REMARK VINA RESULT:")[1].split()[0]) != float(source["saved_score"]):
            raise ValueError("Mismatched archived pose")
        rec = receptor(protein)
        start = time.perf_counter()
        terms = score(rec,query.decode(),first)
        calls.append(dict(stage="pose",seconds=time.perf_counter()-start))
        meta = dict(target=target,ligand=ligand,measurement=meas,receptor=receiver,seed=seed)
        rows.append(make_row(meta,native,scale,terms))
    if records:
        raise ValueError("Incomplete assay series")
    # Same published within-assay aggregation and leave-one-ligand diagnostics.
    old_affinity.ARMS = ALL_ARMS
    summary, edges = old_affinity.summarize([dict(**r,rotors=r["active_sites"]) for r in rows])
    write(out/"affinity_poses.csv",rows)
    write(out/"affinity_summary.csv",summary)
    write(out/"affinity_cost.csv",calls)
    with lzma.open(out/"affinity_pairs.csv.xz","wt",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=edges[0],lineterminator="\n")
        writer.writeheader()
        writer.writerows(edges)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    out = parser.parse_args().output
    out.mkdir(parents=True,exist_ok=True)
    start = time.perf_counter()
    with TemporaryDirectory() as td:
        audit = pose_test(Path(td),out)
    assay_test(out)
    audit["total_wall_s"] = time.perf_counter()-start
    audit["peak_RSS_KiB_linux"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (out/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
