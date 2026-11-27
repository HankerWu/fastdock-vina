"""Stage the frozen inputs, generate matched banks, and run the preregistered search. No reference poses."""
import argparse
import csv
import io
import json
from pathlib import Path
import re
import sys
import tarfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
ROOT = BENCH.parent
TRIPLET = BENCH / "cavity_triplet"
sys.path[:0] = [str(TRIPLET), str(TRIPLET / "acid_pair"), str(ROOT / "tools"), str(HERE.parent / "engine")]
from pose_runs import COUNT_KEYS, digest, invoke, table
from run_acid import chemical_graph
from pocket_triangle import Pocket
from prepare_bank import write_bank

# Query, receiver, prepared state and the exact per-query seed set.
CASES = {
    "82I": ("6YQW", "6Y7L", (91813, 91823)),
    "3FMZ-HA": ("3FMZ", "1RBP", (91813, 93241)),
    "3FMZ-A": ("3FMZ", "1RBP", (91813, 93241)),
    "4GHI": ("4GHI", "5TBM", (93241, 93257)),
    "5TBM": ("5TBM", "4GHI", (93241, 93257)),
}
LIB = {
    "brd": BENCH / "torsion_proposals/fixtures/inputs.tar.xz",
    "human": BENCH / "cubic_grid/human_panel/prepared.tar.xz",
    "acid": BENCH / "cubic_grid/acid_states/acid-pools.tar.xz",
}

def write_json(f, data):
    f.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

def csvrows(f):
    with f.open() as src:
        return list(csv.DictReader(src))

def configs(p):
    configs = json.loads((p / "inputs.json").read_text())["cells"]
    return {name: {k: p / v if k in ("ligand", "receptor", "bonds") else v
                   for k, v in c.items() if k != "sha256"} for name,c in configs.items()}

def prepare(p):
    for sub in ("inputs", "raw", "banks"):
        (p / sub).mkdir(parents=True, exist_ok=True)
    result = {}
    for name, (query, rec, _) in CASES.items():
        pkg = "brd" if name == "82I" else "human"
        r = Path("inputs") / (rec + ".pdbqt")
        l = Path("inputs") / (name + "-ligand.pdbqt")
        g = l.with_suffix(".bonds")
        with tarfile.open(LIB[pkg]) as f:
            (p/r).write_bytes(f.extractfile(rec + "/receptor.pdbqt").read())
            center = json.load(f.extractfile(rec + "/preparation.json"))["center"]
        acid = name.startswith("3FMZ")
        with tarfile.open(LIB["acid" if acid else pkg]) as f:
            member = name + ".pdbqt" if acid else query + "/ligand.pdbqt"
            (p/l).write_bytes(f.extractfile(member).read())
        chemical_graph(p/l, p/g)
        result[name] = dict(ligand=str(l), receptor=str(r), bonds=str(g), center=center,
                            size=24, sha256={str(v):digest(p/v) for v in (l,r,g)})
        if name == "82I" or acid:
            expected = (json.loads((TRIPLET/"results/pose_inputs.json").read_text())["systems"][name]
                        if name == "82I" else json.loads((TRIPLET/"acid_pair/results/inputs.json").read_text())["cells"][name])
            assert all(digest(p / result[name][key]) == expected[key] for key in ("ligand", "receptor", "bonds"))
            assert center == expected["center"]
    write_json(p/"inputs.json", dict(cells=result, archives={str(v.relative_to(ROOT)):digest(v) for v in LIB.values()}))

def interpolation(log, official=False):
    lines = log.read_text()
    if official:
        values = re.findall(r"GRID_CHANNEL derivative=(\d+) n=(\d+)", lines)
    else:
        values = re.findall(r"FIELD_COUNT combined=0 derivative=(\d+) positions=(\d+)", lines)
    counts = {int(d):int(n) for d,n in values}
    assert set(counts) == {0,1}, str(log)
    return dict(scalar_fields=counts[0], gradient_fields=counts[1])

def controls(p, pro, official, counter):
    expected = json.loads((TRIPLET/"results/pose_inputs.json").read_text())["binaries"]["official"]
    assert digest(official) == expected
    costs, extra = [], []
    for name,c in configs(p).items():
        for seed in CASES[name][2]:
            case = f"{name}_{seed}"
            for arm, cap in (("pro",25000), ("pro-high",50000), ("official",25000), ("official-high",50000)):
                is_official = arm.startswith("official")
                stem = f"{case}_{arm}"
                pool = invoke(p, costs, stem, official if is_official else pro, c, seed, cap, official=is_official)
                if is_official:
                    check = invoke(p, extra, stem+"-count", counter, c, seed, cap, official=True)
                    assert check.read_bytes() == pool.read_bytes() and extra[-1]["hunt_evals"] > 0
                    costs[-1].update({k:extra[-1][k] for k in COUNT_KEYS})
                    counts = interpolation(check.with_suffix(".log"), True)
                else:
                    counts = interpolation(pool.with_suffix(".log"))
                costs[-1].update(counts)
    table(p/"control_cost.csv", costs)
    table(p/"official_count_cost.csv", extra)
    write_json(p/"control_binaries.json", {k:digest(v) for k,v in dict(pro=pro,official=official,counter=counter).items()})

def freeze(p, observer, only=None):
    costs, generation, compared = [], [], []
    for name,c in configs(p).items():
        for seed in CASES[name][2]:
            case = f"{name}_{seed}"
            if only is not None and case != only:
                continue
            base = p / "raw" / (case+"_source")
            invoke(p, costs, base.name, observer, c, seed, env=dict(
                FD_SOURCE_STEP=1, FD_SOURCE_TRAJECTORY=1, FD_SOURCE_XYZ=base.with_suffix(".xyz"),
                FD_SOURCE_Q=base.with_suffix(".q")))
            data = base.with_suffix(".xyz").read_text().splitlines()
            serials = np.array(data[0].split()[1:],int)
            assert data[1].split()[:3] == ["source","1","1"]
            x = np.array(data[1].split()[3:],float).reshape(-1,3)
            q = np.loadtxt(base.with_suffix(".q"))[2:]
            dest = p/"banks"/case
            np.savez_compressed(dest.with_suffix(".source.npz"), xyz=x, q=q, serials=serials)
            t0, c0 = time.perf_counter(), time.process_time()
            pocket = Pocket.from_pdbqt(c["receptor"], c["center"], c["size"])
            res = pocket.place(x, seed, arms=("matched",))
            assert len(res["center"]) == 1024
            chosen = int(res["matched_best"])
            bank = dest.with_suffix(".bank")
            poses = write_bank(bank, x, q, res["center"], res["matched_rotation"])
            np.savez_compressed(dest.with_suffix(".npz"), **res)
            dest.with_suffix(".input").write_text("1 1 " + " ".join(f"{v:.17g}" for v in poses[chosen])+"\n")
            generation.append(dict(case=case, wall_s=time.perf_counter()-t0,cpu_s=time.process_time()-c0,
                site_nearest=pocket.n_grid,site_burial=pocket.n_burial,shape_queries=len(x)*1024,
                candidates=1024,triangles=len(pocket.triangles),heavy=len(x),shape_index=chosen,
                bank_sha256=digest(bank)))
            if only is None and (name=="82I" or (name.startswith("3FMZ") and seed==91813)):
                old_name=f"{case}_t1_s1"
                archive = TRIPLET/("candidate_transforms.tar.xz" if name=="82I" else "acid_pair/raw.tar.xz")
                with tarfile.open(archive) as f, np.load(io.BytesIO(f.extractfile(f"candidates/{old_name}.npz").read())) as old:
                    for k in res:
                        np.testing.assert_array_equal(res[k],old[k],err_msg=old_name+"/"+k)
                compared.append(old_name)
    table(p/"source_cost.csv",costs)
    table(p/"generation.csv",generation)
    write_json(p/"freeze.json", dict(observer_sha256=digest(observer), equal_to_previous_bank=compared,
        arrays_checked=list(res),generated_cases=len(generation)))

def dock(p, observer):
    costs = []
    for name,c in configs(p).items():
        for seed in CASES[name][2]:
            case = f"{name}_{seed}"
            base = p/"banks"/case
            for arm in ("shape","true","perm"):
                stem=f"{case}_{arm}"
                root=p/"raw"/stem
                env=dict(FD_RIGID_INTERVENTION=base.with_suffix(".input"), FD_RIGID_LOG=root.with_suffix(".inject"),
                         FD_RIGID_COORDS=root.with_suffix(".stages"), FD_RIGID_HISTORY=root.with_suffix(".history"))
                if arm!="shape":
                    env.update(FD_FIELD_BANK=base.with_suffix(".bank"), FD_FIELD_ARM=arm, FD_FIELD_TRACE=root.with_suffix(".csv"))
                pool = invoke(p,costs,stem,observer,c,seed,env=env)
                costs[-1].update(interpolation(pool.with_suffix(".log")))
                src = np.array(next(l for l in root.with_suffix(".inject").read_text().splitlines() if l.startswith("source ")).split()[1:],float)
                np.testing.assert_array_equal(src,np.load(base.with_suffix(".source.npz"))["q"])
                assert (p/"raw"/f"{case}_source.pdbqt").read_bytes()==(p/"raw"/f"{case}_pro.pdbqt").read_bytes()
    table(p/"dock_cost.csv",costs)
    write_json(p/"rank_build.json", dict(sha256=digest(observer)))

if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("action",choices=("prepare","freeze","controls","dock"))
    a.add_argument("work",type=lambda s:Path(s).resolve())
    for k in ("observer","pro","official","counter"):
        a.add_argument("--"+k,type=lambda s:Path(s).resolve())
    args=a.parse_args()
    if args.action=="prepare": prepare(args.work)
    elif args.action=="freeze": freeze(args.work,args.observer)
    elif args.action=="controls": controls(args.work,args.pro,args.official,args.counter)
    else: dock(args.work,args.observer)
