"""Single-pass 1FKG/1FKH prediction and label-isolated, coordinate-based assessment."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parents[1] / "fragment_anchor/fkbp"
sys.path.insert(0, str(PREVIOUS))
import assay_runner as native
from pro_replay import replay

LIVE = ("edge", "feasible_edge", "coverage", "uniform", "uniform_work")
ARMS = ("official", "base", *LIVE, "high", "official_high", "cost_pro")
SEED = 1047349


def budgets(work):
    path = work / "cost_budget.json"
    return {"cost_pro": json.loads(path.read_text())["budget"]} if path.exists() else {}


def set_budget(work):
    with (work / "ordinary_initial.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    cpu = {a:sum(float(r["cpu_s"]) for r in rows if r["arm"]==a) for a in ("base","high")}
    files = sorted((work / "pools/coverage").glob("*.timing"))
    assert len(files)==16 and len(rows)==32
    target = sum(json.loads(p.read_text())["cpu_s"] for p in files)
    units = 1+(target-cpu["base"])/(cpu["high"]-cpu["base"]) if cpu["high"]>cpu["base"] else target/cpu["base"]
    path = work / "cost_budget.json"
    path.write_text(json.dumps(dict(coverage_cpu_s=target,ordinary_cpu_s=cpu,raw_units_24k=units,
                                    budget=24000*max(1,math.ceil(units))),indent=2)+"\n")


def readouts(work, out):
    import numpy as np
    from check_panel import check
    ordinary = []
    for name in ("ordinary_initial.csv","ordinary_cost.csv"):
        with (work/name).open() as stream: ordinary.extend(csv.DictReader(stream))
    native.table(work/"unmodified_pro.csv",ordinary)
    out.mkdir(parents=True,exist_ok=True)
    check(work,out,arms=ARMS,live=LIVE,seed=SEED,gate_arms=(),bias_arm="coverage",
          bias_controls=("base","edge","feasible_edge","uniform","cost_pro"))
    (out / "cost_budget.json").write_bytes((work/"cost_budget.json").read_bytes())
    ranks = []
    for rec in native.RECEPTORS:
        for n in range(4,12):
            stem = f"{rec}_{n}"
            src = work/"pools/uniform"/stem
            sham = work/"pools/uniform_work"/stem
            for suffix in (".source",".q",".npz",".transform",".stages",".xyz",".history",".pdbqt"):
                assert native.sha(src.with_suffix(suffix))==native.sha(sham.with_suffix(suffix))
            for arm in LIVE:
                p = work/"pools"/arm/stem
                with np.load(p.with_suffix(".npz")) as b:
                    i=np.argmin(b["true_score"])
                    stats=json.loads(p.with_suffix(".rank").read_text())
                    row=dict(arm=arm,receiver=rec,compound=n,bank_cardinality=len(b["triangle_id"]),
                             winner_triangle_id=b["triangle_id"][i],coarse_first_ties=(b["true_score"]==b["true_score"][i]).sum(),
                             **stats)
                    ranks.append(row)
    native.table(out/"gates.csv",ranks)
    with (out/"pose_cells.csv").open() as source:
        cells={(row["receiver"],row["compound"],row["arm"]):row for row in csv.DictReader(source)}
    pairs=[]
    for (rec,n,arm), control in cells.items():
        if arm=="coverage": continue
        tested=cells[rec,n,"coverage"]
        pairs.append(dict(receiver=rec,compound=n,control=arm,
                          **{"delta_"+k+"_A":float(tested[k])-float(control[k]) for k in ("top1","best5","oracle","regret")}))
    native.table(out/"paired_pose.csv",pairs)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_subparsers(dest="action",required=True)
    for task in ("prepare","budget","run","ordinary","live","check","pack","unpack"):
        a=s.add_parser(task)
        a.add_argument("work",type=Path)
        if task=="run":
            for key in ("official","engine","counter"): a.add_argument(key,type=Path)
            a.add_argument("--arms",nargs="+",choices=ARMS,required=True)
            a.add_argument("--resume",action="store_true")
        if task=="live":
            a.add_argument("rec",choices=native.RECEPTORS)
            a.add_argument("compound",type=int)
            a.add_argument("arm",choices=LIVE)
            a.add_argument("engine",type=Path)
        if task=="ordinary":
            a.add_argument("exe",type=Path)
            a.add_argument("--cost",action="store_true")
        if task in ("check","pack","unpack"): a.add_argument("dest",type=Path)
    a=p.parse_args()
    work=a.work.resolve()
    if a.action=="prepare": native.prep(work)
    if a.action=="budget": set_budget(work)
    if a.action=="run":
        assert "cost_pro" not in a.arms or budgets(work)
        native.run(work,a.official.resolve(),a.engine.resolve(),a.counter.resolve(),a.arms,arms=ARMS,
                   live=LIVE,seed=SEED,script=Path(__file__).resolve(),budgets=budgets(work),resume=a.resume)
    if a.action=="live":
        sys.path.insert(0,str(HERE.parent/"engine"))
        from coverage_gate import place
        stem=work/"pools"/a.arm/f"{a.rec}_{a.compound}"
        native.live_one(work,a.rec,a.compound,a.arm,a.engine.resolve(),stem,stem.with_suffix(".conf"),
                        generator=place,seed=SEED)
    if a.action=="ordinary":
        limits=budgets(work) if a.cost else {"base":24000,"high":48000}
        assert limits
        replay(work,a.exe.resolve(),seed=SEED,budgets=limits,csv_name="ordinary_cost.csv" if a.cost else "ordinary_initial.csv")
    if a.action=="check": readouts(work,a.dest.resolve())
    if a.action=="pack":
        from package import create
        create(work,a.dest.resolve(),arms=ARMS,replay_files=[work/f for f in ("ordinary_initial.csv","ordinary_cost.csv","cost_budget.json")])
    if a.action=="unpack":
        from package import unpack
        unpack(work,a.dest.resolve())


if __name__=="__main__":
    main()
