"""Time independent cold-process two-pass predictions, verify their streams and count duplicated work."""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

def worker(p, work, observer, drug, seed, arm):
    from run import configs, csvrows, freeze, interpolation, invoke, write_json
    case = f"{drug}_{seed}"
    (p/"raw").mkdir()
    (p/"banks").mkdir()
    (p/"inputs").symlink_to(work/"inputs",target_is_directory=True)
    (p/"inputs.json").symlink_to(work/"inputs.json")
    freeze(p,observer,only=case)
    c=configs(p)[drug]
    env=dict(FD_RIGID_INTERVENTION=p/"banks"/(case+".input"),FD_RIGID_LOG=p/"raw/inject.log")
    if arm!="shape":
        env.update(FD_FIELD_BANK=p/"banks"/(case+".bank"),FD_FIELD_ARM=arm)
    stats=[]
    pool=invoke(p,stats,"out",observer,c,seed,env=env)
    stem=f"{case}_{arm}"
    assert pool.read_bytes()==(work/"raw"/(stem+".pdbqt")).read_bytes()
    assert (p/"banks"/(case+".bank")).read_bytes()==(work/"banks"/(case+".bank")).read_bytes()
    src=csvrows(p/"source_cost.csv")[0]
    scan=csvrows(p/"generation.csv")[0]
    source_fields=interpolation(p/"raw"/(case+"_source.log"))
    search_fields=interpolation(pool.with_suffix(".log"))
    write_json(p/"measured.json",dict(stem=stem,pool_sha256=stats[0]["pool_sha256"],
        source_wall_s=src["wall_s"], bank_wall_s=scan["wall_s"], search_wall_s=stats[0]["wall_s"],
        source_cpu_s=src["cpu_s"], bank_cpu_s=scan["cpu_s"], search_cpu_s=stats[0]["cpu_s"],
        search_scalar=search_fields["scalar_fields"],search_gradient=search_fields["gradient_fields"],
        source_scalar=source_fields["scalar_fields"],source_gradient=source_fields["gradient_fields"],
        geometry_nearest=int(scan["shape_queries"]),site_nearest=int(scan["site_nearest"]),
        site_burial=int(scan["site_burial"]),**{k:int(stats[0][k])+int(src[k]) for k in
            ("hunt_evals","second_evals","explicit_evals")}))

def run(p,observer,production):
    from run import CASES,configs,csvrows,digest,invoke,table
    elapsed,baseline=[],[]
    # Actual production entry point without diagnostic grid counters or coordinate persistence.
    for name,c in configs(p).items():
        for seed in CASES[name][2]:
            stem=f"{name}_{seed}_pro"
            pool=invoke(p,baseline,stem+"-clock",production,c,seed)
            assert pool.read_bytes()==(p/"raw"/(stem+".pdbqt")).read_bytes()
            for arm in ("shape","true","perm"):
                with tempfile.TemporaryDirectory(dir=p) as temp:
                    folder=Path(temp)
                    args=[sys.executable, __file__, str(p),str(observer),"--one",
                          str(folder),name,str(seed),arm]
                    start=time.perf_counter()
                    before=resource.getrusage(resource.RUSAGE_CHILDREN)
                    subprocess.run(args,check=True)
                    end=time.perf_counter()
                    after=resource.getrusage(resource.RUSAGE_CHILDREN)
                    row=json.loads((folder/"measured.json").read_text())
                    row.update(end_to_end_wall_s=end-start,
                         end_to_end_cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
                    elapsed.append(row)
    table(p/"cold_timing.csv",elapsed)
    table(p/"production_timing.csv",baseline)
    (p/"production.sha256").write_text(digest(production)+"\n")

if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("work",type=Path);a.add_argument("observer",type=Path)
    a.add_argument("--production",type=Path);a.add_argument("--one",nargs=4)
    args=a.parse_args()
    if args.one:
        folder,drug,seed,arm=args.one
        worker(Path(folder),args.work.resolve(),args.observer.resolve(),drug,int(seed),arm)
    else:
        run(args.work.resolve(),args.observer.resolve(),args.production.resolve())
