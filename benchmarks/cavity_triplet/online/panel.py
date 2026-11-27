"""Fresh single-pass engines, geometry arms and unmodified CLI control; compare archived pools."""

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pose_runs import counts, digest, systems, table, SEEDS


def panel(w, engine, out):
    out.mkdir(parents=True,exist_ok=True)
    rows = []
    cfg = systems(w)
    arms = ["pro","matched","sham"]
    for i,(drug,c) in enumerate(cfg.items()):
        for j,seed in enumerate(SEEDS):
            offset = (i+j)%len(arms)
            for arm in arms[offset:]+arms[:offset]:
                stem = f"{drug}_{seed}_{arm}"
                pool = out/(stem+".pdbqt")
                common = ["--ligand",c["ligand"],"--receptor",c["receptor"]]
                for a,v in zip("xyz",c["center"],strict=True):
                    common.extend(["--center_"+a,v,"--size_"+a,c["size"]])
                common += ["--seed",seed,"--exhaustiveness",2,"--cpu",1,"--max_evals",25000,
                           "--num_modes",9,"--energy_range",100,"--out",pool,
                           "--ligand_bonds",c["bonds"],"--pose_selection","late","--search_counters"]
                cmd = [engine] if arm=="pro" else [
                    sys.executable,Path(__file__).with_name("dock.py"),"--engine",engine,"--arm",arm,"--"]
                log = out/(stem+".log")
                env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
                with log.open("w") as stream:
                    before = resource.getrusage(resource.RUSAGE_CHILDREN)
                    t0=time.perf_counter()
                    subprocess.run(list(map(str,cmd+common)),env=env|{"OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1"},
                                   check=True,stdout=stream,stderr=subprocess.STDOUT)
                    wall=time.perf_counter()-t0
                    after=resource.getrusage(resource.RUSAGE_CHILDREN)
                cpu = after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
                expected = w/"raw"/(stem+".pdbqt" if arm=="pro" else f"{drug}_{seed}_t1_s1_{arm}.pdbqt")
                assert pool.read_bytes()==expected.read_bytes()
                counter=counts(log)
                assert counter==counts(expected.with_suffix(".log"))
                note = [l.split("Live cavity placement ")[1] for l in log.read_text().splitlines() if "Live cavity placement " in l]
                assert len(note)==(arm!="pro")
                field=json.loads(note[0]) if note else {}
                rows.append(dict(drug=drug,seed=seed,arm=arm,wall_s=wall,cpu_s=cpu,sha256=digest(pool),
                                 identical_pool=1, atom_queries=field.get("atom_queries",0),
                                 site_nearest=field.get("site_nearest_queries",0),
                                 site_count=field.get("site_count_queries",0),
                                 selected=field.get("selected"),**counter))
    table(out/"jobs.csv",rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("work","engine","out"):
        p.add_argument(k,type=lambda x:Path(x).resolve())
    a=p.parse_args()
    panel(a.work,a.engine,a.out)
