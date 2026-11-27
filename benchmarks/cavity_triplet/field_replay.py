"""Exact query counters on the same production coordinates; separate diagnostic cost."""
import argparse
import csv
import os
import re
import resource
import subprocess
import time
from pathlib import Path

from pose_runs import digest, systems, table


def measure(stem, cmd, directory, prediction, env, output, upstream=False):
    clean = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
    cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    log = subprocess.run(list(map(str,cmd)),cwd=directory,env=clean|env,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True).stdout.decode()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    if output.read_bytes() != prediction.read_bytes():
        raise ValueError(f"field counter changed pool: {stem}")
    counts = {(int(combined),int(deriv)):int(n) for combined,deriv,n in re.findall(
        r"FIELD_COUNT combined=(\d+) derivative=(\d+) positions=(\d+)",log)}
    if upstream:
        channels = dict((int(d),int(n)) for d,n in re.findall(r"GRID_CHANNEL derivative=(\d+) n=(\d+)",log))
        assert len(channels)==2
        combined = stem.startswith("fkbp_")
        divisor = 3 if combined else 1
        if any(n % divisor for n in channels.values()):
            raise ValueError("Unexpected upstream AD4 channel count")
        counts={(i,d):(channels[d]//divisor if i == combined else 0) for i in (0,1) for d in (0,1)}
    assert len(counts)==4
    return dict(cell=stem, pool_sha256=digest(output),
                vina_scalar=counts[0,0],vina_with_gradient=counts[0,1],
                ad4_scalar=counts[1,0],ad4_with_gradient=counts[1,1],
                channel_interpolations=sum((1+2*ad4)*v for (ad4,_),v in counts.items()),
                diagnostic_wall_s=time.perf_counter()-start,
                diagnostic_cpu_s=after.ru_stime+after.ru_utime-cpu.ru_stime-cpu.ru_utime)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("panel","fkbp","binary","work","report"):
        p.add_argument(key,type=lambda s:Path(s).resolve())
    p.add_argument("--upstream",action="store_true")
    a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    rows=[]
    cases=systems(a.panel)
    out=a.work/"replayed.pdbqt"
    for csvfile in ("control_cost.csv","branch_cost.csv"):
        for row in csv.DictReader((a.panel/csvfile).open()):
            stem=row["stem"]
            if ("official" in stem) != a.upstream:
                continue
            name,num,*_=stem.split("_")
            cell=cases[name]
            spec=a.panel/"raw"/(stem+".input")
            env=dict(FD_RIGID_INTERVENTION=str(spec),FD_RIGID_LOG=str(a.work/"injected.log")) if spec.exists() else {}
            cap=50000 if stem.endswith("_double") else 25000
            cmd=[a.binary,"--ligand",cell["ligand"],"--receptor",cell["receptor"]]
            if not a.upstream:
                cmd.extend(["--ligand_bonds",cell["bonds"],"--pose_selection","late","--search_counters"])
            for axis,center in zip("xyz",cell["center"],strict=True):
                cmd.extend([f"--center_{axis}",center,f"--size_{axis}",cell["size"]])
            cmd.extend(["--seed",num,"--exhaustiveness",2,"--cpu",1,"--max_evals",cap,"--num_modes",9,
                        "--energy_range",100,"--out",out])
            rows.append(measure(stem,cmd,a.panel,a.panel/"raw"/(stem+".pdbqt"),env,out,a.upstream))
    for arm in (("official",) if a.upstream else ("base","high","match","shuffle")):
        for original in sorted((a.fkbp/"pools"/arm).glob("*.conf")):
            stem=original.stem
            conf=a.work/"config.txt"
            conf.write_text(re.sub(r"^out\s*=.*$",f"out = {out}",original.read_text(),flags=re.M))
            env={}
            if arm in ("match","shuffle"):
                env=dict(FD_RIGID_INTERVENTION=str(a.fkbp/"plans"/f"{stem}_{arm}.input"),
                         FD_RIGID_LOG=str(a.work/"injected.log"))
            rows.append(measure("fkbp_"+stem+"_"+arm,[a.binary,"--config",conf]+([] if a.upstream else ["--search_counters"]),
                                a.fkbp,original.with_suffix(".pdbqt"),env,out,a.upstream))
    table(a.report,rows)
