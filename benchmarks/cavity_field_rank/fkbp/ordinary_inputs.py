"""Time our two-pass prototype from ordinary inputs to its last reported mode."""
import argparse
import json
import re
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from frozen_assay import ARMS, CASES, HERE, ROOT, RECEIVERS, SEED, cpu_seconds, invoke, rank, run, sha, table
sys.path.insert(0,str(ROOT/"tools"))
from pocket_triangle import Pocket


def from_inputs(root, arm, exe):
    plans, capture = root/"plans", root/"capture"
    plans.mkdir(exist_ok=True)
    capture.mkdir(exist_ok=True)
    collected = []
    for rec,n in CASES:
        p = capture/f"{rec}_{n}"
        collected.append(invoke(root, exe, rec, n, "capture", {
            "FD_SOURCE_STEP":1,"FD_SOURCE_TRAJECTORY":1,
            "FD_SOURCE_XYZ":p.with_suffix(".xyz"),"FD_SOURCE_Q":p.with_suffix(".q")}))
    table(root/"capture_cost.csv",collected)
    for rec in RECEIVERS:
        # The map header carries the exact center used in the original pocket generator.
        with (root/"maps"/rec/"receptor.C.map").open() as f:
            header=[next(f).split() for _ in range(6)]
        assert header[5][0]=="CENTER"
        center=np.array(header[5][1:],float)
        size=np.array(header[4][1:],int)*float(header[3][1])
        pocket=Pocket.from_pdbqt(root/rec/"receptor.pdbqt",center,size)
        for n in range(4,12):
            stem=f"{rec}_{n}"
            lines=(capture/(stem+".xyz")).read_text().splitlines()
            serials=np.array(lines[0][1:].split(),int)
            xyz=np.array(lines[1].split()[3:],float).reshape(-1,3)
            q=np.loadtxt(capture/(stem+".q"))[2:]
            candidates=pocket.place(xyz,SEED,arms=("matched",))
            np.savez_compressed(plans/(stem+".npz"),q=q,xyz=xyz,serials=serials,**candidates)
    rank(root,(arm,))
    run(root,exe,(arm,))


def measure(source, root, arm, exe):
    root.mkdir(parents=True,exist_ok=True)
    (root/"pools/base").mkdir(parents=True,exist_ok=True)
    for name in ("affinity","maps",*RECEIVERS):
        p=root/name
        if not p.exists(): p.symlink_to(source/name,target_is_directory=True)
    # Only ordinary configs enter. The source poses, candidates and historical public pools are absent.
    for rec,n in CASES:
        name=f"pools/base/{rec}_{n}.conf"
        (root/name).write_bytes((source/name).read_bytes())
    start,cpu=time.perf_counter(),cpu_seconds()
    subprocess.run([sys.executable,str(Path(__file__).resolve()),"predict",str(root),arm,str(exe)],check=True)
    elapsed,cput=time.perf_counter()-start,cpu_seconds()-cpu
    cells=[]
    counters={}
    for rec,n in CASES:
        stem=f"{rec}_{n}"
        new,old=root/f"plans/{stem}.npz",source/f"plans/{stem}.npz"
        with np.load(new) as observed,np.load(old) as frozen:
            for key in observed.files:
                np.testing.assert_array_equal(observed[key],frozen[key])
        pool=Path("pools")/arm/(stem+".pdbqt")
        assert (root/pool).read_bytes()==(source/pool).read_bytes()
        assert (root/f"pools/capture/{stem}.pdbqt").read_bytes()==(source/f"pools/base/{stem}.pdbqt").read_bytes()
        for phase in ("capture",arm):
            log=(root/f"pools/{phase}/{stem}.log").read_text()
            queries=re.findall(r"FIELD_COUNT combined=1 derivative=\d+ positions=(\d+)",log)
            assert len(queries)==2
            counters[phase]=counters.get(phase,0)+3*sum(map(int,queries))
        cells.append(dict(receiver=rec,compound=n,pool_sha256=sha(root/pool),
                          source_replay_sha256=sha(root/f"pools/capture/{stem}.pdbqt"),
                          rebuilt_bank_sha256=sha(new)))
    (source/f"ordinary-{arm}.json").write_text(json.dumps(dict(
        arm=arm,wall_s=elapsed,cpu_s=cput,cells=cells,downstream_channels=counters,
        capture_cost=(root/"capture_cost.csv").read_text(),
        rank_cost=(root/"rank_cost.csv").read_text(),
        search_cost=(root/f"runs-{arm}.csv").read_text()),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=("predict","measure"))
    p.add_argument("root",type=Path)
    p.add_argument("arm",choices=ARMS)
    p.add_argument("exe",type=Path)
    p.add_argument("--reference",type=Path)
    a=p.parse_args()
    if a.action=="predict":
        from_inputs(a.root.resolve(),a.arm,a.exe.resolve())
    else:
        measure(a.reference.resolve(),a.root.resolve(),a.arm,a.exe.resolve())
