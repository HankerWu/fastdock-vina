"""Link MC attempts, trajectory container membership, merged minima and published poses."""
import argparse
import hashlib
import json
import lzma
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NAMES = {1:"initial",4:"hunt",5:"accepted",6:"second",8:"merged",9:"explicit_input",
         10:"explicit_output",11:"proposal",12:"trajectory_final"}


def recovered(data, case, name, expected):
    pieces = sorted(data.glob(f"{case}.{name}.xz.*"))
    raw = lzma.decompress(b"".join(p.read_bytes() for p in pieces))
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("Damaged "+name+" stream")
    return raw


def read_stream(raw):
    n = int.from_bytes(raw[:4], "little")
    serial = np.frombuffer(raw, dtype="<u4", count=n, offset=4).tolist()
    dtype = np.dtype([(k,"<u4") for k in ("traj","step","phase","event","index")]
                     + [("energy","<f8"),("xyz","<f8",(n,3))])
    return serial, np.frombuffer(raw,dtype=dtype,offset=4*(n+1))


def released(directory, serials):
    for block in (directory/"dry.pdbqt").read_text().split("ENDMDL")[:-1]:
        positions = {int(l[6:11]): [float(l[k:k+8]) for k in (30,38,46)] for l in block.splitlines()
                     if l.startswith(("ATOM  ","HETATM"))}
        yield [positions[s] for s in serials]


def append_stage(rows, case, stage, source):
    for traj in (0,1,2,3,"all"):
        group = source if traj == "all" else source[source.traj == traj]
        if group.empty:
            continue
        best=group.loc[group.rmsd_A.idxmin()]
        rows.append(dict(receiver=case, stage=stage, trajectory=traj, n=len(group),
                         near_le2=int((group.rmsd_A<=2.).sum()), min_A=best.rmsd_A,
                         best_traj=best.traj, best_step=best.step))


def analyze(results, metadata, dest):
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    runs = json.loads((results/"runs.json").read_text())
    truth = json.loads(metadata.read_text())
    stage_rows, life_rows, links, cells, identities, steps = [], [], [], [], [], []
    for state, manifest in runs.items():
        events = recovered(results,state,"events",manifest["events_sha256"])
        serials, tr = read_stream(events)
        rows = np.arange(len(tr))
        ref = np.array(truth["targets"][state]["truth"])[[truth["serials"].index(s) for s in serials]]
        err = np.sqrt(np.sum((tr["xyz"] - ref)**2,axis=2).mean(axis=1))
        df = pd.DataFrame(dict(row=rows,stage=[NAMES[x] for x in tr["event"]],traj=tr["traj"],
                               step=tr["step"],energy=tr["energy"],rmsd_A=err))
        for stage, group in df.groupby("stage"):
            append_stage(stage_rows,state,stage,group)
        src = HERE.parent/"blind/results/raw"/state/"pro24k_1103503"
        native = (src/"dry.log").read_text()
        for stage, pattern in (("proposal","hunt_calls"),("hunt","hunt_calls"),("accepted","accepted"),("second","second_calls")):
            expected = list(map(int,re.findall(rf"\b{pattern}=(\d+)",native)))
            got = df[df.stage==stage].groupby("traj").size().tolist()
            assert got == expected,(stage,got,expected)
        second = df[df.stage=="second"].set_index(["traj","step"])
        membership = recovered(results,state,"members",manifest["members_sha256"]).decode().splitlines()
        lasts = {t:set() for t in (1,2,3)}
        lifetimes = {}
        for line in membership:
            traj, step, *current = map(int, line.split())
            previous = lasts[traj]
            active = set(current)
            assert len(active)==len(current) and (traj,step) in second.index
            added, removed = active-previous, previous-active
            assert not added or added=={step}, (traj,step,added)
            lifetimes[traj,step] = dict(admitted=step in added, evicted_at=None)
            for gone in removed:
                lifetimes[traj,gone]["evicted_at"] = step
            lasts[traj] = active
        assert set(lifetimes)==set(second.index)
        table = df[df.stage.isin(["proposal", "hunt", "accepted", "second"])].pivot(
            index=["traj", "step"], columns="stage", values="rmsd_A").reset_index()
        table.columns.name = None
        table["receiver"] = state
        table["admitted"] = [lifetimes.get((row.traj, row.step), {}).get("admitted", False)
                             for row in table.itertuples()]
        steps.append(table)
        kept = set((t,s) for t, v in lasts.items() for s in v)
        final = df[df.stage=="trajectory_final"]
        assert kept == set(zip(final.traj, final.step))
        for (t,s), memo in lifetimes.items():
            row = second.loc[(t,s)]
            life_rows.append(dict(receiver=state,traj=t,step=s,rmsd_A=row.rmsd_A,energy=row.energy,**memo,
                                  in_final=(t,s) in kept))
        inserted = second.loc[[k for k,v in lifetimes.items() if v["admitted"]]].reset_index()
        append_stage(stage_rows,state,"admitted",inserted)
        evicted = second.loc[[k for k,v in lifetimes.items() if v["evicted_at"] is not None]].reset_index()
        append_stage(stage_rows,state,"evicted",evicted)
        # For each merged or explicit-input record recover the trajectory-level source.
        parents = {}
        for stage in ("merged","explicit_input"):
            for x in df[df.stage==stage].itertuples():
                matches = final[(final.energy==x.energy) &
                                np.array([np.array_equal(tr["xyz"][i],tr["xyz"][x.row]) for i in final.row])]
                if len(matches)!=1:
                    raise ValueError(f"Ambiguous source for {state} {stage} {x.step}: {len(matches)}")
                parent=matches.iloc[0]
                links.append(dict(receiver=state,stage=stage,ordinal=x.step,traj=parent.traj,step=parent.step,
                                  rmsd_A=x.rmsd_A,rounding_max_A=0))
                if stage=="explicit_input":
                    parents[x.step]=(parent.traj,parent.step)
        pool=list(released(src,serials))
        outputs=df[df.stage=="explicit_output"]
        recovered_steps=set()
        public_rows=[]
        for rank, xyz in enumerate(pool,1):
            delta = np.max(np.abs(tr["xyz"][outputs.row]-xyz),axis=(1,2))
            ids=np.flatnonzero(delta<=.000501)
            if len(ids)!=1:
                raise ValueError(f"{state} published mode {rank} has {len(ids)} origins")
            out=outputs.iloc[ids[0]]
            assert out.step not in recovered_steps
            recovered_steps.add(out.step)
            t,s=parents[out.step]
            rmsd=np.sqrt(np.sum((xyz-ref)**2,axis=1).mean())
            links.append(dict(receiver=state,stage="release",ordinal=rank,traj=t,step=s,rmsd_A=rmsd,
                              rounding_max_A=delta[ids[0]]))
            public_rows.append(dict(traj=t,step=s,rmsd_A=rmsd))
            if rank==1:
                for stage in ("proposal","hunt","accepted","second","trajectory_final"):
                    r=df[(df.traj==t)&(df.step==s)&(df.stage==stage)]
                    if len(r):
                        row=r.iloc[0]
                        identities.append(dict(receiver=state,stage=stage,traj=t,step=s,rmsd_A=row.rmsd_A))
                identities.extend(dict(receiver=state,stage=stage,traj=t,step=s,rmsd_A=e) for stage,e in
                                  (("explicit_input",second.loc[(t,s)].rmsd_A),("explicit_output",out.rmsd_A),("release",rmsd)))
        append_stage(stage_rows,state,"release",pd.DataFrame(public_rows))
        cells.append(dict(state=state,mc_steps=(df.stage=="hunt").sum(),accepted=(df.stage=="accepted").sum(),
                          second=len(second),saved=len(inserted),evicted=len(evicted),trajectory_survivors=len(final),
                          merged=(df.stage=="merged").sum(),explicit=len(outputs),released=len(pool)))
    dest.mkdir(parents=True,exist_ok=True)
    for name, data in (("stages",stage_rows),("trajectory_lifetimes",life_rows),("origins",links),
                       ("cells",cells),("winner_path",identities)):
        pd.DataFrame(data).to_csv(dest/(name+".csv"),index=False)
    pd.concat(steps, ignore_index=True).to_csv(dest/"steps.csv.gz", index=False,
                                              compression={"method":"gzip","mtime":0})
    (dest/"evaluation.json").write_text(json.dumps(dict(wall_s=time.perf_counter()-start_wall,
                                                        cpu_s=time.process_time()-start_cpu),indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("results","metadata","out"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    analyze(a.results,a.metadata,a.out)
