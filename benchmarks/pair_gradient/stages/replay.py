"""Replay four prespecified single-trajectory cells; RMSD is computed outside Vina."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pose"))
from paired import base


def scores(lines, serial, refs):
    indexed = {}
    for field in lines:
        atom, x, y, z = field.split(",")
        indexed[int(atom)] = [float(x), float(y), float(z)]
    if set(indexed) != set(serial):
        raise ValueError("Incomplete trace-to-reference atom correspondence")
    xyz = np.array([indexed[s] for s in serial])
    return float(np.sqrt(((refs-xyz)**2).sum(axis=2).mean(axis=1)).min())


def counters(log):
    prefixes = ("Native merge:", "Search trajectory:", "Hunt chain:", "Search counters:")
    lines = [s for s in log.splitlines() if s.startswith(prefixes)]
    if len(lines) != 4:
        raise ValueError("Incomplete counter log")
    return lines


def dump(path, records):
    with path.open("w") as f:
        w = csv.DictWriter(f, list(records[0]))
        w.writeheader()
        w.writerows(records)


def main(a):
    query, receptor = "3EML", "3UZA"
    qm, rm = [json.loads((a.inputs/p/"preparation.json").read_text()) for p in (query,receptor)]
    rotation, translation, ca_rmsd, n_ca = base.alignment(qm,rm)
    refs = np.asarray(qm["mapped"]) @ rotation + translation
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    a.output.mkdir(parents=True, exist_ok=True)
    all_rows, checks = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for seed in (75467,98251):
            for arm in ("parent","linear"):
                key = f"{query}_{receptor}_{seed}_12000_{arm}"
                record = a.archived / key
                original = record.with_suffix(".pdbqt").read_bytes()
                archived_counts = counters(record.with_suffix(".log").read_text())
                cmd = ["--receptor",a.inputs/receptor/"receptor.pdbqt","--ligand",a.inputs/query/"ligand.pdbqt",
                       "--seed",str(seed),"--max_evals","12000","--num_modes","9","--energy_range","100",
                       "--cpu","1","--exhaustiveness","1","--search_counters"]
                for axis,value in zip("xyz",rm["center"]):
                    cmd += ["--center_"+axis,str(value),"--size_"+axis,"24"]
                order = ("plain","trace") if (seed==75467)==(arm=="parent") else ("trace","plain")
                elapsed = {}
                for kind in order:
                    binary = a.exports/arm/"build/vina" if kind=="plain" else a.traced/arm/"vina"
                    out = Path(tmp) / (kind+".pdbqt")
                    start=time.perf_counter()
                    env = dict(os.environ)
                    for name in ("FD_CAPTURE_HUNT_151", "FD_REPLAY_SINGLE_HUNT"):
                        env.pop(name, None)
                    if kind == "trace" and arm == "linear" and seed == 75467:
                        env["FD_CAPTURE_HUNT_151"] = str(a.output / "hunt151.conf")
                    run=subprocess.run([binary,*cmd,"--out",out],capture_output=True,check=True,env=env)
                    elapsed[kind]=time.perf_counter()-start
                    if original != out.read_bytes() or counters(run.stdout.decode()) != archived_counts:
                        raise ValueError(f"Output/counters differ: {key} {kind}")
                    if kind == "trace":
                        raw = run.stderr
                with (a.output/(key+".trace.gz")).open("wb") as f:
                    with gzip.GzipFile(filename="",mode="wb",fileobj=f,mtime=0) as gz:
                        gz.write(raw)
                rows=[]
                for line in raw.decode().splitlines():
                    if not line.startswith("FDTRACE "):
                        continue
                    _,stage,index,*coords=line.split()
                    rows.append(dict(arm=arm,seed=seed,stage=stage,index=int(index),
                                     rmsd_A=scores(coords,qm["serial"],refs)))
                output = []
                for block in original.decode().split("MODEL ")[1:]:
                    atom_rows = [s for s in block.splitlines() if s.startswith(("ATOM  ","HETATM"))
                                 and int(s[6:11]) in qm["serial"]]
                    fields = [f"{int(s[6:11])},{s[30:38].strip()},{s[38:46].strip()},{s[46:54].strip()}"
                              for s in atom_rows]
                    output.append(dict(arm=arm,seed=seed,stage="output",index=int(block.split()[0]),
                                       rmsd_A=scores(fields,qm["serial"],refs)))
                with (HERE.parent/"pose/results/modes.csv").open() as f:
                    reference = [float(row["rmsd_A"]) for row in csv.DictReader(f) if row["key"] == key]
                if len(reference) != len(output) or max(abs(row["rmsd_A"] - rmsd)
                        for row, rmsd in zip(output, reference)) > 1e-10:
                    raise ValueError("Output disagrees with independently audited pose panel")
                rows.extend(output)
                counts = {stage:sum(r["stage"]==stage for r in rows) for stage in {r["stage"] for r in rows}}
                expected=dict(re.findall(r"(\w+)=(\d+)", " ".join(archived_counts)))
                for keyname,n in (("hunt_calls",counts["proposal"]),("accepted",counts.get("hunt_accepted",0)),
                                   ("second_calls",counts.get("saved_retained",0)+counts.get("saved_dropped",0)),
                                   ("explicit_inputs",counts["explicit_input"]),("explicit_inputs",counts["explicit_output"])):
                    if int(expected[keyname]) != n:
                        raise ValueError(f"Trace count mismatch for {keyname}")
                if counts.get("hunt_rejected",0)+counts["hunt_accepted"] != counts["proposal"]:
                    raise ValueError("Missing MC decision")
                all_rows.extend(rows)
                checks.append(dict(arm=arm,seed=seed,output_sha256=hashlib.sha256(original).hexdigest(),
                                   native_counters=archived_counts,output_byte_identical=True,
                                   plain_s=elapsed["plain"],trace_s=elapsed["trace"],timing_order=list(order),
                                   events=len(rows)-len(output),modes=len(output)))
    dump(a.output/"rmsd.csv",all_rows)
    result = dict(query=query,receptor=receptor,ca_alignment_rmsd_A=ca_rmsd,common_ca=n_ca,
                  symmetry_mappings=len(refs),heavy_atoms=len(qm["serial"]),pinned_cpu=cpu,checks=checks)
    (a.output/"verification.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","archived","exports","traced","output"):
        p.add_argument(name,type=Path)
    main(p.parse_args())
