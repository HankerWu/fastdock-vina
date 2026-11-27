"""Apply the two existing hunt optimizers once to the same saved internal coordinates."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np
from replay import base, scores


def run(inputs, traced, output):
    q, r = "3EML", "3UZA"
    qm, rm = (json.loads((inputs/s/"preparation.json").read_text()) for s in (q,r))
    rot, trans, _, _ = base.alignment(qm, rm)
    refs = np.asarray(qm["mapped"]) @ rot + trans
    conf = output/"hunt151.conf"
    trace = output/"3EML_3UZA_75467_12000_linear.trace.gz"
    with gzip.open(trace, "rt") as f:
        original = {s.split()[1]:s.split()[3:] for s in f if s.startswith(("FDTRACE proposal 151 ",
                                                                            "FDTRACE hunt_accepted 151 "))}
    result = {"input_sha256":hashlib.sha256(conf.read_bytes()).hexdigest(), "arms":{}}
    for arm in ("parent","linear"):
        executable = traced/arm/"vina"
        command = [executable,"--receptor",inputs/r/"receptor.pdbqt","--ligand",inputs/q/"ligand.pdbqt",
                   "--seed","75467","--cpu","1","--exhaustiveness","1","--max_evals","12000","--num_modes","9",
                   "--out",os.devnull]
        for axis,center in zip("xyz",rm["center"]):
            command += ["--center_"+axis,str(center),"--size_"+axis,"24"]
        env = dict(os.environ, FD_REPLAY_SINGLE_HUNT=str(conf))
        env.pop("FD_CAPTURE_HUNT_151",None)
        p = subprocess.run(command,env=env,capture_output=True,text=True,check=True)
        lines = [s for s in p.stderr.splitlines() if s.startswith(("FDTRACE fixed_", "FIXED_HUNT "))]
        (output/(arm+"_fixed.txt")).write_text("\n".join(lines)+"\n")
        coords = {s.split()[1]:s.split()[3:] for s in lines if s.startswith("FDTRACE")}
        stats = next(s.split()[1:] for s in lines if s.startswith("FIXED_HUNT "))
        if coords["fixed_input"] != original["proposal"]:
            raise ValueError("Input restoration changed coordinates")
        if arm == "linear" and coords["fixed_output"] != original["hunt_accepted"]:
            raise ValueError("Linear isolated endpoint differs from predictive trajectory")
        row = dict(local_objective_gradient_calls=int(stats[0]), diagnostic_initial_calls=1,
                   initial_energy=float(stats[1]),final_energy=float(stats[2]),
                   initial_A=scores(coords["fixed_input"],qm["serial"],refs),
                   final_A=scores(coords["fixed_output"],qm["serial"],refs))
        row["source_sha256"]={name:hashlib.sha256((traced/arm/"src/lib"/name).read_bytes()).hexdigest()
                               for name in ("precalculate.h","bfgs.h","quasi_newton.cpp","monte_carlo.cpp")}
        result["arms"][arm]=row
    if result["arms"]["parent"]["initial_energy"] != result["arms"]["linear"]["initial_energy"]:
        raise ValueError("The two controls evaluate different starting objectives")
    (output/"fixed_hunt.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","traced","output"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    run(a.inputs,a.traced,a.output)
