"""Score the first changed proposal and its two local optimizations on the same RNG state."""

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from panel import CASES, SEEDS, holdout


def main(work, always, adaptive):
    dest = work / "replay"
    dest.mkdir(exist_ok=True)
    meta = {p: json.loads((work / p / "preparation.json").read_text()) for p in CASES}
    rows = []
    decisions = []
    for q in ("1FBY", "2LBD"):
        ref = np.array(meta[q]["mapped"])
        for seed in SEEDS:
            for arm, binary in (("always", always), ("adaptive", adaptive)):
                key = f"{q}_{q}_{seed}_12000_{arm}"
                output = dest / (key + ".pdbqt")
                command = holdout.args_for(work, meta, q, q, seed, 12000, output)
                result = subprocess.run(
                    [str(binary.resolve()), *command, "--search_counters"],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                if output.read_bytes() != (work / "runs" / output.name).read_bytes():
                    raise ValueError("observer changed public output: " + key)
                log = result.stdout + result.stderr
                production = (work / "runs" / (key + ".log")).read_text()
                if holdout.paired.counts(log) != holdout.paired.counts(production):
                    raise ValueError("observer changed optimizer work: " + key)
                (dest / (key + ".log")).write_text(log)
                counter = re.search(r"Silent torsions:.*?redirected=(\d+)", log)
                if counter is None:
                    raise ValueError("missing torsion counter: " + key)
                if int(counter[1]) == 0:
                    if "Silent replay step=" in log:
                        raise ValueError("trace without redirection: " + key)
                    continue
                starts = re.findall(
                    r"Silent replay step=(\d+) stage=(\w+) source_rms=(\S+) max_displacement=(\S+)",
                    log,
                )
                if [s[1] for s in starts] != [
                    "source",
                    "native_proposal",
                    "redirected_proposal",
                    "native_hunt",
                    "redirected_hunt",
                ]:
                    raise ValueError(key)
                blocks = re.findall(r"REMARK STAGE [^\n]+(.*?)ENDMDL", log, re.S)
                if len(blocks) != len(starts):
                    raise ValueError((key, len(blocks)))
                for (step, stage, rms, maximum), block in zip(starts, blocks):
                    heavy = holdout.paired.atoms(block.splitlines())
                    xyz = [heavy[i][0] for i in meta[q]["serial"]]
                    rmsd = np.sqrt(np.square(ref - xyz).sum(axis=2).mean(axis=1)).min()
                    rows.append(
                        dict(
                            key=key,
                            step=step,
                            stage=stage,
                            source_rms_A=rms,
                            max_displacement_A=maximum,
                            crystal_rmsd_A=rmsd,
                        )
                    )
                branches = re.findall(
                    r"Silent branch arm=(\d+) initial=(\S+) final=(\S+) reference=(\S+) accepted=(\d+) evals=(\d+)",
                    log,
                )
                if len(branches) != 2:
                    raise ValueError(key)
                for branch, initial, final, reference, accepted, evals in branches:
                    decisions.append(
                        dict(
                            key=key,
                            branch=("native", "redirected")[int(branch)],
                            initial=initial,
                            final=final,
                            reference=reference,
                            accepted=accepted,
                            hunt_evals=evals,
                        )
                    )
    holdout.paired.write(work / "first_proposal.csv", rows)
    holdout.paired.write(work / "first_hunt.csv", decisions)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "always", "adaptive"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    main(a.work.resolve(), a.always, a.adaptive)
