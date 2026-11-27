"""Replay official Vina 1.2.7 and double-budget native controls on the six fixed inputs."""

import argparse
import re
import shlex
import sys
from pathlib import Path

import numpy as np

from transfer import HERE, save
from study import digest, read_csv, timed_process, write_csv
from verify import geometry, point

sys.path.insert(0, str(HERE.parent / "control"))
from benchmark import OFFICIAL


def run(work, output, binary, official):
    import json
    import os

    if "FASTDOCK_RESEARCH_CHPI" in os.environ:
        raise ValueError("Run controls without a research CH-pi field")
    if (
        digest(binary) != json.loads((output / "binary.json").read_text())["sha256"]
        or digest(official) != OFFICIAL
    ):
        raise ValueError("Unexpected binary for the native/official comparison")
    orig = [row for row in read_csv(output / "costs.csv") if row["search"] == "native"]
    records, modes = [], []
    for i, row in enumerate(orig):
        order = [("official", official, 24000), ("native_double", binary, 48000)]
        refs, serials = geometry(work, "GLA", row["query"], row["receptor"])
        for name, engine, cap in order if i % 2 == 0 else reversed(order):
            key = dict(
                query=row["query"],
                receptor=row["receptor"],
                seed=row["seed"],
                arm=name,
                cap=cap,
            )
            args = shlex.split(row["args"])
            args[args.index("--max_evals") + 1] = str(cap)
            if name == "official":
                for flag in ("--pose_selection", "--refine_candidates"):
                    pos = args.index(flag)
                    del args[pos : pos + 2]
                args.remove("--search_counters")
            stem = f"{row['receptor']}_{row['seed']}_{name}"
            pose = output / "pools" / (stem + ".pdbqt")
            log = output / "pools" / (stem + ".log")
            cmd = list(args)
            cmd[cmd.index("--out") + 1] = str(pose)
            timed = timed_process([str(engine), *cmd], work, log)
            errors = []
            for block in (b for b in pose.read_text().split("ENDMDL") if b.strip()):
                xyz = {
                    int(line[6:11]): point(line)
                    for line in block.splitlines()
                    if line.startswith(("ATOM  ", "HETATM"))
                }
                xyz = np.array([xyz[s] for s in serials])
                err = float(np.sqrt(((refs - xyz) ** 2).sum(-1).mean(-1)).min())
                errors.append(err)
                modes.append(dict(**key, rank=len(errors), rmsd=err))
            trace = log.read_text()
            counters = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", trace)
            explicit = re.findall(r"explicit_evals=(\d+)", trace)
            if name == "native_double" and (len(counters), len(explicit)) != (4, 1):
                raise ValueError(f"Missing search counters in {log}")
            records.append(
                dict(
                    **key,
                    top1=errors[0],
                    best5=min(errors[:5]),
                    pool_min=min(errors),
                    regret=errors[0] - min(errors),
                    count=len(errors),
                    **timed,
                    search_evals=sum(int(a) + int(b) for a, b in counters)
                    if counters
                    else "",
                    explicit_evals=sum(map(int, explicit)) if explicit else "",
                    sha256=digest(pose),
                    args=shlex.join(args),
                )
            )
    write_csv(output / "controls.csv", records)
    write_csv(output / "control_modes.csv", modes)
    save(
        output / "control_binaries.json",
        {"official_1.2.7": digest(official), "native_double": digest(binary)},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "output", "binary", "official"):
        parser.add_argument(name, type=lambda p: Path(p).resolve())
    args = parser.parse_args()
    run(args.work, args.output, args.binary, args.official)
