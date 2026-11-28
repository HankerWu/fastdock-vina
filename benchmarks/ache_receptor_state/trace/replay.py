"""Passively rerun the two selected global searches. No affinity or query crystal input."""
import argparse
import hashlib
import json
import lzma
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tarfile
import time

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
COUNTERS = ("Search trajectory:", "Hunt chain:", "Pocket start:", "Search counters:")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def clock():
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.perf_counter(), time.process_time()+children.ru_utime+children.ru_stime


def run(binary, work, out):
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=False)
    blind = HERE.parent/"blind/results"
    hashes = json.loads((blind/"inputs.json").read_text())["hashes"]
    rawhup = (BENCH/"water_rank_envelope/transfer/HUP.pdbqt").read_bytes()
    assert sha(rawhup) == hashes["HUP.pdbqt"]
    receipt = {}
    with tarfile.open(BENCH/"receptor_amide_states/transfer/transfer_evidence.tar.xz") as archive:
        for state in ("4EY6", "4EY7"):
            src = blind/"raw"/state/"pro24k_1103503"
            job = work/state
            job.mkdir()
            (job/"HUP.pdbqt").write_bytes(rawhup)
            rec = archive.extractfile("novel/"+state+"/receptor.pdbqt").read()
            assert sha(rec) == hashes[state+"/receptor.pdbqt"]
            (job/"receptor.pdbqt").write_bytes(rec)
            cmd = json.loads((src/"command.json").read_text())
            env = dict(os.environ, FD_BASIN_TRACE=str(job/"events"), FD_MEMBERS=str(job/"members"))
            start = clock()
            with (job/"trace.log").open("w") as log:
                subprocess.run([str(binary), *cmd[1:]], check=True, cwd=job, env=env, stdout=log, stderr=subprocess.STDOUT)
            end = clock()
            counters = lambda p: [l for l in p.read_text().splitlines() if l.startswith(COUNTERS)]
            result = (job/"dry.pdbqt").read_bytes()
            if result != (src/"dry.pdbqt").read_bytes() or not counters(job/"trace.log") or counters(job/"trace.log") != counters(src/"dry.log"):
                raise ValueError(state+" publication or counters diverged")
            compression = clock()
            chunks = []
            for name in ("events","members"):
                packed = lzma.compress((job/name).read_bytes())
                for i in range(0, len(packed), 1500000):
                    part = f"{state}.{name}.xz.{i//1500000:02d}"
                    (out/part).write_bytes(packed[i:i+1500000])
                    chunks.append(part)
            finish = clock()
            shutil.copyfile(job/"trace.log",out/(state+".log"))
            receipt[state] = dict(output_sha256=sha(result), binary_sha256=sha(binary.read_bytes()),
                                  probe_wall_s=end[0]-start[0],probe_cpu_s=end[1]-start[1],
                                  compression_wall_s=finish[0]-compression[0],compression_cpu_s=finish[1]-compression[1],
                                  events_sha256=sha((job/"events").read_bytes()), members_sha256=sha((job/"members").read_bytes()),
                                  chunks=chunks)
    (out/"runs.json").write_text(json.dumps(receipt,indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("binary","work","out"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    run(a.binary.resolve(),a.work.resolve(),a.out.resolve())
