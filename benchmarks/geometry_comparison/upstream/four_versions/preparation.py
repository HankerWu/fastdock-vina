"""Measure input reconstruction stages on the archived reciprocal depositions."""

import argparse
import json
from pathlib import Path
import resource
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor

# Load the two historical preparation modules in their original import context.
HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[2]
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import run as preparation

sys.path.insert(0, str(BENCH / "cubic_grid/human_panel"))
import panel as human_panel


def cpu():
    return sum(
        r.ru_utime + r.ru_stime
        for r in (
            resource.getrusage(resource.RUSAGE_SELF),
            resource.getrusage(resource.RUSAGE_CHILDREN),
        )
    )


def prepare(job):
    import hashlib

    target, pdb, ligand, ref, work = job
    result = dict(target=target, pdb=pdb, ligand=ligand)

    class Timed:
        def __getattr__(self, name):
            function = getattr(preparation, name)
            if name not in ("ligand", "receptor"):
                return function

            def measured(*args):
                a, b = cpu(), time.perf_counter()
                try:
                    return function(*args)
                finally:
                    result[name + "_cpu_s"] = cpu() - a
                    result[name + "_wall_s"] = time.perf_counter() - b

            return measured

    for f in (work / "raw" / (pdb + ".pdb"), work / "raw" / (ligand + ".cif")):
        if not f.is_file():
            raise FileNotFoundError(f)
    start = cpu(), time.perf_counter()
    metadata = human_panel.prepare_entry((target, pdb, ligand), work, Timed())
    result.update(
        status=metadata["status"],
        error=metadata.get("error"),
        cpu_s=cpu() - start[0],
        wall_s=time.perf_counter() - start[1],
    )
    for name in ["receptor", "ligand"]:
        file = work / pdb / (name + ".pdbqt")
        digest = (
            hashlib.sha256(file.read_bytes()).hexdigest() if file.exists() else None
        )
        result[name + "_sha256"] = digest
        result[name + "_matches_frozen"] = digest == ref[name + ".pdbqt"]
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=False)
    with tarfile.open(BENCH / "cubic_grid/human_panel/depositions.tar.xz") as tf:
        tf.extractall(a.work, filter="data")
    frozen = json.loads((HERE / "inputs.json").read_text())["inputs"]
    jobs = [
        (t, pdb, lig, frozen[pdb], a.work)
        for t, pair in human_panel.PANELS.items()
        for pdb, lig in pair.items()
        if pdb in frozen
    ]
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = list(pool.map(prepare, jobs))
    a.out.write_text(json.dumps(rows, indent=2) + "\n")
