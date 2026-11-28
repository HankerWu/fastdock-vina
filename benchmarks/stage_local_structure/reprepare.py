"""Rebuild the registered inputs from deposited bytes, retaining per-deposit process costs."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import hashlib

HERE = Path(__file__).resolve().parent


def main(inputs, work):
    manifest = json.loads((inputs / "manifest.json").read_text())
    work.mkdir(parents=True, exist_ok=False)
    costs = []
    for row in manifest:
        pdb, lig = row["pdb"], row["ligand"]
        dest = work / pdb
        dest.mkdir()
        for name in (pdb + ".pdb", lig + ".cif"):
            shutil.copyfile(inputs / pdb / name, dest / name)
        argv = [sys.executable, __file__, "--one", str(work), pdb, lig]
        start = time.monotonic()
        with (dest / "stdout.log").open("w") as log:
            proc = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        hashes = {
            n: hashlib.sha256((dest / n).read_bytes()).hexdigest()
            for n in ("metadata.json", "ligand.pdbqt", "receptor.pdbqt")
            if (dest / n).is_file()
        }
        assert (proc.returncode == 0) == (row["status"] == "ready")
        assert all(
            hashes.get(n) == row["sha256"].get(n)
            for n in ("metadata.json", "ligand.pdbqt", "receptor.pdbqt")
        )
        costs.append(
            dict(
                pdb=pdb,
                argv=argv,
                exit=proc.returncode,
                sha256=hashes,
                wall_s=time.monotonic() - start,
                cpu_s=usage.ru_utime + usage.ru_stime,
                peak_rss_kib=usage.ru_maxrss,
            )
        )
    print(json.dumps(costs, indent=2))


if __name__ == "__main__":
    if sys.argv[1] == "--one":
        sys.path.insert(0, str(HERE.parent / "geometry_metric"))
        import panel

        panel.prepare(Path(sys.argv[2]), sys.argv[3], sys.argv[4], "ESR1")
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("inputs", type=Path)
        parser.add_argument("work", type=Path)
        main(**vars(parser.parse_args()))
