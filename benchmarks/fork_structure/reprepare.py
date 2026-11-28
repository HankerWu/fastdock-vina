"""Rebuild from frozen unprocessed deposit and CCD bytes; count per-process memory and CPU."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def main(fixtures, work):
    work.mkdir(parents=True, exist_ok=False)
    rows = []
    for record in json.loads((fixtures / "manifest.json").read_text()):
        pdb, ligand = record["pdb"], record["ligand"]
        folder = work / pdb
        folder.mkdir()
        for name in (pdb + ".pdb", ligand + ".cif"):
            shutil.copyfile(fixtures / pdb / name, folder / name)
        cmd = [sys.executable, __file__, "--one", str(work), pdb, ligand]
        start = time.monotonic()
        with (folder / "stdout.log").open("w") as log:
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        digests = {
            name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
            for name in ("metadata.json", "ligand.pdbqt", "receptor.pdbqt")
            if (folder / name).is_file()
        }
        assert (proc.returncode == 0) == (record["status"] == "ready")
        assert all(
            digests.get(name) == record["sha256"].get(name)
            for name in ("metadata.json", "ligand.pdbqt", "receptor.pdbqt")
        )
        rows.append(
            dict(
                pdb=pdb,
                command=cmd,
                exit=proc.returncode,
                sha256=digests,
                wall_s=time.monotonic() - start,
                cpu_s=usage.ru_utime + usage.ru_stime,
                peak_rss_kib=usage.ru_maxrss,
            )
        )
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    if sys.argv[1] == "--one":
        sys.path.insert(0, str(HERE.parent / "geometry_metric"))
        import panel

        panel.prepare(Path(sys.argv[2]), sys.argv[3], sys.argv[4], "MAPK14")
    else:
        p = argparse.ArgumentParser(description=__doc__)
        p.add_argument("fixtures", type=Path)
        p.add_argument("work", type=Path)
        main(**vars(p.parse_args()))
