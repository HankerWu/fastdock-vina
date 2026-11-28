"""Replay chemistry from deposited bytes; time each molecule including external preparation tools."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
if __name__ == "__main__":
    if sys.argv[1] == "--one":
        sys.path.insert(0, str(HERE.parent / "geometry_metric"))
        import panel

        panel.prepare(Path(sys.argv[2]), sys.argv[3], sys.argv[4], "BRD4")
        sys.exit(0)
    inputs, work = map(Path, sys.argv[1:])
    work.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((inputs / "manifest.json").read_text())
    receipts = []
    for case in manifest:
        pdb, ligand = case["pdb"], case["ligand"]
        folder = work / pdb
        folder.mkdir()
        for name in (pdb + ".pdb", ligand + ".cif"):
            shutil.copyfile(inputs / pdb / name, folder / name)
        start = time.monotonic()
        command = [sys.executable, __file__, "--one", str(work), pdb, ligand]
        with (folder / "stdout.log").open("w") as log:
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        matched = {}
        for name in ("ligand.pdbqt", "receptor.pdbqt", "metadata.json"):
            actual = folder / name
            expected = case["sha256"].get(name)
            matched[name] = (
                hashlib.sha256(actual.read_bytes()).hexdigest() == expected
                if actual.is_file()
                else expected is None
            )
        assert (proc.returncode == 0) == (case["status"] == "ready") and all(
            matched.values()
        )
        receipts.append(
            dict(
                pdb=pdb,
                argv=command,
                status=case["status"],
                exit=proc.returncode,
                offline_wall_s=time.monotonic() - start,
                cpu_s=usage.ru_utime + usage.ru_stime,
                peak_rss_kib=usage.ru_maxrss,
                input_outputs_match=matched,
            )
        )
    print(json.dumps(receipts, indent=2))
