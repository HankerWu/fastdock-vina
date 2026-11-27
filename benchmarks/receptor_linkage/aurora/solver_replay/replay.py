"""Re-run the two affected Pro arms of the frozen AURKA experiment."""

import argparse
import csv
import hashlib
import json
import resource
import shutil
import subprocess
import time
from pathlib import Path


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def replay(frozen, run, binary):
    executable_hash = digest(binary)
    with (frozen / "runs.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 72
    (run / "pools").mkdir(parents=True, exist_ok=True)
    shutil.copytree(frozen / "prepared", run / "prepared", dirs_exist_ok=True)
    with (run / "runs.csv").open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=rows[0].keys())
        wr.writeheader()
        for row in rows:
            method, rec, drug, seed = [
                row[k] for k in ("method", "receptor", "compound", "seed")
            ]
            stem = f"{method}_{rec}_{drug}_{seed}"
            name = Path("pools") / (stem + ".pdbqt")
            log = name.with_suffix(".log")
            assert digest(frozen / name) == row["pool_sha256"]
            if method == "official":
                for rel in (name, log):
                    shutil.copyfile(frozen / rel, run / rel)
            else:
                cap = 25000 if method == "pro" else 50000
                assert int(row["cap"]) == cap
                args = [
                    "--receptor",
                    f"prepared/{rec}.pdbqt",
                    "--ligand",
                    f"prepared/{drug}.pdbqt",
                    "--cpu",
                    "1",
                    "--seed",
                    seed,
                    "--exhaustiveness",
                    "2",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--max_evals",
                    str(cap),
                    "--pose_selection",
                    "late",
                    "--search_counters",
                    "--out",
                    str(name),
                ]
                for axis, val in zip(
                    "xyz",
                    json.loads((run / f"prepared/{rec}.json").read_text())["center"],
                ):
                    args += ["--center_" + axis, str(val), "--size_" + axis, "26"]
                start = time.perf_counter()
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                with (run / log).open("w") as output:
                    subprocess.run(
                        [binary, *args],
                        cwd=run,
                        stdout=output,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                now = resource.getrusage(resource.RUSAGE_CHILDREN)
                scores = [
                    float(l.split()[3])
                    for l in (run / name).read_text().splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                ]
                row.update(
                    wall_s=round(time.perf_counter() - start, 6),
                    cpu_s=round(
                        now.ru_utime + now.ru_stime - before.ru_utime - before.ru_stime,
                        6,
                    ),
                    modes=len(scores),
                    score=scores[0],
                    pool_sha256=digest(run / name),
                )
            wr.writerow(row)
            f.flush()
    data = json.loads((frozen / "binaries.json").read_text())
    if executable_hash != digest(binary):
        raise RuntimeError("The executable changed during the run")
    data["pro"] = executable_hash
    (run / "binaries.json").write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("frozen", "run", "binary"):
        p.add_argument(name, type=lambda x: Path(x).resolve())
    a = p.parse_args()
    replay(a.frozen, a.run, a.binary)
