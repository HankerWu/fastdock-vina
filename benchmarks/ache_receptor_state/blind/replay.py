"""Re-run production pools and test the upstream counting observer against the release."""
import argparse
import csv
import json
import shutil
import tempfile
from pathlib import Path

from dock import load_inputs, run_job, save_json, sha, table, tick


def main(result, out, official, pro, counter):
    out.mkdir(parents=True, exist_ok=False)
    start = tick()
    binaries = {"official": official, "pro": pro, "counter": counter}
    receipt = json.loads((result / "search_receipt.json").read_text())["binary_sha256"]
    hashes = {name: sha(exe.read_bytes()) for name, exe in binaries.items()}
    if any(hashes[name] != value for name, value in receipt.items()):
        raise ValueError("Production binary differs from the one being replayed")
    rows = []
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)
        load_inputs(work)
        for row in csv.DictReader((result / "runs.csv").open()):
            rec, arm, seed = (row[k] for k in ("receiver", "arm", "seed"))
            new = out / rec / f"{arm}_{seed}"
            orig = result / "raw" / rec / f"{arm}_{seed}"
            shutil.copyfile(work / f"{rec}.pdbqt", work / "receptor.pdbqt")
            args = json.loads((orig / "command.json").read_text())[1:]
            for kind in (("official", "counter") if arm == "official24k" else ("pro",)):
                replay = run_job(binaries[kind], args, work, new / kind)
                rows.append(dict(receiver=rec, arm=arm, seed=seed, binary=kind, **replay,
                                 same_output=replay["output_sha256"] == row["output_sha256"]))
                table(out / "replays.csv", rows)
    save_json(out / "receipt.json", dict(binary_sha256=hashes, total_wall_cpu_s=(tick()-start).tolist()))
    if any(r["status"] != "ok" or not r["same_output"] for r in rows if r["binary"] != "counter"):
        raise RuntimeError("Production replay differs; raw evidence retained")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for field in ("result", "out", "official", "pro", "counter"):
        p.add_argument(field, type=Path)
    a = p.parse_args()
    main(a.result, a.out, a.official.resolve(), a.pro.resolve(), a.counter.resolve())
