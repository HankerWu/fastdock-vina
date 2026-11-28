"""Reconstruct stage minima and frequencies from real production-equivalent free trajectories."""

import argparse
import json
import lzma
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

from rescue import HERE, argv, csvfile, digest, fixtures, fold, native, run_panel

sys.path.insert(0, str(HERE.parent / "stage_local_attribution"))
from analyze_stages import observe_trajectory


def main(work, out, binaries, history):
    live = binaries is not None
    if live:
        out.mkdir(parents=True, exist_ok=False)
    results, costs = [], []
    historical = json.loads((HERE / "raw/costs.json").read_text())
    for case, driver, inputs, q, r, seed, _ in fixtures(work, history):
        if driver not in ("stable", "2326430"):
            continue
        ref = run_panel.replay.reciprocal.reference(q, r)
        serials = [
            int(l[6:11])
            for l in (inputs / q["pdb"] / "ligand.pdbqt").read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM"))
            and l[77:].strip() not in {"H", "HD", "HS"}
        ]
        order = [serials.index(i) for i in q["serial"]]

        def error(row):
            xyz = np.fromstring(row[8], dtype=int, sep=" ").reshape(-1, 3)[order] / 1000
            return float(np.sqrt(np.square(xyz[None] - ref).sum(-1).mean(-1).min()))

        for i, arm in enumerate(("stable", "last_trial")):
            key = f"{case}-{arm}"
            trace = work / f"{key}.native"
            packed = out / f"{key}.xz"
            published = out / f"{key}.pdbqt"
            if live:
                for kind, binary in zip(
                    ("stock", "observed"), binaries[2 * i : 2 * i + 2]
                ):
                    for k in ("FD_BFGS_SHADOW", "FD_FORK_COUNTERFACTUAL"):
                        os.environ.pop(k, None)
                    os.environ["FD_BFGS_TRACE"] = str(trace)
                    pose = work / f"{key}-{kind}.pdbqt"
                    cmd = argv(binary, inputs, q, r, seed, pose) + ["--search_counters"]
                    stats = run_panel.execute(cmd, out / f"{key}-{kind}.log", 180)
                    assert stats["status"] == "ok"
                    if kind == "stock":
                        published.write_bytes(pose.read_bytes())
                    else:
                        assert pose.read_bytes() == published.read_bytes()
                        packed.write_bytes(lzma.compress(trace.read_bytes(), preset=3))
                    costs.append(
                        dict(
                            key=key,
                            kind=kind,
                            argv=cmd,
                            binary=run_panel.digest(binary),
                            pose=run_panel.digest(pose),
                            ligand=run_panel.digest(inputs / q["pdb"] / "ligand.pdbqt"),
                            receptor=run_panel.digest(
                                inputs / r["pdb"] / "receptor.pdbqt"
                            ),
                            **stats,
                        )
                    )
                    if kind == "observed":
                        costs[-1]["trace"] = run_panel.digest(packed)
            else:
                receipts = [
                    row
                    for row in json.loads((out / "costs.json").read_text())
                    if row["key"] == key
                ]
                assert len(receipts) == 2 and all(
                    row["pose"] == run_panel.digest(published) for row in receipts
                )
                assert receipts[1]["trace"] == run_panel.digest(packed)
                trace.write_bytes(lzma.decompress(packed.read_bytes()))
            if arm == "stable" and driver == "stable":
                record = next(
                    v for v in historical if v["key"] == f"{case}-stable-free"
                )
                assert run_panel.digest(trace) == record["trace"]
            paths, _ = native(trace)
            allsteps = fold(paths)
            retained = []
            for trajectory in paths:
                terminal, _, _, _ = observe_trajectory(trajectory, error)
                retained.extend(value[1] for value in terminal.values())
            ranks = run_panel.replay.reciprocal.measure(published, q["serial"], ref)
            result = dict(
                case=case,
                arm=arm,
                Top1_A=ranks[0],
                Best5_A=min(ranks[:5]),
                Oracle_A=min(ranks),
                Regret_A=ranks[0] - min(ranks),
            )
            for phase, values in {
                "hunt": [error(step["H"][0]) for step in allsteps],
                "accepted": [error(step["H"][0]) for step in allsteps if step["A"]],
                "post_save": [error(row) for step in allsteps for row in step["D"]],
                "retained": retained,
                "published": ranks,
            }.items():
                result.update(
                    {
                        phase + "_n": len(values),
                        phase + "_near": sum(x <= 2 for x in values),
                        phase + "_min_A": min(values) if values else "",
                    }
                )
            accepted_near = [s for s in allsteps if s["A"] and error(s["H"][0]) <= 2]
            result["accepted_near_triggered_save"] = sum(bool(s["D"]) for s in accepted_near)
            result["accepted_near_saved_near"] = sum(error(row) <= 2 for s in accepted_near for row in s["D"])
            results.append(result)
            print(result, flush=True)
    assert len(results) == 6
    csvfile(out / "stages.csv", results, live)
    if live:
        (out / "costs.json").write_text(json.dumps(costs, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--binaries",
        nargs=4,
        type=Path,
        metavar=("STABLE", "STABLE_OBS", "NEW", "NEW_OBS"),
    )
    p.add_argument("--history", type=Path)
    a = p.parse_args()
    with tempfile.TemporaryDirectory() as work:
        main(Path(work), a.out.resolve(), a.binaries, a.history)
