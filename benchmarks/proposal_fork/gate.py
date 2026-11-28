"""Post hoc replay of the observed default save predicate for the two reverse-ESR1 paths."""

import argparse
import csv
import io
import json
import lzma
from pathlib import Path
import tempfile
import statistics

import numpy as np
from rescue import HERE, fixtures, native, run_panel


def check(work, historical):
    case = "1X7R-1QKT-85231"
    _, _, inputs, q, r, _, _ = next(
        row for row in fixtures(work, historical) if row[0] == case
    )
    ref = run_panel.replay.reciprocal.reference(q, r)
    serials = [
        int(l[6:11])
        for l in (inputs / q["pdb"] / "ligand.pdbqt").read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
        and l[77:].strip() not in {"H", "HD", "HS"}
    ]
    order = [serials.index(i) for i in q["serial"]]
    output, overview = [], {}
    for arm in ("stable", "last_trial"):
        path = work / f"{arm}.native"
        path.write_bytes(
            lzma.decompress((HERE / "free" / f"{case}-{arm}.xz").read_bytes())
        )
        events, _ = native(path)
        accepted, entered = 0, 0
        rejected_near = []
        for t, trajectory in enumerate(events, 1):
            best = float("inf")
            pool_size = 0
            for i, step in enumerate(trajectory):
                if not step["A"]:
                    assert not step["D"]
                    continue
                accepted += 1
                h = step["H"][0]
                energy = float(h[2])
                trigger = energy < best or pool_size < 9
                assert trigger == bool(step["D"]), (arm, t, i, best, pool_size)
                coords = (
                    np.fromstring(h[8], dtype=int, sep=" ").reshape(-1, 3)[order] / 1000
                )
                rms = float(
                    np.sqrt(np.square(coords[None] - ref).sum(-1).mean(-1).min())
                )
                if rms <= 2:
                    gap = energy - best
                    output.append(
                        dict(
                            arm=arm,
                            trajectory=t,
                            step=i,
                            hunt_A=rms,
                            hunt_energy=energy,
                            prior_best_energy=best,
                            hunt_minus_best=gap,
                            pool_size=pool_size,
                            entered=trigger,
                        )
                    )
                    if not trigger:
                        rejected_near.append(gap)
                if trigger:
                    entered += 1
                    pool_size = len(step["M"])
                    # The native default saves into the same tmp object.
                    best = min(best, float(step["D"][0][2]))
        overview[arm] = dict(
            all_accepted_checked=accepted,
            all_save_calls=entered,
            near_blocked=len(rejected_near),
            blocked_energy_gap_min=min(rejected_near),
            blocked_energy_gap_median=statistics.median(rejected_near),
            blocked_energy_gap_max=max(rejected_near),
        )
    return output, overview


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history", type=Path)
    p.add_argument("--write", action="store_true")
    a = p.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        rows, summary = check(Path(tmp), a.history)
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = lzma.compress(buf.getvalue().encode(), preset=3)
    packed = HERE / "free/save_gate.csv.xz"
    report = HERE / "free/save_gate.json"
    if a.write:
        packed.write_bytes(data)
        report.write_text(json.dumps(summary, indent=2) + "\n")
    else:
        assert packed.read_bytes() == data and json.loads(report.read_text()) == summary
    print(summary)
