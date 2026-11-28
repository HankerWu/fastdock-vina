"""Stream raw stage records and recover receiver-frame errors and published endpoints."""

import argparse
import json
import lzma
import struct
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from audit import STUDY, Frame, counters, rmsd, sha


def check(folder, frame):
    manifest = json.loads((folder / "runs.json").read_text())
    output = {}
    for rec, meta in manifest.items():

        def stream(kind, rec=rec, meta=meta):
            raw = lzma.decompress(
                b"".join(p.read_bytes() for p in sorted(folder.glob(f"{rec}.{kind}.xz.*")))
            )
            assert sha(raw) == meta[kind + "_sha256"]
            return raw

        data = stream("events")
        n = struct.unpack_from("<I", data)[0]
        serial = struct.unpack_from(f"<{n}I", data, 4)
        idx = [serial.index(s) for s in frame.order]
        records = defaultdict(dict)
        for traj, step, phase, event, index, e, *xyz in struct.iter_unpack(
            f"<5Id{3 * n}d", data[4 * (1 + n) :]
        ):
            x = np.asarray(xyz).reshape(n, 3)[idx]
            assert (traj, step) not in records[event]
            records[event][traj, step] = (x, e, rmsd(x, frame.targets[rec]))
        src = STUDY / "blind/results/raw" / rec / "pro24k_1103503"
        text = (src / "dry.pdbqt").read_bytes()
        assert sha(text) == meta["output_sha256"]
        assert counters((folder / (rec + ".log")).read_text()) == counters((src / "dry.log").read_text())
        accepted = {}
        for row in stream("members").splitlines():
            t, s, *retained = map(int, row.split())
            assert (t, s) in records[6]
            accepted[t] = retained
        assert {(t, s) for t, kept in accepted.items() for s in kept} == records[12].keys()
        source_ids = {}
        for key, (x, e, d) in records[9].items():
            found = [
                origin
                for origin, (y, energy, _) in records[12].items()
                if e == energy and np.array_equal(x, y)
            ]
            assert len(found) == 1
            source_ids[key] = found[0]
        released = []
        for i, block in enumerate(text.decode().split("ENDMDL")[:-1], 1):
            x = frame.heavy(block)
            ids = [key for key, (y, _, _) in records[10].items() if np.max(np.abs(x - y)) < 0.000501]
            assert len(ids) == 1
            released.append(source_ids[ids[0]])
        all_near = {k: [key for key, (_, _, d) in g.items() if d <= 2] for k, g in records.items()}
        info = {
            "proposals": len(records[11]),
            "near_proposals": len(all_near[11]),
            "near_hunt": len(all_near[4]),
            "near_accepted": len(all_near[5]),
            "near_second": len(all_near[6]),
            "near_final": len(all_near[12]),
            "near_merged": len(all_near[8]),
            "checked_release_origins": released,
            "output_sha256": sha(text),
        }
        if rec == "4EY7":
            assert len(all_near[11]) == 1
            (key,) = all_near[11]
            info["sole_near_proposal"] = {
                "trajectory_step": key,
                "proposal_A": records[11][key][2],
                "hunt_A": records[4][key][2],
                "accepted": key in records[5],
                "has_second": key in records[6],
            }
        output[rec] = info
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    w, c = time.perf_counter(), time.process_time()
    with tempfile.TemporaryDirectory() as tmp:
        results = check(STUDY / "trace/results", Frame(Path(tmp)))
    results.update(wall_s=time.perf_counter() - w, cpu_s=time.process_time() - c)
    args.out.write_text(json.dumps(results, indent=2) + "\n")
