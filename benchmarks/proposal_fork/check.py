"""Recount directional controls directly from trial logs; optionally confirm a separate live run."""

import argparse
import csv
import itertools
import json
import lzma
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent


def check(repeat=None):
    raw = HERE / "raw"
    if repeat:
        for pattern in ("*.xz", "*.csv"):
            for file in raw.glob(pattern):
                assert file.read_bytes() == (repeat / file.name).read_bytes(), file
    pooled = {}
    for file in sorted(raw.glob("*.trials.xz")):
        case = file.name[:4]
        with lzma.open(file, "rt") as stream:
            for _, trials in itertools.groupby(
                csv.reader(stream), key=lambda r: tuple(r[:2])
            ):
                seq = list(trials)
                a = list(map(float, seq[-1]))
                group = pooled.setdefault((case, int(a[-1])), [])
                assert a[10] == a[8] - a[7] and a[9] == a[3] * a[4]
                assert a[-1] == int(a[10] < 1e-4 * a[9])
                assert [int(r[2]) for r in seq] == list(range(len(seq)))
                group.append(a)
    result = {}
    for (case, passed), rows in sorted(pooled.items()):
        result[case + ("-passed" if passed else "-exhausted")] = dict(
            n=len(rows),
            both_negative=sum(r[5] < 0 and r[6] < 0 for r in rows),
            fine_negative=sum(r[6] < 0 for r in rows),
            same_cell_all_atoms=sum(r[15] == 0 for r in rows),
            base_any_outside=sum(r[16] > 0 for r in rows),
            trial_any_outside=sum(r[17] > 0 for r in rows),
            base_inside_to_any_outside=sum(r[16] == 0 and r[17] > 0 for r in rows),
            registered_B=sum(
                r[5] < 0 and r[6] < 0 and (r[14] >= 0.375 or r[15] > 0) for r in rows
            ),
            median_last_RMS_A=statistics.median(r[14] for r in rows),
        )
    assert result["2NNQ-exhausted"]["registered_B"] == 56
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repeat", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    result = check(a.repeat)
    if a.output:
        a.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads((HERE / "controls.json").read_text())
        print("Validated passed/exhausted directional controls.")
