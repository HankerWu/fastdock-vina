"""Check the recorded heavy-atom order in the pose, FKBP12 and acidic-state runs."""

import argparse
import csv
from pathlib import Path
import tarfile

import numpy as np


def check(pose, fkbp, acid, dest):
    def record(series, label, serials):
        ids = np.asarray(serials, dtype=int)
        if len(ids) == 0 or not np.all(np.diff(ids) > 0):
            raise ValueError(f"unsorted/empty raw heavy-atom serial list at {label}")
        return label, series, len(ids), True

    def arrays(path, series):
        for name in sorted(path.glob("*.npz")):
            with np.load(name) as raw:
                yield record(series, name.stem, raw["serials"])

    rows = list(arrays(pose,"pose"))
    with tarfile.open(fkbp) as tar:
        for name in sorted(tar.getnames()):
            if name.startswith("capture/") and name.endswith(".xyz"):
                header = tar.extractfile(name).readline()
                rows.append(record("FKBP12", Path(name).stem, header.split()[1:]))
    rows += list(arrays(acid,"acid_microstate"))
    if len(rows) != 46:
        raise ValueError("expected 24 pose, 16 FKBP and 6 acidic-state sources")
    with dest.open("w", newline="") as stream:
        out = csv.writer(stream)
        out.writerow(("source","panel","heavy_atoms","strictly_increasing_pdbqt_serial"))
        out.writerows(rows)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pose", "fkbp", "acid", "output"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    check(a.pose, a.fkbp, a.acid, a.output)
