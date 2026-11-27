"""Post-hoc label-assisted control: refine the best possible placement of the same input shape."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from audit import read_atoms, reference
from prepare_transfer import BENCH, DIRECTIONS, SYSTEMS
from replay import save

sys.path.insert(0, str(BENCH / "rigid_translation_scan"))
from scan import read_ligand, write_start


def run(base, work, probe):
    meta = {c: json.loads((base / c / "preparation.json").read_text()) for c in SYSTEMS}
    work.mkdir(exist_ok=True, parents=True)
    results = []
    for rec, query in DIRECTIONS:
        if meta[rec]["receptor_status"] != "ready":
            continue
        path = base / query / "ligand.pdbqt"
        topol, atoms, all_xyz, _, _ = read_ligand(probe, path)
        ids = meta[query]["serial"]
        relative = read_atoms(path.read_text(), ids)
        relative -= relative.mean(0)
        refs, _, _ = reference(meta[query], meta[rec])
        choices = [Rotation.align_vectors(ref - ref.mean(0), relative) for ref in refs]
        idx = min(range(len(choices)), key=lambda k: choices[k][1])
        rot, rss = choices[idx]
        start = work / f"{rec}_{query}.pdbqt"
        write_start(start, topol, atoms, rot.apply(all_xyz) + refs[idx].mean(0))
        output = subprocess.check_output([
            str(probe), "refine", str(base / rec / "receptor.pdbqt"), str(path),
            *map(str, meta[rec]["center"]), str(start)], text=True)
        start.with_suffix(".log").write_text(output)
        line = next(l.split() for l in output.splitlines() if l.startswith("result "))
        refined = read_atoms(Path(str(start) + ".refined").read_text(), ids)
        results.append(dict(receiver=rec, query=query, initial_shape_rmsd=rss / np.sqrt(len(ids)),
                            refined_rmsd=np.sqrt(np.square(refined - refs).sum(-1).mean(-1)).min(),
                            explicit_before=float(line[2]), explicit_after=float(line[3])))
    save(work / "fit.csv", results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("base", "work", "probe"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    run(a.base, a.work, a.probe)
