"""Audit the 1FAX cross-docking RMSDs using positional topology and a quaternion alignment."""

import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation


def atoms(path):
    return [
        l
        for l in path.read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM")) and l[16] in " A"
    ]


def coords(line):
    return np.array(list(map(float, (line[30:38], line[38:46], line[46:54]))))


def topology(coords):
    return {
        i: {
            j for j in coords if j != i and np.linalg.norm(coords[j] - coords[i]) < 1.70
        }
        for i in coords
    }


def mappings(src, src_z, dst, dst_z):
    g = topology(src)
    h = topology(dst)
    canon = lambda z: {"A": "C", "OA": "O", "NA": "N"}.get(z, z)

    def sig(i, graph, z):
        return (canon(z[i]), tuple(sorted(canon(z[j]) for j in graph[i])))

    choices = {i: [j for j in dst if sig(i, g, src_z) == sig(j, h, dst_z)] for i in src}

    def visit(m):
        if len(m) == len(g):
            yield m
            return
        i = min(
            (k for k in src if k not in m),
            key=lambda k: (-sum(nei in m for nei in g[k]), len(choices[k]), k),
        )
        assigned = set(m.values())
        for j in choices[i]:
            if j in assigned or any((a in g[i]) != (b in h[j]) for a, b in m.items()):
                continue
            yield from visit({**m, i: j})

    return list(visit({}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "fjs", "p16", "results"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    crys = [l for l in atoms(a.work / "raw" / "1FAX.pdb") if l[21] == "A"]
    xtal = {
        l[12:16].strip(): coords(l)
        for l in crys
        if l[17:20] == "DX9" and l[76:78].strip() != "H"
    }
    xtal_z = {
        l[12:16].strip(): l[76:78].strip()
        for l in crys
        if l[17:20] == "DX9" and l[76:78].strip() != "H"
    }
    prep = [l for l in atoms(a.work / "cross" / "DX9.pdbqt") if l.split()[-1] != "HD"]
    src = {int(l[6:11]): coords(l) for l in prep}
    src_z = {int(l[6:11]): l.split()[-1] for l in prep}
    maps = mappings(src, src_z, xtal, xtal_z)
    assert len(maps) == 8 and all(len(g) == 33 for g in maps)
    ast = next(
        coords(l)
        for l in crys
        if l[17:20] == "ASP" and int(l[22:26]) == 189 and l[12:16].strip() == "CG"
    )
    pool = list(csv.DictReader(a.results.open()))
    discrepancy = 0
    nposes = 0
    for rec, root in (("1FJS", a.fjs), ("2P16", a.p16)):
        ca = {
            l[22:27]: coords(l)
            for l in crys
            if l.startswith("ATOM  ")
            and l[12:16].strip() == "CA"
            and np.linalg.norm(coords(l) - ast) < 15
        }
        ct = {
            l[22:27]: coords(l)
            for l in atoms(root / f"{rec}.pdb")
            if l.startswith("ATOM  ") and l[21] == "A" and l[12:16].strip() == "CA"
        }
        shared = sorted(ca.keys() & ct.keys())
        assert len(shared) == 74
        mobile = np.array([ca[k] for k in shared])
        fixed = np.array([ct[k] for k in shared])
        rot, _ = Rotation.align_vectors(fixed - fixed.mean(0), mobile - mobile.mean(0))
        target = {
            i: rot.apply(v - mobile.mean(0)) + fixed.mean(0) for i, v in xtal.items()
        }
        for row in pool:
            if row["receptor"] != rec:
                continue
            path = a.work / "cross" / f"{rec}_DX9_{row['seed']}_{row['arm']}.pdbqt"
            values = []
            for block in path.read_text().split("MODEL")[1:]:
                xx = {
                    int(l[6:11]): coords(l)
                    for l in block.splitlines()
                    if l.startswith(("ATOM  ", "HETATM"))
                }
                values.append(
                    min(
                        np.sqrt(
                            sum(
                                np.square(xx[i] - target[j]).sum() for i, j in m.items()
                            )
                            / len(m)
                        )
                        for m in maps
                    )
                )
            nposes += len(values)
            discrepancy = max(
                discrepancy,
                abs(values[0] - float(row["rms_top"])),
                abs(min(values) - float(row["rms_best"])),
            )
    assert discrepancy < 1e-9
    print(
        f"{len(pool)} FXa cross-dock pools and {nposes} models. Independent topology and alignment reproduce top/pool heavy-atom RMSD to {discrepancy:.3g} A."
    )


if __name__ == "__main__":
    main()
