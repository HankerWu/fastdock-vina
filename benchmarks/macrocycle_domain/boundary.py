"""Replay real fixed coordinates across box faces; AD4 uses the actual archived map."""
import argparse
import csv
import hashlib
import itertools
import re
import subprocess
import tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
SCORE = re.compile(r"Estimated Free Energy of Binding\s+:\s+([-+\d.]+)")
GLUE = {"G0", "G1", "G2", "G3"}
HYDROGEN = {"HD", "HS", "H"}


def atom_data(block):
    lines = [l for l in block.splitlines() if l.startswith(("ATOM  ", "HETATM"))]
    return [(l.split()[-1], [float(l[i:i+8]) for i in (30, 38, 46)]) for l in lines]


def pose_blocks(data):
    return [block.strip().split("\n", 1)[1] for block in data.split("ENDMDL") if block.strip()]


def call(binary, args):
    run = subprocess.run([binary, "--cpu", "1", "--seed", "1259", *args], text=True, capture_output=True)
    score = SCORE.search(run.stdout)
    if run.returncode:
        if "outside the grid box" not in run.stdout + run.stderr:
            raise RuntimeError(run.stdout + run.stderr)
        return run.returncode, None
    if not score:
        raise RuntimeError(run.stdout + run.stderr)
    return 0, float(score.group(1))


def audit_native_distance(archive, output):
    # Finalized pools are incomplete search histories. Deduplicate identical serialized pools.
    seen = set()
    rows = []
    for name in sorted(n for n in archive.getnames() if n.endswith("_parent.pdbqt")):
        data = archive.extractfile(name).read()
        digest = hashlib.sha256(data).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        blocks = [atom_data(block) for block in pose_blocks(data.decode())]
        for (i, a), (j, b) in itertools.combinations(enumerate(blocks, 1), 2):
            types = [t for t, v in a]
            assert types == [t for t, v in b]
            delta = np.asarray([v for t, v in a]) - np.asarray([v for t, v in b])
            squared = np.sum(delta * delta, axis=1)
            def rmsd(omit):
                return float(np.sqrt(np.mean(squared[[t not in omit for t in types]])))
            all_heavy, chemical = rmsd(HYDROGEN), rmsd(HYDROGEN | GLUE)
            rows.append(dict(source=name, rank_i=i, rank_j=j, with_glue_A=all_heavy, without_glue_A=chemical,
                             disagree_1A=(chemical < 1) != (all_heavy < 1),
                             disagree_2A=(chemical < 2) != (all_heavy < 2)))
    with (output / "archived_distances.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("inputs", "output", "parent", "fixed"):
        parser.add_argument(key, type=Path)
    a = parser.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    results = []
    def evaluate(pdb, source, rank, case, geometry, args):
        for arm in ("parent", "fixed"):
            rc, score = call(getattr(a, arm), args)
            results.append(dict(pdb=pdb, source=source, rank=rank, case=case, geometry=geometry,
                                arm=arm, status=rc, score=score))
    with tarfile.open(BENCH / "concerted_ring/coordinates.tar.xz") as tar:
        audit_native_distance(tar, a.output)
        for pdb in ("3VHA", "6QBH", "8H6P"):
            candidates = sorted(n for n in tar.getnames() if n.startswith(("transfer/raw/", "catd/raw/"))
                                and f"/{pdb}_" in n and n.endswith("_parent.pdbqt"))
            for name in candidates:
                for rank, block in enumerate(pose_blocks(tar.extractfile(name).read().decode()), 1):
                    data = atom_data(block)
                    heavy = np.array([v for t, v in data if t not in GLUE | HYDROGEN])
                    glue = np.array([v for t, v in data if t in GLUE])
                    lo, hi = heavy.min(0), heavy.max(0)
                    extras = np.r_[lo - glue.min(0), glue.max(0) - hi]
                    j = extras.argmax()
                    if extras[j] > 0.05:
                        break
                else:
                    continue
                break
            else:
                results.append(dict(pdb=pdb, source="", rank="", case="no eligible real pose",
                                    geometry="", arm="", status="", score=""))
                continue
            pose = a.output / (pdb + "-fixed.pdbqt")
            pose.write_text(block)
            common = ["--receptor", a.inputs/pdb/"receptor.pdbqt", "--ligand", a.inputs/pdb/"ligand.pdbqt",
                      "--pose_coordinates", pose, "--score_only"]
            for label in ("wide", "marker_outside", "real_outside"):
                d = j % 3
                center = (lo + hi) / 2
                if label != "wide":
                    margin = extras[j] / 2 if label == "marker_outside" else -0.25
                    center[d] = lo[d] - margin + 15 if j < 3 else hi[d] + margin - 15
                if label == "marker_outside":
                    (a.output / (pdb + ".box")).write_text("".join(
                        f"center_{dim} = {v}\nsize_{dim} = 30\n" for dim, v in zip("xyz", center)))
                args = common.copy()
                for axis, value in zip("xyz", center):
                    args += ["--center_"+axis, str(value), "--size_"+axis, "30"]
                evaluate(pdb, name, rank, label, f"axis={'xyz'[d]}, face={'min' if j<3 else 'max'}, "
                         f"marker_extends_beyond_heavy_A={extras[j]:.3f}", args)
    with tarfile.open(BENCH / "cubic_grid/fkbp_ad4/maps-1FKH.tar.xz") as tar:
        tar.extractall(a.output, filter="data")
    with tarfile.open(BENCH / "cubic_grid/fkbp_ad4/pools.tar.xz") as tar:
        name = "pools/10_1FKH_131071_ad4.pdbqt"
        inside = pose_blocks(tar.extractfile(name).read().decode())[0]
    for shift in (0, 100):
        pose = a.output / f"1FKH-{shift}.pdbqt"
        pose.write_text("\n".join(l[:30] + f"{float(l[30:38]) + shift:8.3f}" + l[38:]
                                  if l.startswith(("ATOM  ", "HETATM")) else l for l in inside.splitlines())+"\n")
        for operation in ("score_only", "local_only"):
            args = ["--maps", a.output/"maps/1FKH/receptor", "--ligand", pose, "--scoring", "ad4", "--"+operation]
            if operation == "local_only":
                args += ["--out", a.output/"optimized.pdbqt"]
            evaluate("1FKH", name, 1, operation, f"x_shift_A={shift}", args)
    with (a.output / "boundary.csv").open("w") as out:
        writer = csv.DictWriter(out, results[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    by_case = {(r["pdb"], r["case"], r["geometry"], r["arm"]): (r["status"], r["score"]) for r in results}
    for r in results:
        if r["case"] == "marker_outside":
            wide = by_case[r["pdb"], "wide", r["geometry"], r["arm"]]
            assert wide[0] == 0 and (r["status"] != 0 if r["arm"] == "parent" else r["score"] == wide[1])
        if r["case"] == "real_outside" or (r["pdb"] == "1FKH" and r["geometry"] == "x_shift_A=100" and r["arm"] == "fixed"):
            assert r["status"] != 0


if __name__ == "__main__":
    main()
