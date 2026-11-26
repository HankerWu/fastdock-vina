"""Measure 1FAX cross-docking and deliberately anchor its opposite amidine group."""

import argparse
import itertools
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from dock import ROOT, SEEDS, write, run

sys.path.insert(0, str(ROOT / "benchmarks/directional_ion_pair"))
from cross_dock import alignment, structure_lines
from verify_cross import atoms, coords, mappings
from evaluate import pose_models
from experiment import aspsite
from reversed_head import inverted_anchor


def evaluate(base, out, binary, jobs, parity):
    work = base / "fdv-ionpair"
    crys = work / "raw/1FAX.pdb"
    src = work / "cross/DX9.pdbqt"
    heavy = [l for l in atoms(src) if l.split()[-1] != "HD"]
    x = {int(l[6:11]): coords(l) for l in heavy}
    z = {int(l[6:11]): l.split()[-1] for l in heavy}
    ref = [
        l
        for l in structure_lines(crys)
        if l[17:20] == "DX9" and l[76:78].strip() != "H"
    ]
    xyz = np.asarray([coords(l) for l in ref])
    elements = {i: l[76:78].strip() for i, l in enumerate(ref)}
    maps = mappings(x, z, dict(enumerate(xyz)), elements)
    assert len(maps) == 8 and len(x) == 33
    rows = []
    routing = []
    for rec, folder in (("1FJS", base / "fdv-fxa1"), ("2P16", base / "fdv-fxa")):
        dest = out / "cross_wrong"
        dest.mkdir(parents=True, exist_ok=True)
        site = dest / f"{rec}_DX9.inverted"
        site.write_text(inverted_anchor(src, aspsite(folder / f"{rec}.pdb", "A")[1:]))

        def wrong(seed):
            prior = dest / f"{rec}_DX9_{seed}_inverted.log"
            pure = prior.with_suffix(".pdbqt")
            with prior.open("w") as stream:
                cmd = [
                    binary,
                    "--ligand",
                    src,
                    "--receptor",
                    folder / f"{rec}.pdbqt",
                    "--config",
                    folder / "box.txt",
                    "--seed",
                    str(seed),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "4",
                    "--max_evals",
                    "50000",
                    "--pose_selection",
                    "late",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--search_counters",
                    "--search_anchor",
                    site,
                    "--out",
                    pure,
                ]
                subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, check=True)
            return run(
                (
                    "cross_wrong",
                    rec,
                    "DX9",
                    seed,
                    src,
                    folder / f"{rec}.pdbqt",
                    folder / "box.txt",
                    site,
                    work / "cross" / f"{rec}_DX9_{seed}_native.log",
                    prior,
                ),
                binary,
                out,
                parity,
            )

        with ThreadPoolExecutor(jobs) as pool:
            routing.extend(pool.map(wrong, SEEDS))
        rot, shift, k, ca_rms = alignment(
            crys, folder / f"{rec}.pdb", aspsite(crys, "A")[0]
        )
        assert k == 74
        target = xyz @ rot + shift
        for s, arm in itertools.product(
            SEEDS, ("native", "shell", "mixed", "inverted", "mixed_inverted")
        ):
            directory, name = (
                (work / "cross", arm)
                if arm in ("native", "shell")
                else (
                    (out / "cross", "mixed")
                    if arm == "mixed"
                    else (dest, "mixed" if arm == "mixed_inverted" else "inverted")
                )
            )
            rms = []
            for score, points in pose_models(directory / f"{rec}_DX9_{s}_{name}.pdbqt"):
                rms.append(
                    min(
                        float(
                            np.sqrt(
                                np.mean(
                                    [
                                        np.square(points[i] - target[j]).sum()
                                        for i, j in m.items()
                                    ]
                                )
                            )
                        )
                        for m in maps
                    )
                )
            rows.append(
                dict(
                    receptor=rec,
                    ligand="DX9",
                    seed=s,
                    arm=arm,
                    ca_rms=ca_rms,
                    top=rms[0],
                    oracle=min(rms),
                    modes=len(rms),
                )
            )
    write(out / "routing_cross_wrong.csv", routing)
    write(out / "cross_cells.csv", rows)
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for n in ("base", "binary", "out"):
        p.add_argument(n, type=Path)
    p.add_argument("--parity", choices=("odd", "even"), default="even")
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    result = evaluate(a.base, a.out, a.binary.resolve(), a.jobs, a.parity)
    for arm, rs in itertools.groupby(
        sorted(result, key=lambda r: r["arm"]), key=lambda r: r["arm"]
    ):
        group = list(rs)
        print(
            arm,
            sum(r["top"] <= 2 for r in group),
            sum(r["oracle"] <= 2 for r in group),
            len(group),
        )
