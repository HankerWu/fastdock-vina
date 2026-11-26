"""Independent QMC design, rigid transforms, full-grid evaluations and boundary diagnosis."""

import argparse
import csv
import io
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp
from scipy.stats import qmc
from run import table, write

RT = 0.00198720425864 * 298.15
THETA = math.radians(15)


def nodes(n):
    # Skip the initial zero Sobol point, as in boost::random::sobol.
    sample = (
        qmc.Sobol(d=6, scramble=False).random_base2(int(math.log2(n)) + 1)[1 : n + 1].T
    )

    def axis(u, v):
        az = 2 * math.pi * v
        z = 2 * u - 1
        xy = np.sqrt(1 - z * z)
        return np.stack((xy * np.cos(az), xy * np.sin(az), z), axis=1)

    space = axis(sample[1], sample[2])
    spin = axis(sample[4], sample[5])
    translations = []
    rotations = []
    logweights = []
    vol = lambda t: t - math.sin(t)
    for a in (1, 0.5, 0.25):
        r = a * sample[0] ** (1 / 3)
        t = np.array(
            [
                brentq(lambda x: vol(x) - q * vol(a * THETA), 0, a * THETA)
                for q in sample[3]
            ]
        )
        density = sum(
            ((r <= b) & (t <= b * THETA)) * vol(THETA) / (3 * b**3 * vol(b * THETA))
            for b in (1, 0.5, 0.25)
        )
        translations.extend(r[:, None] * space)
        rotations.extend(t[:, None] * spin)
        logweights.extend(-np.log(density))
    return np.array(translations), np.array(rotations), np.array(logweights)


def compare(case, binary, ref, n, temp):
    receptor, ligand, pool, model = case
    lines = pool.read_text().split("ENDMDL")[model - 1].strip().splitlines()
    ids = [i for i, line in enumerate(lines) if line.startswith(("ATOM", "HETATM"))]
    atomcoords = np.array(
        [[float(lines[i][j : j + 8]) for j in (30, 38, 46)] for i in ids]
    )
    heavy = atomcoords[
        [j for j, i in enumerate(ids) if lines[i].split()[-1] not in ("H", "HD", "HS")]
    ]
    center = heavy.mean(0)
    v = heavy - center
    x = v[np.linalg.norm(v, axis=1).argmax()].copy()
    x /= np.linalg.norm(x)
    proj = v - (v @ x)[:, None] * x
    y = proj[np.linalg.norm(proj, axis=1).argmax()].copy()
    y /= np.linalg.norm(y)
    F = np.stack((x, y, np.cross(x, y)))
    trans, rot, lw = nodes(n)
    transf = temp / "transformed.pdbqt"
    one = temp / "one.pdbqt"
    one.write_text("\n".join(lines) + "\nENDMDL\n")
    with transf.open("w") as out:
        out.write(one.read_text())
        for t, spin in zip(trans, rot, strict=True):
            xyz = (
                Rotation.from_rotvec(spin @ F).apply(atomcoords - center)
                + center
                + t @ F
            )
            new = lines.copy()
            for i, p in zip(ids, xyz, strict=True):
                new[i] = lines[i][:30] + "".join(f"{v:8.3f}" for v in p) + lines[i][54:]
            out.write("\n".join(new) + "\nENDMDL\n")
    direct = np.fromstring(
        subprocess.check_output([ref, receptor, ligand, transf], text=True), sep="\n"
    )
    out = subprocess.check_output(
        [
            binary,
            "--receptor",
            receptor,
            "--ligand",
            ligand,
            "--score_only",
            "--rigid_window_pool",
            one,
            "--rigid_window_samples",
            str(n),
            "--verbosity",
            "0",
            "--cpu",
            "1",
        ],
        text=True,
    )
    row = next(csv.DictReader(io.StringIO(out[out.index("model,inter,intra,") :])))
    delta = -RT * (
        logsumexp(lw - (direct[1:] - direct[0]) * float(row["slope"]) / RT)
        - logsumexp(lw)
    )
    assert abs(float(row["inter"]) - direct[0]) < 2e-8
    if abs(delta - float(row["window_delta"])) > 0.01:
        raise ValueError((pool, model, delta, row))
    return dict(
        pool=pool.name,
        model=model,
        samples=3 * n,
        full_grid_center=direct[0],
        full_grid_delta=delta,
        accelerated_delta=row["window_delta"],
        absolute_error=abs(delta - float(row["window_delta"])),
    )


def boundary(binary, base, results, replayed=None):
    rows = []
    for panel, source, as_run, rec, lig, seed in (
        ("brd4", "fdv-brd4", "fdv-map-brd4", "4OGI", "4OGI", "141421"),
        ("alkyl", "fdv-alkyl", "fdv-alkyl/runs", "A", "4W54", "141421"),
    ):
        inputs = base / source / "prepared"
        pool = base / as_run / f"{rec}_{lig}_{seed}_multi.pdbqt"
        row = next(
            x
            for x in table(results / (panel + ".csv"))
            if (x["receptor"], x["ligand"], x["seed"], x["arm"], x["model"])
            == (rec, lig, seed, "multi", "9")
        )
        with tempfile.TemporaryDirectory() as directory:
            one = Path(directory) / "fixed.pdbqt"
            one.write_text(pool.read_text().split("ENDMDL")[8].strip() + "\nENDMDL\n")
            cmd = [
                binary,
                "--score_only",
                "--ligand",
                inputs / f"{lig}_ligand.pdbqt",
                "--receptor",
                inputs / f"{rec}_receptor.pdbqt",
                "--pose_coordinates",
                one,
                "--cpu",
                "1",
            ]
            box = inputs / ("box.txt" if panel == "alkyl" else rec + "_box.txt")
            regular = subprocess.run(
                cmd + ["--config", box], capture_output=True, text=True
            )
            if regular.returncode and "outside the grid box" not in regular.stderr:
                raise ValueError(regular.stderr)
            vals = []
            for flags in (["--config", box], ["--autobox"]):
                text = subprocess.check_output(
                    cmd + flags + ["--unbound_energy", "0"], text=True
                )
                vals.append(
                    float(re.search(r"Ligand - Receptor\s+:\s+([-0-9.]+)", text)[1])
                )
            repaired = {}
            if replayed is not None:
                change = next(
                    c
                    for c in table(replayed / "changed_models.csv")
                    if c["panel"] == panel
                    and c["pool"] == pool.stem
                    and c["old_rank"] == "9"
                )
                rank = int(change["new_rank"])
                new_block = (
                    (replayed / panel / pool.name).read_text().split("ENDMDL")[rank - 1]
                )
                one.write_text(new_block.strip() + "\nENDMDL\n")
                repaired_text = subprocess.check_output(
                    cmd + ["--config", box], text=True
                )
                reread = float(
                    re.search(r"Ligand - Receptor\s+:\s+([-0-9.]+)", repaired_text)[1]
                )
                scored = float(change["new_inter"])
                if abs(scored - reread) > 0.02:
                    raise ValueError((pool, rank, scored, reread))
                repaired = dict(
                    new_rank=rank, new_saved_inter=scored, new_ordinary_inter=reread
                )
        rows.append(
            dict(
                pool=pool.name,
                model=9,
                saved=row["source_inter"],
                ordinary_score_rejected=int(bool(regular.returncode)),
                original_box=vals[0],
                autobox=vals[1],
                window_center=row["inter"],
                **repaired,
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "binary", "reference", "results"):
        p.add_argument(name, type=Path)
    p.add_argument(
        "--replayed",
        type=Path,
        help="also re-score the repaired models in the original box",
    )
    a = p.parse_args()
    base = a.base
    inputs = base / "fdv-brd4/prepared"
    with tempfile.TemporaryDirectory() as tmp:
        checks = [
            compare(
                (
                    inputs / f"{r}_receptor.pdbqt",
                    inputs / f"{l}_ligand.pdbqt",
                    base / "fdv-map-brd4" / f"{r}_{l}_314159_multi.pdbqt",
                    i,
                ),
                a.binary.resolve(),
                a.reference.resolve(),
                128,
                Path(tmp),
            )
            for r, l, i in (
                ("3MXF", "4O76", 1),
                ("4OGI", "4O77", 2),
                ("4OGJ", "4OGI", 3),
            )
        ]
    write(a.results / "independent.csv", checks)
    write(
        a.results / "boundary.csv",
        boundary(a.binary.resolve(), base, a.results, a.replayed),
    )
    print(checks)
