"""Rescore every candidate with the original covalent graph, movable torsions and a common zero."""

import argparse
import csv
import io
import itertools
import json
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

from run import ROOT, SEEDS, FXA, sources
from amide_rotamers import amides, dihedral_degrees as angle
from select_zinc_pose import read_models, parse_atoms

sys.path.insert(0, str(ROOT / "benchmarks/cutoff_switch/hcaii"))
from coordination_probe import ligand_mapping, error_components
from check_ammonium_states import ccd_vectors, handed


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def normalized_block(text, lines, num):
    # Only move coordinates. Reinstating the original tree reinstates the same
    # scored intraligand pairs and the same conformation-independent slope.
    fields = {
        int(l[6:11]): l[30:54] for l in lines if l.startswith(("ATOM ", "HETATM "))
    }
    block = []
    for l in text.splitlines():
        if l.startswith(("ATOM ", "HETATM ")):
            l = l[:30] + fields.pop(int(l[6:11])) + l[54:]
        block.append(l)
    if fields:
        raise ValueError("Unexpected atoms in returned pose")
    return f"MODEL {num}\n" + "\n".join(block) + "\nENDMDL\n"


def explicit_scores(binary, orig, receptor, bonds, text, models):
    """Score Cartesian poses under one prepared topology and explicit receptor."""
    with tempfile.TemporaryDirectory() as td:
        patched = Path(td) / "pool.pdbqt"
        patched.write_text(
            "".join(normalized_block(text, lines, k) for k, lines in models.items())
        )
        cmd = [
            binary,
            "--score_only",
            "--ligand",
            orig,
            "--ligand_bonds",
            bonds,
            "--receptor",
            receptor,
            "--rigid_window_pool",
            patched,
            "--rigid_window_samples",
            "16",
            "--verbosity",
            "0",
            "--cpu",
            "1",
        ]
        out = subprocess.check_output(cmd, text=True)
    return list(csv.DictReader(io.StringIO(out[out.index("model,inter,intra,") :])))


def evaluate(job):
    binary, base, folder, (panel, rec, lig, states, receptor, box, anchor) = job
    orig = states["orig"]
    text = orig.read_text()
    site = amides(text)[0] if panel != "fxa" or lig in FXA else None
    prep = orig.parent
    expected_angle = None
    if panel == "ca":
        meta = json.loads((prep / "metadata.json").read_text())[lig]
        label = meta["Ki_nM"]
        serials = ligand_mapping(prep, lig)[1]
        ref = meta["references"][rec]["coords"]
        sign_ids, reference_sign = ccd_vectors(prep.parent, lig, meta)
        expected_angle = angle(site, dict(zip(serials, ref, strict=True)))
    elif panel == "brd":
        meta = json.loads((prep / "manifest.json").read_text())[lig]
        label = meta["kd_uM"]
        serials = meta["heavy_serials"]
        refs = np.array(meta["alignments"][rec]["reference_by_match"])
        expected_angle = angle(site, dict(zip(serials, refs[0], strict=True)))
    else:
        with (ROOT / "benchmarks/affinity_partition/fxa_series.csv").open() as f:
            label = next(
                float(r["ki_nm"]) for r in csv.DictReader(f) if r["ligand"] == lig
            )
    allrows = []
    for state, seed, slot, form in itertools.product(
        states, SEEDS, (0, 1), (("free", "near", "far") if site else ("free",))
    ):
        pool = folder / f"{panel}_{rec}_{lig}_{state}_{seed}_{slot}_{form}.pdbqt"
        models = read_models(pool.read_text().splitlines())
        if not models:
            raise ValueError(f"No output poses in {pool}")
        bonds = folder / "inputs" / f"{panel}_{lig}_orig.bonds"
        ev = explicit_scores(binary, orig, receptor, bonds, text, models)
        for (k, lines), v in zip(models.items(), ev, strict=True):
            if int(v["model"]) != k:
                raise ValueError("Model number changed")
            atoms = parse_atoms(lines)
            points = {i: p.xyz for i, p in atoms.items()}
            phi = angle(site, points) if site else None
            rms = None
            hand = None
            if panel == "ca":
                rms = error_components(prep, lig, atoms, ref)["rmsd"]
                hand = int(
                    handed([np.array(points[i]) for i in sign_ids]) == reference_sign
                )
            elif panel == "brd":
                rms = np.sqrt(
                    np.square(refs - np.array([points[i] for i in serials]))
                    .sum(axis=2)
                    .mean(axis=1)
                ).min()
                hand = 1
            allrows.append(
                dict(
                    panel=panel,
                    rec=rec,
                    ligand=lig,
                    state=state,
                    seed=seed,
                    slot=slot,
                    form=form,
                    mode=k,
                    phi=phi,
                    crystal_phi=expected_angle,
                    rmsd=rms,
                    hand=hand,
                    measurement=label,
                    inter=float(v["inter"]),
                    intra=float(v["intra"]),
                    slope=float(v["slope"]),
                )
            )
    return allrows


def summarize(rows, out):
    groups = defaultdict(list)
    for r in rows:
        groups[r["panel"], r["rec"], r["ligand"], r["seed"]].append(r)
    cells = []
    for (p, rec, lig, seed), g in groups.items():
        unbound = min(
            (r for r in g if r["form"] == "free"), key=lambda r: r["inter"] + r["intra"]
        )["intra"]
        for arm in ("free", "near", "far", "both", "reverse"):

            def include(r):
                if p == "fxa" and lig not in FXA:
                    return r["form"] == "free"
                if arm in ("free", "near", "far"):
                    return r["form"] == arm
                return (
                    r["form"]
                    == (("near", "far") if arm == "both" else ("far", "near"))[
                        r["slot"]
                    ]
                )

            pool = [r for r in g if include(r)]
            for barrier in (0, 5, 15):

                def score(r):
                    return (r["inter"] + r["intra"] - unbound) * r["slope"] + (
                        barrier * np.sin(np.radians(r["phi"])) ** 2
                        if r["phi"] is not None
                        else 0.0
                    )

                top = min(pool, key=score)
                structural = top["rmsd"] is not None

                def hit(r):
                    return bool(r["rmsd"] <= 2 and r["hand"])

                cells.append(
                    dict(
                        panel=p,
                        receptor=rec,
                        ligand=lig,
                        seed=seed,
                        arm=arm,
                        barrier=barrier,
                        measurement=top["measurement"],
                        score=score(top),
                        state=top["state"],
                        slot=top["slot"],
                        mode=top["mode"],
                        phi=top["phi"],
                        rmsd=top["rmsd"],
                        strict=int(hit(top)) if structural else None,
                        oracle=int(any(hit(r) for r in pool)) if structural else None,
                        best_rmsd=min(r["rmsd"] for r in pool) if structural else None,
                    )
                )
    summ = []
    edges = []
    for p, arm, bar in sorted({(r["panel"], r["arm"], r["barrier"]) for r in cells}):
        g = [r for r in cells if (r["panel"], r["arm"], r["barrier"]) == (p, arm, bar)]
        prs = defaultdict(list)
        for block in sorted({(r["receptor"], r["seed"]) for r in g}):
            bs = sorted(
                (r for r in g if (r["receptor"], r["seed"]) == block),
                key=lambda r: r["measurement"],
            )
            for a, b in itertools.combinations(bs, 2):
                if b["measurement"] >= 3 * a["measurement"]:
                    gap = b["score"] - a["score"]
                    prs[a["ligand"], b["ligand"]].append(gap)
                    edges.append(
                        dict(
                            panel=p,
                            arm=arm,
                            barrier=bar,
                            receptor=block[0],
                            seed=block[1],
                            strong=a["ligand"],
                            weak=b["ligand"],
                            gap=gap,
                        )
                    )
        summ.append(
            dict(
                panel=p,
                arm=arm,
                barrier=bar,
                cells=len(g),
                top=sum(r["strict"] for r in g) if p != "fxa" else None,
                oracle=sum(r["oracle"] for r in g) if p != "fxa" else None,
                unique=sum((np.sign(np.median(v)) + 1) / 2 for v in prs.values()),
                pairs=len(prs),
            )
        )
    write(out / "cells.csv", cells)
    write(out / "pairs.csv", edges)
    write(out / "summary.csv", summ)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("binary", "base", "folder", "out"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    a = parser.parse_args()
    tasks = [(a.binary.resolve(), a.base, a.folder, case) for case in sources(a.base)]
    with ThreadPoolExecutor(a.jobs) as w:
        rows = list(itertools.chain.from_iterable(w.map(evaluate, tasks)))
    a.out.mkdir(exist_ok=True, parents=True)
    write(a.out / "modes.csv", rows)
    summarize(rows, a.out)
    print(len(rows), "common-topology model scores")
