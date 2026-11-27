"""Offline AD4 score/assay audit and read-only localization from the FKBP12 coordinates."""

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from place_assay import OLD, RECEIVERS, rows, table

sys.path.insert(0, str(OLD))
import evaluate as historical


def coordinates(block, serial):
    return historical.xyz(block, list(serial))


def rms(x, references):
    return np.min(
        np.sqrt(
            np.mean(np.sum((x[..., None, :, :] - references) ** 2, axis=-1), axis=-1)
        ),
        axis=-1,
    )


def check(work, out):
    historical.inspect(
        work,
        out,
        arms=("official", "base", "high", "match", "shuffle"),
        comparison="match",
    )
    refs = {
        r: historical.reference(work, r, c)
        for r, c in zip(RECEIVERS, ("SB3", "SBX"), strict=True)
    }
    tracks, coordinates_audit, proximity = [], [], []
    descriptors = rows(out / "descriptors.csv")
    file_scores = rows(out / "scores.csv")
    scores = {
        (r["arm"], int(r["compound"])): float(r["top1_score"])
        for r in file_scores
        if r["receiver"] == "mean"
    }
    contrast = []
    for arm, ref in (("match", "base"), ("shuffle", "base"), ("match", "shuffle")):
        changes = [
            scores[arm, int(d["compound"])] - scores[ref, int(d["compound"])]
            for d in descriptors
        ]
        contrast.append(
            dict(
                arm=arm,
                ref=ref,
                n=len(changes),
                rho_heavy_atoms=spearmanr(
                    [float(d["heavy_atoms"]) for d in descriptors], changes
                ).statistic,
                rho_torsdof=spearmanr(
                    [float(d["torsdof"]) for d in descriptors], changes
                ).statistic,
                distinct_formal_charges=len({d["formal_charge"] for d in descriptors}),
            )
        )
    table(out / "bias_contrasts.csv", contrast)
    for rec in RECEIVERS:
        atoms = [
            l
            for l in (work / rec / "receptor.pdbqt").read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM"))
            and l.split()[-1] not in ("HD", "H", "HS")
        ]
        tree = cKDTree([[float(l[i : i + 8]) for i in (30, 38, 46)] for l in atoms])
        for n in range(4, 12):
            ligand = (work / "affinity" / f"{n}.pdbqt").read_text()
            mol, serials, _ = historical.molecule(ligand)
            stem = f"{rec}_{n}"
            crystal = {8: "1FKG", 9: "1FKH"}.get(n)
            with np.load(work / "plans" / (stem + ".npz")) as data:
                order = [data["serials"].tolist().index(i) for i in serials]
                ref = (
                    historical.crystal_positions(refs[crystal], refs[rec], mol)
                    if crystal
                    else None
                )
                for arm, key in (("match", "matched"), ("shuffle", "sham")):
                    p = work / "pools" / arm / stem
                    blocks = [
                        b
                        for b in p.with_suffix(".pdbqt").read_text().split("ENDMDL")
                        if b.strip()
                    ]
                    released = np.array([coordinates(b, serials) for b in blocks])
                    distances = tree.query(released)[0].min(axis=1)
                    proximity.extend(
                        dict(receiver=rec, compound=n, arm=arm, rank=k, minimum_A=d)
                        for k, d in enumerate(distances, 1)
                    )
                    if not crystal:
                        continue
                    local = data["xyz"] - data["xyz"].mean(0)
                    candidates = (
                        local @ data[key + "_rotation"].transpose(0, 2, 1)
                        + data["center"][:, None, :]
                    )
                    distances = rms(candidates[:, order], ref)
                    selected = int(data[key + "_best"])
                    lines = p.with_suffix(".xyz").read_text().splitlines()
                    assert list(map(int, lines[0][1:].split())) == list(data["serials"])
                    hist = [
                        list(map(int, r.split()[4:]))
                        for r in p.with_suffix(".history").read_text().splitlines()
                    ]
                    steps = [
                        int(r.split()[0])
                        for r in p.with_suffix(".history").read_text().splitlines()
                    ]
                    events = {}
                    for line in lines[1:]:
                        tag, tr, step, *points = line.split()
                        x = np.array(points, float).reshape(-1, 3)[order]
                        coordinates_audit.append(
                            dict(
                                receiver=rec,
                                compound=n,
                                arm=arm,
                                stage=tag,
                                scope=tr,
                                step_or_rank=step,
                                rmsd_A=float(rms(x, ref)),
                            )
                        )
                        if tag in ("selected", "hunt", "second"):
                            events[tag] = x
                    np.testing.assert_allclose(
                        events["selected"],
                        candidates[selected, order],
                        rtol=0,
                        atol=1e-10,
                    )
                    assert 1 in hist[0] and steps[0] == 1
                    # Output coordinates are rounded to 0.001 A. AD4 does no late explicit refinement.
                    same = (
                        np.max(abs(released - events["second"]), axis=(1, 2)) <= 0.00051
                    )
                    rank = np.flatnonzero(same) + 1
                    tracks.append(
                        dict(
                            receiver=rec,
                            compound=n,
                            arm=arm,
                            candidates=len(distances),
                            source_A=float(rms(data["xyz"][order], ref)),
                            best_generated_A=min(distances),
                            generated_lt2=int(sum(distances < 2)),
                            generated_lt3=int(sum(distances < 3)),
                            winning_index=selected,
                            selected_A=distances[selected],
                            hunt_A=float(rms(events["hunt"], ref)),
                            second_A=float(rms(events["second"], ref)),
                            insert_step=1,
                            eviction_step=next(
                                (s for s, h in zip(steps, hist) if 1 not in h), ""
                            ),
                            preserved_exact_at_output_rank=" ".join(map(str, rank)),
                        )
                    )
    table(out / "stages.csv", tracks)
    table(out / "stage_coordinates.csv", coordinates_audit)
    table(out / "new_receptor_contact.csv", proximity)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    check(a.work.resolve(), a.out.resolve())
