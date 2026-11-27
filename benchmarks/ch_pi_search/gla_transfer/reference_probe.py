"""Check representability by rigidly placing the prepared conformer on the reference."""

import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from transfer import CASES, QUERY, save
from study import read_csv, timed_process, write_csv
from verify import geometry, point


def score_placed(binary, args, coords, work, log):
    timing = timed_process(
        [
            str(binary),
            *args,
            "--score_only",
            "--pose_coordinates",
            str(coords),
            "--unbound_energy",
            "0",
        ],
        work,
        log,
    )
    score = float(
        re.search(r"Estimated Free Energy of Binding +: +([-\d.]+)", log.read_text())[1]
    )
    return score, timing


def probe(work, output, binary):
    source = work / QUERY / "protonated.pdbqt"
    original = source.read_text().splitlines()
    atom_rows = [i for i, line in enumerate(original) if line.startswith("ATOM  ")]
    serials = [int(original[i][6:11]) for i in atom_rows]
    points = np.array([point(original[i]) for i in atom_rows])
    stem = output / "probe"
    stem.mkdir(exist_ok=True)
    results = []
    meta = json.loads((work / "metadata.json").read_text())
    compared = []
    for receptor in CASES:
        refs, order = geometry(work, "GLA", QUERY, receptor)
        (reference,) = refs
        idx = [serials.index(n) for n in order]
        mobile = points[idx]
        rotation = Rotation.align_vectors(
            reference - reference.mean(0), mobile - mobile.mean(0)
        )[0]
        placed = rotation.apply(points - mobile.mean(0)) + reference.mean(0)
        lines = original.copy()
        for i, xyz in zip(atom_rows, placed):
            lines[i] = lines[i][:30] + "".join(f"{x:8.3f}" for x in xyz) + lines[i][54:]
        pose = stem / (receptor + "_placed.pdbqt")
        pose.write_text("\n".join(lines) + "\n")
        args = [
            "--cpu",
            "1",
            "--seed",
            "43103",
            "--ligand",
            str(pose),
            "--receptor",
            str(work / receptor / "receptor.pdbqt"),
        ]
        for axis, center in zip("xyz", meta[receptor]["center"]):
            args += ["--center_" + axis, str(center), "--size_" + axis, "24"]
        local = stem / (receptor + "_relaxed.pdbqt")
        timing = timed_process(
            [str(binary), *args, "--local_only", "--out", str(local)],
            work,
            stem / (receptor + "_relax.log"),
        )
        stage_scores = {}
        for stage, coordinates in (("placed", pose), ("relaxed", local)):
            atoms = {
                int(line[6:11]): point(line)
                for line in coordinates.read_text().splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            }
            rmsd = np.sqrt(
                np.square([np.array(atoms[n]) - r for n, r in zip(order, reference)])
                .sum(1)
                .mean()
            )
            log = stem / f"{receptor}_{stage}_score.log"
            score, score_time = score_placed(binary, args, coordinates, work, log)
            stage_scores[stage] = score
            results.append(
                dict(
                    receptor=receptor,
                    stage=stage,
                    rmsd_A=rmsd,
                    zero_reference_score=score,
                    score_cost=score_time,
                    relaxation_cost=timing if stage == "relaxed" else None,
                )
            )
        for row in read_csv(output / "costs.csv"):
            if row["receptor"] != receptor:
                continue
            tag = "_".join(
                str(row[k])
                for k in ("panel", "query", "receptor", "seed", "cap", "search")
            )
            pool = (output / "pools" / (tag + ".pdbqt")).read_text()
            first = stem / (tag + "_first.pdbqt")
            first.write_text(pool.split("ENDMDL")[0] + "ENDMDL\n")
            log = stem / (tag + "_first.log")
            dry, score_time = score_placed(binary, args, first, work, log)
            compared.append(
                dict(
                    receptor=receptor,
                    seed=row["seed"],
                    search=row["search"],
                    predicted_first_zero=dry,
                    reference_local_zero=stage_scores["relaxed"],
                    predicted_minus_reference=dry - stage_scores["relaxed"],
                    wall_s=score_time["wall_s"],
                )
            )
    write_csv(output / "probe_competition.csv", compared)
    save(output / "reference_probe.json", results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "output", "binary"):
        p.add_argument(name, type=lambda x: Path(x).resolve())
    a = p.parse_args()
    probe(a.work, a.output, a.binary)
