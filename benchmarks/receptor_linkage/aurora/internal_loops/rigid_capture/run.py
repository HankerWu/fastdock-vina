"""Crystal-informed common-state rigid transforms and full-production continuations, diagnostic only."""

import argparse
import csv
import hashlib
import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_rows(path, rows):
    with path.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def rms(x, refs):
    return np.sqrt(np.mean(np.sum((x - refs) ** 2, axis=-1), axis=-1)).min()


def best_rotation(x, refs):
    u, _, vt = np.linalg.svd(
        (x - x.mean(axis=0)).T @ (refs - refs.mean(axis=1, keepdims=True))
    )
    u[:, :, -1] *= np.linalg.det(u @ vt)[:, None]
    row = u @ vt
    residual = np.mean(
        np.sum(
            ((x - x.mean(axis=0)) @ row - refs + refs.mean(axis=1, keepdims=True)) ** 2,
            axis=-1,
        ),
        axis=-1,
    )
    return Rotation.from_matrix(row[np.argmin(residual)].T)


def placements(x, refs):
    c, goal = x.mean(axis=0), refs[0].mean(axis=0)
    rotate = best_rotation(x, refs)
    yield "native", Rotation.identity(), c
    yield "center_only", Rotation.identity(), goal
    yield "orientation_only", rotate, c
    yield "fit", rotate, goal
    for axis, unit in zip("xyz", np.eye(3), strict=True):
        for sign, suffix in ((-1, "-"), (1, "+")):
            extra = Rotation.from_rotvec(unit * np.deg2rad(sign * 30))
            yield f"shift_{axis}{suffix}", rotate, goal + sign * 2 * unit
            yield f"turn_{axis}{suffix}", extra * rotate, goal
            yield f"both_{axis}{suffix}", extra * rotate, goal + sign * 2 * unit


def transformed(q, x, rot, center, tag):
    if tag == "native":
        return q.copy(), x.copy()
    centroid = x.mean(axis=0)
    origin = rot.apply(q[:3] - centroid) + center
    qrot = rot * Rotation.from_quat(q[[4, 5, 6, 3]])
    return np.r_[origin, qrot.as_quat()[[3, 0, 1, 2]], q[7:]], rot.apply(
        x - centroid
    ) + center


def main(work, binary, out):
    out.mkdir(exist_ok=False, parents=True)
    raw = out / "raw"
    raw.mkdir()
    ligand, bonds, query_order, refs = replay.parameters(work, "at9283")
    receptor = work / replay.PROTEINS["1OL7_released"]
    atoms = [
        line
        for line in receptor.read_text().splitlines()
        if line.startswith(("ATOM ", "HETATM"))
        and line.split()[-1] not in ("HD", "H", "HS")
    ]
    protein = np.array(
        [[float(line[k : k + 8]) for k in (30, 38, 46)] for line in atoms]
    )
    box = [
        s
        for axis, value in zip(
            "xyz",
            json.loads((work / "prepared/2X6E.json").read_text())["center"],
            strict=True,
        )
        for s in (f"--center_{axis}", str(value), f"--size_{axis}", "26")
    ]
    common = [
        str(binary),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--ligand_bonds",
        str(bonds),
        *box,
        "--cpu",
        "1",
        "--exhaustiveness",
        "2",
        "--seed",
        "61231",
        "--max_evals",
        "25000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--pose_selection",
        "late",
        "--search_counters",
    ]
    costs = []

    def run(stem, variables):
        p = raw / stem
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        env = {
            k: v
            for k, v in os.environ.items()
            if not k.startswith(("FD_RIGID_", "FD_DUMP_", "FD_QUENCH_"))
        }
        with p.with_suffix(".log").open("w") as log:
            subprocess.run(
                [*common, "--out", str(p.with_suffix(".pdbqt"))],
                env=env | {k: str(v) for k, v in variables.items()},
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        costs.append(
            dict(
                process=stem,
                wall_s=time.perf_counter() - start,
                cpu_s=after.ru_utime
                + after.ru_stime
                - before.ru_utime
                - before.ru_stime,
                main_evals=sum(
                    map(
                        int,
                        re.findall(
                            r"\b(?:hunt|second|explicit)_evals=(\d+)",
                            p.with_suffix(".log").read_text(),
                        ),
                    )
                ),
            )
        )

    base = HERE.parent / "proposal_crossover/results"
    input_paths = []
    metadata, general, intended, original, expected, end_native = [], [], [], [], [], []
    for source in ("complete", "inferred"):
        path = base / (source + "_on_complete.npz")
        input_paths.append(path)
        with np.load(path) as data:
            serials = data["serials"].tolist()
            ordered_refs = refs[:, [query_order.index(s) for s in serials]]
            for i, (tr, step) in enumerate(
                zip(data["trajectory"], data["step"], strict=True)
            ):
                if source == "complete" and int(step) != 1 and step % 100:
                    continue
                if source == "inferred" and (tr, step) != (2, 706):
                    continue
                x, q = data["start_xyz"][i], data["generalized"][i]
                for tag, rotation, destination in placements(x, ordered_refs):
                    g, cart = transformed(q, x, rotation, destination, tag)
                    metadata.append(
                        dict(
                            source=source,
                            trajectory=int(tr),
                            step=int(step),
                            placement=tag,
                        )
                    )
                    general.append(g)
                    intended.append(cart)
                    original.append(q)
                    end_native.append(data["return_xyz"][i])
                    expected.append(float(data["energy"][i]))

    qpath, qout = raw / "quench.q", raw / "quench.xyz"
    seq = np.arange(len(metadata))
    np.savetxt(
        qpath,
        np.c_[np.ones_like(seq), seq + 1, general],
        fmt=["%d", "%d"] + ["%.17g"] * len(general[0]),
    )
    run("quench", {"FD_QUENCH_INPUT": qpath, "FD_QUENCH_OUTPUT": qout})
    with qout.open() as f:
        assert list(map(int, f.readline()[1:].split())) == serials
    record = np.loadtxt(qout)
    assert np.array_equal(record[:, 1], seq + 1)
    n = 3 * len(serials)
    entry = record[:, 3 : 3 + n].reshape(-1, len(serials), 3)
    finals = record[:, 5 + n :].reshape(-1, len(serials), 3)
    evals = record[:, 3 + n].astype(int)
    np.testing.assert_allclose(entry, intended, atol=1e-8, rtol=0)
    for i, row in enumerate(metadata):
        if row["placement"] == "native":
            np.testing.assert_array_equal(finals[i], end_native[i])
            assert record[i, 4 + n] == expected[i]
    result = []
    for i, row in enumerate(metadata):
        values = {}
        for tag, xyz in (("start", entry[i]), ("return", finals[i])):
            d = cdist(xyz, protein)
            values.update(
                {
                    f"{tag}_A": rms(xyz, ordered_refs),
                    f"{tag}_nearest_heavy_A": d.min(),
                    f"{tag}_heavy_pairs_lt_2A": int((d < 2).sum()),
                }
            )
        result.append(
            dict(
                **row,
                **values,
                entry_score=record[i, 2],
                return_score=record[i, 4 + n],
                hunt_evals=int(evals[i]),
            )
        )
    write_rows(out / "quench.csv", result)
    np.savez_compressed(
        out / "quench.npz",
        start_xyz=entry,
        return_xyz=finals,
        serials=serials,
        generalized=general,
        entry_score=record[:, 2],
        return_score=record[:, 4 + n],
        hunt_evals=evals,
    )
    # The compact archive contains the same 17-digit round-tripped binary64 values.
    qout.unlink()
    qpath.unlink()
    branches, modes, stage_rows = [], [], []
    for i, meta in enumerate(metadata):
        if meta["source"] != "complete" or meta["placement"] != "fit":
            continue
        tr, step = meta["trajectory"], meta["step"]
        stem = f"fit_t{tr}_s{step}"
        spec, log = raw / (stem + ".input"), raw / (stem + ".stages")
        spec.write_text(
            " ".join([str(tr), str(step), *[f"{v:.17g}" for v in general[i][:7]]])
            + "\n"
        )
        stagepath = raw / (stem + ".xyz")
        run(
            stem,
            {
                "FD_RIGID_INTERVENTION": spec,
                "FD_RIGID_LOG": log,
                "FD_RIGID_COORDS": stagepath,
            },
        )
        stage_lines = stagepath.read_text().splitlines()
        assert list(map(int, stage_lines[0][1:].split())) == serials
        stage_minima = {}
        for line in stage_lines[1:]:
            stage, scope, index, *points = line.split()
            xyz = np.array(points, float).reshape(len(serials), 3)
            err = rms(xyz, ordered_refs)
            stage_rows.append(
                dict(
                    trajectory=tr,
                    step=step,
                    stage=stage,
                    scope=int(scope),
                    index=int(index),
                    rmsd_A=err,
                )
            )
            if stage != "end_trajectory" or int(scope) == tr:
                stage_minima[stage] = min(err, stage_minima.get(stage, np.inf))
        lines = log.read_text().splitlines()
        np.testing.assert_array_equal(
            np.array(lines[0].split()[1:], float), original[i]
        )
        assert float(lines[1].split()[1]) == result[i]["return_score"]
        accepted = "accepted" in lines
        second = next((v for v in lines if v.startswith("second ")), None)
        pdbqt = raw / (stem + ".pdbqt")
        errors = [rms(p, refs) for p in replay.geom.mode_coords(pdbqt, query_order)]
        energies = [
            float(line.split()[3])
            for line in pdbqt.read_text().splitlines()
            if line.startswith("REMARK VINA RESULT:")
        ]
        for rank, (err, score) in enumerate(zip(errors, energies, strict=True), 1):
            modes.append(
                dict(trajectory=tr, step=step, rank=rank, rmsd_A=err, score=score)
            )
        branches.append(
            dict(
                trajectory=tr,
                step=step,
                start_A=result[i]["start_A"],
                after_hunt_A=result[i]["return_A"],
                accepted=int(accepted),
                second=int(second is not None),
                inserted=int(second.split()[-1]) if second else 0,
                **{
                    f"{stage}_min_A": stage_minima.get(stage, "")
                    for stage in (
                        "second",
                        "end_trajectory",
                        "pooled",
                        "fin_input",
                        "fin_output",
                    )
                },
                n=len(errors),
                top1_A=errors[0],
                best5_A=min(errors[:5]),
                oracle_A=min(errors),
                regret_A=errors[0] - min(errors),
                hash=digest(pdbqt),
            )
        )
    write_rows(out / "branches.csv", branches)
    write_rows(out / "modes.csv", modes)
    write_rows(out / "stages.csv", stage_rows)
    write_rows(out / "costs.csv", costs)
    manifest = {
        k: digest(p)
        for k, p in {
            "binary": binary,
            "receptor": receptor,
            "ligand": ligand,
            "bonds": bonds,
            "box": work / "prepared/2X6E.json",
            "ligand_sdf": work / "prepared/at9283.sdf",
            "queries": work / "2W1G.pdb",
            "receiver": work / "raw/2X6E.pdb",
            "ccd": work / "L0G.cif",
            "receiver_archive": replay.ARCHIVE,
            **{p.name: p for p in input_paths},
        }.items()
    }
    (out / "manifest.json").write_text(
        json.dumps(
            dict(
                sha256=manifest,
                quench_rows=len(evals),
                entry_observer_evals=len(evals),
                extra_hunt_evals=int(evals.sum()),
                max_start_cartesian_error_A=float(np.max(np.abs(entry - intended))),
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "binary", "out"):
        p.add_argument(k, type=Path)
    main(*(v.resolve() for v in vars(p.parse_args()).values()))
