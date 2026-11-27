"""Diagnose changes in receptor-frame error on captured coordinates; no fitting of scores."""

import argparse
import collections
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tarfile
import tempfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
PANELS = ("rigid_pose_coverage", "fragment_anchor")
EPS = 1e-8
PARTS = ("centroid", "shape", "orientation")


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def by(rows, key):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[row[key]].append(row)
    return groups


def manifest(path):
    return {name: digest for digest, name in
            (line.split(maxsplit=1) for line in path.read_text().splitlines())}


def verify(path, expected):
    with path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != expected:
            raise ValueError(f"SHA-256 mismatch: {path}")


def unpack(parent, work):
    # The two parent experiments publish both compressed hashes and member hashes.
    archives = manifest(parent / "ARCHIVE.sha256")
    for name, expected in archives.items():
        verify(parent / name, expected)
        with tarfile.open(parent / name) as archive:
            archive.extractall(work, filter="data")
    members = manifest(parent / "results/raw_files.sha256")
    for name, expected in members.items():
        verify(work / name, expected)
    return len(archives), len(members)


def reference_reader(panel, parent, work):
    # Reuse the registered receptor alignment, CCD mappings and microstate equivalences.
    evaluator = BENCH / "cavity_field_rank/pose/evaluate.py"
    spec = importlib.util.spec_from_file_location("previous_reference", evaluator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cells = read_csv(parent / "results/cells.csv")
    module.panel.CASES = {r["drug"]: (r["query"], r["receiver"], ()) for r in cells}
    if panel == "rigid_pose_coverage":
        module.COMPONENTS = {"3ANR": "HRM", "7AKL": "RK5", "3FMZ": "2T1", "1RBP": "RTL"}
    else:
        module.COMPONENTS.update({"2NNQ": "T4B", "3FR2": "8CA"})
    historical = json.loads((parent / "results/references.json").read_text())
    if "archives" in historical:
        for name, expected in historical["archives"].items():
            verify(BENCH.parent / name, expected)
    else:
        verify(BENCH / "cubic_grid/human_panel/depositions.tar.xz", historical["raw_sha256"])
    prepared = module.reference_inputs(work)
    for name, (cell, *_) in prepared.items():
        saved = historical["serials"][name]
        assert tuple(cell["serials"]) == tuple(saved["serials"])
        assert len(cell["refs"]) == saved.get("chemical_mappings", saved.get("equivalent_frames"))
    return {name: v[0] for name, v in prepared.items()}, module.xyz


def horn_fit_squared(cross, xnorm, ynorm, n):
    # Maximum over SO(3) as the largest eigenvalue of Horn's quadratic form.
    tr = np.trace(cross, axis1=-2, axis2=-1)
    antisym = cross - cross.transpose(0, 2, 1)
    z = antisym[:, [1, 2, 0], [2, 0, 1]]
    h = np.zeros((len(cross), 4, 4))
    h[:, 0, 0] = tr
    h[:, 1:, 0] = h[:, 0, 1:] = z
    h[:, 1:, 1:] = cross + cross.transpose(0, 2, 1) - tr[:, None, None] * np.eye(3)
    return (xnorm + ynorm - 2 * np.linalg.eigvalsh(h)[:, -1]) / n


def components(x, refs):
    """Proper-rotation Procrustes floor and squared receptor-frame identity (Å, Å²)."""
    refs = np.asarray(refs)
    if not np.isfinite(x).all() or not np.isfinite(refs).all():
        raise ValueError("Nonfinite coordinates")
    ybar = refs.mean(axis=1)
    np.testing.assert_allclose(ybar, np.broadcast_to(ybar[0], ybar.shape), rtol=0, atol=1e-10)
    xbar = x.mean(axis=0)
    local, targets = x - xbar, refs - ybar[:, None]
    d2 = np.square(refs - x).sum(axis=-1).mean(axis=-1)
    t2 = np.sum((xbar - ybar[0]) ** 2)
    c2 = np.square(targets - local).sum(axis=-1).mean(axis=-1)
    np.testing.assert_allclose(d2, c2 + t2, rtol=0, atol=1e-10)
    cross = local.T @ targets
    u, _, vt = np.linalg.svd(cross)
    u[:, :, -1] *= np.linalg.det(u @ vt)[:, None]
    rotations = u @ vt
    fitted = np.square(local @ rotations - targets).sum(axis=-1).mean(axis=-1)
    check = horn_fit_squared(cross, np.square(local).sum(), np.square(targets).sum(axis=(1, 2)), len(x))
    np.testing.assert_allclose(fitted, check, rtol=0, atol=1e-10)
    direct = int(d2.argmin())
    optimal = int(fitted.argmin())
    orientation = c2[direct] - fitted[optimal]
    if orientation < -EPS:
        raise ValueError("Negative orientation residual")
    return dict(rmsd_A=np.sqrt(d2[direct]), centroid_A=np.sqrt(t2),
                shape_A=np.sqrt(fitted[optimal]), orientation_A=np.sqrt(max(orientation, 0)),
                direct_map=direct, shape_map=optimal, shape_direct_map_A=np.sqrt(fitted[direct]),
                map_release_A2=fitted[direct] - fitted[optimal],
                rotation_direct_map_A2=max(c2[direct] - fitted[direct], 0),
                horn_max_abs_error_A2=np.max(np.abs(fitted - check)))


def check_synthetic():
    x = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]], dtype=float)
    y = np.array([x])
    rot = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    for a in (x, x + [3, 4, 0], (x - x.mean(0)) @ rot + x.mean(0)):
        r = components(a, y)
        assert r["shape_A"] < 1e-12
        assert abs(r["rmsd_A"]**2 - sum(r[k+"_A"]**2 for k in PARTS)) < EPS
    assert abs(components(x + [3, 4, 0], y)["centroid_A"] - 5) < 1e-12
    assert components(x * [-1, 1, 1], y)["shape_A"] > 0.1  # a reflection is forbidden
    y = np.stack([x, x[[0, 2, 1, 3]]])
    assert components(y[1], y)["rmsd_A"] == 0


def difference(stem, contrast, old, new):
    row = dict(stem=stem, contrast=contrast, from_id=old["id"], to_id=new["id"],
               from_A=old["rmsd_A"], to_A=new["rmsd_A"],
               delta_total_A2=new["rmsd_A"]**2 - old["rmsd_A"]**2)
    deltas = {p: new[p+"_A"]**2 - old[p+"_A"]**2 for p in PARTS}
    if abs(row["delta_total_A2"] - sum(deltas.values())) > EPS:
        raise ValueError(f"Squared error identity failed in {stem}")
    return row | {"delta_"+p+"_A2": d for p, d in deltas.items()} | dict(
        mixed=int(max(deltas.values()) > EPS and min(deltas.values()) < -EPS),
        flat=int(max(map(abs, deltas.values())) <= EPS))


def published_modes(panel, parent, work):
    previous_modes = read_csv(parent / "results/modes.csv")
    files = {r["stem"]: r["stem"] + ("-count" if r["arm"] == "pro-cost" else "")
             for r in previous_modes}
    if panel == "rigid_pose_coverage":
        selected = {r["case"]: r for r in read_csv(work / "extra_trajectory_selected.csv")}
        for r in read_csv(parent / "results/extra_trajectory_modes.csv"):
            drug, seed = r["case"].split("_")
            stem = r["case"] + "_pro-extra"
            files[stem] = selected[r["case"]]["selected_stem"]
            previous_modes.append(dict(stem=stem, drug=drug, seed=seed, arm="pro-extra",
                                       rank=r["rank"], rmsd_A=r["rmsd_A"]))
    return previous_modes, files


def panel_rows(panel, work, out):
    parent = BENCH / panel / "pose"
    n_archives, n_members = unpack(parent, work)
    references, xyz = reference_reader(panel, parent, work)
    previous_modes, files = published_modes(panel, parent, work)

    modes, stages, pairs, pools = [], [], [], []
    stage_records = [r for r in read_csv(parent / "results/stages.csv")
                     if r["scope"] == "1" and r["index"] == "1"
                     and r["stage"] in ("selected", "hunt", "second")]
    by_stage = by(stage_records, "stem")

    def measure(key, x, cell, label, parent_rmsd):
        r = components(x, cell["refs"])
        err = abs(r["rmsd_A"] - float(parent_rmsd))
        if err >= 1e-7:
            raise ValueError(f"Historical RMSD mismatch: {key}, {label}, {err}")
        return key | dict(id=label, **r, parent_abs_error_A=err)

    for stem, records in by(previous_modes, "stem").items():
        key = {c: records[0][c] for c in ("stem", "drug", "seed", "arm")}
        cell = references[key["drug"]]
        published = (work / "raw" / (files[stem]+".pdbqt")).read_text().split("ENDMDL")[:-1]
        actual = []
        for block, old in zip(published, records, strict=True):
            actual.append(measure(key, xyz(block, cell["serials"]), cell, int(old["rank"]), old["rmsd_A"]))
        modes.extend(actual)
        best = min(actual, key=lambda row: (row["rmsd_A"], row["id"]))
        shape = min(actual, key=lambda row: (row["shape_A"], row["id"]))
        pairs.append(difference(stem, "top_to_oracle", actual[0], best))
        pools.append(key | dict(n=len(actual), top_A=actual[0]["rmsd_A"],
                                best5_A=min(r["rmsd_A"] for r in actual[:5]),
                                oracle_rank=best["id"], oracle_A=best["rmsd_A"],
                                regret_A=actual[0]["rmsd_A"]-best["rmsd_A"],
                                best_shape_rank=shape["id"], best_shape_A=shape["shape_A"],
                                best_shape_actual_A=shape["rmsd_A"],
                                top_shape_A=actual[0]["shape_A"], oracle_shape_A=best["shape_A"]))
        if stem not in by_stage:
            continue
        lines = (work / "dump" / stem / "stages.xyz").read_text().splitlines()
        ids = list(map(int, lines[0].split()[1:]))
        order = [ids.index(s) for s in cell["serials"]]
        checks = {r["stage"]: r["rmsd_A"] for r in by_stage[stem]}
        measured = {}
        for line in lines[1:]:
            tag, scope, index, values = line.split(maxsplit=3)
            if tag in checks and scope == "1" and index == "1":
                x = np.fromstring(values, sep=" ").reshape(len(ids), 3)[order]
                measured[tag] = measure(key, x, cell, tag, checks[tag])
        if measured.keys() != checks.keys():
            raise ValueError(f"Missing stage in {stem}")
        stages.extend(measured.values())
        pairs.append(difference(stem, "selected_to_hunt", measured["selected"], measured["hunt"]))
        if "second" in measured:
            pairs.append(difference(stem, "hunt_to_second", measured["hunt"], measured["second"]))

    # Keep missing or extra modes/stages explicit, including any historical empty cells.
    cell_rows = read_csv(parent / "results/cells.csv")
    if any(int(c["n"]) == 0 for c in cell_rows):
        raise ValueError("Unexpected empty historical pool; report its missing geometry")
    expected = len(cell_rows) + (sum(pool["arm"] == "pro-extra" for pool in pools))
    assert len(pools) == expected and len(stages) == len(stage_records)
    for label, rows in (("modes", modes), ("stages", stages), ("paired", pairs), ("pools", pools)):
        table(out / (panel+"_"+label+".csv"), rows)

    summaries = []
    for contrast, group in by(pairs, "contrast").items():
        delta = np.array([r["delta_total_A2"] for r in group])
        summaries.append(dict(panel=panel, contrast=contrast, n=len(group),
                              improve=int((delta < -EPS).sum()), worse=int((delta > EPS).sum()),
                              same=int((abs(delta) <= EPS).sum()), mixed=sum(r["mixed"] for r in group),
                              all_components_fixed=sum(r["flat"] for r in group)))
    return summaries, dict(panel=panel, archives=n_archives, archive_members=n_members,
            pools=len(pools), modes=len(modes), stages=len(stages),
            max_parent_RMSD_error_A=max(r["parent_abs_error_A"] for r in modes+stages),
            max_Horn_SVD_difference_A2=max(r["horn_max_abs_error_A2"] for r in modes+stages))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("out", type=Path)
    args = cli.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    check_synthetic()
    checks, summaries = [], []
    for panel in PANELS:
        with tempfile.TemporaryDirectory(prefix="pose_geometry_") as tmp:
            tic, cpu = time.perf_counter(), time.process_time()
            summary, receipt = panel_rows(panel, Path(tmp), args.out)
            summaries.extend(summary)
            checks.append(receipt | dict(wall_s=time.perf_counter()-tic, cpu_s=time.process_time()-cpu))
    table(args.out/"summary.csv", summaries)
    table(args.out/"verification.csv", checks)


if __name__ == "__main__":
    main()
