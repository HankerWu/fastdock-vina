"""Compare equal-node interpolants of the native receptor field near observed ligand atoms."""
import argparse
import csv
import gzip
import io
import json
import subprocess
from pathlib import Path

import numpy as np


PDB_IDS = ("1VJ5", "2B8L", "2E9P", "2YWP", "2ZJM", "3ANS", "6Y7L", "6YQW")
CAP = 10.0


def weights(t, cubic):
    if not cubic:
        return np.array([0., 1-t, t, 0.]), np.array([0., -1., 1., 0.])
    # Hermite end values and centered (Catmull–Rom) end slopes.
    tangents = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [-.5, 0, .5, 0], [0, -.5, 0, .5]])
    h = np.array([1-3*t*t+2*t**3, 3*t*t-2*t**3, t-2*t*t+t**3, -t*t+t**3])
    dh = np.array([-6*t+6*t*t, 6*t-6*t*t, 1-4*t+3*t*t, -2*t+3*t*t])
    return h @ tangents, dh @ tangents


def interp(nodes, fraction, cubic, spacing):
    w, dw = zip(*(weights(t, cubic) for t in fraction))

    def tensor(vectors):
        return float(np.einsum("i,j,k,ijk", *vectors, nodes))

    energy = tensor(w)
    grad = np.array([tensor([dw[j] if j == axis else w[j] for j in range(3)])
                     for axis in range(3)]) / spacing
    if energy > 0:
        saturation = CAP / (CAP + energy)
        return energy * saturation, grad * saturation**2, energy
    return energy, grad, energy


def run(probe, inputs, pdb_index, spacing):
    meta = json.loads((inputs / "preparation.json").read_text())
    rng = np.random.default_rng(260126 + pdb_index)
    origin = np.asarray(meta["center"]) - 12
    queries = []
    for serial, xyz in zip(meta["serial"], meta["mapped"][0]):
        for radius in (0., .25, .75, 1.5):
            displacement = np.zeros(3) if radius == 0 else rng.normal(size=3)
            if radius:
                displacement *= radius / np.linalg.norm(displacement)
            point = np.asarray(xyz) + displacement
            # The production cubic boundary rule is irrelevant to these interior probes.
            index = np.floor((point-origin)/spacing)
            if (index < 1).any() or (index > round(24/spacing)-2).any():
                raise ValueError("Probe outside the full interior stencil")
            queries.append((serial, radius, point))
    def record(vals):
        return " ".join(map(str, vals)) + "\n"
    stdin = record([*origin, spacing]) + "".join(record([s, *q]) for s, _, q in queries)
    output = subprocess.run([probe, inputs/"receptor.pdbqt", inputs/"ligand.pdbqt"],
                            input=stdin, text=True, capture_output=True, check=True)
    values = np.loadtxt(io.StringIO(output.stdout), ndmin=2)
    if values.shape != (len(queries), 79):
        raise ValueError(f"Wrong probe output shape {values.shape}")
    rows = []
    for (serial, radius, q), raw in zip(queries, values):
        nodes = raw[15:].reshape((4, 4, 4))
        fraction = np.mod((q-origin)/spacing, 1)
        linear, lg, _ = interp(nodes, fraction, False, spacing)
        cubic, cg, cubic_raw = interp(nodes, fraction, True, spacing)
        if max(abs(linear - raw[7]), max(abs(lg - raw[8:11]))) > 2e-11:
            raise ValueError("Independent trilinear arithmetic differs from native grid")
        if max(abs(cubic - raw[11]), max(abs(cg - raw[12:15]))) > 2e-11:
            raise ValueError("Independent Hermite arithmetic differs from engine cubic grid")
        f = nodes[1:3, 1:3, 1:3]
        row = dict(target=meta["target"], pdb=meta["pdb"], serial=serial, radius_A=radius,
                   x=q[0], y=q[1], z=q[2], xs=int(raw[0]), cutoff_contact=int(raw[1]),
                   explicit_fast=raw[2], explicit_smooth=raw[3], linear=linear, cubic=cubic,
                   signed_overshoot=max(0., cubic_raw-f.max()) - max(0., f.min()-cubic_raw))
        for axis, ref, lin, cub in zip("xyz", raw[4:7], lg, cg):
            for k, v in (("ref", ref), ("linear", lin), ("cubic", cub)):
                row[k+"_g"+axis] = v
        rows.append(row)
    return rows


def summarise(rows):
    linear = np.array([r["linear"] for r in rows])
    cubic = np.array([r["cubic"] for r in rows])
    d = {k: np.array([r[k] for r in rows]) for k in ("explicit_fast", "explicit_smooth", "signed_overshoot")}
    gradients = {k: np.array([[r[k+"_g"+axis] for axis in "xyz"] for r in rows])
                 for k in ("ref", "linear", "cubic")}
    summary = {"n": len(rows)}
    for name, pred in (("linear", linear), ("cubic", cubic)):
        for ref in ("explicit_fast", "explicit_smooth"):
            summary[f"{name}_mae_{ref}"] = float(np.abs(pred - d[ref]).mean())
        gerr = np.linalg.norm(gradients[name] - gradients["ref"], axis=1)
        summary[f"{name}_mean_gradient_error"] = float(gerr.mean())
        summary[f"{name}_p95_gradient_error"] = float(np.percentile(gerr, 95))
    summary["cubic_overshoots"] = int((np.abs(d["signed_overshoot"]) > 1e-10).sum())
    summary["max_overshoot"] = float(np.max(np.abs(d["signed_overshoot"])))
    return summary


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("probe", type=Path)
    args.add_argument("inputs", type=Path, help="extracted torsion_proposals/fixtures/inputs.tar.xz")
    args.add_argument("results", type=Path)
    args.add_argument("--spacing", type=float, choices=(0.375, 0.25, 0.1875), default=0.375)
    a = args.parse_args()
    a.results.mkdir(exist_ok=True, parents=True)
    rows = []
    for i, pdb in enumerate(PDB_IDS):
        rows.extend(run(a.probe.resolve(), a.inputs / pdb, i, a.spacing))
    with (a.results/"field_probes.csv.gz").open("wb") as file:
        with gzip.GzipFile(filename="", fileobj=file, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                w = csv.DictWriter(stream, list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    groups = [("all", rows)]
    for key in ("target", "pdb", "radius_A"):
        for value in sorted({r[key] for r in rows}):
            groups.append((f"{key}={value}", [r for r in rows if r[key] == value]))
    groups.extend([("near_cutoff", [r for r in rows if r["cutoff_contact"]]),
                   ("away_from_cutoff", [r for r in rows if not r["cutoff_contact"]])])
    result = {name: summarise(group) for name, group in groups}
    (a.results/"field_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
