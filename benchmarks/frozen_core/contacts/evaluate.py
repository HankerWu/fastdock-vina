"""Evaluate contact transfer on archived coordinates; write one score per pose and one outcome per run."""

import argparse
import csv
import hashlib
import json
import sys
import tarfile
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
GEOM = BENCH / "geometry_comparison"
sys.path.insert(0, str(GEOM / "upstream"))
from verify import references


@lru_cache(maxsize=12)
def archive(path):
    # Small compressed evidence bundles. Keep only text inputs used by this experiment.
    with tarfile.open(path) as stream:
        return {
            m.name: stream.extractfile(m).read().decode()
            for m in stream
            if m.isfile() and m.name.endswith((".pdbqt", ".json"))
        }


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def atoms(pdbqt):
    elements, points, serial = [], [], []
    aliases = {"A": "C", "OA": "O", "NA": "N", "SA": "S"}
    for line in pdbqt.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = line.split()[-1]
            if atom in {"H", "HD", "HS", "W", "G0", "G1", "G2", "G3"}:
                continue
            elements.append(aliases.get(atom, atom))
            points.append(tuple(float(line[i : i + 8]) for i in (30, 38, 46)))
            serial.append(int(line[6:11]))
    return np.array(points), np.array(elements), serial


def profile(rec, ligand):
    r, re, _ = rec
    x, xe, _ = ligand
    d = cdist(r, x)
    general = np.clip((6.0 - d.min(axis=1)) / 2.0, 0, 1)
    rp, xp = np.isin(re, ["N", "O"]), np.isin(xe, ["N", "O"])
    polar = (
        np.clip((4.5 - d[rp][:, xp].min(axis=1)) / 1.5, 0, 1)
        if xp.any()
        else np.zeros(rp.sum())
    )
    return general, polar


def tanimoto(a, b):
    dot = a @ b
    norm = a @ a + b @ b - dot
    return dot / norm if norm else None


def shape(ligand, reference):
    x, xe, _ = ligand
    y, ye, _ = reference
    d = cdist(x, y, metric="sqeuclidean")
    d[xe[:, None] != ye] = np.inf
    return (
        np.exp(-d.min(axis=1) / 4.5).mean() + np.exp(-d.min(axis=0) / 4.5).mean()
    ) / 2


def native_components(block):
    def value(label):
        return float(
            next(l for l in block.splitlines() if l.startswith("REMARK " + label))
            .split(":")[1]
            .split()[0]
        )

    score = value("VINA RESULT:")
    raw = value("INTER + INTRA:") - value("UNBOUND:")
    # For default Vina the model-derived rotor count is integral or half-integral.
    ntors = round((raw / score - 1) / 0.05846 * 2) / 2
    divisor = 1 + 0.05846 * ntors
    if abs(score * divisor - raw) > 0.003:
        raise ValueError("Vina torsion scaling reconstruction failed")
    return score, divisor


def input_data(name, code):
    if name == "paired":
        f = archive(BENCH / "paired_crossdocking" / "fixtures.tar.xz")
        return (
            json.loads(f["metadata.json"])[code],
            f[code + "/ligand.pdbqt"],
            f[code + "/receptor.pdbqt"],
        )
    if name == "old":
        f = archive(BENCH / "cubic_grid" / "human_panel" / "prepared.tar.xz")
        meta = "preparation.json"
    elif name == "brdt":
        f = archive(GEOM / "water_transfer" / "inputs.tar.xz")
        code = "inputs/" + code
        meta = "metadata.json"
    elif code in ("7RWP", "7RWQ"):
        f = archive(GEOM / "bptf_water" / "raw.tar.xz")
        code = "external/" + code
        meta = "metadata.json"
    else:
        f = archive(GEOM / "pose" / "inputs.tar.xz")
        meta = "metadata.json"
    return (
        json.loads(f[code + "/" + meta]),
        f[code + "/ligand.pdbqt"],
        f[code + "/receptor.pdbqt"],
    )


def distinct_template(name, q, r):
    qm, _, _ = input_data(name, q)
    rm, ligand, receptor = input_data(name, r)
    if q == r or qm["canonical_smiles"] == rm["canonical_smiles"]:
        raise ValueError("A different reference compound is required")
    _, elements, serial = atoms(ligand)
    ids = rm["serial"]
    donor = (
        np.asarray(rm["mapped"][0]),
        np.array([elements[serial.index(i)] for i in ids]),
        ids,
    )
    return atoms(receptor), donor


def cases():
    root = GEOM / "upstream" / "reciprocal"
    for stage in ("fixed", "matched", "per_complex"):
        d = root if stage == "fixed" else root / stage
        for row in rows(d / "cells.csv"):
            if row["arm"] != "native":
                continue
            q, r, seed = [row[k] for k in ("query", "receiver", "seed")]
            raw = archive(d / "raw.tar.xz")[f"{q}-{r}-{seed}-native.pdbqt"]
            assert hashlib.sha256(raw.encode()).hexdigest() == row["pool_sha256"]
            qm, rm = input_data("old", q)[0], input_data("old", r)[0]
            ref = references(qm, rm)
            rms = []
            for block in raw.split("ENDMDL"):
                if "MODEL" in block:
                    xyz, _, ids = atoms(block)
                    xyz = xyz[[ids.index(i) for i in qm["serial"]]]
                    rms.append(
                        float(
                            np.sqrt(
                                np.mean(np.sum((ref - xyz) ** 2, axis=-1), axis=-1)
                            ).min()
                        )
                    )
            assert abs(rms[0] - float(row["first"])) < 1e-8
            yield (stage, row["target"], q, r, seed), raw, rms, "old"
    panels = [
        ("bptf", "bptf_water", "confirm_modes.csv", "raw.tar.xz", "dry/confirm", "dry"),
        ("brdt", "water_transfer", "modes.csv", "outputs.tar.xz", "dry/poses", "dry"),
        (
            "bcl6_120k",
            "pose",
            "high_modes.csv",
            "high_outputs.tar.xz",
            "poses",
            "armijo",
        ),
        ("bcl6_12k", "pose", "low_modes.csv", "low_outputs.tar.xz", "poses", "armijo"),
    ]
    for stage, folder, table, bundle, prefix, arm in panels:
        groups = defaultdict(list)
        for row in rows(GEOM / folder / table):
            if (
                row["arm"] != arm
                or row["query"] == row["receiver"]
                or (stage.startswith("bcl6") and row["target"] != "BCL6")
            ):
                continue
            groups[
                tuple(row[k] for k in ("target", "query", "receiver", "seed"))
            ].append(row)
        for (target, q, r, seed), modes in groups.items():
            raw = archive(GEOM / folder / bundle)[
                f"{prefix}/{q}-{r}-{seed}-{arm}.pdbqt"
            ]
            yield (
                (stage, target, q, r, seed),
                raw,
                [float(m["rmsd_A"]) for m in modes],
                stage,
            )
    d = BENCH / "paired_crossdocking"
    for cap in (3000, 12000):
        groups = defaultdict(list)
        for row in rows(d / "results" / f"poses_{cap}.csv"):
            if row["arm"] == "native" and row["ligand"] != row["receiver"]:
                groups[
                    tuple(row[k] for k in ("target", "ligand", "receiver", "seed"))
                ].append(row)
        for (target, q, r, seed), modes in groups.items():
            raw = archive(d / f"poses-{cap}.tar.xz")[
                f"runs_{cap}/{q}_{r}_{seed}_native.pdbqt"
            ]
            yield (
                (f"paired_{cap}", target, q, r, seed),
                raw,
                [float(m["rmsd"]) for m in modes],
                "paired",
            )


def evaluate():
    result, cells = [], []
    for key, raw, rms, name in cases():
        stage, target, q, r, seed = key
        rec, donor = distinct_template(name, q, r)
        y = profile(rec, donor)
        pool_rows = []
        for idx, block in enumerate(b for b in raw.split("ENDMDL") if "MODEL" in b):
            ligand = atoms(block)
            score, divisor = native_components(block)
            x = profile(rec, ligand)
            g, p = [tanimoto(a, b) for a, b in zip(x, y)]
            overlaps = [z for z in (g, p) if z is not None]
            overlap = np.mean(overlaps) if overlaps else 0.0
            rank = {
                "native": score,
                "element": score - 6 * shape(ligand, donor) / divisor,
                "contact": score - 6 * overlap / divisor,
                "general": score - 6 * (g or 0.0) / divisor,
                "contact_only": -overlap,
            }
            line = dict(
                stage=stage,
                target=target,
                query=q,
                receiver=r,
                seed=seed,
                model=idx + 1,
                rmsd=rms[idx],
                divisor=divisor,
                **rank,
            )
            pool_rows.append(line)
        assert len(rms) == len(pool_rows)
        cell = {
            "stage": stage,
            "target": target,
            "query": q,
            "receiver": r,
            "seed": seed,
            "pool_sha256": hashlib.sha256(raw.encode()).hexdigest(),
            "modes": len(pool_rows),
            "oracle": min(rms),
        }
        for arm in rank:
            best = min(pool_rows, key=lambda x: (x[arm], x["model"]))
            cell[arm] = best["rmsd"]
            cell[arm + "_rank"] = best["model"]
        result.extend(pool_rows)
        cells.append(cell)
    return result, cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    out = parser.parse_args().output
    out.mkdir(parents=True, exist_ok=True)
    modes, cells = evaluate()
    for tag, data in (("modes", modes), ("cells", cells)):
        with (out / (tag + ".csv")).open("w", newline="") as fh:
            csv_out = csv.DictWriter(fh, data[0].keys())
            csv_out.writeheader()
            csv_out.writerows(data)
    for stage in dict.fromkeys(c["stage"] for c in cells):
        batch = [c for c in cells if c["stage"] == stage]
        print(
            stage,
            len(batch),
            *(
                sum(c[a] <= 2 for c in batch)
                for a in (
                    "native",
                    "element",
                    "contact",
                    "general",
                    "contact_only",
                    "oracle",
                )
            ),
        )


if __name__ == "__main__":
    main()
