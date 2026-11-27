"""Independent reconstruction from deposited structures, PDBQT and printed energies."""

import argparse
import hashlib
import json
import shlex
import re
import subprocess
import tempfile
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

K = ["panel", "query", "receptor", "seed", "cap", "search"]
A = ["native_native", "field_native", "native_field", "field_field"]


def point(line):
    return tuple(float(line[i : i + 8]) for i in (30, 38, 46))


@lru_cache(maxsize=None)
def typed(receptor, ligand, sidecar):
    carbon_names = {
        "PHE": "CG CD1 CE1 CZ CE2 CD2",
        "TYR": "CG CD1 CE1 CZ CE2 CD2",
        "TRP": "CG CD1 CD2 CE2 CE3 CZ2 CZ3 CH2",
        **{h: "CG CD2 CE1" for h in ("HIS", "HIP", "HID", "HIE", "HSD", "HSE", "HSP")},
    }
    atom_centers = [
        point(line)
        for line in receptor.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
        and line[12:16].strip() in carbon_names.get(line[17:20], "").split()
    ]
    spec = sidecar.read_text().splitlines()
    assert np.allclose([float(v) for v in spec[0].split()], [0.078, 4.49, 0.75, 8.0])
    assert Counter(tuple(map(float, l.split())) for l in spec[2:]) == Counter(
        atom_centers
    )
    graph = ligand.read_text().splitlines()
    (chem_smi,) = [
        l.split()[2]
        for l in graph
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX ")
    ]
    mol = Chem.MolFromSmiles(chem_smi)
    assert mol is not None and not mol.HasSubstructMatch(Chem.MolFromSmarts("[c]"))
    hits = mol.GetSubstructMatches(Chem.MolFromSmarts("[C;X4;!H0]"))
    mapping = {}
    for line in graph:
        if line.startswith("REMARK SMILES IDX "):
            ints = list(map(int, line.split()[3:]))
            mapping.update(zip(ints[::2], ints[1::2], strict=True))
    serials = [mapping[i + 1] for (i,) in hits]
    assert set(serials) == set(map(int, spec[1].split()))
    return serials, cKDTree(atom_centers)


@lru_cache(maxsize=None)
def raw_atoms(work, pdb):
    return (work / "raw" / f"{pdb}.pdb").read_text().splitlines()


def ca_positions(raw, record):
    positions = {}
    for line in raw:
        if (
            line.startswith("ATOM  ")
            and line[12:16].strip() == "CA"
            and line[21] == record["chain"]
        ):
            key = line[22:27].strip() + ":" + line[17:20].strip()
            if key in record["ca"] and np.isclose(point(line), record["ca"][key]).all():
                positions[key] = np.array(point(line))
    return positions


@lru_cache(maxsize=None)
def geometry(work, panel, query, rec):
    meta = json.loads((work / "metadata.json").read_text())
    sugar = panel in ("GCK", "GLUT3")
    a = meta[f"{rec}-{query}"] if sugar else meta[query]
    raw = raw_atoms(work, rec if sugar else query)
    refs = np.array(a["mapped"])
    deposited = {
        point(l)
        for l in raw
        if l.startswith("HETATM")
        and l[17:20] == a["comp"]
        and (sugar or l[21] == a["chain"])
    }
    if not all(tuple(p) in deposited for ref in refs for p in ref):
        raise ValueError("Reference atom missing from original experimental ligand")
    if not sugar and query != rec:
        b = meta[rec]
        q_ca, r_ca = ca_positions(raw, a), ca_positions(raw_atoms(work, rec), b)
        keys = sorted(
            k
            for k in q_ca.keys() & r_ca.keys()
            if np.linalg.norm(q_ca[k] - a["center"]) < 14
            and np.linalg.norm(r_ca[k] - b["center"]) < 14
        )
        assert len(keys) >= 12
        q, r = (np.array([positions[k] for k in keys]) for positions in (q_ca, r_ca))
        qc, rc = q.mean(axis=0), r.mean(axis=0)
        rotate = Rotation.align_vectors(r - rc, q - qc)[0]
        refs = np.array([rotate.apply(ref - qc) + rc for ref in refs])
    return refs, a["serial" if sugar else "protonated_serial"]


def rescore_close(binary, pdbqt, work, options, contenders):
    """On archived coordinates use a common zero unbound reference for every candidate."""
    blocks = pdbqt.read_text().split("ENDMDL")
    out = []
    with tempfile.TemporaryDirectory() as temp:
        ligand = Path(temp) / "pose.pdbqt"
        for p in contenders:
            block = blocks[p["rank"] - 1]
            ligand.write_text(
                "\n".join(
                    line for line in block.splitlines() if not line.startswith("MODEL")
                )
                + "\n"
            )
            cmd = [
                str(binary.resolve()),
                "--score_only",
                "--unbound_energy",
                "0",
                "--cpu",
                "1",
                "--ligand",
                str(ligand),
                "--receptor",
                str((work / options[options.index("--receptor") + 1]).resolve()),
            ]
            for name in (
                "center_x",
                "center_y",
                "center_z",
                "size_x",
                "size_y",
                "size_z",
            ):
                flag = "--" + name
                cmd.extend([flag, options[options.index(flag) + 1]])
            text = subprocess.check_output(cmd, text=True)
            dry = float(
                re.search(r"Estimated Free Energy of Binding +: +([-0-9.]+)", text)[1]
            )
            out.append(
                dict(
                    rank=p["rank"],
                    common_reference_dry=dry,
                    corrected=dry + 0.3 * p["factor"] * p["raw_chpi"],
                )
            )
    return out


def reconstruct(result, workspaces, binary):
    modes, cells, sensitivity, changes = [], [], [], []
    costs = pd.read_csv(result / "costs.csv")
    for row in costs.itertuples():
        keys = {k: getattr(row, k) for k in K}
        work = workspaces[row.panel]
        options = shlex.split(row.args)
        accepts, rec_field = typed(
            work / options[options.index("--receptor") + 1],
            work / options[options.index("--ligand") + 1],
            result / "typing" / f"{row.panel}_{row.query}_{row.receptor}.txt",
        )
        tag = "_".join(str(keys[k]) for k in K)
        path = result / "pools" / (tag + ".pdbqt")
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row.sha256
        refs, order = geometry(work, row.panel, row.query, row.receptor)
        pool = []
        torsions = set(range(65))
        for block in path.read_text().split("ENDMDL"):
            if not block.strip():
                continue
            lines = block.strip().splitlines()
            rank = int(lines[0].split()[1])
            xyz = {
                int(l[6:11]): point(l)
                for l in lines
                if l.startswith(("ATOM  ", "HETATM"))
            }
            remarks = {
                l.split(":")[0]: l.split(":")[1].strip().split()[0]
                for l in lines
                if l.startswith("REMARK ") and ":" in l
            }
            native = float(remarks["REMARK VINA RESULT"])
            unscaled = float(remarks["REMARK INTER + INTRA"]) - float(
                remarks["REMARK UNBOUND"]
            )
            torsions &= {
                n
                for n in torsions
                if abs(native * (1 + 0.05846 * n) - unscaled)
                <= 0.0005 * (1 + 0.05846 * n) + 0.0010001
            }
            sites = np.asarray([xyz[i] for i in accepts])
            raw = 0.0
            for x, neighbors in zip(
                sites, rec_field.query_ball_point(sites, 8), strict=True
            ):
                d = np.linalg.norm(rec_field.data[neighbors] - x, axis=1)
                raw -= (
                    float(np.exp(-0.5 * ((d[d < 8] - 4.49) / 0.75) ** 2).sum()) * 0.26
                )
            rms = float(
                np.sqrt(
                    np.square(np.array([xyz[i] for i in order]) - refs).sum(-1).mean(-1)
                ).min()
            )
            pool.append(dict(**keys, rank=rank, rmsd=rms, native=native, raw_chpi=raw))
        if len(torsions) != 1:
            raise ValueError(
                f"Rounded energies do not determine a unique rotor count for {tag}: {torsions}"
            )
        f = 1 / (1 + 0.05846 * torsions.pop())
        for p in pool:
            p.update(factor=f, corrected=p["native"] + 0.3 * f * p["raw_chpi"])
        modes.extend(pool)
        for ranker, col in (("native", "rank"), ("field", "corrected")):
            seq = sorted(pool, key=lambda p: (p[col], p["rank"]))
            winner = seq[0]
            cells.append(
                dict(
                    **keys,
                    rerank=ranker,
                    output_count=len(seq),
                    winner=winner["rank"],
                    top1=winner["rmsd"],
                    best5=min(x["rmsd"] for x in seq[:5]),
                    pool_min=min(x["rmsd"] for x in seq),
                    regret=winner["rmsd"] - min(x["rmsd"] for x in seq),
                )
            )
            if ranker == "field" and winner["rank"] != 1:
                changes.append(
                    dict(
                        **keys,
                        winner=winner["rank"],
                        before=pool[0]["rmsd"],
                        after=winner["rmsd"],
                        printed_dry_gap=winner["native"] - pool[0]["native"],
                        pair_advantage=0.3
                        * f
                        * (pool[0]["raw_chpi"] - winner["raw_chpi"]),
                    )
                )
            if ranker == "field":
                contenders = [
                    p
                    for p in seq
                    if p["corrected"] - winner["corrected"] <= 0.0010000001
                ]
                if len(contenders) > 1:
                    sensitivity.append(
                        dict(
                            **keys,
                            best_rank=winner["rank"],
                            best_rmsd=winner["rmsd"],
                            score_gap=seq[1]["corrected"] - winner["corrected"],
                            runner_up_rank=seq[1]["rank"],
                            runner_up_rmsd=seq[1]["rmsd"],
                            coordinate_step_angstrom=0.001,
                            printed_score_step=0.001,
                            score_only_common_reference=(
                                rescore_close(binary, path, work, options, contenders)
                                if binary
                                else None
                            ),
                            min_plausible=min(p["rmsd"] for p in contenders),
                            max_plausible=max(p["rmsd"] for p in contenders),
                        )
                    )
    cells, modes = pd.DataFrame(cells), pd.DataFrame(modes)
    max_errors = {}
    for frame, name, key in (
        (modes, "modes", K + ["rank"]),
        (cells, "cells", K + ["rerank"]),
    ):
        reference = pd.read_csv(result / (name + ".csv")).set_index(key).sort_index()
        frame = frame.set_index(key).sort_index()
        assert len(reference) == len(frame) and reference.index.equals(frame.index)
        err = (frame[reference.columns] - reference).abs().max().max()
        assert err < 1e-8, (name, err)
        max_errors[name] = float(err)
    cells["arm"] = cells.search + "_" + cells.rerank
    output = []
    for (panel, cap), group in cells.groupby(["panel", "cap"]):
        pivot = group.pivot(
            index=["query", "receptor", "seed"], columns="arm", values="top1"
        )
        means = pivot[A].mean()
        bounds = {a: [float(means[a]), float(means[a])] for a in A}
        for s in sensitivity:
            if s["panel"] == panel and s["cap"] == cap:
                arm = s["search"] + "_field"
                bounds[arm][0] += (s["min_plausible"] - s["best_rmsd"]) / len(pivot)
                bounds[arm][1] += (s["max_plausible"] - s["best_rmsd"]) / len(pivot)
        interaction_min = (
            bounds[A[3]][0] - bounds[A[2]][1] - bounds[A[1]][1] + bounds[A[0]][0]
        )
        interaction_max = (
            bounds[A[3]][1] - bounds[A[2]][0] - bounds[A[1]][0] + bounds[A[0]][1]
        )
        output.append(
            dict(
                panel=panel,
                cap=int(cap),
                conditions=len(pivot),
                field_field_within_printed_score_tolerance=bounds[A[3]],
                interaction_within_printed_score_tolerance=[
                    interaction_min,
                    interaction_max,
                ],
                **{a: float(means[a]) for a in A},
                interaction=float(
                    means[A[3]] - means[A[2]] - means[A[1]] + means[A[0]]
                ),
            )
        )
    return dict(
        pools=len(costs),
        poses=len(modes),
        cells=len(cells),
        max_errors=max_errors,
        means=output,
        switched=len(changes),
        largest_reorder=sorted(changes, key=lambda p: -abs(p["before"] - p["after"]))[
            :3
        ],
        close_corrected_winners=sensitivity,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--gba", type=Path)
    parser.add_argument("--glucose", type=Path)
    parser.add_argument("--gaa", type=Path)
    parser.add_argument(
        "--binary", type=Path, help="Optional Vina CLI for rounded coordinate replay"
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    results = reconstruct(
        args.results,
        dict(GBA1=args.gba, GAA=args.gaa, GCK=args.glucose, GLUT3=args.glucose),
        args.binary,
    )
    args.out.write_text(
        json.dumps(results, indent=2, default=lambda x: x.item()) + "\n"
    )
