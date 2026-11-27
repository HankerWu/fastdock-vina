"""Overlay seven observed Aurora A structures on frozen 2X6E poses without a loop graft."""

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
DONORS = ("1OL5", "1OL7", "1MQ4", "2J4Z", "2J50", "2W1G", "6CPG")
NEW = {"1OL5", "1OL7", "1MQ4"}
SEGMENTS = {"281-291": tuple(range(281, 292)), "303-305": (303, 304, 305)}
BACKBONE = {"N", "CA", "C", "O"}
SEEDS = ("witness", "60913", "60943", "61211", "61231")


def xyz(l):
    return [float(l[k : k + 8]) for k in (30, 38, 46)]


def record(work, pdb):
    if pdb in NEW:
        raw = gzip.decompress((HERE / "donors" / f"{pdb}.pdb.gz").read_bytes())
    else:
        path = (work / "raw" if pdb in {"2J4Z", "2X6E"} else work) / f"{pdb}.pdb"
        raw = path.read_bytes()
    lines = raw.decode().splitlines()
    if not any(
        l.startswith("DBREF ") and l[12] == "A" and "O14965" in l for l in lines
    ):
        raise ValueError(f"Unknown Aurora A chain A mapping for {pdb}")
    return lines, hashlib.sha256(raw).hexdigest()


def atoms(lines):
    values = {}
    for l in lines:
        if (
            (l.startswith("ATOM  ") or (l.startswith("HETATM") and l[17:20] == "TPO"))
            and l[21] == "A"
            and l[16] in " A"
            and l[76:78].strip() not in {"H", "D"}
        ):
            assert l[26] == " "
            key = int(l[22:26]), l[17:20], l[12:16].strip()
            if key in values:
                raise ValueError(("duplicate heavy atom", key))
            values[key] = np.array(xyz(l))
    return values


def ca(at):
    return {
        (i, "THR" if res == "TPO" else res): pos
        for (i, res, name), pos in at.items()
        if name == "CA"
    }


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(
                {k: f"{v:.6f}" if isinstance(v, float) else v for k, v in row.items()}
            )


def ligands(work):
    result = []
    for drug in ("at9283", "danusertib"):
        for idx, seed in enumerate(SEEDS, 1):
            lines = (
                (work / "coordinates" / f"2X6E_{drug}_2X6E_pocket_{idx}.pdbqt")
                .read_text()
                .splitlines()
            )
            rows = [
                l
                for l in lines
                if l.startswith(("ATOM", "HETATM"))
                and l.split()[-1] not in {"H", "HD", "HS"}
            ]
            result.append(
                (
                    drug,
                    seed,
                    [int(l[6:11]) for l in rows],
                    np.array([xyz(l) for l in rows]),
                )
            )
    return result


def rmsd(x, y):
    return np.sqrt(np.square(x - y).sum(axis=-1).mean()).item()


def alignment(src, dest):
    rotation, _ = Rotation.align_vectors(dest - dest.mean(0), src - src.mean(0))
    origin, destination = src.mean(0), dest.mean(0)

    def apply(points):
        return rotation.apply(np.array(points) - origin) + destination

    return apply


def replay(work, out, model_paths):
    receiver, _ = record(work, "2X6E")
    observed = atoms(receiver)
    dest_ca = ca(observed)
    ym4 = [
        xyz(l)
        for l in receiver
        if l.startswith("HETATM")
        and l[21] == "A"
        and l[17:20] == "YM4"
        and l[16] in " A"
        and l[76:78].strip() not in {"H", "D"}
    ]
    missing = {
        int(i): r
        for r, i in re.findall(
            r"^REMARK 465\s+([A-Z]{3}) A\s+(\d+)\s*$", "\n".join(receiver), re.MULTILINE
        )
    }
    pocket = {k for k, v in dest_ca.items() if cdist([v], ym4).min() <= 12}
    poses = ligands(work)
    donors, provenance = {}, {}
    for pdb in DONORS:
        lines, digest = record(work, pdb)
        donors[pdb] = atoms(lines)
        provenance[pdb] = {
            "url": f"https://files.rcsb.org/download/{pdb}.pdb",
            "sha256": digest,
            "chain": "A",
            "resolution": next(
                l[10:].strip() for l in lines if l.startswith("REMARK   2 RESOLUTION.")
            ),
            "title": " ".join(l[10:].strip() for l in lines if l.startswith("TITLE")),
        }
    shared = set(dest_ca).intersection(*(ca(a).keys() for a in donors.values()))
    built = {
        path.as_posix(): atoms((work / path).read_text().splitlines())
        for path in model_paths
    }
    frames, distances, models = [], [], []
    for name, keys in (("pocket", sorted(pocket & shared)), ("whole", sorted(shared))):
        dest = np.array([dest_ca[k] for k in keys])
        for pdb, structure in donors.items():
            original_ca = ca(structure)
            src = np.array([original_ca[k] for k in keys])
            move = alignment(src, dest)

            anchors = [
                (n, res)
                for n, res in dest_ca
                if n in (280, 292, 302, 306) and (n, res) in original_ca
            ]
            anchor_dist = {
                i: float(np.linalg.norm(move(original_ca[i, res]) - dest_ca[i, res]))
                for i, res in anchors
            }
            common = {"donor": pdb, "frame": name}
            frames.append(
                {
                    **common,
                    "n_ca": len(keys),
                    "fit_rmsd_A": rmsd(move(src), dest),
                    **{
                        f"anchor_{i}_A": anchor_dist.get(i, "")
                        for i in (280, 292, 302, 306)
                    },
                }
            )
            for segment, span in SEGMENTS.items():
                labels = [k for k in structure if k[0] in span]
                numbers = {k[0] for k in labels}
                for i, res, _ in labels:
                    assert ("THR" if res == "TPO" else res) == missing[i]
                header = {
                    **common,
                    "segment": segment,
                    "deposited_residues": " ".join(map(str, sorted(numbers))),
                    "n_deposited_heavy": len(labels),
                }
                xyz_loop = (
                    move([structure[k] for k in labels]) if labels else np.empty((0, 3))
                )
                for drug, seed, serials, pts in poses:
                    near = {
                        "min_A": "",
                        "nearest_donor_atom": "",
                        "ligand_serial": "",
                        "min_backbone_A": "",
                    }
                    if len(xyz_loop):
                        ds = cdist(xyz_loop, pts)
                        i, j = np.unravel_index(np.argmin(ds), ds.shape)
                        nr, res, atom = labels[i]
                        near = {
                            "min_A": float(ds[i, j]),
                            "nearest_donor_atom": f"{res}{nr}.{atom}",
                            "ligand_serial": serials[j],
                            "min_backbone_A": float(
                                ds[[a in BACKBONE for _, _, a in labels]].min()
                            ),
                        }
                    distances.append({**header, "drug": drug, "pose": seed, **near})
                for path, model in built.items():
                    active = [
                        (i, (nr, "THR" if res == "TPO" else res, atom))
                        for i, (nr, res, atom) in enumerate(labels)
                        if atom in BACKBONE
                    ]
                    value = (
                        rmsd(
                            np.array([model[key] for _, key in active]),
                            xyz_loop[[i for i, _ in active]],
                        )
                        if active
                        else ""
                    )
                    models.append(
                        {
                            **header,
                            "model": path,
                            "shared_backbone_atoms": len(active),
                            "backbone_rmsd_A": value,
                        }
                    )
    out.mkdir(parents=True, exist_ok=True)
    write(out / "donor_frames.csv", frames)
    write(out / "donor_distances.csv", distances)
    if models:
        write(out / "donor_model_rmsd.csv", models)
    (out / "donor_sources.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument(
        "--models",
        type=Path,
        nargs="+",
        default=[],
        help="Optional completed receptor paths relative to WORK",
    )
    args = ap.parse_args()
    replay(args.work, args.out, args.models)
