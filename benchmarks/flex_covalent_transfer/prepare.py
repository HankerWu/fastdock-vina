"""Prepare the predeclared human cognate complexes and unpack the frozen assay."""

import argparse
import csv
import json
import subprocess
import sys
import tarfile
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "geometry_metric"))
import panel

STRUCTURES = {
    "5DS3": ("09L", "PARP1"),
    "1F5K": ("BEN", "PLAU"),
    "1F5L": ("AMR", "PLAU"),
}
STRUCTURES.update(
    {
        "3P5O": ("EAM", "BRD4"),
        "5EHR": ("5OD", "PTPN11"),
        "6CMS": ("5OD", "PTPN11"),
        "3OCS": ("746", "BTK"),
    }
)
TRANSFER = HERE.parent / "archive_overlap" / "transfer"
MEASUREMENTS = TRANSFER / "fabp4" / "assay" / "measured.csv"


def records(lines):
    return Counter(
        l[12:26] + l[30:54] + l[70:]
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
    )


def prepare_pose(work, pdb, code, target):
    root = work / pdb
    meta = panel.prepare(work, pdb, code, target)
    observed = np.asarray(meta["mapped"][0])
    prot = (root / "receptor.pdbqt").read_text().splitlines()
    atoms = [
        l
        for l in prot
        if l.startswith("ATOM")
        and l[77:].strip() != "HD"
        and l[17:20] not in ("ALA", "GLY", "PRO", "CYS")
        and l[12:16].strip() not in ("O", "OXT", "N", "CA", "C")
    ]
    dist, atom = min(
        (float(np.linalg.norm(observed - panel.base.coords(l), axis=1).min()), l)
        for l in atoms
    )
    res_id = atom[22:27].strip()
    with (root / "flex_prepare.log").open("w") as out:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                str(root / "chain-A.pdb"),
                "-f",
                "A:" + res_id,
                "-p",
                str(root / "partition"),
            ],
            check=True,
            stdout=out,
            stderr=subprocess.STDOUT,
        )
    rigid, side = [
        records((root / f"partition_{name}.pdbqt").read_text().splitlines())
        for name in ("rigid", "flex")
    ]
    if rigid & side or rigid + side != records(prot):
        raise RuntimeError("incomplete receptor partition")
    return dict(
        pdb=pdb,
        kind="pose",
        target=target,
        ligand=code,
        status="ready",
        site=atom[17:20] + res_id,
        site_distance=dist,
        maps=len(meta["mapped"]),
    )


def unpack(archive, names, root):
    with tarfile.open(archive) as tar:
        for member in names:
            if member in tar.getnames():
                dst = root / names[member]
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_bytes(tar.extractfile(member).read())


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("work", type=Path)
    w = arg.parse_args().work.resolve()
    w.mkdir(parents=True, exist_ok=True)
    manifest = []
    for pdb, (code, target) in STRUCTURES.items():
        try:
            manifest.append(prepare_pose(w, pdb, code, target))
        except (
            RuntimeError,
            ValueError,
            KeyError,
            subprocess.CalledProcessError,
        ) as exc:
            manifest.append(
                dict(
                    pdb=pdb,
                    kind="pose",
                    target=target,
                    ligand=code,
                    status="excluded",
                    reason=str(exc).replace(str(w), "WORK"),
                )
            )
        print(manifest[-1], flush=True)
    labels = list(csv.DictReader(MEASUREMENTS.open()))
    names = {r["pdb"] + "/ligand.pdbqt": "assay/" + r["id"] + ".pdbqt" for r in labels}
    for p in ("2NNQ", "3FR2", "1RBP", "3FMZ"):
        names[p + "/ligand.pdbqt"] = p + "/ligand.pdbqt"
        names[p + "/preparation.json"] = p + "/metadata.json"
    for a in (
        HERE.parent / "cubic_grid/human_panel/prepared.tar.xz",
        TRANSFER / "inputs.tar.xz",
        TRANSFER / "fabp4/assay/new_inputs.tar.xz",
    ):
        unpack(a, names, w)
    parts = {
        p + "_" + k + ".pdbqt": p + "/partition_" + k + ".pdbqt"
        for p in ("2NNQ", "3FR2", "1RBP", "3FMZ")
        for k in ("flex", "rigid")
    }
    unpack(HERE.parent / "cubic_grid/flex_replay/flex-pools.tar.xz", parts, w)
    for rel in list(names.values()) + list(parts.values()):
        if not (w / rel).is_file():
            raise ValueError(f"Missing archive member {rel}")
    (w / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
