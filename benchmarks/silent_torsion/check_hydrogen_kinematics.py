"""Pair two CLIs on an off-axis hydroxyl H, including a hydrogen-sensitive AD4 score."""

import argparse
import csv
import hashlib
import re
import subprocess
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
ATOM = ("ATOM", "HETATM")
AFF = re.compile(r"Estimated Free Energy of Binding\s+:\s+([-\d.]+)")
REFINE = re.compile(
    r"Local refinement: optimizer=\w+ gradient_evals=(\d+) objective_kcal_mol=(\S+)"
)


def write_ligand(path, origin, sign):
    def atom(ident, element, atomtype, offset):
        xyz = tuple(a + b for a, b in zip(origin, offset))
        return (
            f"ATOM  {ident:5} {element:^4} UNL     1    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00     0.000 {atomtype}\n"
        )

    path.write_text(
        "ROOT\n"
        + atom(1, "C", "C", (-1.5, 0, 0))
        + atom(2, "C", "C", (0, 0, 0))
        + "ENDROOT\nBRANCH   2   3\n"
        + atom(3, "O", "OA", (1.43, 0, 0))
        + atom(4, "H", "HD", (1.87, 0.82 * sign, 0.20 * sign))
        + "ENDBRANCH   2   3\nTORSDOF 1\n"
    )


def coords(text):
    return [
        tuple(float(line[i : i + 8]) for i in (30, 38, 46))
        for line in text.splitlines()
        if line.startswith(ATOM)
    ]


def compare(parent, candidate, w, ligand, receptor, sf, mode, flags):
    states = []
    for binary in (parent, candidate):
        pose = w / "out.pdbqt"
        pose.unlink(missing_ok=True)
        command = [
            binary,
            "--scoring",
            sf,
            "--ligand",
            str(ligand),
            "--cpu",
            "1",
            "--seed",
            "37",
            "--verbosity",
            "0",
            "--out",
            str(pose),
            *receptor,
            *flags,
        ]
        process = subprocess.run(command, capture_output=True, text=True, check=True)
        affinity = float(AFF.search(process.stdout).group(1))
        local = REFINE.search(process.stdout)
        output = pose.read_bytes() if pose.exists() else b""
        states.append((affinity, local.groups() if local else (), output))
    if states[0] != states[1]:
        raise ValueError(
            f"Score, accepted objective/count, or complete PDBQT changed: {sf}, {mode}"
        )
    score, detail, output = states[0]
    row = {
        "scoring": sf,
        "mode": mode,
        "score": score,
        "gradients": detail[0] if detail else None,
        "objective": detail[1] if detail else None,
        "pdbqt_sha256": hashlib.sha256(output).hexdigest() if output else None,
    }
    return row, output


def run(parent, candidate, out):
    basic = ROOT / "examples/basic"
    center = coords((basic / "1iep_ligand.pdbqt").read_text())[0]
    rows = []
    frozen = {}
    with TemporaryDirectory() as directory:
        w = Path(directory)
        header = (
            "GRID_PARAMETER_FILE e.gpf\nGRID_DATA_FILE e.fld\nMACROMOLECULE e.pdbqt\n"
            "SPACING 1.0\nNELEMENTS 12 12 12\nCENTER 0.0 0.0 0.0\n"
        )
        for suffix in ("C", "OA", "HD", "d", "e"):
            field = [
                str(y if suffix == "HD" else 0)
                for z in range(-6, 7)
                for y in range(-6, 7)
                for x in range(-6, 7)
            ]
            (w / f"e.{suffix}.map").write_text(header + "\n".join(field) + "\n")
        for sf in ("vina", "vinardo", "ad4"):
            receptor = (
                ["--maps", str(w / "e")]
                if sf == "ad4"
                else ["--receptor", str(basic / "1iep_receptor.pdbqt"), "--autobox"]
            )
            for sign in (1, -1):
                ligand = w / "input.pdbqt"
                write_ligand(ligand, (0, 0, 0) if sf == "ad4" else center, sign)
                options = [("fixed", ["--score_only"]), ("local", ["--local_only"])]
                if sf == "vina":
                    for alg in ("armijo", "geometry", "rigid_metric"):
                        for cap in (1, 7, 40):
                            options.append(
                                (
                                    f"{alg}-{cap}",
                                    [
                                        "--local_only",
                                        "--refinement_optimizer",
                                        alg,
                                        "--max_evals",
                                        str(cap),
                                    ],
                                )
                            )
                    options.append(("map-local", ["--local_only", "--no_refine"]))
                for mode, flags in options:
                    row, output = compare(
                        parent, candidate, w, ligand, receptor, sf, mode, flags
                    )
                    if mode == "fixed":
                        frozen[sf, sign] = row["score"]
                    if mode.endswith("-1") and coords(output.decode()) != coords(
                        ligand.read_text()
                    ):
                        raise ValueError("One-evaluation run moved some atoms")
                    if (
                        mode == "armijo-40"
                        and coords(output.decode())[-1]
                        == coords(ligand.read_text())[-1]
                    ):
                        raise ValueError(
                            "Fixture did not move its H atom on strict local optimization"
                        )
                    rows.append({"fixture": "ethanol", "sign": sign, **row})
        with tarfile.open(ROOT / "benchmarks/galectin_affinity/evidence.tar.xz") as arc:
            receptor = w / "gal3.pdbqt"
            receptor.write_bytes(arc.extractfile("prepared/1KJL.pdbqt").read())
            source = arc.extractfile("pools/1KJL_5_141421_native.pdbqt").read().decode()
        ligand = w / "gal3_5.pdbqt"
        first = source.split("ENDMDL", 1)[0].splitlines()[1:]
        ligand.write_text("\n".join(first) + "\n")
        opts = [
            ("fixed", ["--score_only"]),
            ("native", ["--local_only"]),
            ("map-local", ["--local_only", "--no_refine"]),
        ]
        for algorithm in ("armijo", "geometry", "rigid_metric"):
            for cap in (1, 7, 40):
                opts.append(
                    (
                        f"{algorithm}-{cap}",
                        [
                            "--local_only",
                            "--refinement_optimizer",
                            algorithm,
                            "--max_evals",
                            str(cap),
                        ],
                    )
                )
        for mode, flags in opts:
            row, _ = compare(
                parent,
                candidate,
                w,
                ligand,
                ["--receptor", str(receptor), "--autobox"],
                "vina",
                mode,
                flags,
            )
            rows.append({"fixture": "1KJL_5_141421", "sign": None, **row})
    if (
        frozen["vina", 1] != frozen["vina", -1]
        or frozen["vinardo", 1] != frozen["vinardo", -1]
    ):
        raise ValueError("XS score depends on ignored H orientation")
    if abs(frozen["ad4", 1] - frozen["ad4", -1]) < 1:
        raise ValueError("HD-specific map control did not respond to H orientation")
    with out.open("w") as stream:
        table = csv.DictWriter(stream, rows[0].keys(), lineterminator="\n")
        table.writeheader()
        table.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("parent", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    run(a.parent.resolve(), a.candidate.resolve(), a.out)
