"""Apply both cutoff treatments to identical saved poses, preserving ligand topology."""

import argparse
import csv
import itertools
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory


def rescore(binary, prepared, pool, receptor, ligand, seed):
    # Use this existing diagnostic for its full-precision ordinary inter/intra
    # terms. None of the polar corrections enter the energy below.
    common = [
        str(binary),
        "--receptor",
        str(prepared / f"{receptor}_receptor.pdbqt"),
        "--ligand",
        str(prepared / f"{ligand}_ligand.pdbqt"),
        "--score_only",
        "--net_polar",
        "native",
        "--cpu",
        "1",
        "--verbosity",
        "0",
    ]
    rows = []
    blocks = re.findall(r"^MODEL +\d+\s*\n.*?^ENDMDL$", pool.read_text(), re.M | re.S)
    if not blocks:
        raise ValueError(f"No complete models in {pool}")
    with TemporaryDirectory() as tmp:
        coordinates = Path(tmp) / "fixed.pdbqt"
        for number, block in enumerate(blocks, 1):
            coordinates.write_text(block + "\n")
            row = dict(receptor=receptor, ligand=ligand, seed=seed, model=number)
            for arm, flags in (("native", []), ("smooth", ["--smooth_cutoff"])):
                report = subprocess.check_output(
                    common + ["--pose_coordinates", str(coordinates)] + flags, text=True
                )
                inter, intra, zero = map(
                    float,
                    re.search(r"^Net polar native: (.+)$", report, re.M)[1].split(),
                )
                scale = float(
                    re.search(r"^Net polar normalization: (.+)$", report, re.M)[1]
                )
                # 'zero' selects the best pose within a ligand (shared reference).
                # Its own intra energy is then the reference for its affinity.
                row.update(
                    {
                        arm + "_inter": inter,
                        arm + "_intra": intra,
                        arm + "_zero": zero,
                        arm + "_top_affinity": round(scale * inter, 9),
                    }
                )
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("series", choices=["t4l", "hcaii"])
    parser.add_argument(
        "work", type=Path, help="Series directory containing prepared/ and runs/"
    )
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    inputs = args.work / "prepared"
    ligands = sorted(p.name.split("_")[0] for p in inputs.glob("*_ligand.pdbqt"))
    receptors = ("1L90", "186L") if args.series == "t4l" else ("6WQ4", "6WQ8")
    expected = 9 if args.series == "t4l" else 5
    if len(ligands) != expected:
        raise ValueError(f"Expected {expected} prepared ligands, found {len(ligands)}")

    def run(cell):
        rec, lig, seed = cell
        stem = f"{lig}_{rec}" if args.series == "t4l" else f"{rec}_{lig}"
        pool = args.work / "runs" / f"{stem}_{seed}_native.pdbqt"
        return rescore(args.binary.resolve(), inputs, pool, rec, lig, seed)

    cells = itertools.product(receptors, ligands, (42, 1729))
    with ThreadPoolExecutor(args.workers) as executor:
        rows = list(itertools.chain.from_iterable(executor.map(run, cells)))
    with args.output.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
