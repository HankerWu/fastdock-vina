"""Fixed top-mode torsion integration on the two external Factor Xa receptors."""

import argparse
import csv
import io
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from fxa_holdout import HERE, read, write


def evaluate(case, probe, shared, angles):
    row, work = case
    rec, lig, seed = (row[k] for k in ("receptor", "ligand", "seed"))
    text = (work / f"{lig}_{seed}.pdbqt").read_text()
    top, sep, _ = text.partition("ENDMDL")
    if not sep:
        raise ValueError("Incomplete docking pool")
    path = shared / f"{rec}_{lig}_{seed}.pdbqt"
    path.write_text(top + "ENDMDL\n")
    output = subprocess.run(
        [
            probe,
            work / f"{rec}.pdbqt",
            work / f"{lig}.pdbqt",
            path,
            str(angles),
            "298.15",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    (result,) = csv.DictReader(io.StringIO(output.stdout))
    terms = {k: float(result[k]) for k in ("inter", "intra", "work", "half_grid_work")}
    return {
        **row,
        **terms,
        "rotors": result["rotors"],
        "rotor": float(row["no_divisor"]) + terms["work"],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("probe", type=Path)
    p.add_argument("first", type=Path, help="2P16 work directory")
    p.add_argument("second", type=Path, help="1FJS work directory")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--angles", type=int, default=72)
    a = p.parse_args()
    labels = [
        (r, work)
        for rec, work in [("2P16", a.first), ("1FJS", a.second)]
        for r in read(HERE / f"results/fxa_{rec}_cells.csv")
    ]
    with tempfile.TemporaryDirectory() as shared, ThreadPoolExecutor(a.jobs) as pool:
        results = list(
            pool.map(
                lambda case: evaluate(case, a.probe.resolve(), Path(shared), a.angles),
                labels,
            )
        )
    suffix = "" if a.angles == 72 else f"_{a.angles}"
    write(HERE / f"results/fxa_rotor_cells{suffix}.csv", results)


if __name__ == "__main__":
    main()
