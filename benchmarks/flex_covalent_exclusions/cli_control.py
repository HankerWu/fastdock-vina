"""Score a bonded backbone atom and a four-bond neighbor on identical pose geometries."""
import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
FLEX = (HERE.parent / "fixtures/flex-residue.pdbqt").read_text()
LIGAND = """ROOT
ATOM      1    C UNL     1      11.100  50.700  15.300  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("fixed", type=Path)
    ap.add_argument("csv", type=Path)
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as temporary, args.csv.open("w", newline="") as out:
        temp = Path(temporary)
        ligand, flex = temp / "ligand.pdbqt", temp / "flex.pdbqt"
        ligand.write_text(LIGAND)
        rows = csv.writer(out)
        rows.writerow(("arm", "field", "angle_rad", "flex_both", "flex_four_bonds_only",
                       "flex_difference", "ligand_difference"))
        for angle in (0, .71, -1.29):
            xyz = "".join(f"{v:8.3f}" for v in (12, 50 + 1.5 * math.cos(angle), 14 + 1.5 * math.sin(angle)))
            flex.write_text("\n".join(line[:30] + xyz + line[54:] if line.startswith("ATOM") and
                                      line[12:16].strip() == "CG" else line for line in FLEX.splitlines()) + "\n")
            for name, binary in (("baseline", args.baseline), ("fixed", args.fixed)):
                for field in ("explicit", "grid"):
                    values = []
                    for rigid in ("covalent", "four_bonds"):
                        cmd = [binary.resolve(), "--ligand", ligand, "--receptor", HERE / f"{rigid}.pdbqt",
                               "--flex", flex, "--score_only", "--center_x", "11", "--center_y", "50", "--center_z",
                               "14", "--size_x", "10", "--size_y", "10", "--size_z", "10", "--spacing", ".3125",
                               "--seed", "921", "--cpu", "1"]
                        if field == "grid":
                            cmd.append("--no_refine")
                        log = subprocess.check_output(cmd, text=True)
                        values.append([float(re.search(pattern + r"\s*:\s*([-+0-9.]+)", log)[1])
                                       for pattern in (r"Flex\s+- Receptor", r"Ligand - Receptor")])
                    rows.writerow((name, field, angle, values[0][0], values[1][0],
                                   round(values[0][0] - values[1][0], 3), round(values[0][1] - values[1][1], 3)))


if __name__ == "__main__":
    main()
