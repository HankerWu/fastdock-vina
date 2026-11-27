"""Replay the published A2A dry branch and water-network differences against SPR Kd."""
import argparse
import csv
import math
from pathlib import Path

RT = 0.00198720425864083 * 298.15


def main(source, output):
    with source.open() as f:
        rows = list(csv.DictReader(f))
    results = []
    for row in rows:
        ratio = float(row["to_kd_M"]) / float(row["from_kd_M"])
        measured = RT * math.log(ratio)
        dry = float(row["dry_kcal_mol"])
        corrected = dry + float(row["network_difference_kcal_mol"])
        results.append(dict(from_ligand=row["from_ligand"], to_ligand=row["to_ligand"],
                            kd_ratio=ratio, measured_kcal_mol=measured,
                            dry_kcal_mol=dry, cycle_estimate_kcal_mol=corrected,
                            dry_abs_error=abs(dry - measured),
                            cycle_abs_error=abs(corrected - measured)))
    with output.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=results[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.source, args.output)
