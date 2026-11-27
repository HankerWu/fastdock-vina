"""Measure energies and gradients on identical Cartesian coordinates and topology."""
import argparse
import csv
import io
import itertools
import subprocess
import tempfile
from pathlib import Path


def observations(binary, ligand):
    lines = ligand.read_text().splitlines(keepends=True)
    start, end = lines.index("ROOT\n") + 1, lines.index("ENDROOT\n")
    if end - start != 4:
        raise ValueError("Expected the four-membered baricitinib torsion-tree root")
    root_atoms = lines[start:end]
    with tempfile.TemporaryDirectory() as temp:
        variants = []
        for atoms in itertools.permutations(root_atoms):
            order = "".join(line[6:11].strip() for line in atoms)
            path = Path(temp) / (order + ".pdbqt")
            path.write_text("".join(lines[:start] + list(atoms) + lines[end:]))
            variants.append(path)
        def rows(*args):
            result = subprocess.run([binary, *args, *variants], capture_output=True, text=True, check=True)
            return list(csv.DictReader(io.StringIO(result.stdout), delimiter="\t"))
        counts = {row["path"]: row for row in rows()}
        data = rows("--forces")
    baseline = {r["serial"]: r for r in data if Path(r["path"]).stem == "1234"}
    for path in variants:
        current = [r for r in data if r["path"] == str(path)]
        count = counts[str(path)]
        def max_delta(fields):
            return max(abs(float(r[f]) - float(baseline[r["serial"]][f])) for r in current for f in fields)
        yield {
            "order": path.stem, "pairs": int(count["pairs"]), "spurious": int(count["spurious"]),
            "energy": float(current[0]["energy"]), "spurious_energy": float(current[0]["excess"]),
            "delta_energy": max_delta(["energy"]), "max_component_gradient_change": max_delta(["gx", "gy", "gz"]),
            "max_spurious_gradient_component": max(abs(float(r[f])) for r in current for f in ("excess_gx","excess_gy","excess_gz"))
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit_binary", type=Path)
    parser.add_argument("--ligand", type=Path, default=Path(__file__).with_name("baricitinib.pdbqt"))
    parser.add_argument("--require-valid", action="store_true")
    args = parser.parse_args()
    rows = list(observations(args.audit_binary.resolve(), args.ligand))
    import sys
    table = csv.DictWriter(sys.stdout, rows[0].keys(), lineterminator="\n")
    table.writeheader()
    table.writerows(rows)
    if args.require_valid:
        for r in rows:
            if (r["pairs"] != 202 or r["spurious"] or r["max_spurious_gradient_component"]
                or r["delta_energy"] > 1e-11 or r["max_component_gradient_change"] > 1e-11):
                raise ValueError(f"Graph-distance exclusion or atom-order invariance failed: {r}")


if __name__ == "__main__":
    main()
