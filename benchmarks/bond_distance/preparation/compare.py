"""Compare all ordered pair entries, then interleave parser timings by arm."""
import argparse
import csv
import hashlib
import itertools
import subprocess
import tempfile
from pathlib import Path


def csv_file(path, fields, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("parent", "bfs", "release", "macro_inputs", "repo", "out"):
        p.add_argument("--" + name.replace("_", "-"), required=True, type=Path)
    p.add_argument("--repeats", type=int, default=192)
    args = p.parse_args()
    macro = args.macro_inputs
    ligands = {
        "1IEP": args.repo / "examples/basic/1iep_ligand.pdbqt",
        "baricitinib": args.repo / "benchmarks/bond_distance/baricitinib.pdbqt",
        "5QCT": macro / "bace/raw/single-0.pdbqt",
        "1YET": macro / "hsp/1YET/ligand.pdbqt",
        "1FKB": macro / "extra/1FKB/ligand.pdbqt",
        "1FKF": macro / "extra/1FKF/ligand.pdbqt",
    }
    args.out.mkdir(parents=True, exist_ok=True)
    arms = {"parent": args.parent.resolve(), "bfs": args.bfs.resolve(), "release": args.release.resolve()}
    tables, times = [], []
    with tempfile.TemporaryDirectory() as temp:
        for code, ligand in ligands.items():
            content = ligand.read_text()
            variants = {"native": ligand}
            if "CG0" in content:
                for i in range(1,4):
                    path = Path(temp) / f"{code}-G{i}.pdbqt"
                    path.write_text(content.replace("G0", "G" + str(i)))
                    variants[f"G{i}"] = path
            for label, path in variants.items():
                for typing in ("XS", "AD"):
                    packed = {}
                    for name, binary in arms.items():
                        data = subprocess.check_output([binary, "table", typing, "1", path])
                        packed[name] = data
                        tables.append(dict(ligand=code, label=label, typing=typing, arm=name,
                                           pairs=data.count(b"\n"), sha256=hashlib.sha256(data).hexdigest()))
                    if packed["bfs"] != packed["release"]:
                        raise ValueError(f"Cached and direct lists differ: {code} {label} {typing}")
    # All six possible arm orders once per ligand, counterbalanced by rotating
    # the molecule order by round. Use process CPU as well as elapsed wall time.
    for block, order in enumerate(itertools.permutations(arms)):
        cases = list(ligands)
        for code in cases[block:] + cases[:block]:
            for position, name in enumerate(order):
                output = subprocess.check_output([arms[name], "time", "XS", str(args.repeats), ligands[code]], text=True)
                count, cpu, wall = output.split()
                times.append(dict(block=block + 1, position=position + 1, ligand=code, arm=name,
                                  parses=args.repeats, sum_internal_pairs=int(count),
                                  cpu_us=1e6*float(cpu), wall_us=1e6*float(wall)))
    csv_file(args.out / "tables.csv", tables[0], tables)
    csv_file(args.out / "timing.csv", times[0], times)


if __name__ == "__main__":
    main()
