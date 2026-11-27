"""Independent CCD-graph control on three registered clinical inhibitor scaffolds."""

import argparse
import csv
import hashlib
import io
import subprocess
import sys
from pathlib import Path
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "paired_crossdocking"))
import run as base


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--parent", type=Path, required=True)
    p.add_argument("--fixed", type=Path, required=True)
    p.add_argument("--parent-graph", type=Path, required=True)
    p.add_argument("--fixed-graph", type=Path, required=True)
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    outputs = []
    for pdb, code in [
        ("6WTO", "3JW"),
        ("1XMU", "ROF"),
        ("3WZD", "LEV"),
        ("4LMN", "EUI"),
    ]:
        pdbfile, ccd = args.work / (pdb + ".pdb"), args.work / (code + ".cif")
        p_sha = base.fetch(pdbfile, "https://files.rcsb.org/download/" + pdb + ".pdb")
        c_sha = base.fetch(
            ccd, "https://files.rcsb.org/ligands/download/" + code + ".cif"
        )
        ligand = args.work / (pdb + ".pdbqt")
        data = base.ligand(
            ccd, base.dep_records(pdbfile.read_text().splitlines(), code), ligand
        )
        mol = Chem.MolFromSmiles(data["canonical_smiles"])
        serial = data["serial"]
        ref_bonds = {
            tuple(sorted((serial[b.GetBeginAtomIdx()], serial[b.GetEndAtomIdx()])))
            for b in mol.GetBonds()
        }
        chemical_dist = Chem.GetDistanceMatrix(mol)
        atom_index = {s: i for i, s in enumerate(serial)}
        for name, binary, graph in [
            ("parent", args.parent, args.parent_graph),
            ("fixed", args.fixed, args.fixed_graph),
        ]:
            edges = {"bond": set(), "pair": set()}
            for line in subprocess.check_output(
                [graph, ligand], text=True
            ).splitlines():
                kind, left, right = line.split()
                edges[kind].add((int(left), int(right)))
            short = sorted(
                (a, b)
                for a, b in edges["pair"]
                if chemical_dist[atom_index[a], atom_index[b]] <= 3
            )

            def labels(pairs):
                return " ".join(f"{a}-{b}" for a, b in sorted(pairs))

            def run(*args):
                text = subprocess.check_output([binary, *args, ligand], text=True)
                return list(csv.DictReader(io.StringIO(text), delimiter="\t"))

            pairs = run()
            forces = run("--forces")
            outputs.append(
                dict(
                    pdb=pdb,
                    ccd=code,
                    pdb_sha256=p_sha,
                    ccd_sha256=c_sha,
                    prepared_sha256=hashlib.sha256(ligand.read_bytes()).hexdigest(),
                    arm=name,
                    smiles=data["canonical_smiles"],
                    accepted_pairs=pairs[0]["pairs"],
                    spurious_pairs=pairs[0]["spurious"],
                    spurious_serials=" ".join(
                        row["serial_a"] + "-" + row["serial_b"]
                        for row in pairs
                        if row["serial_a"]
                    ),
                    missing_ccd_bonds=labels(ref_bonds - edges["bond"]),
                    extra_bonds=labels(edges["bond"] - ref_bonds),
                    spurious_on_ccd_graph=labels(short),
                    intraligand_energy=forces[0]["energy"],
                    spurious_energy=forces[0]["excess"],
                )
            )
    with (args.work / "results.csv").open("w") as f:
        w = csv.DictWriter(f, outputs[0], lineterminator="\n")
        w.writeheader()
        w.writerows(outputs)


if __name__ == "__main__":
    main()
