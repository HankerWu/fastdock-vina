"""Reproduce the native/mixed/soft comparison from the frozen OpenFF archive and binaries."""

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
from urllib.request import urlopen

from rdkit import Chem

SHA = "fd88824f9114244f95a14b485e6d6c96c1de716d"
URL = (
    f"https://raw.githubusercontent.com/openforcefield/protein-ligand-benchmark/{SHA}/"
)


def load_prior(root):
    spec = importlib.util.spec_from_file_location(
        "paired_panel", root / "benchmarks/mixed_vinardo/paired_panel.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("work", "prior_source", "mixed_vina", "soft_vina"):
        p.add_argument(arg, type=Path)
    args = p.parse_args()
    prior = load_prior(args.prior_source)
    work = args.work.resolve()
    inputs = work / "openff"
    hashes = {}
    expected = json.loads(Path(__file__).with_name("input_hashes.json").read_text())
    for target, anchor in prior.AFFINITY.items():
        root = work / target
        root.mkdir(parents=True, exist_ok=True)
        for f in (
            "00_data/ligands.yml",
            "01_protein/crd/protein.pdb",
            "02_ligands/ligands.sdf",
        ):
            rel = "data/" + target + "/" + f
            path = inputs / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                path.write_bytes(urlopen(URL + rel, timeout=60).read())
            hashes[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
            if expected.get(rel) != hashes[rel]:
                raise ValueError("Input checksum differs: " + rel)
        a = SimpleNamespace(
            binary=args.mixed_vina.resolve(),
            data=inputs,
            output=root,
            subset=target,
            exhaustiveness=4,
            cap=3000,
        )
        prior.run_affinity(a)
        rows = list(csv.DictReader((root / "affinity_rows.csv").open()))
        (ref,) = [
            m
            for m in Chem.SDMolSupplier(
                str(inputs / "data" / target / "02_ligands/ligands.sdf"), removeHs=False
            )
            if m and m.GetProp("_Name") == anchor
        ]
        prior.ARMS["soft"] = ("--soft_receptor_search",)
        try:
            for row in list(rows):
                if row["arm"] != "native":
                    continue
                if row["status"] != "ok":
                    rows.append(dict(row, arm="soft"))
                    continue
                folder = root / target / row["ligand"]
                _, _, metrics = prior.docking(
                    args.soft_vina.resolve(),
                    root / target / "receptor.pdbqt",
                    folder / "ligand.pdbqt",
                    prior.coords(ref).mean(axis=0),
                    folder,
                    int(row["seed"]),
                    "soft",
                    4,
                    3000,
                )
                rows.append(dict(row, arm="soft", **metrics))
        finally:
            del prior.ARMS["soft"]
        prior.write_rows(root / "three_arms.csv", rows)
    (work / "input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main()
