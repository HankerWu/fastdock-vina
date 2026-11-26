"""Export the crystallographic query in each receiver frame for post-docking evaluation."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from paired import aligned_receivers
from rdkit import Chem


def main(repo, work):
    sys.path.insert(0, str(repo.resolve() / "benchmarks/paired_crossdocking"))
    import run

    metadata = json.loads((work / "metadata.json").read_text())
    directory = work / "frames"
    directory.mkdir(exist_ok=True)
    for query, source in metadata.items():
        for receptor, target in metadata.items():
            if source["target"] != target["target"]:
                continue
            rotation, translation, _, _ = run.alignment(
                *aligned_receivers(source, target)
            )
            coordinates = np.array(source["mapped"][0]) @ rotation + translation
            mol = Chem.MolFromSmiles(source["canonical_smiles"])
            conformer = Chem.Conformer(mol.GetNumAtoms())
            for atom, point in enumerate(coordinates):
                conformer.SetAtomPosition(atom, point)
            mol.AddConformer(conformer)
            with Chem.SDWriter(str(directory / f"{query}_{receptor}.sdf")) as writer:
                writer.write(mol)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    main(args.repo, args.work)
