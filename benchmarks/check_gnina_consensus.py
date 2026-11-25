#!/usr/bin/env python3
"""Validate the GNINA input contract using one existing real scored pool."""

import argparse
import json
from math import isclose, log
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from check_late_pose_selection import ROOT, require

sys.path.insert(0, str(ROOT / "tools"))
from score_physical_pool import export_gnina_sdf, score_gnina_pool


def check_macrocycle_export(work):
    # Meeko opens the 12-membered ring for docking and marks two real
    # carbons and two glue sites. Its hydrated ligand also contains W sites.
    molecule = Chem.AddHs(Chem.MolFromSmiles("OC1CCCCCCCCCCC1"))
    require(AllChem.EmbedMolecule(molecule, randomSeed=42) == 0, "Could not embed cyclododecanol")
    [prepared] = MoleculePreparation(hydrate=True).prepare(molecule)
    pdbqt, success, message = PDBQTWriterLegacy.write_string(prepared)
    require(success, message)
    types = [line.split()[-1] for line in pdbqt.splitlines() if line.startswith("ATOM ")]
    require(types.count("CG0") == types.count("G0") == 2 and "W" in types,
            "The real macrocycle did not exercise both pseudo-atom kinds")

    def block(model, ligand):
        return f"MODEL {model}\nREMARK VINA RESULT: -5 0 0\n{ligand}ENDMDL\n"

    pool = work / "macrocycle.pdbqt"
    destination = work / "macrocycle.sdf"
    cli = [sys.executable, str(ROOT / "tools/score_physical_pool.py"), "--export-gnina-sdf",
           "--pool", str(pool), "--output", str(destination)]
    pool.write_text(block(1, pdbqt) + block(2, pdbqt.replace(" CG0\n", " OA\n", 1)))
    run = subprocess.run(cli, capture_output=True, text=True, timeout=30)
    require(run.returncode != 0 and "heavy-atom elements or coordinates differ for MODEL 2" in run.stderr,
            run.stderr or "Element mismatch was accepted for the second MODEL")
    require(not destination.exists() and not list(work.glob(".gnina-export-*")),
            "Failed export left a new output or staging directory")

    pool.write_text(block(1, pdbqt) + block(2, pdbqt))
    subprocess.run(cli, capture_output=True, text=True, timeout=30, check=True)
    expected = destination.read_bytes()
    exported = list(Chem.SDMolSupplier(str(destination), removeHs=False))
    require([m.GetIntProp("native_model") for m in exported] == [1, 2],
            "Wrong macrocycle MODEL metadata")
    for mol in exported:
        observed = sorted(a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
        require(observed == [6] * 12 + [8], "Exported physical macrocycle has the wrong elements")
    try:
        export_gnina_sdf(pool, destination)
    except FileExistsError:
        require(destination.read_bytes() == expected, "The existing SDF was modified")
    else:
        raise RuntimeError("Export accepted an existing destination")


def check(pool, sdf, work):
    output = work / "scores.json"
    subprocess.run([sys.executable, str(ROOT / "tools/score_physical_pool.py"),
                    "--pool", str(pool), "--gnina-scores-sdf", str(sdf), "--output", str(output)],
                   capture_output=True, text=True, timeout=30, check=True)
    result = json.loads(output.read_text())
    molecules = list(Chem.SDMolSupplier(str(sdf), sanitize=False, removeHs=False))
    cnn_scores = {int(mol.GetProp("native_model")): float(mol.GetProp("CNNscore")) for mol in molecules}
    require(len(molecules) >= 3 and len(molecules) == len(result["poses"]), "Expected a complete multi-pose fixture")
    for row in result["poses"]:
        cnn = cnn_scores[row["model"]]
        c = max(1e-6, min(1 - 1e-6, cnn))
        require(cnn == row["gnina"] and isclose(row["cnn_consensus"], row["native"] - 0.592 * log(c / (1 - c))),
                "Wrong raw or consensus score")
    for arm, sign in (("gnina", -1), ("cnn_consensus", 1)):
        expected_order = sorted(result["poses"], key=lambda row: (sign * row[arm], row["native_rank"]))
        require(result["orders"][arm] == [row["model"] for row in expected_order], f"Wrong order for {arm}")

    modified = work / "modified.sdf"

    def write(records):
        with Chem.SDWriter(str(modified)) as writer:
            for record in records:
                writer.write(record)

    def rejects(records, message):
        write(records)
        try:
            score_gnina_pool(pool, modified)
        except ValueError as error:
            require(message in str(error), str(error))
        else:
            raise RuntimeError(f"Invalid GNINA input accepted ({message})")

    first, *others = molecules
    rejects(others, "missing native_model")
    rejects([first, first, *others], "duplicate")
    for invalid in ("nan", "1.02"):
        changed_score = Chem.Mol(first)
        changed_score.SetProp("CNNscore", invalid)
        rejects([changed_score, *others], "outside [0, 1]")

    moved = Chem.Mol(first)
    heavy_index = next(a.GetIdx() for a in moved.GetAtoms() if a.GetAtomicNum() > 1)
    positions = moved.GetConformer()
    point = positions.GetAtomPosition(heavy_index)
    point.x += 0.0001
    positions.SetAtomPosition(heavy_index, point)
    rejects([moved, *others], "elements or coordinates differ")
    changed_element = Chem.Mol(first)
    atom = changed_element.GetAtomWithIdx(heavy_index)
    atom.SetAtomicNum(8 if atom.GetAtomicNum() != 8 else 6)
    rejects([changed_element, *others], "elements or coordinates differ")

    # Equal raw probabilities use native ranks even if the SDF records run
    # backwards. Exact endpoint values remain available in the raw score.
    for i, molecule in enumerate(molecules):
        molecule.SetProp("CNNscore", "0" if i == 0 else "1")
    write(reversed(molecules))
    result = score_gnina_pool(pool, modified)
    lowest, *tied = result["poses"]
    require(lowest["gnina"] == 0 and all(row["gnina"] == 1 for row in tied), "Raw endpoints mutated")
    ordered = sorted(tied, key=lambda row: row["native_rank"]) + [lowest]
    require(result["orders"]["gnina"] == [row["model"] for row in ordered], "Ties ignored native ranks")
    for row in result["poses"]:
        c = 1 - 1e-6 if row["gnina"] else 1e-6
        require(isclose(row["cnn_consensus"], row["native"] - 0.592 * log(c / (1 - c))),
                "Endpoint clipping differs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pool", type=Path)
    parser.add_argument("gnina_scores", type=Path)
    args = parser.parse_args()
    with TemporaryDirectory(prefix="vina-gnina-check-") as directory:
        work = Path(directory)
        check(args.pool, args.gnina_scores, work)
        check_macrocycle_export(work)
    print("GNINA acceptance passed: fixed-coordinate ranking, malformed input rejection and macrocycle export.")


if __name__ == "__main__":
    main()
