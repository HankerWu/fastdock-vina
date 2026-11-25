"""Verify fixed-topology scoring on one drifting and one unchanged real pool."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from rdkit import Chem

TOOL = Path(__file__).resolve().parents[1] / "tools/score_physical_pool.py"


def run(arguments, output, *, expected_error=None):
    process = subprocess.run([sys.executable, str(TOOL), *map(str, arguments),
                              "--output", str(output)], text=True, capture_output=True, timeout=60)
    if expected_error is not None:
        assert process.returncode and expected_error in process.stderr and not output.exists(), process.stderr
        return None
    assert process.returncode == 0, process.stderr
    return json.loads(output.read_text()) if output.suffix == ".json" else output.read_text()


def fixture(label, binary, prepared, receptor, pool, work):
    arguments = ["--pool", pool, "--binary", binary, "--receptor", receptor]
    export = ["--prepared-ligand", prepared, "--export-ligand-bonds"]
    bonds = work / f"{label}.txt"
    graph = run(export, bonds)
    edges = [tuple(map(int, line.split())) for line in graph.splitlines()]
    assert edges == sorted(set(edges))
    default = run(arguments, work / f"{label}-default.json")
    stable = run([*arguments, "--prepared-ligand", prepared, "--ligand-bonds", bonds],
                 work / f"{label}-stable.json")
    assert "fixed_topology" not in default and "saved_native" not in default["poses"][0]
    assert stable["fixed_topology"]["saved_native_order"] == default["orders"]["native"]
    assert [row["saved_native"] for row in stable["poses"]] == [row["native"] for row in default["poses"]]
    for row in stable["poses"]:
        engine = row["engine"]
        expect = engine["zero_reference"] - engine["scale"]*stable["fixed_topology"]["unbound"]
        assert abs(row["native"] - expect) < 1e-10
    if label == "9FA3":
        assert len(edges) == 16 and stable["orders"]["native"][0] == 3
    else:
        for first, second in zip(default["poses"], stable["poses"], strict=True):
            for name in ("inter", "intra", "zero_reference", "scale"):
                assert abs(first["engine"][name] - second["engine"][name]) < 1e-8
    return stable, arguments, bonds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--bound", type=Path, nargs=3, required=True, metavar=("PREPARED", "RECEPTOR", "POOL"),
                        help="Broad-validation 9FA3 seed 42; the expected reordered model is fixture-specific")
    parser.add_argument("--unchanged", type=Path, nargs=3, required=True, metavar=("PREPARED", "RECEPTOR", "POOL"),
                        help="Post-reveal T4 HIE benzene seed 43, which has no inferred bond changes")
    parser.add_argument("--box", required=True, type=Path)
    parser.add_argument("--public", required=True, type=Path)
    args = parser.parse_args()
    with TemporaryDirectory(prefix="fixed-topology-check-") as temporary:
        work = Path(temporary)
        stable, command, bonds = fixture("9FA3", args.binary, *args.bound, work)
        fixture("unchanged", args.binary, *args.unchanged, work)
        prepared, receptor, pool = args.bound
        exported = work / "same-pool.sdf"
        run(["--pool", pool, "--export-gnina-sdf"], exported)
        scored = work / "same-pool-cnn.sdf"
        with Chem.SDWriter(str(scored)) as writer:
            for molecule in Chem.SDMolSupplier(str(exported), sanitize=False, removeHs=False):
                molecule.SetDoubleProp("CNNscore", 0.5)
                writer.write(molecule)
        baseline = run([*command, "--gnina-scores-sdf", scored], work / "legacy-cnn.json")
        assert baseline["orders"]["cnn_consensus"] == baseline["orders"]["native"]
        cnn = run([*command, "--prepared-ligand", prepared, "--ligand-bonds", bonds,
                   "--gnina-scores-sdf", scored], work / "stable-cnn.json")
        assert cnn["orders"]["cnn_consensus"] == stable["orders"]["native"]
        assert [row["cnn_consensus"] for row in cnn["poses"]] == [row["native"] for row in cnn["poses"]]
        vinardo = ["--pool", pool, "--prepared-ligand", prepared, "--ligand-bonds", bonds,
                   "--receptor", receptor, "--vinardo-box", args.box, "--vinardo-binary"]
        new = run([*vinardo, args.binary], work / "vinardo-stable.json")
        assert not new["vinardo_scoring"]["failures"] and len(new["orders"]["vinardo"]) == 9
        old = run([*vinardo, args.public], work / "vinardo-unsupported.json")
        assert len(old["vinardo_scoring"]["failures"]) == 9
        assert old["orders"]["vinardo"] == old["orders"]["native"] and "vinardo" in old["abstentions"]
        broken = work / "bad-bonds.txt"
        broken.write_text("1 999\n")
        run([*command, "--prepared-ligand", prepared, "--ligand-bonds", broken], work / "bad.json",
            expected_error="differs from the prepared Meeko graph")
        bad_pose = work / "bad-ligand.pdbqt"
        bad_pose.write_text(prepared.read_text().replace("REMARK SMILES IDX 4 1", "REMARK SMILES IDX 2 1", 1))
        run(["--prepared-ligand", bad_pose, "--export-ligand-bonds"], work / "broken.txt",
            expected_error="must map each encoded atom one-to-one")
        run([*command, "--prepared-ligand", bad_pose], work / "bad-pool.json",
            expected_error="differs from the prepared PDBQT outside coordinates")
    print("9FA3, unchanged topology, bond validation, and Vinardo fallback checks passed")


if __name__ == "__main__":
    main()
