"""Check a serialized pose with the fixed-coordinate CLI and reuse its annotated ligand."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREFIX = "REMARK VINA CONTACT ROTOR SCORE:"


def test(binary):
    checked = []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        with tarfile.open(
            HERE.parents[1] / "geometry_comparison/affinity/inputs.tar.xz"
        ) as inputs:
            for name in (
                "4DJW.box",
                "4DJW.pdbqt",
                "CHEMBL2151138.pdbqt",
                "CHEMBL2151136.pdbqt",
            ):
                (work / name).write_bytes(inputs.extractfile("prepared/" + name).read())
        with tarfile.open(HERE / "outputs.tar.xz") as pools:
            for ligand in ("CHEMBL2151138", "CHEMBL2151136"):
                key = f"pools/4DJW_{ligand}_96431.pdbqt"
                block = pools.extractfile(key).read().decode().split("ENDMDL")[0]
                block = (
                    "\n".join(
                        line
                        for line in block.splitlines()
                        if not line.startswith("MODEL")
                    )
                    + "\n"
                )
                (work / "fixed.pdbqt").write_text(block)
                before = next(
                    line[len(PREFIX) :].split()
                    for line in block.splitlines()
                    if line.startswith(PREFIX)
                )
                unbound = next(
                    line.split()[-1]
                    for line in block.splitlines()
                    if line.startswith("REMARK UNBOUND:")
                )
                cmd = [
                    str(binary),
                    "--config",
                    "4DJW.box",
                    "--receptor",
                    "4DJW.pdbqt",
                    "--ligand",
                    ligand + ".pdbqt",
                    "--score_only",
                    "--pose_coordinates",
                    "fixed.pdbqt",
                    "--unbound_energy",
                    unbound,
                    "--contact_rotor_score",
                ]
                out = subprocess.run(
                    cmd, cwd=work, check=True, text=True, capture_output=True
                ).stdout
                after = next(
                    line[len(PREFIX) :].split()
                    for line in out.splitlines()
                    if line.startswith(PREFIX)
                )
                assert (
                    before[1:] == after[1:]
                    and abs(float(before[0]) - float(after[0])) < 0.02
                )
                # A model split from an annotated file does not carry stale diagnostics into a new native output.
                cmd = [
                    str(binary),
                    "--config",
                    "4DJW.box",
                    "--receptor",
                    "4DJW.pdbqt",
                    "--ligand",
                    "fixed.pdbqt",
                    "--local_only",
                    "--out",
                    "native.pdbqt",
                    "--cpu",
                    "1",
                ]
                subprocess.run(cmd, cwd=work, check=True, capture_output=True)
                assert PREFIX not in (work / "native.pdbqt").read_text()
                checked.append({"ligand": ligand, "stored": before, "fixed": after})
    print(json.dumps(checked, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    test(p.parse_args().binary.resolve())
