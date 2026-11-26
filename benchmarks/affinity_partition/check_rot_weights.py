"""Exercise Vina/Vinardo torsion weights through CLI and config-file entry points."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

BASIC = Path(__file__).resolve().parents[2] / "examples/basic"


def score(binary, mode, *options):
    args = [
        binary,
        "--scoring",
        mode,
        "--score_only",
        "--cpu",
        "1",
        "--receptor",
        BASIC / "1iep_receptor.pdbqt",
        "--ligand",
        BASIC / "1iep_ligand.pdbqt",
        *options,
    ]
    run = subprocess.run(args, check=True, text=True, capture_output=True)

    def term(label):
        return float(
            re.search(re.escape(label) + r"\s*:\s*([-\d.]+)", run.stdout).group(1)
        )

    return (
        term("Estimated Free Energy of Binding"),
        term("Final Intermolecular Energy"),
        term("Torsional Free Energy"),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument(
        "--historical", type=Path, help="Compare defaults to a pre-fix binary"
    )
    a = p.parse_args()
    binary = a.binary.resolve()
    for mode, active, other in (
        ("vina", "weight_rot", "weight_vinardo_rot"),
        ("vinardo", "weight_vinardo_rot", "weight_rot"),
    ):
        default = score(binary, mode)
        disabled = score(binary, mode, "--" + active, "0")
        assert disabled[2] == 0 and disabled[0] == disabled[1]
        assert default[1] == disabled[1] and default[2] > 0
        assert score(binary, mode, "--" + other, "0") == default
        assert score(binary, mode, "--" + active, "0", "--" + other, "0.2") == disabled
        with tempfile.NamedTemporaryFile(mode="w") as config:
            config.write(active + " = 0\n")
            config.flush()
            assert score(binary, mode, "--config", config.name) == disabled
        if a.historical:
            assert score(a.historical.resolve(), mode) == default
        print(mode, "default", default, "zero torsion", disabled)


if __name__ == "__main__":
    main()
