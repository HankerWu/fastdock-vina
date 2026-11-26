"""Regression checks for pool delimiters, error atomicity, strict topology and an empty interaction field."""

import argparse
import csv
import io
import subprocess
import tempfile
from pathlib import Path


def atom(x, y, z, atom_type="C"):
    return (
        f"ATOM      1  C   UNL     1    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00    +0.000 {atom_type}\n"
    )


def run(binary):
    with tempfile.TemporaryDirectory() as directory:
        base = Path(directory)
        rec = base / "receptor.pdbqt"
        lig = base / "ligand.pdbqt"
        pool = base / "pool.pdbqt"
        rec.write_text(atom(0, 0, 0))
        ligand = "ROOT\n" + atom(20, 0, 0) + "ENDROOT\nTORSDOF 0\n"
        lig.write_text(ligand)
        args = [
            str(binary),
            "--score_only",
            "--ligand",
            str(lig),
            "--receptor",
            str(rec),
            "--cpu",
            "1",
            "--verbosity",
            "0",
        ]
        flags = ["--rigid_window_pool", str(pool), "--rigid_window_samples", "16"]
        first = "MODEL\t3  \n" + ligand + "ENDMDL      \n"
        second = "MODEL 19\n" + ligand + "ENDMDL\n"
        pool.write_bytes((first + second).replace("\n", "\r\n").encode())

        def call(extra):
            return subprocess.run(
                args + extra, text=True, capture_output=True, check=False
            )

        ok = call(flags)
        if ok.returncode:
            raise AssertionError(ok.stderr)
        start = ok.stdout.index("model,inter,intra,")
        rows = list(csv.DictReader(io.StringIO(ok.stdout[start:])))
        assert [r["model"] for r in rows] == ["3", "19"]
        for r in rows:
            assert float(r["slope"]) == 1
            assert 0 < float(r["effective_samples"]) <= 48
            for name in (
                "inter",
                "intra",
                "window_delta",
                "half_delta",
                "uniform_delta",
            ):
                assert abs(float(r[name])) < 1e-12
        for bad in (
            "",
            first + first,
            first + "MODEL 2\n",
            first + "garbage\n",
            first + "MODEL 2\nENDMDL\n",
            first + second.replace("ROOT\n", "MODEL 20\n", 1),
            first + second.replace("ENDMDL", "ENDMDL extra"),
            first + second.replace("+0.000 C", "+0.000 N"),
            "MODEL -3\n" + ligand + "ENDMDL\n",
        ):
            pool.write_text(bad)
            result = call(flags)
            assert result.returncode != 0 and "model,inter,intra," not in result.stdout
        pool.write_text(second)
        for bad_flags in (
            flags[:-1] + ["17"],
            ["--rigid_window_samples", "16"],
            flags + ["--unbound_energy", "0"],
        ):
            result = call(bad_flags)
            assert result.returncode != 0 and "model,inter,intra," not in result.stdout
        plain = call([])
        assert plain.returncode == 0 and "model,inter,intra," not in plain.stdout
        print(
            "PASS: two spaced/CRLF model delimiters, degenerate one-atom far field, "
            "nine invalid pools with no partial rows, incompatible flags and ordinary score-only entry"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    run(parser.parse_args().binary.resolve())
