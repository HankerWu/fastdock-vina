"""Confirm byte parity and invariants on the global docking entry point."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main(binary, baseline):
    root = Path(__file__).resolve().parents[2]
    case = root / "examples/basic"
    ligand = case / "1iep_ligand.pdbqt"

    def dock(exe, out, lig, flags=(), cpu=1, n=4):
        return subprocess.run(
            [
                exe,
                "--receptor",
                case / "1iep_receptor.pdbqt",
                "--ligand",
                lig,
                "--config",
                case / "1iep_receptor.box.txt",
                "--pose_selection",
                "late",
                "--num_modes",
                "6",
                "--exhaustiveness",
                str(n),
                "--max_evals",
                "5000",
                "--energy_range",
                "100",
                "--seed",
                "59311",
                "--cpu",
                str(cpu),
                "--search_counters",
                "--out",
                out,
                *flags,
            ],
            capture_output=True,
            text=True,
        )

    with tempfile.TemporaryDirectory() as temp:
        out = Path(temp) / "poses.pdbqt"

        def run(exe, lig, flags=(), cpu=1):
            p = dock(exe, out, lig, flags, cpu)
            p.check_returncode()
            tracks = p.stdout.splitlines()
            return out.read_bytes(), [
                l for l in tracks if l.startswith("Search trajectory:")
            ]

        ng = run(baseline, ligand)
        native = run(binary, ligand)
        assert ng == native
        for flag in ("--prepared_torsion_starts", "--balanced_torsion_search"):
            changed = run(binary, ligand, [flag])
            assert changed == run(binary, ligand, [flag], 2)
            assert changed[0] != native[0]
            if flag == "--prepared_torsion_starts":
                assert len(changed[1]) == 4 and changed[1][::2] == native[1][::2]
                assert dock(binary, out, ligand, [flag], n=1).returncode != 0
        # The rigid root alone is a real, fixed heavy-atom fragment. Drop
        # atom-mapping remarks as the software invariant only requires coordinates.
        fragment = Path(temp) / "rigid.pdbqt"
        original = ligand.read_text().splitlines()
        fragment.write_text(
            "\n".join(original[original.index("ROOT") : original.index("ENDROOT") + 1])
            + "\nTORSDOF 0\n"
        )
        rigid = run(binary, fragment)
        for flag in ("--prepared_torsion_starts", "--balanced_torsion_search"):
            assert run(binary, fragment, [flag]) == rigid
        print(
            "Default byte parity, threads 1/2, odd trajectories, rigid-fragment null and one-trajectory rejection passed."
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("baseline", type=Path)
    a = p.parse_args()
    main(a.binary.resolve(), a.baseline.resolve())
