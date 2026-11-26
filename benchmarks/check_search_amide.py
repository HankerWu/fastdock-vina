"""Check search guidance and final amide ranking on 6WQ8 in receptor 6WQ4."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path, help="CLI from the parent commit")
    parser.add_argument(
        "fixtures", type=Path, help="extracted amide_search fixture directory"
    )
    args = parser.parse_args()
    prepared = args.fixtures / "fdv-second/prepared"
    inputs = [
        "--receptor",
        prepared / "6WQ4_receptor.pdbqt",
        "--config",
        prepared / "6WQ4_box.txt",
        "--ligand",
        prepared / "6WQ8_ligand.pdbqt",
    ]
    opts = [
        "--seed",
        "314159",
        "--exhaustiveness",
        "2",
        "--max_evals",
        "5500",
        "--num_modes",
        "4",
        "--search_counters",
    ]
    anchor = [
        "--search_anchor",
        args.fixtures / "fdv-map-states-declared/inputs/6WQ4_6WQ8.anchor",
    ]
    late = ["--pose_selection", "late"]
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        out, site, pose = (
            directory / name for name in ("pool.pdbqt", "amide.txt", "pose.pdbqt")
        )
        amide = ["--search_amide", site]

        def dock(
            *flags, binary=args.binary, cpu=1, window=100, ligand=None, table=False
        ):
            proc = subprocess.run(
                [
                    binary,
                    *inputs[:-1],
                    ligand or inputs[-1],
                    *opts,
                    *flags,
                    "--cpu",
                    str(cpu),
                    "--energy_range",
                    str(window),
                    "--out",
                    out,
                ],
                capture_output=True,
                check=True,
            )
            traces = re.findall(rb"^Search trajectory: .+$", proc.stdout, re.MULTILINE)
            assert len(traces) == 2
            return (
                (out.read_bytes(), proc.stdout) if table else (out.read_bytes(), traces)
            )

        for flags in ([], late, anchor):
            assert dock(*flags) == dock(*flags, binary=args.baseline)
        site.write_text("2 1 3 4 17 0\n")
        assert dock(*amide) == dock(*late)
        assert dock(*anchor, *amide) == dock(*anchor)
        site.write_text("2 1 3 4 17 5\n")
        amide_only = dock(*amide)
        combined = dock(*amide, *anchor)
        assert amide_only[0] != dock(*late)[0]
        assert combined[0] != dock(*anchor)[0]
        assert amide_only == dock(*amide, cpu=2)
        assert combined == dock(*amide, *anchor, cpu=2)
        assert dock(*amide, *anchor, "--search_anchor_trajectories", "even")[1] == [
            amide_only[1][0],
            combined[1][1],
        ]
        assert dock(*amide, *anchor, "--search_anchor_trajectories", "odd")[1] == [
            combined[1][0],
            amide_only[1][1],
        ]
        rejected = subprocess.run(
            [args.binary, *inputs, "--score_only", *amide, "--cpu", "1"],
            capture_output=True,
            check=False,
        )
        assert rejected.returncode and b"search_amide" in rejected.stderr

        models = re.findall(
            rb"^MODEL .*?^ENDMDL$",
            amide_only[0] + combined[0],
            re.MULTILINE | re.DOTALL,
        )
        assert models
        worst = 0.0

        def remark(model, field):
            return float(
                re.search(rb"^REMARK " + field + rb":\s+(\S+)", model, re.MULTILINE)[1]
            )

        for model in models:
            pose.write_bytes(model + b"\n")
            score = subprocess.check_output(
                [
                    args.binary,
                    *inputs,
                    "--score_only",
                    "--pose_coordinates",
                    pose,
                    "--unbound_energy",
                    str(remark(model, b"UNBOUND")),
                    "--cpu",
                    "1",
                ],
                text=True,
                stderr=subprocess.PIPE,
            )
            for label, value in (
                ("Ligand - Receptor", remark(model, b"INTER")),
                ("Ligand", remark(model, b"INTRA")),
                ("Estimated Free Energy of Binding", remark(model, b"VINA RESULT")),
            ):
                measured = float(
                    re.search(
                        r"^\s*" + label + r"\s*:\s+([-+0-9.]+)", score, re.MULTILINE
                    )[1]
                )
                worst = max(worst, abs(measured - value))
        # Written atom coordinates are rounded to 0.001 Angstrom.
        assert worst < 0.01, worst
        # Final ranking must preserve the entire searched and refined pool, its
        # native reference and original full-precision affinities.
        rank = ["--rank_amide", site]
        site.write_text("2 1 3 4 17 0\n")
        assert dock(*rank) == dock(*late)
        site.write_text("2 1 3 4 17 15\n")

        def blocks(content):
            return re.findall(rb"^MODEL .*?^ENDMDL$", content, re.MULTILINE | re.DOTALL)

        def numbers(block, label):
            return [
                float(s)
                for s in re.search(
                    rb"^REMARK VINA " + label + rb": (.*)$", block, re.MULTILINE
                )[1].split()
            ]

        def atoms(block):
            return [
                line
                for line in block.splitlines()
                if line.startswith((b"ATOM ", b"HETATM "))
            ]

        ranked, table = dock(*anchor, *rank, table=True)
        originals = blocks(dock(*anchor)[0])
        reordered = blocks(ranked)
        assert len(reordered) == len(originals)
        new_scores, native_ranks = [], []
        for model in reordered:
            energy, native_rank = numbers(model, b"NATIVE")
            native_ranks.append(int(native_rank))
            original = originals[native_ranks[-1] - 1]
            assert atoms(model) == atoms(original)
            assert numbers(model, b"RESULT")[0] == numbers(original, b"RESULT")[0]
            assert abs(numbers(model, b"RESULT")[0] - energy) <= 0.000501
            assert re.findall(
                rb"^REMARK (?:INTER|INTRA|UNBOUND).*?$", model, re.MULTILINE
            ) == re.findall(
                rb"^REMARK (?:INTER|INTRA|UNBOUND).*?$", original, re.MULTILINE
            )
            value = numbers(model, b"AMIDE RESULT")[0]
            (amount,) = numbers(model, b"AMIDE PENALTY")
            assert amount >= 0 and abs(value - energy - amount) < 1e-12
            new_scores.append(value)
        assert new_scores == sorted(new_scores)
        assert sorted(native_ranks) == list(range(1, len(originals) + 1))
        assert native_ranks[0] != 1, (
            "fixture must distinguish native- and amide-ranked winners"
        )
        assert b"amide score" in table
        table_scores = [
            float(s)
            for s in re.findall(
                rb"^ *[0-9]+ +([-+.0-9eE]+) +[-+.0-9eE]+ +[-+.0-9eE]+$",
                table,
                re.MULTILINE,
            )
        ]
        assert len(table_scores) == len(new_scores)
        assert all(abs(a - b) < 0.0051 for a, b in zip(table_scores, new_scores))
        narrow = blocks(dock(*anchor, *rank, window=0)[0])
        assert len(narrow) == 1
        assert atoms(narrow[0]) == atoms(originals[0])
        assert numbers(narrow[0], b"RESULT")[1:] == [0, 0]
        assert numbers(narrow[0], b"AMIDE RESULT")[1:] == [0, 0]
        both_stages = dock(*anchor, *amide, *rank)
        assert both_stages == dock(*anchor, *amide, *rank, cpu=2)
        # A previous full score block is discarded on reading even for default
        # redocking (where score headers from older protocols are preserved).
        old = directory / "prior-ranked.pdbqt"
        old.write_bytes(
            b"REMARK test note\n" + b"\n".join(reordered[0].splitlines()[1:-1]) + b"\n"
        )
        redocked = blocks(dock(ligand=old)[0])
        assert redocked
        assert all(
            b"REMARK VINA AMIDE" not in m
            and m.count(b"REMARK VINA RESULT:") == 1
            and b"REMARK test note" in m
            for m in redocked
        )
        failure = subprocess.run(
            [args.binary, *inputs, "--score_only", *rank],
            capture_output=True,
            check=False,
        )
        assert failure.returncode and b"rank_amide" in failure.stderr
        print(
            f"Final ranking preserved {len(reordered)} original modes, native-window RMSD, "
            "raw components, score table and redocking metadata"
        )
        print(
            f"Native/default, zero barrier, CPU and trajectory parity passed; "
            f"{len(models)} released poses rescored, max discrepancy {worst:.4f} kcal/mol"
        )


if __name__ == "__main__":
    main()
