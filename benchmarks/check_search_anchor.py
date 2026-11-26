"""Check search-only guidance, determinism, unmodified scores and invalid inputs."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

BASIC = Path(__file__).resolve().parents[1] / "examples" / "basic"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    inputs = [
        "--receptor",
        BASIC / "1iep_receptor.pdbqt",
        "--config",
        BASIC / "1iep_receptor.box.txt",
        "--ligand",
        BASIC / "1iep_ligand.pdbqt",
    ]
    opts = [
        "--seed",
        "314159",
        "--exhaustiveness",
        "3",
        "--max_evals",
        "2500",
        "--num_modes",
        "4",
        "--energy_range",
        "100",
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        anchor = tmp / "sites.txt"
        out = tmp / "run.pdbqt"

        def run(binary, *flags):
            result = subprocess.run(
                [binary, *inputs, *opts, *flags, "--out", out],
                capture_output=True,
                check=True,
            )
            return out.read_bytes(), re.findall(
                rb"^Search trajectory: .*$", result.stdout, re.MULTILINE
            )

        for select in ("native", "late"):
            flags = ["--pose_selection", select, "--cpu", "1"]
            assert run(args.binary, *flags) == run(args.baseline, *flags)
        anchor.write_text("3 16 52 15\n39 18 53 18 shell 2\n")
        one, _ = run(args.binary, "--search_anchor", anchor, "--cpu", "1")
        assert one == run(args.binary, "--search_anchor", anchor, "--cpu", "2")[0]
        assert one != run(args.binary, "--pose_selection", "late", "--cpu", "1")[0]
        original = anchor.read_text()
        anchor.write_text(
            original + "OR\n" + "\n".join(reversed(original.splitlines())) + "\n"
        )
        assert one == run(args.binary, "--search_anchor", anchor, "--cpu", "1")[0]
        anchor.write_text(original + "OR\n3 18 53 18\n39 16 52 15 shell 2\n")
        ensemble, guided = run(
            args.binary, "--search_anchor", anchor, "--cpu", "1", "--search_counters"
        )
        assert ensemble == run(args.binary, "--search_anchor", anchor, "--cpu", "2")[0]
        _, native = run(
            args.binary, "--pose_selection", "late", "--cpu", "1", "--search_counters"
        )
        flags = [
            "--search_anchor",
            anchor,
            "--search_anchor_trajectories",
            "even",
            "--search_counters",
        ]
        mixed, traces = run(args.binary, *flags, "--cpu", "1")
        assert (mixed, traces) == run(args.binary, *flags, "--cpu", "2")
        assert len(guided) == len(native) == 3 and guided[1] != native[1]
        assert traces == [native[0], guided[1], native[2]]
        odd, reverse = run(
            args.binary,
            "--search_anchor",
            anchor,
            "--search_anchor_trajectories",
            "odd",
            "--search_counters",
            "--cpu",
            "1",
        )
        assert reverse == [guided[0], native[1], guided[2]]
        assert (
            run(
                args.binary,
                "--search_anchor",
                anchor,
                "--search_anchor_trajectories",
                "all",
                "--cpu",
                "1",
            )[0]
            == ensemble
        )
        for flags in (
            ["--search_anchor", anchor, "--search_anchor_trajectories", "other"],
            ["--search_anchor_trajectories", "even"],
            [
                "--search_anchor",
                anchor,
                "--search_anchor_trajectories",
                "even",
                "--exhaustiveness",
                "1",
            ],
        ):
            failed = subprocess.run(
                [args.binary, *inputs, *flags, "--out", out],
                capture_output=True,
                text=True,
            )
            assert failed.returncode and "search_anchor_trajectories" in failed.stderr
        for model in re.findall(
            rb"^MODEL .*?^ENDMDL$",
            one + ensemble + mixed + odd,
            re.MULTILINE | re.DOTALL,
        ):
            fixed = tmp / "fixed.pdbqt"
            fixed.write_bytes(model + b"\n")
            score = subprocess.check_output(
                [
                    args.binary,
                    *inputs,
                    "--pose_coordinates",
                    fixed,
                    "--score_only",
                    "--unbound_energy",
                    "0",
                    "--cpu",
                    "1",
                ],
                text=True,
            )
            inter = float(re.search(r"Ligand - Receptor\s+:\s+(\S+)", score)[1])
            saved = float(re.search(rb"^REMARK INTER:\s+(\S+)", model, re.MULTILINE)[1])
            assert abs(inter - saved) < 0.06
        for bad in (
            "",
            "OR\n3 16 52 15",
            "3 16 52 15\nOR",
            "3 16 52 15\nOR\nOR\n3 16 52 15",
            "3 16 52 15\nOR extra\n39 16 52 15",
            "3 0 0 0",
            "99999 16 52 15",
            "4 16 52 15",
            "3 16 52 15\n3 17 53 16",
            "3 16 52 15 shell -1",
            "3 16 52 15 unknown 1",
            "3 16 52 15 shell 100",
            "3 16 52 15 shell 2 extra",
            "3 nan 52 15",
            "3 16 52",
        ):
            anchor.write_text(bad)
            failed = subprocess.run(
                [
                    args.binary,
                    *inputs,
                    *opts,
                    "--search_anchor",
                    anchor,
                    "--cpu",
                    "1",
                    "--out",
                    out,
                ],
                capture_output=True,
                check=False,
                text=True,
            )
            assert failed.returncode and "search_anchor" in failed.stderr, (
                bad,
                failed.stderr,
            )
        print(
            "Native/Late compatibility, deterministic workers, released native scores and invalid inputs passed"
        )


if __name__ == "__main__":
    main()
