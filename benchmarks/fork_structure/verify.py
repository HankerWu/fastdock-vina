"""Recompute the structural endpoints, fixed comparisons, CPU allocations and identities."""

import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import evidence
import holdout

D = Path(__file__).resolve().parent
B = D.parent
P = evidence.panel


def identical(file, other):
    if file.read_bytes() != other.read_bytes():
        raise ValueError("reconstruction differs: " + str(other))


def summaries(path, temp):
    temp.mkdir()
    (temp / "cells.csv").symlink_to((path / "cells.csv").resolve())
    evidence.evaluate.summaries(temp)
    subprocess.run(
        [sys.executable, B / "stage_local_structure/directions.py", temp], check=True
    )
    for name in ("summary.csv", "paired.csv", "paired_summary.csv", "by_direction.csv"):
        identical(temp / name, path / name)


def main():
    versions = json.loads((D / "builds.json").read_text())["versions"]
    identities = {k: v["binary_sha256"] for k, v in versions.items()}
    assert (D / "fixed/fixed/binary.sha256").read_text().strip() == identities[
        "candidate"
    ]
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        with tarfile.open(D / "external/inputs.tar.xz") as archive:
            archive.extractall(work / "inputs", filter="data")
        holdout.structures(work / "inputs")
        for part in ("external", "external/repeat"):
            folder = D / part
            assert json.loads((folder / "binaries.json").read_text()) == identities
            check = holdout.external.verify(work / "inputs", folder)
            assert check == json.loads((folder / "verification.json").read_text())
            print(part, check)

        def pools(folder):
            return {
                r["key"]: tuple(
                    r[k]
                    for k in (
                        "pool_sha256",
                        "status",
                        "modes",
                        "hunt_evals",
                        "second_evals",
                        "explicit_calls",
                        "explicit_evals",
                        "late_prescore",
                    )
                )
                for r in evidence.evaluate.read(folder / "cells.csv")
            }

        assert pools(D / "external") == pools(D / "external/repeat")
        identical(D / "external/ranks.csv", D / "external/repeat/ranks.csv")
        summaries(D / "external", work / "ext")
        fixed = D / "fixed"
        assert P.verify(fixed / "fixed", 65) == json.loads(
            (fixed / "fixed/verification.json").read_text()
        )
        P.verify_receipts(fixed, "fixed", {})
        summaries(fixed, work / "summary")
        old = evidence.evaluate.read(B / "stage_local_structure/fixed/cells.csv")
        expected = [
            r
            for r in old
            if r["arm"] != "candidate"
            and r["cohort"] in ("reciprocal", "cognate", "default", "cost")
        ]
        assert evidence.evaluate.read(
            fixed / "cells.csv"
        ) == expected + evidence.evaluate.read(fixed / "fixed/cells.csv")
        evidence.report_actual(work / "summary", "cost")
        for name in ("actual_cost_summary.csv", "actual_cost_pairs.csv"):
            identical(work / "summary" / name, fixed / name)
        for stage in ("probes", "validation", "matched"):
            path = D / "timing" / stage
            if (path / "cells.csv").exists():
                assert json.loads((path / "binaries.json").read_text()) == identities
                cert = evidence.check_cost(path, stage)
                assert cert == json.loads((path / "verification.json").read_text())
                print(stage, cert)
        evidence.verify_allocations(D / "timing")
        evidence.engine.csvfile(
            work / "repeats.csv",
            evidence.cpu_repeats(evidence.evaluate.read(D / "timing/probes/cells.csv")),
        )
        identical(work / "repeats.csv", D / "timing/repeats.csv")
    src = versions["candidate"]["revision"]
    assert (
        subprocess.check_output(
            ["git", "-C", D, "rev-parse", src + ":src"], text=True
        ).strip()
        == versions["candidate"]["src_tree"]
    )
    print("All archived four-source structural results reproduced.")


if __name__ == "__main__":
    main()
