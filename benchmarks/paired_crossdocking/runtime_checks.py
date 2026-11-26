"""Controls on one real crystal, including the exact zero-repulsion limiting case."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(root, binary, baseline, work):
    work.mkdir(parents=True, exist_ok=True)
    # This uses the already completed 1A52/3ERD screen cell for a compiler/output parity check.
    meta = json.loads((root / "metadata.json").read_text())
    s, r = "1A52", "3ERD"
    arg = [
        "--ligand",
        str(root / s / "ligand.pdbqt"),
        "--receptor",
        str(root / r / "receptor.pdbqt"),
        "--pose_selection",
        "late",
        "--exhaustiveness",
        "4",
        "--max_evals",
        "3000",
        "--num_modes",
        "10",
        "--refine_candidates",
        "40",
        "--energy_range",
        "100",
        "--cpu",
        "1",
        "--search_counters",
        "--seed",
        "43103",
    ]
    for axis, pos in zip("xyz", meta[r]["center"]):
        arg += ["--center_" + axis, str(pos), "--size_" + axis, "24"]
    specs = {
        "baseline": (baseline, []),
        "patch_native": (binary, []),
        "patch_soft": (binary, ["--soft_receptor_search"]),
        "zero_native": (binary, ["--weight_repulsion", "0"]),
        "zero_soft": (binary, ["--weight_repulsion", "0", "--soft_receptor_search"]),
    }
    out = {}
    for label, (exe, flags) in specs.items():
        dest = work / (label + ".pdbqt")
        subprocess.run(
            [str(exe), *arg, *flags, "--out", str(dest)],
            capture_output=True,
            text=True,
            check=True,
        )
        out[label] = sha(dest)
    assert (
        out["baseline"]
        == out["patch_native"]
        == sha(root / "runs_3000" / "1A52_3ERD_43103_native.pdbqt")
    )
    assert out["zero_native"] == out["zero_soft"]
    assert out["patch_soft"] == sha(root / "runs_3000" / "1A52_3ERD_43103_soft.pdbqt")
    assert out["patch_soft"] != out["patch_native"]
    threaded = arg.copy()
    threaded[threaded.index("--cpu") + 1] = "2"
    dest = work / "threaded.pdbqt"
    subprocess.run(
        [str(binary), *threaded, "--soft_receptor_search", "--out", str(dest)],
        capture_output=True,
        check=True,
    )
    out["two_threads"] = sha(dest)
    assert out["two_threads"] == out["patch_soft"]
    for incompatible in (
        ["--mixed_cutoff_search"],
        ["--smooth_cutoff"],
        ["--score_only"],
    ):
        result = subprocess.run(
            [
                str(binary),
                *arg,
                "--soft_receptor_search",
                *incompatible,
                "--out",
                str(work / "invalid.pdbqt"),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert (
            result.returncode
            and "soft_receptor_search" in result.stderr + result.stdout
        )
    (work / "binary_checks.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for k in ("root", "binary", "baseline", "work"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    main(*(x.resolve() for x in (a.root, a.binary, a.baseline, a.work)))
