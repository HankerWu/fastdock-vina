"""Construct the two single-stage control source trees from the fixed code and its parent."""

import argparse
import subprocess
from pathlib import Path

BASE = "21714700baecee21e617f0c472b1136e4be33b4a"


def patch(root, arm):
    old = subprocess.check_output(
        ["git", "show", BASE + ":src/lib/coords.cpp"], text=True
    )
    body = old[
        old.index("void add_to_output_container(") : old.index(
            "output_container select_distinct_poses("
        )
    ]
    files = ("parallel_mc.cpp", "vina.cpp") if arm == "save" else ("monte_carlo.cpp",)
    for name in files:
        path = root / "src/lib" / name
        text = path.read_text()
        if text.count("add_to_output_container(") != 1:
            raise ValueError(f"Unexpected call sites in {path}")
        text = text.replace("add_to_output_container(", "baseline_insert(")
        body_local = body.replace(
            "void add_to_output_container(", "static void baseline_insert("
        )
        path.write_text('#include "coords.h"\n' + body_local + text)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("checkout", type=Path)
    p.add_argument("arm", choices=["save", "merge"])
    args = p.parse_args()
    patch(args.checkout, args.arm)
