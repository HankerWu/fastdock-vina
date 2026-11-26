"""Check the public MMFF affinity annotation against frozen BACE1 results."""

import argparse
import csv
import json
import lzma
import subprocess
import sys
import tempfile
from pathlib import Path

from pose import GEOM, HERE, ROOT, unpack


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    a = p.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        w = Path(temp)
        unpack(GEOM / "affinity/inputs.tar.xz", w)
        unpack(GEOM / "affinity/cap60000/outputs.tar.xz", w)
        comp = "CHEMBL2151135"
        pool = w / "pools" / f"4DJW_{comp}_89609_armijo.pdbqt"
        command = [
            sys.executable,
            str(ROOT / "tools/score_physical_pool.py"),
            "--binary",
            str(a.binary.resolve()),
            "--receptor",
            str(w / "prepared/4DJW.pdbqt"),
            "--prepared-ligand",
            str(w / "prepared" / (comp + ".pdbqt")),
            "--pool",
            str(pool),
        ]

        def run(label, flags):
            out = w / (label + ".json")
            subprocess.run(command + flags + ["--output", str(out)], check=True)
            return json.loads(out.read_text())

        default = run("default", [])
        annotated = run("annotation", ["--mmff-affinity"])
        if "mmff_affinity" in default or default["orders"] != annotated["orders"]:
            raise ValueError("Legacy orders or default output changed")
        annot = annotated["mmff_affinity"]
        with lzma.open(HERE / "affinity_modes.csv.xz", "rt") as stream:
            saved = [
                r
                for r in csv.DictReader(stream)
                if r["cap"] == "60000"
                and r["compound"] == comp
                and r["receptor"] == "4DJW"
                and r["seed"] == "89609"
            ]
        for v, r in zip(annot["per_model"], saved, strict=True):
            if abs(v["affinity_score"] - float(r["full_unit"])) > 1e-8 or v[
                "model"
            ] != int(r["rank"]):
                raise ValueError("MMFF score differs from archived calculation")
        assert (
            annot["native_rank_one"]["model"] == 1
            and annot["minimum_across_emitted_models"]["model"] == 5
        )
        print(
            "BACE1 native rank 1 and pool minimum:",
            annot["native_rank_one"],
            annot["minimum_across_emitted_models"],
        )


if __name__ == "__main__":
    main()
