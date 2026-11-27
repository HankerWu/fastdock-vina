"""Compare each Pro pool at identical coordinates, including an independent first-pose affinity score."""

import argparse
import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def main(a):
    sys.path.insert(0, str(a.repo.resolve() / "benchmarks/bond_distance"))
    from fixed import TERMS
    from run import read, write

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        pose = Path(tmp) / "model.pdbqt"
        for r in read(a.out / "runs.csv"):
            if r["arm"] != "parent":
                continue
            pool = (a.out / "raw" / (r["key"] + ".pdbqt")).read_text()
            blocks = re.findall(r"^MODEL .*?^ENDMDL", pool, re.M | re.S)
            if len(blocks) != int(r["modes"]):
                raise ValueError(r["key"])
            for rank, block in enumerate(blocks, 1):
                pose.write_text(block + "\n")
                for reference in ("zero", "self") if rank == 1 else ("zero",):
                    common = [
                        "--receptor",
                        a.work / r["receiver"] / "receptor.pdbqt",
                        "--ligand",
                        a.work / r["query"] / "ligand.pdbqt",
                        "--config",
                        a.work / r["receiver"] / "assay.box",
                        "--pose_coordinates",
                        pose,
                        "--score_only",
                        "--cpu",
                        "1",
                        "--seed",
                        r["seed"],
                    ]
                    if reference == "zero":
                        common += ["--unbound_energy", "0"]
                    pair = {}
                    for arm, binary in (
                        ("parent", a.parent),
                        ("ccd", a.ccd),
                        ("lock", a.lock),
                    ):
                        result = subprocess.check_output([binary, *common], text=True)
                        values = re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", result)
                        if len(values) != len(TERMS):
                            raise ValueError(result)
                        pair[arm] = values
                    if any(terms != pair["parent"] for terms in pair.values()):
                        raise ValueError(
                            f"Fixed score mismatch {r['key']}, rank {rank}"
                        )
                    rows.append(
                        dict(
                            key=r["key"],
                            rank=rank,
                            ref=reference,
                            serialized_pose_sha256=hashlib.sha256(
                                pose.read_bytes()
                            ).hexdigest(),
                            identical=True,
                            **dict(zip(TERMS, pair["parent"])),
                        )
                    )
    write(a.out / "fixed.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("repo", "work", "out", "parent", "ccd", "lock"):
        p.add_argument(k, type=Path)
    main(p.parse_args())
