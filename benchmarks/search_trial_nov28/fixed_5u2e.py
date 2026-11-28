"""Replay the preregistered post-unblinding two-pose, two-executable score diagnostic."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "fixed_score_domain"))
from publication import run
from run import score_blocks


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fixtures", "runs", "stable", "candidate", "output"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    fixture = a.fixtures.resolve() / "5U2E"
    raw = a.runs.resolve() / "raw" / "cognate-5U2E-5U2E-85237-candidate.pdbqt"
    work = a.output.resolve()
    work.mkdir(parents=True, exist_ok=False)
    meta = json.loads((fixture / "metadata.json").read_text())
    manifest = json.loads((a.fixtures / "manifest.json").read_text())
    inputs = next(r for r in manifest if r["pdb"] == "5U2E")["sha256"]
    for key in ("receptor.pdbqt", "ligand.pdbqt", "metadata.json"):
        assert sha(fixture / key) == inputs[key]
    box = []
    for axis, val in zip("xyz", meta["center"]):
        box.extend(["--center_" + axis, str(val), "--size_" + axis, "24"])
    blocks = raw.read_text().split("ENDMDL")[:2]
    binaries = {"stable": a.stable.resolve(), "candidate": a.candidate.resolve()}
    rows = []
    for rank, block in enumerate(blocks, 1):
        pose = work / f"rank{rank}.pdbqt"
        lines = [
            l for l in block.splitlines() if l.strip() and not l.startswith("MODEL")
        ]
        pose.write_text("\n".join(lines) + "\n")
        published = float(
            next(l for l in lines if l.startswith("REMARK VINA RESULT:")).split()[3]
        )
        row = dict(rank=rank, published=published, pose_sha256=sha(pose), scores={})
        for arm, binary in binaries.items():
            log = work / f"{arm}_rank{rank}.log"
            cmd = [
                binary,
                "--receptor",
                fixture / "receptor.pdbqt",
                "--ligand",
                fixture / "ligand.pdbqt",
                "--pose_coordinates",
                pose,
                "--score_only",
                "--unbound_energy",
                "0",
                "--spacing",
                "0.375",
                "--cpu",
                "1",
                "--seed",
                "85237",
                *box,
            ]
            receipt = run(cmd, log)
            receipt["terms"] = score_blocks(log.read_text())
            row["scores"][arm] = receipt
        rows.append(row)
    report = dict(
        source_pool=raw.name,
        source_sha256=sha(raw),
        inputs={k: inputs[k] for k in ("ligand.pdbqt", "receptor.pdbqt")},
        executables={arm: sha(b) for arm, b in binaries.items()},
        rows=rows,
    )
    comparable = all(
        len(r["scores"]["stable"]["terms"]) == 1
        and r["scores"]["stable"]["terms"] == r["scores"]["candidate"]["terms"]
        and all(s["exit"] == 0 for s in r["scores"].values())
        and len(r["scores"]["stable"]["terms"][0]) == 10
        for r in rows
    )
    report["comparable"] = comparable
    if comparable:
        far, near = [r["scores"]["candidate"]["terms"][0] for r in rows]
        diff = {k: round(near[k] - far[k], 6) for k in far}
        published_gap = round(rows[1]["published"] - rows[0]["published"], 6)
        fixed_gap = diff["Estimated Free Energy of Binding"]
        report.update(
            near_minus_far=diff,
            published_gap=published_gap,
            gap_discrepancy=round(fixed_gap - published_gap, 6),
            passes_prediction=fixed_gap > 0 and abs(fixed_gap - published_gap) <= 0.03,
        )
    (work / "receipt.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
