"""Compare serialized release modes on identical fixed-budget queries by rank."""

import argparse
import math
from pathlib import Path
from run import cells
from evaluate import write


def compare(work, out):
    summary = []
    for baseline in ("official", "parent", "stable"):
        rows = []
        for state, metadata, seed in cells("nominal"):
            name = f"{state}_{metadata['compound']}_{seed}.pdbqt"
            base, final = [
                work / "outputs" / engine / "nominal" / name for engine in (baseline, "candidate")
            ]
            if not base.is_file() or not final.is_file():
                rows.append(
                    dict(
                        receptor=state,
                        compound=metadata["compound"],
                        seed=seed,
                        coverage="missing",
                        baseline=baseline,
                        rank=None,
                        same_block=None,
                        score_delta=None,
                        max_atom_shift=None,
                    )
                )
                continue
            old, new = [p.read_bytes().split(b"ENDMDL\n") for p in (base, final)]
            assert old.pop() == new.pop() == b""
            for i in range(max(len(old), len(new))):
                both = i < min(len(old), len(new))
                row = dict(
                    receptor=state,
                    compound=metadata["compound"],
                    seed=seed,
                    baseline=baseline,
                    rank=i + 1,
                    coverage="both" if both else "baseline_only" if i < len(old) else "candidate_only",
                    same_block=None,
                    score_delta=None,
                    max_atom_shift=None,
                )
                if both:
                    lines = [mode.splitlines() for mode in (old[i], new[i])]
                    scores = [
                        float(next(l for l in block if l.startswith(b"REMARK VINA RESULT:")).split()[3])
                        for block in lines
                    ]
                    coords = [
                        [
                            [float(l[j : j + 8]) for j in (30, 38, 46)]
                            for l in block
                            if l.startswith((b"ATOM", b"HETATM"))
                        ]
                        for block in lines
                    ]
                    assert len(coords[0]) == len(coords[1])
                    row.update(
                        same_block=old[i] == new[i],
                        score_delta=round(scores[1] - scores[0], 4),
                        max_atom_shift=max(math.dist(a, b) for a, b in zip(*coords)),
                    )
                rows.append(row)
        write(out / f"pool_changes_{baseline}.csv", rows)
        summary.append(
            dict(
                baseline=baseline,
                shared_modes=sum(r["coverage"] == "both" for r in rows),
                identical_modes=sum(r["same_block"] is True for r in rows),
                identical_top=sum(r["rank"] == 1 and r["same_block"] is True for r in rows),
                changed_printed_score=sum(
                    r["score_delta"] is not None and r["score_delta"] != 0 for r in rows
                ),
            )
        )
    write(out / "pool_changes.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    compare(a.work, a.out)
