"""Check every assay projection against archived CLI modes and all trajectory counters."""

import argparse
import itertools
import json
import re
from pathlib import Path

from assay import MEASURED, OLD, read, run_panel


def main(original, new, results):
    cells = read(results / "all_cells.csv")
    assert len(cells) == 60 and len({r["key"] for r in cells}) == 60
    rows = read(results / "fixed_top1.csv")
    assert len(rows) == 60 and len({(r["key"], r["arm"]) for r in rows}) == 60
    for cell in cells:
        d = original if cell["query"] in OLD else new
        raw = d / "runs" / (cell["key"] + ".pdbqt")
        assert run_panel.sha(raw) == cell["output_sha256"]
        top = float(
            re.search(r"^REMARK VINA RESULT:\s*(\S+)", raw.read_text(), re.M).group(1)
        )
        assert abs(top - float(cell["top1_score"])) < 1e-9
        log = (d / "runs" / (cell["key"] + ".counts")).read_text()
        for k in (
            "hunt_calls",
            "hunt_evals",
            "second_calls",
            "second_evals",
            "explicit_calls",
            "explicit_evals",
            "explicit_inputs",
        ):
            vals = re.findall(r"\b" + k + r"=(\d+)", log)
            assert len(vals) == (1 if k.startswith("explicit") else 4) and sum(
                map(int, vals)
            ) == int(cell[k])
    for key in {r["key"] for r in rows}:
        fields = ("score", "inter", "intra", "torsional", "unbound", "pose_sha256")
        assert len({tuple(r[k] for k in fields) for r in rows if r["key"] == key}) == 1
    w = read(results / "ligands.csv")
    lg = {i["id"]: i for i in w if i["arm"] == "parent" and i["receiver"] == "both"}
    ok = tie = large_ok = 0
    for a, b in itertools.combinations(
        sorted(lg, key=lambda i: float(lg[i]["Kd_uM"])), 2
    ):
        na, nb = int(lg[a]["N_heavy"]), int(lg[b]["N_heavy"])
        ok += na > nb
        tie += na == nb
        if float(lg[b]["Kd_uM"]) / float(lg[a]["Kd_uM"]) >= 3:
            large_ok += na > nb
    print(
        json.dumps(
            dict(
                cells=len(cells),
                raw_poses_and_counters=True,
                fixed_scores=60,
                fixed_coordinates=20,
                all_fixed_terms_identical=True,
                size_baseline_correct=ok,
                size_baseline_tied=tie,
                size_baseline_large_correct=large_ok,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for a in ("original", "new", "results"):
        p.add_argument(a, type=Path)
    v = p.parse_args()
    main(v.original, v.new, v.results)
