"""Run MAPK14 using the historical independent dock/evaluate paths and frozen identities."""

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stage_local_structure"))
import external


def structures(fixtures):
    manifest = json.loads((fixtures / "manifest.json").read_text())
    expected = [("1BL7", "SB4"), ("1DI9", "MSQ")]
    assert [(r["pdb"], r["ligand"]) for r in manifest] == expected
    cases, exclusions = [], []
    for row in manifest:
        folder = fixtures / row["pdb"]
        for name, digest in row["sha256"].items():
            assert external.run_panel.digest(folder / name) == digest
        if row["status"] != "ready":
            exclusions.append(dict(query=row["pdb"], reason=row["reason"]))
            continue
        meta = json.loads((folder / "metadata.json").read_text())
        lines = (folder / (row["pdb"] + ".pdb")).read_text().splitlines()
        assert meta["target"] == "MAPK14" and meta["chain"] == "A"
        assert any(
            l.startswith("DBREF ")
            and l[12] == "A"
            and l[33:39] == "Q16539"
            and int(l[14:18]) == int(l[55:60])
            and int(l[20:24]) == int(l[62:67])
            for l in lines
        )
        ca = {}
        for line in lines:
            if (
                line.startswith("ATOM  ")
                and line[21] == "A"
                and line[12:16].strip() == "CA"
                and line[16] in " A"
            ):
                key = line[22:27].strip() + ":" + line[17:20]
                assert key not in ca
                ca[key] = [float(line[k : k + 8]) for k in (30, 38, 46)]
        cases.append(dict(meta, ca=ca))
    return cases, exclusions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("run", "collect", "check"):
        cli = sub.add_parser(action)
        cli.add_argument("fixtures", type=Path)
        if action != "check":
            cli.add_argument("work", type=Path)
        if action == "run":
            cli.add_argument("arm", choices=external.run_panel.ARMS)
            cli.add_argument("binary", type=Path)
        else:
            cli.add_argument("output", type=Path)
    args = vars(parser.parse_args())
    action = args.pop("action")
    external.structures = structures
    if action == "check":
        print(json.dumps(external.verify(**args), indent=2))
    else:
        {"run": external.run, "collect": external.collect}[action](**args)
