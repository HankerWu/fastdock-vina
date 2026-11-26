"""Check Kd labels, raw rank-one energies, and reported pair counts independently."""
import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def audit(work):
    identities = {r["compound"]: r for r in json.loads((ROOT / "identities.json").read_text())["ligands"]}
    # Ycas Table 1: unlabelled BPTF SPR 0.200 and 9.0 micromolar.
    kd = {"TP-238": 200, "Bromosporine": 9000}
    chembl = {"CHEMBL4872300": "BZ1", "CHEMBL4861491": "21", "CHEMBL4848701": "22", "CHEMBL4863631": "24"}
    for row in json.loads((ROOT / "binding_records.json").read_text()):
        assert (row["assay_chembl_id"], row["standard_type"], row["standard_units"], row["standard_relation"]) == (
            "CHEMBL4828237", "Kd", "nM", "=")
        kd[chembl[row["molecule_chembl_id"]]] = float(row["standard_value"])
    assert kd == {k: v["kd_nM"] for k, v in identities.items()}
    observed = {}
    for p in (work / "pools").glob("*.pdbqt"):
        name, state, cap, seed, arm = p.stem.split("_")
        assert p.read_bytes().endswith(b"ENDMDL\n")
        score = float(p.read_text().splitlines()[1].split()[3])
        observed.setdefault((name, state, cap, arm), []).append((int(seed), score))
    count = sum(map(len, observed.values()))
    assert count == len(read(work / "cells.csv"))
    for cell in read(work / "summary.csv"):
        arm_names = [x for x in cell if x not in {"assay", "receptor", "cap", "scope", "pairs"}]
        for arm in arm_names:
            pred = {}
            for n, entry in identities.items():
                if entry["assay"] != cell["assay"] or (cell["scope"] == "without-24" and n == "24"):
                    continue
                run = observed[n, cell["receptor"], cell["cap"], arm]
                assert len({seed for seed, _ in run}) == 3
                pred[n] = sorted(score for _, score in run)[1]
            good, total = 0.0, 0
            for x, y in combinations(sorted(pred, key=kd.get), 2):
                if kd[y] < 3 * kd[x]:
                    continue
                total += 1
                good += (pred[x] < pred[y]) + 0.5 * (pred[x] == pred[y])
            assert (good, total) == (float(cell[arm]), int(cell["pairs"]))
    print(f"Verified {count} docking outputs, all Kd records and {len(read(work / 'summary.csv'))} summary rows.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    audit(parser.parse_args().work)
