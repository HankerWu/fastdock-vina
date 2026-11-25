"""Check side-chain transfer geometry and verify results of an extraction run."""

import argparse
import json
from math import isclose
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from extract_independent_receptor_evidence import (  # noqa: E402
    TargetSite, angle_difference, chis, move_sidechain, water_sites,
)


def check_synthetic():
    position_rows = [
        ("N", (0, 0, 1), "N"), ("CA", (0, 0, 0), "C"),
        ("C", (-1, 0, 0), "C"), ("O", (-2, 0, 0), "OA"),
        ("CB", (1.5, 0, 0), "C"), ("CG", (2, 1.2, 0), "C"),
        ("OD1", (3.2, 1.2, 0.2), "OA"), ("ND2", (1.4, 2.1, -0.1), "N"),
        ("HD21", (1.3, 3, -0.05), "HD"),
    ]
    rows = {name: (i, np.array(xyz, dtype=float), "ASN", kind)
            for i, (name, xyz, kind) in enumerate(position_rows)}
    before = {name: row[1] for name, row in rows.items()}
    baseline = chis("ASN", before)
    target = (baseline[0] + 120, baseline[1] - 65)
    changed = move_sidechain(TargetSite((61, ""), "ASN", rows, baseline, 0), target)
    after = {name: changed.get(row[0], row[1]) for name, row in rows.items()}
    assert np.allclose([angle_difference("ASN", i, a, b) for i, (a, b)
                        in enumerate(zip(target, chis("ASN", after)))], [0, 0], atol=1e-10)
    assert all(np.array_equal(before[name], after[name]) for name in ("N", "CA", "C", "O", "CB"))
    assert not np.array_equal(before["HD21"], after["HD21"])
    for a, b in (("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2"), ("ND2", "HD21")):
        assert isclose(np.linalg.norm(before[a] - before[b]), np.linalg.norm(after[a] - after[b]),
                       abs_tol=1e-10)

    hits = [{"pdb_id": pdb_id, "chain": "A", "occupancy": occ,
             "source_resseq": i, "source_icode": "", "xyz": [x, 0, 0]}
            for i, (pdb_id, occ, x) in enumerate([("1AAA", 1.0, 0), ("1AAA", 0.6, 0.1),
                                                  ("1AAB", 1.0, 0.5), ("1AAC", 1.0, 3)])]
    first, second = water_sites(hits)
    assert first["xyz"] == [0.25, 0, 0] and first["independent_pdb_count"] == 2
    assert len(first["same_pdb_repeats"]) == 1 and second["independent_pdb_count"] == 1


def verify_run(manifest, directory):
    targets = json.loads(manifest.read_text())["targets"]
    summary = json.loads((directory / "summary.json").read_text())
    assert len(targets) == len(summary["targets"])
    n_alternates = 0
    for source, reported in zip(targets, summary["targets"]):
        assert source["id"] == reported["id"]
        target = directory / source["id"]
        provenance = json.loads((target / "provenance.json").read_text())
        coordinates = np.array([list(map(float, line.split())) for line in (target / "waters.xyz").read_text().splitlines()])
        sidechains = json.loads((target / "sidechains.json").read_text())
        alternates = sidechains["alternates"]
        assert isclose(sidechains["baseline_population"] + sum(x["population"] for x in alternates), 1.0)
        assert len(coordinates) == reported["water_sites"] and len(alternates) == reported["alternate_receptors"]
        if "donors" in provenance:
            accepted = {donor["pdb_id"] for donor in provenance["donors"] if donor["status"] == "accepted"}
            assert source["query"]["pdb_id"].upper() not in accepted and len(accepted) == reported["accepted_donors"]
            for point, water in zip(coordinates, provenance["water_sites"], strict=True):
                origins = [hit["pdb_id"] for hit in water["hits"]]
                assert len(set(origins)) == len(origins) == water["independent_pdb_count"]
                assert set(origins) <= accepted
                assert np.allclose(point, np.mean([hit["xyz"] for hit in water["hits"]], axis=0),
                                   atol=5.1e-7, rtol=0)
            assert all(set(alt["source_pdbs"]) <= accepted for alt in alternates)
        if alternates:
            original = (manifest.parent / source["query"]["receptor"]).read_text().splitlines()
            for alt in alternates:
                updated = (target / alt["receptor"]).read_text().splitlines()
                for before, after in zip(original, updated, strict=True):
                    if before.startswith(("ATOM  ", "HETATM")):
                        assert before[:30] == after[:30] and before[54:] == after[54:]
                        if before[12:16].strip() in {"N", "CA", "C", "O"}:
                            assert before == after
                    else:
                        assert before == after
                assert original != updated
            n_alternates += len(alternates)
    return n_alternates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, nargs="?")
    parser.add_argument("output", type=Path, nargs="?")
    args = parser.parse_args()
    if (args.manifest is None) != (args.output is None):
        parser.error("provide both manifest and evidence output, or neither")
    check_synthetic()
    if args.manifest:
        count = verify_run(args.manifest, args.output)
        print(f"Verified populated output and same-topology coordinates for {count} alternate receptors.")
    print("Side-chain chi transfer and explicit hydrogen, and cross-PDB water de-duplication passed.")


if __name__ == "__main__":
    main()
