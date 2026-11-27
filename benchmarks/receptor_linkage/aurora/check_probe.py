"""Reconstruct reference frames, stored-pose RMSD and common-zero contrasts."""

import argparse
import csv
import json
import re
import tarfile
from pathlib import Path

import gemmi
from frame_check import ca, read_query, reference_graph
from geometry import RECEPTORS, REFERENCES
from rdkit.Chem import rdMolAlign


def check(run, refs, archive, out):
    raw = out / "raw"
    with (out / "probe_stages.csv").open() as f:
        stages = list(csv.DictReader(f))
    with (out / "probe_pairs.csv").open() as f:
        pairs = list(csv.DictReader(f))
    frames = {}
    with tarfile.open(archive) as f:
        for rid, comp in RECEPTORS.items():
            s = gemmi.read_pdb_string(f.extractfile(f"raw/{rid}.pdb").read())
            s.merge_chain_parts()
            chain = s[0]["A"]
            lig = [
                a.pos
                for r in chain
                if r.name == comp
                for a in r
                if not a.element.is_hydrogen and a.altloc in ("\x00", "A")
            ]
            frames[rid] = {
                key: v
                for key, v in ca(chain).items()
                if min(v.dist(x) for x in lig) <= 12
            }

    def read_score(name):
        content = (raw / (name + "_score.log")).read_text()
        lines = [l for l in content.splitlines() if "Unbound System's Energy" in l]
        assert len(lines) == 1 and float(lines[0].split(":")[1].split()[0]) == 0
        return float(
            re.search(r"Estimated Free Energy of Binding\s*:\s*(\S+)", content)[1]
        )

    rmsd_diffs = []
    energy_diffs = []
    for drug, refs_of_drug in REFERENCES.items():
        probe, serials = read_query(run, drug)
        refid, comp, _ = refs_of_drug[0]
        pdb = refs / (refid + ".pdb")
        atoms = ca(gemmi.read_structure(str(pdb))[0]["A"])
        for rec, frame in frames.items():
            anchors = sorted(atoms.keys() & frame.keys())
            fit = gemmi.superpose_positions(
                [frame[k] for k in anchors], [atoms[k] for k in anchors]
            )
            ref = reference_graph(pdb, comp, "A")
            conf = ref.GetConformer()
            for i, v in enumerate(conf.GetPositions()):
                conf.SetAtomPosition(i, tuple(fit.transform.apply(gemmi.Position(*v))))
            for row in (
                r for r in stages if r["compound"] == drug and r["receptor"] == rec
            ):
                stem = "_".join(
                    row[k]
                    for k in ("receptor", "compound", "graph", "witness", "stage")
                )
                coords = {
                    int(l[6:11]): [float(l[k : k + 8]) for k in (30, 38, 46)]
                    for l in (raw / (stem + ".pdbqt")).read_text().splitlines()
                    if l.startswith(("ATOM  ", "HETATM"))
                }
                for i, serial in enumerate(serials):
                    probe.GetConformer().SetAtomPosition(i, coords[serial])
                rms = rdMolAlign.CalcRMS(
                    probe, ref, symmetrizeConjugatedTerminalGroups=False
                )
                rmsd_diffs.append(abs(rms - float(row["rmsd_A"])))
                score = read_score(stem)
                energy_diffs.append(abs(score - float(row["zero_reference_score"])))
    for row in pairs:
        winner = (
            "_".join(row[k] for k in ("receptor", "compound", "graph", "seed"))
            + "_winner"
        )
        witness = "_".join(
            row[k] for k in ("receptor", "compound", "graph", "witness", "stage")
        )
        gap = read_score(winner) - read_score(witness)
        energy_diffs.append(abs(gap - float(row["winner_minus_witness"])))
    assert (
        len(rmsd_diffs) == len(stages)
        and max(rmsd_diffs) < 1e-9
        and max(energy_diffs) < 1e-12
    )
    result = {
        "stage_rows": len(stages),
        "checked_pairs": len(pairs),
        "max_rmsd_difference_A": max(rmsd_diffs),
        "max_score_difference": max(energy_diffs),
        "all_score_logs_common_unbound": 0,
    }
    (out / "probe_independent.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "references", "archive", "output"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    check(a.run, a.references, a.archive, a.output)
