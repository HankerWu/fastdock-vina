"""Replace the four archived modes of the one known malformed FXa ligand."""

import argparse
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

from solvent_strain import ARMS, FROZEN, decompose, pair_rows, score, write


def main(out, platform):
    old = FROZEN / "strain_transfer_audit"
    ident = "CHEMBL369480"
    assay = "CHEMBL661804"
    meta = {
        r["compound"]: r
        for r in decompose.read(FROZEN / "strain_components/scores.csv")
    }
    rows = [r for r in decompose.read(out / "modes.csv.xz") if r["assay"] == assay]
    ledger = decompose.read(old / "redock_cells.csv")
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        with tarfile.open(old / "redock_raw.tar.xz") as archive:
            archive.extractall(work, filter="data")
        cases = []
        for r in (r for r in ledger if r["arm"] == "random_coords"):
            p = work / f"random_coords_{r['receptor']}_{r['seed']}.pdbqt"
            cases.append((decompose.read_first(p), r))
        assert len(cases) == 4
        replacements, _ = score(
            {"compound": ident, "assay": assay}, cases[0][0], cases, platform
        )
    write(out / "corrected_modes.csv", replacements)
    rows = [r for r in rows if r["compound"] != ident] + replacements
    cells = {}
    for r in rows:
        cells.setdefault(r["compound"], {})[r["pdb"], r["seed"]] = r
    assert len(cells) == 28 and all(len(v) == 4 for v in cells.values())
    expected = {
        (r["receptor"], r["subset"]): r
        for r in decompose.read(old / "redock_summary.csv")
        if r["arm"] == "random_coords"
    }
    results = []
    for receptor in ("1FJS", "2P16", "both"):
        states = {
            c: {k: r for k, r in v.items() if receptor == "both" or k[0] == receptor}
            for c, v in cells.items()
        }
        edges = list(pair_rows(assay, receptor, states, meta, ARMS))
        for subset in ("all", "incident", "remaining"):
            group = [
                e
                for e in edges
                if subset == "all"
                or ((ident in (e["stronger"], e["weaker"])) == (subset == "incident"))
            ]
            summary = {arm: sum(e[arm] for e in group) for arm in ARMS}
            if any(
                summary[a] != float(expected[receptor, subset][a])
                for a in ("native", "full")
            ):
                raise ValueError("Historical input replacement diverged")
            results.append(
                dict(receptor=receptor, subset=subset, pairs=len(group), **summary)
            )
    write(out / "corrected_summary.csv", results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--platform", choices=("Reference", "CPU"), default="Reference")
    a = p.parse_args()
    main(a.out, a.platform)
