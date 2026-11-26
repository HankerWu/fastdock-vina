"""Filter declared-graph donor assignments by the actual geometry of each input ammonium state."""

import argparse
import itertools
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "tools"), str(ROOT / "benchmarks/cutoff_switch/hcaii")]
from ammonium_states import (
    LIGANDS,
    RECEPTORS,
    SEEDS,
    dock,
    evaluate,
    evaluate_all,
    mirror_state,
    signed_volume,
    xyzs,
)
from coordination_probe import write_csv
from make_scaffold_anchors import query_graph
from rank_anchor_pool import read_sites
from rank_state_pools import rank


def signature(text):
    graph, serial = query_graph(text)
    (center,) = [
        a
        for a in graph.GetAtoms()
        if a.GetAtomicNum() == 7
        and a.GetFormalCharge() == 1
        and a.GetDegree() == 3
        and a.GetTotalNumHs() == 1
    ]
    if center.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
        raise ValueError("Expected input with unspecified amine configuration")
    return (
        serial[center.GetIdx()],
        *sorted(serial[a.GetIdx()] for a in center.GetNeighbors()),
    )


def prepare(work, anchors, out, relax):
    out.mkdir(parents=True, exist_ok=False)
    inputs = out / "inputs"
    inputs.mkdir()
    graphs = json.loads((anchors / "graphs.json").read_text())
    rows = []
    for lig in LIGANDS:
        donor = graphs[lig]["donor"]
        if donor == lig:
            raise ValueError("Source must be a different compound")
        original = (work / "prepared" / f"{lig}_ligand.pdbqt").read_text()
        sig = signature(original)
        mirror, _ = mirror_state(work / "prepared", lig)
        (inputs / f"{lig}_reflected.pdbqt").write_text(mirror)
        for rec, state in itertools.product(RECEPTORS, ("original", "reflected")):
            text = (anchors / f"{rec}_{lig}.multi").read_text()
            blocks = text.rstrip("\n").split("\nOR\n")
            if len(blocks) != len(graphs[lig]["maps"]):
                raise ValueError("Map count mismatch")
            inputvol = signed_volume(
                xyzs(original if state == "original" else mirror), sig
            )
            if abs(inputvol) < 1e-4:
                raise ValueError("Degenerate input geometry")
            accepted = []
            decisions = []
            for block in blocks:
                sites = read_sites(block)
                if not set(sig) <= sites.keys():
                    decisions.append("partial")
                    accepted.append(block)
                    continue
                pts = {s: sites[s][0] for s in sig}
                if any(sites[s][1] != 0 for s in sig):
                    raise ValueError("Expected point anchors")
                v = signed_volume(pts, sig)
                if abs(v) < 1e-4:
                    raise ValueError("Degenerate donor geometry")
                decisions.append("compatible" if inputvol * v > 0 else "opposite")
                if decisions[-1] == "compatible":
                    accepted.append(block)
            fallback = not accepted
            if fallback:
                if relax:
                    # Predefined second experiment: release the center and first carbon shell.
                    accepted = [
                        "\n".join(
                            line
                            for line in b.splitlines()
                            if int(line.split()[0]) not in sig
                        )
                        for b in blocks
                    ]
                    accepted = list(dict.fromkeys(accepted))
                    if any(not block for block in accepted):
                        raise ValueError("No distal landmarks remain")
                else:
                    accepted = blocks
            (inputs / f"{rec}_{lig}_{state}.anchor").write_text(
                "\nOR\n".join(accepted) + "\n"
            )
            rows.append(
                {
                    "receptor": rec,
                    "ligand": lig,
                    "donor": donor,
                    "state": state,
                    "center": sig[0],
                    "neighbors": " ".join(map(str, sig[1:])),
                    "input_volume": inputvol,
                    "n_input": len(blocks),
                    "compatible": decisions.count("compatible"),
                    "opposite": decisions.count("opposite"),
                    "partial": decisions.count("partial"),
                    "fallback": int(fallback),
                    "n_output": len(accepted),
                    "decisions": " ".join(decisions),
                }
            )
    write_csv(out / "filter.csv", rows)
    return rows


def run(work, anchors, out, binary, jobs, relax):
    prepare(work, anchors, out, relax)
    items = list(
        itertools.product(
            RECEPTORS, LIGANDS, SEEDS, ("common", "second_original", "second_reflected")
        )
    )

    def run_dock(item):
        rec, lig, _, slot = item
        state = "reflected" if slot == "second_reflected" else "original"
        dock(
            item,
            binary,
            work,
            out,
            anchor_file=out / "inputs" / f"{rec}_{lig}_{state}.anchor",
        )

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(run_dock, items))
    metadata = json.loads((work / "prepared/metadata.json").read_text())
    with ThreadPoolExecutor(jobs) as pool:
        rows = [
            r
            for part in pool.map(
                lambda k: evaluate(k, binary, work, out, metadata), items
            )
            for r in part
        ]
    cells, pairs = evaluate_all(rows, metadata)
    for name, rs in (("candidates", rows), ("cells", cells), ("pairs", pairs)):
        write_csv(out / f"{name}.csv", rs)
    rank(out, metadata, state_anchors=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("anchors", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--release-conflicting-neighborhood", action="store_true")
    a = p.parse_args()
    run(
        a.work.resolve(),
        a.anchors.resolve(),
        a.out.resolve(),
        a.binary.resolve(),
        a.jobs,
        a.release_conflicting_neighborhood,
    )
