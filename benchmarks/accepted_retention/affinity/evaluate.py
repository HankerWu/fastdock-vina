"""Recompute scientific endpoints directly from saved first poses and complete pool logs."""

import argparse
import hashlib
import itertools
import json
import re
from decimal import Decimal
from pathlib import Path

import numpy as np
from rdkit import Chem
from run import BENCH, COUNT, RECS, SEED, cells, read, write
from scipy.stats import kendalltau, spearmanr

SCORE = re.compile(r"^REMARK VINA RESULT:\s+(\S+)", re.MULTILINE)


def label_table(work, panel):
    if panel == "Gal3":
        return {
            r["compound"]: (Decimal(r["kd_uM"]), int(r["heavy"]), int(r["charge"]))
            for r in read(BENCH / "pair_gradient/galectin_external/manifest.csv")
        }
    data = json.loads((work / "ca/prepared/metadata.json").read_text())
    result = {}
    for ident, meta in data.items():
        mol = Chem.MolFromSmiles(meta["smiles"])
        result[ident] = (
            Decimal(str(meta["Ki_nM"])),
            mol.GetNumHeavyAtoms(),
            Chem.GetFormalCharge(mol),
        )
    return result


def recompute(work, panel):
    observed = []
    for stage in ("upstream", "fixed", "matched"):
        path = work / f"{panel}_{stage}.csv"
        if not path.exists():
            if stage == "matched":
                continue
            raise FileNotFoundError(path)
        for r in read(path):
            prefix = work / "raw" / r["key"]
            pdbqt = prefix.with_suffix(".pdbqt").read_text()
            values = SCORE.findall(pdbqt)
            assert values[0] == r["first_score"] and len(values) == int(r["modes"])
            contents = pdbqt.split("ENDMDL", 1)[0]
            coords = "\n".join(
                l for l in contents.splitlines() if l.startswith(("ATOM", "HETATM"))
            )
            log = prefix.with_suffix(
                ".counts" if stage == "upstream" else ".log"
            ).read_text()
            for k in COUNT:
                actual = re.findall(r"\b" + k + r"=(\d+)", log)
                assert len(actual) == 1 and actual[0] == r[k], (prefix, k, actual)
            assert int(r["joint"]) == sum(int(r[k]) for k in COUNT[:3])
            deferred = next(
                (line for line in log.splitlines() if "Deferred hunt:" in line), ""
            )
            fields = dict(re.findall(r"(\w+)=(\S+)", deferred))
            observed.append(
                dict(
                    **r,
                    method="matched" if stage == "matched" else r["arm"],
                    explicit_inputs=int(re.search(r"\bexplicit_inputs=(\d+)", log)[1]),
                    pose_sha256=hashlib.sha256(pdbqt.encode()).hexdigest(),
                    first_coordinates_sha256=hashlib.sha256(
                        coords.encode()
                    ).hexdigest(),
                    **{
                        "deferred_" + k: fields.get(k, "")
                        for k in (
                            "skipped",
                            "checks",
                            "records",
                            "eligible",
                            "replaced",
                            "slots",
                            "slot",
                            "pending_energy",
                            "previous_energy",
                        )
                    },
                )
            )
    return observed


def evaluate(work, output):
    output.mkdir(parents=True, exist_ok=True)
    rows, ranks, pairs, uncertainty, costs, changes, inputs = [], [], [], [], [], [], []
    for panel in RECS:
        observed = recompute(work, panel)
        rows.extend(observed)
        labels = label_table(work, panel)
        ids = sorted(labels)
        n = len(ids)
        assert n == {"Gal3": 47, "CAII": 5}[panel]
        recs = RECS[panel]
        methods = sorted({r["method"] for r in observed})
        index = {(r["method"], r["receptor"], r["compound"]): r for r in observed}
        assert len(index) == len(observed) == n * len(recs) * len(methods)
        for method in set(methods) & {"retention", "ungated"}:
            for receiver, compound in itertools.product(recs, ids):
                reference = index["parent", receiver, compound]
                experiment = index[method, receiver, compound]
                for k in (
                    "hunt_evals",
                    "second_evals",
                    "hunt_calls",
                    "second_calls",
                    "explicit_inputs",
                ):
                    assert reference[k] == experiment[k], (
                        panel,
                        method,
                        receiver,
                        compound,
                        k,
                    )
        for _, rec, c, protein, box, lig in cells(work, panel):
            kd, heavy, q = labels[c]
            inputs.append(
                {
                    "panel": panel,
                    "receptor": rec,
                    "compound": c,
                    "k": kd,
                    "unit": "uM" if panel == "Gal3" else "nM",
                    "heavy": heavy,
                    "formal_charge": q,
                    "receptor_sha256": hashlib.sha256(protein.read_bytes()).hexdigest(),
                    "box_sha256": hashlib.sha256(box.read_bytes()).hexdigest(),
                    "ligand_sha256": hashlib.sha256(lig.read_bytes()).hexdigest(),
                }
            )
        draws = np.random.default_rng(SEED).multinomial(
            n, np.repeat(1 / n, n), size=10000
        )
        for rec in (*recs, "mean"):
            group = recs if rec == "mean" else (rec,)
            # Integer sum preserves exact ordering by one score or by the receiver mean.
            score = {}
            for m in methods:
                score[m] = np.array(
                    [
                        sum(
                            int(1000 * Decimal(index[m, s, c]["first_score"]))
                            for s in group
                        )
                        for c in ids
                    ]
                )
            for m in methods:
                li = [index[m, s, c] for s in group for c in ids]
                par = [index["parent", s, c] for s in group for c in ids]
                paired = list(zip(li, par, strict=True))

                totals = {
                    k: sum(int(r[k]) for r in li)
                    for k in ("joint", *COUNT, "explicit_inputs", "modes", "elapsed_ns")
                }

                costs.append(
                    dict(
                        panel=panel,
                        receiver=rec,
                        method=m,
                        cells=len(li),
                        elapsed_s=round(totals["elapsed_ns"] / 1e9, 8),
                        **{
                            k: totals[k]
                            for k in ("joint", *COUNT, "explicit_inputs", "modes")
                        },
                        replacements=sum(int(r["deferred_replaced"] or 0) for r in li),
                        first_score_different=sum(
                            l["first_score"] != p["first_score"] for l, p in paired
                        ),
                        first_coordinates_different=sum(
                            l["first_coordinates_sha256"]
                            != p["first_coordinates_sha256"]
                            for l, p in paired
                        ),
                        output_different=sum(
                            l["pose_sha256"] != p["pose_sha256"] for l, p in paired
                        ),
                    )
                )
            score["size_only"] = np.array([-labels[c][1] for c in ids])
            for fold, subset in (
                (3, "all"),
                (3, "equal_heavy_charge"),
                (10, "all"),
                (10, "equal_heavy_charge"),
            ):
                edges = [
                    (i, j)
                    for i, j in itertools.permutations(range(n), 2)
                    if labels[ids[j]][0] >= fold * labels[ids[i]][0]
                    and (subset == "all" or labels[ids[i]][1:] == labels[ids[j]][1:])
                ]
                strong = np.array([i for i, _ in edges], dtype=int)
                weak = np.array([j for _, j in edges], dtype=int)
                weights = draws[:, strong] * draws[:, weak]
                den = weights.sum(axis=1)
                valid = den > 0
                credits = {
                    m: 1 + np.sign(v[weak] - v[strong]) for m, v in score.items()
                }
                for m, v in score.items():
                    k = [float(labels[c][0]) for c in ids]
                    ranks.append(
                        {
                            "panel": panel,
                            "receiver": rec,
                            "method": m,
                            "fold": fold,
                            "subset": subset,
                            "n_ligands": n,
                            "n_pairs": len(edges),
                            "credit": sum(credits[m]) / 2,
                            "ties": np.count_nonzero(credits[m] == 1),
                            "rho": round(float(spearmanr(k, v).statistic), 6),
                            "tau_b": round(float(kendalltau(k, v).statistic), 6),
                        }
                    )
                    if m == "parent" or not edges:
                        continue
                    delta = credits[m] - credits["parent"]
                    shifts = 50 * (weights @ delta)[valid] / den[valid]
                    q = np.percentile(shifts, [2.5, 97.5])
                    uncertainty.append(
                        {
                            "panel": panel,
                            "receiver": rec,
                            "relative_to_parent": m,
                            "fold": fold,
                            "subset": subset,
                            "n_ligands": n,
                            "in_edges": len(set(strong) | set(weak)),
                            "eligible": len(edges),
                            "valid_draws": valid.sum(),
                            "delta_pp": round(float(50 * np.mean(delta)), 6),
                            "lo_pp": round(float(q[0]), 6),
                            "hi_pp": round(float(q[1]), 6),
                        }
                    )
                if fold == 3 and subset == "all":
                    for h, (i, j) in enumerate(edges):
                        pairs.append(
                            dict(
                                panel=panel,
                                receiver=rec,
                                strong=ids[i],
                                weak=ids[j],
                                k_strong=labels[ids[i]][0],
                                k_weak=labels[ids[j]][0],
                                size_charge_equal=labels[ids[i]][1:]
                                == labels[ids[j]][1:],
                                fold10=labels[ids[j]][0] >= 10 * labels[ids[i]][0],
                                **{
                                    m: (credits[m][h] / 2 if m in credits else "")
                                    for m in (
                                        "size_only",
                                        "upstream",
                                        "parent",
                                        "retention",
                                        "ungated",
                                        "matched",
                                    )
                                },
                            )
                        )
            if rec != "mean":
                for c in ids:
                    p = index["parent", rec, c]
                    for m in methods:
                        if m == "parent":
                            continue
                        r = index[m, rec, c]
                        changes.append(
                            {
                                "panel": panel,
                                "receptor": rec,
                                "compound": c,
                                "method": m,
                                "first_parent": p["first_score"],
                                "first_method": r["first_score"],
                                "coords_differ": r["first_coordinates_sha256"]
                                != p["first_coordinates_sha256"],
                                "emitted_differs": r["pose_sha256"] != p["pose_sha256"],
                                "modes_parent": p["modes"],
                                "modes_method": r["modes"],
                                "replaced": r["deferred_replaced"],
                            }
                        )
    manifest = read(BENCH / "pair_gradient/galectin_external/manifest.csv")
    scaffold = {r["compound"]: r["scaffold"] for r in manifest}
    group_summary = []
    for rec in (*RECS["Gal3"], "mean"):
        eligible = [r for r in pairs if r["panel"] == "Gal3" and r["receiver"] == rec]
        for g in [*sorted(set(scaffold.values())), "between_scaffolds"]:
            if g == "between_scaffolds":
                e = [
                    r for r in eligible if scaffold[r["strong"]] != scaffold[r["weak"]]
                ]
                members = len(scaffold)
            else:
                e = [
                    r
                    for r in eligible
                    if scaffold[r["strong"]] == scaffold[r["weak"]] == g
                ]
                members = sum(a == g for a in scaffold.values())
            group_summary.append(
                {
                    "receiver": rec,
                    "group": g,
                    "members": members,
                    "n_pairs": len(e),
                    **{
                        m: sum(r[m] for r in e)
                        for m in ("size_only", "upstream", "parent", "retention")
                    },
                }
            )
    for name, out in (
        ("galectin_scaffolds", group_summary),
        ("cells", rows),
        ("ranks", ranks),
        ("pairs", pairs),
        ("bootstrap", uncertainty),
        ("cost", costs),
        ("changes", changes),
        ("inputs", inputs),
    ):
        write(output / (name + ".csv"), out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    evaluate(args.work, args.output)
