"""Independent graph, receptor-frame and common-precision energy audit for GBA1."""

import argparse
import hashlib
import json
import re
import subprocess
import sys
from functools import cache
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from diagnose import decomposition

SOURCES = {"2NSX": ("B", "IFM"), "2V3D": ("A", "NBV"), "2V3E": ("A", "NND")}


def coords(row):
    return np.fromstring(row[30:54], sep=" ")


@cache
def source(w, code):
    chain, ligand = SOURCES[code]
    rows = (w / "raw" / (code + ".pdb")).read_text().splitlines()
    ca = {
        r[17:20] + r[22:27]: coords(r)
        for r in rows
        if r.startswith("ATOM  ")
        and r[21] == chain
        and r[12:16].strip() == "CA"
        and r[16] in " A"
    }
    lig = {
        r[12:16].strip(): coords(r)
        for r in rows
        if r.startswith("HETATM")
        and r[17:20] == ligand
        and r[21] == chain
        and r[16] in " A"
        and r[76:78].strip() not in ("H", "D")
    }
    pdbqt = (w / code / "receptor.pdbqt").read_text().splitlines()
    prepared = {
        r[17:20] + r[22:27]: coords(r)
        for r in pdbqt
        if r.startswith("ATOM  ") and r[21] == "A" and r[12:16].strip() == "CA"
    }
    # Prepared coordinates preserve the selected deposited chain (B in 2NSX).
    intersection = sorted(prepared.keys() & ca.keys())
    if (
        len(intersection) < 450
        or max(np.linalg.norm(prepared[k] - ca[k]) for k in intersection) > 0.002
    ):
        raise ValueError(("receptor frame changed", code))
    return ca, lig, np.mean(list(lig.values()), axis=0)


@cache
def fit(w, src, dst):
    a, _, c = source(w, src)
    b, _, d = source(w, dst)
    keys = [
        k
        for k in a.keys() & b.keys()
        if np.linalg.norm(a[k] - c) < 14 and np.linalg.norm(b[k] - d) < 14
    ]
    x = np.asarray([a[k] for k in sorted(keys)])
    y = np.asarray([b[k] for k in sorted(keys)])
    xc, yc = x.mean(0), y.mean(0)
    rot = Rotation.align_vectors(y - yc, x - xc)[0]
    return rot, xc, yc


@cache
def reference(w, code, dst, state):
    _, lig, _ = source(w, code)
    ccd = gemmi.cif.read_file(str(w / "raw" / (SOURCES[code][1] + ".cif"))).sole_block()
    graph, named_elements = {}, {}
    for name, elt in ccd.find("_chem_comp_atom.", ["atom_id", "type_symbol"]):
        if elt not in ("H", "D"):
            graph[name] = set()
            named_elements[name] = elt
    for s, t, order in ccd.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        if s in graph and t in graph:
            assert order == "SING"
            graph[s].add(t)
            graph[t].add(s)
    if set(lig) != set(graph):
        raise ValueError("incomplete deposited ligand")
    tmpl = (w / code / (state + ".pdbqt")).read_text().splitlines()
    (smiles,) = [
        r.split()[2]
        for r in tmpl
        if r.startswith("REMARK SMILES ") and " IDX " not in r
    ]
    mol = Chem.MolFromSmiles(smiles)
    if Chem.GetFormalCharge(mol) != int(state == "protonated"):
        raise ValueError("incorrect graph charge")
    smi = {a.GetIdx(): {b.GetIdx() for b in a.GetNeighbors()} for a in mol.GetAtoms()}
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("unexpected ligand bond order")
    allowed = {
        i: [
            k
            for k in graph
            if named_elements[k] == mol.GetAtomWithIdx(i).GetSymbol()
            and len(graph[k]) == len(smi[i])
        ]
        for i in smi
    }
    order = sorted(smi, key=lambda i: (len(allowed[i]), -len(smi[i]), i))
    mappings = []

    def extend(mapping):
        if len(mapping) == len(order):
            mappings.append(dict(mapping))
            return
        i = order[len(mapping)]
        for candidate in allowed[i]:
            if candidate in mapping.values():
                continue
            if all(
                (j in smi[i]) == (k in graph[candidate]) for j, k in mapping.items()
            ):
                mapping[i] = candidate
                extend(mapping)
                del mapping[i]

    extend({})
    if len(mappings) != 1:
        raise ValueError(
            "GBA ligands expected to have unique element-labelled correspondence"
        )
    (matching,) = mappings
    ref = np.array([lig[matching[a.GetIdx()]] for a in mol.GetAtoms()])
    b = Chem.Mol(mol)
    conf = Chem.Conformer(len(ref))
    for k, x in enumerate(ref):
        conf.SetAtomPosition(k, x)
    b.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(b, replaceExistingTags=True)
    desired = dict(
        Chem.FindMolChiralCenters(
            mol, includeUnassigned=False, useLegacyImplementation=False
        )
    )
    actual = dict(
        Chem.FindMolChiralCenters(
            b, includeUnassigned=False, useLegacyImplementation=False
        )
    )
    for k, cip in desired.items():
        if actual.get(k) != cip:
            raise ValueError(
                ("deposited stereochemistry inconsistent", code, k, cip, actual.get(k))
            )
    mapping = {}
    for s in tmpl:
        if s.startswith("REMARK SMILES IDX"):
            a = list(map(int, s.split()[3:]))
            mapping.update(zip(a[::2], a[1::2]))
    ids = [mapping[i + 1] for i in range(mol.GetNumAtoms())]
    rot, x, y = fit(w, code, dst)
    return ids, rot.apply(ref - x) + y


def unwrapped(body):
    return (
        "\n".join(
            l for l in body.splitlines() if l and not l.startswith(("MODEL", "ENDMDL"))
        )
        + "\n"
    )


def rescore(binary, w, dst, body, tmp):
    path = tmp / "single.pdbqt"
    path.write_text(unwrapped(body))
    args = [
        str(binary),
        "--ligand",
        str(path),
        "--receptor",
        str(w / dst / "receptor.pdbqt"),
        "--cpu",
        "1",
        "--seed",
        "43103",
        "--score_only",
    ]
    for ax, c in zip("xyz", source(w, dst)[2]):
        args.extend(["--center_" + ax, str(c), "--size_" + ax, "24"])
    report = subprocess.check_output(args, text=True)

    def value(key):
        return float(re.search(key + r"\s*:\s*([-0-9.]+)", report)[1])

    return value("Final Intermolecular Energy") + value("Final Total Internal Energy")


def verify(w, binary, out):
    all_deltas, err = [], []
    results = []
    frame_max = 0.0
    archived = pd.read_csv(HERE / "results/crystal_diagnostic.csv")
    with TemporaryDirectory() as scratch:
        tmp = Path(scratch)
        for cap in (12000, 12100, 24000):
            raw = w if cap == 12000 else w / f"budget_{cap}"
            root = HERE / "results" if cap == 12000 else HERE / "results" / f"cap{cap}"
            modes = pd.read_csv(root / "modes.csv")
            if cap != 12000:
                modes = modes[(modes.state == "protonated") & (modes.engine == "pro")]
            paired = pd.read_csv(root / "paired.csv").set_index(
                ["state", "query", "receptor", "seed"]
            )
            for (src, dst, state, engine, seed), group in modes.groupby(
                ["query", "receptor", "state", "engine", "seed"]
            ):
                ids, ref = reference(w, src, dst, state)
                rpath = w / dst / "receptor.pdbqt"
                lpath = w / src / (state + ".pdbqt")
                path = raw / "runs" / f"{src}_{dst}_{state}_{engine}_{seed}.pdbqt"
                blocks = path.read_text().strip().split("ENDMDL")[:-1]
                rms, chp, replay = [], [], []
                for body, row in zip(blocks, group.itertuples(), strict=True):
                    atoms = {
                        int(l[6:11]): coords(l)
                        for l in body.splitlines()
                        if l.startswith(("ATOM  ", "HETATM"))
                    }
                    e = np.sqrt(
                        np.square(np.array([atoms[i] for i in ids]) - ref).sum()
                        / len(ids)
                    )
                    rms.append(e)
                    err.append(abs(e - row.rmsd))
                    if engine != "pro":
                        continue
                    comp = sum(decomposition(rpath, lpath, body, set()))
                    (score,) = [
                        float(l.split()[3])
                        for l in body.splitlines()
                        if l.startswith("REMARK VINA RESULT:")
                    ]
                    val = score + 0.3 * row.rotor_scale * comp
                    if abs(val - row.chpi) > 1e-7:
                        raise ValueError(
                            ("score mismatch", src, dst, seed, row.rank, val, row.chpi)
                        )
                    chp.append(val)
                    if cap == 12000 and state == "protonated":
                        replay.append((rescore(binary, w, dst, body, tmp), comp))
                if engine == "pro":
                    idx = np.argmin(chp)
                    diff = rms[idx] - rms[0]
                    report = paired.loc[state, src, dst, seed]
                    if abs(diff - report.delta) > 1e-7 or idx + 1 != report.chosen_rank:
                        raise ValueError(("re-ranking mismatch", src, dst, state, seed))
                    all_deltas.append(
                        {
                            "cap": cap,
                            "src": src,
                            "dst": dst,
                            "state": state,
                            "seed": seed,
                            "delta": diff,
                            "rank": int(idx) + 1,
                        }
                    )
                if not replay:
                    continue
                for stage in ["deposited", "relaxed"]:
                    probe = (
                        w / "probes" / f"{src}_{dst}_{stage}_single.pdbqt"
                    ).read_text()
                    atoms = {
                        int(l[6:11]): coords(l)
                        for l in probe.splitlines()
                        if l.startswith(("ATOM  ", "HETATM"))
                    }
                    frame = np.sqrt(
                        np.square(np.array([atoms[i] for i in ids]) - ref).sum()
                        / len(ids)
                    )
                    if stage == "deposited":
                        frame_max = max(frame_max, frame)
                    ep = rescore(binary, w, dst, probe, tmp)
                    corr = sum(decomposition(rpath, lpath, probe, set()))
                    a, c = np.array(replay).T
                    native = group.iloc[0].rotor_scale * (ep - a.min())
                    new = group.iloc[0].rotor_scale * (
                        ep + 0.3 * corr - (a + 0.3 * c).min()
                    )
                    old = archived[
                        (archived["query"] == src)
                        & (archived["receptor"] == dst)
                        & (archived.stage == stage)
                        & (archived.seed == seed)
                    ].iloc[0]
                    results.append(
                        {
                            "query": src,
                            "receptor": dst,
                            "seed": seed,
                            "stage": stage,
                            "probe_rmsd": frame,
                            "ref_inter_plus_intra": ep,
                            "rescored_top_native_min": float(a.min()),
                            "ranked_logged_native_winner": int(np.argmin(a)) + 1,
                            "gap_native": native,
                            "gap_chpi": new,
                            "change_from_archived_native": native - old.gap_native,
                            "change_from_archived_chpi": new - old.gap_chpi,
                        }
                    )
    if max(err) > 1e-7 or frame_max > 0.001:
        raise ValueError(("geometry discrepancy", max(err), frame_max))
    df = pd.DataFrame(all_deltas)
    summary = []
    for (cap, state), g in df.groupby(["cap", "state"]):
        summary.append(
            {
                "cap": int(cap),
                "state": state,
                "mean": g.delta.mean(),
                "better": int((g.delta < -1e-10).sum()),
                "worse": int((g.delta > 1e-10).sum()),
                "n": len(g),
            }
        )
    details = pd.DataFrame(results)
    out.mkdir(exist_ok=True, parents=True)
    details.to_csv(out / "independent_placement.csv", index=False, float_format="%.9f")
    (out / "independent.json").write_text(
        json.dumps(
            {
                "modes": len(err),
                "max_RMSD_discrepancy_A": max(err),
                "max_deposited_placement_roundtrip_A": frame_max,
                "common_precision_placement_rows": len(results),
                "worst_native_gap_delta": float(
                    details.change_from_archived_native.abs().max()
                ),
                "worst_chpi_gap_delta": float(
                    details.change_from_archived_chpi.abs().max()
                ),
                "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                "summaries": summary,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "binary", "output"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    verify(a.work.resolve(), a.binary.resolve(), a.output)
