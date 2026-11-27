"""Independent coordinate, scalar-score and mixture checks for the two-receptor study."""
import argparse
import csv
import io
import itertools
import json
import re
import subprocess
import sys
import tarfile
from decimal import Decimal, localcontext
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import rdMolAlign

AURORA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(AURORA))
from frame_check import ca, positions, read_query, reference_graph

STATES = {"2J4Z": "626", "2X6E": "YM4"}
REF = {"at9283": ("2W1G", "L0G"), "danusertib": ("2J50", "627")}
LABELS = ("witness", "60913", "60943", "61211", "61231")
RT = Decimal("0.00198720425864") * Decimal("298.15")


def read(path):
    with path.open() as inp:
        return list(csv.DictReader(inp))


def cli(args):
    return subprocess.run(list(map(str, args)), text=True, capture_output=True, check=True).stdout


def extract_inputs(work):
    for name in ("docking.tar.xz", "probe_raw.tar.xz", "solver_replay/outputs.tar.xz", "references.tar.xz"):
        with tarfile.open(AURORA / name) as src:
            src.extractall(work, filter="data")
    frames = {}
    with tarfile.open(AURORA.parents[1] / "torsion_kinase_affinity/evidence.tar.xz") as archive:
        for rec, ligand in STATES.items():
            model = gemmi.read_pdb_string(archive.extractfile(f"raw/{rec}.pdb").read())
            model.merge_chain_parts()
            chain = model[0]["A"]
            atoms = ca(chain)
            box = [a.pos for r in chain if r.name == ligand for a in r
                   if not a.element.is_hydrogen and a.altloc in ("\x00", "A")]
            pocket = {k: x for k, x in atoms.items() if min(x.dist(v) for v in box) <= 12}
            frames[rec] = {"whole": atoms, "pocket": pocket}
    return frames


def fit(source, target):
    keys = sorted(source.keys() & target.keys())
    return gemmi.superpose_positions([target[k] for k in keys], [source[k] for k in keys])


def first_body(path):
    first = path.read_text().split("ENDMDL")[0].strip().splitlines()
    return [l for l in first if not l.startswith("MODEL")]


def pool_file(bodies, transform, path):
    parts = []
    for i, body in enumerate(bodies, 1):
        lines = []
        for line in body:
            if transform is not None and line.startswith(("ATOM  ", "HETATM")):
                pos = gemmi.Position(*(float(line[j:j+8]) for j in (30, 38, 46)))
                new = transform.apply(pos)
                line = line[:30] + "".join(f"{v:8.3f}" for v in new) + line[54:]
            lines.append(line)
        parts += [f"MODEL {i}", *lines, "ENDMDL"]
    path.write_text("\n".join(parts) + "\n")


def geometry(work, drug, source_ca, pool):
    pdb, comp = REF[drug]
    path = work / f"{pdb}.pdb"
    structure = gemmi.read_structure(str(path))
    structure.merge_chain_parts()
    align = fit(ca(structure[0]["A"]), source_ca).transform
    ref = reference_graph(path, comp, "A")
    conf = ref.GetConformer()
    for j, coord in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(j, tuple(align.apply(gemmi.Position(*coord))))
    query, ids = read_query(work, drug)
    result, posed_queries = [], []
    for pose in positions(pool):
        for i, serial in enumerate(ids):
            query.GetConformer().SetAtomPosition(i, pose[serial])
        result.append(rdMolAlign.CalcRMS(query, ref, symmetrizeConjugatedTerminalGroups=False))
        posed_queries.append(Chem.Mol(query))
    separations = [rdMolAlign.CalcRMS(a, b, symmetrizeConjugatedTerminalGroups=False)
                   for a, b in itertools.combinations(posed_queries[1:], 2)]
    return result, {"minimum_A": min(separations), "maximum_A": max(separations)}


def reconstruct(exe, work, records, centers):
    frames = extract_inputs(work)
    nodes = sorted({int(r["nodes_per_scale"]) for r in records})
    expected = {(r["origin"], r["ligand"], r["state"], r["registration"],
                 int(r["nodes_per_scale"]), r["pose"]): r for r in records}
    references = {(r["origin"], r["ligand"], r["pose"]): float(r["rmsd_A"]) for r in centers}
    checked = set()
    largest = {"rmsd_A": 0., "center": 0., "window": 0., "printed_score": 0.}
    fits, decoy_spacing = {}, {}
    for src, drug in itertools.product(STATES, REF):
        ligand = work / "prepared" / f"{drug}.pdbqt"
        graph = "chemical" if drug == "at9283" else "native"
        chem = work / "raw" / f"{drug}_{graph}.bonds"
        raw = work / "raw" / f"{src}_{drug}_{'sdf' if drug == 'at9283' else 'native'}_legal_placed.pdbqt"
        near = work / "localized.pdbqt"
        center = json.loads((work / "prepared" / f"{src}.json").read_text())["center"]
        box = [v for axis, val in zip("xyz", center, strict=True)
               for v in (f"--center_{axis}", val, f"--size_{axis}", 26)]
        cli([exe, "--receptor", work / "prepared" / f"{src}.pdbqt", "--ligand", raw,
             "--cpu", 1, "--seed", 60913, "--ligand_bonds", chem, *box, "--local_only", "--out", near])
        bodies = [first_body(near)] + [first_body(work / "pools" / f"pro_{src}_{drug}_{seed}.pdbqt")
                                       for seed in LABELS[1:]]
        for dst in STATES:
            for frame in (("pocket",) if dst == src else ("pocket", "whole")):
                transform = fit(frames[src][frame], frames[dst][frame])
                key = f"{src}_to_{dst}_{frame}"
                fits[key] = {"atoms": transform.count, "rmsd_A": transform.rmsd}
                pool = work / "transported.pdbqt"
                pool_file(bodies, None if src == dst else transform.transform, pool)
                if src == dst:
                    rmsds, decoy_spacing[f"{src}_{drug}"] = geometry(work, drug, frames[src]["pocket"], pool)
                    for label, val in zip(LABELS, rmsds, strict=True):
                        largest["rmsd_A"] = max(largest["rmsd_A"], abs(val - references[src, drug, label]))
                args = [exe, "--receptor", work / "prepared" / f"{dst}.pdbqt", "--cpu", 1,
                        "--ligand", ligand, "--ligand_bonds", chem, "--score_only"]
                scores = []
                for block in pool.read_text().split("ENDMDL")[:-1]:
                    single = work / "single.pdbqt"
                    single.write_text(block.strip() + "\nENDMDL\n")
                    msg = cli([*args, "--pose_coordinates", single, "--unbound_energy", 0, "--autobox"])
                    val = re.search(r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)", msg)
                    scores.append(float(val[1]))
                for n in nodes:
                    msg = cli([*args, "--rigid_window_samples", n, "--rigid_window_pool", pool, "--verbosity", 0])
                    rows = csv.DictReader(io.StringIO(msg[msg.index("model,inter,intra,"):]))
                    for label, score, row in zip(LABELS, scores, rows, strict=True):
                        key = (src, drug, dst, frame, n, label)
                        ref = expected[key]
                        e = (float(row["inter"]) + float(row["intra"])) * float(row["slope"])
                        f = e + float(row["window_delta"])
                        for field, diff in (("center", e - float(ref["center"])),
                                            ("window", f - float(ref["window"])), ("printed_score", e - score)):
                            largest[field] = max(largest[field], abs(diff))
                        checked.add(key)
    if checked != expected.keys() or max(largest[k] for k in ("rmsd_A", "center", "window")) > 1e-8:
        raise ValueError(("Coordinate reconstruction disagrees", len(checked), largest))
    if largest["printed_score"] > 0.000501:
        raise ValueError(("Independent score_only reading disagrees", largest))
    return {"rows_rebuilt": len(checked), "geometry_centers": len(references), "max_errors": largest,
            "registrations": fits, "within_origin_blind_rmsd": decoy_spacing}


def mixture(rows, bounds, records):
    """Directly exponentiate decimal energies. Solve affine inequalities in p, including endpoints."""
    certify = []
    raw = {(r["registration"],r["origin"],r["ligand"],r["pose"],r["state"],r["nodes_per_scale"]):r
           for r in records}
    with localcontext() as ctx:
        ctx.prec = 90
        for r in bounds:
            keys = ("registration", "nodes", "model")
            subset = [s for s in rows if all(r[k] == s[k] for k in keys)]
            pairs, intersections = [], set()
            max_margin_error = max_odds_error = 0.
            for s in subset:
                for name, pose in (("near","witness"), ("blind",s["seed"])):
                    for state in STATES:
                        reg = "pocket" if s["origin"] == state else s["registration"]
                        row = raw[reg,s["origin"],s["ligand"],pose,state,s["nodes"]]
                        if Decimal(s[f"{name}_{state}"]) != Decimal(row[s["model"]]):
                            raise ValueError(("Source energies differ", s))
                values = [[(-Decimal(s[f"{which}_{state}"]) / RT).exp() for state in STATES]
                          for which in ("near", "blind")]
                diff = [v - b for v, b in zip(*values, strict=True)]
                pairs.append(diff)
                for state, d in zip(STATES, diff, strict=True):
                    if (int(d > 0)-int(d < 0)) != int(s[f"sign_{state}"]):
                        raise ValueError(("Wrong sign convention", s))
                cross = s["crossover_log_odds_2X6E"]
                if diff[0] * diff[1] < 0:
                    p = diff[0] / (diff[0] - diff[1])
                    intersections.add(p)
                    odds = (p / (1-p)).ln()
                    max_odds_error = max(max_odds_error, abs(float(odds)-float(cross)))
                elif cross:
                    raise ValueError(("Unexpected crossover", s))
                margin = RT * (sum(values[0]) / sum(values[1])).ln()
                max_margin_error = max(max_margin_error, abs(float(margin)-float(s["margin_half"])))
            breaks = sorted({Decimal(0), Decimal(1), *intersections})
            points = {(a+b)/2 for a, b in itertools.pairwise(breaks)} | set(breaks)
            def count(p):
                return sum(a * (1-p) + b * p > 0 for a, b in pairs)
            checks = {"wins_2J4Z": count(0), "wins_2X6E": count(1), "wins_half": count(Decimal("0.5")),
                      "max_shared": max(map(count, points)), "max_per_pair": sum(max(x) > 0 for x in pairs)}
            if max(max_margin_error,max_odds_error) > 1e-8 or not all(int(r[k]) == v for k,v in checks.items()):
                raise ValueError((r, checks, max_margin_error, max_odds_error))
            certify.append(dict(**{k:r[k] for k in keys}, max_margin_error=max_margin_error,
                                max_crossover_log_odds_error=max_odds_error))
    return certify


def sensitivity(records, pairs):
    levels = sorted({int(r["nodes_per_scale"]) for r in records})
    by_pose = {(r["registration"], r["origin"], r["ligand"], r["pose"], r["state"], int(r["nodes_per_scale"])): r
               for r in records}
    extremes = []
    for lo, hi in itertools.pairwise(levels):
        differences = [(abs(float(row["window"]) - float(by_pose[(*key[:-1],lo)]["window"])), key[:-1])
                       for key,row in by_pose.items() if key[-1] == hi]
        d,key = max(differences)
        extremes.append({"coarse": lo, "fine": hi, "max_window_difference": d, "case": list(key)})
    w = [r for r in pairs if r["model"] == "window"]
    grouped = {}
    for r in w:
        for s in STATES:
            k = (r["registration"], r["origin"], r["ligand"], r["seed"], s)
            grouped.setdefault(k,{})[int(r["nodes"])] = float(r[f"blind_{s}"]) - float(r[f"near_{s}"])
    switches = [{"case": list(k), "blind_minus_witness": v} for k,v in grouped.items()
                if len({x > 0 for x in v.values()}) > 1]
    fin = [r for r in records if int(r["nodes_per_scale"]) == levels[-1]]
    ess = min(fin,key=lambda r:float(r["effective_samples"]))
    return {"levels": levels, "max_energy_changes": extremes, "sign_changes_any_level": switches,
            "sign_changes_8192_to_finest": sum((v[8192]>0)!=(v[levels[-1]]>0) for v in grouped.values())
            if 8192 in levels else None,
            "min_effective_sample_size_at_finest": float(ess["effective_samples"]),
            "min_ess_case": [ess[k] for k in ("origin","ligand","pose","state","registration")],
            "reliability": "Deterministic Sobol level agreement checks sign sensitivity only. "
                           "Low effective sample sizes and near-zero margins leave some crossings unresolved."}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("binary","work","results"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    records, centers, pairs, bounds = [read(a.results / f"{name}.csv")
                                      for name in ("state_scores", "centers", "pairs", "bounds")]
    report = reconstruct(a.binary.resolve(), a.work.resolve(), records, centers)
    report["decimal_mixture_checks"] = mixture(pairs, bounds, records)
    report["quadrature"] = sensitivity(records, pairs)
    (a.results / "verification.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
