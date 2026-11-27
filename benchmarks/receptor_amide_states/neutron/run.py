"""Compare a protein-only amide prior with deposited neutron N-D assignments."""

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from coupled import CoupledAmides
from model import AmideReceptor, RT, TERMINI, atoms


def save(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def xyz(line):
    return np.array([float(line[i : i + 8]) for i in (30, 38, 46)])


def reference_status(site, raw):
    oname, nname, _ = TERMINI[site[:3]]
    terminals = {n: [l for l in raw if l[12:16].strip() == n] for n in [oname, nname]}
    if any(
        len(ls) != 1 or ls[0][16] != " " or float(ls[0][54:60]) < 1.0
        for ls in terminals.values()
    ):
        return "ambiguous_terminal", []
    n, o = xyz(terminals[nname][0]), xyz(terminals[oname][0])
    # Each D must be at covalent N distance and closer to N than O.
    candidates = sorted(
        (float(l[54:60]), l)
        for l in raw
        if l[76:78].strip() == "D"
        and 0.7 <= np.linalg.norm(xyz(l) - n) <= 1.25
        and np.linalg.norm(xyz(l) - n) < np.linalg.norm(xyz(l) - o)
    )
    distinct = []
    for occ, line in reversed(candidates):
        if all(np.linalg.norm(xyz(line) - xyz(prior)) > 0.65 for _, prior in distinct):
            distinct.append((occ, line))
    if len(distinct) != 2:
        return "missing_or_ambiguous_ND", distinct
    return (
        "eligible" if min(v[0] for v in distinct) >= 0.5 else "low_D_occupancy"
    ), distinct


def audit(entry, work, mk):
    raw = gzip.decompress((HERE / (entry + ".pdb.gz")).read_bytes())
    lines = raw.decode().splitlines()
    if not any(l.startswith("EXPDTA") and "NEUTRON" in l for l in lines):
        raise ValueError(f"{entry}: no neutron experiment declared")
    residues = defaultdict(list)
    for l in lines:
        if l.startswith("ATOM  ") and l[21] == "A":
            residues[l[17:27]].append(l)
    # Preserve deposited protein heavy coordinates. Rebuild all H with templates.
    selected = [
        l[:16] + " " + l[17:]
        for group in residues.values()
        for l in group
        if l[76:78].strip() not in {"H", "D"} and l[16] in {" ", "A"}
    ]
    pdb, prepared = work / (entry + ".pdb"), work / (entry + ".pdbqt")
    pdb.write_text("\n".join(selected) + "\nEND\n")
    subprocess.run(
        [mk, "--read_pdb", str(pdb), "-p", str(prepared)],
        check=True,
        capture_output=True,
    )
    text = prepared.read_text()
    parsed = {
        (a.residue, a.name): a
        for a in atoms(text).values()
        if a.kind not in {"H", "HD"}
    }
    geom = {key: a.xyz for key, a in parsed.items()}
    source = {(l[17:27], l[12:16].strip()): xyz(l) for l in selected}
    # Meeko can interchange the two equivalent terminal carboxylate O atom names.
    # Permit only an exact coordinate swap of two OAs; no amide permutation is allowed.
    terminal_swaps = 0
    for residue in residues:
        o, xt = (residue, "O"), (residue, "OXT")
        if {o, xt} <= source.keys() & geom.keys() and parsed[o].kind == parsed[
            xt
        ].kind == "OA":
            if np.array_equal(geom[o], source[xt]) and np.array_equal(
                geom[xt], source[o]
            ):
                geom[o], geom[xt] = geom[xt], geom[o]
                terminal_swaps += 1
    missing = source.keys() - geom.keys()
    max_error = max(
        np.linalg.norm(source[k] - geom[k]) for k in source.keys() & geom.keys()
    )
    if missing or max_error > 1e-6:
        raise ValueError(
            f"{entry}: missing observed heavy atoms {missing}, displacement {max_error}"
        )
    receptor = AmideReceptor(text)
    graph = CoupledAmides(receptor)
    joint = {}
    for component in graph.components:
        for i, site_index in enumerate(component.nodes):
            clamped = component.field.copy()
            clamped[i] = np.inf
            joint[receptor.sites[site_index].residue] = np.exp(
                component.log_partition(clamped) - component.apo
            )
    scores = {s.residue: s.apo_delta for s in receptor.sites}
    rows = []
    for site, group in residues.items():
        if site[:3] not in TERMINI:
            continue
        status, ds = reference_status(site, group)
        if site not in scores:
            status = "unscorable:" + status
        e = scores.get(site, np.nan)
        rows.append(
            dict(
                entry=entry,
                residue=site.strip(),
                status=status,
                nD=len(ds),
                D_occupancies=" ".join(f"{d[0]:.2f}" for d in ds),
                delta_apo=e,
                p_deposited=expit(e / RT),
                p_deposited_joint=joint.get(site, np.nan),
            )
        )
    valid = [r for r in rows if r["status"] == "eligible"]
    summary = dict(
        entry=entry,
        total=len(rows),
        eligible=len(valid),
        raw_sha256=hashlib.sha256(raw).hexdigest(),
        original_heavy=len(source),
        added_heavy=len(geom.keys() - source.keys()),
        max_observed_shift_A=max_error,
        equivalent_terminal_name_swaps=terminal_swaps,
    )
    for key in ["p_deposited", "p_deposited_joint"]:
        ps = np.array([r[key] for r in valid])
        summary[key + "_preferred"] = int(np.sum(ps > 0.5))
        summary[key + "_ties"] = int(np.sum(ps == 0.5))
        summary[key + "_mean"] = np.mean(ps) if len(ps) else np.nan
    return rows, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--meeko", default="mk_prepare_receptor.py")
    args = parser.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    rows, summaries = [], []
    with TemporaryDirectory() as tmp:
        for code in ["3TMJ", "4G0C", "2R24", "5CE4"]:
            rs, s = audit(code, Path(tmp), args.meeko)
            rows.extend(rs)
            summaries.append(s)
    save(args.output / "sites.csv", rows)
    save(args.output / "summary.csv", summaries)
