"""Audit fixed HUP start scores and decompose matched receptor-state changes."""
import argparse
import hashlib
import json
import re
import resource
import subprocess
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

TERMS = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]
RADII = np.array([1.9, 1.9, 1.8, 1.8, 1.8, 1.8, 1.7, 1.7, 1.7, 1.7,
                  2., 2.1, 1.5, 1.8, 2., 2.2, 2.2, 2.3, 1.2])
ORACLE = {"pro24k": 1, "pro48k": 2, "official24k": 4}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clock():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return np.array([time.perf_counter(), time.process_time() + child.ru_utime + child.ru_stime])


def pdb_atoms(path):
    return [(line, tuple(float(line[c:c+8]) for c in (30, 38, 46)))
            for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]


def resname(line):
    return f"{line[21]}:{line[17:20].strip()}{line[22:27].strip()}"


def deposited_ids(state, receptor_path, meta):
    archive = Path(__file__).resolve().parents[2]/"receptor_amide_states/transfer/transfer_evidence.tar.xz"
    if digest(archive) != meta["sources_archive_sha256"]:
        raise ValueError("Source archive digest mismatch")
    with tarfile.open(archive) as tar:
        data = tar.extractfile(f"novel/raw/{state}.pdb").read()
    if hashlib.sha256(data).hexdigest() != meta["hashes"][state+".pdb"]:
        raise ValueError("Deposited PDB digest mismatch")
    ca = {}
    for line in data.decode().splitlines():
        if line.startswith("ATOM  ") and line[21] == "A" and line[16] in " A" and line[12:16].strip() == "CA":
            key = line[17:20], tuple(float(line[c:c+8]) for c in (30, 38, 46))
            if key in ca:
                raise ValueError("Duplicate C-alpha")
            ca[key] = resname(line)
    result = {}
    for line, position in pdb_atoms(receptor_path):
        if line[12:16].strip() == "CA":
            result[resname(line)] = ca[line[17:20], position]
    if len(result) != len(set(result.values())):
        raise ValueError("Many-to-one residue mapping")
    return result


def atom_motion(raw, frame, mappings):
    shared = []
    for state in ("4EY6", "4EY7"):
        records = {}
        for line, xyz in pdb_atoms(raw/state/"receptor.pdbqt"):
            if line.split()[-1] in ("H", "HS", "HD"):
                continue
            residue = mappings[state][resname(line)]
            key = residue, line[12:16].strip()
            if key in records:
                raise ValueError(f"Duplicate identity {state} {key}")
            records[key] = np.asarray(xyz)
        shared.append(records)
    first, second = shared
    rotation = np.asarray(frame["from_4ey6_rotation"])
    shift = np.asarray(frame["from_4ey6_translation"])
    for residue, atom in sorted(first.keys() & second.keys()):
        delta = first[residue, atom] @ rotation + shift - second[residue, atom]
        yield dict(residue=residue, receptor_atom=atom, displacement_A=np.linalg.norm(delta))


def probe(binary, kind, path):
    output = subprocess.check_output([binary, kind, path], text=True)
    return [line.split() for line in output.splitlines()]


def receptor(binary, path, names):
    records = {pos: line for line, pos in pdb_atoms(path) if line.split()[-1] not in ("H", "HD", "HS")}
    positions, labels, serials, types = [], [], [], []
    for tag, atom_type, x, y, z in probe(binary, "receptor", path):
        point = tuple(map(float, (x, y, z)))
        line = records.pop(point)
        labels.append(names[resname(line)])
        positions.append(point)
        serials.append(int(line[6:11]))
        types.append(int(atom_type))
        if tag != "R":
            raise ValueError(line)
    if records:
        raise ValueError("Unmatched receptor atoms")
    return np.array(positions), np.array(labels), np.array(serials), np.array(types)


def pair_terms(r2, a, b):
    """Both native smooth samples, then their mean; radii in Å."""
    bins = np.floor(32. * r2)
    radius = RADII[a] + RADII[b]
    hydrophobic = [0, 12, 13, 14, 15]
    donor, acceptor = [3, 5, 7, 9, 18], [4, 5, 8, 9]
    hp = np.isin(a, hydrophobic) & np.isin(b, hydrophobic)
    hb = ((np.isin(a, donor) & np.isin(b, acceptor))
          | (np.isin(b, donor) & np.isin(a, acceptor)))
    result = np.zeros((len(r2), len(TERMS)))
    for node in (bins, bins + 1):
        d = np.sqrt(node / 32.) - radius
        values = np.column_stack((-.035579*np.exp(-4*d*d),
                                  -.005156*np.exp(-((d-3.)/2.)**2),
                                  .840245*np.minimum(d, 0.)**2,
                                  -.035069*hp*np.clip(1.5-d, 0., 1.),
                                  -.587439*hb*np.clip(-d/.7, 0., 1.)))
        result += .5 * (node < 2048)[:, None] * values
    return result


def score(binary, path, rec, tree, center):
    topology = probe(binary, "ligand", path)
    active_torsions = int(topology[0][1])
    atoms = {int(line[6:11]): pos for line, pos in pdb_atoms(path)}
    ligand_types = {int(row[1]): int(row[2]) for row in topology if row[0] == "L"}
    serials = np.array(list(ligand_types))
    points = np.array([atoms[s] for s in serials])
    box_margin = 18. - np.abs(points - center).max()
    if box_margin < .001:
        raise ValueError("Coordinates beyond or too close to printed 36 Å box")
    receiver, residues, receiver_serial, receiver_types = rec
    triples = [(i, j) for i, js in enumerate(tree.query_ball_point(points, 8.)) for j in js]
    i, j = np.array(triples, dtype=int).T
    delta = points[i] - receiver[j]
    r2 = (delta * delta).sum(axis=1)
    i, j, r2 = i[r2 < 64.], j[r2 < 64.], r2[r2 < 64.]
    vals = pair_terms(r2, np.array(list(ligand_types.values()))[i], receiver_types[j])
    atom_energy = np.bincount(i, weights=vals.sum(axis=1), minlength=len(points))
    curled = atom_energy / (1. + np.maximum(atom_energy, 0.) / 1000.)
    cap_remainder = (curled - atom_energy).sum()
    internal = 0.
    pair_count = 0
    for row in topology:
        if row[0] == "P":
            s, t = map(int, row[1:])
            d = np.subtract(atoms[s], atoms[t])
            if d @ d < 64.:
                v = pair_terms(np.array([d @ d]), ligand_types[s], ligand_types[t]).sum()
                internal += v / (1. + max(v, 0.) / 1000.)
                pair_count += 1
    divisor = 1. + .05846 * active_torsions
    energy = (curled.sum() + internal) / divisor
    contacts = pd.DataFrame(vals / divisor, columns=TERMS)
    contacts.insert(0, "r_A", np.sqrt(r2))
    contacts.insert(0, "receptor_serial", receiver_serial[j])
    contacts.insert(0, "ligand_serial", serials[i])
    contacts.insert(0, "residue", residues[j])
    stats = dict(score=energy, raw_pair_sum=vals.sum()/divisor,
                 cap_remainder=cap_remainder/divisor, internal=internal/divisor,
                 active_torsions=active_torsions, internal_pairs=pair_count,
                 inter_pairs=len(r2), box_margin_A=box_margin)
    return stats, contacts


def csv_gz(frame, path):
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def main(raw, binary, out, inputs):
    started = clock()
    out.mkdir(parents=True, exist_ok=True)
    receipts, failures, pair_tables, contrasts, totals = [], [], [], [], []
    meta = json.loads(inputs.read_text())
    if meta["oracle_mode"] != ORACLE:
        raise ValueError("Original Oracle identities changed")
    checksums = {"probe": digest(binary), "inputs": digest(inputs)}
    mappings = {}
    for state in ("4EY5", "4EY6", "4EY7"):
        folder = raw/state
        recpath = folder/"receptor.pdbqt"
        if state not in meta["included"]:
            failures.append(dict(state=state, id="", reason=meta["cognate_guard"][state]["excluded"]))
            continue
        if not recpath.exists():
            failures.append(dict(state=state, id="", reason="No prepared receptor in registered archive"))
            continue
        if digest(recpath) != meta["hashes"][f"raw/{state}/receptor.pdbqt"]:
            raise ValueError(f"Substituted receptor in {state}")
        mappings[state] = deposited_ids(state, recpath, meta)
        rec = receptor(binary, recpath, mappings[state])
        tree = cKDTree(rec[0])
        checksums[f"{state}/receptor.pdbqt"] = digest(recpath)
        summaries = {}
        coordinate_ids = sorted({"projection", *(arm+"_01" for arm in ORACLE),
                                 *(f"{arm}_{rank:02d}" for arm, rank in ORACLE.items())})
        for case in coordinate_ids:
            job = folder/case
            pose = job/"input.pdbqt"
            log = job/"initial.log"
            if not pose.exists() or not log.exists():
                failures.append(dict(state=state, id=case, reason="Missing start or native score"))
                continue
            content = log.read_text()
            center = [float(v) for v in re.search(r"Grid center:\s*X ([^ ]+) Y ([^ ]+) Z ([^\n ]+)", content).groups()]
            reference = float(re.search(r"Estimated Free Energy of Binding\s*:\s*(\S+)", content)[1])
            vals, pairs = score(binary, pose, rec, tree, center)
            dev = vals["score"] - reference
            ok = abs(dev) <= .001
            receipts.append(dict(state=state, id=case, input_sha256=digest(pose),
                                 native=reference, discrepancy=dev, arithmetic_ok=ok, **vals))
            if not ok:
                failures.append(dict(state=state, id=case, reason=f"Score disagreement {dev:.7f}"))
                continue
            sums = pairs.groupby("residue")[TERMS].sum()
            sums["sum"] = sums.sum(axis=1)
            summaries[case] = vals, sums, reference
            pair_tables.append(pairs.assign(state=state, id=case))
        for arm, rank in ORACLE.items():
            baseline = arm + "_01"
            for test, case in (("oracle", f"{arm}_{rank:02d}"), ("projection", "projection")):
                if baseline not in summaries or case not in summaries:
                    continue
                x, p, ref1 = summaries[case]
                y, q, ref0 = summaries[baseline]
                frame = p.subtract(q, fill_value=0.)
                contrasts.append(frame.reset_index().assign(state=state, pool=arm, contrast=test))
                chem = {k: x[k]-y[k] for k in ("score", "cap_remainder", "internal")}
                totals.append(dict(state=state, pool=arm, contrast=test, test=case,
                                   baseline=baseline, native_gap=ref1-ref0, **chem,
                                   **frame[TERMS].sum().to_dict()))
                if abs(chem["score"]-chem["cap_remainder"]-chem["internal"]-frame["sum"].sum()) > 1e-8:
                    raise ValueError("Residue contrast does not close")
    pd.DataFrame(receipts).to_csv(out/"scores.csv", index=False)
    pd.DataFrame(failures, columns=["state","id","reason"]).to_csv(out/"failures.csv", index=False)
    pd.DataFrame(totals).to_csv(out/"contrasts.csv", index=False)
    csv_gz(pd.concat(pair_tables, ignore_index=True), out/"pairs.csv.gz")
    residue = pd.concat(contrasts, ignore_index=True)
    csv_gz(residue, out/"residues.csv.gz")
    map_rows = [dict(state=state, prepared=prepared, deposited=deposited)
                for state, lookup in mappings.items() for prepared, deposited in lookup.items()]
    csv_gz(pd.DataFrame(map_rows), out/"numbering.csv.gz")
    changes = []
    if {"4EY6", "4EY7"} <= set(meta["included"]):
        motion = pd.DataFrame(atom_motion(raw, meta["targets"]["4EY7"], mappings))
        motion = motion[motion.residue.isin(residue.residue)]
        csv_gz(motion, out/"receptor_motion.csv.gz")
    for (pool, contrast), sub in residue.groupby(["pool","contrast"]):
        frames = {state: rows.set_index("residue")[TERMS+["sum"]] for state, rows in sub.groupby("state")}
        for state in ["4EY5", "4EY7"]:
            if state in frames and "4EY6" in frames:
                change = frames[state].subtract(frames["4EY6"], fill_value=0)
                changes.append(change.reset_index().assign(pool=pool, contrast=contrast,
                                                           receptor_delta=state+"-4EY6"))
    csv_gz(pd.concat(changes, ignore_index=True), out/"state_changes.csv.gz")
    cost = clock() - started
    manifest = dict(hashes=checksums, start_evaluations=len(receipts), gradients=0,
                    receptor_ligand_pairs=sum(row["inter_pairs"] for row in receipts),
                    requested_contrasts=18, verified_contrasts=len(totals),
                    cpu_s=cost[1], wall_s=cost[0])
    (out/"receipt.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("raw", "binary", "out", "inputs"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    main(args.raw, args.binary.resolve(), args.out, args.inputs)
