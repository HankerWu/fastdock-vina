"""Replicate the fixed BPTF water protocol on human BRDT BD1, Asn109 pocket."""
import argparse
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bptf_water"))
from water import (align_donor, aligned_waters, align_with_numbering, fixed_water, load_meta,
                   measure, structure_chain, water_sites, write_csv)
from registration import match_observed

RECEIVERS = {"7MRH": "ZMJ", "7MRG": "ZMM"}
DONORS = ("4FLP", "9YCQ")
ARMS = ("dry", "polar", "supported", "registered")
SEEDS = (95701, 95713, 95717)


def sites_and_fits(input_dir, pdb, donors):
    entry = input_dir / pdb
    _, protein = structure_chain(entry / (pdb + ".pdb"), pdb, "A")
    key = (109, "")
    if protein[key].name != "ASN":
        raise ValueError("Recognition Asn109 is absent")
    anchor = protein[key].atoms["ND2"]
    waters, fits = [], {}
    for donor in DONORS:
        path = donors / (donor + ".pdb")
        if not path.exists():
            path.write_bytes(urlopen(f"https://files.rcsb.org/download/{donor}.pdb", timeout=45).read())
        model, chain = structure_chain(path, donor, "A")
        report = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        mapping, rotation, shift = align_donor(protein, chain, anchor, report)
        if mapping[key] != key or chain[key].name != "ASN":
            raise ValueError("Asn109 does not map between the two BRDT structures")
        waters.extend(w for w in aligned_waters(model, donor, anchor, np.full(3, 16), rotation, shift)
                      if w["chain"] == "A" and w["occupancy"] >= 0.8)
        fits[donor] = report
    dry = entry / "receptor.pdbqt"
    protein_xyz = np.array([[float(l[i:i+8]) for i in (30,38,46)] for l in dry.read_text().splitlines()
                            if l.startswith(("ATOM  ", "HETATM")) and l[77:].strip() not in ("H","HD")])
    sites = [s for s in water_sites(waters) if s["independent_pdb_count"] == 2
             and np.linalg.norm(s["xyz"]-anchor) <= 8
             and np.linalg.norm(protein_xyz-np.asarray(s["xyz"]), axis=1).min() >= 2.5]
    matched = match_observed(sites, entry / (pdb+".pdb"), dry)
    return dict(anchor="ASN109:ND2", fits=fits, sites=sites, matched=matched)


def build_receptors(input_dir, work):
    evidence = {}
    donors = work / "donors"
    donors.mkdir(exist_ok=True)
    for pdb in RECEIVERS:
        evidence[pdb] = e = sites_and_fits(input_dir, pdb, donors)
        for arm in ARMS:
            mapped = e["sites"] if arm == "polar" else e["matched"] if arm != "dry" else []
            key = "observed_xyz" if arm == "registered" else "xyz"
            positions = [np.asarray(site[key]) for site in mapped]
            if any(np.linalg.norm(a-b) < 2.4 for i,a in enumerate(positions) for b in positions[:i]):
                raise ValueError("Adjacent water oxygens are closer than 2.4 A")
            directory = work / arm / pdb
            directory.mkdir(parents=True, exist_ok=True)
            (directory/"ligand.pdbqt").write_bytes((input_dir/pdb/"ligand.pdbqt").read_bytes())
            text = (input_dir/pdb/"receptor.pdbqt").read_text().rstrip() + "\n"
            for i,xyz in enumerate(positions):
                text += "\n".join(fixed_water({"xyz":xyz}, 9001+3*i, 800+i, True)) + "\n"
            (directory/"receptor.pdbqt").write_text(text)
    (work/"evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path, help="Prepared inputs and donor PDB files in inputs/ and donors/")
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    work = a.work.resolve()
    inputs = work / "inputs"
    build_receptors(inputs, work)
    meta = {pdb:load_meta(inputs, pdb) for pdb in RECEIVERS}
    jobs = []
    for query,qm in meta.items():
        for receiver,rm in meta.items():
            rot,tr,rms,n,off = np.eye(3),np.zeros(3),0.,0,0
            if query != receiver:
                rot,tr,rms,n,off = align_with_numbering(qm,rm)
            mapped = np.asarray(qm["mapped"])@rot+tr
            for i,seed in enumerate(SEEDS):
                for arm in ARMS[i:]+ARMS[:i]:
                    out = work/arm/"poses"
                    out.mkdir(exist_ok=True)
                    jobs.append((work/arm,out,a.binary.resolve(),qm,rm,mapped,rms,n,off,seed,arm,120000))
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        results = list(pool.map(measure, jobs))
    write_csv(work/"cells.csv", [r for r,m in results])
    write_csv(work/"modes.csv", [p for r,m in results for p in m])


if __name__ == "__main__":
    main()
