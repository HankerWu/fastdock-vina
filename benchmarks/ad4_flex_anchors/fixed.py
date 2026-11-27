"""Counterfactual rescoring of the exact same rounded, released parent poses.

AD4 parameters below cover the observed ligand types; the two axis atoms are
carbon. Energies are kcal/mol, positions Angstrom, charges in elementary units.
"""
import argparse
import csv
import itertools
import math
from pathlib import Path
import re
import subprocess
import tempfile

from panel import atoms, coordinates, write_csv

# AD4_parameters.dat: half equilibrium diameter, epsilon, solpar, volume.
PARAMETERS = {
    "C": (2, .15, -.00143, 33.5103), "A": (2, .15, -.00052, 33.5103),
    "N": (1.75, .16, -.00162, 22.4493), "NA": (1.75, .16, -.00162, 22.4493),
    "OA": (1.60, .20, -.00251, 17.1573), "HD": (1, .02, .00051, 0),
    "F": (1.545, .08, -.00110, 15.448)
}


def correction(ligand, side):
    """Independent AD4 midpoint table followed by per-pair positive-energy cap."""
    anchor_atoms = [line for line in atoms(side) if line[12:16].strip() in ("CA", "CB")]
    assert len(anchor_atoms) == 2 and all(line[77:].strip() == "C" for line in anchor_atoms)
    terms = []
    for a, b in itertools.product(atoms(ligand), anchor_atoms):
        d2 = sum((i-j)**2 for i,j in zip(coordinates(a), coordinates(b)))
        if d2 >= 20.48**2:
            continue
        ta, tb = PARAMETERS[a[77:].strip()], PARAMETERS["C"]
        qa, qb = float(a[70:76]), float(b[70:76])
        def node(r):
            if r >= 20.48:
                return 0.
            rij = ta[0]+tb[0]
            radial = min(max(rij, r-.25), r+.25)
            vdw = math.sqrt(ta[1]*tb[1])*((rij/radial)**12-2*(rij/radial)**6)
            dielectric = -8.5525 + 86.9525/(1+7.7839*math.exp(-86.9525*.003627*r))
            electro = 332*qa*qb*(min(100.,1/(r*dielectric)) if r else 100./dielectric)
            solv = ((ta[2]+.01097*abs(qa))*tb[3]+(tb[2]+.01097*abs(qb))*ta[3])*math.exp(-.5*(r/3.6)**2)
            return .1662*min(100000.,vdw)*(r<8) + .1406*electro + .1322*solv
        index = math.floor(32*d2)
        e = (node(math.sqrt(index/32))+node(math.sqrt((index+1)/32)))/2
        terms.append(e if e <= 0 else 1000*e/(1000+e))
    return sum(terms)


def score(binary, prefix, lig, flex):
    res = subprocess.run([binary, "--scoring", "ad4", "--maps",prefix,"--ligand",lig,"--flex",flex,
                          "--score_only","--cpu","1"],capture_output=True,text=True)
    if res.returncode:
        raise RuntimeError(res.stderr)
    return {key:float(re.search(pat+r"\s*:\s*([-+.\d]+)",res.stdout)[1]) for key,pat in (
        ("score", "Estimated Free Energy of Binding"), ("side", "Ligand - Flex side chains"),
        ("grid","Ligand - Receptor"),("internal", "Final Total Internal Energy"),("torsion","Torsional Free Energy"))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("experiment","old","new"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    out = a.experiment.resolve()
    labels = {(r["query"],r["receiver"],r["seed"],r["rank"]):r for r in csv.DictReader((out/"modes.csv").open())
              if r["arm"]=="old"}
    rows = []
    with tempfile.TemporaryDirectory(dir=out) as scratch:
        lig, flex = [Path(scratch)/name for name in ("ligand.pdbqt","flex.pdbqt")]
        for pool in sorted((out/"pools").glob("*-old.pdbqt")):
            q, r, seed, _ = pool.stem.split("-")
            for rank, block in enumerate(pool.read_text().split("ENDMDL")[:-1],1):
                body = block.split("\nROOT\n",1)[1]
                first, rest = body.split("BEGIN_RES",1)
                ligand, side = "ROOT\n"+first, "BEGIN_RES"+rest
                lig.write_text(ligand)
                flex.write_text(side)
                pred = correction(ligand, side)
                old, new = [score(binary.resolve(),out/"maps"/r/"receptor",lig,flex) for binary in (a.old,a.new)]
                err = max(abs(new[k]-old[k]-pred) for k in ("side","score"))
                if err > .00101 or any(new[k]!=old[k] for k in ("grid","internal","torsion")):
                    raise ValueError(f"{pool} {rank}: paired scores disagree with formula ({err})")
                label = labels[q,r,seed,str(rank)]
                rows.append(dict(query=q,receiver=r,seed=seed,rank=rank,rmsd=label["rmsd"],
                                 min_anchor_A=label["min_anchor_A"],analytic_delta=pred,
                                 rounding_error=err,**{"old_"+k:v for k,v in old.items()},
                                 **{"new_"+k:v for k,v in new.items()}))
    write_csv(out/"fixed.csv",rows)

if __name__ == "__main__":
    main()
