"""Verify membership by energy-priority set operations and exact cartesian crystal labels."""
import csv
import json
import lzma
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stage_local_hunt"))
from coverage import specs


def verify(base, output):
    metadata = list(csv.DictReader((base / "results/arrivals.csv").open()))
    counts = dict(cells=0, arrivals=0, reference_coordinates=0, max_label_difference_A=0.0)
    with tempfile.TemporaryDirectory() as tmp:
        for _, job in specs(Path(tmp)):
            _, inputs, q, r, ref, seed, _ = job
            serials = [int(s[6:11]) for s in (inputs/q["pdb"]/"ligand.pdbqt").read_text().splitlines()
                       if s.startswith(("ATOM  ","HETATM")) and s[77:].strip() not in {"H","HD","HS"}]
            order = [serials.index(i) for i in q["serial"]]
            for arm in ("stable","candidate"):
                key = f"reciprocal-{q['pdb']}-{r['pdb']}-{seed}-{arm}"
                live = lzma.decompress((base / "raw" / (key+".tsv.xz")).read_bytes()).decode().splitlines()
                old = lzma.decompress((HERE.parent/"stage_local_hunt/coverage"/(key+".xz")).read_bytes())
                old_fields = (line.split("\t") for line in old.decode().splitlines())
                saved = {f[3]: np.array(f[8].split(), dtype=float).reshape(-1,3)/1000.
                         for f in old_fields if f[0]=="D"}
                ranked = [r for r in metadata if r["key"]==key]
                radius, cap = (float(x) for x in live[0].split("\t")[1:])
                coords, energy, residents = [], [], set()
                for line in live[1:]:
                    ev,*fields = line.split("\t")
                    if ev == "C":
                        i = len(coords)
                        xyz = np.array(fields[4:],dtype=float).reshape(-1,3)
                        coords.append(xyz)
                        energy.append(float(fields[2]))
                        rounded = saved[fields[3]]
                        # Recover the original atom permutation by unambiguous matching
                        # of exact positions to the old rounded coordinates of this state.
                        matches = np.max(np.abs(rounded[:,None] - xyz), axis=2) <= .00050000001
                        assert np.all(matches.sum(axis=1)==1) and np.all(matches.sum(axis=0)==1)
                        ordered = xyz[matches.argmax(axis=1)][order]
                        exact = ((np.linalg.norm(ordered-ref,axis=2)**2).mean(axis=1)**.5).min()
                        delta = abs(exact-float(ranked[i]["crystal_A"]))
                        assert delta < .00087 and (exact<=2)==(float(ranked[i]["crystal_A"])<=2)
                        counts["reference_coordinates"] += len(ordered)
                        counts["max_label_difference_A"] = max(counts["max_label_difference_A"],delta)
                        if key == "reciprocal-1RBP-3FMZ-85237-candidate" and i in (18,20):
                            counts[f"primary_ID{i}_exact_crystal_A"] = exact
                        # Independent set formulation with the native lexicographic
                        # priority (lower energy, earlier arrival) for exact ties.
                        better = {j for j in residents if (energy[j],j)<(energy[i],i)}
                        nearby = {j for j in residents if
                                  np.linalg.norm(xyz-coords[j])/len(xyz)**.5 < radius}
                        if len(better)<cap and not (better & nearby):
                            residents -= nearby
                            residents.add(i)
                            if len(residents)>cap:
                                residents.remove(max(residents,key=lambda j:(energy[j],j)))
                        counts["arrivals"] += 1
                    elif ev in {"M","I"}:
                        assert list(map(int,fields)) == sorted(residents,key=lambda j:(energy[j],j))
                    elif ev=="G":
                        assert len(set(fields))==len(fields)
                counts["cells"] += 1
    output.write_text(json.dumps(counts,indent=2)+"\n")


if __name__ == "__main__":
    verify(HERE,Path(sys.argv[1]))
