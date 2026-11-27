"""Recheck 9GIO using archived CCD atom-name matches, Gemmi registration and images."""
import argparse
import csv
import json
import sys
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from verify import transform, coords


def check(work, results):
    old = {p: json.loads((work / p / "metadata.json").read_text()) for p in ("4W9H", "5LLI")}
    structs = {p: gemmi.read_structure(str(work / "raw" / (p+".pdb"))) for p in old}
    st = gemmi.read_structure(str(HERE / "9gio.cif.gz"))
    st.setup_entities()
    cn = st[0]["C"]
    lig = next(r for r in cn if r.name == "3JF")
    ligand = {a.name: tuple(a.pos) for a in lig if a.occ and not a.element.is_hydrogen}
    raw_ligand = next(r for r in structs["4W9H"][0]["C"] if r.name == "3JF")
    old_names = {tuple(a.pos): a.name for a in raw_ligand}
    maps = np.array([[ligand[old_names[tuple(p)]] for p in eq] for eq in old["4W9H"]["crystal"]])
    ca = {"ca": {}}
    for res in cn:
        for atom in res:
            k = f"C{res.seqid.num:4}{res.seqid.icode}{res.name}"
            if k in old["4W9H"]["ca"] and k in old["5LLI"]["ca"] and atom.name == "CA":
                assert k not in ca["ca"]
                ca["ca"][k] = list(atom.pos)
    assert len(ca["ca"]) == 55
    center = gemmi.Position(*np.array(list(ligand.values())).mean(0))
    # Only the C-centred trimer is present in this asymmetric unit. Query a
    # sphere smaller than its shortest plane spacing / 2 to avoid duplicates.
    radius = 27.
    reciprocal = st.cell.reciprocal()
    assert 2 * radius < 1 / max(reciprocal.a, reciprocal.b, reciprocal.c)
    ns = gemmi.NeighborSearch(st[0], st.cell, radius).populate(include_h=False)
    cloud = []
    for mark in ns.find_atoms(center, radius=radius):
        cra = mark.to_cra(st[0])
        if cra.residue.entity_type != gemmi.EntityType.Polymer or not cra.atom.occ:
            continue
        im = st.cell.find_nearest_pbc_image(center, cra.atom.pos, mark.image_idx)
        if not im.same_asu():
            cloud.append(list(st.cell.find_nearest_pbc_position(center, cra.atom.pos, mark.image_idx)))
    from scipy.spatial.distance import cdist
    with (results / "forms.csv").open() as f:
        rows = list(csv.DictReader(f))
    errors, geom = [], []
    for row in rows:
        if row["group"] not in ("blind", "witness", "explicit_graph", "deposited_9GIO"):
            continue
        if row["group"] == "deposited_9GIO":
            xyz = maps[0]
        else:
            p = coords(dict(row, query="4W9H", receiver="5LLI", frame="receiver"), work, old, structs)
            if row["group"] in ("blind","explicit_graph","witness"):
                # coords reads atoms in PDBQT serial order. Switch to query SMILES order.
                serials = old["4W9H"]["serial"]
                order = {s:i for i,s in enumerate(sorted(serials))}
                p = p[[order[s] for s in serials]]
            xyz = transform(p, old["5LLI"], ca)
        e = np.sqrt(np.square(maps-xyz).sum(axis=-1).mean(axis=-1)).min()
        geom.append(abs(e-float(row["rmsd_9GIO_A"])))
        dist = cdist(xyz, cloud)
        assert dist.min() + np.linalg.norm(xyz-np.array(list(center)), axis=1).max() < radius
        errors.append(abs(dist.min()-float(row["closest_9GIO_mate_A"])))
        assert int((dist<2.3).sum()) == int(row["pairs_lt_2p3"])
    assert max(errors) < 1e-9 and max(geom) < 1e-9
    print(json.dumps({"checked_rows": len(errors), "max_RMSD_difference_A": max(geom),
                      "max_nearest_distance_difference_A": max(errors)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("results", type=Path)
    a = parser.parse_args()
    check(a.work, a.results)
