"""Re-score a nested molecular pair on its archived FKBP12 coordinates."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import tempfile

from rdkit import Chem

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
LEGACY = BENCH / "cubic_grid/fkbp_ad4"
HUMAN = BENCH / "cubic_grid/human_panel"
COVER = BENCH / "rigid_pose_coverage/fkbp"
sys.path.insert(0, str(LEGACY))
from analyze import calculate_components, read_map
from run_official import RECEPTORS, SEEDS, TYPES, sha, table


def digest(data):
    return hashlib.sha256(data).hexdigest()


def files(archive, archive_sha, requested):
    if sha(archive) != archive_sha:
        raise ValueError(f"Archive digest mismatch: {archive}")
    result = {}
    with tarfile.open(archive, "r|xz") as tar:
        for member in tar:
            if member.name in requested:
                data = tar.extractfile(member).read()
                if digest(data) != requested[member.name]:
                    raise ValueError(f"Member digest mismatch: {member.name}")
                result[member.name] = data
    if result.keys() != requested.keys():
        raise ValueError(f"Missing archive members: {requested.keys() - result.keys()}")
    return result


def atoms(lines):
    return [line for line in lines if line.startswith(("ATOM", "HETATM"))]


def indexed(lines):
    smiles = next(line[14:] for line in lines if line.startswith("REMARK SMILES "))
    fields = [int(value) for line in lines if line.startswith("REMARK SMILES IDX ")
              for value in line.split()[3:]]
    idx = dict(zip(fields[::2], fields[1::2], strict=True))
    records = {int(line[6:11]): line for line in atoms(lines)}
    assert sorted(idx.keys()) == list(range(1, len(records) + 1))
    assert sorted(idx.values()) == sorted(records)
    return Chem.MolFromSmiles(smiles), [records[idx[i]] for i in sorted(idx)]


def torsion(lines):
    return next(line for line in lines if line.startswith("TORSDOF"))


def models(data):
    return [block.strip().splitlines() for block in data.decode().split("ENDMDL") if block.strip()]


def header_score(lines):
    values = {s.split(":")[0]: float(s.split(":")[1].split()[0]) for s in lines
              if s.startswith("REMARK ") and ":" in s}
    assert abs(values["REMARK INTRA"] - values["REMARK UNBOUND"]) < 1e-10
    score = values["REMARK VINA RESULT"]
    assert abs(score - values["REMARK INTER"] - values["REMARK CONF_INDEPENDENT"]) < .0011
    return score


def audited_maps(receiver, inputs, work):
    names = [f"maps/{receiver}/receptor.{kind}.map" for kind in (*TYPES, "e", "d")]
    archive = LEGACY / f"maps-{receiver}.tar.xz"
    data = files(archive, inputs["archives"][str(archive.relative_to(BENCH.parent))],
                 {name: inputs["files"][name] for name in names})
    fields = {}
    for name, payload in data.items():
        path = work / Path(name).name
        path.write_bytes(payload)
        fields[path.name.split(".")[1]] = read_map(path)
    return fields


def pool_sets():
    manifest = json.loads((COVER / "raw/manifest.json").read_text())
    for filename, item in manifest.items():
        if filename == "replay.tar.xz":
            continue
        arm = filename.split(".")[0]
        names = {f"pools/{arm}/{r}_{n}.pdbqt":(r, n) for r in RECEPTORS for n in (5, 6)}
        data = files(COVER / "raw" / filename, item["sha256"], {n:item["files"][n] for n in names})
        for r in RECEPTORS:
            yield ("coverage", arm, 1047349, r), [data[name] for name, key in names.items() if key[0] == r]
    rows = list(csv.DictReader((LEGACY / "results/cells.csv").read_text().splitlines()))
    hashes = {f"pools/{x['compound']}_{x['receiver']}_{x['seed']}_ad4.pdbqt":x["sha256"]
              for x in rows if x["arm"] == "ad4" and x["compound"] in ("5", "6")}
    archives = {key:value for value,key in [s.split() for s in (LEGACY / "archives.sha256").read_text().splitlines()]}
    legacy = files(LEGACY / "pools.tar.xz", archives["pools.tar.xz"], hashes)
    for seed in SEEDS:
        for r in RECEPTORS:
            yield ("official_legacy", "official", seed, r), [legacy[f"pools/{n}_{r}_{seed}_ad4.pdbqt"] for n in (5,6)]


def match(prep5, prep6):
    m5, a5 = indexed(prep5)
    m6, a6 = indexed(prep6)
    mappings = m6.GetSubstructMatches(m5, useChirality=True, uniquify=False)
    assert len(mappings) == 4
    omitted = sorted(set(range(len(a6))) - set(mappings[0]))
    assert len(a5) == 27 and len(a6) == 33 and len(omitted) == 6
    assert all(set(index) == set(mappings[0]) for index in mappings)
    # Every mapping assigns the same actual type and charge on the same parent sites.
    param = lambda index: {j: a5[i][70:].strip() for i, j in enumerate(index)}
    assert all(param(index) == param(mappings[0]) for index in mappings)
    return a5, a6, mappings[0], omitted


def calculate_cell(key, data, templates, fields):
    panel, arm, seed, receiver = key
    reference5, reference6, indices, omitted = templates
    rows = []
    five, six = map(models, data)
    assert len(five) == len(six) == 9
    for rank, (l5, l6) in enumerate(zip(five, six, strict=True), 1):
        assert int(l5[0].split()[1]) == int(l6[0].split()[1]) == rank
        _, a5 = indexed(l5)
        _, a6 = indexed(l6)
        for actual, ref in ((a5, reference5), (a6, reference6)):
            assert [line[70:] for line in actual] == [line[70:] for line in ref]
        projected = [line[:30] + a6[j][30:54] + line[54:] for line, j in zip(a5, indices, strict=True)]
        operator = "separate" if arm in ("official", "official_high") else "combined"
        def components(lines):
            return calculate_components(lines, fields, operator=operator, uncapped_terms=True)

        score5 = components(l5)
        score6 = components(l6)
        scorep = components([*projected, torsion(l5)])
        added = components([*(a6[j] for j in omitted), "TORSDOF 0"])
        core = components([*(a6[j] for j in indices), "TORSDOF 0"])
        s5, s6, sp = (sum(c.values()) for c in (score5, score6, scorep))
        if rank == 1:
            first5 = s5
        row = dict(panel=panel, arm=arm, seed=seed, receiver=receiver, rank=rank,
                   s5_native=s5, s6_native=s6, s5_projected=sp,
                   substitution=s6-sp, placement=sp-first5, net=s6-first5)
        terms = {}
        for full_key, abbrev in (("type_maps","type"), ("electrostatic","electrostatic"),
                                 ("charge_desolvation","charge_desolvation"), ("cap_adjustment","cap_adjustment")):
            terms["added_"+abbrev] = added[full_key]
            terms["shared_"+abbrev] = core[full_key] - scorep[full_key]
        terms["delta_torsion"] = score6["torsion"] - scorep["torsion"]
        assert abs(sum(terms.values()) - row["substitution"]) < 1e-10
        row.update(terms, reported5=header_score(l5), reported6=header_score(l6),
                   s5_official_operator=sum(calculate_components(l5, fields).values()),
                   s6_official_operator=sum(calculate_components(l6, fields).values()))
        assert max(abs(s5-row["reported5"]), abs(s6-row["reported6"])) < .01
        rows.append(row)
    return rows


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    inputs = json.loads((COVER / "results/inputs.json").read_text())
    lig_archive = HUMAN / "affinity-pools.tar.xz"
    names = [f"affinity/{i}.pdbqt" for i in (5,6)]
    preps = files(lig_archive, inputs["archives"][str(lig_archive.relative_to(BENCH.parent))],
                  {n:inputs["files"][n] for n in names})
    templates = match(*(preps[n].decode().splitlines() for n in names))
    label_path = HUMAN / "affinity/ligands.csv"
    assert sha(label_path) == inputs["evaluation"][str(label_path.relative_to(BENCH.parent))]
    labels = {int(row["compound"]):row for row in csv.DictReader(label_path.read_text().splitlines())}
    for n, expected in ((5,"110"),(6,"12")):
        assert labels[n]["ki_low_nm"] == labels[n]["ki_high_nm"] == expected
        assert ("REMARK SMILES " + labels[n]["smiles"]) in preps[f"affinity/{n}.pdbqt"].decode().splitlines()
    all_rows, firsts = [], []
    with tempfile.TemporaryDirectory() as tmp:
        fields = {r:audited_maps(r, inputs, Path(tmp)) for r in RECEPTORS}
        for key, pools in pool_sets():
            rows = calculate_cell(key, pools, templates, fields[key[-1]])
            first = rows[0]
            firsts.append(dict(**first, favorable_in_6_pool=sum(r["substitution"] < 0 for r in rows),
                               min_5_projection=min(r["s5_projected"] for r in rows),
                               min_5_native=min(r["s5_native"] for r in rows)))
            all_rows.extend(rows)
    assert len(all_rows) == 26*9
    table(out / "poses.csv", all_rows)
    table(out / "first.csv", firsts)
    summary = []
    for panel in ("coverage","official_legacy"):
        for r in RECEPTORS:
            cells = [x for x in firsts if x["panel"] == panel and x["receiver"] == r]
            wrong = [x for x in cells if x["net"] >= 0]
            summary.append(dict(panel=panel, receiver=r, cells=len(cells), correct=len(cells)-len(wrong),
                                wrong_with_favorable_substitution=sum(x["substitution"] < 0 for x in wrong),
                                wrong_with_unfavorable_substitution=sum(x["substitution"] >= 0 for x in wrong),
                                all_favorable_substitutions=sum(x["favorable_in_6_pool"] for x in cells),
                                projection_improves_5_pool=sum(x["min_5_projection"] < x["min_5_native"] for x in cells),
                                n_poses=9*len(cells)))
    table(out / "summary.csv", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    run(parser.parse_args().out)
