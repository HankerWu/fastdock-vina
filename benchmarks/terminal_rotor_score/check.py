"""Check terminal-axis recodings, closure atoms and unchanged complete native searches."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile
from tempfile import TemporaryDirectory

import numpy as np

from assay import BENCH, archive_files, count, save


def scored_field(output, name):
    match = re.search(re.escape(name) + r"\s+:\s+(-?[\d.]+)", output)
    if match is None:
        raise ValueError(name)
    return float(match[1])


def controls(args, work):
    assets = archive_files(BENCH / "rigid_flip/inputs.tar.xz", r"1SJ0/.+\.pdbqt")
    archived = archive_files(BENCH / "silent_torsion/primary/raw.tar.xz", r"runs/1SJ0_1SJ0_202711_3000_parent\.pdbqt")
    receptor = work / "receptor.pdbqt"
    receptor.write_bytes(assets["1SJ0/receptor.pdbqt"])
    first = next(iter(archived.values())).decode().split("ENDMDL")[0].splitlines()[1:]
    center = np.mean([[float(s[q:q + 8]) for q in (30, 38, 46)] for s in first if s.startswith("ATOM")], axis=0)
    header = "GRID_PARAMETER_FILE r\nGRID_DATA_FILE f\nMACROMOLECULE r\nSPACING 10\nNELEMENTS 4 4 4\nCENTER " + " ".join(map(str, center)) + "\n"
    types = {s.split()[-1] for s in first if s.startswith("ATOM")}
    for atom_type in types | {"e", "d"}:
        (work / f"test.{atom_type}.map").write_text(header + "0\n" * 125)
    rows = []
    for removed in range(3):
        axes = {("18", "23"), ("14", "25")} if removed == 2 else ({("18", "23")} if removed else set())
        pose = work / f"frozen{removed}.pdbqt"
        output = []
        for line in first:
            fields = line.split()
            if fields[0] in ("BRANCH", "ENDBRANCH") and tuple(fields[1:]) in axes:
                continue
            output.append(f"TORSDOF {8 - removed}" if fields[0] == "TORSDOF" else line)
        pose.write_text("\n".join(output) + "\n")
        for sf in ("vina", "vinardo", "ad4"):
            inputs = ["--maps", work / "test"] if sf == "ad4" else ["--receptor", receptor]
            for version, binary, probe in (("old", args.parent, args.parent_probe), ("fixed", args.fixed, args.fixed_probe)):
                response = subprocess.check_output([binary, "--ligand", pose, "--cpu", "1", "--seed", "1",
                                                    "--score_only", "--scoring", sf, *inputs], text=True)
                rows.append(dict(sf=sf, arms=version, frozen_axes=removed, nrot=count(probe, pose),
                                 inter=scored_field(response, "(1) Final Intermolecular Energy"),
                                 intra=scored_field(response, "(2) Final Total Internal Energy"),
                                 score=scored_field(response, "Estimated Free Energy of Binding")))
    for sf in ("vina", "vinardo", "ad4"):
        old = [r for r in rows if r["sf"] == sf and r["arms"] == "old"]
        fixed = [r for r in rows if r["sf"] == sf and r["arms"] == "fixed"]
        assert [r["nrot"] for r in old] == [7, 6.5, 6]
        assert [r["nrot"] for r in fixed] == [6, 6, 6]
        if sf == "ad4":
            assert [r["score"] for r in old] == [r["score"] for r in fixed]
            assert all(abs(r["score"] - .2983 * (8 - r["frozen_axes"])) < .00051 for r in fixed)
        else:
            assert len({(r["inter"], r["intra"]) for r in old + fixed}) == 1
            assert len({r["score"] for r in old}) == 3 and len({r["score"] for r in fixed}) == 1
    save(args.out / "recoding.csv", rows)
    ring = BENCH / "fixtures/cyclododecane.pdbqt"
    assert count(args.parent_probe, ring) == count(args.fixed_probe, ring) == 11
    repaired = archive_files(BENCH / "frozen_core/strain_transfer_audit/redock_raw.tar.xz",
                             r"(original|random_coords)\.pdbqt")
    geometry = []
    for filename, content in sorted(repaired.items()):
        path = work / filename
        path.write_bytes(content)
        geometry.append(dict(input=filename, old=count(args.parent_probe, path), fixed=count(args.fixed_probe, path)))
    assert [(r["old"], r["fixed"]) for r in geometry] == [(5, 7), (7, 7)]
    save(args.out / "distorted_input.csv", geometry)


def rmsds(path, meta):
    out = []
    for block in path.read_text().split("ENDMDL")[:-1]:
        xyz = {int(s[6:11]): [float(s[q:q + 8]) for q in (30, 38, 46)]
               for s in block.splitlines() if s.startswith(("ATOM ", "HETATM"))}
        pos = np.array([xyz[n] for n in meta["serial"]])
        out.append(float(np.sqrt(np.min(np.mean(np.sum((np.array(meta["mapped"]) - pos)**2, axis=2), axis=1)))))
    return out


def searches(args, work):
    structure = archive_files(BENCH / "rigid_flip/inputs.tar.xz", r"[^/]+/(?:receptor.pdbqt|ligand.pdbqt|preparation.json)")
    raw = work / "raw"
    raw.mkdir()
    rows = []
    for pdb in ("1SJ0", "1T48", "1T49", "3ERT", "4UND"):
        model = json.loads(structure[f"{pdb}/preparation.json"])
        flags = []
        for key in ("ligand", "receptor"):
            path = work / (key + ".pdbqt")
            path.write_bytes(structure[f"{pdb}/{key}.pdbqt"])
            flags.extend(["--" + key, path])
        for axis, v in zip("xyz", model["center"]):
            flags.extend(["--center_" + axis, str(v), "--size_" + axis, "24"])
        for seed in (202711, 202727):
            for cap in (3000, 12000):
                snapshots = []
                for version, binary in (("old", args.parent), ("fixed", args.fixed), ("upstream", args.upstream)):
                    stem = raw / f"{pdb}-{seed}-{cap}-{version}"
                    path = stem.with_suffix(".pdbqt")
                    log, cost = stem.with_suffix(".log"), stem.with_suffix(".json")
                    cli = [binary, *flags, "--cpu", "1", "--exhaustiveness", "1", "--seed", str(seed),
                           "--max_evals", str(cap), "--num_modes", "9", "--energy_range", "100", "--out", path]
                    if version != "upstream":
                        cli.append("--search_counters")
                    subprocess.run([sys.executable, "-S", BENCH / "rigid_flip/process_cost.py", log, cost, *cli], check=True)
                    coords = [s for s in path.read_bytes().splitlines() if s.startswith((b"ATOM ", b"HETATM", b"ENDMDL"))]
                    counters = [s for s in log.read_text().splitlines() if s.startswith(("Search trajectory:", "Hunt chain:", "Pocket start:", "Native merge:", "Search counters:"))]
                    if version != "upstream":
                        assert counters, stem
                        snapshots.append((coords, counters))
                    default_count = ""
                    if version in ("old", "fixed"):
                        # Compare the actual default 3 kcal/mol consumer window. A common
                        # rescaling can change eligibility without changing the search.
                        default = stem.with_suffix(".default.pdbqt")
                        default_cli = list(cli)
                        i = default_cli.index("--energy_range")
                        del default_cli[i:i + 2]
                        default_cli[default_cli.index("--out") + 1] = default
                        subprocess.run(default_cli, stdout=subprocess.DEVNULL, check=True)
                        default_ds = rmsds(default, model)
                        default_count = len(default_ds)
                        default_coords = [s for s in default.read_bytes().splitlines() if s.startswith((b"ATOM ", b"HETATM", b"ENDMDL"))]
                        assert coords[:len(default_coords)] == default_coords
                    ds = rmsds(path, model)
                    rows.append(dict(pdb=pdb, target=model["target"], seed=seed, cap=cap, arm=version,
                                     modes=len(ds), default_window_modes=default_count, rmsd_top1=ds[0], rmsd_best5=min(ds[:5]), rmsd_oracle=min(ds),
                                     coordinate_sha256=hashlib.sha256(b"\n".join(coords)).hexdigest(),
                                     **json.loads(cost.read_text())))
                assert snapshots[0] == snapshots[1], (pdb, seed, cap)
    save(args.out / "redocking.csv", rows)
    with tarfile.open(args.out / "redocking.tar.xz", "w:xz") as dest:
        dest.add(raw, arcname="raw")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("parent", "fixed", "parent_probe", "fixed_probe", "upstream", "out"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        controls(a, Path(tmp))
        searches(a, Path(tmp))
