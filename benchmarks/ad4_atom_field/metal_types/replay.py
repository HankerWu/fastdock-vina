"""Reproduce the AD4 probe-type alias at two deposited glucokinase ion positions."""
import argparse
import csv
from pathlib import Path
import re
import subprocess
import tarfile


ARCHIVE = Path(__file__).resolve().parents[2] / "ligand_hydroxyl/glucose/evidence.tar.xz"
SCORE = re.compile(r"Estimated Free Energy of Binding\s+:\s+([-\d.]+)")


def field_center(directory, channel):
    path = directory / f"site.{channel}.map"
    lines = path.read_text().splitlines()
    assert lines[4].split() == ["NELEMENTS", "16", "16", "16"]
    values = lines[6:]
    assert len(values) == 17 ** 3
    return float(values[len(values) // 2])


def map_inputs(directory, pdb, element):
    with tarfile.open(ARCHIVE) as archive:
        deposited = archive.extractfile(f"raw/{pdb}.pdb").read().decode()
        receptor = archive.extractfile(f"{pdb}/receptor.pdbqt").read()
    line, = [s for s in deposited.splitlines() if s.startswith("HETATM")
             and s[21:26] == "A 600" and s[76:78].strip().upper() == element.upper()]
    (directory / "receptor.pdbqt").write_bytes(receptor)
    types = sorted({s.split()[-1] for s in receptor.decode().splitlines()
                    if s.startswith(("HETATM", "ATOM  "))})
    assert element not in types
    center = " ".join(line[30:54].split())
    (directory / "site.gpf").write_text(
        f"npts 16 16 16\ngridfld site.maps.fld\nspacing 0.375\nreceptor_types {' '.join(types)}\n"
        f"ligand_types C Fe Zn\nreceptor receptor.pdbqt\ngridcenter {center}\nsmooth 0.5\n"
        "map site.C.map\nmap site.Fe.map\nmap site.Zn.map\nelecmap site.e.map\n"
        "dsolvmap site.d.map\ndielectric -42\n")
    return line


def run(cli, args, log):
    r = subprocess.run([cli, *args, "--cpu", "1", "--seed", "331", "--score_only"],
                       capture_output=True, text=True)
    log.write_text(r.stdout + r.stderr)
    match = SCORE.search(r.stdout)
    return r.returncode, float(match[1]) if match else None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("parent", "fixed", "autogrid", "out"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rows = []
    entry_points = []
    for pdb, element in (("1V4S", "Na"), ("3IDH", "K")):
        directory = a.out / pdb
        directory.mkdir(exist_ok=True)
        line = map_inputs(directory, pdb, element)
        subprocess.run([a.autogrid.resolve(), "-p", "site.gpf", "-l", "site.glg"],
                       cwd=directory, check=True, capture_output=True)
        for probe, charge in ((element, 1), (element, 0), ("C", 0), ("Fe", 2), ("Zn", 2)):
            stem = f"{probe}_{charge}"
            ligand = directory / f"{stem}.pdbqt"
            ligand.write_text(f"ROOT\n{line[:66]}    {charge:6.3f} {probe}\nENDROOT\nTORSDOF 0\n")
            typed = field_center(directory, probe if probe != element else "e")
            phi, dsolv = (field_center(directory, t) for t in ("e", "d"))
            raw = typed + charge * phi + abs(charge) * dsolv
            expected = raw if raw <= 0 else raw / (1 + raw / 1000)
            args = ["--ligand", ligand, "--scoring", "ad4", "--maps", directory / "site"]
            old_rc, old = run(a.parent, args, directory / f"{stem}-parent.log")
            new_log = directory / f"{stem}-fixed.log"
            new_rc, new = run(a.fixed, args, new_log)
            assert old_rc == 0 and abs(old - expected) < 0.00051
            if probe == element:
                assert new_rc != 0 and f"No AD4 atom parameters for {probe}." in new_log.read_text()
            else:
                assert new_rc == 0 and old == new
            rows.append(dict(pdb=pdb, site="A600", probe=probe, charge=charge, scoring="ad4",
                             atom_map_source="e" if probe == element else probe, phi=phi, D=dsolv,
                             typed=typed, reconstructed=round(expected, 6), parent_score=old,
                             fixed_score=new, fixed_exit=new_rc))
        args = ["--ligand", directory / f"{element}_1.pdbqt", "--scoring", "ad4",
                "--maps", directory / "site", "--cpu", "1", "--seed", "331",
                "--exhaustiveness", "1", "--max_evals", "500"]
        for mode, options in (("local", ["--local_only"]), ("dock", [])):
            for label, engine in (("parent", a.parent), ("fixed", a.fixed)):
                stem = directory / f"{mode}-{label}"
                r = subprocess.run([engine, *args, *options, "--out", stem.with_suffix(".pdbqt")],
                                   capture_output=True, text=True)
                stem.with_suffix(".log").write_text(r.stdout + r.stderr)
                if label == "parent":
                    assert r.returncode == 0
                else:
                    assert r.returncode != 0 and f"No AD4 atom parameters for {element}." in r.stderr
                entry_points.append(dict(pdb=pdb, probe=element, charge=1, mode=mode,
                                         arm=label, exit_code=r.returncode))
        # Same original +1 ion under the empirical XS metal type is a different valid API path.
        for sf in ("vina", "vinardo"):
            args = ["--ligand", directory / f"{element}_1.pdbqt", "--scoring", sf,
                    "--receptor", directory / "receptor.pdbqt"]
            old_rc, old = run(a.parent, args, directory / f"{sf}-parent.log")
            new_rc, new = run(a.fixed, args, directory / f"{sf}-fixed.log")
            assert old_rc == new_rc == 0 and old == new
            rows.append(dict(pdb=pdb, site="A600", probe=element, charge=1, scoring=sf,
                             parent_score=old, fixed_score=new, fixed_exit=new_rc))
    for name, records in (("scores", rows), ("entry_points", entry_points)):
        with (a.out / f"{name}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, records[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)


if __name__ == "__main__":
    main()
