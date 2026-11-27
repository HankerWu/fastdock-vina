"""Independently exercise cut sides, rooting and kinematics on the compiled C++ model."""

import argparse
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))


def run(probe, *paths):
    count, nlig, *movement = subprocess.check_output([probe, *paths], text=True).split()
    return float(count), int(nlig), list(map(float, movement))


def atom(serial, typ, point):
    x, y, z = point
    return (
        f"ATOM  {serial:5d}  {typ[0]:<2}  UNL     1    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00     0.000 {typ}"
    )


def synthetic(work, probe):
    cases = {}
    for name, elt, x2, y2, gap, count, visible in [
        ("normal", "C", 1.4, 1.2, 1.4, 1, True),
        ("stretched_branch", "C", 3.5, 1.2, 3.5, 1, True),
        ("terminal_hydrogen", "HD", 1.4, 1.2, 1.4, 0, False),
        ("terminal_water_proxy", "W", 1.4, 1.2, 1.4, 0, False),
        ("collinear_tip", "NA", 2.55, 0, 1.4, 1, False),
    ]:
        points = [[0.0, 1.35, 0.0], [0.0, 0.0, 0.0], [gap, 0.0, 0.0], [x2, y2, 0.0]]
        atos = [
            atom(i + 1, t, p)
            for i, (t, p) in enumerate(zip(["C", "C", "C", elt], points))
        ]
        for root, axes, order in [
            ("forward", (2, 3), [0, 1, 2, 3]),
            ("reversed", (3, 2), [3, 2, 1, 0]),
        ]:
            axis = f"{axes[0]} {axes[1]}"
            txt = [
                "ROOT",
                atos[order[0]],
                atos[order[1]],
                "ENDROOT",
                "BRANCH " + axis,
                atos[order[2]],
                atos[order[3]],
                "ENDBRANCH " + axis,
                "TORSDOF 1",
            ]
            path = work / (name + "_" + root + ".pdbqt")
            path.write_text("\n".join(txt) + "\n")
            observed, n, move = run(probe, path)
            assert observed == count and n == 1 and len(move) == 1
            assert (move[0] > 1e-6) == visible, (name, root, move)
            cases[name, root] = (
                path,
                dict(
                    case=name,
                    root=root,
                    torsions=observed,
                    max_internal_distance_change=move[0],
                ),
            )
    files = [
        cases["normal", "forward"][0],
        cases["terminal_hydrogen", "reversed"][0],
        cases["collinear_tip", "forward"][0],
    ]
    # Overlapping atom numbering is valid across separately parsed ligand models.
    assert run(probe, *files)[:2] == (2.0, 3)
    assert run(probe, *reversed(files))[:2] == (2.0, 3)
    return [x[1] for x in cases.values()]


def independently_count(pdbqt):
    """Use the flat file's branch intervals. Leave geometric radii out of the count."""
    atoms = {}
    intervals = []
    open_branches = []
    for ln in pdbqt.decode().splitlines():
        fields = ln.split()
        if not fields:
            continue
        if fields[0] == "BRANCH":
            begin = len(atoms)
            open_branches.append((int(fields[1]), int(fields[2]), begin))
        elif fields[0] == "ENDBRANCH":
            a, b, begin = open_branches.pop()
            assert (a, b) == tuple(map(int, fields[1:]))
            intervals.append((a, b, begin, len(atoms)))
        elif fields[0] in ("ATOM", "HETATM"):
            serial = int(ln[6:11])
            assert serial not in atoms
            atoms[serial] = fields[-1]
    assert not open_branches

    def site(t):
        return t not in ("H", "HD", "W")

    site_flags = [site(typ) for typ in atoms.values()]
    sites = sum(site_flags)
    axes = 0
    for a, b, i, j in intervals:
        component = sum(site_flags[i:j])
        clean_ends = all(
            site(atoms[s]) and atoms[s] not in ("G0", "G1", "G2", "G3") for s in (a, b)
        )
        axes += int(clean_ends and min(component, sites - component) >= 2)
    return axes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("probe", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--synthetic-only",
        action="store_true",
        help="Run the lightweight CI controls without the archived assay data",
    )
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        rows = synthetic(work, a.probe)
        with (a.out / "controls.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        ring = HERE.parents[1] / "fixtures/cyclododecane.pdbqt"
        count, n, distances = run(a.probe, ring)
        assert (count, n) == (11, 1) and sum(d > 1e-10 for d in distances) == 9
        if a.synthetic_only:
            print(
                "10 rooted controls, both multi-ligand orders and the macrocycle passed"
            )
            return
        from assay import assay_inputs, save

        files = {(t, key, rec): lig for t, key, _, rec, _, _, lig, _ in assay_inputs()}
        files["macro", "cyclododecane", "fixture"] = ring.read_bytes()
        full = []
        for (target, lig, rec), contents in sorted(files.items()):
            path = work / "query.pdbqt"
            path.write_bytes(contents)
            expected = independently_count(contents)
            observed, n, finite = run(a.probe, path)
            assert observed == expected and n == 1, (
                target,
                lig,
                rec,
                observed,
                expected,
            )
            full.append(
                dict(
                    target=target,
                    ligand=lig,
                    receiver=rec,
                    pdbqt_count=expected,
                    cpp_count=observed,
                    tree_coordinates=len(finite),
                    exact_fixed_internal_distances=sum(x < 1e-10 for x in finite),
                    max_changes_angstrom=" ".join(f"{x:.7g}" for x in finite),
                )
            )
        save(a.out / "prepared.csv", full)
    print(f"{len(full)} coordinate files matched PDBQT counts")


if __name__ == "__main__":
    main()
