"""Instrument only the first map evaluation; capture exact coordinates and forces then stop."""

import argparse
import csv
import os
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory

MARKER = "AD4FIRST"
TRACE = r"""
    std::cout << std::setprecision(17) << "AD4FIRST total " << ENERGY << '\n';
    VINA_FOR(i, m.num_movable_atoms()) {
        const auto& a = m.atoms[i];
        std::cout << "AD4FIRST atom " << i << ' ' << a.ad << ' ' << a.charge;
        for (sz j=0; j<3; ++j) std::cout << ' ' << m.coords[i][j];
        for (sz j=0; j<3; ++j) std::cout << ' ' << m.minus_forces[i][j];
        std::cout << '\n';
    }
    std::exit(0);
"""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "objects", "work", "out"):
        p.add_argument(name, type=Path)
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--arm", required=True)
    a = p.parse_args()
    source = a.source.resolve()
    a.out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = (source / "src/lib/ad4cache.cpp").read_text()
        start = src.index("fl ad4cache::eval_deriv(")
        end = src.index("bool ad4cache::is_in_grid", start)
        block = src[start:end]
        match = re.search(r"\breturn (e|energy);\s*\}", block)
        if not match:
            raise ValueError("unrecognized map derivative routine")
        block = block[: match.start()] + TRACE.replace("ENERGY", match[1]) + "}\n\n"
        unit = tmp / "ad4cache.cpp"
        unit.write_text("#include <iomanip>\n" + src[:start] + block + src[end:])
        cxx = os.environ.get("CXX", "c++")
        inc = ["-I" + str(source / "src/lib")]
        libs = []
        if a.boost_prefix:
            inc.append("-I" + str(a.boost_prefix / "include"))
            libs.extend(
                [
                    "-L" + str(a.boost_prefix / "lib"),
                    "-Wl,-rpath," + str(a.boost_prefix / "lib"),
                ]
            )
        libs.extend(
            [
                "-lboost_" + key
                for key in ("thread", "serialization", "filesystem", "program_options")
            ]
        )
        obj = tmp / "ad4cache.o"
        subprocess.run(
            [
                cxx,
                "-O3",
                "-DNDEBUG",
                "-std=c++14",
                "-pthread",
                *inc,
                "-c",
                unit,
                "-o",
                obj,
            ],
            check=True,
        )
        exe = tmp / "probe"
        subprocess.run(
            [
                cxx,
                "-pthread",
                obj,
                *[o for o in sorted(a.objects.glob("*.o")) if o.name != "ad4cache.o"],
                *libs,
                "-o",
                exe,
            ],
            check=True,
        )
        for pdb in ("1FKG", "1FKH", "2BTS", "2PRG"):
            for seed in (14591, 30553):
                proc = subprocess.run(
                    [
                        exe,
                        "--scoring",
                        "ad4",
                        "--maps",
                        a.work / "maps" / pdb / "receptor",
                        "--ligand",
                        a.work / pdb / "ligand.pdbqt",
                        "--cpu",
                        "1",
                        "--exhaustiveness",
                        "1",
                        "--max_evals",
                        "3000",
                        "--num_modes",
                        "5",
                        "--seed",
                        str(seed),
                        "--out",
                        tmp / "unused.pdbqt",
                    ],
                    text=True,
                    capture_output=True,
                    check=True,
                )
                lines = [
                    l[l.index(MARKER) :].split()
                    for l in proc.stdout.splitlines()
                    if MARKER in l
                ]
                if not lines or lines[0][1] != "total":
                    raise ValueError(proc.stdout)
                with (a.out / f"{pdb}-{seed}-{a.arm}.csv").open("w", newline="") as f:
                    w = csv.writer(f, lineterminator="\n")
                    w.writerow("atom type charge x y z gx gy gz map_total".split())
                    for l in lines[1:]:
                        assert l[1] == "atom"
                        w.writerow([*l[2:], lines[0][2]])


if __name__ == "__main__":
    main()
