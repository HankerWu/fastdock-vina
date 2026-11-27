"""Build a passive one-worker search observer from an already-built arm."""

import argparse
from pathlib import Path
import subprocess


def replace(source, anchor, replacement):
    if source.count(anchor) != 1:
        raise ValueError("Expected one matching source location: " + anchor)
    return source.replace(anchor, replacement)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "objects", "out"):
        p.add_argument(name, type=Path)
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--cxx", default="c++")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    lib = a.source / "src/lib"
    mc = (lib / "monte_carlo.cpp").read_text()
    anchor = '#include "monte_carlo.h"'
    mc = replace(
        mc,
        anchor,
        anchor
        + """
#include <cstdlib>
#include <fstream>
static void record(model& m, const conf& c, const char* stage, sz index) {
    const char* prefix = std::getenv("FASTDOCK_RING_TRACE");
    VINA_CHECK(prefix);
    m.set(c);
    std::ofstream out(std::string(prefix) + "." + stage + ".pdbqt", std::ios::app);
    m.write_model(out, index + 1, "");
}
""",
    )
    anchor = "mixed_reorientation, separate_hunt_chain, centered_rotation, balanced_torsion);"
    if anchor not in mc:
        anchor = "ring ? &ring->allowed_moves() : nullptr);"
    mc = replace(
        mc,
        anchor,
        anchor + '\n        if (step == 0) record(m, candidate.c, "raw", step);',
    )
    anchor = "const int before_hunt = evalcount;"
    mc = replace(
        mc,
        anchor,
        'if (step == 0) record(m, candidate.c, "proposal", step);\n        ' + anchor,
    )
    anchor = "quasi_newton_par(m, p, ig, candidate, hunt_cap, evalcount);"
    mc = replace(mc, anchor, anchor + '\n        record(m, candidate.c, "hunt", step);')
    anchor = "saved.coords = m.get_heavy_atom_movable_coords();"
    mc = replace(
        mc, anchor, 'record(m, saved.c, "second", step);\n                ' + anchor
    )
    anchor = "archive->native_evaluations = evalcount;"
    mc = replace(
        mc,
        anchor,
        anchor
        + '\n    for (sz i = 0; i < out.size(); ++i) record(m, out[i].c, "saved", i);',
    )
    src = a.out / "observed.cpp"
    src.write_text(mc)
    args = [
        a.cxx,
        "-O3",
        "-DNDEBUG",
        "-std=c++14",
        "-pthread",
        "-fvisibility=hidden",
        "-I" + str(lib),
    ]
    linker = []
    if a.boost_prefix:
        args += ["-I" + str(a.boost_prefix / "include")]
        boost_lib = str(a.boost_prefix / "lib")
        linker += ["-L" + boost_lib, "-Wl,-rpath," + boost_lib]
    obj = a.out / "monte_carlo.o"
    subprocess.run(args + ["-c", src, "-o", obj], check=True)
    others = sorted(
        f for f in a.objects.glob("*.o") if f.name not in ("monte_carlo.o", "split.o")
    )
    subprocess.run(
        args
        + [
            obj,
            *others,
            *linker,
            "-lboost_thread",
            "-lboost_serialization",
            "-lboost_filesystem",
            "-lboost_program_options",
            "-o",
            a.out / "vina",
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
