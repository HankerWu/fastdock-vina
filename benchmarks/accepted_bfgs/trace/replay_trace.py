"""Replay the 2X6E/AT9283 line search from archived inputs using an instrumented local binary."""

import argparse
import json
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path


def extract(archive, member, directory):
    target = directory / Path(member).name
    with tarfile.open(archive) as contents, contents.extractfile(member) as stream:
        target.write_bytes(stream.read())
    return target


def instrument(source, directory):
    for name in ("quasi_newton.cpp", "quasi_newton.h", "bfgs.h"):
        shutil.copyfile(source / name, directory / name)
    header = directory / "bfgs.h"
    code = "#include <cstdio>\n" + header.read_text()
    additions = {
        "++evalcount;\n\t\tif (StrictArmijo)": (
            '++evalcount;\n\t\tstd::fprintf(stderr, "trial,%d,%d,%.17g,%.17g,%.17g,%.17g\\n", '
            "int(StrictArmijo), evalcount, alpha, pg, f0, f1);\n\t\tif (StrictArmijo)"
        ),
        "const fl sy = alpha * yp;": (
            'const fl sy = alpha * yp;\n\t\t\tstd::fprintf(stderr, "curvature,%d,%u,%.17g,%.17g\\n", '
            "int(StrictArmijo), unsigned(step), alpha, sy);"
        ),
        "std::swap(x, original);": (
            'std::fprintf(stderr, "RESTORE,%d,%.17g,%.17g\\n", int(StrictArmijo), energy, initial_energy);'
            "\n\t\t\t\tstd::swap(x, original);"
        ),
    }
    for anchor, replacement in additions.items():
        if code.count(anchor) != 1:
            raise ValueError(f"Unexpected number of source matches for {anchor!r}")
        code = code.replace(anchor, replacement)
    header.write_text(code)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("source", "build", "boost", "docking", "probe", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--cxx", default="c++")
    args = parser.parse_args()
    source = args.source.resolve() / "src/lib"
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        instrument(source, directory)
        lib = args.boost.resolve() / "lib"
        flags = [
            args.cxx,
            "-O3",
            "-DNDEBUG",
            "-std=c++14",
            "-pthread",
            "-fvisibility=hidden",
            "-I" + str(source),
            "-I" + str(args.boost.resolve() / "include"),
        ]
        patched_object = directory / "quasi_newton.o"
        subprocess.run(
            [*flags, "-c", directory / "quasi_newton.cpp", "-o", patched_object],
            check=True,
        )
        objects = [
            p
            for p in sorted((args.build / "obj").glob("*.o"))
            if p.name not in ("quasi_newton.o", "split.o", "vina_split.o")
        ]
        libraries = [
            "-lboost_" + name
            for name in ("thread", "serialization", "filesystem", "program_options")
        ]
        if list(lib.glob("libboost_system.*")):
            libraries.append("-lboost_system")
        engine = directory / "trace_vina"
        subprocess.run(
            [
                *flags,
                *objects,
                patched_object,
                "-L" + str(lib),
                "-Wl,-rpath," + str(lib),
                *libraries,
                "-o",
                engine,
            ],
            check=True,
        )
        receptor = extract(args.docking, "prepared/2X6E.pdbqt", directory)
        center = json.loads(
            extract(args.docking, "prepared/2X6E.json", directory).read_text()
        )["center"]
        ligand = extract(
            args.probe, "raw/2X6E_at9283_sdf_legal_placed.pdbqt", directory
        )
        bonds = extract(args.probe, "raw/at9283_chemical.bonds", directory)
        cmd = [
            "--receptor",
            receptor,
            "--ligand",
            ligand,
            "--ligand_bonds",
            bonds,
            "--cpu",
            "1",
            "--seed",
            "60913",
            "--local_only",
        ]
        for axis, value in zip("xyz", center):
            cmd += ["--center_" + axis, str(value), "--size_" + axis, "26"]
        for mode in ("native", "armijo"):
            options = [] if mode == "native" else ["--refinement_optimizer", mode]
            with (args.output / (mode + ".trace")).open("w") as log:
                output = args.output / (mode + ".pdbqt")
                subprocess.run(
                    [engine, *cmd, *options, "--out", output],
                    stderr=log,
                    stdout=subprocess.DEVNULL,
                    check=True,
                )
                reference = directory / "production.pdbqt"
                subprocess.run(
                    [args.build / "vina", *cmd, *options, "--out", reference],
                    stdout=subprocess.DEVNULL,
                    check=True,
                )
                if output.read_bytes() != reference.read_bytes():
                    raise ValueError(
                        "Instrumentation changed the production PDBQT for " + mode
                    )


if __name__ == "__main__":
    main()
