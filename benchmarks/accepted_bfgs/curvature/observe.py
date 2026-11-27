"""Instrument the optimizer on the archived cells and check exported coordinates."""

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import replay


def replace_once(code, match, update):
    if code.count(match) != 1:
        raise ValueError(f"Unexpected source layout at {match}")
    return code.replace(match, update)


def observer(source, build, boost, scratch):
    libsrc = source / "src/lib"
    for file in ("quasi_newton.cpp", "quasi_newton.h", "bfgs.h"):
        shutil.copyfile(libsrc / file, scratch / file)
    header = scratch / "bfgs.h"
    code = "#include <cstdio>\n" + header.read_text()
    # Only used with one CPU. Number minimizations per thread. Interpret call numbers together with the per-trajectory evaluation offset.
    code = replace_once(
        code,
        "original = x;",
        """original = x;
        static thread_local unsigned ordinal = 0;
        const unsigned call = ++ordinal;
        bool negative_first = false;
        const int begin_fg = evalcount;""",
    )
    code = replace_once(
        code,
        "const fl sy = alpha * yp;",
        """const fl sy = alpha * yp;
            if (step == 0 && sy < epsilon_fl) {
                negative_first = true;
                std::fprintf(stderr, "first,%u,%d,%.17g,%.17g,%.17g,%.17g,%.17g\\n",
                    call, begin_fg, initial_energy, energy, alpha, sy, scalar_product(gradient,gradient));
            }""",
    )
    code = replace_once(
        code,
        "if (alpha == 0) {",
        """if (negative_first && step == 1)
                std::fprintf(stderr, "next,%u,%.17g,%.17g,%.17g\\n",
                    call, scalar_product(direction,gradient), alpha, trial_energy);
            if (alpha == 0) {""",
    )
    code = replace_once(
        code,
        "return {energy, restore_coordinates};",
        """if (negative_first)
            std::fprintf(stderr, "finish,%u,%d,%.17g\\n", call, evalcount - begin_fg, energy);
        return {energy, restore_coordinates};""",
    )
    header.write_text(code)
    flags = [
        "c++",
        "-O3",
        "-DNDEBUG",
        "-std=c++14",
        "-pthread",
        "-fvisibility=hidden",
        "-I" + str(libsrc),
        "-I" + str(boost / "include"),
    ]
    obj = scratch / "quasi_newton.o"
    subprocess.run([*flags, "-c", scratch / "quasi_newton.cpp", "-o", obj], check=True)
    others = [
        p
        for p in sorted((build / "obj").glob("*.o"))
        if p.name not in ("quasi_newton.o", "split.o", "vina_split.o")
    ]
    lib = boost / "lib"
    names = ["thread", "serialization", "filesystem", "program_options"]
    if list(lib.glob("libboost_system.*")):
        names.append("system")
    engine = scratch / "trace_vina"
    subprocess.run(
        [
            *flags,
            obj,
            *others,
            "-L" + str(lib),
            "-Wl,-rpath," + str(lib),
            *["-lboost_" + n for n in names],
            "-o",
            engine,
        ],
        check=True,
    )
    return engine


def run(source, build, boost, work, arm):
    meta = json.loads((work / "metadata.json").read_text())
    trace_dir = work / "traces"
    trace_dir.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        temp = Path(tmp).resolve()
        engine = observer(source, build, boost, temp)
        for case, entry in sorted(meta.items()):
            for cap in replay.CAPS:
                for seed in replay.SEEDS:
                    tag = f"{case}_{cap}_{seed}_{arm}"
                    args = replay.dock_args(
                        case, cap, seed, entry, temp / "observed.pdbqt"
                    )
                    with (trace_dir / (tag + ".trace")).open("w") as trace:
                        subprocess.run(
                            [str(engine), *args],
                            cwd=work,
                            check=True,
                            stderr=trace,
                            stdout=subprocess.DEVNULL,
                        )
                    if (temp / "observed.pdbqt").read_bytes() != (
                        work / "runs" / (tag + ".pdbqt")
                    ).read_bytes():
                        raise ValueError(
                            f"Observation changed production output for {tag}"
                        )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "build", "boost", "work"):
        cli.add_argument(name, type=Path)
    cli.add_argument("arm", choices=("parent", "fixed"))
    a = cli.parse_args()
    run(
        a.source.resolve(),
        a.build.resolve(),
        a.boost.resolve(),
        a.work.resolve(),
        a.arm,
    )
