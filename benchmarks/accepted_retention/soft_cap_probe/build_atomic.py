"""Build the two diagnostic translation units against an exported candidate arm."""

import argparse
import concurrent.futures
import os
import shutil
import subprocess
from pathlib import Path

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("source", "objects", "output"):
        p.add_argument(k, type=Path)
    p.add_argument("--cxx", default="c++")
    p.add_argument("--boost-prefix", type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True)
    lib = out / "src/lib"
    shutil.copytree(a.source.resolve() / "src/lib", lib)
    subprocess.run(
        ["git", "apply", Path(__file__).with_name("explicit_grid.patch").resolve()],
        cwd=out,
        env={**os.environ, "GIT_CEILING_DIRECTORIES": str(out.parent)},
        check=True,
    )
    flags = [
        "-O3",
        "-DNDEBUG",
        "-std=c++14",
        "-pthread",
        "-fvisibility=hidden",
        '-DVERSION="v1.2.7-27-g3c65c0b"',
        f"-I{lib}",
    ]
    link = []
    if a.boost_prefix:
        prefix = a.boost_prefix.resolve()
        flags.append(f"-I{prefix / 'include'}")
        link = [f"-L{prefix / 'lib'}", f"-Wl,-rpath,{prefix / 'lib'}"]

    def compile_one(name):
        subprocess.run(
            [a.cxx, *flags, "-c", lib / (name + ".cpp"), "-o", out / (name + ".o")],
            check=True,
        )

    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        list(pool.map(compile_one, ["monte_carlo", "vina"]))
    objs = [
        p
        for p in sorted(a.objects.glob("*.o"))
        if p.stem not in ("monte_carlo", "vina", "split")
    ]
    if not any(p.stem == "main" for p in objs):
        raise ValueError("CLI object set is required")
    libs = [
        "-lboost_" + x
        for x in ("thread", "serialization", "filesystem", "program_options")
    ]
    subprocess.run(
        [
            a.cxx,
            "-pthread",
            *objs,
            out / "monte_carlo.o",
            out / "vina.o",
            *link,
            *libs,
            "-o",
            out / "observer",
        ],
        check=True,
    )
