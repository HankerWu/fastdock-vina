"""Apply the existing passive refinement trace in a private copy and relink one object."""

import argparse
import os
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("src", "build", "output"):
        p.add_argument(key, type=Path)
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--cxx", default="c++")
    a = p.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True)
    shutil.copytree(a.src.resolve() / "src", out / "src")
    patch = HERE.parents[1] / "pair_gradient/holdout/stage_trace.patch"
    subprocess.run(
        ["git", "apply", patch],
        cwd=out,
        check=True,
        env={**os.environ, "GIT_CEILING_DIRECTORIES": str(out.parent)},
    )
    flags = [
        "-O3",
        "-DNDEBUG",
        "-std=c++14",
        "-pthread",
        "-fvisibility=hidden",
        '-DVERSION="v1.2.7-27-g3c65c0b"',
        f"-I{out / 'src/lib'}",
    ]
    libs = [
        "-lboost_" + n
        for n in ("thread", "serialization", "filesystem", "program_options")
    ]
    if a.boost_prefix:
        prefix = a.boost_prefix.resolve()
        flags.append(f"-I{prefix / 'include'}")
        libs[:0] = [f"-L{prefix / 'lib'}", f"-Wl,-rpath,{prefix / 'lib'}"]
    subprocess.run(
        [a.cxx, *flags, "-c", out / "src/lib/vina.cpp", "-o", out / "vina.o"],
        check=True,
    )
    objects = [
        p
        for p in sorted((a.build.resolve() / "obj").glob("*.o"))
        if p.stem not in ("vina", "split")
    ]
    subprocess.run(
        [a.cxx, *flags, *objects, out / "vina.o", *libs, "-o", out / "vina"], check=True
    )
