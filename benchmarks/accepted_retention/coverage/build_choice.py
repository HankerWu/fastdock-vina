"""Recompile the passive per-skip observer against the final selector's own ABI."""

import argparse
import os
import shutil
import subprocess
from pathlib import Path


def build(a):
    dest = a.output.resolve()
    dest.mkdir(parents=True)
    lib = dest / "src/lib"
    lib.mkdir(parents=True)
    name = "monte_carlo.cpp"
    source = a.source.resolve() / "src/lib"
    shutil.copy2(source / name, lib / name)
    patch = a.patch.resolve()
    subprocess.run(
        ["git", "apply", patch],
        cwd=dest,
        check=True,
        env=dict(os.environ, GIT_CEILING_DIRECTORIES=str(dest.parent)),
    )
    flags = [
        "-O3",
        "-DNDEBUG",
        "-std=c++14",
        "-pthread",
        "-fvisibility=hidden",
        f"-I{source}",
    ]
    link = []
    if a.boost_prefix:
        path = a.boost_prefix.resolve()
        flags.append(f"-I{path / 'include'}")
        link = [f"-L{path / 'lib'}", f"-Wl,-rpath,{path / 'lib'}"]
    obj = dest / "monte_carlo.o"
    subprocess.run([a.cxx, *flags, "-c", lib / name, "-o", obj], check=True)
    objects = [
        p
        for p in sorted(a.objects.glob("*.o"))
        if p.stem not in ("monte_carlo", "split")
    ]
    if not any(p.stem == "main" for p in objects):
        raise ValueError("CLI object set required")
    libs = [
        f"-lboost_{x}"
        for x in ("thread", "serialization", "filesystem", "program_options")
    ]
    subprocess.run(
        [a.cxx, "-pthread", obj, *objects, *link, *libs, "-o", dest / "observer"],
        check=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "objects", "output"):
        p.add_argument(name, type=Path)
    p.add_argument(
        "--patch", type=Path, default=Path(__file__).with_name("choice_observer.patch")
    )
    p.add_argument("--cxx", default="c++")
    p.add_argument("--boost-prefix", type=Path)
    build(p.parse_args())
