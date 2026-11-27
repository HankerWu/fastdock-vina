"""Recompile only the two traced translation units of pinned parent/linear builds."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent


def build(exports, arm, output, prefix, cxx, overlay=None):
    base = exports / arm
    dest = output / arm
    dest.mkdir(parents=True)
    shutil.copytree(base / "src", dest / "src")
    subprocess.run(["git", "apply", HERE / "mc_trace.patch"], cwd=dest,
                   env={**os.environ, "GIT_CEILING_DIRECTORIES": str(dest.parent)}, check=True)
    if overlay:
        subprocess.run(["git", "apply", overlay.resolve()], cwd=dest, check=True,
                       env={**os.environ, "GIT_CEILING_DIRECTORIES": str(dest.parent)})
    flags = ["-O3", "-DNDEBUG", "-std=c++14", "-pthread", "-fvisibility=hidden",
             '-DVERSION="v1.2.7-27-g3c65c0b"', f"-I{dest/'src/lib'}"]
    if arm == "linear":
        flags.append("-DFASTDOCK_EXACT_PAIR_GRADIENT")
    if prefix:
        flags.append(f"-I{prefix/'include'}")
    objects = []
    for obj in sorted((base / "build/obj").glob("*.o")):
        if obj.stem in ("monte_carlo", "vina"):
            source = dest / "src/lib" / (obj.stem + ".cpp")
            obj = dest / obj.name
            subprocess.run([cxx, *flags, "-c", source, "-o", obj], check=True)
        if obj.stem != "split":
            objects.append(obj)
    libs = [f"-lboost_{name}" for name in ("thread","serialization","filesystem","program_options")]
    if prefix:
        libs[:0] = [f"-L{prefix/'lib'}",f"-Wl,-rpath,{prefix/'lib'}"]
    subprocess.run([cxx,*flags,*objects,*libs,"-o",dest/"vina"],check=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("exports", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--cxx", default=os.environ.get("CXX","c++"))
    p.add_argument("--overlay", type=Path, help="Optional additional diagnostic patch after the stage observer")
    p.add_argument("--arm", choices=("parent","linear"), required=True)
    a = p.parse_args()
    build(a.exports.resolve(), a.arm, a.output.resolve(), a.boost_prefix, a.cxx, a.overlay)
