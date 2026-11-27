"""Build the archived research commit and optional fixed-start observers."""
import argparse
import hashlib
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1] / "pair_gradient"))
from build_arms import export_src, git

BASE = "b26101b9028b24e488718024c8ba01e3c680d813"


def write_wrappers(root):
    for arm in ("native", "halving", "quadratic"):
        wrapper = root / arm / "vina"
        wrapper.parent.mkdir()
        wrapper.write_text(
            "#!/bin/sh\nexport FASTDOCK_RESEARCH_BACKTRACK=" + arm
            + "\nexec " + shlex.quote(str(root / "build/vina")) + ' "$@"\n'
        )
        wrapper.chmod(0o755)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--first-hunt-replay", action="store_true", help="build the fixed-start diagnostic executable")
    p.add_argument("--profile-exhaustion", action="store_true", help="build fixed-start diagnostic with rejected-line probes")
    p.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    p.add_argument("--jobs", default=2, type=int)
    args = p.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True)
    tree = root / "tree"
    tree.mkdir()
    export_src(REPO, BASE, tree)
    (tree / "tools").mkdir()
    builder = tree / "tools/build.py"
    builder.write_bytes(git(REPO, "show", BASE + ":tools/build.py"))
    patches = []
    if args.first_hunt_replay or args.profile_exhaustion:
        patches += [HERE / "first_hunt.patch", HERE / "quarter_control.patch"]
    if args.profile_exhaustion:
        patches.append(HERE / "exhaustion_profile.patch")
    for item in patches:
        subprocess.run(["git", "apply", "--whitespace=error", item], cwd=tree, check=True,
                       env={**os.environ, "GIT_CEILING_DIRECTORIES": str(root)})
    command = [sys.executable, str(builder), "--output-dir", str(root / "build"),
               "--jobs", str(args.jobs), "--cxx", args.cxx]
    includes = ["-I" + str(tree / "src/lib")]
    if args.boost_prefix:
        command += ["--boost-prefix", str(args.boost_prefix.resolve())]
        includes += ["-I" + str(args.boost_prefix.resolve() / "include")]
    subprocess.run(command, check=True)
    unit = root / "check"
    subprocess.run([args.cxx, "-std=c++14", "-O2", *includes, str(HERE / "check.cpp"),
                    *[str(tree / ("src/lib/" + name)) for name in ("quaternion.cpp", "random.cpp")],
                    "-o", str(unit)], check=True)
    subprocess.run([unit], check=True)
    write_wrappers(root)
    print("source", BASE)
    for path in (*patches, root / "build/vina"):
        print(hashlib.sha256(path.read_bytes()).hexdigest(), path.name)


if __name__ == "__main__":
    main()
