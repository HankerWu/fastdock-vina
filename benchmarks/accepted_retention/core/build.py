"""Build a pinned parent or one-slot research arm in a separate directory."""

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "pair_gradient"))
from build_arms import export_src, git

PARENT = "31d073fdfb4baf7859b90a926f5debc71de26528"
REPO = HERE.parents[2]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arm", choices=("parent", "candidate", "ungated", "recurrence"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--boost-prefix", type=Path)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True)
    export_src(REPO, PARENT, root)
    (root / "tools").mkdir()
    builder = root / "tools/build.py"
    builder.write_bytes(git(REPO, "show", f"{PARENT}:tools/build.py"))
    patches = [] if args.arm == "parent" else [HERE / "retention.patch"]
    if args.arm == "recurrence":
        patches.append(HERE.parent / "persistence/recurrence.patch")
    for patch in patches:
        subprocess.run(
            ["git", "apply", "--whitespace=error", patch],
            cwd=root,
            check=True,
            env={**os.environ, "GIT_CEILING_DIRECTORIES": str(root.parent)},
        )
    compiler = args.cxx
    if args.arm == "ungated":
        wrapper = root / "research_cxx"
        wrapper.write_text(
            "#!/bin/sh\nexec "
            + shlex.quote(compiler)
            + ' -DFASTDOCK_UNGATED_DEFERRED "$@"\n'
        )
        wrapper.chmod(0o755)
        compiler = str(wrapper)
    command = [
        sys.executable,
        builder,
        "--output-dir",
        root / "build",
        "--cxx",
        compiler,
        "--jobs",
        str(args.jobs),
    ]
    if args.boost_prefix:
        command += ["--boost-prefix", args.boost_prefix.resolve()]
    subprocess.run(command, check=True)
