"""Export pinned FastDock/Vina sources and independently compile each research arm."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PARENT = "ac0a2f6e40cac8f62748848227f298d44b7b764e"
UPSTREAM = "3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645"
EXACT_HEADER = "e325e67b7b593bdea34a6970361d31c3bfafbb38"
HEADER = Path("src/lib/precalculate.h")
ARMS = ("parent", "default", "linear", "hermite", "shifted", "upstream", "upstream-counts")


def git(root, *args):
    return subprocess.check_output(["git", "-C", root, *args])


def export_src(repo, revision, dest):
    archive = git(repo, "archive", "--format=tar", revision, "src")
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        tar.extractall(dest, filter="data")


def build_arm(arm, args, upstream):
    root = args.output / arm
    root.mkdir()  # never overwrite an existing source or build
    is_upstream = arm.startswith("upstream")
    revision = UPSTREAM if is_upstream else PARENT
    export_src(upstream if is_upstream else REPO, revision, root)
    (root / "tools").mkdir()
    build = root / "tools/build.py"
    build.write_bytes(git(REPO, "show", f"{PARENT}:tools/build.py"))
    if arm in ("default", "linear", "shifted"):
        (root / HEADER).write_bytes(git(REPO, "cat-file", "blob", EXACT_HEADER))
    patches = {"upstream-counts": "upstream_counts.patch", "hermite": "monotone.patch",
               "shifted": "shifted_explicit.patch"}
    if arm in patches:
        # Also allow an export under the calling checkout's build directory.
        # Git must treat root as a standalone directory when applying paths.
        patch_env = dict(os.environ, GIT_CEILING_DIRECTORIES=str(root.parent))
        subprocess.run(["git", "apply", "--whitespace=error", HERE / patches[arm]],
                       cwd=root, env=patch_env, check=True)
    compiler = args.cxx
    exact = arm in ("linear", "shifted")
    if exact:
        wrapper = root / "compiler"
        wrapper.write_text("#!/bin/sh\nexec " + shlex.quote(compiler) +
                           ' -DFASTDOCK_EXACT_PAIR_GRADIENT "$@"\n')
        wrapper.chmod(0o755)
        compiler = str(wrapper)
    command = [sys.executable, str(build), "--output-dir", str(root / "build"), "--jobs",
               str(args.jobs), "--cxx", compiler]
    if args.boost_prefix:
        command += ["--boost-prefix", str(args.boost_prefix)]
    with (root / "build.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    manifest = dict(revision=revision, exact_macro=exact, build_script_sha256=sha(build),
                    patch=patches.get(arm), patch_sha256=(sha(HERE / patches[arm]) if arm in patches else None),
                    compiler=subprocess.check_output([args.cxx, "--version"], text=True).splitlines()[0],
                    build=command, binary_sha256=sha(root / "build/vina"),
                    header_sha256=sha(root / HEADER))
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(arm, manifest["binary_sha256"], flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path, help="Exports are written to new ARM subdirectories")
    p.add_argument("--arms", choices=ARMS, nargs="+", default=ARMS)
    p.add_argument("--upstream", type=Path, help="Optional existing official AutoDock-Vina Git checkout")
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--cxx", default="c++")
    p.add_argument("--jobs", type=int, default=2)
    args = p.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.boost_prefix:
        args.boost_prefix = args.boost_prefix.resolve()
    if args.jobs < 1 or len(set(args.arms)) != len(args.arms):
        p.error("Require positive jobs and distinct arms")
    if any((args.output / arm).exists() for arm in args.arms):
        p.error("Use a fresh output directory for each requested arm")
    if not shutil.which(args.cxx):
        p.error("Compiler unavailable")
    upstream = args.upstream
    if any(arm.startswith("upstream") for arm in args.arms):
        if upstream is None:
            upstream = args.output / "AutoDock-Vina.git"
            if not upstream.exists():
                subprocess.run(["git", "init", "--bare", upstream], check=True)
            if subprocess.run(["git", "-C", upstream, "cat-file", "-e", UPSTREAM],
                              stderr=subprocess.DEVNULL).returncode:
                subprocess.run(["git", "-C", upstream, "fetch", "--depth=1",
                                "https://github.com/ccsb-scripps/AutoDock-Vina.git", UPSTREAM], check=True)
        if git(upstream, "rev-parse", UPSTREAM + "^{commit}").decode().strip() != UPSTREAM:
            p.error("Pinned upstream commit unavailable")
    for arm in args.arms:
        build_arm(arm, args, upstream)


if __name__ == "__main__":
    main()
