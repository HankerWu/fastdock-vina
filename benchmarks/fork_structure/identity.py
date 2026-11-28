"""Bind measured executables to real Git revisions, declared defaults and runtime dependencies."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
START = "0449b1d0bc3adec768dd9a9f243f3ce5a29cef58"


def git(repo, ref):
    return subprocess.check_output(
        ["git", "-C", repo, "rev-parse", ref], text=True
    ).strip()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def identity(upstream, binaries, revision):
    prior = json.loads(
        (HERE.parent / "stage_local_structure/fixed/builds.json").read_text()
    )
    source = {
        "official": git(upstream, "v1.2.7"),
        "parent": git(ROOT, "5f3fc4b^"),
        "stable": START,
    }
    if "candidate" in binaries:
        source["candidate"] = git(ROOT, revision)
    result = {}
    for arm, binary in binaries.items():
        tree = git(upstream if arm == "official" else ROOT, source[arm] + ":src")
        digest = sha(binary.read_bytes())
        version = subprocess.check_output([binary, "--version"], text=True).strip()
        help_bytes, advanced = [
            subprocess.check_output([binary, opt])
            for opt in ("--help", "--help_advanced")
        ]
        defaults = dict(
            re.findall(
                r"^\s*--(\w+) arg \(=([^)]*)\)", (help_bytes + advanced).decode(), re.M
            )
        )
        dependencies = {
            name: sha(Path(file).read_bytes())
            for name, file in re.findall(
                r"^\s*(libboost_\S+)\s+=> (\S+) ",
                subprocess.check_output(["ldd", binary], text=True),
                re.M,
            )
        }
        if arm != "candidate":
            expected = prior["versions"][arm]
            assert tree == expected["src_tree"] and digest == expected["binary_sha256"]
            assert (
                version == expected["printed_version"]
                and sha(help_bytes) == expected["help_sha256"]
            )
            assert arm == "stable" or source[arm] == expected["revision"]
        assert (
            dependencies == prior["versions"]["stable"]["runtime_boost_sha256"]
            or dependencies == prior["versions"]["official"]["runtime_boost_sha256"]
        )
        result[arm] = dict(
            revision=source[arm],
            src_tree=tree,
            binary_sha256=digest,
            bytes=binary.stat().st_size,
            version=version,
            help_sha256=sha(help_bytes),
            help_advanced_sha256=sha(advanced),
            option_defaults=defaults,
            runtime_boost_sha256=dependencies,
        )
    compiler = Path(shutil.which("c++")).resolve()
    assert sha(compiler.read_bytes()) == prior["compiler_sha256"]
    return dict(
        versions=result,
        upstream_build_receipts="benchmarks/stage_local_structure/restored/upstream_rebuild.json",
        stable_build_receipt="benchmarks/fork_structure/stable_build.json",
        compiler=subprocess.check_output(
            [compiler, "--version"], text=True
        ).splitlines()[0],
        compiler_sha256=prior["compiler_sha256"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream", type=Path)
    parser.add_argument("output", type=Path)
    for arm in ("official", "parent", "stable", "candidate"):
        parser.add_argument("--" + arm, type=Path, required=arm != "candidate")
    parser.add_argument("--revision")
    a = parser.parse_args()
    binaries = {
        arm: getattr(a, arm)
        for arm in ("official", "parent", "stable", "candidate")
        if getattr(a, arm)
    }
    if bool(a.candidate) != bool(a.revision):
        parser.error("supply both a candidate executable and its revision")
    a.output.write_text(
        json.dumps(identity(a.upstream, binaries, a.revision), indent=2) + "\n"
    )
