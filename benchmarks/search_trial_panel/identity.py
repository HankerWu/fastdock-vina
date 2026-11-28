"""Bind tested binaries and their source identity to immutable source trees."""

import argparse
import hashlib
import json
import platform
import re
import subprocess
from pathlib import Path

from panel import HERE, UP
from run_panel import digest


def output(args):
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT).strip()


def git(repo, expression):
    return output(["git", "-C", repo, "rev-parse", expression])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("upstream", type=Path, help="unmodified official git repository")
    ap.add_argument(
        "candidate_revision", help="public commit containing the trial source"
    )
    ap.add_argument("output", type=Path)
    for role in ("official", "parent", "stable", "candidate"):
        ap.add_argument("--" + role, type=Path, required=True)
    a = ap.parse_args()
    repo = HERE.parents[1]
    archived = json.loads((UP / "builds.json").read_text())
    roles = {}
    revisions = {
        "official": archived["versions"]["official"]["revision"],
        "parent": archived["versions"]["parent"]["revision"],
        "stable": "ed0d42d50dd3f3de011b99b57ba76438929aae07",
        "candidate": a.candidate_revision,
    }
    for name, revision in revisions.items():
        source = a.upstream if name == "official" else repo
        binary = getattr(a, name)
        full = git(source, revision)
        tree = git(source, full + ":src")
        dependencies = {}
        for lib, file in re.findall(
            r"^\s*(libboost_\S+)\s+=> (\S+) ", output(["ldd", binary]), re.M
        ):
            dependencies[lib] = digest(Path(file))
        h = digest(binary)
        expected = (
            archived["versions"]["candidate" if name == "stable" else name]
            if name != "candidate"
            else None
        )
        if expected and (
            expected["src_tree"] != tree or expected["binary_sha256"] != h
        ):
            raise ValueError(f"archived identity mismatch: {name}")
        if not dependencies or any(
            archived["boost_soname_sha256"].get(lib) != v
            for lib, v in dependencies.items()
        ):
            raise ValueError(f"runtime Boost mismatch: {name}")
        help_bytes = subprocess.check_output([binary, "--help"])
        defaults = dict(
            re.findall(r"^\s*--(\w+) arg \(=([^)]*)\)", help_bytes.decode(), re.M)
        )
        roles[name] = {
            "revision": full,
            "src_tree": tree,
            "binary_sha256": h,
            "printed_version": output([binary, "--version"]),
            "help_sha256": hashlib.sha256(help_bytes).hexdigest(),
            "printed_option_defaults": defaults,
            "runtime_boost_sha256": dependencies,
        }
    if (
        git(repo, "5f3fc4b^") != roles["parent"]["revision"]
        or git(a.upstream, "v1.2.7") != roles["official"]["revision"]
    ):
        raise ValueError("source baseline identity mismatch")
    cxx = Path(output(["which", "c++"])).resolve()
    a.output.write_text(
        json.dumps(
            {
                "system": platform.system(),
                "arch": platform.machine(),
                "compiler": output([cxx, "--version"]).splitlines()[0],
                "compiler_sha256": digest(cxx),
                "versions": roles,
            },
            indent=2,
        )
        + "\n"
    )
