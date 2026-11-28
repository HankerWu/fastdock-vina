"""Verify the upstream commit, root tree, archived source tree, license and three local revisions."""

import base64
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
ROOT = BENCH.parent


def git(repo, *args, data=None):
    return (
        subprocess.check_output(["git", "-C", repo, *args], input=data).decode().strip()
    )


def main():
    upstream = json.loads((HERE / "upstream.json").read_text())
    builds = json.loads(
        (BENCH / "stage_local_structure/fixed/builds.json").read_text()
    )["versions"]
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(HERE / "vina-1.2.7-src.tar.xz") as archive:
            assert archive.pax_headers["comment"] == upstream["revision"]
            archive.extractall(tmp, filter="data")
        git(tmp, "init", "-q")
        root = git(
            tmp,
            "hash-object",
            "-w",
            "-t",
            "tree",
            "--stdin",
            data=base64.b64decode(upstream["root_tree_base64"]),
        )
        revision = git(
            tmp,
            "hash-object",
            "-w",
            "-t",
            "commit",
            "--stdin",
            data=upstream["commit"].encode(),
        )
        assert revision == upstream["revision"] == builds["official"]["revision"]
        assert git(tmp, "rev-parse", revision + "^{tree}") == root
        git(tmp, "add", "src", "LICENSE")
        assert git(tmp, "rev-parse", ":LICENSE") == git(
            tmp, "rev-parse", root + ":LICENSE"
        )
        actual = git(tmp, "write-tree", "--prefix=src/")
        assert (
            actual
            == git(tmp, "rev-parse", root + ":src")
            == builds["official"]["src_tree"]
        )
    for role in ("parent", "stable", "candidate"):
        item = builds[role]
        assert git(ROOT, "rev-parse", item["revision"] + ":src") == item["src_tree"]
    assert git(ROOT, "rev-parse", "5f3fc4b^") == builds["parent"]["revision"]
    print(
        "Verified four immutable source trees, upstream original commit and license, and true pre-Pro parent"
    )


if __name__ == "__main__":
    main()
