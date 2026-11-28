"""Rebuild a registered Pro revision and compare the entire measured executable."""

import argparse
import hashlib
import json
import os
import re
import shlex
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
START = "0449b1d0bc3adec768dd9a9f243f3ce5a29cef58"


def build(work, baseline, boost, output, round_candidate=False):
    out = work.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = out / "src-tree"
    source.mkdir()
    role = "candidate" if round_candidate else "stable"
    registered = json.loads((HERE / "builds.json").read_text())["versions"][role]
    revision = registered["revision"]
    if not round_candidate:
        assert revision == START
    payload = subprocess.check_output(
        ["git", "-C", ROOT, "archive", revision, "src", "tools/build.py"]
    )
    subprocess.run(["tar", "-x", "-C", source], input=payload, check=True)
    historical = json.loads(
        (
            HERE.parent / "stage_local_structure/restored/path_reproducibility.json"
        ).read_text()
    )
    cxx = out / "cxx"
    original = historical["old_source_directory"]
    if round_candidate:
        command_line = (HERE / "candidate_build.log").read_text().splitlines()[0]
        original = re.search(r"-I([^ ]+)/src/lib ", command_line)[1]
    remap = "-ffile-prefix-map=" + str(source) + "=" + original
    cxx.write_text("#!/bin/sh\nexec /usr/bin/c++ " + shlex.quote(remap) + ' "$@"\n')
    cxx.chmod(0o755)
    command = [
        sys.executable,
        str(source / "tools/build.py"),
        "--boost-prefix",
        str(boost.resolve()),
        "--output-dir",
        str(out / "bin"),
        "--cxx",
        str(cxx),
        "--jobs",
        "2",
    ]
    start = time.monotonic()
    with (out / "build.log").open("w") as stream:
        proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    if proc.returncode:
        raise RuntimeError("build failed: " + str(out / "build.log"))
    binary = (out / "bin/vina").read_bytes()
    digest = hashlib.sha256(binary).hexdigest()
    assert digest == registered["binary_sha256"] and binary == baseline.read_bytes()
    result = dict(
        revision=revision,
        src_tree=subprocess.check_output(
            ["git", "-C", ROOT, "rev-parse", revision + ":src"], text=True
        ).strip(),
        source_prefix_map=remap,
        command=command,
        cpu_s=usage.ru_utime + usage.ru_stime,
        wall_s=time.monotonic() - start,
        peak_rss_kib=usage.ru_maxrss,
        binary_sha256=digest,
        bytes=len(binary),
        byte_identical=True,
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    output.with_suffix(".log").write_bytes((out / "build.log").read_bytes())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("work", "baseline", "boost", "output"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--round-candidate", action="store_true")
    build(**vars(parser.parse_args()))
