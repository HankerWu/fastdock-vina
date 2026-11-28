"""Verify a new stable rebuild against the byte-reproducible historical executable."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
START = "0449b1d0bc3adec768dd9a9f243f3ce5a29cef58"


def build(work, baseline, boost, output):
    out = work.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = out / "src-tree"
    source.mkdir()
    payload = subprocess.check_output(
        ["git", "-C", ROOT, "archive", START, "src", "tools/build.py"]
    )
    subprocess.run(["tar", "-x", "-C", source], input=payload, check=True)
    historical = json.loads(
        (
            HERE.parent / "stage_local_structure/restored/path_reproducibility.json"
        ).read_text()
    )
    cxx = out / "cxx"
    import shlex

    remap = (
        "-ffile-prefix-map=" + str(source) + "=" + historical["old_source_directory"]
    )
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
    assert digest == historical["archived_sha256"] and binary == baseline.read_bytes()
    result = dict(
        revision=START,
        src_tree=subprocess.check_output(
            ["git", "-C", ROOT, "rev-parse", START + ":src"], text=True
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
    build(**vars(parser.parse_args()))
