"""Verify four real executable/source identities, including the exact round-start source tree."""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
START = "666f49d2570c047470d3c281df30d18c5ebd3ce1"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("upstream", type=Path)
parser.add_argument("candidate_revision")
parser.add_argument("output", type=Path)
roles = ("official", "parent", "stable", "candidate")
for role in roles:
    parser.add_argument("--" + role, required=True, type=Path)
a = parser.parse_args()
legacy = HERE.parent / "search_trial_panel/identity.py"
subprocess.run([sys.executable, legacy, *sys.argv[1:]], check=True)
manifest = json.loads(a.output.read_text())
stable = manifest["versions"]["stable"]
tree = subprocess.check_output(
    ["git", "-C", ROOT, "rev-parse", START + ":src"], text=True
).strip()
assert tree == stable["src_tree"]
stable["same_source_archived_build_revision"] = stable["revision"]
stable["revision"] = START
for role in roles:
    binary = getattr(a, role)
    assert (
        hashlib.sha256(binary.read_bytes()).hexdigest()
        == manifest["versions"][role]["binary_sha256"]
    )
manifest["build_recipe"] = (
    "benchmarks/geometry_comparison/upstream/four_versions/build.sh for controls; archive CANDIDATE src tools/build.py, run build.py --boost-prefix BOOST --output-dir OUT --jobs 2"
)
a.output.write_text(json.dumps(manifest, indent=2) + "\n")
