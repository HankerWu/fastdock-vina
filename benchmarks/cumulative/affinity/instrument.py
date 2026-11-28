"""Apply existing upstream passive BFGS telemetry to disposable exact-source copies."""

import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
PATCH = ROOT / "benchmarks/pair_gradient/upstream_counts.patch"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, help="fresh disposable official or parent source")
    args = p.parse_args()
    source = args.source.resolve()
    patch = PATCH.read_text()
    # The parent reuses a candidate buffer. Adjust only the context of the
    # existing instrumentation patch, with no candidate or RNG mutations.
    if "candidate.c = tmp.c;" in (source / "src/lib/monte_carlo.cpp").read_text():
        patch = patch.replace(" \t\toutput_type candidate = tmp;", " \t\tcandidate.c = tmp.c;")
        patch = patch.replace(
            " \t\t\ttmp = candidate;",
            " \t\t\tstd::swap(tmp.c, candidate.c);\n \t\t\ttmp.e = candidate.e;",
        )
    subprocess.run(["git", "apply", "--recount", "-"], input=patch, text=True, cwd=source, check=True)


if __name__ == "__main__":
    main()
