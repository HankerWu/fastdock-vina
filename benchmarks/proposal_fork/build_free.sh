#!/usr/bin/env bash
set -eu
root=$(git rev-parse --show-toplevel)
src=$(realpath -m "$1")
boost=$(realpath "$2")
mkdir "$src"
git -C "$root" archive 1c672a8e2ce12c7a16e20172a3315ea5a19c23d8 src tools/build.py | tar -xC "$src"
python "$root/benchmarks/bfgs_source_replay/instrument.py" "$src"
python "$src/tools/build.py" --boost-prefix "$boost" --output-dir "$src/build" --jobs "${JOBS:-3}"
