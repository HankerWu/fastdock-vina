#!/bin/bash
set -eu
D=$(cd "$(dirname "$0")" && pwd)
R=$(git -C "$D" rev-parse --show-toplevel)
S=$(realpath -m "$1")
B=$(realpath "$2")
mkdir "$S"
git -C "$R" archive 0bb827ebfe69a32504b74488730f80a21495eccf src tools/build.py | tar -xC "$S"
python "$D/instrument.py" "$S"
python "$S/tools/build.py" --boost-prefix "$B" --output-dir "$S/build" --jobs "${JOBS:-2}"
