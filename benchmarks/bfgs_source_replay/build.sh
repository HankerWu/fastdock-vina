#!/usr/bin/env bash
# Build four diagnostic-only executables from pinned trees and their build scripts.
set -eu
if [ "$#" -ne 2 ]; then
    echo 'usage: bash build.sh FRESH_BUILD_DIRECTORY BOOST_PREFIX' >&2
    exit 2
fi
base=$(realpath -m "$1")
boost=$(realpath "$2")
protocol=$(cd "$(dirname "$0")" && pwd)
root=$(git -C "$protocol" rev-parse --show-toplevel)
mkdir "$base"
for rev in e7031e9aabc0b31fd371370763bc3810748a4cfa 232643049094dee7fece70498f77f5271d621949 ed0d42d50dd3f3de011b99b57ba76438929aae07 0bb827ebfe69a32504b74488730f80a21495eccf; do
    arm=${rev:0:7}
    source="$base/$arm-source"
    mkdir "$source"
    git -C "$root" archive "$rev" src tools/build.py | tar -xC "$source"
    case "$arm" in
        e7031e9|2326430) python "$protocol/instrument.py" "$source" ;;
        *) python "$protocol/instrument_candidate.py" "$source" ;;
    esac
    python "$source/tools/build.py" --boost-prefix "$boost" --jobs "${JOBS:-2}" --output-dir "$base/$arm"
done
