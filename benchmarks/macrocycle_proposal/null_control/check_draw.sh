#!/usr/bin/env bash
# Usage: check_draw.sh NULL_SOURCE NULL_BUILD BOOST_ROOT BACE_INPUT MER_INPUT
set -eu
here=$(cd "$(dirname "$0")" && pwd)
d=$(mktemp -d)
trap 'rm -rf "$d"' EXIT
objects=()
for f in "$2"/obj/*.o; do
    if [[ $(basename "$f") != main.o ]]; then objects+=("$f"); fi
done
c++ -O2 -DNDEBUG -std=c++14 -pthread -I"$1/src/lib" -I"$3/include" "$here/check_draw.cpp" \
    "${objects[@]}" -L"$3/lib" -Wl,-rpath,"$3/lib" \
    -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -o "$d/check"
printf 'pdb,all_entities,allowed,locked,same_nonnull_draws,same_null_draws\n'
printf '5QCT,'
"$d/check" "$4"
printf '5K0X,'
"$d/check" "$5"
