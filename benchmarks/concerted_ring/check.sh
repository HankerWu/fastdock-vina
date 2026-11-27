#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$1 build=$2 inputs=$3 dest=$4
shift 4
flags=(-O3 -DNDEBUG -std=c++14 -pthread -I"$src/src/lib")
libs=(-lboost_thread -lboost_filesystem -lboost_program_options -lboost_serialization)
if [[ -n ${BOOST_PREFIX:-} ]]; then
    flags+=(-I"$BOOST_PREFIX/include" -L"$BOOST_PREFIX/lib" -Wl,-rpath,"$BOOST_PREFIX/lib")
fi
objects=()
for obj in "$build"/obj/*.o; do
    case "${obj##*/}" in main.o|split.o) continue;; esac
    objects+=("$obj")
done
c++ "${flags[@]}" "$here/check.cpp" "${objects[@]}" "${libs[@]}" -o "$dest"
"$dest" "$inputs/5T9W/ligand.pdbqt" "$inputs/7MBO/ligand.pdbqt" "$inputs/5K0X/ligand.pdbqt" "$@"
