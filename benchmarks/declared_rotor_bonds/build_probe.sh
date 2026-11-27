#!/usr/bin/env bash
# Link a diagnostic against the same engine objects produced by tools/build.py.
# BOOST_PREFIX is optional when Boost is installed in the compiler search paths.
set -eu
tree=$1
build=$2
source=$3
output=$4
flags=(-O2 -DNDEBUG -std=c++14 -pthread "-I$tree/src/lib")
libs=(-lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options)
if [[ -n ${BOOST_PREFIX:-} ]]; then
    flags+=("-I$BOOST_PREFIX/include" "-L$BOOST_PREFIX/lib" "-Wl,-rpath,$BOOST_PREFIX/lib")
fi
objs=()
for obj in "$build"/obj/*.o; do
    [[ $obj == */main.o || $obj == */split.o ]] || objs+=("$obj")
done
"${CXX:-c++}" "${flags[@]}" "$source" "${objs[@]}" "${libs[@]}" -o "$output"
