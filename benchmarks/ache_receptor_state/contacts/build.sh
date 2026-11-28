#!/usr/bin/env bash
# Link the typing-only probe against the same objects as the ordinary Pro.
set -eu
build=$1
boost=$2
output=$3
objects=()
for obj in "$build"/obj/*.o; do
    case ${obj##*/} in main.o|split.o) continue ;; esac
    objects+=("$obj")
done
c++ -std=c++14 -O3 -DNDEBUG -pthread -Isrc/lib -I"$boost/include" \
    benchmarks/ache_receptor_state/contacts/topology.cpp "${objects[@]}" \
    -L"$boost/lib" -Wl,-rpath,"$boost/lib" -lboost_thread -lboost_filesystem \
    -lboost_serialization -lboost_program_options -o "$output"
