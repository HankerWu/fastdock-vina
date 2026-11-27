#!/bin/sh
set -eu
# Run from repository root; use the objects of the same source tree.
objects=$1
output=$2
shift 2
set -- "$@" -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options
for o in "$objects"/*.o; do
    case "${o##*/}" in main.o|split.o) continue ;; esac
    set -- "$o" "$@"
done
c++ -std=c++14 -DNDEBUG -pthread -Isrc/lib benchmarks/cavity_field_rank/audit/types.cpp "$@" -o "$output"
