#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 3 ]]; then echo 'usage: build.sh VINA_OBJECT_DIR BOOST_PREFIX OUTPUT' >&2; exit 2; fi
here=$(cd "$(dirname "$0")" && pwd)
root=$(git -C "$here" rev-parse --show-toplevel)
objs=()
for obj in "$1"/*.o; do
    case "$(basename "$obj")" in main.o|split.o) continue ;; esac
    objs+=("$obj")
done
"${CXX:-c++}" -std=c++14 -O3 -DNDEBUG -pthread -fvisibility=hidden -I"$root/src/lib" -I"$2/include" "$here/fixed.cpp" "${objs[@]}" -L"$2/lib" -Wl,-rpath,"$2/lib" -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -o "$3"
