#!/usr/bin/env bash
# Build immutable source versions with one toolchain. Run from the repository root.
set -eu
if [ "$#" -ne 2 ]; then
    echo 'usage: bash build.sh NEW_BUILD_DIRECTORY BOOST_PREFIX' >&2
    exit 2
fi
work=$(realpath -m "$1")
boost=$(realpath "$2")
mkdir "$work"
root=$(git rev-parse --show-toplevel)
official=8eb40404f4f45608acb3b01427587ac049f27c1f
parent=05ec6ed580999d9c541664673c6dfe46b89be0f5
stable=f46aa3774f9c77503b98654ec3903fbdb1c72fc0
candidate=${CANDIDATE_COMMIT:-8d24ebebbd272d7ef05f15cd5ec9643460192f66}
git clone --bare --depth 1 --branch v1.2.7 https://github.com/forlilab/AutoDock-Vina "$work/official.git"
test "$(git --git-dir="$work/official.git" rev-parse v1.2.7)" = "$official"
for arm in official parent stable candidate; do
    src="$work/$arm-source"
    mkdir -p "$src/tools"
    if [ "$arm" = official ]; then
        git --git-dir="$work/official.git" archive "$official" src | tar -x -C "$src"
        git -C "$root" show "$stable:tools/build.py" | sed 's/v1.2.7-27-g3c65c0b/v1.2.7/' > "$src/tools/build.py"
    else
        case "$arm" in parent) commit=$parent ;; stable) commit=$stable ;; candidate) commit=$candidate ;; esac
        git -C "$root" archive "$commit" src tools/build.py | tar -x -C "$src"
    fi
    python "$src/tools/build.py" --boost-prefix "$boost" --output-dir "$work/$arm" --jobs "${JOBS:-2}"
done
