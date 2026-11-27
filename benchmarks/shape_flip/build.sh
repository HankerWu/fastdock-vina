#!/usr/bin/env bash
# Build the two research revisions outside the production checkout.
set -euo pipefail
if (( $# < 1 )); then echo "usage: build.sh NEW_WORKDIR [tools/build.py options]" >&2; exit 2; fi
work=$(realpath -m "$1")
shift
repo=$(git rev-parse --show-toplevel)
rev=${SHAPE_SOURCE:-b8f9aad6bd5d757c83d6bd6a4e0ad2ed4fa7bf95}
compiler=$(command -v "${CXX:-c++}")
mkdir -p "$work"
mkdir "$work/shape-source"
git -C "$repo" archive "$rev" | tar -x -C "$work/shape-source"
cat > "$work/always-cxx" <<'WRAPPER'
#!/usr/bin/env bash
exec "$SHAPE_CXX" -DSHAPE_FLIP_ALWAYS=1 "$@"
WRAPPER
chmod +x "$work/always-cxx"
export SHAPE_CXX="$compiler"
python3 "$repo/tools/build.py" --output-dir "$work/parent" "$@"
for variant in adaptive always observer-adaptive observer-always; do
    if [[ "$variant" == observer-adaptive ]]; then
        (cd "$work/shape-source" && git apply --unidiff-zero "$repo/benchmarks/shape_flip/observer.patch")
    fi
    cxx="$compiler"
    if [[ "$variant" == *always ]]; then cxx="$work/always-cxx"; fi
    python3 "$work/shape-source/tools/build.py" --output-dir "$work/$variant" "$@" --cxx "$cxx"
done
