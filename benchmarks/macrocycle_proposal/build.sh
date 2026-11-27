#!/usr/bin/env bash
# Build the four recorded C++ experiments and both reference engines.
set -euo pipefail
root=$(git rev-parse --show-toplevel)
dest=$(realpath -m "$1")
shift
mkdir -p "$dest"
while read -r name rev; do
    git -C "$root" worktree add --quiet --detach "$dest/src-$name" "$rev"
    python3 "$dest/src-$name/tools/build.py" --output-dir "$dest/bin-$name" "$@"
done <<'ARMS'
parent 5000d1e00b348f5cb5593d0f322e6511004a86b3
ccd 6d2f76ff15eed4a6144c6fca8b1bbc6e4bb1a1b1
reset 9280bbbe26d02f87a3d69c6007e6994c3eff277e
lock fda379306a4f64ed008d046b2be5630c925877ad
null fed12226584e21899bc9edcc8d6c789a15baebdf
ARMS
python3 "$root/benchmarks/pair_gradient/build_arms.py" "$dest/up" --arms upstream upstream-counts "$@"
