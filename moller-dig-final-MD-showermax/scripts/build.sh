#!/usr/bin/env bash
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${REMOLL_DIR:?Set REMOLL_DIR to your remoll installation containing include and lib64/libremoll.so}"
cmake -S "$project_dir" -B "$project_dir/build" -DREMOLL_DIR="$REMOLL_DIR"
cmake --build "$project_dir/build" -j "${BUILD_JOBS:-4}"
