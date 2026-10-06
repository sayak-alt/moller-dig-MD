#!/usr/bin/env bash
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
if (( $# < 2 || $# > 3 )); then
  echo 'Usage: scripts/run.sh CONFIG.dat FILELIST.txt [Nevents]' >&2
  exit 1
fi
exec "$project_dir/build/mollerdig" "$@"
