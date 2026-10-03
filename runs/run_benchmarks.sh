#!/usr/bin/env bash
# Run the benchmark parameter files of models/parameter-files/attractor{E,T}_*.in, one after another.
# Usage (from anywhere):  runs/run_benchmarks.sh [name ...] [-- key=value ...]
#   names are file stems without the "attractor" prefix, e.g. E_k6_a1 T_k10_a1 (default: all 14)
#   keys after "--" override the .in file for every run, e.g. -- tMax=300 withGWs=true
# Each run goes to runs/benchmarks/<name>/ (run.in = copy of the benchmark file + overrides).
# Runs that already finished (file "done" present) are skipped. Never run two batches at once.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PF="$ROOT/models/parameter-files"

names=(); overrides=()
while [[ $# -gt 0 ]]; do
    if [[ "$1" == "--" ]]; then shift; overrides=("$@"); break; fi
    names+=("$1"); shift
done
if [[ ${#names[@]} -eq 0 ]]; then
    for f in "$PF"/attractor[ET]_*.in; do b="$(basename "$f" .in)"; names+=("${b#attractor}"); done
fi

# threads = physical cores is fastest
THREADS="${THREADS:-$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)}"

for name in "${names[@]}"; do
    model="attractor${name:0:1}"
    exe="$ROOT/build_$model/$model"
    out="$ROOT/runs/benchmarks/$name"
    [[ -x "$exe" ]] || { echo "missing $exe: build it first (cmake .. -DMODEL=$model)"; exit 1; }
    if [[ -f "$out/done" ]]; then echo "skip $name (done)"; continue; fi
    mkdir -p "$out"
    cp "$PF/attractor$name.in" "$out/run.in"
    for kv in "${overrides[@]}"; do
        key="${kv%%=*}"
        if grep -qE "^$key[[:space:]]*=" "$out/run.in"; then
            sed -i -E "s|^$key[[:space:]]*=.*|$kv # override|" "$out/run.in"
        else
            echo "$kv # override" >> "$out/run.in"
        fi
    done
    echo "$(date +%T) start $name"
    start=$SECONDS
    (cd "$out" && OMP_NUM_THREADS="$THREADS" OMP_PROC_BIND=false "$exe" input=run.in overwriteFiles=true > log.txt 2>&1)
    echo "$((SECONDS - start)) s" > "$out/done"
    echo "$(date +%T) done  $name ($((SECONDS - start)) s)"
done
