#!/bin/bash
# run.sh E02_firm_dynamics E04_size_gradient ...: runs the named extension files side by side, one Julia process each,
# with logs in logs/<name>.txt; prints each exit code when all are done. Run from the FablePreSolve folder (Git Bash).
cd "$(dirname "$0")"
JL="$HOME/.julia/juliaup/julia-1.12.7+0.x64.w64.mingw32/bin/julia.exe"
export OPENBLAS_NUM_THREADS=1
mkdir -p logs
pids=()
for f in "$@"; do
    "$JL" --project=../../.. "$f.jl" > "logs/$f.txt" 2>&1 &
    pids+=($!)
done
i=0
for f in "$@"; do
    wait "${pids[$i]}"; echo "$f: exit $?"
    i=$((i+1))
done
