#!/bin/bash
# SPDX-License-Identifier: MIT
# Live demo on a SLURM cluster: a 3 x 3 grid + base station runs on compute nodes for
# 20 s; the base station prints every congestion alert and its decision, then the summary.
#   bash scripts/demo_cluster.sh
# Load your cluster's Open MPI module first (it must wrap a compiler with OpenMP).
set -u
echo "user=$(whoami) host=$(hostname -s) date=$(date '+%F %T') dir=$(pwd)"
md5sum bin/wsn
OUT=results/caas/demo_$(date +%Y%m%d_%H%M%S)
srun -p defq --time=00:03:00 --mem=2G -N 1 -n 10 \
    bin/wsn 3 3 -c 500 -t 20 -a 0.7 -d 0.2 -W 2 -o "$OUT"
echo
echo "---- base station log (first alerts) : $OUT/base_log.txt"
sed -n '1,14p' "$OUT/base_log.txt"
echo "...."
grep -A12 "SIMULATION SUMMARY" "$OUT/base_log.txt"
echo
echo "---- node 4 (centre) availability history : $OUT/node_004.log"
grep -A8 "Availability history" "$OUT/node_004.log"
