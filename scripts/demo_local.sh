#!/bin/bash
# SPDX-License-Identifier: MIT
# Local demo: a 3 x 3 grid of charging stations + base station for 20 s.
# The base station prints every congestion alert and its decision, then the summary.
#   make demo            (or: LAUNCHER="mpirun --oversubscribe" bash scripts/demo_local.sh)
set -u
LAUNCH=${LAUNCHER:-mpirun --oversubscribe}
OUT=out/demo_$(date +%Y%m%d_%H%M%S)
$LAUNCH -n 10 bin/wsn 3 3 -c 500 -t 20 -a 0.7 -d 0.2 -W 2 -o "$OUT"
echo
echo "---- base station log (first alerts) : $OUT/base_log.txt"
sed -n '1,14p' "$OUT/base_log.txt"
echo "...."
grep -A12 "SIMULATION SUMMARY" "$OUT/base_log.txt"
echo
echo "---- node 4 (centre) availability history : $OUT/node_004.log"
grep -A8 "Availability history" "$OUT/node_004.log"
