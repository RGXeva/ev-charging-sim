#!/bin/bash
# SPDX-License-Identifier: MIT
# Run every line of a case file; results go to results/raw/<experiment>/<label>_r<rep>/.
# usage:   bash scripts/run_cases.sh jobs/cases/A_grid_1.txt
# env:     LAUNCHER  (default "srun"; locally use "mpirun --oversubscribe")
#          BUDGET    seconds available for runs (default 540, CAAS jobs are 10 min)
#          RAW       output root (default results/raw)
# A run whose summary.csv already exists is skipped, so a case file can simply be
# resubmitted after a job ran out of time.
set -u
CASES=$1
LAUNCH=${LAUNCHER:-srun}
BUDGET=${BUDGET:-540}
RAW=${RAW:-results/raw}
SECONDS=0

echo "ENV job=${SLURM_JOB_ID:-local} nodelist=${SLURM_JOB_NODELIST:-$(hostname -s)} host=$(hostname -s) start=$(date '+%F %T') cases=$CASES launcher='$LAUNCH'"
(md5sum bin/wsn 2>/dev/null || md5 -r bin/wsn) | sed 's/^/MD5 /'

n_ok=0; n_fail=0; n_skip=0
while read -r exp label rep m n args; do
    [[ -z "${exp:-}" || "$exp" == \#* ]] && continue
    np=$((m * n + 1))
    out="$RAW/$exp/${label}_r$rep"
    if [[ -f "$out/summary.csv" ]]; then
        echo "DONE $exp $label r$rep (already)"; continue
    fi
    t=$(sed -n 's/.*-t \([0-9.]*\).*/\1/p' <<< "$args"); t=${t%.*}
    est=$(( ${t:-10} + 6 ))
    if (( SECONDS + est > BUDGET )); then
        echo "SKIP $exp $label r$rep (time budget: ${SECONDS}s used)"; n_skip=$((n_skip + 1)); continue
    fi
    rm -rf "$out"; mkdir -p "$out"
    extra=""
    if [[ "$LAUNCH" == srun* ]]; then
        extra="-N $(( (np + 15) / 16 )) --distribution=block"   # fill one 16-core node first
    fi
    t0=$SECONDS
    # shellcheck disable=SC2086
    $LAUNCH -n $np $extra bin/wsn $m $n $args -s $((1000 + rep)) -o "$out" < /dev/null > "$out/stdout.txt" 2>&1
    rc=$?
    if [[ $rc -eq 0 && -f "$out/summary.csv" ]]; then
        n_ok=$((n_ok + 1)); status=OK
    else
        n_fail=$((n_fail + 1)); status=FAIL
    fi
    echo "$status $exp $label r$rep np=$np rc=$rc ${SECONDS}s(+$((SECONDS - t0)))"
done < "$CASES"
echo "END ok=$n_ok fail=$n_fail skip=$n_skip elapsed=${SECONDS}s $(date '+%F %T')"
