#!/bin/sh -e
zstd -qf -d -q "${DATADIR}/rfam_lookup_target_clan.tsv.zst" -o "${RESULTS}/rfam_lookup_target_clan.tsv"
zstd -f -d -q "${DATADIR}/target.tsv.zst" -o "${RESULTS}/target.tsv"
export DATADIR="${RESULTS}"

QUERY="${EXAMPLEDIR}/QUERY.fasta"
QUERYDB="${RESULTS}/query"
"${RIBOSEEK}" createdb "${QUERY}" "${QUERYDB}"

TARGET="${EXAMPLEDIR}/DB.fasta"
TARGETDB="${RESULTS}/target"
"${RIBOSEEK}" createdb "${TARGET}" "${RESULTS}/target_tmp"
"${RIBOSEEK}" makepaddedseqdb "${RESULTS}/target_tmp" "${TARGETDB}"
"${RIBOSEEK}" createindex "${TARGETDB}" "${RESULTS}/index_tmp" --index-subset 10

"${RIBOSEEK}" gpuserver "${TARGETDB}" --db-load-mode 0 --gpu-server-namespace first \
    > "${RESULTS}/gpuserver-first.log" 2>&1 &
GPUSERVER_FIRST_PID=$!
"${RIBOSEEK}" gpuserver "${TARGETDB}" --db-load-mode 0 --gpu-server-namespace second \
    > "${RESULTS}/gpuserver-second.log" 2>&1 &
GPUSERVER_SECOND_PID=$!
trap 'kill "${GPUSERVER_FIRST_PID}" "${GPUSERVER_SECOND_PID}" 2>/dev/null' EXIT INT TERM

wait_for_server() {
    LOG=$1
    PID=$2
    ATTEMPT=0
    while [ "$ATTEMPT" -lt 120 ] && ! grep -Eq '^[0-9]+$' "$LOG"; do
        kill -0 "$PID"
        sleep 0.5
        ATTEMPT=$((ATTEMPT + 1))
    done
    grep -Eq '^[0-9]+$' "$LOG"
}
wait_for_server "${RESULTS}/gpuserver-first.log" "${GPUSERVER_FIRST_PID}"
wait_for_server "${RESULTS}/gpuserver-second.log" "${GPUSERVER_SECOND_PID}"

# A duplicate owner must fail instead of reopening and corrupting the live
# server's shared-memory segment.
if "${RIBOSEEK}" gpuserver "${TARGETDB}" --db-load-mode 0 --gpu-server-namespace first \
    > "${RESULTS}/gpuserver-duplicate.log" 2>&1; then
    echo "Duplicate gpuserver unexpectedly succeeded" >&2
    exit 1
fi
grep -q "already exists" "${RESULTS}/gpuserver-duplicate.log"

"${RIBOSEEK}" search "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/tmp" \
    --prefilter-mode 1 --gpu 1 --gpu-server 1 --gpu-server-namespace first --db-load-mode 2

"${RIBOSEEK}" search "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results-second" "${RESULTS}/tmp-second" \
    --prefilter-mode 1 --gpu 1 --gpu-server 1 --gpu-server-namespace second --db-load-mode 2

"${RIBOSEEK}" convertalis "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/results.m8"

"${EVALUATE}" "${RESULTS}/results.m8" "${RESULTS}/roc1_auc.tsv" | tee "${RESULTS}/evaluation.log"
ACTUAL=$(grep "^ROC1-AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f2 | cut -d"," -f1)

EXPECTED="0.473903"
awk -v actual="$ACTUAL" -v expected="$EXPECTED" \
    'BEGIN { print (actual == expected) ? "GOOD" : "BAD"; print "Expected: ", expected; print "Actual: ", actual; }' \
    > "${RESULTS}.report"
