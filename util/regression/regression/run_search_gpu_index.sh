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

"${RIBOSEEK}" gpuserver "${TARGETDB}" --db-load-mode 0 > "${RESULTS}/gpuserver.log" 2>&1 &
GPUSERVER_PID=$!
trap 'kill "${GPUSERVER_PID}" 2>/dev/null' EXIT INT TERM

"${RIBOSEEK}" search "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/tmp" \
    --prefilter-mode 1 --gpu 1 --gpu-server 1 --db-load-mode 2

"${RIBOSEEK}" convertalis "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/results.m8"

"${EVALUATE}" "${RESULTS}/results.m8" "${RESULTS}/roc1_auc.tsv" | tee "${RESULTS}/evaluation.log"
ACTUAL=$(grep "^ROC1-AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f2 | cut -d"," -f1)

EXPECTED="0.473903"
awk -v actual="$ACTUAL" -v expected="$EXPECTED" \
    'BEGIN { print (actual == expected) ? "GOOD" : "BAD"; print "Expected: ", expected; print "Actual: ", actual; }' \
    > "${RESULTS}.report"
