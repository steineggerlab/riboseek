#!/bin/sh -e
zstd -qf -d -q "${DATADIR}/rfam_lookup_target_clan.tsv.zst" -o "${RESULTS}/rfam_lookup_target_clan.tsv"
zstd -f -d -q "${DATADIR}/target.tsv.zst" -o "${RESULTS}/target.tsv"
export DATADIR="${RESULTS}"

QUERY="${EXAMPLEDIR}/QUERY.fasta"
QUERYDB="${RESULTS}/query"
"${RIBOSEEK}" createdb "${QUERY}" "${QUERYDB}"

TARGET="${EXAMPLEDIR}/DB.fasta"
TARGETDB="${RESULTS}/target"
"${RIBOSEEK}" createdb "${TARGET}" "${TARGETDB}"
"${RIBOSEEK}" rnasearch "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/tmp"

"${RIBOSEEK}" convertalis "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/results.m8"

"${EVALUATE}" "${RESULTS}/results.m8" "${RESULTS}/roc1_auc.tsv" | tee "${RESULTS}/evaluation.log"
ACTUAL=$(grep "^ROC1-AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f2 | cut -d"," -f1)
EXPECTED="0.485652"
awk -v actual="$ACTUAL" -v expected="$EXPECTED" \
    'BEGIN { print (actual == expected) ? "GOOD" : "BAD"; print "Expected: ", expected; print "Actual: ", actual; }' \
    > "${RESULTS}.report"
