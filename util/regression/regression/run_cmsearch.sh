#!/bin/sh -e
QUERY="${RESULTS}/query.fasta"
awk '/^>/{h=$0; getline s; if(length(s)<=500 && n<15){print h; print s; n++}}' \
    "${EXAMPLEDIR}/QUERY.fasta" > "${QUERY}"
QUERYDB="${RESULTS}/query"
"${RIBOSEEK}" createdb "${QUERY}" "${QUERYDB}"

TARGET="${EXAMPLEDIR}/DB.fasta"
TARGETDB="${RESULTS}/target"
"${RIBOSEEK}" createdb "${TARGET}" "${TARGETDB}"
"${RIBOSEEK}" search "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/tmp"
"${RIBOSEEK}" cmbuild "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/cmDB"
"${RIBOSEEK}" cmscan "${RESULTS}/cmDB" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/cmAln"

"${RIBOSEEK}" convertalis "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results" "${RESULTS}/search.m8"
"${RIBOSEEK}" convertalis "${QUERYDB}" "${TARGETDB}" "${RESULTS}/cmAln" "${RESULTS}/cm.m8"

# total aligned columns
SEARCH=$(awk '{s+=$4} END {print s}' "${RESULTS}/search.m8")
ACTUAL=$(awk '{s+=$4} END {print s}' "${RESULTS}/cm.m8")
EXPECTED="48686"
awk -v actual="$ACTUAL" -v expected="$EXPECTED" -v search="$SEARCH" \
    'BEGIN { print (actual == expected && actual > search) ? "GOOD" : "BAD"; print "Expected: ", expected; print "Actual: ", actual; print "Search baseline: ", search; }' \
    > "${RESULTS}.report"
