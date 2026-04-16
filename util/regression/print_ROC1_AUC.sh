#!/bin/sh -e
in_m8="$1"
out_roc1_auc_file="$2"

BASE="$(dirname "$(cd "$(dirname "$0")"; pwd)/$(basename "$0")")"
LOOKUP="${DATADIR}/rfam_lookup_target_clan.tsv"
TARGET_STRAND="${DATADIR}/target.tsv"
ROC1_AWK="${BASE}/roc1.awk"
TMP_M8="${out_roc1_auc_file}.tmp.m8"

sort -k1,1 -k11,11g "$in_m8" \
    | awk 'BEGIN {OFS="\t"} !($1$2 in a) {print $1, $2, $(NF-1); a[$1$2]}' \
    | awk 'BEGIN {OFS="\t"; argind=0} FNR==1 {argind++} argind==1 {a[$1]=$2;next;} argind==2 {b[$1]=$2;next;} argind==3 {print $0, b[$1], b[$2], a[$2]}' "$TARGET_STRAND" "$LOOKUP" - \
    > "$TMP_M8"

"$ROC1_AWK" "$LOOKUP" "$TMP_M8" | sort -k4,4rn > "$out_roc1_auc_file"
rm -f "$TMP_M8"

printf "ROC1-AUC: %s, number of queries: %s\n" \
    "$(awk '{sum+=$4} END {print sum/NR}' "$out_roc1_auc_file")" \
    "$(wc -l < "$out_roc1_auc_file")"
