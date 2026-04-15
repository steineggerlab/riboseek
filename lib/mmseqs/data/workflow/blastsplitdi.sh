#!/bin/sh -e
# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# # check number of input variables
# [ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
cur=1
last=$(( $# - 2 ))
chk=""
while [ "$cur" -le "$last" ]; do
    eval "chk=\${$cur}"
    [ ! -f "${chk}.dbtype" ] &&  echo "${chk}.dbtype not found!" && exit 1;
    cur=$(( cur + 1 ))
done
# [ ! -f "$1.dbtype" ] &&  echo "$1.dbtype not found!" && exit 1;
# [ ! -f "$2.dbtype" ] &&  echo "$2.dbtype not found!" && exit 1;
# [   -f "$3.dbtype" ] &&  echo "$3.dbtype exists already!" && exit 1;
# [ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";
i=$(( $# - 1 ))
eval "chk=\${$i}"
[ -f "${chk}.dbtype" ] &&  echo "${chk}.dbtype exists already!" && exit 1;
i=$#
eval "chk=\${$i}"
[ ! -d "${chk}" ] &&  echo "tmp directory ${chk} not found!" && mkdir -p "${chk}";


QUERY="$1"
# TARGET="$2"
TMP_PATH="${chk}" # tmp directory

if [ -n "$EXTRACTQUERYPROFILES" ]; then
    if notExists "$TMP_PATH/query_prof.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractqueryprofiles "$QUERY" "$TMP_PATH/query_prof" ${EXTRACT_QUERY_PROFILES_PAR}  \
            || fail "Extractqueryprofiles died"
    fi
    QUERY="$TMP_PATH/query_prof"
fi

if [ -n "$NEEDQUERYSPLIT" ]; then
    if notExists "$TMP_PATH/query_seqs_split.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" splitsequence "$QUERY" "$TMP_PATH/query_seqs_split" ${SPLITSEQUENCE_PAR}  \
        || fail "Split sequence died"
    fi
    QUERY="$TMP_PATH/query_seqs_split"
fi

# Iterate through targets
i=0
cur=2
last=$(( $# - 2 ))
res_base=$(( $# - 1 ))
eval "res_base=\${$res_base}"
while [ "$cur" -le "$last" ]; do
    eval "TARGET=\${$cur}"
    ori_target="$TARGET"
    if [ -n "$NEEDTARGETSPLIT" ]; then
        if notExists "$TMP_PATH/target_seqs_split.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" splitsequence "$TARGET" "$TMP_PATH/target_seqs_split" ${SPLITSEQUENCE_PAR}  \
                || fail "Split sequence died"
        fi
        TARGET="$TMP_PATH/target_seqs_split"
    fi

    mkdir -p "$TMP_PATH/search"
    if notExists "$TMP_PATH/aln_${i}.dbtype"; then
        # search does not need a parameter because the environment variables will be set by the workflow
        # shellcheck disable=SC2086
        "$SEARCH" "${QUERY}" "${TARGET}" "$TMP_PATH/aln_${i}" "$TMP_PATH/search"  \
            || fail "Search step died"
    fi

    if notExists "${res_base}_${i}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" offsetalignment "$1" "${QUERY}" "$ori_target" "${TARGET}" "$TMP_PATH/aln_${i}"  "${res_base}_${i}" ${OFFSETALIGNMENT_PAR} \
            || fail "Offset step died"
    fi
    i=$((i + 1))
    cur=$(( cur + 1 ))
done


if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/q_orfs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/q_orfs_aa" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/t_orfs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/t_orfs_aa" ${VERBOSITY}
fi

