#!/bin/sh -e
# Iterative sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "$(abspath "${TDB}.index")" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(wc -c < "${TDB}.index")"
    awk -v size="$INDEX_SIZE" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
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
chk=""
eval "chk=\${$i}"
[ -f "${chk}.dbtype" ] &&  echo "${chk}.dbtype exists already!" && exit 1;
i=$#
eval "chk=\${$i}"
[ ! -d "${chk}" ] &&  echo "tmp directory ${chk} not found!" && mkdir -p "${chk}";

QUERYDB="$1"
ORIQUERYDB="$1"
TARGETDB="$2"
TMP_PATH="${chk}" # tmp directory
MERGE_INPUT=""

# if [ -n "$NEEDTARGETSPLIT" ]; then
#     if notExists "$TMP_PATH/target_seqs_split.dbtype"; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" splitsequence "$TARGETDB" "$TMP_PATH/target_seqs_split" ${SPLITSEQUENCE_PAR}  \
#             || fail "Split sequence died"
#     fi
#     TARGETDB="$TMP_PATH/target_seqs_split"
# fi

STEP=0
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ "$STEP" -lt "$NUM_IT" ]; do
    # Iterate through targets
    i=0
    cur=2
    last=$(( $# - 2 ))
    res_base=$(( $# - 1 ))
    eval "res_base=\${$res_base}"

    TARGETDBS=""
    ALIGNDBS=""

    if [ -n "$EXTRACTQUERYPROFILES" ]; then
        if notExists "$TMP_PATH/prof_${STEP}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" extractqueryprofiles "$QUERYDB" "$TMP_PATH/prof_${STEP}" ${EXTRACT_QUERY_PROFILES_PAR}  \
                || fail "Extractqueryprofiles died"
        fi
        QUERYDB="$TMP_PATH/prof_${STEP}"
    fi

    if [ -n "$NEEDQUERYSPLIT" ]; then
        if notExists "$TMP_PATH/query_seqs_split_${STEP}_${i}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" splitsequence "$QUERYDB" "$TMP_PATH/query_seqs_split_${STEP}" ${SPLITSEQUENCE_PAR}  \
            || fail "Split sequence died"
        fi
        QUERYDB="$TMP_PATH/query_seqs_split_${STEP}"
    fi

    while [ "$cur" -le "$last" ]; do
        eval "TARGETDB=\${$cur}"
        ori_target="$TARGETDB"

        if [ $cur -eq 2 ]; then
            cp ${QUERYDB}.index ${QUERYDB}.index.tmp
            awk '(FNR % 2 == 1)' ${QUERYDB}.index.tmp > ${QUERYDB}.index
        fi

        PREF_TMP=""
        # call prefilter module
        if notExists "$TMP_PATH/pref_tmp_${STEP}_${i}.done"; then
            if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
                TMP=""
                PREF="fake_pref"
            elif [ "$PREFMODE" = "UNGAPPED" ]; then
                PARAM=""
                if [ $cur -eq 2 ]; then
                    PARAM="UNGAPPEDPREFILTER_STR_PAR_$STEP"
                else
                    PARAM="UNGAPPEDPREFILTER_PAR_$STEP"
                fi
                eval TMP="\$$PARAM"
                PREF="${MMSEQS} ungappedprefilter"
            else
                PARAM="PREFILTER_PAR_$STEP"
                eval TMP="\$$PARAM"
                PREF="${MMSEQS} prefilter"
            fi

            if [ "$STEP" -eq 0 ]; then
                # shellcheck disable=SC2086
                $RUNNER $PREF "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_${STEP}_${i}" ${TMP} \
                    || fail "Prefilter died"
                PREF_TMP="$TMP_PATH/pref_${STEP}_${i}"
            else
                # shellcheck disable=SC2086
                $RUNNER $PREF "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_tmp_${STEP}_${i}" ${TMP} \
                    || fail "Prefilter died"
                PREF_TMP="$TMP_PATH/pref_tmp_${STEP}_${i}"
            fi
            touch "$TMP_PATH/pref_tmp_${STEP}_${i}.done"
        fi

        if [ $cur -eq 2 ]; then
            mv ${QUERYDB}.index.tmp ${QUERYDB}.index
            # fake pref for reverse strands
            awk 'BEGIN {OFS="\t"} {print;print $1+1, $2+$3-1, 0}' ${PREF_TMP}.index > ${PREF_TMP}.index.tmp
            mv ${PREF_TMP}.index.tmp ${PREF_TMP}.index
        fi

        if [ "$STEP" -ge 1 ]; then
            if notExists "$TMP_PATH/pref_${STEP}_${i}.done"; then
                STEPONE=$((STEP-1))
                # shellcheck disable=SC2086
                "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_${STEP}_${i}" "$TMP_PATH/aln_double_${STEPONE}_${i}" "$TMP_PATH/pref_${STEP}_${i}" $SUBSTRACT_PAR \
                    || fail "Substract died"
                # # shellcheck disable=SC2086
                # "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_${STEP}_${i}" "$TMP_PATH/aln_${STEP}_${i}ONE" "$TMP_PATH/pref_${STEP}_${i}" $SUBSTRACT_PAR \
                #     || fail "Substract died"
                "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_${STEP}_${i}"
            fi
            touch "$TMP_PATH/pref_${STEP}_${i}.done"
        fi

        # call alignment module
        if notExists "$TMP_PATH/aln_tmp_${STEP}_${i}.done"; then
            PARAM=""
            if [ $cur -eq 2 ]; then
                PARAM="ALIGNMENT_STR_PAR_$STEP"
            else
                PARAM="ALIGNMENT_PAR_$STEP"
            fi
            eval TMP="\$$PARAM"

            if [ "$STEP" -eq 0 ]; then
                MERGE_INPUT="$TMP_PATH/aln_double_${STEP}_${i}"
                # shellcheck disable=SC2086
                $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_${STEP}_${i}" "$TMP_PATH/aln_double_${STEP}_${i}" ${TMP} \
                    || fail "Alignment died"
            else
                MERGE_INPUT="$TMP_PATH/aln_tmp_double_${STEP}_${i}"
                # shellcheck disable=SC2086
                $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_${STEP}_${i}" "$TMP_PATH/aln_tmp_double_${STEP}_${i}" ${TMP} \
                    || fail "Alignment died"
            fi
            touch "$TMP_PATH/aln_tmp_${STEP}_${i}.done"
        fi

        # merge alignments for the next round
        if [ "$STEP" -gt 0 ]; then
            STEPONE=$((STEP-1))
            if notExists "$TMP_PATH/merge_aln_double_${STEP}_${i}.done"; then
                if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
                    # merge alignments before offset
                    "$MMSEQS" mergedbs "$QUERYDB" "$TMP_PATH/aln_double_${STEP}_${i}" "$TMP_PATH/aln_double_${STEPONE}_${i}" "$MERGE_INPUT" \
                        || fail "Alignment died"
                    "$MMSEQS" rmdb "$TMP_PATH/aln_double_${STEPONE}_${i}"
                else
                    "$MMSEQS" rmdb "$TMP_PATH/aln_double_${STEPONE}_${i}"
                fi
                touch "$TMP_PATH/merge_aln_double_${STEP}_${i}.done"
            fi
        fi

        # offset alignment
        if notExists "$TMP_PATH/aln_offset_${STEP}_${i}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" offsetalignment "$ORIQUERYDB" "$QUERYDB" "$ori_target" "$TARGETDB" "$MERGE_INPUT" "$TMP_PATH/aln_offset_${STEP}_${i}" ${OFFSETALIGNMENT_PAR} \
                || fail "Offset step died"
            # replace alignment with offset alignment
            if [ "$STEP" -eq 0 ]; then
                "$MMSEQS" mvdb "$TMP_PATH/aln_offset_${STEP}_${i}" "$TMP_PATH/aln_${STEP}_${i}"
            else
                "$MMSEQS" mvdb "$TMP_PATH/aln_offset_${STEP}_${i}" "$TMP_PATH/aln_tmp_${STEP}_${i}"
            fi
            touch "$TMP_PATH/aln_offset_${STEP}_${i}.done"
        fi

        if [ "$STEP" -gt 0 ]; then
            if notExists "$TMP_PATH/aln_${STEP}_${i}.done"; then
                STEPONE=$((STEP-1))
                if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
                    "$MMSEQS" mergedbs "$ORIQUERYDB" "$TMP_PATH/aln_${STEP}_${i}" "$TMP_PATH/aln_${STEPONE}_${i}" "$TMP_PATH/aln_tmp_${STEP}_${i}" \
                        || fail "Alignment died"
                else
                    "$MMSEQS" mergedbs "$ORIQUERYDB" "${res_base}_${i}" "$TMP_PATH/aln_${STEPONE}_${i}" "$TMP_PATH/aln_tmp_${STEP}_${i}" \
                        || fail "Alignment died"
                fi
                "$MMSEQS" rmdb "$TMP_PATH/aln_${STEPONE}_${i}"
                "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_${STEP}_${i}"
                touch "$TMP_PATH/aln_${STEP}_${i}.done"
            fi
        fi

        TARGETDBS="$TARGETDBS $ori_target"
        ALIGNDBS="$ALIGNDBS $TMP_PATH/aln_${STEP}_${i}"

        i=$((i + 1))
        cur=$(( cur + 1 ))
    done

    # create profiles
    if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
        if notExists "$TMP_PATH/profile_${STEP}_${i}.dbtype"; then
            PARAM="PROFILE_PAR_$STEP"
            eval TMP="\$$PARAM"
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" results2profile "$ORIQUERYDB" $TARGETDBS $ALIGNDBS "$TMP_PATH/profile_${STEP}" ${TMP} \
                || fail "Create profile died"
        fi
    fi
    QUERYDB="$TMP_PATH/profile_${STEP}"
	STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/query_seqs_split" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/query_seqs_split_h" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/target_seqs_split" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/target_seqs_split_h" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/query_prof" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/query_prof_h" ${VERBOSITY}
    while [ "$STEP" -lt "$NUM_IT" ]; do
        if [ "$STEP" -gt 0 ]; then
            rm -f -- "$TMP_PATH/aln_$STEP.done" "$TMP_PATH/pref_$STEP.done"
        fi
        rm -f -- "$TMP_PATH/aln_tmp_$STEP.done" "$TMP_PATH/pref_tmp_${STEP}.done" "$TMP_PATH/aln_offset_$STEP.done"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_tmp_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/query_seqs_split_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/query_seqs_split_${STEP}_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/target_seqs_split_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/target_seqs_split_${STEP}_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/query_prof_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/query_prof_${STEP}_h" ${VERBOSITY}

        STEP=$((STEP+1))
    done
    rm -f "$TMP_PATH/blastpgp.sh"
fi

