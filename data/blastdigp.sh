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
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERYDB="$1"
ORIQUERYDB="$1"
# offsetalignment needs a nucleotide-typed source DB to build the ORF-key lookup
# (lib/mmseqs's offsetalignment gates the lookup block on queryDbType == NUCLEOTIDES).
# The profile we overwrite ORIQUERYDB with between iterations is HMM_PROFILE, so
# we keep a separate handle to the original nucleotide query DB for offsetalignment.
OFFSET_QUERYDB="$1"
TARGETDB="$2"
TMP_PATH="$4"
MERGE_INPUT=""

# if [ -n "$NEEDTARGETSPLIT" ]; then
#     if notExists "$TMP_PATH/target_seqs_split.dbtype"; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" splitsequence "$TARGETDB" "$TMP_PATH/target_seqs_split" ${SPLITSEQUENCE_PAR}  \
#             || fail "Split sequence died"
#     fi
#     TARGETDB="$TMP_PATH/target_seqs_split"
# fi

# Target DB stores raw nucleotide characters. Set DINUCLEOTIDE flag so
# Sequence::mapSequence uses dinucleotide pair encoding (via auxRegistry).
# Only do if the workflow doesn't use GPU
# Also make soft-link to make the target DB dinucleotide-aware without altering the original target DB
if [ -n "$SPLITSTRAND" ] && [ -z "$GPU" ]; then
    if notExists "$TMP_PATH/target_db_dinuc.dbtype"; then
        ln -s "$(abspath "$2")" "$TMP_PATH/target_db_dinuc"
        ln -s "$(abspath "$2.index")" "$TMP_PATH/target_db_dinuc.index"
        ln -s "$(abspath "${2}_h")" "$TMP_PATH/target_db_dinuc_h" 2>/dev/null || true
        ln -s "$(abspath "${2}_h.index")" "$TMP_PATH/target_db_dinuc_h.index" 2>/dev/null || true
        ln -s "$(abspath "${2}_h.dbtype")" "$TMP_PATH/target_db_dinuc_h.dbtype" 2>/dev/null || true
        ln -s "$(abspath "${2}.lookup")" "$TMP_PATH/target_db_dinuc.lookup" 2>/dev/null || true
        ln -s "$(abspath "${2}.source")" "$TMP_PATH/target_db_dinuc.source" 2>/dev/null || true
        awk 'BEGIN { printf("%c%c%c%c",0,0,64,0); exit; }' > "$TMP_PATH/target_db_dinuc.dbtype"
    fi
    TARGETDB="$TMP_PATH/target_db_dinuc"
fi

STEP=0
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ "$STEP" -lt "$NUM_IT" ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_tmp_${STEP}.done"; then
        if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
            TMP=""
            PREF="fake_pref"
        elif [ "$PREFMODE" = "UNGAPPED" ]; then
            PARAM="UNGAPPEDPREFILTER_PAR_$STEP"
            eval TMP="\$$PARAM"
            PREF="${MMSEQS} ungappedprefilter"
        else
            PARAM="PREFILTER_PAR_$STEP"
            eval TMP="\$$PARAM"
            PREF="${MMSEQS} prefilter"
        fi

        if [ "$STEP" -eq 0 ]; then
            if [ -n "$SPLITSTRAND" ]; then
                if notExists "$TMP_PATH/query_strand.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" splitstrand "$QUERYDB" "$TMP_PATH/query_strand" ${SPLIT_STRAND_PAR} \
                        || fail "Split strand died"
                fi
                QUERYDB="$TMP_PATH/query_strand"
            fi
            if [ -n "$NEEDQUERYSPLIT" ]; then
                if notExists "$TMP_PATH/query_seqs_split.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" splitsequence "$QUERYDB" "$TMP_PATH/query_seqs_split" ${SPLITSEQUENCE_PAR}  \
                    || fail "Split sequence died"
                fi
                QUERYDB="$TMP_PATH/query_seqs_split"
            fi
            # Force query to HMM_PROFILE + DINUCLEOTIDE + SRC_SEQUENCE + STRAND_SPLIT
            # so prefilter treats sequences as profiles (on-the-fly profile from submat)
            # HMM_PROFILE=2, DINUCLEOTIDE=64, SRC_SEQUENCE=32, STRAND_SPLIT=128 => ext=0xE0
            if [ -n "$SPLITSTRAND" ]; then
                if [ -L "$QUERYDB.dbtype" ]; then
                    rm -f "$QUERYDB.dbtype"
                fi
                LC_ALL=C awk 'BEGIN { printf("%c%c%c%c",2,0,224,0); exit; }' > "$QUERYDB.dbtype"
            fi
            # shellcheck disable=SC2086
            $RUNNER $PREF "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_$STEP" ${TMP} \
                || fail "Prefilter died"
        else
            if [ -n "$SPLITSTRAND" ]; then
                if notExists "$TMP_PATH/query_strand_$STEP.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" splitstrand "$QUERYDB" "$TMP_PATH/query_strand_$STEP" ${SPLIT_STRAND_PAR} \
                        || fail "Split strand died"
                fi
                QUERYDB="$TMP_PATH/query_strand_${STEP}"
            fi
            if [ -n "$NEEDQUERYSPLIT" ]; then
                if notExists "$TMP_PATH/query_seqs_split_$STEP.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" splitsequence "$QUERYDB" "$TMP_PATH/query_seqs_split_$STEP" ${SPLITSEQUENCE_PAR}  \
                    || fail "Split sequence died"
                fi
                QUERYDB="$TMP_PATH/query_seqs_split_${STEP}"
            fi
            # shellcheck disable=SC2086
            $RUNNER $PREF "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_tmp_$STEP" ${TMP} \
                || fail "Prefilter died"
        fi
        touch "$TMP_PATH/pref_tmp_${STEP}.done"
    fi

    if [ "$STEP" -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.done"; then
            STEPONE=$((STEP-1))
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_double_$STEPONE" "$TMP_PATH/pref_$STEP" $SUBSTRACT_PAR \
                || fail "Substract died"
            # # shellcheck disable=SC2086
            # "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/pref_$STEP" $SUBSTRACT_PAR \
            #     || fail "Substract died"
            "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_$STEP"
        fi
        touch "$TMP_PATH/pref_$STEP.done"
    fi

	# call alignment module
	if notExists "$TMP_PATH/aln_tmp_$STEP.done"; then
	    PARAM="ALIGNMENT_PAR_$STEP"
        eval TMP="\$$PARAM"

        if [ "$STEP" -eq 0 ]; then
            MERGE_INPUT="$TMP_PATH/aln_double_$STEP"
        else
            MERGE_INPUT="$TMP_PATH/aln_tmp_double_$STEP"
        fi
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_$STEP" "$MERGE_INPUT" ${TMP} \
            || fail "Alignment died"
        touch "$TMP_PATH/aln_tmp_$STEP.done"
    fi

    # merge alignments for the next round
    if [ "$STEP" -gt 0 ]; then
        STEPONE=$((STEP-1))
        if notExists "$TMP_PATH/merge_aln_double_$STEP.done"; then
            if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
                # merge alignments before offset
                "$MMSEQS" mergedbs "$QUERYDB" "$TMP_PATH/aln_double_$STEP" "$TMP_PATH/aln_double_$STEPONE" "$MERGE_INPUT" \
                    || fail "Alignment died"
                "$MMSEQS" rmdb "$TMP_PATH/aln_double_$STEPONE"
            else
                "$MMSEQS" rmdb "$TMP_PATH/aln_double_$STEPONE"
            fi
            touch "$TMP_PATH/merge_aln_double_$STEP.done"
        fi
    fi

    # offset alignment
    if notExists "$TMP_PATH/aln_offset_$STEP.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" offsetalignment "$OFFSET_QUERYDB" "$QUERYDB" "$2" "$TARGETDB" "$MERGE_INPUT" "$TMP_PATH/aln_offset_$STEP" ${OFFSETALIGNMENT_PAR} \
            || fail "Offset step died"
        # replace alignment with offset alignment
        if [ "$STEP" -eq 0 ]; then
            "$MMSEQS" mvdb "$TMP_PATH/aln_offset_$STEP" "$TMP_PATH/aln_$STEP"
        else
            "$MMSEQS" mvdb "$TMP_PATH/aln_offset_$STEP" "$TMP_PATH/aln_tmp_$STEP"
        fi
        touch "$TMP_PATH/aln_offset_$STEP.done"
    fi

    if [ "$STEP" -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.done"; then
            STEPONE=$((STEP-1))
            if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
                "$MMSEQS" mergedbs "$ORIQUERYDB" "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/aln_tmp_$STEP" \
                    || fail "Alignment died"
            else
                "$MMSEQS" mergedbs "$ORIQUERYDB" "$3" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/aln_tmp_$STEP" \
                    || fail "Alignment died"
            fi
            "$MMSEQS" rmdb "$TMP_PATH/aln_$STEPONE"
            "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_$STEP"
            touch "$TMP_PATH/aln_$STEP.done"
        fi
    fi

# create profiles
    if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
        if notExists "$TMP_PATH/profile_$STEP.dbtype"; then
            PARAM="PROFILE_PAR_$STEP"
            eval TMP="\$$PARAM"

            # For dinucleotide search, result2profile must read target sequences
            # in dinucleotide space so the PSSM is built over dinucleotide codes.
            R2P_TARGET="$2"
            R2P_QUERY="$ORIQUERYDB"
            if [ -n "$SPLITSTRAND" ] && [ -z "$GPU" ]; then
                # Create a dinuc-flagged symlink of the original target DB (once)
                # offsetalignment maps split coords back to original, so r2p needs the original.
                if notExists "$TMP_PATH/target_r2p.dbtype"; then
                    ln -s "$(abspath "$2")" "$TMP_PATH/target_r2p"
                    ln -s "$(abspath "$2.index")" "$TMP_PATH/target_r2p.index"
                    ln -s "$(abspath "${2}_h")" "$TMP_PATH/target_r2p_h" 2>/dev/null || true
                    ln -s "$(abspath "${2}_h.index")" "$TMP_PATH/target_r2p_h.index" 2>/dev/null || true
                    ln -s "$(abspath "${2}_h.dbtype")" "$TMP_PATH/target_r2p_h.dbtype" 2>/dev/null || true
                    ln -s "$(abspath "${2}.lookup")" "$TMP_PATH/target_r2p.lookup" 2>/dev/null || true
                    ln -s "$(abspath "${2}.source")" "$TMP_PATH/target_r2p.source" 2>/dev/null || true
                    awk 'BEGIN { printf("%c%c%c%c",0,0,64,0); exit; }' > "$TMP_PATH/target_r2p.dbtype"
                fi
                R2P_TARGET="$TMP_PATH/target_r2p"
                # For step 0, create a dinuc-flagged symlink of the original query
                if [ "$STEP" -eq 0 ]; then
                    if notExists "$TMP_PATH/query_r2p.dbtype"; then
                        ln -s "$(abspath "$1")" "$TMP_PATH/query_r2p"
                        ln -s "$(abspath "$1.index")" "$TMP_PATH/query_r2p.index"
                        ln -s "$(abspath "${1}_h")" "$TMP_PATH/query_r2p_h" 2>/dev/null || true
                        ln -s "$(abspath "${1}_h.index")" "$TMP_PATH/query_r2p_h.index" 2>/dev/null || true
                        ln -s "$(abspath "${1}_h.dbtype")" "$TMP_PATH/query_r2p_h.dbtype" 2>/dev/null || true
                        ln -s "$(abspath "${1}.lookup")" "$TMP_PATH/query_r2p.lookup" 2>/dev/null || true
                        ln -s "$(abspath "${1}.source")" "$TMP_PATH/query_r2p.source" 2>/dev/null || true
                        awk 'BEGIN { printf("%c%c%c%c",0,0,64,0); exit; }' > "$TMP_PATH/query_r2p.dbtype"
                    fi
                    R2P_QUERY="$TMP_PATH/query_r2p"
                fi
            else
                # For step 0, create a dinuc-flagged symlink of the original query
                if [ "$STEP" -eq 0 ]; then
                    if notExists "$TMP_PATH/query_r2p.dbtype"; then
                        ln -s "$(abspath "$1")" "$TMP_PATH/query_r2p"
                        ln -s "$(abspath "$1.index")" "$TMP_PATH/query_r2p.index"
                        ln -s "$(abspath "${1}_h")" "$TMP_PATH/query_r2p_h" 2>/dev/null || true
                        ln -s "$(abspath "${1}_h.index")" "$TMP_PATH/query_r2p_h.index" 2>/dev/null || true
                        ln -s "$(abspath "${1}_h.dbtype")" "$TMP_PATH/query_r2p_h.dbtype" 2>/dev/null || true
                        ln -s "$(abspath "${1}.lookup")" "$TMP_PATH/query_r2p.lookup" 2>/dev/null || true
                        ln -s "$(abspath "${1}.source")" "$TMP_PATH/query_r2p.source" 2>/dev/null || true
                        awk 'BEGIN { printf("%c%c%c%c",0,0,64,0); exit; }' > "$TMP_PATH/query_r2p.dbtype"
                    fi
                    R2P_QUERY="$TMP_PATH/query_r2p"
                fi
            fi

            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" result2profile "$R2P_QUERY" "$R2P_TARGET" "$TMP_PATH/aln_$STEP" "$TMP_PATH/profile_$STEP" ${TMP} \
                || fail "Create profile died"

            # Plain HMM_PROFILE type: splitstrand applies the per-column dinucleotide
            # revcomp permutation inline for the reverse-strand entry, so downstream
            # uses standard mapProfile (matching the forked extractqueryprofiles).
            # RnaAlign forces the DINUCLEOTIDE flag itself for alignment.
            if [ -n "$SPLITSTRAND" ]; then
                awk 'BEGIN { printf("%c%c%c%c",2,0,0,0); exit; }' > "$TMP_PATH/profile_$STEP.dbtype"
            fi
        fi
    fi
    ORIQUERYDB="$TMP_PATH/profile_$STEP"
    QUERYDB="$TMP_PATH/profile_$STEP"
	STEP=$((STEP+1))
done

# For single iteration (NUM_IT=1), result is in aln_0; move to output
if [ "$NUM_IT" -eq 1 ]; then
    "$MMSEQS" mvdb "$TMP_PATH/aln_0" "$3"
fi

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
    "$MMSEQS" rmdb "$TMP_PATH/query_strand" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$TMP_PATH/query_strand_h" ${VERBOSITY}
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
        "$MMSEQS" rmdb "$TMP_PATH/query_strand_${STEP}" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$TMP_PATH/query_strand_${STEP}_h" ${VERBOSITY}

        STEP=$((STEP+1))
    done
    rm -f "$TMP_PATH/blastpgp.sh"
fi

