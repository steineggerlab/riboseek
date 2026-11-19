#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
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


# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

mkdir -p "${TMP_PATH}/linclust"
if notExists "${TMP_PATH}/clu_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" linclust "$INPUT" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/linclust" ${LINCLUST_PAR} \
        || fail "linclust died"
fi

if notExists "${TMP_PATH}/input_step_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
        || fail "createsubdb died"
fi

if [ "$CLUSTER_MODULE" = "cluster2" ]; then
    INPUT_PRECLUST="${TMP_PATH}/input_step_redundancy"
    if notExists "${TMP_PATH}/pref_rep.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "$INPUT_PRECLUST" "$INPUT_PRECLUST" "${TMP_PATH}/pref_rep" ${PREFILTER_PAR} \
            || fail "Prefilter died"
    fi

    if notExists "${TMP_PATH}/clu2_rep.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" align2clust "$INPUT_PRECLUST" "${TMP_PATH}/pref_rep" "${TMP_PATH}/clu2_rep" ${ALIGN2CLUST_PAR} --filter-cludb-file "${TMP_PATH}/clu_redundancy" --filter-seqdb-file "$INPUT" \
            || fail "align2clust step died"
    fi

    if notExists "$2.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergeclusters "$SOURCE" "$2" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/clu2_rep" $MERGECLU_PAR \
            || fail "Merging of clusters has died"
    fi
elif [ "$CLUSTER_MODULE" = "cluster1" ]; then
    INPUT="${TMP_PATH}/input_step_redundancy"
    STEP=0
    STEPS=${STEPS:-1}
    CLUSTER_STR=""
    while [ "$STEP" -lt "$STEPS" ]; do
        # riboseek cluster --strand 2 (STRAND_DB: key 2k forward, 2k+1 reverse complement): the queries are both strands
        # of the step's sequences, the targets their forward strands; mergestrandaln keeps the better strand per pair
        # for clust and addclustrand writes each member's strand (mergeclusters composes it across the steps)
        QUERY="$INPUT"
        if [ -n "$STRAND_DB" ]; then
            QUERY="${TMP_PATH}/input_step${STEP}_both"
            if notExists "$QUERY.dbtype"; then
                awk '{ printf("%.0f\n%.0f\n", 2 * $1, 2 * $1 + 1) }' "$INPUT.index" > "$QUERY.keys"
                # shellcheck disable=SC2086
                "$MMSEQS" createsubdb "$QUERY.keys" "$STRAND_DB" "$QUERY" ${VERBOSITY} --subdb-mode 1 \
                    || fail "Strand step $STEP died"
            fi
        fi
        PARAM=PREFILTER${STEP}_PAR
        eval TMP="\$$PARAM"
        if notExists "${TMP_PATH}/pref_step$STEP.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "$QUERY" "$INPUT" "${TMP_PATH}/pref_step$STEP" ${TMP} \
                || fail "Prefilter step $STEP died"
        fi
        PARAM=ALIGNMENT${STEP}_PAR
        eval TMP="\$$PARAM"
        if [ -n "$STRAND_DB" ]; then
            if notExists "${TMP_PATH}/aln_step${STEP}_strand.dbtype"; then
                # shellcheck disable=SC2086
                $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERY" "$INPUT" "${TMP_PATH}/pref_step$STEP" "${TMP_PATH}/aln_step${STEP}_strand" ${TMP} \
                    || fail "Alignment step $STEP died"
            fi
            if notExists "${TMP_PATH}/aln_step$STEP.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" mergestrandaln "$INPUT" "${TMP_PATH}/aln_step${STEP}_strand" "${TMP_PATH}/aln_step$STEP" ${TMP} \
                    || fail "Strand merge step $STEP died"
            fi
        elif notExists "${TMP_PATH}/aln_step$STEP.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$INPUT" "${TMP_PATH}/pref_step$STEP" "${TMP_PATH}/aln_step$STEP" ${TMP} \
                || fail "Alignment step $STEP died"
        fi
        PARAM=CLUSTER${STEP}_PAR
        eval TMP="\$$PARAM"
        if [ -n "$STRAND_DB" ]; then
            if notExists "${TMP_PATH}/clu_step${STEP}_nostrand.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln_step$STEP" "${TMP_PATH}/clu_step${STEP}_nostrand" ${TMP} \
                    || fail "Clustering step $STEP died"
            fi
            if notExists "${TMP_PATH}/clu_step$STEP.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" addclustrand "${TMP_PATH}/clu_step${STEP}_nostrand" "${TMP_PATH}/aln_step${STEP}_strand" "${TMP_PATH}/clu_step$STEP" ${MERGECLU_PAR} \
                    || fail "Strand clustering step $STEP died"
            fi
        elif notExists "${TMP_PATH}/clu_step$STEP.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln_step$STEP" "${TMP_PATH}/clu_step$STEP" ${TMP} \
                || fail "Clustering step $STEP died"
        fi

        # FIXME: This won't work if paths contain spaces
        CLUSTER_STR="${CLUSTER_STR} ${TMP_PATH}/clu_step$STEP"
        NEXTINPUT="${TMP_PATH}/input_step$((STEP+1))"
        if [ "$STEP" -eq "$((STEPS-1))" ]; then
        if [ -n "$REASSIGN" ]; then
            if notExists "${TMP_PATH}/clu.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} \
                || fail "Merging of clusters has died"
            fi
        else
                # shellcheck disable=SC2086
                "$MMSEQS" mergeclusters "$SOURCE" "$2" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} $MERGECLU_PAR \
                || fail "Merging of clusters has died"
        fi
        else
            if notExists "$NEXTINPUT.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "$INPUT" "$NEXTINPUT" ${VERBOSITY} --subdb-mode 1 \
                    || fail "Order step $STEP died"
            fi
        fi

        INPUT="$NEXTINPUT"
        STEP=$((STEP+1))
    done

    if [ -n "$REASSIGN" ]; then
        STEP=$((STEP-1))
        PARAM=ALIGNMENT${STEP}_PAR
        eval ALIGNMENT_PAR="\$$PARAM"
        CLU="${TMP_PATH}/clu"
        if [ -n "$STRAND_DB" ]; then
            # riboseek --strand 2: every member is aligned to its representative in both orientations (the
            # representative's strand entries 2r and 2r+1 against the member's forward strand), the better one kept;
            # the rest runs on the cluster db without its strand column, and members that fail are searched again on
            # both strands
            if notExists "${TMP_PATH}/clu_strandq.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" splitclustrand "${TMP_PATH}/clu" "${TMP_PATH}/clu_strandq" ${THREADSANDCOMPRESS} \
                    || fail "splitclustrand reassign died"
            fi
            if notExists "${TMP_PATH}/aln_strand.dbtype"; then
                # shellcheck disable=SC2086
                $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$STRAND_DB" "$SOURCE" "${TMP_PATH}/clu_strandq" "${TMP_PATH}/aln_strand" ${ALIGNMENT_REASSIGN_PAR} \
                    || fail "align1 reassign died"
            fi
            if notExists "${TMP_PATH}/aln.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" mergestrandaln "$SOURCE" "${TMP_PATH}/aln_strand" "${TMP_PATH}/aln" ${ALIGNMENT_REASSIGN_PAR} \
                    || fail "mergestrandaln reassign died"
            fi
            if notExists "${TMP_PATH}/clu_plain.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" filterdb "${TMP_PATH}/clu" "${TMP_PATH}/clu_plain" --trim-to-one-column ${THREADSANDCOMPRESS} \
                    || fail "filterdb strand reassign died"
            fi
            CLU="${TMP_PATH}/clu_plain"
        # align to cluster sequences
        elif notExists "${TMP_PATH}/aln.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$SOURCE" "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/aln" ${ALIGNMENT_REASSIGN_PAR} \
                    || fail "align1 reassign died"
        fi
        # create file of cluster that do not align based on given criteria
        if notExists "${TMP_PATH}/clu_not_accepted.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "$CLU" "${TMP_PATH}/aln" "${TMP_PATH}/clu_not_accepted" --e-profile 100000000 -e 100000000 ${THREADSANDCOMPRESS} \
                    || fail "subtractdbs1 reassign died"
        fi
        if notExists "${TMP_PATH}/clu_not_accepted_swap.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" swapdb "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_not_accepted_swap" ${THREADSANDCOMPRESS} \
                    || fail "swapdb1 reassign died"
        fi
        # short circuit if nothing can be reassigned
        if [ ! -s "${TMP_PATH}/clu_not_accepted_swap.index" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" mvdb "${TMP_PATH}/clu" "$2" ${VERBOSITY}
            if [ -n "$REMOVE_TMP" ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted_swap" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
                if [ -n "$STRAND_DB" ]; then
                    for db in clu_strandq aln_strand clu_plain; do
                        # shellcheck disable=SC2086
                        "$MMSEQS" rmdb "${TMP_PATH}/$db" ${VERBOSITY}
                    done
                fi
            fi
        else
            # create file of cluster that do align based on given criteria
            if notExists "${TMP_PATH}/clu_accepted.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" subtractdbs "$CLU" "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_accepted" --e-profile 100000000 -e 100000000 ${THREADSANDCOMPRESS} \
                        || fail "subtractdbs2 reassign died"
            fi
            # create sequences database that were wrong assigned
            if notExists "${TMP_PATH}/seq_wrong_assigned.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" createsubdb "${TMP_PATH}/clu_not_accepted_swap" "$SOURCE" "${TMP_PATH}/seq_wrong_assigned" ${VERBOSITY} \
                        || fail "createsubdb1 reassign died"
            fi
            # build seed sequences
            if notExists "${TMP_PATH}/seq_seeds.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" createsubdb "$CLU" "$SOURCE" "${TMP_PATH}/seq_seeds" ${VERBOSITY} \
                        || fail "createsubdb2 reassign died"
            fi
            PARAM=PREFILTER${STEP}_PAR
            eval PREFILTER_PAR="\$$PARAM"
            # riboseek --strand 2: the wrongly assigned sequences are searched on both strands
            WRONG="${TMP_PATH}/seq_wrong_assigned"
            if [ -n "$STRAND_DB" ]; then
                WRONG="${TMP_PATH}/seq_wrong_assigned_both"
                if notExists "$WRONG.dbtype"; then
                    awk '{ printf("%.0f\n%.0f\n", 2 * $1, 2 * $1 + 1) }' "${TMP_PATH}/seq_wrong_assigned.index" > "$WRONG.keys"
                    # shellcheck disable=SC2086
                    "$MMSEQS" createsubdb "$WRONG.keys" "$STRAND_DB" "$WRONG" ${VERBOSITY} --subdb-mode 1 \
                        || fail "createsubdb strand reassign died"
                fi
            fi
            # try to find best matching centroid sequences for prev. wrong assigned sequences
            if notExists "${TMP_PATH}/seq_wrong_assigned_pref.dbtype"; then
                if notExists "${TMP_PATH}/seq_seeds.merged.dbtype"; then
                    # combine seq dbs
                    MAXOFFSET=$(awk '($2+$3) > max { max = $2+$3 } END { printf("%.0f\n", max); }' "${TMP_PATH}/seq_seeds.index")
                    awk -v OFFSET="${MAXOFFSET}" 'FNR == NR { print $0; next } { printf("%s\t%.0f\t%s\n", $1, $2+OFFSET, $3); }' "${TMP_PATH}/seq_seeds.index" \
                        "${TMP_PATH}/seq_wrong_assigned.index" > "${TMP_PATH}/seq_seeds.merged.index"
                    ln -s "$(abspath "${TMP_PATH}/seq_seeds")" "${TMP_PATH}/seq_seeds.merged.0"
                    ln -s "$(abspath "${TMP_PATH}/seq_wrong_assigned")" "${TMP_PATH}/seq_seeds.merged.1"
                    cp "${TMP_PATH}/seq_seeds.dbtype" "${TMP_PATH}/seq_seeds.merged.dbtype"
                fi
                # shellcheck disable=SC2086
                $RUNNER "$MMSEQS" prefilter "$WRONG" "${TMP_PATH}/seq_seeds.merged" "${TMP_PATH}/seq_wrong_assigned_pref" ${PREFILTER_REASSIGN_PAR} \
                        || fail "Prefilter reassign died"
            fi
            if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref" "${TMP_PATH}/seq_wrong_assigned_pref_swaped" ${THREADSANDCOMPRESS} \
                        || fail "swapdb2 reassign died"
            fi
            if [ -n "$STRAND_DB" ]; then
                # targets are strand entries 2m / 2m+1: keyed by them (swapdb), merged to one hit per sequence pair
                # (mergestrandaln) and keyed back by the seeds; the strand-keyed alignments are kept for the strands
                if notExists "${TMP_PATH}/seq_wrong_assigned_strand_aln.dbtype"; then
                    # shellcheck disable=SC2086
                    $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${TMP_PATH}/seq_seeds.merged" "$WRONG" \
                                                        "${TMP_PATH}/seq_wrong_assigned_pref_swaped" "${TMP_PATH}/seq_wrong_assigned_strand_aln" ${ALIGNMENT_REASSIGN_PAR} \
                            || fail "align2 reassign died"
                fi
                if notExists "${TMP_PATH}/seq_wrong_assigned_strand_aln_swap.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_strand_aln" "${TMP_PATH}/seq_wrong_assigned_strand_aln_swap" ${THREADSANDCOMPRESS} \
                            || fail "swapdb3 reassign died"
                fi
                if notExists "${TMP_PATH}/seq_wrong_assigned_strand_aln_merged.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" mergestrandaln "${TMP_PATH}/seq_wrong_assigned" "${TMP_PATH}/seq_wrong_assigned_strand_aln_swap" \
                                             "${TMP_PATH}/seq_wrong_assigned_strand_aln_merged" ${ALIGNMENT_REASSIGN_PAR} \
                            || fail "mergestrandaln2 reassign died"
                fi
                if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_strand_aln_merged" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${THREADSANDCOMPRESS} \
                            || fail "swapdb4 reassign died"
                fi
            elif notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln.dbtype"; then
                # shellcheck disable=SC2086
                $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${TMP_PATH}/seq_seeds.merged" "${TMP_PATH}/seq_wrong_assigned" \
                                                    "${TMP_PATH}/seq_wrong_assigned_pref_swaped" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${ALIGNMENT_REASSIGN_PAR} \
                        || fail "align2 reassign died"
            fi

            if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" filterdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol" --trim-to-one-column ${THREADSANDCOMPRESS} \
                            || fail "filterdb2 reassign died"
            fi

            if notExists "${TMP_PATH}/clu_accepted_plus_wrong.dbtype"; then
                # combine clusters
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "${TMP_PATH}/seq_seeds.merged" "${TMP_PATH}/clu_accepted_plus_wrong" "${TMP_PATH}/clu_accepted" \
                                "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol" ${MERGEDBS_PAR} \
                                    || fail "mergedbs reassign died"
            fi

            if notExists "${TMP_PATH}/missing.single.seqs.db.dbtype"; then
                awk 'FNR==NR{if($3 > 1){ f[$1]=1; }next} !($1 in f){print $1"\t"$1}' "${TMP_PATH}/clu_accepted_plus_wrong.index" "${SOURCE}.index" > "${TMP_PATH}/missing.single.seqs"
                # shellcheck disable=SC2086
                "$MMSEQS" tsv2db "${TMP_PATH}/missing.single.seqs" "${TMP_PATH}/missing.single.seqs.db" --output-dbtype 6 ${VERBCOMPRESS} \
                                    || fail "tsv2db reassign died"
            fi

            if notExists "${TMP_PATH}/clu_accepted_plus_wrong_plus_single.dbtype"; then
                # combine clusters
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "${SOURCE}" "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" "${TMP_PATH}/clu_accepted_plus_wrong" \
                                "${TMP_PATH}/missing.single.seqs.db" ${MERGEDBS_PAR} \
                                    || fail "mergedbs2 reassign died"
            fi

            PARAM=CLUSTER${STEP}_PAR
            eval TMP="\$$PARAM"
            if [ -n "$STRAND_DB" ]; then
                # strands from the alignments behind every final pair: members kept (aln_strand, keyed by the
                # representative's strand entry) and members reassigned (keyed by their own strand entry)
                if notExists "${TMP_PATH}/clu_reassign_nostrand.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" clust "${SOURCE}" "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" "${TMP_PATH}/clu_reassign_nostrand" ${TMP} \
                            || fail "Clustering step $STEP died"
                fi
                if notExists "${TMP_PATH}/strand_evidence.dbtype"; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" mergedbs "$STRAND_DB" "${TMP_PATH}/strand_evidence" "${TMP_PATH}/aln_strand" \
                                       "${TMP_PATH}/seq_wrong_assigned_strand_aln_swap" ${MERGEDBS_PAR} \
                            || fail "mergedbs strand reassign died"
                fi
                # shellcheck disable=SC2086
                "$MMSEQS" addclustrand "${TMP_PATH}/clu_reassign_nostrand" "${TMP_PATH}/strand_evidence" "${2}" ${THREADSANDCOMPRESS} \
                        || fail "addclustrand reassign died"
            else
                # shellcheck disable=SC2086
                "$MMSEQS" clust "${SOURCE}" "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" "${2}" ${TMP} \
                        || fail "Clustering step $STEP died"
            fi

            if [ -n "$REMOVE_TMP" ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted_swap" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_seeds" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_seeds.merged" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_ocol" ${VERBOSITY}
                rm -f "${TMP_PATH}/missing.single.seqs"
                rm -f "${TMP_PATH}/clu_accepted_plus_wrong.tsv"
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/missing.single.seqs.db" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted_plus_wrong" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted_plus_wrong_plus_single" ${VERBOSITY}
                if [ -n "$STRAND_DB" ]; then
                    for db in clu_strandq aln_strand clu_plain seq_wrong_assigned_both seq_wrong_assigned_strand_aln \
                              seq_wrong_assigned_strand_aln_swap seq_wrong_assigned_strand_aln_merged clu_reassign_nostrand \
                              strand_evidence; do
                        # shellcheck disable=SC2086
                        "$MMSEQS" rmdb "${TMP_PATH}/$db" ${VERBOSITY}
                    done
                    rm -f "${TMP_PATH}/seq_wrong_assigned_both.keys"
                fi
            fi
        fi
    fi
fi


if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy_h" ${VERBOSITY}

    if [ "$CLUSTER_MODULE" = "cluster2" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_rep" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/clu2_rep" ${VERBOSITY}
    elif [ "$CLUSTER_MODULE" = "cluster1" ]; then

        STEP=0
        while [ "$STEP" -lt "$STEPS" ]; do
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/pref_step$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aln_step$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_step$STEP" ${VERBOSITY}
            if [ -n "$STRAND_DB" ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aln_step${STEP}_strand" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/clu_step${STEP}_nostrand" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/input_step${STEP}_both" ${VERBOSITY}
                rm -f "${TMP_PATH}/input_step${STEP}_both.keys"
            fi
            STEP=$((STEP+1))
        done

        STEP=1
        while [ "$STEP" -lt "$STEPS" ]; do
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_step$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_step${STEP}_h" ${VERBOSITY}
            STEP=$((STEP+1))
        done
    fi

    rm -f "${TMP_PATH}/cascaded_clustering.sh"
fi