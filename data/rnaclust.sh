#!/bin/sh -e
# Clustering in dinucleotide space: encode the nucleotide database into the 25 letter
# dinucleotide alphabet, then run the MMseqs2 protein clustering workflow on it. The
# encoding keeps the database keys, so the resulting cluster database can be used with
# the original nucleotide database.
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"

if [ -n "$CONVERT" ]; then
    if notExists "${TMP_PATH}/input_dinuc.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" dinucdb "$INPUT" "${TMP_PATH}/input_dinuc" ${DINUCDB_PAR} \
            || fail "dinucdb died"
    fi
    INPUT="${TMP_PATH}/input_dinuc"
fi

# --strand 2: both strands of every sequence (key 2k forward, 2k+1 reverse complement); the cascade searches them
# against the forward strands (cascaded_clustering.sh)
if [ -n "$STRAND_BOTH" ]; then
    if notExists "${TMP_PATH}/input_dinuc_both.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" dinucdb "$1" "${TMP_PATH}/input_dinuc_both" ${DINUCDB_BOTH_PAR} \
            || fail "dinucdb --strand 2 died"
    fi
    STRAND_DB="${TMP_PATH}/input_dinuc_both"
    export STRAND_DB
else
    unset STRAND_DB
fi

# no $RUNNER here: the called workflow picks the MPI runner up from the environment and
# applies it to the modules that support it
mkdir -p "${TMP_PATH}/clust"
if notExists "$2.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${CLUST_MODULE}" "$INPUT" "$2" "${TMP_PATH}/clust" ${CLUST_PAR} \
        || fail "${CLUST_MODULE} died"
fi

if [ -n "$REMOVE_TMP" ]; then
    if [ -n "$CONVERT" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_dinuc" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_dinuc_h" ${VERBOSITY}
    fi
    if [ -n "$STRAND_BOTH" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_dinuc_both" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/clust"
    rm -f "${TMP_PATH}/rnaclust.sh"
fi
