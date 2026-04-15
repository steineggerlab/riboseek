#ifndef RIBOSEEK_DINUCLEOTIDEMAPPING_H
#define RIBOSEEK_DINUCLEOTIDEMAPPING_H

class Sequence;

// Register dinucleotide mapSequence callback via SeqAuxInfo
void registerDinucleotideMapping();

// Get the num2outputnum lookup table (maps dinuc encoding -> second nucleotide, numeric)
// Populated after registerDinucleotideMapping() + matrix setup
const unsigned char* getDinucToNucTable();

// Re-encode a Sequence with reverse-shifted dinucleotide pairing.
// Call after mapSequence() on target sequences when the query is on the reverse strand.
// Uses the raw sequence data (seq->seqData) to re-encode with (prev, curr) pairing.
void dinucEncodeReverse(Sequence *seq);

#endif
