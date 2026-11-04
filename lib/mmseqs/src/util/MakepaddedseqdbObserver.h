#ifndef MAKEPADDEDSEQDB_OBSERVER_H
#define MAKEPADDEDSEQDB_OBSERVER_H

#include <cstddef>

struct MakepaddedseqdbObserver {
    virtual ~MakepaddedseqdbObserver() {}
    // called once, after parameters are parsed and before any sequence is seen
    virtual void begin(unsigned int threads) = 0;
    // the source entry's raw data, from the thread that pads it
    virtual void addSequence(unsigned int threadIdx, const char *seq, size_t len) = 0;
};

// NULL unless a caller (riboseek's makepaddedseqdb wrapper) sets it around the makepaddedseqdb call.
extern MakepaddedseqdbObserver *makepaddedseqdbObserver;

#endif
