#ifndef CREATEDB_OBSERVER_H
#define CREATEDB_OBSERVER_H

#include <cstddef>

struct CreatedbObserver {
    virtual ~CreatedbObserver() {}
    virtual void reset() = 0;
    virtual void addSequence(const char *seq, size_t len) = 0;
};

// NULL unless a caller (riboseek's createdb wrapper) sets it around the createdb call.
extern CreatedbObserver *createdbObserver;

#endif
