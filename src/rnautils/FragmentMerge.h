#ifndef RIBOSEEK_FRAGMENTMERGE_H
#define RIBOSEEK_FRAGMENTMERGE_H

// Install the offsetalignment hook that stitches alignments fragmented by a
// splitsequence'd target back together via windowed realignment on the original target.
void registerFragmentMerger();

#endif
