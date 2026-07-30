#!/usr/bin/env python3
# Smoke / regression test: fold a stratified sample of ArchiveII sequences
# (28-282 nt) and check the MFE *and* the traceback structure against
# precomputed golden values in data/archiveii_smoke.tsv.
#
# The golden values are tornadofold's own deterministic output; every MFE was
# verified byte-exact against RNAfold -d2 when the file was generated (see
# util/gen note below). Structure is checked against tornadofold's golden rather
# than RNAfold's because co-optimal (equal-energy) structures are resolved by
# implementation-specific traceback tie-breaking. The sample avoids hairpin
# loops > 30 nt so the golden reproduces exactly across platforms/libm.
#
# Usage: python3 tests/smoke_test.py   (set TORNADOFOLD=path to override the binary)
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("TORNADOFOLD", os.path.join(HERE, "..", "tornadofold"))
DATA = os.path.join(HERE, "data", "archiveii_smoke.tsv")


def main():
    if not os.path.exists(BIN):
        print(f"error: tornadofold binary not found at {BIN} (build it first: `make`)")
        return 2
    rows = [l.rstrip("\n").split("\t") for l in open(DATA)][1:]
    fails = 0
    for name, seq, ref, mfe, db in rows:
        out = subprocess.run([BIN], input=seq + "\n", capture_output=True,
                             text=True).stdout.strip().split("\n")
        got = out[1].rsplit(" ", 1)
        got_db, got_mfe = got[0], got[1].strip("()")
        if got_mfe != mfe or got_db != db:
            fails += 1
            print(f"FAIL {name} (len {len(seq)}): "
                  f"MFE {got_mfe} vs {mfe}; structure {'ok' if got_db == db else 'DIFFERS'}")
    print(f"smoke: {len(rows) - fails}/{len(rows)} sequences reproduce golden MFE + structure")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
