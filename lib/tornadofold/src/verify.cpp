// Self-consistency test: for random sequences the DP optimum must equal the
// energy of its own traceback structure re-scored by the energy model.
#include "tornadofold.h"
#include <cstdio>
#include <random>

int main() {
    std::mt19937 rng(42);
    const char* B = "ACGU";
    std::uniform_int_distribution<int> d(0, 3);
    int fails = 0, tot = 0;
    for (int L : {20, 35, 50, 80, 120, 200}) {
        for (int t = 0; t < 40; ++t) {
            std::string s;
            for (int i = 0; i < L; ++i) {
                s += B[d(rng)];
            }
            tornadofold::TornadoFold f;
            int e = f.fold(s);
            std::string db = f.traceback(e);
            int e2 = f.em.evalStructure(db); // re-score the traceback
            tot++;
            if (e != e2) {
                if (fails < 8) {
                    printf("MISMATCH L=%d dp=%d eval=%d\n%s\n%s\n", L, e, e2, s.c_str(),
                           db.c_str());
                }
                fails++;
            }
        }
    }
    printf("verify: %d/%d consistent (%d mismatches)\n", tot - fails, tot, fails);
    return fails ? 1 : 0;
}
