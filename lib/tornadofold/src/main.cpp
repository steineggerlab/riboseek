// tornadofold driver: reads FASTA (or one seq per line) from stdin/file, prints
// dot-bracket and MFE. Scalar reference path (fold.h).
#include "tornadofold.h"
#include <cstdio>
#include <iostream>
#include <chrono>

int main(int argc, char** argv) {
    bool timing = false, evalMode = false, verbose = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-t") {
            timing = true;
        } else if (a == "-e") {
            evalMode = true;
        } else if (a == "-v") {
            verbose = true;
        }
    }
    if (evalMode) {
        // read seq then structure lines, print my model's energy
        std::string seqL, dbL;
        while (std::getline(std::cin, seqL)) {
            if (seqL.empty() || seqL[0] == '>') {
                continue;
            }
            if (!std::getline(std::cin, dbL)) {
                break;
            }
            std::string s;
            for (char c : seqL) {
                if (isspace((unsigned char)c)) {
                    continue;
                }
                char u = toupper(c);
                if (u == 'T') {
                    u = 'U';
                }
                s += u;
            }
            std::vector<int> bb(s.size());
            for (size_t i = 0; i < s.size(); ++i) {
                char c = s[i];
                bb[i] = (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : ((c == 'U') ? 3 : 4);
            }
            en::EM em;
            em.set(s, bb.data(), (int)s.size());
            int e = em.evalStructure(dbL, verbose);
            printf("%s (%.2f)\n", dbL.c_str(), e / 100.0);
        }
        return 0;
    }

    std::string line, name;
    std::vector<std::pair<std::string, std::string>> seqs;
    while (std::getline(std::cin, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            name = line.substr(1);
            continue;
        }
        // preserve length: uppercase, T->U, unknown kept as-is (treated as N)
        std::string s;
        for (char c : line) {
            if (isspace((unsigned char)c)) {
                continue;
            }
            char u = toupper(c);
            if (u == 'T') {
                u = 'U';
            }
            s += u;
        }
        if (!s.empty()) {
            seqs.push_back({name, s});
            name.clear();
        }
    }

    // Fold the batch (optionally across cores via OpenMP), preserving order.
    int N = (int)seqs.size();
    std::vector<std::string> dbs(N);
    std::vector<int> ens(N);
    std::vector<double> tms(N);
    auto tw0 = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        // one folder per thread, reused across sequences (buffers amortise)
        tornadofold::TornadoFold f;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int k = 0; k < N; ++k) {
            auto t0 = std::chrono::high_resolution_clock::now();
            int e = f.fold(seqs[k].second);
            auto t1 = std::chrono::high_resolution_clock::now();
            tms[k] = std::chrono::duration<double, std::milli>(t1 - t0).count();
            ens[k] = e;
            dbs[k] = f.traceback(e);
        }
    }
    auto tw1 = std::chrono::high_resolution_clock::now();
    double wall = std::chrono::duration<double, std::milli>(tw1 - tw0).count();
    double total = 0;
    for (int k = 0; k < N; ++k) {
        total += tms[k];
        printf("%s\n%s (%.2f)", seqs[k].second.c_str(), dbs[k].c_str(), ens[k] / 100.0);
        if (timing) {
            printf("  [%.2f ms]", tms[k]);
        }
        printf("\n");
    }
    if (timing) {
        fprintf(stderr, "sum fold time: %.2f ms | wall: %.2f ms\n", total, wall);
    }
    return 0;
}
