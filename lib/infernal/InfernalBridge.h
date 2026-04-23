#ifndef MMSEQS_INFERNALBRIDGE_H
#define MMSEQS_INFERNALBRIDGE_H

#include <string>

namespace InfernalBridge {

struct WorkerPool; // opaque

bool isConfigured();

// Create N pre-forked Infernal worker processes. Fork them BEFORE allocating
// large memory (DBReader mmap etc.) so each worker stays small — this avoids
// the kernel-mm_struct contention that serializes concurrent forks from a fat
// parent. Returns NULL on failure. Thread-safe to use from OMP regions.
WorkerPool *startWorkerPool(int nWorkers);
void stopWorkerPool(WorkerPool *pool);

// If calibrate is true, worker runs infernal cmcalibrate after cmbuild
// (slow, ~1-2s/query) to fit exp-tail E-value parameters. If false, returns
// the uncalibrated CM and downstream consumers must supply their own E-value
// mapping.
bool buildCmFromStockholmText(WorkerPool *pool,
                              const std::string &stockholmText,
                              std::string &cmText,
                              std::string &error,
                              bool calibrate = false);

} // namespace InfernalBridge

#endif
