#ifndef RIBOSEEK_CMBUILDSCAN_H
#define RIBOSEEK_CMBUILDSCAN_H

#include "DBReader.h"
#include "SubstitutionMatrix.h"
#include "MsaFilter.h"
#include "Sequence.h"
#include <string>
#include <vector>

class LocalParameters;
struct cm_s;
typedef struct cm_s CM_t;

// per-thread build context
struct CmBuildCtx {
    DBReader<unsigned int>* qDbr = nullptr;
    DBReader<unsigned int>* tDbr = nullptr;
    DBReader<unsigned int>* resultReader = nullptr;
    SubstitutionMatrix* subMat = nullptr;
    MsaFilter* msaFilter = nullptr;
    Sequence* targetMapper = nullptr;
    std::vector<int> qid_vec;
    float minColCoverage = 0.30f;
    bool doMsaFilter = false;
    bool decodeTargetDinuc = false;
    bool targetGpuDb = false;
    int targetSeqType = 0;
};

extern bool buildQueryCmText(
    LocalParameters &par, CmBuildCtx &ctx, size_t id,
    unsigned int thread_idx, std::string &cmText, std::string &err,
    bool forScan = false, bool useInside = true, bool useLocal = false,
    CM_t **ret_cm = nullptr
);

#endif
