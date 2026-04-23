#ifndef RIBOSEEK_LOCALPARAMETERS_H
#define RIBOSEEK_LOCALPARAMETERS_H

#include "Parameters.h"

class LocalParameters : public Parameters {
public:
    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initParameterSingleton();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    LocalParameters();

    static void initInstance() {
        new LocalParameters;
    }

    static const unsigned int DBTYPE_EXTENDED_DINUCLEOTIDE = 64;
    static const unsigned int DBTYPE_EXTENDED_STRAND_SPLIT = 128;
    static const int CM_MODE_INSIDE = 1;

    float cmRegionFlanking;
    int cmMode;
    size_t dbSize;
    bool calibrateCm;

    PARAMETER(PARAM_CM_REGION)
    PARAMETER(PARAM_CM_MODE)
    PARAMETER(PARAM_DB_SIZE)
    PARAMETER(PARAM_CALIBRATE_CM)

    std::vector<MMseqsParameter*> splitstrand;
    std::vector<MMseqsParameter*> rnaalign;
    std::vector<MMseqsParameter*> cmbuild;
};

#endif
