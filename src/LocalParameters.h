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
    int cmAlign;
    int cmAlignBanded;
    int cmLocal;
    size_t dbSize;
    double cmliteMsaEvalThr;
    double cmbuildEre;
    double cmbuildSymfrac;
    int cmbuildNoss;
    int localMin;
    int chainBand;
    int chainGap;
    int envMargin;
    int envMaxIter;
    int chainWord;
    int chainXdrop;
    int chainMinHsp;
    int chainBridge;
    int chainMaxHsp;
    int chainExtend;
    double chainExtendEval;
    std::string kmerAlph;

    PARAMETER(PARAM_CM_REGION)
    PARAMETER(PARAM_CM_MODE)
    PARAMETER(PARAM_CM_ALIGN)
    PARAMETER(PARAM_CM_ALIGN_BANDED)
    PARAMETER(PARAM_CM_LOCAL)
    PARAMETER(PARAM_DB_SIZE)
    PARAMETER(PARAM_CMLITE_MSA_EVAL)
    PARAMETER(PARAM_CMBUILD_ERE)
    PARAMETER(PARAM_CMBUILD_SYMFRAC)
    PARAMETER(PARAM_CMBUILD_NOSS)
    PARAMETER(PARAM_LOCAL_MIN)
    PARAMETER(PARAM_CHAIN_BAND)
    PARAMETER(PARAM_CHAIN_GAP)
    PARAMETER(PARAM_ENV_MARGIN)
    PARAMETER(PARAM_ENV_MAX_ITER)
    PARAMETER(PARAM_CHAIN_WORD)
    PARAMETER(PARAM_CHAIN_XDROP)
    PARAMETER(PARAM_CHAIN_MIN_HSP)
    PARAMETER(PARAM_CHAIN_BRIDGE)
    PARAMETER(PARAM_CHAIN_MAX_HSP)
    PARAMETER(PARAM_CHAIN_EXTEND)
    PARAMETER(PARAM_CHAIN_EXTEND_EVAL)
    PARAMETER(PARAM_KMER_ALPH)

    std::vector<MMseqsParameter*> splitstrand;
    std::vector<MMseqsParameter*> rnaalign;
    std::vector<MMseqsParameter*> cmbuild;
    std::vector<MMseqsParameter*> cmscan;
    std::vector<MMseqsParameter*> dinucdb;
    std::vector<MMseqsParameter*> envalign;
    std::vector<MMseqsParameter*> chainalign;
    std::vector<MMseqsParameter*> redkmermatcher;
};

#endif
