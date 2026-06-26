#ifndef RIBOSEEK_RNAFOLD_H
#define RIBOSEEK_RNAFOLD_H

#include <string>

namespace RiboseekRnaFold {

bool predictDotBracket(const std::string &rna, std::string &dotBracket, double *scoreOut = 0);

}

#endif
