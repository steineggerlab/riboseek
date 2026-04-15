#ifndef MMSEQS_INFERNALBRIDGE_H
#define MMSEQS_INFERNALBRIDGE_H

#include <string>

namespace InfernalBridge {

bool isConfigured();
bool buildCmFromStockholmText(const std::string &stockholmText, std::string &cmText, std::string &error);

} // namespace InfernalBridge

#endif
