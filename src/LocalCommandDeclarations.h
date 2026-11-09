#ifndef RIBOSEEK_LOCALCOMMANDDECLARATIONS_H
#define RIBOSEEK_LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int splitstrand(int argc, const char **argv, const Command& command);
extern int rnaalign(int argc, const char **argv, const Command& command);
extern int riboseekSearch(int argc, const char **argv, const Command& command);
extern int riboseekCluster(int argc, const char **argv, const Command& command);
extern int riboseekLinclust(int argc, const char **argv, const Command& command);
extern int envalign(int argc, const char **argv, const Command& command);
extern int xdropalign(int argc, const char **argv, const Command &command);
extern int chainalign(int argc, const char **argv, const Command &command);
extern int redkmermatcher(int argc, const char **argv, const Command& command);
extern int riboseekKmermatcher(int argc, const char **argv, const Command& command);
extern int dinucdb(int argc, const char **argv, const Command& command);
extern int generatecm(int argc, const char **argv, const Command& command);
extern int cmbuild(int argc, const char **argv, const Command& command);
extern int cmsearch(int argc, const char **argv, const Command& command);
extern int cmscan(int argc, const char **argv, const Command& command);

#endif
