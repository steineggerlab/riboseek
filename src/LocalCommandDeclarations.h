#ifndef RIBOSEEK_LOCALCOMMANDDECLARATIONS_H
#define RIBOSEEK_LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int splitstrand(int argc, const char **argv, const Command& command);
extern int rnaalign(int argc, const char **argv, const Command& command);
extern int rnasearch(int argc, const char **argv, const Command& command);
extern int generatecm(int argc, const char **argv, const Command& command);
extern int cmbuild(int argc, const char **argv, const Command& command);
extern int cmsearch(int argc, const char **argv, const Command& command);
extern int cmscan(int argc, const char **argv, const Command& command);

#endif
