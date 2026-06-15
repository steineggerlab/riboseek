#include <cstddef>
#include <cstring>
#include "Command.h"
#include "LocalParameters.h"
#include "DinucleotideMapping.h"

const char* binary_name = "riboseek";
const char* tool_name = "riboseek";
const char* tool_introduction = "Riboseek enables fast RNA homology search using covariance models and structure-aware alignment.\n";
const char* main_author = "";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
const char* index_version_compatible = "rs1";
bool hide_base_commands = false;
bool hide_base_downloads = false;

extern std::vector<Command> baseCommands;
extern std::vector<Command> riboseekCommands;
static void removeBaseCommand(const char *name) {
    std::vector<Command> filtered;
    filtered.reserve(baseCommands.size());
    for (std::vector<Command>::const_iterator it = baseCommands.begin(); it != baseCommands.end(); ++it) {
        if (std::strcmp(it->cmd, name) != 0) {
            filtered.push_back(*it);
        }
    }
    baseCommands.swap(filtered);
}

void init() {
    removeBaseCommand("search");
    registerCommands(&baseCommands);
    registerCommands(&riboseekCommands);
    registerDinucleotideMapping();
    #ifdef HAVE_CUDA
    registerDinucleotideFilterConfig();
    #endif
}
void (*initCommands)(void) = init;
void initParameterSingleton() { new LocalParameters; }
