#include <cstddef>
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
void init() {
    registerCommands(&baseCommands);
    registerCommands(&riboseekCommands);
    registerDinucleotideMapping();
}
void (*initCommands)(void) = init;
void initParameterSingleton() { new LocalParameters; }
