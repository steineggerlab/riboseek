// WebAssembly entry point: expose the folder to JavaScript via embind.
// Built with emcc (see `make wasm`); the NEON block in tornadofold.h is inactive
// on wasm, so this is the exact-int32 scalar path.
#include "tornadofold.h"
#include <emscripten/bind.h>
#include <cctype>
#include <string>

// Fold a sequence and return "<dot-bracket>\t<mfe centikcal>". The input is
// sanitized to match the CLI: whitespace dropped, lowercase upper-cased, T->U,
// anything else treated as N (kept, so positions still line up).
static std::string foldSeq(const std::string& seq) {
    std::string s;
    for (char c : seq) {
        if (std::isspace((unsigned char)c)) {
            continue;
        }
        char u = (char)std::toupper((unsigned char)c);
        if (u == 'T') {
            u = 'U';
        }
        s += u;
    }
    if (s.empty()) {
        return std::string("\t0");
    }
    tornadofold::TornadoFold f;
    int e = f.fold(s);
    std::string db = f.traceback(e);
    return db + "\t" + std::to_string(e);
}

EMSCRIPTEN_BINDINGS(tornadofold) {
    emscripten::function("foldSeq", &foldSeq);
}
