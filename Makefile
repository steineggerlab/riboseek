CXX ?= clang++
ifneq (,$(filter x86_64 amd64,$(shell uname -m)))
SIMDFLAGS ?= -mavx2
endif
CXXFLAGS ?= -O3 -std=c++11 -Wall $(SIMDFLAGS)
# OpenMP: plain -fopenmp on Linux (clang/gcc); Apple clang needs the libomp form.
ifeq ($(shell uname -s),Darwin)
OMPFLAGS ?= -Xpreprocessor -fopenmp -lomp
else
OMPFLAGS ?= -fopenmp
endif

HDRS = src/ribossfold.h src/energy.h src/t2004.h

all: ribossfold verify

ribossfold: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

# OpenMP batch build (needs libomp: brew install libomp)
ribossfold_omp: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ src/main.cpp

# NEON disabled, isolates the scalar constant factor
ribossfold_noneon: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -DDISABLE_NEON -o $@ src/main.cpp

verify: src/verify.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/verify.cpp

EMCC ?= em++
EMFLAGS ?= -O3 -std=c++17 -Isrc -lembind \
	-sMODULARIZE=1 -sEXPORT_NAME=ribossfold -sENVIRONMENT=web,worker \
	-sINITIAL_MEMORY=268435456

wasm: web/ribossfold.js web/ribossfold-simd.js
web/ribossfold.js: web/ribossfold_wasm.cpp $(HDRS)
	$(EMCC) $(EMFLAGS) -o $@ web/ribossfold_wasm.cpp
web/ribossfold-simd.js: web/ribossfold_wasm.cpp $(HDRS)
	$(EMCC) $(EMFLAGS) -msimd128 -o $@ web/ribossfold_wasm.cpp


NPMFLAGS ?= $(EMFLAGS) -sSINGLE_FILE=1

npm: npm/ribossfold.js npm/ribossfold-simd.js npm/LICENSE npm/README.md
npm/ribossfold.js: web/ribossfold_wasm.cpp $(HDRS)
	$(EMCC) $(NPMFLAGS) -o $@ web/ribossfold_wasm.cpp
npm/ribossfold-simd.js: web/ribossfold_wasm.cpp $(HDRS)
	$(EMCC) $(NPMFLAGS) -msimd128 -o $@ web/ribossfold_wasm.cpp
npm/LICENSE: LICENSE
	cp LICENSE $@
npm/README.md: README.md
	cp README.md $@

clean:
	rm -f ribossfold ribossfold_omp ribossfold_noneon verify \
	  web/ribossfold.js web/ribossfold.wasm web/ribossfold-simd.js web/ribossfold-simd.wasm \
	  npm/ribossfold.js npm/ribossfold-simd.js npm/LICENSE npm/README.md npm/*.tgz

.PHONY: all clean wasm npm
