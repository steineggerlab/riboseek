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

HDRS = src/tornadofold.h src/energy.h src/t2004.h

all: tornadofold verify

tornadofold: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

# OpenMP batch build (needs libomp: brew install libomp)
tornadofold_omp: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ src/main.cpp

# NEON disabled, isolates the scalar constant factor
tornadofold_noneon: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -DDISABLE_NEON -o $@ src/main.cpp

verify: src/verify.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/verify.cpp

EMCC ?= emcc
EMFLAGS ?= -O3 -std=c++17 -Isrc -lembind \
	-sMODULARIZE=1 -sEXPORT_NAME=tornadofold -sENVIRONMENT=web,worker \
	-sINITIAL_MEMORY=268435456

wasm: web/tornadofold.js web/tornadofold-simd.js
web/tornadofold.js: web/tornadofold_wasm.cpp $(HDRS)
	$(EMCC) $(EMFLAGS) -o $@ web/tornadofold_wasm.cpp
web/tornadofold-simd.js: web/tornadofold_wasm.cpp $(HDRS)
	$(EMCC) $(EMFLAGS) -msimd128 -o $@ web/tornadofold_wasm.cpp

clean:
	rm -f tornadofold tornadofold_omp tornadofold_noneon verify \
	  web/tornadofold.js web/tornadofold.wasm web/tornadofold-simd.js web/tornadofold-simd.wasm

.PHONY: all clean wasm
