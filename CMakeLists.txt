cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
if(CMAKE_VERSION VERSION_LESS "4.0")
    cmake_policy(SET CMP0060 OLD)
endif()
cmake_policy(SET CMP0074 NEW)
project(riboseek C CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

if(NOT XCODE AND NOT MSVC AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)
target_compile_definitions(mmseqs-framework PUBLIC RIBOSEEK=1)

# Hack: apply static linker flags until the scope is fixed in mmseqs
if(FORCE_STATIC_DEPS)
    if(ENABLE_CUDA)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    endif()
endif()


set(RIBOSEEK_FRAMEWORK_ONLY 0 CACHE BOOL "Framework mode (don't create riboseek executable)")
if (RIBOSEEK_FRAMEWORK_ONLY)
    set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
endif()

# Infernal bridge for cmbuild
set(ENABLE_INFERNAL_BRIDGE 1 CACHE BOOL "Enable embedded Infernal-backed CM build bridge")
set(INFERNAL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/lib/infernal/vendor" CACHE PATH "Path to vendored Infernal source/build tree")
if (ENABLE_INFERNAL_BRIDGE)
    if (NOT EXISTS ${INFERNAL_ROOT}/src/cmbuild.c)
        message(FATAL_ERROR "ENABLE_INFERNAL_BRIDGE=ON requires INFERNAL_ROOT with src/cmbuild.c")
    endif ()
    add_subdirectory(lib/infernal)
endif ()

# Optional ViennaRNA reference library for development parity checks.
# Production riboseek RNA folding uses src/rnautils/RnaFold.cpp.
option(ENABLE_VIENNARNA_REFERENCE "Build vendored ViennaRNA reference target" OFF)
if (ENABLE_VIENNARNA_REFERENCE)
    add_subdirectory(lib/viennarna)
endif()

add_subdirectory(data)
add_subdirectory(src)
