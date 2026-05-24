cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
if(CMAKE_VERSION VERSION_LESS "4.0")
    cmake_policy(SET CMP0060 OLD)
endif()
cmake_policy(SET CMP0074 NEW)
project(riboseek C CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")
set(RIBOSEEK_LINEAR_TURBOFOLD_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/lib/LinearTurboFold")
set(RIBOSEEK_LINEAR_TURBOFOLD_EXPECTED_COMMIT "445808976857df01c57e5d9302a0bbee10515bb5")
set(RIBOSEEK_LINEAR_TURBOFOLD_PATCH
    "${CMAKE_CURRENT_SOURCE_DIR}/patches/linear_turbofold_linearpartition_optimizations.patch")
if (EXISTS "${RIBOSEEK_LINEAR_TURBOFOLD_ROOT}/src/LinearPartition/src/LinearPartition.cpp"
    AND EXISTS "${RIBOSEEK_LINEAR_TURBOFOLD_PATCH}")
    find_package(Git QUIET)
    if (Git_FOUND)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${RIBOSEEK_LINEAR_TURBOFOLD_ROOT}" rev-parse HEAD
            OUTPUT_VARIABLE RIBOSEEK_LINEAR_TURBOFOLD_HEAD
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${RIBOSEEK_LINEAR_TURBOFOLD_ROOT}" status --short
            OUTPUT_VARIABLE RIBOSEEK_LINEAR_TURBOFOLD_STATUS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${RIBOSEEK_LINEAR_TURBOFOLD_ROOT}"
                    apply --check "${RIBOSEEK_LINEAR_TURBOFOLD_PATCH}"
            RESULT_VARIABLE RIBOSEEK_LINEAR_TURBOFOLD_PATCH_CAN_APPLY
            OUTPUT_QUIET
            ERROR_QUIET)
        if (RIBOSEEK_LINEAR_TURBOFOLD_PATCH_CAN_APPLY EQUAL 0)
            message(STATUS "Applying riboseek LinearPartition optimization patch")
            execute_process(
                COMMAND "${GIT_EXECUTABLE}" -C "${RIBOSEEK_LINEAR_TURBOFOLD_ROOT}"
                        apply "${RIBOSEEK_LINEAR_TURBOFOLD_PATCH}"
                RESULT_VARIABLE RIBOSEEK_LINEAR_TURBOFOLD_PATCH_APPLY_STATUS)
            if (NOT RIBOSEEK_LINEAR_TURBOFOLD_PATCH_APPLY_STATUS EQUAL 0)
                message(FATAL_ERROR "Failed to apply LinearPartition optimization patch")
            endif ()
        else ()
            execute_process(
                COMMAND "${GIT_EXECUTABLE}" -C "${RIBOSEEK_LINEAR_TURBOFOLD_ROOT}"
                        apply --reverse --check "${RIBOSEEK_LINEAR_TURBOFOLD_PATCH}"
                RESULT_VARIABLE RIBOSEEK_LINEAR_TURBOFOLD_PATCH_REVERSES
                OUTPUT_QUIET
                ERROR_QUIET)
            if (NOT RIBOSEEK_LINEAR_TURBOFOLD_PATCH_REVERSES EQUAL 0)
                message(FATAL_ERROR
                        "LinearPartition optimization patch is neither applicable nor already applied.\n"
                        "LinearTurboFold HEAD: ${RIBOSEEK_LINEAR_TURBOFOLD_HEAD}\n"
                        "Expected HEAD: ${RIBOSEEK_LINEAR_TURBOFOLD_EXPECTED_COMMIT}\n"
                        "LinearTurboFold status:\n${RIBOSEEK_LINEAR_TURBOFOLD_STATUS}\n"
                        "Fix: run `git submodule update --init --recursive lib/LinearTurboFold` from the riboseek repository root, then rerun CMake. If the submodule has local changes, stash or discard them first.")
            endif ()
        endif ()
    else ()
        message(FATAL_ERROR "Git is required to apply the LinearPartition optimization patch")
    endif ()
endif ()

if(NOT XCODE AND NOT MSVC AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)
target_compile_definitions(mmseqs-framework PUBLIC RIBOSEEK=1)

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

# ViennaRNA minimal library (hand-rolled CMake wrapper, autotools-free)
add_subdirectory(lib/viennarna)

add_subdirectory(data)
add_subdirectory(src)
