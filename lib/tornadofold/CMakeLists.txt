cmake_minimum_required(VERSION 3.14)
project(tornadofold LANGUAGES CXX)

# tornadofold is a header-only Turner 2004 SIMD MFE folder (everything lives in
# src/*.h). Consumers embed it with:
#
#   add_subdirectory(path/to/rna EXCLUDE_FROM_ALL)   # headers only, no binaries
#   target_link_libraries(myapp PRIVATE tornadofold::tornadofold)
#
# EXCLUDE_FROM_ALL keeps the CLI drivers out of the parent's default build; the
# INTERFACE library has no build artifact of its own, so nothing compiles until
# a consumer actually #includes and uses the headers.

add_library(tornadofold INTERFACE)
add_library(tornadofold::tornadofold ALIAS tornadofold)
target_include_directories(tornadofold INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>)
target_compile_features(tornadofold INTERFACE cxx_std_11)

# Build the CLI drivers only when tornadofold is the top-level project. Included via
# add_subdirectory() it stays header-only regardless of EXCLUDE_FROM_ALL; set
# -DTORNADOFOLD_BUILD_TOOLS=ON to build them as a subproject anyway.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_tornadofold_tools_default ON)
else()
    set(_tornadofold_tools_default OFF)
endif()
option(TORNADOFOLD_BUILD_TOOLS "Build the tornadofold CLI and verify binaries" ${_tornadofold_tools_default})

if(TORNADOFOLD_BUILD_TOOLS)
    add_executable(tornadofold_cli src/main.cpp)
    set_target_properties(tornadofold_cli PROPERTIES OUTPUT_NAME tornadofold)
    target_link_libraries(tornadofold_cli PRIVATE tornadofold::tornadofold)

    add_executable(tornadofold_verify src/verify.cpp)
    set_target_properties(tornadofold_verify PROPERTIES OUTPUT_NAME verify)
    target_link_libraries(tornadofold_verify PRIVATE tornadofold::tornadofold)

    # NEON-disabled build (scalar constant-factor baseline).
    add_executable(tornadofold_noneon src/main.cpp)
    target_compile_definitions(tornadofold_noneon PRIVATE DISABLE_NEON)
    target_link_libraries(tornadofold_noneon PRIVATE tornadofold::tornadofold)

    # Optional OpenMP batch build (parallel across sequences).
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        add_executable(tornadofold_omp src/main.cpp)
        target_link_libraries(tornadofold_omp PRIVATE tornadofold::tornadofold OpenMP::OpenMP_CXX)
    endif()
endif()
