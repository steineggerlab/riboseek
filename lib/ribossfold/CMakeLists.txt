cmake_minimum_required(VERSION 3.14)
project(ribossfold LANGUAGES CXX)

# ribossfold is a header-only Turner 2004 SIMD MFE folder.
#
#   add_subdirectory(path/to/rna EXCLUDE_FROM_ALL)   # headers only, no binaries
#   target_link_libraries(myapp PRIVATE ribossfold::ribossfold)
#
# EXCLUDE_FROM_ALL keeps the CLI drivers out of the parent's default build; the
# INTERFACE library has no build artifact of its own, so nothing compiles until
# a consumer actually #includes and uses the headers.

add_library(ribossfold INTERFACE)
add_library(ribossfold::ribossfold ALIAS ribossfold)
target_include_directories(ribossfold INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>)
target_compile_features(ribossfold INTERFACE cxx_std_11)

# Build the CLI drivers only when ribossfold is the top-level project. Included via
# add_subdirectory() it stays header-only regardless of EXCLUDE_FROM_ALL; set
# -DRIBOSSFOLD_BUILD_TOOLS=ON to build them as a subproject anyway.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_ribossfold_tools_default ON)
else()
    set(_ribossfold_tools_default OFF)
endif()
option(RIBOSSFOLD_BUILD_TOOLS "Build the ribossfold CLI and verify binaries" ${_ribossfold_tools_default})

if(RIBOSSFOLD_BUILD_TOOLS)
    add_executable(ribossfold_cli src/main.cpp)
    set_target_properties(ribossfold_cli PROPERTIES OUTPUT_NAME ribossfold)
    target_link_libraries(ribossfold_cli PRIVATE ribossfold::ribossfold)

    add_executable(ribossfold_verify src/verify.cpp)
    set_target_properties(ribossfold_verify PROPERTIES OUTPUT_NAME verify)
    target_link_libraries(ribossfold_verify PRIVATE ribossfold::ribossfold)

    # NEON-disabled build (scalar constant-factor baseline).
    add_executable(ribossfold_noneon src/main.cpp)
    target_compile_definitions(ribossfold_noneon PRIVATE DISABLE_NEON)
    target_link_libraries(ribossfold_noneon PRIVATE ribossfold::ribossfold)

    # Optional OpenMP batch build (parallel across sequences).
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        add_executable(ribossfold_omp src/main.cpp)
        target_link_libraries(ribossfold_omp PRIVATE ribossfold::ribossfold OpenMP::OpenMP_CXX)
    endif()
endif()
