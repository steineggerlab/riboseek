include_directories(.)
include_directories(commons)
add_subdirectory(commons)
add_subdirectory(rnautils)
add_subdirectory(workflow)
add_subdirectory(version)

add_library(riboseek-framework
    ${commons_source_files}
    ${rnautils_source_files}
    ${workflow_source_files}
    LocalCommandDeclarations.h
    LocalParameters.h
    LocalParameters.cpp
    RiboseekBase.cpp
)
mmseqs_setup_derived_target(riboseek-framework)
# Put riboseek's generated/ BEFORE mmseqs's so our blastdigp.sh.h (etc.) shadows
# the stock MMseqs2 copy when RnaSearch.cpp does #include "blastdigp.sh.h".
target_include_directories(riboseek-framework BEFORE PUBLIC ${CMAKE_BINARY_DIR}/generated)
add_dependencies(riboseek-framework generated_riboseek)
target_include_directories(riboseek-framework PUBLIC ${CMAKE_SOURCE_DIR}/lib)

# Infernal bridge
if (TARGET infernal-bridge)
    target_link_libraries(riboseek-framework infernal-bridge)
    target_compile_definitions(riboseek-framework PUBLIC -DHAVE_INFERNAL_BRIDGE=1)
endif()


if (NOT RIBOSEEK_FRAMEWORK_ONLY)
    add_executable(riboseek riboseek.cpp)
    mmseqs_setup_derived_target(riboseek riboseek-framework)
    target_link_libraries(mmseqs-framework version)
    target_link_libraries(riboseek version)
    install(TARGETS riboseek DESTINATION bin)
endif()
