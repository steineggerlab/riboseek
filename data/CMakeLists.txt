include(MMseqsResourceCompiler)

set(GENERATED_MATRICES
        dinuc.out
        )

set(GENERATED_WORKFLOWS
        blastdi.sh
        blastdigp.sh
        blastsplitdi.sh
        blastsplitdigp.sh
        blastsplitdigp2.sh
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${GENERATED_MATRICES} ${GENERATED_WORKFLOWS})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(generated_riboseek ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
