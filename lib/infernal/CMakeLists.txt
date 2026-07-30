set(INFERNAL_SRC_DIR ${INFERNAL_ROOT}/src)
set(INFERNAL_HMMER_DIR ${INFERNAL_ROOT}/hmmer)
set(INFERNAL_HMMER_SRC_DIR ${INFERNAL_ROOT}/hmmer/src)
set(INFERNAL_EASEL_DIR ${INFERNAL_ROOT}/easel)

if (NOT EXISTS ${INFERNAL_SRC_DIR}/config.h)
    message(FATAL_ERROR "vendored infernal requires ${INFERNAL_SRC_DIR}/config.h")
endif ()
if (NOT EXISTS ${INFERNAL_CONFIG_DIR}/p7_config.h OR NOT EXISTS ${INFERNAL_CONFIG_DIR}/esl_config.h)
    message(FATAL_ERROR "missing config headers in ${INFERNAL_CONFIG_DIR}")
endif ()

# Build Infernal/HMMER/Easel with the same SIMD/arch flags as riboseek (MMSEQS_ARCH).
if (MMSEQS_ARCH)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MMSEQS_ARCH}")
endif ()

set(EASEL_SRCS
        ${INFERNAL_EASEL_DIR}/easel.c
        # ${INFERNAL_EASEL_DIR}/esl_alloc.c
        ${INFERNAL_EASEL_DIR}/esl_alphabet.c
        ${INFERNAL_EASEL_DIR}/esl_arr2.c
        ${INFERNAL_EASEL_DIR}/esl_arr3.c
        ${INFERNAL_EASEL_DIR}/esl_bitfield.c
        ${INFERNAL_EASEL_DIR}/esl_buffer.c
        ${INFERNAL_EASEL_DIR}/esl_cluster.c
        # ${INFERNAL_EASEL_DIR}/esl_composition.c
        # ${INFERNAL_EASEL_DIR}/esl_cpu.c
        ${INFERNAL_EASEL_DIR}/esl_dirichlet.c
        ${INFERNAL_EASEL_DIR}/esl_distance.c
        ${INFERNAL_EASEL_DIR}/esl_dmatrix.c
        # ${INFERNAL_EASEL_DIR}/esl_dsqdata.c
        ${INFERNAL_EASEL_DIR}/esl_exponential.c
        ${INFERNAL_EASEL_DIR}/esl_fileparser.c
        # ${INFERNAL_EASEL_DIR}/esl_gamma.c
        # ${INFERNAL_EASEL_DIR}/esl_gencode.c
        ${INFERNAL_EASEL_DIR}/esl_getopts.c
        # ${INFERNAL_EASEL_DIR}/esl_gev.c
        ${INFERNAL_EASEL_DIR}/esl_graph.c
        ${INFERNAL_EASEL_DIR}/esl_gumbel.c
        # ${INFERNAL_EASEL_DIR}/esl_heap.c
        # ${INFERNAL_EASEL_DIR}/esl_histogram.c
        ${INFERNAL_EASEL_DIR}/esl_hmm.c
        # ${INFERNAL_EASEL_DIR}/esl_huffman.c
        # ${INFERNAL_EASEL_DIR}/esl_hyperexp.c
        # ${INFERNAL_EASEL_DIR}/esl_iset.c
        # ${INFERNAL_EASEL_DIR}/esl_json.c
        ${INFERNAL_EASEL_DIR}/esl_keyhash.c
        # ${INFERNAL_EASEL_DIR}/esl_lognormal.c
        ${INFERNAL_EASEL_DIR}/esl_matrixops.c
        ${INFERNAL_EASEL_DIR}/esl_mem.c
        ${INFERNAL_EASEL_DIR}/esl_minimizer.c
        ${INFERNAL_EASEL_DIR}/esl_mixdchlet.c
        # ${INFERNAL_EASEL_DIR}/esl_mixgev.c
        # ${INFERNAL_EASEL_DIR}/esl_mpi.c
        ${INFERNAL_EASEL_DIR}/esl_msa.c
        ${INFERNAL_EASEL_DIR}/esl_msacluster.c
        ${INFERNAL_EASEL_DIR}/esl_msafile.c
        # ${INFERNAL_EASEL_DIR}/esl_msafile2.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_a2m.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_afa.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_clustal.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_phylip.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_psiblast.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_selex.c
        ${INFERNAL_EASEL_DIR}/esl_msafile_stockholm.c
        # ${INFERNAL_EASEL_DIR}/esl_msashuffle.c
        ${INFERNAL_EASEL_DIR}/esl_msaweight.c
        # ${INFERNAL_EASEL_DIR}/esl_normal.c
        # ${INFERNAL_EASEL_DIR}/esl_paml.c
        ${INFERNAL_EASEL_DIR}/esl_quicksort.c
        ${INFERNAL_EASEL_DIR}/esl_random.c
        ${INFERNAL_EASEL_DIR}/esl_rand64.c
        ${INFERNAL_EASEL_DIR}/esl_randomseq.c
        # ${INFERNAL_EASEL_DIR}/esl_ratematrix.c
        # ${INFERNAL_EASEL_DIR}/esl_recorder.c
        # ${INFERNAL_EASEL_DIR}/esl_red_black.c
        # ${INFERNAL_EASEL_DIR}/esl_regexp.c
        ${INFERNAL_EASEL_DIR}/esl_rootfinder.c
        # ${INFERNAL_EASEL_DIR}/esl_scorematrix.c
        ${INFERNAL_EASEL_DIR}/esl_sq.c
        ${INFERNAL_EASEL_DIR}/esl_sqio.c
        ${INFERNAL_EASEL_DIR}/esl_sqio_ascii.c
        ${INFERNAL_EASEL_DIR}/esl_sqio_ncbi.c
        ${INFERNAL_EASEL_DIR}/esl_ssi.c
        ${INFERNAL_EASEL_DIR}/esl_stack.c
        ${INFERNAL_EASEL_DIR}/esl_stats.c
        ${INFERNAL_EASEL_DIR}/esl_stopwatch.c
        # ${INFERNAL_EASEL_DIR}/esl_stretchexp.c
        # ${INFERNAL_EASEL_DIR}/esl_subcmd.c
        # ${INFERNAL_EASEL_DIR}/esl_threads.c
        ${INFERNAL_EASEL_DIR}/esl_tree.c
        # ${INFERNAL_EASEL_DIR}/esl_varint.c
        ${INFERNAL_EASEL_DIR}/esl_vectorops.c
        # ${INFERNAL_EASEL_DIR}/esl_weibull.c
        # ${INFERNAL_EASEL_DIR}/esl_workqueue.c
        ${INFERNAL_EASEL_DIR}/esl_wuss.c
)

set(HMMER_SRCS
        # ${INFERNAL_HMMER_SRC_DIR}/build.c
        # ${INFERNAL_HMMER_SRC_DIR}/cachedb.c
        # ${INFERNAL_HMMER_SRC_DIR}/cachedb_shard.c
        ${INFERNAL_HMMER_SRC_DIR}/emit.c
        ${INFERNAL_HMMER_SRC_DIR}/errors.c
        ${INFERNAL_HMMER_SRC_DIR}/evalues.c
        ${INFERNAL_HMMER_SRC_DIR}/eweight.c
        ${INFERNAL_HMMER_SRC_DIR}/generic_decoding.c
        ${INFERNAL_HMMER_SRC_DIR}/generic_fwdback.c
        # ${INFERNAL_HMMER_SRC_DIR}/generic_fwdback_chk.c
        # ${INFERNAL_HMMER_SRC_DIR}/generic_fwdback_banded.c
        ${INFERNAL_HMMER_SRC_DIR}/generic_null2.c
        # ${INFERNAL_HMMER_SRC_DIR}/generic_msv.c
        ${INFERNAL_HMMER_SRC_DIR}/generic_optacc.c
        ${INFERNAL_HMMER_SRC_DIR}/generic_stotrace.c
        # ${INFERNAL_HMMER_SRC_DIR}/generic_viterbi.c
        # ${INFERNAL_HMMER_SRC_DIR}/generic_vtrace.c
        # ${INFERNAL_HMMER_SRC_DIR}/h2_io.c
        # ${INFERNAL_HMMER_SRC_DIR}/heatmap.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmlogo.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmdmstr.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmdmstr_shard.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmd_search_status.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmdwrkr.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmdwrkr_shard.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmdutils.c
        ${INFERNAL_HMMER_SRC_DIR}/hmmer.c
        ${INFERNAL_HMMER_SRC_DIR}/logsum.c
        ${INFERNAL_HMMER_SRC_DIR}/modelconfig.c
        ${INFERNAL_HMMER_SRC_DIR}/modelstats.c
        # ${INFERNAL_HMMER_SRC_DIR}/mpisupport.c
        # ${INFERNAL_HMMER_SRC_DIR}/seqmodel.c
        ${INFERNAL_HMMER_SRC_DIR}/tracealign.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_alidisplay.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_bg.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_builder.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_domain.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_domaindef.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_gbands.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_gmx.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_gmxb.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_gmxchk.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_hit.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_hmm.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_hmmcache.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_hmmd_search_stats.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_hmmfile.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_hmmwindow.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_pipeline.c
        # ${INFERNAL_HMMER_SRC_DIR}/p7_prior.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_profile.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_spensemble.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_tophits.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_trace.c
        ${INFERNAL_HMMER_SRC_DIR}/p7_scoredata.c
        # ${INFERNAL_HMMER_SRC_DIR}/hmmpgmd2msa.c
        ${INFERNAL_HMMER_SRC_DIR}/fm_alphabet.c
        ${INFERNAL_HMMER_SRC_DIR}/fm_general.c
        ${INFERNAL_HMMER_SRC_DIR}/fm_sse.c
        ${INFERNAL_HMMER_SRC_DIR}/fm_ssv.c
)

if (ARM)
    set(HMMER_IMPL_DIR ${INFERNAL_HMMER_SRC_DIR}/impl_neon)
elseif (X86 OR X64)
    set(HMMER_IMPL_DIR ${INFERNAL_HMMER_SRC_DIR}/impl_sse)
else ()
    message(FATAL_ERROR "Unsupported architecture: neither ARM nor X86/X64 is set")
endif ()

set(HMMER_IMPL_SRCS
        ${HMMER_IMPL_DIR}/decoding.c
        ${HMMER_IMPL_DIR}/fwdback.c
        ${HMMER_IMPL_DIR}/io.c
        ${HMMER_IMPL_DIR}/ssvfilter.c
        ${HMMER_IMPL_DIR}/msvfilter.c
        ${HMMER_IMPL_DIR}/null2.c
        ${HMMER_IMPL_DIR}/optacc.c
        ${HMMER_IMPL_DIR}/stotrace.c
        ${HMMER_IMPL_DIR}/vitfilter.c
        ${HMMER_IMPL_DIR}/p7_omx.c
        ${HMMER_IMPL_DIR}/p7_oprofile.c
        # ${HMMER_IMPL_DIR}/mpi.c
)
if (ARM)
    list(APPEND EASEL_SRCS ${INFERNAL_EASEL_DIR}/esl_neon.c)
elseif (X86 OR X64)
    list(APPEND EASEL_SRCS ${INFERNAL_EASEL_DIR}/esl_sse.c)
endif ()

set(INFERNAL_CORE_SRCS
        ${INFERNAL_SRC_DIR}/cm.c
        ${INFERNAL_SRC_DIR}/cm_alidisplay.c
        ${INFERNAL_SRC_DIR}/cm_alndata.c
        ${INFERNAL_SRC_DIR}/cm_dpalign.c
        ${INFERNAL_SRC_DIR}/cm_dpalign_trunc.c
        ${INFERNAL_SRC_DIR}/cm_dpsearch.c
        ${INFERNAL_SRC_DIR}/cm_dpsearch_trunc.c
        ${INFERNAL_SRC_DIR}/cm_dpsmall.c
        ${INFERNAL_SRC_DIR}/cm_file.c
        ${INFERNAL_SRC_DIR}/cm_modelconfig.c
        ${INFERNAL_SRC_DIR}/cm_modelmaker.c
        ${INFERNAL_SRC_DIR}/cm_mx.c
        ${INFERNAL_SRC_DIR}/cm_parsetree.c
        ${INFERNAL_SRC_DIR}/cm_pipeline.c
        ${INFERNAL_SRC_DIR}/cm_qdband.c
        ${INFERNAL_SRC_DIR}/cm_submodel.c
        ${INFERNAL_SRC_DIR}/cm_tophits.c
        ${INFERNAL_SRC_DIR}/cm_trunc.c
        # ${INFERNAL_SRC_DIR}/cm_p7_band.c
        ${INFERNAL_SRC_DIR}/cm_p7_domaindef.c
        ${INFERNAL_SRC_DIR}/cm_p7_modelconfig_trunc.c
        ${INFERNAL_SRC_DIR}/cm_p7_modelmaker.c
        ${INFERNAL_SRC_DIR}/cp9.c
        ${INFERNAL_SRC_DIR}/cp9_dp.c
        ${INFERNAL_SRC_DIR}/cp9_modelmaker.c
        ${INFERNAL_SRC_DIR}/cp9_mx.c
        ${INFERNAL_SRC_DIR}/cp9_trace.c
        ${INFERNAL_SRC_DIR}/alphabet.c
        ${INFERNAL_SRC_DIR}/display.c
        ${INFERNAL_SRC_DIR}/errors.c
        ${INFERNAL_SRC_DIR}/eweight.c
        ${INFERNAL_SRC_DIR}/hmmband.c
        ${INFERNAL_SRC_DIR}/logsum.c
        # ${INFERNAL_SRC_DIR}/mpisupport.c
        ${INFERNAL_SRC_DIR}/prior.c
        ${INFERNAL_SRC_DIR}/rnamat.c
        ${INFERNAL_SRC_DIR}/stats.c
        ${INFERNAL_SRC_DIR}/truncyk.c
)

add_library(infernal-vendor-easel STATIC ${EASEL_SRCS})
target_compile_definitions(infernal-vendor-easel PRIVATE HAVE_CONFIG_H=1)
target_include_directories(infernal-vendor-easel PRIVATE
        ${INFERNAL_CONFIG_DIR}
        ${INFERNAL_EASEL_DIR}
)
target_link_libraries(infernal-vendor-easel PUBLIC m)

add_library(infernal-vendor-hmmer STATIC ${HMMER_SRCS} ${HMMER_IMPL_SRCS})
target_compile_definitions(infernal-vendor-hmmer PRIVATE HAVE_CONFIG_H=1)
target_include_directories(infernal-vendor-hmmer PRIVATE
        ${INFERNAL_CONFIG_DIR}
        ${INFERNAL_HMMER_SRC_DIR}
        ${HMMER_IMPL_DIR}
        ${INFERNAL_EASEL_DIR}
)
target_link_libraries(infernal-vendor-hmmer PUBLIC infernal-vendor-easel m)

add_library(infernal-vendor-core STATIC ${INFERNAL_CORE_SRCS})
target_compile_definitions(infernal-vendor-core PRIVATE HAVE_CONFIG_H=1)
target_include_directories(infernal-vendor-core PRIVATE
        ${INFERNAL_CONFIG_DIR}
        ${INFERNAL_SRC_DIR}
        ${INFERNAL_HMMER_SRC_DIR}
        ${INFERNAL_EASEL_DIR}
)
target_link_libraries(infernal-vendor-core PUBLIC infernal-vendor-hmmer infernal-vendor-easel m)
