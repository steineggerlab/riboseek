# ViennaRNA minimal static library for MFE folding
# Built from the ViennaRNA package source tree

# ViennaRNA source: look relative to the riboseek root, then fall back to env/cache
set(VRNA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ViennaRNA/src/ViennaRNA")

if(NOT EXISTS "${VRNA_SRC_DIR}/fold_compound.c")
    message(WARNING "ViennaRNA source not found at ${VRNA_SRC_DIR}, skipping viennarna library")
    return()
endif()

# Core sources (from Makefile.am libRNA_conv_la_SOURCES)
set(VRNA_CONV_SOURCES
    ${VRNA_SRC_DIR}/fold_compound.c
    ${VRNA_SRC_DIR}/dist_vars.c
    ${VRNA_SRC_DIR}/treedist.c
    ${VRNA_SRC_DIR}/ProfileDist.c
    ${VRNA_SRC_DIR}/RNAstruct.c
    ${VRNA_SRC_DIR}/stringdist.c
    ${VRNA_SRC_DIR}/ProfileAln.c
    ${VRNA_SRC_DIR}/duplex.c
    ${VRNA_SRC_DIR}/mm.c
    ${VRNA_SRC_DIR}/2Dfold.c
    ${VRNA_SRC_DIR}/2Dpfold.c
    ${VRNA_SRC_DIR}/plex_functions.c
    ${VRNA_SRC_DIR}/ali_plex.c
    ${VRNA_SRC_DIR}/c_plex.c
    ${VRNA_SRC_DIR}/plex.c
    ${VRNA_SRC_DIR}/snofold.c
    ${VRNA_SRC_DIR}/snoop.c
    ${VRNA_SRC_DIR}/perturbation_fold.c
    ${VRNA_SRC_DIR}/model.c
    ${VRNA_SRC_DIR}/unstructured_domains.c
    ${VRNA_SRC_DIR}/heat_capacity.c
)

# JSON (bundled with ViennaRNA)
set(VRNA_JSON_SOURCES
    ${VRNA_SRC_DIR}/../json/json.c
)

# Cephes math library
set(VRNA_CEPHES_SOURCES
    ${VRNA_SRC_DIR}/../cephes/kn.c
    ${VRNA_SRC_DIR}/../cephes/const.c
    ${VRNA_SRC_DIR}/../cephes/mtherr.c
    ${VRNA_SRC_DIR}/../cephes/expn.c
)

# Eval
set(VRNA_EVAL_SOURCES
    ${VRNA_SRC_DIR}/eval/eval_exterior.c
    ${VRNA_SRC_DIR}/eval/eval_hairpin.c
    ${VRNA_SRC_DIR}/eval/eval_internal.c
    ${VRNA_SRC_DIR}/eval/eval_multibranch.c
    ${VRNA_SRC_DIR}/eval/eval_gquad.c
    ${VRNA_SRC_DIR}/eval/exp_eval_exterior.c
    ${VRNA_SRC_DIR}/eval/exp_eval_hairpin.c
    ${VRNA_SRC_DIR}/eval/exp_eval_internal.c
    ${VRNA_SRC_DIR}/eval/exp_eval_multibranch.c
    ${VRNA_SRC_DIR}/eval/eval_structures.c
    ${VRNA_SRC_DIR}/eval/eval_covariance.c
    ${VRNA_SRC_DIR}/eval/eval_wrappers.c
)

# Constraints
set(VRNA_CONSTRAINTS_SOURCES
    ${VRNA_SRC_DIR}/constraints/constraints.c
    ${VRNA_SRC_DIR}/constraints/hard.c
    ${VRNA_SRC_DIR}/constraints/soft.c
    ${VRNA_SRC_DIR}/constraints/soft_cb_multi.c
    ${VRNA_SRC_DIR}/constraints/sc_cb_mod.c
    ${VRNA_SRC_DIR}/constraints/sc_cb_mod_parser.c
    ${VRNA_SRC_DIR}/constraints/sc_cb_mod_wrappers.c
    ${VRNA_SRC_DIR}/constraints/ligand.c
)

# Probing
set(VRNA_PROBING_SOURCES
    ${VRNA_SRC_DIR}/probing/probing.c
    ${VRNA_SRC_DIR}/probing/strategy_chain.c
    ${VRNA_SRC_DIR}/probing/strategy_deigan.c
    ${VRNA_SRC_DIR}/probing/strategy_eddy.c
    ${VRNA_SRC_DIR}/probing/strategy_nlogp.c
    ${VRNA_SRC_DIR}/probing/strategy_zarringhalam.c
    ${VRNA_SRC_DIR}/probing/SHAPE.c
)

# Data
set(VRNA_DATA_SOURCES
    ${VRNA_SRC_DIR}/data/transform.c
)

# Math
set(VRNA_MATH_SOURCES
    ${VRNA_SRC_DIR}/math/fun_bin.c
    ${VRNA_SRC_DIR}/math/fun_chain.c
    ${VRNA_SRC_DIR}/math/fun_gaussian.c
    ${VRNA_SRC_DIR}/math/fun_kde.c
    ${VRNA_SRC_DIR}/math/fun_linear.c
    ${VRNA_SRC_DIR}/math/fun_log.c
    ${VRNA_SRC_DIR}/math/fun_logistic.c
)

# Structures
set(VRNA_STRUCTURES_SOURCES
    ${VRNA_SRC_DIR}/structures/structure_benchmark.c
    ${VRNA_SRC_DIR}/structures/centroid.c
    ${VRNA_SRC_DIR}/structures/structure_dotbracket.c
    ${VRNA_SRC_DIR}/structures/structure_helix.c
    ${VRNA_SRC_DIR}/structures/mea.c
    ${VRNA_SRC_DIR}/structures/structure_metrics.c
    ${VRNA_SRC_DIR}/structures/structure_pairtable.c
    ${VRNA_SRC_DIR}/structures/structure_problist.c
    ${VRNA_SRC_DIR}/structures/structure_shapes.c
    ${VRNA_SRC_DIR}/structures/structure_tree.c
    ${VRNA_SRC_DIR}/structures/structure_utils.c
)

# Sequences
set(VRNA_SEQUENCES_SOURCES
    ${VRNA_SRC_DIR}/sequences/alphabet.c
    ${VRNA_SRC_DIR}/sequences/sequence.c
    ${VRNA_SRC_DIR}/sequences/msa.c
    ${VRNA_SRC_DIR}/sequences/seq_utils.c
)

# Utils
set(VRNA_UTILS_SOURCES
    ${VRNA_SRC_DIR}/utils/utils.c
    ${VRNA_SRC_DIR}/utils/string_utils.c
    ${VRNA_SRC_DIR}/utils/higher_order_functions.c
    ${VRNA_SRC_DIR}/utils/cpu.c
    ${VRNA_SRC_DIR}/utils/units.c
    ${VRNA_SRC_DIR}/utils/log.c
    ${VRNA_SRC_DIR}/io/io_utils.c
    ${VRNA_SRC_DIR}/io/file_formats.c
    ${VRNA_SRC_DIR}/io/file_formats_msa.c
    ${VRNA_SRC_DIR}/search/BoyerMoore.c
    ${VRNA_SRC_DIR}/io/commands.c
)

# Plotting
set(VRNA_PLOTTING_SOURCES
    ${VRNA_SRC_DIR}/plotting/alignments.c
    ${VRNA_SRC_DIR}/plotting/layouts.c
    ${VRNA_SRC_DIR}/plotting/probabilities.c
    ${VRNA_SRC_DIR}/plotting/structures.c
    ${VRNA_SRC_DIR}/plotting/plot_utils.c
    ${VRNA_SRC_DIR}/plotting/svg.c
    ${VRNA_SRC_DIR}/plotting/eps.c
    ${VRNA_SRC_DIR}/plotting/gml.c
    ${VRNA_SRC_DIR}/plotting/xrna.c
    ${VRNA_SRC_DIR}/plotting/ssv.c
    ${VRNA_SRC_DIR}/plotting/RNApuzzler/RNApuzzler.c
    ${VRNA_SRC_DIR}/plotting/RNApuzzler/RNAturtle.c
)

# Naview layout (optional, always include)
set(VRNA_NAVIEW_SOURCES
    ${VRNA_SRC_DIR}/plotting/naview/naview.c
)

# MFE
set(VRNA_MFE_SOURCES
    ${VRNA_SRC_DIR}/mfe/mfe.c
    ${VRNA_SRC_DIR}/mfe/mfe_window.c
    ${VRNA_SRC_DIR}/mfe/mfe_wrappers.c
    ${VRNA_SRC_DIR}/mfe/mfe_window_wrappers.c
    ${VRNA_SRC_DIR}/mfe/mfe_exterior.c
    ${VRNA_SRC_DIR}/mfe/mfe_exterior_window.c
    ${VRNA_SRC_DIR}/mfe/mfe_internal.c
    ${VRNA_SRC_DIR}/mfe/mfe_multibranch.c
    ${VRNA_SRC_DIR}/mfe/mfe_gquad.c
    ${VRNA_SRC_DIR}/mfe/alifold.c
    ${VRNA_SRC_DIR}/mfe/cofold.c
    ${VRNA_SRC_DIR}/mfe/Lfold.c
    ${VRNA_SRC_DIR}/mfe/fold.c
)

# Backtrack
set(VRNA_BACKTRACK_SOURCES
    ${VRNA_SRC_DIR}/backtrack/bt_exterior.c
    ${VRNA_SRC_DIR}/backtrack/bt_exterior_f5.c
    ${VRNA_SRC_DIR}/backtrack/bt_exterior_f3.c
    ${VRNA_SRC_DIR}/backtrack/bt_hairpin.c
    ${VRNA_SRC_DIR}/backtrack/bt_internal.c
    ${VRNA_SRC_DIR}/backtrack/bt_multibranch.c
    ${VRNA_SRC_DIR}/backtrack/bt_gquad.c
)

# Partition function
set(VRNA_PARTFUNC_SOURCES
    ${VRNA_SRC_DIR}/partfunc/partfunc.c
    ${VRNA_SRC_DIR}/partfunc/pf_window.c
    ${VRNA_SRC_DIR}/partfunc/pf_wrappers.c
    ${VRNA_SRC_DIR}/partfunc/pf_multifold.c
    ${VRNA_SRC_DIR}/partfunc/pf_exterior.c
    ${VRNA_SRC_DIR}/partfunc/pf_internal.c
    ${VRNA_SRC_DIR}/partfunc/pf_multibranch.c
    ${VRNA_SRC_DIR}/partfunc/pf_gquad.c
    ${VRNA_SRC_DIR}/partfunc/alipfold.c
    ${VRNA_SRC_DIR}/partfunc/part_func_up.c
    ${VRNA_SRC_DIR}/partfunc/part_func_co.c
    ${VRNA_SRC_DIR}/partfunc/pf_fold.c
)

# Probabilities
set(VRNA_PROBABILITIES_SOURCES
    ${VRNA_SRC_DIR}/probabilities/equilibrium_probs.c
    ${VRNA_SRC_DIR}/probabilities/probs_structures.c
)

# Subopt
set(VRNA_SUBOPT_SOURCES
    ${VRNA_SRC_DIR}/subopt/subopt.c
    ${VRNA_SRC_DIR}/subopt/subopt_zuker.c
    ${VRNA_SRC_DIR}/subopt/subopt_gquad.c
)

# Sampling
set(VRNA_SAMPLING_SOURCES
    ${VRNA_SRC_DIR}/sampling/boltzmann_sampling.c
    ${VRNA_SRC_DIR}/sampling/bs_wrappers.c
)

# Params
set(VRNA_PARAMS_SOURCES
    ${VRNA_SRC_DIR}/params/io.c
    ${VRNA_SRC_DIR}/params/default.c
    ${VRNA_SRC_DIR}/params/params.c
    ${VRNA_SRC_DIR}/params/convert.c
    ${VRNA_SRC_DIR}/params/salt.c
    ${VRNA_SRC_DIR}/params/ribosum.c
)

# Inverse
set(VRNA_INVERSE_SOURCES
    ${VRNA_SRC_DIR}/inverse/inverse.c
)

# Combinatorics
set(VRNA_COMBINATORICS_SOURCES
    ${VRNA_SRC_DIR}/combinatorics/combinatorics.c
)

# Grammar
set(VRNA_GRAMMAR_SOURCES
    ${VRNA_SRC_DIR}/grammar/grammar.c
    ${VRNA_SRC_DIR}/grammar/gr_extension_mfe.c
    ${VRNA_SRC_DIR}/grammar/gr_extension_pf.c
)

# Datastructures
set(VRNA_DATASTRUCTURES_SOURCES
    ${VRNA_SRC_DIR}/datastructures/array.c
    ${VRNA_SRC_DIR}/datastructures/basic_datastructures.c
    ${VRNA_SRC_DIR}/datastructures/lists.c
    ${VRNA_SRC_DIR}/datastructures/char_stream.c
    ${VRNA_SRC_DIR}/datastructures/stream_output.c
    ${VRNA_SRC_DIR}/datastructures/string.c
    ${VRNA_SRC_DIR}/datastructures/hash_tables.c
    ${VRNA_SRC_DIR}/datastructures/heap.c
    ${VRNA_SRC_DIR}/datastructures/sparse_mx.c
    ${VRNA_SRC_DIR}/datastructures/dp_matrices.c
)

# Landscape
set(VRNA_LANDSCAPE_SOURCES
    ${VRNA_SRC_DIR}/move_set.c
    ${VRNA_SRC_DIR}/landscape/move.c
    ${VRNA_SRC_DIR}/landscape/findpath.c
    ${VRNA_SRC_DIR}/landscape/neighbor.c
    ${VRNA_SRC_DIR}/landscape/local_neighbors.c
    ${VRNA_SRC_DIR}/landscape/walk.c
)

# Special const + sanitize
set(VRNA_SPECIAL_SOURCES
    ${VRNA_SRC_DIR}/params/special_const.c
    ${VRNA_SRC_DIR}/io/sanitize.c
)

add_library(viennarna OBJECT
    ${VRNA_CONV_SOURCES}
    ${VRNA_JSON_SOURCES}
    ${VRNA_CEPHES_SOURCES}
    ${VRNA_EVAL_SOURCES}
    ${VRNA_CONSTRAINTS_SOURCES}
    ${VRNA_PROBING_SOURCES}
    ${VRNA_DATA_SOURCES}
    ${VRNA_MATH_SOURCES}
    ${VRNA_STRUCTURES_SOURCES}
    ${VRNA_SEQUENCES_SOURCES}
    ${VRNA_UTILS_SOURCES}
    ${VRNA_PLOTTING_SOURCES}
    ${VRNA_NAVIEW_SOURCES}
    ${VRNA_MFE_SOURCES}
    ${VRNA_BACKTRACK_SOURCES}
    ${VRNA_PARTFUNC_SOURCES}
    ${VRNA_PROBABILITIES_SOURCES}
    ${VRNA_SUBOPT_SOURCES}
    ${VRNA_SAMPLING_SOURCES}
    ${VRNA_PARAMS_SOURCES}
    ${VRNA_INVERSE_SOURCES}
    ${VRNA_COMBINATORICS_SOURCES}
    ${VRNA_GRAMMAR_SOURCES}
    ${VRNA_DATASTRUCTURES_SOURCES}
    ${VRNA_LANDSCAPE_SOURCES}
    ${VRNA_SPECIAL_SOURCES}
)

# ViennaRNA sources include headers as <ViennaRNA/...> so we need the parent dir
# Also add VRNA_SRC_DIR itself for relative .inc includes (e.g. "landscape/local_neighbors.inc")
target_include_directories(viennarna PUBLIC
    ${VRNA_SRC_DIR}/..
    ${VRNA_SRC_DIR}/../json
    ${VRNA_SRC_DIR}/../cephes
    ${CMAKE_CURRENT_SOURCE_DIR}  # for vrna_config.h as "config.h"
    ${VRNA_SRC_DIR}              # for relative .inc includes
)

# Provide config.h — we don't use autotools, so define HAVE_CONFIG_H
# and provide a minimal config.h in this directory
target_compile_definitions(viennarna PRIVATE
    HAVE_CONFIG_H
    VRNA_WITH_NAVIEW_LAYOUT
)
