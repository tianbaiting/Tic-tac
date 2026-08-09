# C++ unit tests for 3NF components
# Each test is a standalone executable that returns 0 on success, non-zero on failure.

# Inherit parent include directories and link libraries
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/interactions
    ${CMAKE_SOURCE_DIR}/src/config
    ${CMAKE_SOURCE_DIR}/src/core/state_space
)

# Shared source files needed by tests (subset of main build)
set(TEST_UTIL_SOURCES
    ${CMAKE_SOURCE_DIR}/src/utils/coupling_coefficients.cpp
)

# L1: Coupling coefficient tests
add_executable(test_coupling_coefficients
    test_coupling_coefficients.cpp
    ${TEST_UTIL_SOURCES}
)
target_link_libraries(test_coupling_coefficients ${GSL_LIBRARIES})
add_test(NAME coupling_coefficients COMMAND test_coupling_coefficients)

# L2: Spin-isospin algebra tests
add_executable(test_spin_isospin_algebra
    test_spin_isospin_algebra.cpp
    ${TEST_UTIL_SOURCES}
    ${CMAKE_SOURCE_DIR}/src/utils/spin_isospin_algebra.cpp
)
target_link_libraries(test_spin_isospin_algebra ${GSL_LIBRARIES})
add_test(NAME spin_isospin_algebra COMMAND test_spin_isospin_algebra)

# L2: 3NF matrix element tests
set(TEST_3NF_SOURCES
    ${TEST_UTIL_SOURCES}
    ${CMAKE_SOURCE_DIR}/src/utils/spin_isospin_algebra.cpp
    ${CMAKE_SOURCE_DIR}/src/core/state_space/make_pw_symm_states.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/error_management.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/chiral_3nf_recoupling.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/gauss_legendre.cpp
)
add_executable(test_3nf_matrix_elements
    test_3nf_matrix_elements.cpp
    ${TEST_3NF_SOURCES}
)
target_link_libraries(test_3nf_matrix_elements ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
add_test(NAME 3nf_matrix_elements COMMAND test_3nf_matrix_elements)

# L3: 3NF symmetry tests
add_executable(test_3nf_symmetries
    test_3nf_symmetries.cpp
    ${TEST_3NF_SOURCES}
)
target_link_libraries(test_3nf_symmetries ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
add_test(NAME 3nf_symmetries COMMAND test_3nf_symmetries)

# L2: Chiral 3NF spin-isospin recoupling helpers (scalar + rank-2)
add_executable(test_chiral_3nf_recoupling
    test_chiral_3nf_recoupling.cpp
    ${TEST_UTIL_SOURCES}
    ${CMAKE_SOURCE_DIR}/src/utils/chiral_3nf_recoupling.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/spin_isospin_algebra.cpp
)
target_link_libraries(test_chiral_3nf_recoupling ${GSL_LIBRARIES})
add_test(NAME chiral_3nf_recoupling COMMAND test_chiral_3nf_recoupling)

# L4: Faddeev/AGS operator-ordering dense test (3NF audit B4)
add_executable(test_faddeev_operator_order
    test_faddeev_operator_order.cpp
)
target_link_libraries(test_faddeev_operator_order)
add_test(NAME faddeev_operator_order COMMAND test_faddeev_operator_order)

# L0: Phase 0 operator-level oracle (fix/3nf-physics-contract).
# Exercises the REAL calculate_CPVC_col with mock P,V,W1,C and verifies the
# algebra equals C^T·(P·V + W^(1)·(1+P))·C. See docs/three_nf_equation_contract.md §6.
add_executable(test_3nf_operator_oracle
    test_3nf_operator_oracle.cpp
    ${CMAKE_SOURCE_DIR}/src/core/faddeev_solver/cpvc_kernel.cpp
    ${CMAKE_SOURCE_DIR}/src/core/faddeev_solver/coupled_channel_transform.cpp
    ${CMAKE_SOURCE_DIR}/src/core/faddeev_solver/pade_approximant.cpp
    ${CMAKE_SOURCE_DIR}/src/interactions/w1_pw_cache.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/matrix_routines.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/gauss_legendre.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_keys.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_manifest.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_p123.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_w1.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_layer.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party/sha256.cpp
)
target_include_directories(test_3nf_operator_oracle PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/core/faddeev_solver
    ${CMAKE_SOURCE_DIR}/src/interactions
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party
    ${CMAKE_SOURCE_DIR}/include
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(test_3nf_operator_oracle
    OpenMP::OpenMP_CXX
    ${GSL_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${LAPACKE_LIB}
    ${LAPACK_LIB}
    ${BLAS_LIB}
)
add_test(NAME 3nf_operator_oracle COMMAND test_3nf_operator_oracle)

# Phase 3 oracle driver: prints production W1_element for the independent
# Python angular oracle (tools/3nf_oracle/angular_oracle.py). NOT a test — a
# thin wrapper around chiral_N2LO_3NF::W1_element so the oracle can compare
# against the real production code without re-transcribing the formula.
add_executable(print_w1_element
    ${CMAKE_SOURCE_DIR}/tools/3nf_oracle/print_w1_element.cpp
    ${CMAKE_SOURCE_DIR}/src/interactions/three_nucleon_force_model.cpp
    ${CMAKE_SOURCE_DIR}/src/core/state_space/make_pw_symm_states.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/coupling_coefficients.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/spin_isospin_algebra.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/chiral_3nf_recoupling.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/gauss_legendre.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/error_management.cpp
)
target_include_directories(print_w1_element PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/interactions
    ${CMAKE_SOURCE_DIR}/src/core/state_space
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(print_w1_element ${GSL_LIBRARIES})
set_target_properties(print_w1_element PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools/3nf_oracle)

# Resolvent analytic Im parts
add_executable(resolvent_im_test
    resolvent_im_test.cpp
    ${CMAKE_SOURCE_DIR}/src/core/resolvent/make_resolvent.cpp
    ${CMAKE_SOURCE_DIR}/src/core/state_space/make_pw_symm_states.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/error_management.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/coupling_coefficients.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/spin_isospin_algebra.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/kinetic_conversion.cpp
)
target_include_directories(resolvent_im_test PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/core/resolvent
    ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(resolvent_im_test ${GSL_LIBRARIES})
add_test(NAME resolvent_im_test COMMAND resolvent_im_test)

# Cache layer C++ unit tests
add_executable(cache_layer_test
    cache_layer_test.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_keys.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_manifest.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_p123.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_w1.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_layer.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party/sha256.cpp
)
target_include_directories(cache_layer_test PRIVATE
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party
    ${CMAKE_SOURCE_DIR}/include
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(cache_layer_test PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
add_test(NAME cache_layer_test COMMAND cache_layer_test)
