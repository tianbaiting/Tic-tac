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

# Cache layer C++ unit tests
# Note: cache_keys.cpp etc. are added in later tasks. We include them only if they exist
# so CMake configure succeeds even before those files are introduced. The test executable
# may fail to link until all sources are present, which is acceptable for early bootstrap.
set(CACHE_LAYER_TEST_SOURCES cache_layer_test.cpp)
foreach(_src
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_keys.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_manifest.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_layer.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_p123.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_w1.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party/sha256.cpp)
    if(EXISTS ${_src})
        list(APPEND CACHE_LAYER_TEST_SOURCES ${_src})
    endif()
endforeach()
add_executable(cache_layer_test ${CACHE_LAYER_TEST_SOURCES})
target_include_directories(cache_layer_test PRIVATE
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party
    ${CMAKE_SOURCE_DIR}/include
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(cache_layer_test PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
add_test(NAME cache_layer_test COMMAND cache_layer_test)
