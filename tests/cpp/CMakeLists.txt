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
)
add_executable(test_3nf_matrix_elements
    test_3nf_matrix_elements.cpp
    ${TEST_3NF_SOURCES}
)
target_link_libraries(test_3nf_matrix_elements ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
add_test(NAME 3nf_matrix_elements COMMAND test_3nf_matrix_elements)

# L3: 3NF symmetry tests (added in Task 7)
