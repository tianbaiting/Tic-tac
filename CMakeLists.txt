# Tic-tac 三体核物理计算程序
# CMake 构建系统 - 重构版本

cmake_minimum_required(VERSION 3.12)

# 生成 compile_commands.json 用于 clangd/LSP
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Generate compile_commands.json")

project(Tic-tac
    VERSION 1.0.0
    DESCRIPTION "Three-body nuclear physics Faddeev equation solver"
    LANGUAGES CXX Fortran)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置构建类型
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-unused-variable")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Fortran编译选项
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fdefault-real-8 -fdefault-double-8 -fPIC")

# 查找必要的包
find_package(OpenMP COMPONENTS CXX Fortran)
find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

# Cache layer rollback switch
option(TICTAC_USE_NEW_CACHE_LAYER "Use the hash-keyed cache layer for P123/W1" ON)
if(TICTAC_USE_NEW_CACHE_LAYER)
    add_compile_definitions(TICTAC_USE_NEW_CACHE_LAYER=1)
else()
    add_compile_definitions(TICTAC_USE_NEW_CACHE_LAYER=0)
endif()

# BLAS和LAPACK - 使用更简单的查找方式
find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)
find_library(LAPACKE_LIB lapacke)

# 包含目录
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/app
    ${CMAKE_SOURCE_DIR}/src/config
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/core/faddeev_solver
    ${CMAKE_SOURCE_DIR}/src/core/state_space
    ${CMAKE_SOURCE_DIR}/src/core/potential
    ${CMAKE_SOURCE_DIR}/src/core/resolvent
    ${CMAKE_SOURCE_DIR}/src/interactions
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/utils
    ${HDF5_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
)

# 收集所有源文件
file(GLOB_RECURSE SOURCES_APP "${CMAKE_SOURCE_DIR}/src/app/*.cpp")
file(GLOB_RECURSE SOURCES_CONFIG "${CMAKE_SOURCE_DIR}/src/config/*.cpp")
file(GLOB_RECURSE SOURCES_CORE "${CMAKE_SOURCE_DIR}/src/core/*/*.cpp")
file(GLOB_RECURSE SOURCES_INTERACTIONS "${CMAKE_SOURCE_DIR}/src/interactions/*.cpp")
file(GLOB_RECURSE SOURCES_IO "${CMAKE_SOURCE_DIR}/src/io/*.cpp")
file(GLOB_RECURSE SOURCES_UTILS "${CMAKE_SOURCE_DIR}/src/utils/*.cpp")

# Fortran源文件
file(GLOB_RECURSE SOURCES_FORTRAN_F90 "${CMAKE_SOURCE_DIR}/src/interactions/**/*.f90")
file(GLOB_RECURSE SOURCES_FORTRAN_F77 "${CMAKE_SOURCE_DIR}/src/interactions/**/*.f")

foreach(source_file ${SOURCES_FORTRAN_F90})
    set_source_files_properties(
        ${source_file}
        PROPERTIES
        Fortran_FORMAT FREE
        COMPILE_OPTIONS "-cpp;-ffree-line-length-1000;-std=legacy"
    )
endforeach()

foreach(source_file ${SOURCES_FORTRAN_F77})
    set_source_files_properties(
        ${source_file}
        PROPERTIES
        Fortran_FORMAT FIXED
        COMPILE_OPTIONS "-std=legacy;-ffixed-line-length-none"
    )
endforeach()

# 主程序源文件
set(SOURCES_MAIN "${CMAKE_SOURCE_DIR}/src/main.cpp")

# 所有源文件
set(ALL_SOURCES
    ${SOURCES_MAIN}
    ${SOURCES_APP}
    ${SOURCES_CONFIG}
    ${SOURCES_CORE}
    ${SOURCES_INTERACTIONS}
    ${SOURCES_IO}
    ${SOURCES_UTILS}
    ${SOURCES_FORTRAN_F90}
    ${SOURCES_FORTRAN_F77}
)

# 创建可执行文件
add_executable(${PROJECT_NAME} 
    ${ALL_SOURCES}
)

# 链接库
target_link_libraries(${PROJECT_NAME}
    OpenMP::OpenMP_CXX
    ${GSL_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${LAPACKE_LIB}
    ${LAPACK_LIB}
    ${BLAS_LIB}
    pthread
    m
    dl
    gfortran
    stdc++fs
    gomp
)

# 设置输出目录
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Tests
enable_testing()
add_subdirectory(tests/cpp)

# W1 block-database worker (Phase D/E/F): plan/build/verify exact W1 blocks
# without running the Faddeev solve. Embarrassingly parallel via --shard K/N.
add_executable(w1_worker
    ${CMAKE_SOURCE_DIR}/tools/w1_worker/w1_worker.cpp
    ${CMAKE_SOURCE_DIR}/src/config/set_run_parameters.cpp
    ${CMAKE_SOURCE_DIR}/src/config/run_organizer.cpp
    ${CMAKE_SOURCE_DIR}/src/core/state_space/make_pw_symm_states.cpp
    ${CMAKE_SOURCE_DIR}/src/core/state_space/make_wp_states.cpp
    ${CMAKE_SOURCE_DIR}/src/io/disk_io_routines.cpp
    ${CMAKE_SOURCE_DIR}/src/interactions/w1_pw_cache.cpp
    ${CMAKE_SOURCE_DIR}/src/interactions/chiral_N2LO_3NF_factorized.cpp
    ${CMAKE_SOURCE_DIR}/src/interactions/chiral_N2LO_3NF_full_reference.cpp
    ${CMAKE_SOURCE_DIR}/src/interactions/three_nucleon_force_model.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/coupling_coefficients.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/spin_isospin_algebra.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/chiral_3nf_recoupling.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/gauss_legendre.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/error_management.cpp
    ${CMAKE_SOURCE_DIR}/src/utils/kinetic_conversion.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_keys.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_manifest.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_p123.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_io_w1.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/cache_layer.cpp
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party/sha256.cpp
)
target_include_directories(w1_worker PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/config
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/core/state_space
    ${CMAKE_SOURCE_DIR}/src/interactions
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer
    ${CMAKE_SOURCE_DIR}/src/io/cache_layer/third_party
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(w1_worker
    OpenMP::OpenMP_CXX
    ${GSL_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
)
set_target_properties(w1_worker PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# 安装规则
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)

# 打印配置信息
message(STATUS "=== Tic-tac 构建配置 ===")
message(STATUS "构建类型: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ 编译器: ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran 编译器: ${CMAKE_Fortran_COMPILER}")
message(STATUS "C++ 标准: ${CMAKE_CXX_STANDARD}")
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "GSL: ${GSL_FOUND}")
message(STATUS "HDF5: ${HDF5_FOUND}")
message(STATUS "========================")

# CMAKE_EXPORT_COMPILE_COMMANDS 已在文件开头设置
