# Tic-tac 三体核物理计算程序
# CMake 构建系统 - 重构版本

cmake_minimum_required(VERSION 3.12)
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
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fdefault-real-8 -fdefault-double-8")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-form -ffree-line-length-1000 -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

# 查找必要的包
find_package(OpenMP COMPONENTS CXX)
find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

# BLAS和LAPACK - 使用更简单的查找方式
find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)
find_library(LAPACKE_LIB lapacke)

# 包含目录
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
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
file(GLOB_RECURSE SOURCES_CONFIG "${CMAKE_SOURCE_DIR}/src/config/*.cpp")
file(GLOB_RECURSE SOURCES_CORE "${CMAKE_SOURCE_DIR}/src/core/*/*.cpp")
file(GLOB_RECURSE SOURCES_INTERACTIONS "${CMAKE_SOURCE_DIR}/src/interactions/*.cpp")
file(GLOB_RECURSE SOURCES_IO "${CMAKE_SOURCE_DIR}/src/io/*.cpp")
file(GLOB_RECURSE SOURCES_UTILS "${CMAKE_SOURCE_DIR}/src/utils/*.cpp")

# Fortran源文件
file(GLOB_RECURSE SOURCES_FORTRAN 
    "${CMAKE_SOURCE_DIR}/src/interactions/chp/*.f90"
#    "${CMAKE_SOURCE_DIR}/src/interactions/nijmegen/*.f"  # 暂时跳过有问题的Fortran 77文件
)

# 使用CPP目录中已编译的nijmegen对象文件
set(NIJMEGEN_OBJECTS
    "${CMAKE_SOURCE_DIR}/CPP/Interactions/nijmegen/nijmegen_interface.o"
    "${CMAKE_SOURCE_DIR}/CPP/Interactions/nijmegen/pnijm.o"
)

# 主程序源文件
set(SOURCES_MAIN "${CMAKE_SOURCE_DIR}/src/main.cpp")

# 所有源文件
set(ALL_SOURCES
    ${SOURCES_MAIN}
    ${SOURCES_CONFIG}
    ${SOURCES_CORE}
    ${SOURCES_INTERACTIONS}
    ${SOURCES_IO}
    ${SOURCES_UTILS}
    ${SOURCES_FORTRAN}
)

# 创建可执行文件
add_executable(${PROJECT_NAME} 
    ${ALL_SOURCES}
    ${NIJMEGEN_OBJECTS}  # 添加预编译的nijmegen对象文件
)

# 链接库
target_link_libraries(${PROJECT_NAME}
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

# 可选：创建编译命令数据库（用于IDE和代码补全）
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
