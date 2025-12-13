cmake_minimum_required(VERSION 3.10)
project(PoissonCpp)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Optional SuiteSparse / CHOLMOD
set(SUITESPARSE_ROOT "/opt/homebrew/opt/suite-sparse")
find_path(CHOLMOD_INCLUDE_DIR cholmod.h
    HINTS ${SUITESPARSE_ROOT}/include/suitesparse
          ${SUITESPARSE_ROOT}/include
          /usr/local/include/suitesparse
          /usr/include/suitesparse)
find_library(CHOLMOD_LIBRARY cholmod
    HINTS ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(AMD_LIBRARY amd
    HINTS ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(COLAMD_LIBRARY colamd
    HINTS ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(CAMD_LIBRARY camd
    HINTS ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(CCOLAMD_LIBRARY ccolamd
    HINTS ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig
    HINTS ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(METIS_LIBRARY metis
    HINTS ${SUITESPARSE_ROOT}/lib /opt/homebrew/opt/metis/lib /usr/local/lib)

if(CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARY)
    set(SUITESPARSE_FOUND TRUE)
    add_definitions(-DEIGEN_USE_CHOLMOD)
endif()

set(COMMON_DIR ${CMAKE_CURRENT_LIST_DIR}/src/common)
set(IPDG_DIR   ${CMAKE_CURRENT_LIST_DIR}/src/ipdg)
set(MIXED_DIR  ${CMAKE_CURRENT_LIST_DIR}/src/mixed)

add_library(poisson_common STATIC
    ${COMMON_DIR}/Mesh.cpp
    ${COMMON_DIR}/FEM.cpp
    ${COMMON_DIR}/Quadrature.cpp)

target_include_directories(poisson_common PUBLIC
    ${COMMON_DIR}
    ${EIGEN3_INCLUDE_DIR})

add_library(dg_assembly STATIC
    ${IPDG_DIR}/DGAssembly.cpp)
target_include_directories(dg_assembly PUBLIC
    ${IPDG_DIR}
    ${COMMON_DIR}
    ${EIGEN3_INCLUDE_DIR})
target_link_libraries(dg_assembly PUBLIC poisson_common)

add_library(hdg STATIC
    ${MIXED_DIR}/VecFEM.cpp
    ${MIXED_DIR}/AssemblyHDG.cpp
    ${MIXED_DIR}/PostProcess.cpp)
target_include_directories(hdg PUBLIC
    ${MIXED_DIR}
    ${COMMON_DIR}
    ${IPDG_DIR}
    ${EIGEN3_INCLUDE_DIR})
target_link_libraries(hdg PUBLIC poisson_common dg_assembly)

add_executable(poisson_ipdg
    ${IPDG_DIR}/ipdg_main.cpp)
target_link_libraries(poisson_ipdg PRIVATE dg_assembly poisson_common)

add_executable(poisson_mixed_hdg
    ${MIXED_DIR}/mixed_main.cpp)
target_link_libraries(poisson_mixed_hdg PRIVATE hdg dg_assembly poisson_common)

if(SUITESPARSE_FOUND)
    target_include_directories(poisson_common PUBLIC ${CHOLMOD_INCLUDE_DIR})
    foreach(tgt poisson_ipdg poisson_mixed_hdg)
        target_link_libraries(${tgt} PRIVATE
            ${CHOLMOD_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY}
            ${CAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${SUITESPARSECONFIG_LIBRARY})
        if(METIS_LIBRARY)
            target_link_libraries(${tgt} PRIVATE ${METIS_LIBRARY})
        endif()
        if(APPLE)
            find_library(ACCELERATE_FRAMEWORK Accelerate)
            if(ACCELERATE_FRAMEWORK)
                target_link_libraries(${tgt} PRIVATE ${ACCELERATE_FRAMEWORK})
            endif()
        endif()
    endforeach()
endif()
