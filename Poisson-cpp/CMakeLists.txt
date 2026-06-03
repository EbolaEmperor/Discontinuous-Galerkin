cmake_minimum_required(VERSION 3.10)
project(PoissonCpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# Eigen: prefer a 3.x package (the codebase targets the Eigen 3 API), but accept
# any Eigen whose CMake config is found. Targets link the imported Eigen3::Eigen
# target (which always carries the correct include path), rather than the legacy
# EIGEN3_INCLUDE_DIR variable that newer Eigen configs no longer set.
find_package(Eigen3 3.3 NO_MODULE QUIET)
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 REQUIRED NO_MODULE)
endif()
message(STATUS "Using Eigen ${Eigen3_VERSION}")
find_package(OpenMP)

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
set(CG_DIR     ${CMAKE_CURRENT_LIST_DIR}/src/cg)
set(CH_DIR     ${CMAKE_CURRENT_LIST_DIR}/src/cahn_hilliard)

add_library(poisson_common STATIC
    ${COMMON_DIR}/Mesh.cpp
    ${COMMON_DIR}/FEM.cpp
    ${COMMON_DIR}/Quadrature.cpp)

target_include_directories(poisson_common PUBLIC
    ${COMMON_DIR})
# PUBLIC so every dependent target inherits Eigen's include path portably.
target_link_libraries(poisson_common PUBLIC Eigen3::Eigen)

add_library(dg_assembly STATIC
    ${IPDG_DIR}/DGAssembly.cpp)
target_include_directories(dg_assembly PUBLIC
    ${IPDG_DIR}
    ${COMMON_DIR})
target_link_libraries(dg_assembly PUBLIC poisson_common)
if(OpenMP_CXX_FOUND)
    target_link_libraries(dg_assembly PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP not found; assembly will run single-threaded.")
endif()

add_library(hdg STATIC
    ${MIXED_DIR}/VecFEM.cpp
    ${MIXED_DIR}/AssemblyHDG.cpp
    ${MIXED_DIR}/PostProcess.cpp)
target_include_directories(hdg PUBLIC
    ${MIXED_DIR}
    ${COMMON_DIR}
    ${IPDG_DIR})
target_link_libraries(hdg PUBLIC poisson_common dg_assembly)

add_executable(poisson_ipdg
    ${IPDG_DIR}/ipdg_main.cpp)
target_link_libraries(poisson_ipdg PRIVATE dg_assembly poisson_common)

add_executable(poisson_mixed_hdg
    ${MIXED_DIR}/mixed_main.cpp)
target_link_libraries(poisson_mixed_hdg PRIVATE hdg dg_assembly poisson_common)

add_executable(poisson_pk
    ${CG_DIR}/pk_main.cpp)
target_link_libraries(poisson_pk PRIVATE dg_assembly poisson_common)

add_executable(biharmonic_ipcg
    ${IPDG_DIR}/biharmonic_ipcg_main.cpp)
target_link_libraries(biharmonic_ipcg PRIVATE dg_assembly poisson_common)

# Cahn-Hilliard: mixed-DG IMEX time-stepping demo that renders an evolution movie.
add_executable(cahn_hilliard
    ${CH_DIR}/ch_main.cpp
    ${CH_DIR}/CahnHilliard.cpp)
target_include_directories(cahn_hilliard PRIVATE
    ${CH_DIR} ${COMMON_DIR} ${IPDG_DIR})
target_link_libraries(cahn_hilliard PRIVATE dg_assembly poisson_common)
# ch_main.cpp uses <filesystem>; GCC < 9.1 / libc++ < 9 need an extra runtime lib.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(cahn_hilliard PRIVATE stdc++fs)
endif()

# Cahn-Hilliard space-time convergence test (MMS, no video).
add_executable(cahn_hilliard_convergence
    ${CH_DIR}/ch_convergence_main.cpp
    ${CH_DIR}/CahnHilliard.cpp)
target_include_directories(cahn_hilliard_convergence PRIVATE
    ${CH_DIR} ${COMMON_DIR} ${IPDG_DIR})
target_link_libraries(cahn_hilliard_convergence PRIVATE dg_assembly poisson_common)

if(SUITESPARSE_FOUND)
    target_include_directories(poisson_common PUBLIC ${CHOLMOD_INCLUDE_DIR})
    foreach(tgt poisson_ipdg poisson_mixed_hdg poisson_pk biharmonic_ipcg cahn_hilliard cahn_hilliard_convergence)
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
