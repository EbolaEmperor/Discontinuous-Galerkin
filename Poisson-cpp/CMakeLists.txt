cmake_minimum_required(VERSION 3.10)
project(PoissonCpp)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
endif()

# Eigen (and DGAssembly) emit OpenMP calls; make sure the OpenMP runtime is on
# every executable's link line.  The OpenMP imported target does not always
# propagate its link flag under the bundled MinGW toolchain, so add it directly.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
endif()

# Local prefix holding bundled dependencies (Eigen headers + prebuilt
# SuiteSparse/CHOLMOD for the MinGW toolchain). See external/ (git-ignored).
set(POISSON_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external)
set(SS_PREFIX ${POISSON_EXTERNAL}/ssprefix/mingw64)

# ---- Eigen3 (header-only) ----
# Prefer the bundled copy (the system find_package can be hijacked by a stale
# user CMake package registry pointing at another project's build tree).
if(EXISTS ${POISSON_EXTERNAL}/eigen3/Eigen/Dense)
    set(EIGEN3_INCLUDE_DIR ${POISSON_EXTERNAL}/eigen3)
    message(STATUS "Using bundled Eigen at ${EIGEN3_INCLUDE_DIR}")
else()
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif()

find_package(OpenMP)

# Optional SuiteSparse / CHOLMOD. The bundled MinGW prefix (SS_PREFIX) is
# searched first, then the usual macOS/Linux locations.
set(SUITESPARSE_ROOT "/opt/homebrew/opt/suite-sparse")
find_path(CHOLMOD_INCLUDE_DIR cholmod.h
    HINTS ${SS_PREFIX}/include/suitesparse
          ${SS_PREFIX}/include
          ${SUITESPARSE_ROOT}/include/suitesparse
          ${SUITESPARSE_ROOT}/include
          /usr/local/include/suitesparse
          /usr/include/suitesparse)
find_library(CHOLMOD_LIBRARY cholmod
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(AMD_LIBRARY amd
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(COLAMD_LIBRARY colamd
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(CAMD_LIBRARY camd
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(CCOLAMD_LIBRARY ccolamd
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(OPENBLAS_LIBRARY NAMES openblas
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /usr/local/lib)
find_library(METIS_LIBRARY metis
    HINTS ${SS_PREFIX}/lib ${SUITESPARSE_ROOT}/lib /opt/homebrew/opt/metis/lib /usr/local/lib)

if(CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARY)
    set(SUITESPARSE_FOUND TRUE)
    add_definitions(-DEIGEN_USE_CHOLMOD)
endif()

set(COMMON_DIR  ${CMAKE_CURRENT_LIST_DIR}/src/common)
set(IPDG_DIR    ${CMAKE_CURRENT_LIST_DIR}/src/ipdg)
set(MIXED_DIR   ${CMAKE_CURRENT_LIST_DIR}/src/mixed)
set(CG_DIR      ${CMAKE_CURRENT_LIST_DIR}/src/cg)
set(ARGYRIS_DIR ${CMAKE_CURRENT_LIST_DIR}/src/argyris)

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
if(OpenMP_CXX_FOUND)
    target_link_libraries(dg_assembly PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP not found; assembly will run single-threaded.")
endif()

add_library(hdg STATIC
    ${MIXED_DIR}/VecFEM.cpp
    ${MIXED_DIR}/AssemblyHDG.cpp
    ${MIXED_DIR}/PostProcess.cpp
    ${MIXED_DIR}/HybridHDG.cpp)
target_include_directories(hdg PUBLIC
    ${MIXED_DIR}
    ${COMMON_DIR}
    ${IPDG_DIR}
    ${EIGEN3_INCLUDE_DIR})
target_link_libraries(hdg PUBLIC poisson_common dg_assembly)

add_library(argyris STATIC
    ${ARGYRIS_DIR}/Argyris.cpp)
target_include_directories(argyris PUBLIC
    ${ARGYRIS_DIR}
    ${COMMON_DIR}
    ${EIGEN3_INCLUDE_DIR})
target_link_libraries(argyris PUBLIC poisson_common)

add_executable(poisson_ipdg
    ${IPDG_DIR}/ipdg_main.cpp)
target_link_libraries(poisson_ipdg PRIVATE dg_assembly poisson_common)

add_executable(poisson_mixed_hdg
    ${MIXED_DIR}/mixed_main.cpp)
target_link_libraries(poisson_mixed_hdg PRIVATE hdg dg_assembly poisson_common)

add_executable(poisson_hdg_hybrid
    ${MIXED_DIR}/hybrid_main.cpp)
target_link_libraries(poisson_hdg_hybrid PRIVATE hdg dg_assembly poisson_common)

add_executable(poisson_pk
    ${CG_DIR}/pk_main.cpp)
target_link_libraries(poisson_pk PRIVATE dg_assembly poisson_common)

add_executable(biharmonic_ipcg
    ${IPDG_DIR}/biharmonic_ipcg_main.cpp)
target_link_libraries(biharmonic_ipcg PRIVATE dg_assembly poisson_common)

add_executable(biharmonic_argyris
    ${ARGYRIS_DIR}/argyris_main.cpp)
target_link_libraries(biharmonic_argyris PRIVATE argyris poisson_common)

if(SUITESPARSE_FOUND)
    target_include_directories(poisson_common PUBLIC ${CHOLMOD_INCLUDE_DIR})
    set(POISSON_EXES poisson_ipdg poisson_mixed_hdg poisson_pk biharmonic_ipcg poisson_hdg_hybrid biharmonic_argyris)
    foreach(tgt ${POISSON_EXES})
        target_link_libraries(${tgt} PRIVATE
            ${CHOLMOD_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY}
            ${CAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${SUITESPARSECONFIG_LIBRARY})
        if(OPENBLAS_LIBRARY)
            target_link_libraries(${tgt} PRIVATE ${OPENBLAS_LIBRARY})
        endif()
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

    # On Windows the prebuilt CHOLMOD is a set of DLLs; copy the runtime
    # closure next to each executable so they run from the build tree.
    if(WIN32 AND EXISTS ${SS_PREFIX}/bin)
        set(SS_RUNTIME_DLLS
            libcholmod libamd libcamd libccolamd libcolamd libsuitesparseconfig
            libopenblas libgfortran-5 libquadmath-0 libgomp-1 libgcc_s_seh-1
            libwinpthread-1 libstdc++-6)
        foreach(tgt ${POISSON_EXES})
            foreach(dll ${SS_RUNTIME_DLLS})
                if(EXISTS ${SS_PREFIX}/bin/${dll}.dll)
                    add_custom_command(TARGET ${tgt} POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                                ${SS_PREFIX}/bin/${dll}.dll $<TARGET_FILE_DIR:${tgt}>)
                endif()
            endforeach()
        endforeach()
    endif()
endif()
