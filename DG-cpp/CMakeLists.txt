cmake_minimum_required(VERSION 3.10)
project(PoissonCpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Default Release optimisation is -Ofast (= -O3 + -ffast-math + aggressive FP).
# -ffast-math relaxes IEEE semantics (reassociation, finite-math-only); it does
# NOT change the discretisation order, but for a clean convergence / round-off
# study reconfigure with -DDG_FAST_MATH=OFF to fall back to plain -O3.
option(DG_FAST_MATH "Use -Ofast (-ffast-math) for Release builds" ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    if(DG_FAST_MATH)
        set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
    else()
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    endif()
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

# ---- Eigen3 ----
# Prefer the bundled header-only copy (the system find_package can be hijacked by
# a stale user CMake package registry pointing at another project's build tree).
# Either way we expose the imported target Eigen3::Eigen so every dependent
# target links Eigen uniformly and portably.
if(EXISTS ${POISSON_EXTERNAL}/eigen3/Eigen/Dense)
    set(EIGEN3_INCLUDE_DIR ${POISSON_EXTERNAL}/eigen3)
    if(NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        set_target_properties(Eigen3::Eigen PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    endif()
    message(STATUS "Using bundled Eigen at ${EIGEN3_INCLUDE_DIR}")
else()
    find_package(Eigen3 3.3 NO_MODULE QUIET)
    if(NOT TARGET Eigen3::Eigen)
        find_package(Eigen3 REQUIRED NO_MODULE)
    endif()
    message(STATUS "Using Eigen ${Eigen3_VERSION}")
endif()

find_package(OpenMP)

# Robust geometry backend for DistMesh retriangulation.  The internal
# Bowyer-Watson fallback is useful for portability, but production body-fitted
# ALE remeshing should use CGAL's exact-predicate Delaunay triangulation.
option(DG_USE_CGAL "Use CGAL for 2-D Delaunay triangulation in MeshGen" ON)
if(DG_USE_CGAL)
    if(APPLE)
        foreach(_cgal_prefix /opt/homebrew/opt/cgal /usr/local/opt/cgal)
            if(EXISTS "${_cgal_prefix}/lib/cmake/CGAL/CGALConfig.cmake")
                list(APPEND CMAKE_PREFIX_PATH "${_cgal_prefix}")
            endif()
        endforeach()
    endif()
    find_package(CGAL CONFIG QUIET)
    if(NOT CGAL_FOUND)
        message(FATAL_ERROR
            "DG_USE_CGAL=ON requires CGAL. Install it before configuring, e.g. "
            "`brew install cgal` on macOS or `apt install libcgal-dev` on Debian/Ubuntu. "
            "Pass -DDG_USE_CGAL=OFF only for the less robust internal fallback.")
    endif()
    message(STATUS "Using CGAL ${CGAL_VERSION} for MeshGen Delaunay triangulation")
endif()

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
set(CH_DIR      ${CMAKE_CURRENT_LIST_DIR}/src/cahn_hilliard)
set(ARGYRIS_DIR ${CMAKE_CURRENT_LIST_DIR}/src/argyris)
set(NS_DIR      ${CMAKE_CURRENT_LIST_DIR}/src/navier_stokes)

add_library(poisson_common STATIC
    ${COMMON_DIR}/Mesh.cpp
    ${COMMON_DIR}/MeshGen.cpp
    ${COMMON_DIR}/FEM.cpp
    ${COMMON_DIR}/Quadrature.cpp)

target_include_directories(poisson_common PUBLIC
    ${COMMON_DIR})
# PUBLIC so every dependent target inherits Eigen's include path portably.
target_link_libraries(poisson_common PUBLIC Eigen3::Eigen)
if(DG_USE_CGAL)
    target_compile_definitions(poisson_common PUBLIC DG_USE_CGAL)
    if(TARGET CGAL::CGAL)
        target_link_libraries(poisson_common PUBLIC CGAL::CGAL)
    else()
        target_include_directories(poisson_common PUBLIC ${CGAL_INCLUDE_DIRS})
        target_link_libraries(poisson_common PUBLIC ${CGAL_LIBRARIES})
    endif()
endif()

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
    ${MIXED_DIR}/PostProcess.cpp
    ${MIXED_DIR}/HybridHDG.cpp)
target_include_directories(hdg PUBLIC
    ${MIXED_DIR}
    ${COMMON_DIR}
    ${IPDG_DIR})
target_link_libraries(hdg PUBLIC poisson_common dg_assembly)

add_library(argyris STATIC
    ${ARGYRIS_DIR}/Argyris.cpp)
target_include_directories(argyris PUBLIC
    ${ARGYRIS_DIR}
    ${COMMON_DIR})
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

# Navier-Stokes: DG / high-order-splitting incompressible solver.  The cylinder
# demo renders a von Karman vortex-street movie; the convergence test validates
# the spatial (k+1) and temporal (2nd) order on the Taylor-Green exact solution.
add_library(navier_stokes STATIC
    ${NS_DIR}/NavierStokes.cpp
    ${NS_DIR}/CosseratFilament.cpp
    ${NS_DIR}/IBCoupler.cpp
    ${NS_DIR}/ElasticTadpole.cpp
    ${NS_DIR}/Tadpole.cpp
    ${NS_DIR}/Spoon.cpp)
target_include_directories(navier_stokes PUBLIC
    ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR})
target_link_libraries(navier_stokes PUBLIC dg_assembly poisson_common)

add_executable(navier_stokes_cylinder
    ${NS_DIR}/ns_main.cpp)
target_include_directories(navier_stokes_cylinder PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR} ${CH_DIR})
target_link_libraries(navier_stokes_cylinder PRIVATE navier_stokes dg_assembly poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(navier_stokes_cylinder PRIVATE stdc++fs)  # <filesystem>
endif()

add_executable(navier_stokes_convergence
    ${NS_DIR}/ns_convergence_main.cpp)
target_include_directories(navier_stokes_convergence PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR})
target_link_libraries(navier_stokes_convergence PRIVATE navier_stokes dg_assembly poisson_common)

add_executable(ns_temporal_diag
    ${NS_DIR}/ns_temporal_diag.cpp)
target_include_directories(ns_temporal_diag PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR})
target_link_libraries(ns_temporal_diag PRIVATE navier_stokes dg_assembly poisson_common)

# Standalone test for the Cosserat-rod stepper: clamped cantilever, plucked
# at the tip; verifies the discrete period vs the analytic Euler-Bernoulli
# first-mode period and bounded energy drift.
add_executable(cosserat_test
    ${NS_DIR}/cosserat_test_main.cpp)
target_include_directories(cosserat_test PRIVATE ${NS_DIR} ${COMMON_DIR})
target_link_libraries(cosserat_test PRIVATE navier_stokes poisson_common)

# Cylinder + flexible filament FSI driver (Cosserat rod + immersed boundary).
add_executable(navier_stokes_filament
    ${NS_DIR}/ns_filament_main.cpp)
target_include_directories(navier_stokes_filament PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR} ${CH_DIR})
target_link_libraries(navier_stokes_filament PRIVATE navier_stokes dg_assembly poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(navier_stokes_filament PRIVATE stdc++fs)
endif()

# Cylinder + tadpole drift driver (rigid head + tail, station-keeping in the wake).
add_executable(navier_stokes_tadpole
    ${NS_DIR}/ns_tadpole_main.cpp)
target_include_directories(navier_stokes_tadpole PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR} ${CH_DIR})
target_link_libraries(navier_stokes_tadpole PRIVATE navier_stokes dg_assembly poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(navier_stokes_tadpole PRIVATE stdc++fs)
endif()

# Cylinder + tadpole drift driver with an immersed elastic tail.
add_executable(navier_stokes_tadpole_elastic
    ${NS_DIR}/ns_tadpole_elastic_main.cpp)
target_include_directories(navier_stokes_tadpole_elastic PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR} ${CH_DIR})
target_link_libraries(navier_stokes_tadpole_elastic PRIVATE navier_stokes dg_assembly poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(navier_stokes_tadpole_elastic PRIVATE stdc++fs)
endif()

# Closed circular bowl + a prescribed-motion spoon blade that stirs once,
# shedding a vortex dipole (immersed-boundary direct forcing).
add_executable(navier_stokes_spoon
    ${NS_DIR}/ns_spoon_main.cpp)
target_include_directories(navier_stokes_spoon PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR} ${CH_DIR})
target_link_libraries(navier_stokes_spoon PRIVATE navier_stokes dg_assembly poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(navier_stokes_spoon PRIVATE stdc++fs)
endif()

# Euler: DG (dP_k) + 2nd-order IMEX-RK solver for the 2-D compressible Euler
# equations.
set(EULER_DIR ${CMAKE_CURRENT_LIST_DIR}/src/euler)
set(EULER_CORE_DIR ${EULER_DIR}/core)
set(EULER_CASE_DIR ${EULER_DIR}/cases)
set(ALE_CORE_DIR ${EULER_DIR}/ALE_core)
set(ALE_CASE_DIR ${EULER_DIR}/ALE_cases)

add_library(euler STATIC
    ${EULER_CORE_DIR}/DG.cpp
    ${EULER_CORE_DIR}/DGUtils.cpp
    ${EULER_CORE_DIR}/DGRender.cpp
    ${EULER_CORE_DIR}/DGLTS.cpp
    ${EULER_CORE_DIR}/Render.cpp
    ${EULER_CORE_DIR}/AMR.cpp)
target_include_directories(euler PUBLIC
    ${EULER_CORE_DIR} ${EULER_DIR} ${COMMON_DIR})
target_link_libraries(euler PUBLIC poisson_common)
if(OpenMP_CXX_FOUND)
    target_link_libraries(euler PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(euler_convergence
    ${EULER_CASE_DIR}/convergence.cpp)
target_include_directories(euler_convergence PRIVATE ${EULER_CORE_DIR} ${COMMON_DIR})
target_link_libraries(euler_convergence PRIVATE euler poisson_common)

add_executable(euler_dmr
    ${EULER_CASE_DIR}/dmr.cpp)
target_include_directories(euler_dmr PRIVATE ${EULER_CORE_DIR} ${COMMON_DIR} ${CH_DIR})
target_link_libraries(euler_dmr PRIVATE euler poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(euler_dmr PRIVATE stdc++fs)
endif()

add_executable(euler_dmr_amr
    ${EULER_CASE_DIR}/dmr_amr.cpp)
target_include_directories(euler_dmr_amr PRIVATE ${EULER_CORE_DIR} ${COMMON_DIR} ${CH_DIR})
target_link_libraries(euler_dmr_amr PRIVATE euler poisson_common)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(euler_dmr_amr PRIVATE stdc++fs)
endif()

set(EULER_SCENE_TARGETS
    euler_riemann:riemann.cpp
    euler_shock_bubble:shock_bubble.cpp
    euler_corner_diffraction:corner.cpp
    euler_rmi_amr:rmi_amr.cpp)
foreach(scene_spec ${EULER_SCENE_TARGETS})
    string(REPLACE ":" ";" scene_parts ${scene_spec})
    list(GET scene_parts 0 scene_target)
    list(GET scene_parts 1 scene_source)
    add_executable(${scene_target} ${EULER_CASE_DIR}/${scene_source})
    target_include_directories(${scene_target} PRIVATE ${EULER_CORE_DIR} ${COMMON_DIR} ${CH_DIR})
    target_link_libraries(${scene_target} PRIVATE euler poisson_common)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
        target_link_libraries(${scene_target} PRIVATE stdc++fs)
    endif()
endforeach()

# ALE moving-mesh infrastructure
add_library(ale_core STATIC
    ${ALE_CORE_DIR}/Core.cpp
    ${ALE_CORE_DIR}/Checkpoint.cpp
    ${ALE_CORE_DIR}/DGState.cpp
    ${ALE_CORE_DIR}/Movie.cpp)
target_include_directories(ale_core PUBLIC
    ${ALE_CORE_DIR}
    ${ALE_CASE_DIR}
    ${EULER_CORE_DIR}
    ${EULER_DIR}
    ${COMMON_DIR})
target_link_libraries(ale_core PUBLIC euler poisson_common)

add_library(elastic_solid STATIC
    ${ALE_CORE_DIR}/ElasticSolid.cpp
    ${ALE_CORE_DIR}/SolidALEMap.cpp
    ${ALE_CORE_DIR}/NeoHookeanSolid.cpp)
target_link_libraries(elastic_solid PUBLIC ale_core)

add_library(body_mesh STATIC
    ${ALE_CORE_DIR}/BodyFittedMesh.cpp)
target_link_libraries(body_mesh PUBLIC elastic_solid)

add_library(piston_model STATIC
    ${ALE_CORE_DIR}/PistonModel.cpp)
target_link_libraries(piston_model PUBLIC ale_core)

add_executable(ale_convergence
    ${ALE_CASE_DIR}/Convergence.cpp)
target_link_libraries(ale_convergence PRIVATE elastic_solid)

add_executable(piston
    ${ALE_CASE_DIR}/Piston.cpp)
target_link_libraries(piston PRIVATE piston_model elastic_solid)

add_executable(fsi
    ${ALE_CASE_DIR}/FSI.cpp)
target_link_libraries(fsi PRIVATE piston_model elastic_solid)

add_executable(panel
    ${ALE_CASE_DIR}/Panel.cpp)
target_link_libraries(panel PRIVATE elastic_solid)

add_executable(turek_flag
    ${ALE_CASE_DIR}/TurekFlag.cpp)
target_link_libraries(turek_flag PRIVATE elastic_solid)

add_executable(shock_disk
    ${ALE_CASE_DIR}/ShockDisk.cpp)
target_link_libraries(shock_disk PRIVATE body_mesh)

add_executable(blast_rod
    ${ALE_CASE_DIR}/BlastRod/Case.cpp
    ${ALE_CASE_DIR}/BlastRod/Coupling.cpp
    ${ALE_CASE_DIR}/BlastRod/Domain.cpp
    ${ALE_CASE_DIR}/BlastRod/Geometry.cpp)
target_include_directories(blast_rod PRIVATE ${ALE_CASE_DIR}/BlastRod ${ALE_CASE_DIR} ${CH_DIR})
target_link_libraries(blast_rod PRIVATE body_mesh)

add_executable(blast_rod_cg
    ${ALE_CASE_DIR}/BlastRod/CGCase.cpp
    ${ALE_CASE_DIR}/BlastRod/Domain.cpp
    ${ALE_CASE_DIR}/BlastRod/Geometry.cpp)
target_include_directories(blast_rod_cg PRIVATE ${ALE_CASE_DIR}/BlastRod ${ALE_CASE_DIR})
target_link_libraries(blast_rod_cg PRIVATE body_mesh)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    foreach(tgt piston fsi panel turek_flag shock_disk blast_rod blast_rod_cg)
        target_link_libraries(${tgt} PRIVATE stdc++fs)
    endforeach()
endif()

if(SUITESPARSE_FOUND)
    target_include_directories(poisson_common PUBLIC ${CHOLMOD_INCLUDE_DIR})
    set(POISSON_EXES poisson_ipdg poisson_mixed_hdg poisson_hdg_hybrid poisson_pk
                     biharmonic_ipcg biharmonic_argyris cahn_hilliard cahn_hilliard_convergence)
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

add_executable(ns_cyl_order
    ${NS_DIR}/ns_cyl_order.cpp)
target_include_directories(ns_cyl_order PRIVATE ${NS_DIR} ${COMMON_DIR} ${IPDG_DIR})
target_link_libraries(ns_cyl_order PRIVATE navier_stokes dg_assembly poisson_common)
