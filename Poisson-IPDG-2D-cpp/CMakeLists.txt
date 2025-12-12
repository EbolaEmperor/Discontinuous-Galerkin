cmake_minimum_required(VERSION 3.10)
project(Poisson-IPDG-2D-cpp)

set(CMAKE_CXX_STANDARD 17)

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
endif()

# Find Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find SuiteSparse (CHOLMOD)
# On macOS with Homebrew, paths usually need hints
set(SUITESPARSE_ROOT "/opt/homebrew/opt/suite-sparse")
find_path(CHOLMOD_INCLUDE_DIR cholmod.h
    HINTS ${SUITESPARSE_ROOT}/include/suitesparse
    ${SUITESPARSE_ROOT}/include
    /usr/local/include/suitesparse
    /usr/local/include
    /usr/include/suitesparse
)

find_library(CHOLMOD_LIBRARY cholmod
    HINTS ${SUITESPARSE_ROOT}/lib
    /usr/local/lib
)
find_library(AMD_LIBRARY amd
    HINTS ${SUITESPARSE_ROOT}/lib
    /usr/local/lib
)
find_library(COLAMD_LIBRARY colamd
    HINTS ${SUITESPARSE_ROOT}/lib
    /usr/local/lib
)
find_library(CAMD_LIBRARY camd
    HINTS ${SUITESPARSE_ROOT}/lib
    /usr/local/lib
)
find_library(CCOLAMD_LIBRARY ccolamd
    HINTS ${SUITESPARSE_ROOT}/lib
    /usr/local/lib
)
find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig
    HINTS ${SUITESPARSE_ROOT}/lib
    /usr/local/lib
)
find_library(METIS_LIBRARY metis
    HINTS ${SUITESPARSE_ROOT}/lib
    /opt/homebrew/opt/metis/lib
    /usr/local/lib
)

if(CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARY)
    set(SUITESPARSE_FOUND TRUE)
    message(STATUS "Found CHOLMOD: ${CHOLMOD_LIBRARY}")
    # Enable CHOLMOD support in Eigen
    add_definitions(-DEIGEN_USE_CHOLMOD)
    
    # On macOS, linking Accelerate framework is often needed for SuiteSparse dependencies
    if(APPLE)
        find_library(ACCELERATE_FRAMEWORK Accelerate)
        if(ACCELERATE_FRAMEWORK)
            message(STATUS "Found Accelerate framework: ${ACCELERATE_FRAMEWORK}")
            list(APPEND EXTRA_LIBS ${ACCELERATE_FRAMEWORK})
        endif()
    endif()
else()
    message(WARNING "CHOLMOD not found. SuiteSparse solver will not be available.")
endif()

include_directories(${EIGEN3_INCLUDE_DIR})
if(SUITESPARSE_FOUND)
    include_directories(${CHOLMOD_INCLUDE_DIR})
endif()
include_directories(src)

file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.h")

add_executable(main ${SOURCES} ${HEADERS})

if(SUITESPARSE_FOUND)
    target_link_libraries(main ${CHOLMOD_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY} ${CAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${SUITESPARSECONFIG_LIBRARY} ${EXTRA_LIBS})
    if(METIS_LIBRARY)
        target_link_libraries(main ${METIS_LIBRARY})
    endif()
endif()

