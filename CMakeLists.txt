# Ensure that cmake has fetch_content commands with minimum of 3.13
CMAKE_MINIMUM_REQUIRED( VERSION 3.13 )

if(NOT DEFINED CPPDATAPIPELINEREF)
    SET(CPPDATAPIPELINEREF "heads/main")
endif()

project(BEEPmbp
    VERSION 0.1.0
    DESCRIPTION "Bayesian Estimation of Epidemic Parameters using Model-Based Proposals"
    HOMEPAGE_URL "https://github.com/ScottishCovidResponse/BEEPmbp"
    LANGUAGES CXX C
)

set (BEEPMBP BEEPmbp)

# Set Output Directories to avoid issues on multi release compilers
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set C++ Standard to 11
SET( CMAKE_CXX_STANDARD 11 )
SET( CMAKE_CXX_STANDARD_REQUIRED ON )

# Set Compiler flag to add available instructions (for SIMD)
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g  -mfma -O3 -std=c++11 -fmax-errors=3" )

if (MSVC)
    # supres warnings
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /w")
else()
    # supress warnings
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
endif()

# Require MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Add the src dirctory compiling it's CMakeLists.txt
ADD_SUBDIRECTORY( src )