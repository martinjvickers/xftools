# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (sRNA_mismatch_aligner CXX)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBRARIES OFF)
set(CMAKE_EXE_LINKER_FLAGS "-static")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB REQUIRED)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ")

# Add executable and link against SeqAn dependencies.
add_executable (sRNA_mismatch_aligner sRNA_mismatch_aligner.cpp)
target_link_libraries (sRNA_mismatch_aligner ${ZLIB_LIBRARIES} ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
