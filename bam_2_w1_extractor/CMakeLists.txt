# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (bam_2_w1_extractor CXX)

# Make static binary
set (CMAKE_CXX_FLAGS  "-static ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wl,--whole-archive -Wl,-no-whole-archive")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Load zlib as a dependency for SeqAn and BAM and fail if not found.
find_package (ZLIB REQUIRED)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (bam_2_w1_extractor bam_2_w1_extractor.cpp)
target_link_libraries (bam_2_w1_extractor ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
