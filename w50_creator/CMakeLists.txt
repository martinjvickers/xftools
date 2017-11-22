# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (w50_creator_trim CXX)

set(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage -static ${CMAKE_CXX_FLAGS} -Wl,--whole-archive -Wl,-no-whole-archive")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
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
add_definitions ( -DSEQAN_DISABLE_VERSION_CHECK )

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ")

# Add executable and link against SeqAn dependencies.
add_executable (w50_tests w50_tests.cpp w50_functions.cpp)
target_link_libraries (w50_tests ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (w50_creator w50_creator.cpp w50_functions.cpp)
target_link_libraries (w50_creator ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

