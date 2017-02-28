# Minimummapper.cpp cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (mapper CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

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
add_executable (mapper mapper.cpp)
target_link_libraries (mapper ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (build_reference build_reference.cpp)
target_link_libraries (build_reference ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (mapper2 mapper2.cpp)
target_link_libraries (mapper2 ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (bamwrite_test bamwrite_test.cpp)
target_link_libraries (bamwrite_test ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
