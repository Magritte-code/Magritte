# Set minimun required CMake version
cmake_minimum_required (VERSION 3.11)

# Set project title
project (Tests)

set (PARACABS_USE_MESSAGE_PASSING true)
set (PARACABS_USE_MPI             true)

set (PARACABS_USE_MULTI_THREADING true)
set (PARACABS_USE_OPENMP          true)


# Include Paracabs directory
include (../CMakeLists.txt)
include_directories (../)

# Set directory to put executables in
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY tests)

# Create executables
#add_executable (test_1 test_1.cpp)
add_executable (test_message_passing  test_message_passing.cpp)

if    (PARACABS_USE_MPI)
    target_link_libraries (test_message_passing ${MPI_C_LIBRARIES})
endif (PARACABS_USE_MPI)
