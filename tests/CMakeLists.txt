# Set minimun required CMake version
cmake_minimum_required (VERSION 3.11)

# Set project title
project (Tests)

set (PARACABS_USE_MESSAGE_PASSING true)
set (PARACABS_USE_MPI             true)

set (PARACABS_USE_MULTI_THREADING true)
set (PARACABS_USE_OPENMP          true)

set (PARACABS_USE_ACCELERATOR     true)
set (PARACABS_USE_CUDA            true)

set (PARACABS_USE_SIMD            false)

find_package (CUDA REQUIRED)
set (CUDA_NVCC_FLAGS "-x")

# Include Paracabs directory
include (../CMakeLists.txt)
include_directories (../Paracabs)

# Set directory to put executables in
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY tests)

# Create executables
#     add_executable (test_message_passing test_message_passing.cpp)
#     add_executable (test_multi_threading test_multi_threading.cpp)
cuda_add_executable (test_accelerator     test_accelerator.cpp    )
cuda_add_executable (test_datatypes       test_datatypes.cpp      )

if    (PARACABS_USE_MPI)
#    target_link_libraries (test_message_passing ${MPI_C_LIBRARIES})
#    target_link_libraries (test_multi_threading ${MPI_C_LIBRARIES})
endif (PARACABS_USE_MPI)
