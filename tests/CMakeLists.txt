# Set minimun required CMake version
cmake_minimum_required (VERSION 3.11)

# Set project title
project (Tests)

#set (CMAKE_BUILD_TYPE Release)
set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_CUDA_FLAGS " -g -G -x cu -O5") WARNING -G cfr issue:#112
#set (CMAKE_CUDA_FLAGS " -g")
#set (CMAKE_CUDA_FLAGS "-O5 -use_fast_math")
set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda --use_fast_math --expt-relaxed-constexpr")

set (PARACABS_USE_MESSAGE_PASSING true)
set (PARACABS_USE_MPI             true)

set (PARACABS_USE_MULTI_THREADING true)
set (PARACABS_USE_OPENMP          true)

set (PARACABS_USE_ACCELERATOR     false)
set (PARACABS_USE_CUDA            false)

set (PARACABS_USE_SIMD            false)

set (SOURCE_FILES
        test_datatypes.cpp
)


if    (PARACABS_USE_CUDA)
    enable_language (CUDA)
    # Specify the accelerator architecture (currently we don't)
    set (CMAKE_CUDA_ARCHITECTURES OFF)
    # Set source file language to CUDA
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (PARACABS_USE_CUDA)


# Include Paracabs directory
include             (../CMakeLists.txt)
include_directories (../Paracabs      )

# Set directory to put executables in
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY tests)

# Create executables
#      add_executable (test_message_passing test_message_passing.cpp)
#      add_executable (test_multi_threading test_multi_threading.cpp)
#      add_executable (test_datatypes       test_datatypes.cpp      )
#cuda_add_executable(test test.cpp)

#cuda_add_executable (test_accelerator test_accelerator.cu)
#cuda_add_executable (test_allocator   test_allocator.cu  )
add_executable (test_datatypes   test_datatypes.cpp  )
#cuda_add_executable (test             test.cu            )

if    (PARACABS_USE_MPI)
#    target_link_libraries (test_message_passing ${MPI_C_LIBRARIES})
#    target_link_libraries (test_multi_threading ${MPI_C_LIBRARIES})
#    target_link_libraries (test_accelerator     ${MPI_C_LIBRARIES})
#    target_link_libraries (test_datatypes       ${MPI_C_LIBRARIES})
endif (PARACABS_USE_MPI)
