# Set project title
project ("Magritte")

# Set version numbers
set (MAGRITTE_VERSION_MAJOR 1)
set (MAGRITTE_VERSION_MINOR 0)

# Set build type
#set (CMAKE_BUILD_TYPE Release)
set (CMAKE_BUILD_TYPE Debug)
#set (CUDA_NVCC_FLAGS " -g")
#set (CUDA_NVCC_FLAGS "-O5 -use_fast_math")

# Set minimun required CMake version
cmake_minimum_required (VERSION 3.11)

# Write all binary files to the bin directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Build options
option (PYTHON_IO        "Handling all io through python (slow)" OFF)
option (PYTHON_BINDINGS  "Python front-end?"                     OFF)
option (OMP_PARALLEL     "OpenMP parallelisation required?"      OFF)
option (MPI_PARALLEL     "MPI parallelisation required?"         OFF)
option (GPU_ACCELERATION "Build using the CUDA GPU solver"       OFF)

# Convert options to bools for configuration file (MUST BE A BETTER WAY!)
if    (PYTHON_IO)
    set (MAGRITTE_PYTHON_IO true)
else  (PYTHON_IO)
    set (MAGRITTE_PYTHON_IO false)
endif (PYTHON_IO)
if    (OMP_PARALLEL)
    set (MAGRITTE_OMP_PARALLEL true)
else  (OMP_PARALLEL)
    set (MAGRITTE_OMP_PARALLEL false)
endif (OMP_PARALLEL)
if    (MPI_PARALLEL)
    set (MAGRITTE_MPI_PARALLEL true)
else  (MPI_PARALLEL)
    set (MAGRITTE_MPI_PARALLEL false)
endif (MPI_PARALLEL)
if    (GPU_ACCELERATION)
    set (MAGRITTE_GPU_ACCELERATION true)
else  (GPU_ACCELERATION)
    set (MAGRITTE_GPU_ACCELERATION false)
endif (GPU_ACCELERATION)

# Write configuration file
configure_file (${CMAKE_SOURCE_DIR}/src/configure.hpp.in
                ${CMAKE_SOURCE_DIR}/src/configure.hpp   )

# Set Paracabs options
set (PARACABS_USE_MESSAGE_PASSING ${MAGRITTE_MPI_PARALLEL})
set (PARACABS_USE_MPI             ${MAGRITTE_MPI_PARALLEL})
set (PARACABS_USE_MULTI_THREADING ${MAGRITTE_OMP_PARALLEL})
set (PARACABS_USE_OPENMP          ${MAGRITTE_OMP_PARALLEL})
set (PARACABS_USE_ACCELERATOR     ${MAGRITTE_GPU_ACCELERATION})
set (PARACABS_USE_CUDA            ${MAGRITTE_GPU_ACCELERATION})
set (PARACABS_USE_SIMD            false)

# Include Paracabs directory
include             (~/Dropbox/GitHub/Paracabs/CMakeLists.txt)
include_directories (~/Dropbox/GitHub/Paracabs/Paracabs      )

# Link to pybind11 directory if required
if    (PYTHON_BINDINGS OR PYTHON_IO)
    add_subdirectory (dependencies/pybind11 pybind11)
endif (PYTHON_BINDINGS OR PYTHON_IO)

# Add source directory
add_subdirectory (src)

# Add tests directory
add_subdirectory (tests)