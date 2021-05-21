# Set minimun required CMake version
cmake_minimum_required (VERSION 3.18)

# Set CMake version policy
cmake_policy (SET CMP0048 NEW)
cmake_policy (SET CMP0104 OLD)

# Set project title, version, description, and url
project ("Magritte"
    VERSION      0.0.1
    DESCRIPTION  "A modern software library for simulating radiation transport."
    HOMEPAGE_URL "https://github.com/Magritte-code/Magritte"
)

# Set build type
set (CMAKE_BUILD_TYPE Release)
# set (CMAKE_BUILD_TYPE Debug)

# Write all binary files to the bin directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Use C++11
set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          -std=c++11")
set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -std=c++11")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++11")

# Build options
option (PYTHON_IO        "Handling all io through python (slow)" ON)
option (PYTHON_BINDINGS  "Python front end?"                     ON)
option (OMP_PARALLEL     "OpenMP parallelisation required?"      ON)
option (MPI_PARALLEL     "MPI parallelisation required?"         OFF)
option (GPU_ACCELERATION "Use the GPU solver"                    OFF)
option (GPU_CUDA         "Use Paracabs CUDA implementation"      OFF)
option (GPU_SYCL         "Usa Paracabs SYCL implementation"      OFF)

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
    if    (GPU_SYCL)
        set (MAGRITTE_GPU_SYCL true)
        set (MAGRITTE_GPU_CUDA false)
    else  (GPU_SYCL)
        set (MAGRITTE_GPU_SYCL false)
        set (MAGRITTE_GPU_CUDA true)
    endif (GPU_SYCL)
else  (GPU_ACCELERATION)
    set (MAGRITTE_GPU_ACCELERATION false)
    set (MAGRITTE_GPU_CUDA         false)
    set (MAGRITTE_GPU_SYCL         false)
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
set (PARACABS_USE_CUDA            ${MAGRITTE_GPU_CUDA})
set (PARACABS_USE_SYCL            ${MAGRITTE_GPU_SYCL})
set (PARACABS_USE_SIMD            false)

# Enable CUDA if required
if (GPU_ACCELERATION)
    enable_language (CUDA)
    # Specify the accelerator architecture (currently we don't)
	set (CMAKE_CUDA_ARCHITECTURES OFF)
    set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda --use_fast_math --expt-relaxed-constexpr")
endif ()

# Include Paracabs directory
include             (dependencies/Paracabs/CMakeLists.txt)
include_directories (dependencies/Paracabs/Paracabs      )

# Link to pybind11 directory if required
if (PYTHON_BINDINGS OR PYTHON_IO)
    add_subdirectory (dependencies/pybind11 pybind11)
endif ()

# Add source directory
add_subdirectory (src)

# Add tests directory
add_subdirectory (tests)
