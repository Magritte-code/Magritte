# Set minimun required CMake version
cmake_minimum_required (VERSION 3.18)

# Set CMake version policy
cmake_policy (SET CMP0048 NEW)

# Set project title, version, description, and url
project ("Paracabs"
    VERSION      0.0.1
    DESCRIPTION  "Parallelisation and Acceleration abstractions."
    HOMEPAGE_URL "https://github.com/Magritte-code/Paracabs"
)

# Set build type to Release to set all optimisation flags
#set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_BUILD_TYPE Relwithdebinfo)
#set (CMAKE_BUILD_TYPE Debug)
#set (CUDA_NVCC_FLAGS " -g -G -x cu -O5") WARNING -G cfr issue:#112
#set (CUDA_NVCC_FLAGS " -g")
#set (CUDA_NVCC_FLAGS "-O5 -use_fast_math")

# Use C++11
set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          -std=c++11")
set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -std=c++11")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++11")

# set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          -std=c++17 -fsycl")
# set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -std=c++17 -fsycl")
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++17 -fsycl")

#set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          -std=c++11 -fsycl --gcc-toolchain=/cosma/local/gcc/7.3.0/")
#set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -std=c++11 -fsycl --gcc-toolchain=/cosma/local/gcc/7.3.0/")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++11 -fsycl --gcc-toolchain=/cosma/local/gcc/7.3.0/")


if (PARACABS_USE_OPENMP)
    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C                "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS          "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES      "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY   ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY  ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -latomic")
    else ()
        find_package (OpenMP REQUIRED)
        # Set compiler and linker flags
        set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX             "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS       "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES   "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY   ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY  ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -latomic")
    else ()
        find_package (OpenMP REQUIRED)
        # Set compiler and linker flags
        set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ()

if (PARACABS_USE_MPI)
    find_package (MPI REQUIRED)
    # Include headers
    include_directories (SYSTEM ${MPI_INCLUDE_PATH})
endif ()

if (PARACABS_USE_CUDA)
    find_package (CUDA REQUIRED)
endif ()


# Get current directory path
get_filename_component (CURRENT_DIRECTORY ${CMAKE_CURRENT_LIST_FILE} DIRECTORY)

# Write configuration file
configure_file (${CURRENT_DIRECTORY}/Paracabs/configure.hpp.in
                ${CURRENT_DIRECTORY}/Paracabs/configure.hpp   )

# Includes for pybind11
include_directories (dependencies/pybind11)
include_directories (dependencies/pybind11/include)
