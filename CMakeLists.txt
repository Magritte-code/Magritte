# Set project title
project ("Paracabs")

# Set build type to Release to set all optimisation flags
#set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_BUILD_TYPE Relwithdebinfo)
#set (CMAKE_BUILD_TYPE Debug)
#set (CUDA_NVCC_FLAGS " -g -G -x cu -O5") WARNING -G cfr issue:#112
#set (CUDA_NVCC_FLAGS " -g")
#set (CUDA_NVCC_FLAGS "-O5 -use_fast_math")

set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          -std=c++17 -fsycl")
set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -std=c++17 -fsycl")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++17 -fsycl")

#set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          -std=c++11 -fsycl --gcc-toolchain=/cosma/local/gcc/7.3.0/")
#set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -std=c++11 -fsycl --gcc-toolchain=/cosma/local/gcc/7.3.0/")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -std=c++11 -fsycl --gcc-toolchain=/cosma/local/gcc/7.3.0/")

#set (OpenMP_C_FLAGS   "-qopenmp")
#set (OpenMP_CXX_FLAGS "-qopenmp")

# Set minimun required CMake version
cmake_minimum_required (VERSION 3.11)

if    (PARACABS_USE_OPENMP)
    find_package (OpenMP REQUIRED)
    # Set compiler and linker flags
    set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (PARACABS_USE_OPENMP)
if    (PARACABS_USE_MPI)
    find_package (MPI REQUIRED)
    # Include headers
    include_directories (SYSTEM ${MPI_INCLUDE_PATH})
endif (PARACABS_USE_MPI)
if    (PARACABS_USE_CUDA)
    find_package (CUDA REQUIRED)
endif (PARACABS_USE_CUDA)

#include (CMakePackageConfigHelpers)
#
#write_basic_package_version_file (
#    "${PROJECT_NAME}ConfigVersion.cmake"
#    VERSION 0.0.1
#    COMPATIBILITY SameMajorVersion
#)
#
#configure_package_config_file (
#    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
#    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
#    INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
#)
#
#install (
#    EXPORT ${PROJECT_NAME}_Targets
#    FILE ${PROJECT_NAME}Targets.cmake
#    NAMESPACE ${PROJECT_NAME}::
#    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
#)
#
#install (
#    FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
#    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
#    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
#)
#
#install (
#    DIRECTORY ${PROJECT_SOURCE_DIR}/include/SI
#    DESTINATION include
#)

# Get current directory path
get_filename_component (CURRENT_DIRECTORY ${CMAKE_CURRENT_LIST_FILE} DIRECTORY)

# Write configuration file
configure_file (${CURRENT_DIRECTORY}/Paracabs/configure.hpp.in
                ${CURRENT_DIRECTORY}/Paracabs/configure.hpp   )

# Create Paracabs library
#add_library (${PROJECT_NAME} INTERFACE)

# Add include directory, containing the headers
#target_include_directories (${PROJECT_NAME} INTERFACE Paracabs/)
