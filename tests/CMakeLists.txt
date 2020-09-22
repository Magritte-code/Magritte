# Include directory containing source
include_directories (../src/)
include_directories (../dependencies)

set (SOURCE_FILES
    test_raytracer.cpp
    test_multigrid.cpp
)


if    (OMP_PARALLEL)
    find_package (OpenMP)
endif (OMP_PARALLEL)

if    (GPU_ACCELERATION)
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (GPU_ACCELERATION)

add_executable        (test_raytracer test_raytracer.cpp)
target_link_libraries (test_raytracer Magritte)

add_executable        (test_multigrid test_multigrid.cpp)
target_link_libraries (test_multigrid Magritte)

if (OpenMP_CXX_FOUND)
    target_link_libraries (test_raytracer OpenMP::OpenMP_CXX)
endif()
