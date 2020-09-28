# Include directory containing source
include_directories (../src/)
include_directories (../dependencies)

set (SOURCE_FILES
    test_raytracer.cpp
    test_multigrid.cpp
    test_0th_short.cpp
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

add_executable        (test_0th_short test_0th_short.cpp)
target_link_libraries (test_0th_short Magritte)

if (OpenMP_CXX_FOUND)
    target_link_libraries (test_raytracer OpenMP::OpenMP_CXX)
    target_link_libraries (test_multigrid OpenMP::OpenMP_CXX)
    target_link_libraries (test_0th_short OpenMP::OpenMP_CXX)
endif()
