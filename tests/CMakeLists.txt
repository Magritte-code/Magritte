# Include directory containing source
include_directories (../src/)
include_directories (../dependencies)

set (SOURCE_FILES
    test_raytracer.cpp
)

if    (GPU_ACCELERATION)
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (GPU_ACCELERATION)

add_executable        (test_raytracer test_raytracer.cpp)
target_link_libraries (test_raytracer Magritte)
