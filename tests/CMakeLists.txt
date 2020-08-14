# Include directory containing source
include_directories (../src/)
include_directories (../dependencies)

find_package (CUDA REQUIRED)
cuda_add_executable   (test_raytracer test_raytracer.cu)
target_link_libraries (test_raytracer Magritte)
