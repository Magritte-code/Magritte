# Write the python library to the magritte python package
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/magritte)

# Include directories
include_directories (../)
include_directories (../../dependencies)
include_directories (../../dependencies/pybind11/include)

set (SOURCE_FILES
    pybindings.cpp
    ../io/cpp/io_cpp_text.cpp
    ../io/python/io_python.cpp
    ../model/model.cpp
    ../model/parameters/parameters.cpp
    ../model/geometry/geometry.cpp
    ../model/geometry/points/points.cpp
    ../model/geometry/rays/rays.cpp
    ../model/geometry/boundary/boundary.cpp
#    ../model/chemistry/chemistry.cpp
#    ../model/chemistry/species/species.cpp
#    ../model/thermodynamics/thermodynamics.cpp
#    ../model/thermodynamics/temperature/temperature.cpp
#    ../model/thermodynamics/turbulence/turbulence.cpp
#    ../model/lines/lines.cpp
#    ../model/lines/lineProducingSpecies/lineProducingSpecies.cpp
#    ../model/lines/lineProducingSpecies/linedata/linedata.cpp
#    ../model/lines/lineProducingSpecies/quadrature/quadrature.cpp
    ../solver/solver.cpp
)

if    (MPI_PARALLEL)
  find_package (MPI REQUIRED)
  include_directories (SYSTEM ${MPI_INCLUDE_PATH})
endif (MPI_PARALLEL)
if    (OMP_PARALLEL)
  find_package (OpenMP REQUIRED)
  set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (OMP_PARALLEL)
if    (GPU_ACCELERATION)
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (GPU_ACCELERATION)


if    (GPU_ACCELERATION)
    # Add Python package separately
    find_package (PythonLibs REQUIRED)
    # Add python dirs separately
    include_directories (${PYTHON_INCLUDE_DIRS})
    # Add CUDA library
    cuda_add_library (core SHARED ${SOURCE_FILES})
    # Link Python libraries
    target_link_libraries (core ${PYTHON_LIBRARIES})
else  (GPU_ACCELERATION)
    # pybind11 command to create the python module in bin directory
    pybind11_add_module (core SHARED ${SOURCE_FILES})
    # Link MPI libraries
endif (GPU_ACCELERATION)

# Link MPI libraries
if    (MPI_PARALLEL)
    target_link_libraries (core PRIVATE ${MPI_C_LIBRARIES})
endif (MPI_PARALLEL)

# Set library properties
set_target_properties (core PROPERTIES PREFIX "")
set_target_properties (core PROPERTIES SUFFIX ".so")