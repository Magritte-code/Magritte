# Include directories
include_directories (.)
include_directories (../dependencies)
include_directories (../dependencies/eigen)
include_directories (../dependencies/pybind11)
include_directories (../dependencies/pybind11/include)

set (SOURCE_FILES
    io/cpp/io_cpp_text.cpp
    model/model.cpp
    model/parameters/parameters.cpp
    model/geometry/geometry.cpp
    model/geometry/points/points.cpp
    model/geometry/rays/rays.cpp
    model/geometry/boundary/boundary.cpp
    model/chemistry/chemistry.cpp
    model/chemistry/species/species.cpp
    model/thermodynamics/thermodynamics.cpp
    model/thermodynamics/temperature/temperature.cpp
    model/thermodynamics/turbulence/turbulence.cpp
    model/lines/lines.cpp
    model/lines/lineProducingSpecies/lineProducingSpecies.cpp
    model/lines/lineProducingSpecies/linedata/linedata.cpp
    model/lines/lineProducingSpecies/linedata/collisionPartner/collisionPartner.cpp
    model/lines/lineProducingSpecies/quadrature/quadrature.cpp
    model/radiation/radiation.cpp
    model/radiation/frequencies/frequencies.cpp
    model/image/image.cpp
    solver/solver.cpp
)

if    (MPI_PARALLEL)
    find_package (MPI REQUIRED)
    include_directories (SYSTEM ${MPI_INCLUDE_PATH})
endif (MPI_PARALLEL)
if    (OMP_PARALLEL)
    find_package (OpenMP REQUIRED)
endif (OMP_PARALLEL)
if    (GPU_ACCELERATION)
    set_source_files_properties(${SOURCE_FILES}         PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(io/python/io_python.cpp PROPERTIES LANGUAGE CUDA)
endif (GPU_ACCELERATION)

# Add python dirs separately
include_directories (${PYTHON_INCLUDE_DIRS})

# Create Magritte library
add_library (Magritte ${SOURCE_FILES})

if    (MPI_PARALLEL)
    target_link_libraries (Magritte ${MPI_C_LIBRARIES})
endif (MPI_PARALLEL)

if    (PYTHON_IO)
    # Create library for python io
    add_library (PyIo io/python/io_python.cpp)
    # Link embedding and link to Magritte
    target_link_libraries (PyIo PRIVATE pybind11::embed)
    target_link_libraries (Magritte PyIo)
endif (PYTHON_IO)

if    (PYTHON_BINDINGS)
    add_subdirectory (bindings)
endif (PYTHON_BINDINGS)
