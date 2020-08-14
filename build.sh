#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Enter build directory
mkdir build; cd build

# Get Python interpreter
PYTHON_EXECUTABLE=$(which python)

# Get compilers
CC_FLAG=$(which gcc)
CXX_FLAG=$(which g++)

# Run cmake
CC=$CC_FLAG CXX=$CXX_FLAG                         \
cmake                                             \
  -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE \
  -DPYTHON_IO=ON                                  \
  -DPYTHON_BINDINGS=ON                            \
  -DOMP_PARALLEL=ON                               \
  -DMPI_PARALLEL=OFF                              \
  -DGPU_ACCELERATION=OFF                          \
  $DIR

# Run make
make -j4