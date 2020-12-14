#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Enter build directory
mkdir build; cd build

# Get Python interpreter
PYTHON_EXECUTABLE=$(which python)

# Get compilers
COMPILER_CC=$(which gcc)
COMPILER_CXX=$(which g++)

# Run cmake
CC=$COMPILER_CC                                     \
CXX=$COMPILER_CXX                                   \
cmake                                               \
  -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_EXECUTABLE   \
  -DOMP_PARALLEL=ON                                 \
  -DMPI_PARALLEL=OFF                                \
  -DGPU_ACCELERATION=OFF                            \
  # -DCMAKE_CXX_FLAGS=-pg                             \
  # -DCMAKE_EXE_LINKER_FLAGS=-pg                      \
  # -DCMAKE_SHARED_LINKER_FLAGS=-pg                   \
  $DIR

# Run make
# note: this is just such that i have 1 thread left on my pc
make -j3
