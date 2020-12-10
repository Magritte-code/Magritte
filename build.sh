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
  $DIR

# Run make
make -j4

# Go to Magritte root directory
cd $DIR

echo "Trying to uninstall magritte (to avoid duplication)..."
pip uninstall magritte

echo "Installing magritte python package..."
python setup.py install
