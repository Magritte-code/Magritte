#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Enter build directory
mkdir $DIR/build; cd $DIR/build

# Get Python interpreter
PYTHON_EXECUTABLE=$(which python)

# Get compilers
# COMPILER_CC=$(which gcc)
# COMPILER_CXX=$(which g++)
# if the compilers are not set, use the default gcc and g++ compilers
COMPILER_CC=${CC:-$(which gcc)}
COMPILER_CXX=${CXX:-$(which g++)}
echo "C compiler: " $COMPILER_CC
# COMPILER_CC=$(which gcc-11)
# COMPILER_CXX=$(which g++-11)
#COMPILER_CC=$(which clang)
#COMPILER_CXX=$(which clang++)

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
