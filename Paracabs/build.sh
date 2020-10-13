#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Enter build directory
mkdir build; cd build

# Get compilers
COMPILER_CC=$(which gcc)
COMPILER_CXX=$(which g++)

# Run cmake
CC=$COMPILER_CC   \
CXX=$COMPILER_CXX \
cmake $DIR/tests

# Run make
make -j4
