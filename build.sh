#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Enter build directory
mkdir build; cd build

# Get compilers
CC_FLAG=$(which dpcpp)
CXX_FLAG=$(which dpcpp)

# Run cmake
#CC=$CC_FLAG CXX=$CXX_FLAG \
cmake $DIR/tests

# Run make
make -j4