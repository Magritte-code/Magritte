#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

echo "WARNING: This installer assumes a Linux x84_64 system!"

# Download CMake
wget https://github.com/Kitware/CMake/releases/download/v3.18.4/cmake-3.18.4-Linux-x86_64.tar.gz
# Untar the file
tar -zxvf cmake-3.18.4-Linux-x86_64.tar.gz
# Rename the folder
mv cmake-3.18.4-Linux-x86_64 cmake
# Remove the tar ball
rm cmake-3.18.4-Linux-x86_64.tar.gz
