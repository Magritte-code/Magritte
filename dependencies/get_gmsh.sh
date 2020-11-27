#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

echo "WARNING: This installer assumes a Linux x84_64 system!"

# Download gmsh
wget https://gmsh.info/bin/Linux/gmsh-4.6.0-Linux64.tgz
# Untar the file
tar -zxvf gmsh-4.6.0-Linux64
# Rename the folder
mv gmsh-4.6.0-Linux64 gmsh
# Remove the tar ball
rm gmsh-4.6.0-Linux64
