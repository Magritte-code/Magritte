#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Remove everything
rm -rf Eigen
rm -rf pybind11
rm -rf Paracabs
rm -rf googletest

# Get latest Eigen (there is no stable release yet that works with latest CUDA)
wget https://github.com/eigenteam/eigen-git-mirror/archive/master.zip
# Extract only the Eigen headers
unzip master.zip 'eigen-git-mirror-master/Eigen/*'
# Rename folder
mv eigen-git-mirror-master/Eigen/ Eigen/
# Remove old folder
rm -r eigen-git-mirror-master/
# Remove zip file
rm master.zip

# Get pybind11
wget https://github.com/pybind/pybind11/archive/v2.2.4.tar.gz
# Extract whole directory
tar -zxvf v2.2.4.tar.gz
# Rename the folder
mv pybind11-2.2.4 pybind11
# Remove tar ball
rm v2.2.4.tar.gz

# Get Paracabs
wget https://github.com/FredDeCeuster/Paracabs/archive/master.zip
# Extract the whole directory
unzip master.zip
# Rename the folder
mv Paracabs-master Paracabs
# Remove zip file
rm master.zip

# Get googletest
wget https://github.com/google/googletest/archive/release-1.10.0.tar.gz
# Extract whole directory
tar -zxvf release-1.10.0.tar.gz
# Rename the folder
mv googletest-release-1.10.0 googletest
# Remove tar ball
rm release-1.10.0.tar.gz
