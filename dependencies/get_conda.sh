#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

echo "WARNING: This installer assumes a Linux x84_64 system!"

# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda.sh
# Run the install script
bash miniconda.sh -b -p $DIR/miniconda3
# Remove the installer
rm miniconda.sh
# Configure conda
source "$DIR/miniconda3/etc/profile.d/conda.sh"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
# Print all conda info (for debugging)
conda info -a
