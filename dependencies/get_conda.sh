#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

echo "WARNING: This installer assumes a Linux x84_64 system!"

# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda3.sh
# Set permissions
chmod +x miniconda3.sh
# Run the install script
./miniconda3.sh -b
# Remove the installer
rm miniconda3.sh
ls
# Configure conda
export PATH=$DIR/miniconda/bin:$PATH
echo $PATH
echo "1"
# source "$DIR/miniconda3/etc/profile.d/conda.sh"
echo "2"
# hash -r
echo "3"
# conda config --set always_yes yes --set changeps1 no
echo "4"
conda update --yes conda
echo "5"
# Print all conda info (for debugging)
# conda info -a
# echo "6"
