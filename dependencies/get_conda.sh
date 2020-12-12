#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Check operating system and download appropriate installation file
case "`uname -s`" in
    Linux*)
        echo "Recognized Linux as OS (assuming x86_64), installing..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
    ;;
    Darwin*)
        echo "Recognized macOS as OS (assuming x86_64), installing..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda3.sh
    ;;
    *)
        echo "Could not recognize OS. Aborting."
        return
esac

# Run the install script
bash miniconda3.sh -b -p $DIR/miniconda3

# Remove the installer
rm miniconda3.sh

# Export path to conda
export PATH=$DIR/miniconda3/bin:$PATH

# Remove old program locations
hash -r

# Configure and update conda
conda config --set always_yes yes --set changeps1 no
conda update --yes conda

# Create conda environment
conda env create -f conda_env.yml

# Activate the magritte conda environment
source activate magritte
