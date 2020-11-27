#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Create conda environment
conda env create -f conda_env.yml
# Activate the magritte conda environment
conda activate magritte
