#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Create directories to store models and results
mkdir models
mkdir results

# Unit tests
echo "Running unit tests..."

# Integration tests
echo "Running integration tests..."

python3 benchmarks/analytic/all_constant_1D.py
python3 benchmarks/analytic/density_distribution_1D.py
python3 benchmarks/numeric/vanZadelhoff_1_1D.py
