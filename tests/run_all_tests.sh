#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR

# Create directories to store results
mkdir results

# Unit tests
echo "Running unit tests..."

# Integration tests
echo "Running integration tests..."

cd $DIR/benchmarks/analytic
python3 all_constant_single_ray.py         nosave
python3 density_distribution_single_ray.py nosave

cd $DIR/benchmarks/numeric
python3 vanZadelhoff_1_1D.py               nosave
