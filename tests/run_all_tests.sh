#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Go to this directory
cd $DIR


mkdir models
mkdir results

# Unit tests
echo "Running unit tests..."

# Integration tests
echo "Running integration tests..."

python benchmarks/analytic/all_constant_1D.py
python benchmarks/numeric/vanZadelhoff_1_1D.py
