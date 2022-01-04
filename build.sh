#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Enter the directory this script is in
cd $DIR

# echo "Trying to uninstall magritte (to avoid duplication)..."
pip uninstall magritte

# echo "Installing magritte python package..."
python setup.py install
