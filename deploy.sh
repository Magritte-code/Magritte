#! /bin/bash

# Get directory this script is in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export PYTHONDONTWRITEBYTECODE=0

# Build package for conda (using meta.yaml)
conda-build --channel conda-forge --channel anaconda $DIR 