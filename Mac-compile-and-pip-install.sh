#!/bin/bash

# I brew-installed gcc-15 and openmp, and then
# can use that in an isolated way with this script

eval "$(/opt/homebrew/bin/brew shellenv)"
export CC=/opt/homebrew/bin/gcc-15
export CXX=/opt/homebrew/bin/g++-15
export CFLAGS="-Xpreprocessor -fopenmp"
export CXXFLAGS="-Xpreprocessor -fopenmp"
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

python setup.py build_ext --inplace

pip install -e .

