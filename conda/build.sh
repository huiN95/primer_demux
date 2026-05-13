#!/bin/bash

set -ex

# Set Rust compilation environment variables
export CARGO_HOME=$BUILD_PREFIX/cargo
mkdir -p $CARGO_HOME

# Force linking to external htslib provided by Conda
export HTSLIB_INCLUDE_DIR=$PREFIX/include
export HTSLIB_LIBRARY_DIR=$PREFIX/lib

# Compile and install to the bin directory of the Conda environment
cargo install --locked --root $PREFIX --path .

# Clean up unnecessary cargo metadata
rm -f $PREFIX/.crates.toml $PREFIX/.crates2.json
