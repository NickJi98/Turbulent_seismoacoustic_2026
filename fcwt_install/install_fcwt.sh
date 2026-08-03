#!/bin/bash
#
# Install the fCWT package from source with the patched files in this directory.
# Automates the manual steps described in the top-level README.txt.
#
# Usage (after activating the conda environment from environment.yml):
#   bash fcwt_install/install_fcwt.sh
#
# Supported platforms: macOS (tested on Apple M3) and Linux (tested on HPC cluster).
# On Linux clusters you may need to load fftw first (e.g. `module load fftw/3.3.10`).

set -e

# Directory containing this script (and the patched files)
patch_dir="$(cd "$(dirname "$0")" && pwd)"

# Check required commands
for cmd in git pip python; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: '$cmd' not found. Activate the conda environment first."; exit 1; }
done

# Build in a temporary directory
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT

echo "==> Cloning fCWT into $build_dir"
git clone https://github.com/fastlib/fCWT.git "$build_dir/fCWT"

echo "==> Applying patched files from $patch_dir"
cp "$patch_dir/setup.py"        "$build_dir/fCWT/setup.py"
cp "$patch_dir/boilerplate.py"  "$build_dir/fCWT/src/fcwt/boilerplate.py"
cp "$patch_dir/fcwt.cpp"        "$build_dir/fCWT/src/fcwt/fcwt.cpp"

if [ "$(uname)" = "Darwin" ]; then
    echo "==> macOS detected: adding libomp.dylib"
    mkdir -p "$build_dir/fCWT/libs"
    cp "$patch_dir/libomp.dylib" "$build_dir/fCWT/libs/"
fi

echo "==> Installing fCWT"
pip install "$build_dir/fCWT"

echo "==> Testing import"
python -c "import fcwt; print('fCWT installed successfully')"
