#!/bin/bash
#
# This script prepares CM1 runs with dense 1 Hz outputs,
# after the main runs have finished.
# For each parameter pair, it:
#   1. Creates a directory named V{V}_Cd{Cd×1000, zero-padded}_dense.
#   2. Copies the base input files (input_sounding and namelist.input_dense).
#   3. Edits key parameters (hurr_vg, cnstce, cnstcd) in the copied namelist.
#   4. Links the restart files from the corresponding main run directory.

# Parameter ranges to search
V_arr=$(seq 38 2 46)
Cd_arr=$(seq 0.018 0.002 0.026)  

# Input files to be copied
file1="input_sounding"
file2="namelist.input_dense"

# Check input file
if [[ ! -f $file1 || ! -f $file2 ]]; then
    echo "Error: Input files '$file1' and '$file2' not found."
    exit 1
fi

# Loop over parameter combinations
for Cd in $Cd_arr; do
    for V in $V_arr; do

        # Create directory
        dir_name="V${V}_Cd$(printf "%02d" $(echo "$Cd * 1000 / 1" | bc))_dense"
        dir1="V${V}_Cd$(printf "%02d" $(echo "$Cd * 1000 / 1" | bc))"
        mkdir -p "$dir_name"

        # Copy input files
        cp "$file1" "$dir_name/input_sounding"

        # Modify input parameters
        sed -e "s/hurr_vg\s*=\s*[0-9.]*,/hurr_vg       =      $(printf "%.1f" $V),/" \
            -e "s/cnstce\s*=\s*[0-9.]*,/cnstce     =   ${Cd},/" \
            -e "s/cnstcd\s*=\s*[0-9.]*,/cnstcd     =   ${Cd},/" \
            "$file2" > "${dir_name}/namelist.input"
            
        # Restart files
        for fl in "$dir1"/cm1rst_000006*; do
            [ -e "$fl" ] || continue
            ln -fs "$(realpath "$fl")" "$dir_name/"
        done
    done
done
