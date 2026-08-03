#!/bin/bash

# Target directory for gathering all files
main_dir=$(pwd)
database_dir="$main_dir/param_search"

# Create database directory if it doesn't exist
mkdir -p "$database_dir"

# Loop over each directory
for rs_dir in "$main_dir"/V*_Cd*[02468]; do
    # Skip if not a directory
    [[ -d "$rs_dir" ]] || continue

    # Extract the name (e.g., V*_Cd*)
    rs_name=$(basename "$rs_dir")

    # Create corresponding subdirectory in database
    target_dir="$database_dir/$rs_name"
    mkdir -p "$target_dir"

    echo "Gathering files from $rs_name ..."

    # Copy files from the root directory
    cp "$rs_dir/namelist.input" "$target_dir/"
    cp "$rs_dir/input_sounding" "$target_dir/"

    # Copy files from /Data/
    data_dir="$rs_dir/Data"
    if [[ -d "$data_dir" ]]; then
        cp "$data_dir"/* "$target_dir/"
    fi

done

echo "All files gathered under $database_dir"
