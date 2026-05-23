#!/bin/bash
print_valid_arguments() {
  echo "Valid arguments are 'True' or 'False'."
}

if [ "$#" -ne 1 ]; then
    echo "Error: This script requires exactly one argument."
    print_valid_arguments
    exit 1
fi

case "$1" in
  True)
    files="single_trees_True_root"
    ;;
  False)
    files="single_trees_False_root"
    ;;
  *)
    echo "Error: Invalid argument."
    print_valid_arguments
    exit 1
    ;;
esac

# Find the most recent folder in ../tmp
latest_folder=$(ls -td ../tmp/*/ | head -n 1)

result_file="${latest_folder}/${files}"

if [ ! -f "${result_file}" ]; then
    echo "Error: ${result_file} file not found!"
    exit 1
fi
ncks -O -L9 -h "${result_file}" output/${files}
