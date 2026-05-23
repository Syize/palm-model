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
    files="lcz_dcep_True_root"
    ;;
  False)
    files="lcz_dcep_False_root"
    ;;
  *)
    echo "Error: Invalid argument."
    print_valid_arguments
    exit 1
    ;;
esac

for f in $files; do
    if [ ! -f "$f" ]; then
        echo "Error: $f file not found!"
        exit 1
    fi
    ncks -O -L9 -h "$f" output/"$f"
done
