#!/bin/bash

print_valid_arguments() {
  echo "Valid arguments are 'root', 'N02', 'root_np', 'N02_np', 'root_nt', 'N02_nt', 'trees_N02', 'rootgr', 'N02gr'."
}

if [ "$#" -ne 1 ]; then
    echo "Error: This script requires exactly one argument."
    print_valid_arguments
    exit 1
fi

diffveg=false
case "$1" in
  root)
    files="berlin_tiergarten_root"
    ;;
  N02)
    files="berlin_tiergarten_N02"
    ;;
  root_np)
    files="berlin_tiergarten_no_patches_root"
    diffveg=true
    ;;
  N02_np)
    files="berlin_tiergarten_no_patches_N02"
    diffveg=true
    ;;
  root_nt)
    files="berlin_tiergarten_no_trees_root"
    diffveg=true
    ;;
  N02_nt)
    files="berlin_tiergarten_no_trees_N02"
    diffveg=true
    ;;
  trees_N02)
    files="berlin_tiergarten_trees_N02"
    diffveg=true
    ;;
  rootgr)
    files="berlin_tiergarten_geo_referenced_root"
    ;;
  N02gr)
    files="berlin_tiergarten_geo_referenced_N02"
    ;;
  *)
    echo "Error: Invalid argument."
    print_valid_arguments
    exit 1
    ;;
esac

# Find the most recent folder in ../tmp
latest_folder=$(ls -td ../tmp/*/ | head -n 1)

if [ "$diffveg" = true ]; then

  for f in $files; do
    if [ ! -f "$latest_folder/$f" ]; then
        echo "Error: $f file not found!"
        exit 1
    fi
    ncks -O -L9 -v lad,bad,tree_id,tree_type,vegetation_pars,vegetation_type,zlad "$latest_folder/$f" output/diff_$f
  done

else

  for f in $files; do
      if [ ! -f "$latest_folder/$f" ]; then
          echo "Error: $f file not found!"
          exit 1
      fi
      ncks -O -L9 -h "$latest_folder/$f" output/$f
  done

fi
