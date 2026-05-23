# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2025  Leibniz Universitaet Hannover
# Copyright 2022-2025  Technische Universitaet Berlin

"""Tools for testing."""

from pathlib import Path
from pprint import pprint
from typing import Dict, List, Optional, Sequence, Tuple, Union

import rasterio as rio
import yaml
from deepdiff import DeepDiff
from netCDF4 import Dataset
from numpy import ma


def ncdf_equal(
    ncdf_ref: Union[Path, str],
    ncdf_com: Union[Path, str],
    check_metadata: bool = True,
    check_fields: bool = True,
    metadata_significant_digits: Optional[int] = None,
    metadata_exclude_regex_paths: Optional[List[str]] = None,
    fields_only: Optional[List[str]] = None,
    fields_exclude: Optional[List[str]] = None,
) -> bool:
    """Test if two netCDF files are equal in terms of values and metadata."""
    nc_data_ref = Dataset(ncdf_ref, "r")
    nc_data_com = Dataset(ncdf_com, "r")

    diff = {}
    if check_metadata:
        # metadata diff, a Dataset does not yet load the values of the fields
        # exclude "_grpid" and "_varid" because they are expected to differ and
        #  do not have an impact on the content of a netCDF file
        exclude = ["_grpid", "_varid", "_dimid", "_grp"]
        if metadata_exclude_regex_paths is not None:
            exclude.extend(metadata_exclude_regex_paths)
        if fields_exclude is not None:
            exclude.extend(fields_exclude)

        metadata_diff = DeepDiff(
            nc_data_ref.__dict__ | nc_data_ref.variables | nc_data_ref.dimensions,
            nc_data_com.__dict__ | nc_data_com.variables | nc_data_com.dimensions,
            exclude_regex_paths=exclude,
            significant_digits=metadata_significant_digits,
            number_format_notation="e",
        )
        diff.update(metadata_diff.to_dict())

    if check_fields:
        # value diff
        value_diff = {}
        keys_ref = nc_data_ref.variables.keys()
        keys_com = nc_data_com.variables.keys()

        if fields_only is None and keys_ref != keys_com:
            diff.update({"Dataset": "variables differ"})

        for variable in keys_ref & keys_com:
            if fields_only is not None and variable not in fields_only:
                continue
            if fields_exclude is not None and variable in fields_exclude:
                continue

            field_ref = nc_data_ref.variables[variable][:]
            field_com = nc_data_com.variables[variable][:]
            if field_ref.dtype == object or field_com.dtype == object:
                if (field_ref != field_com).any():
                    value_diff.update({variable: "field values differ"})
            else:
                if not ma.allclose(field_ref, field_com):
                    value_diff.update({variable: "field values differ"})
            if not ma.allequal(ma.getmaskarray(field_ref), ma.getmaskarray(field_com)):
                value_diff.update({variable: "field masks differ"})
        diff.update(value_diff)

    nc_data_ref.close()
    nc_data_com.close()

    if diff:
        pprint(diff)
        return False

    return True


def geotiff_equal(geotiff_ref: Union[Path, str], geotiff_com: Union[Path, str]) -> bool:
    """Compare two GeoTIFF files in terms of metadata and data arrays.

    Args:
        geotiff_ref: Reference GeoTIFF file.
        geotiff_com: Comparison GeoTIFF file.

    Returns:
        True if the GeoTIFF files are equal, False otherwise.
    """
    with rio.open(geotiff_ref) as ref, rio.open(geotiff_com) as com:
        # Compare metadata
        if ref.meta != com.meta:
            print(f"Metadata of {geotiff_ref} and {geotiff_com} do not match.")
            return False

        # Read data
        data_ref = ref.read(1, masked=True)
        data_com = com.read(1, masked=True)

        # Compare data arrays
        if not ma.allclose(data_ref, data_com) or not ma.allequal(
            ma.getmaskarray(data_ref), ma.getmaskarray(data_com)
        ):
            print(f"Data arrays of {geotiff_ref} and {geotiff_com} do not match.")
            return False

    return True


def modify_configuration(
    config_in: Union[str, Path],
    config_out: Optional[Union[str, Path]] = None,
    to_delete: Optional[Sequence[Sequence[str]]] = None,
    to_set: Optional[Sequence[Tuple[Sequence[str], Union[str, float]]]] = None,
    to_replace: Optional[Sequence[Tuple[Sequence[str], str, str]]] = None,
    to_rename: Optional[Sequence[Tuple[Sequence[str], str]]] = None,
) -> Dict:
    """Update a configuration YAML with deleted, set, replaced and renamed configuration elements.

    If config_out is not None, the resulting dictionary is stored in a file.

    Args:
        config_in: Input configuration file.
        config_out: Output configuration file. Defaults to None.
        to_delete: Elements to delete. Defaults to None.
        to_set: Elements to set. Defaults to None.
        to_replace: Elements to replaces. Defaults to None.
        to_rename: Elements to rename. Defaults to None.

    Returns:
        Updated dictionary.
    """
    with open(config_in, "r", encoding="utf-8") as file:
        complete_dict = yaml.safe_load(file)

    def nested_del(dic: Dict, keys: Sequence[str]):
        for key in keys[:-1]:
            dic = dic[key]
        del dic[keys[-1]]

    def nested_set(dic: Dict, keys: Sequence[str], value: Union[str, float]):
        for key in keys[:-1]:
            dic = dic.setdefault(key, {})
        dic[keys[-1]] = value

    def nested_replace(dic: Dict, keys: Sequence[str], value: str, replace: str):
        # def update(dic: Any, value: str, replace: str):
        #     if isinstance(dic, dict):
        #         for key in dic:
        #             update(dic[key], value, replace)
        #     elif isinstance(dic, str):
        #         dic = dic.replace(value, replace)

        if len(keys) > 1:
            for key in keys[:-1]:
                dic = dic[key]

        if isinstance(dic[keys[-1]], dict):
            for k in dic[keys[-1]]:
                # if isinstance(v, str):
                #     dic[keys[-1]][k] = v.replace(value, replace)
                # elif isinstance(v, dict):
                nested_replace(dic[keys[-1]], [k], value, replace)
        elif isinstance(dic[keys[-1]], str):
            dic[keys[-1]] = dic[keys[-1]].replace(value, replace)

    def nested_rename(dic: Dict, keys: Sequence[str], value: str):
        for key in keys[:-1]:
            dic = dic[key]
        dic[value] = dic[keys[-1]]
        del dic[keys[-1]]

    if to_delete is not None:
        for td in to_delete:
            nested_del(complete_dict, td)
    if to_set is not None:
        for ts in to_set:
            nested_set(complete_dict, ts[0], ts[1])
    if to_replace is not None:
        for tr in to_replace:
            nested_replace(complete_dict, tr[0], tr[1], tr[2])
    if to_rename is not None:
        for trn in to_rename:
            nested_rename(complete_dict, trn[0], trn[1])

    if config_out is not None:
        with open(config_out, "w", encoding="utf-8") as file:
            yaml.dump(complete_dict, file)

    return complete_dict


def modify_configuration_output(
    config_in: Union[str, Path], config_out: Union[str, Path], output_path: Path
) -> Dict:
    """Update a configuration YAML for output files.

    Args:
        config_in: Input configuration file.
        config_out: Output configuration file.
        output_path: Path for output files.

    Returns:
        Updated dictionary.
    """
    return modify_configuration(
        config_in, config_out, to_set=[(["output", "path"], str(output_path))]
    )


def add_to_stem(file_path: Path, addition: str) -> Path:
    """Add a string to the stem of a file path.

    Args:
        file_path: Original file path.
        addition: String to add to the stem.

    Returns:
        New file path with the added string in the stem.
    """
    return file_path.with_stem(file_path.stem + addition)


def add_root_n02(file_path: Path) -> Tuple[Path, Path]:
    """Add '_root' and '_N02' to the stem of a file path.

    Args:
        file_path: Original file path.

    Returns:
        Tuple containing the new file path with '_root' and '_N02' added to the stem.
    """
    return (
        add_to_stem(file_path, "_root"),
        add_to_stem(file_path, "_N02"),
    )
