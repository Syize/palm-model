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

"""Run the vegetation patch application with combinations of input data."""

import csv
import itertools
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np
import numpy.ma as ma
import numpy.ma.core as ma_core
import pytest
import rasterio as rio
from netCDF4 import Dataset  # type: ignore

from palm_csd.create_driver import create_driver
from tests.tools import add_to_stem, modify_configuration, modify_configuration_output

test_folder = Path("tests/97_vegetation_patch_application/")
test_folder_ref = test_folder / "output"

width = 6
height = 7
dz = 3.0


@pytest.fixture
def create_input_data(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Tuple[Path, Path, Dict[str, ma.MaskedArray]]:
    """Create input data and configuration for the vegetation patch application test.

    The data in input_data.csv is read. The palm_csd input data is stored as GeoTIFF, the other
    variables are used for result comparison. The configuration is adjusted accordingly.

    Args:
        request: Parametrized input data. request.param includes the data with request.param[0] as
            the overhanging_trees setting and request.param[1] as the high_vegetation setting.
        tmp_path: Temporary path for storing input data and configuration.

    Returns:
        Configuration path, output file path, and variables for comparison.
    """
    # Run parameters
    overhanging_trees = request.param[0]
    replace_high_vegetation = request.param[1]
    estimate_lai_from_vegetation_height = request.param[2]

    # General strings for OH, RHV, ELH being True or False in input_data.csv
    setting_string_oh = f"OH{'T' if overhanging_trees else 'F'}"
    setting_string_rhv = f"RHV{'T' if replace_high_vegetation else 'F'}"
    setting_string_elh = f"ELH{'T' if estimate_lai_from_vegetation_height else 'F'}"
    setting_string = f"{setting_string_oh}_{setting_string_rhv}_{setting_string_elh}"

    # Variable-specific strings in input_data.csv
    patch_setting_string = f"patch_{setting_string_oh}_{setting_string_rhv}"
    vegetation_type_setting_string = f"vegetation_type_{setting_string_oh}_{setting_string_rhv}"
    lai_lsm_setting_string = f"lai_lsm_{setting_string_oh}_{setting_string_rhv}"
    lai_patch_setting_string = f"lai_patch_{setting_string_elh}"

    test_input = tmp_path / setting_string

    os.makedirs(test_input, exist_ok=True)

    # Read input data from CSV file.
    dtype = {
        "lai": np.float32,  # LAI input
        "patch_height": np.float32,  # Patch height input
        "patch_type": np.uint8,  # Patch type input
        "vegetation_type": np.uint8,  # Vegetation type input
        "pavement_type": np.uint8,  # Pavement type input
        patch_setting_string: np.bool_,  # Identified vegetation patches
        lai_patch_setting_string: np.float32,  # LAI used for vegetation patches
        "patch_height_filled": np.float32,  # Patch height used for vegetation patches
        "patch_type_filled": np.uint8,  # Patch type used for vegetation patches
        vegetation_type_setting_string: np.uint8,  # Output vegetation type
        lai_lsm_setting_string: np.float32,  # Output LAI for LSM
    }

    fillvalue = {
        "lai": -9999.0,
        "patch_height": -9999.0,
        "patch_type": 255,
        "vegetation_type": 255,
        "pavement_type": 255,
        patch_setting_string: False,
        lai_patch_setting_string: -9999.0,
        "patch_height_filled": -9999.0,
        "patch_type_filled": 255,
        vegetation_type_setting_string: 255,
        lai_lsm_setting_string: -9999.0,
    }

    variables: Dict[str, ma.MaskedArray] = {}
    for variable in dtype:
        variables[variable] = ma.empty(
            width * height, dtype=dtype[variable], fill_value=fillvalue[variable]
        )

    with open(test_folder / "input_data.csv") as input_data_csv:
        reader = csv.DictReader(input_data_csv)
        i = 0
        for row in reader:
            # Process each row. Each row is a dict of strings.
            # remove whitespaces
            row = {key.strip(): value.strip() for key, value in row.items()}
            for variable in variables:
                if variable == patch_setting_string:
                    variables[variable][i] = row[variable] == "True"
                else:
                    variables[variable][i] = (
                        ma.masked if row[variable] == "" else dtype[variable](row[variable])
                    )
            i += 1

    # Rename setting_string elements to general ones.
    variables["patch"] = variables.pop(patch_setting_string)
    variables["vegetation_type_adjusted"] = variables.pop(vegetation_type_setting_string)
    variables["lai_lsm"] = variables.pop(lai_lsm_setting_string)
    variables["lai_patch"] = variables.pop(lai_patch_setting_string)

    variables["zt"] = ma_core.zeros_like(variables["patch_height"])
    variables["zt"].mask = False

    for variable in variables:
        variables[variable] = variables[variable].reshape((height, width))

    # Variables needed for result checks.
    variables_return: Dict[str, ma.MaskedArray] = {
        key: ma.MaskedArray(np.flipud(variables[key]))
        for key in [
            "patch",
            "lai_patch",
            "lai_lsm",
            "patch_height_filled",
            "patch_type_filled",
            "vegetation_type_adjusted",
        ]
    }
    # Variables needed for palm_csd input.
    variables_save = {
        key: variables[key] for key in variables.keys() if key not in variables_return.keys()
    }

    for variable in variables_save:
        with rio.open(
            test_input / f"{variable}.tif",
            "w",
            driver="GTiff",
            height=height,
            width=width,
            dtype=variables[variable].dtype,
            nodata=variables[variable].fill_value,
            count=1,
        ) as output_file:
            output_file.write(variables[variable], 1)

    # Set up configuration.
    file_out = f"vegetation_patch_application_{setting_string}"

    to_set: List[Tuple[List[str], Union[str, float]]] = [(["input", "path"], str(test_input))]
    to_set.extend([(["output", "file_out"], file_out)])
    to_set.extend([(["domain", "dz"], dz)])
    to_set.extend([(["domain", "nx"], width - 1)])
    to_set.extend([(["domain", "ny"], height - 1)])
    to_set.extend([(["domain", "overhanging_trees"], overhanging_trees)])
    to_set.extend([(["domain", "replace_high_vegetation_types"], replace_high_vegetation)])
    to_set.extend(
        [(["domain", "estimate_lai_from_vegetation_height"], estimate_lai_from_vegetation_height)]
    )

    config_in = test_folder / "vegetation_patch.yml"
    config_out = tmp_path / f"vegetation_patch_{setting_string}.yml"

    modify_configuration(config_in, config_out, to_set=to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    return config_out, tmp_path / file_out, variables_return


overhanging_trees = [True, False]
replace_high_vegetation = [True, False]
estimate_lai_from_vegetation_height = [True, False]
combinations = list(
    itertools.product(
        overhanging_trees, replace_high_vegetation, estimate_lai_from_vegetation_height
    )
)
names = [f"OH {x} RHV {y} ELH {z}" for x, y, z in combinations]


@pytest.mark.parametrize(
    "create_input_data",
    combinations,
    ids=names,
    indirect=True,
)
@pytest.mark.usefixtures("config_counters")
def test_vegetation_patches(create_input_data: Tuple[Path, Path, Dict[str, ma.MaskedArray]]):
    """Test vegetation patches.

    Args:
        create_input_data: Configuration path, output file path, reference file path, and variables
        for comparison.
    """
    config = create_input_data[0]
    file_out = add_to_stem(create_input_data[1], "_root")
    variables = create_input_data[2]
    create_driver(config)

    with Dataset(file_out, "r") as nc_result:
        lai = nc_result.variables["vegetation_pars"][1, :, :]
        lad = nc_result.variables["lad"][:]
        vegetation_type = nc_result.variables["vegetation_type"][:]

    # Check if patches are identified correctly.
    assert np.array_equal(lad.any(axis=0).filled(False), variables["patch"])
    # Check if LAD sums up to LAI.
    assert ma.allclose(lad.sum(axis=0) * dz, variables["lai_patch"])

    # Check if the height is correct.
    # Find the highest LAD output index of dimension 0 that is not masked.
    highest_indices = np.argmax(~lad.mask[::-1, :, :], axis=0)
    highest_indices = ma.masked_where(
        ~np.any(~lad.mask, axis=0), lad.shape[0] - 1 - highest_indices
    )
    # Highest expected index according to patch height and patch height rules.
    highest_indices_expected = np.ceil(variables["patch_height_filled"] / dz).astype(int)
    assert ma.allclose(highest_indices, highest_indices_expected)

    # Check if the vegetation type and LAI are correct.
    assert ma.allequal(vegetation_type, variables["vegetation_type_adjusted"])
    assert np.array_equal(vegetation_type.mask, variables["vegetation_type_adjusted"].mask)
    assert ma.allequal(lai, variables["lai_lsm"])
    assert np.array_equal(lai.mask, variables["lai_lsm"].mask)
