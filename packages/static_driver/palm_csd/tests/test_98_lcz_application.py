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

"""Run the LCZ test cases."""

from pathlib import Path
from typing import Tuple

import numpy as np
import pytest
import rasterio as rio
import rasterio.transform as riotf
from numpy import ma

from palm_csd.create_driver import create_driver
from palm_csd.geo_converter import GeoConverter
from palm_csd.lcz import LCZTypes
from tests.tools import add_to_stem, modify_configuration, modify_configuration_output, ncdf_equal

test_folder = Path("tests/98_lcz_application/")
test_folder_ref = test_folder / "output/"


@pytest.fixture
def configuration_dcep(request: pytest.FixtureRequest, tmp_path: Path) -> Tuple[Path, Path, Path]:
    """Generate a configuration file with domain_root.dcep True or False."""
    config_in = test_folder / "full_lcz.yml"
    config_out = tmp_path / f"full_lcz_dcep_{request.param}.yml"
    file_out = f"lcz_dcep_{request.param}"
    file_ref = test_folder_ref / f"lcz_dcep_{request.param}"

    to_set = [(["output", "file_out"], file_out)]
    to_set.extend([(["domain", "dcep"], request.param)])
    to_set.extend([(["input", "path"], str(tmp_path))])

    modify_configuration(config_in, config_out, to_set=to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    # Create input data
    # LCZ values
    lcz_index: ma.MaskedArray = ma.MaskedArray(np.full((4, 10), 17, dtype=np.uint8))
    lcz_index[1:3, 1:9] = np.arange(1, 17, 1, dtype=np.uint8).reshape((2, 8))
    lcz_types = LCZTypes("summer", True)
    lcz_rgb = lcz_types.lcz_index_to_rgb(lcz_index)

    lcz_output = tmp_path / "lcz.tif"
    with rio.open(
        lcz_output,
        "w",
        driver="GTiff",
        height=lcz_rgb.shape[1],
        width=lcz_rgb.shape[2],
        count=3,
        dtype=lcz_rgb.dtype,
        crs=GeoConverter.crs_wgs84,
        transform=riotf.from_origin(13.0, 52.0, 0.005, 0.005),
    ) as output_file:
        output_file.write(lcz_rgb)

    # orography values
    zt = ma.zeros((1, 80, 150))
    for iy, ix in np.ndindex(zt.shape[1:]):
        zt[0, iy, ix] += iy * 0.5 + ix * 1.0

    zt_output = tmp_path / "zt.tif"
    with rio.open(
        zt_output,
        "w",
        driver="GTiff",
        height=zt.shape[1],
        width=zt.shape[2],
        count=1,
        dtype=zt.dtype,
        crs=GeoConverter.crs_wgs84,
        transform=riotf.from_origin(12.99, 52.01, 0.0005, 0.0005),
    ) as output_file:
        output_file.write(zt)

    return config_out, tmp_path / file_out, file_ref


@pytest.mark.parametrize("configuration_dcep", [True, False], indirect=True)
@pytest.mark.usefixtures("config_counters")
def test_lcz_run(configuration_dcep: Tuple[Path, Path, Path]):
    """Run the LCZ test case and compare with correct output."""
    create_driver(configuration_dcep[0], verbose={"gis": True})

    output_root = add_to_stem(configuration_dcep[1], "_root")
    output_root_ref = add_to_stem(configuration_dcep[2], "_root")

    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=10,
    ), "Root driver does not comply with reference"
