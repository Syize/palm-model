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

"""Run netcdf_data tests."""

from pathlib import Path

import pytest
from netCDF4 import Dataset
from numpy import ma

from palm_csd.netcdf_data import NCDFDimension, NCDFVariable, remove_existing_file
from tests.tools import ncdf_equal


def test_variables(tmp_path: Path):
    """Test dimensions and variables."""
    # with defined values
    x_dim = NCDFDimension(
        name="x",
        values=ma.MaskedArray([1.0, 3.0, 5.0]),
        data_type="f4",
        standard_name="projection_x_coordinate",
        long_name="x",
        units="m",
    )
    assert x_dim.size == 3
    assert len(x_dim) == 3

    # without defined values
    y_dim = NCDFDimension(
        name="y",
        data_type="f4",
        standard_name="projection_y_coordinate",
        long_name="y",
        units="m",
    )
    with pytest.raises(ValueError):
        print(y_dim.size)
    with pytest.raises(ValueError):
        len(y_dim)

    # write to file
    to_file_dimension = tmp_path / "dimension.nc"
    remove_existing_file(to_file_dimension)
    nc_data = Dataset(to_file_dimension, "a", format="NETCDF4")
    x_dim.to_dataset(nc_data)
    # fail because no values
    with pytest.raises(ValueError):
        y_dim.to_dataset(nc_data)

    # assign values
    y_dim.values = ma.MaskedArray([2.0, 4.0, 6.0, 8.0])
    assert y_dim.size == 4
    assert len(y_dim) == 4
    y_dim.to_dataset(nc_data)

    nc_data.close()

    # create variable
    buildings_var = NCDFVariable(
        name="buildings_2d",
        dimensions=(y_dim, x_dim),
        data_type="f4",
        fill_value=-9999.0,
        long_name="buildings",
        units="m",
        lod=1,
        coordinates="E_UTM N_UTM lon lat",
        grid_mapping="crs",
    )
    to_file_variable = tmp_path / "variable.nc"
    remove_existing_file(to_file_variable)

    # no defined values
    with pytest.raises(ValueError):
        buildings_var.to_nc(file=to_file_dimension)

    buildings_var.values = ma.arange(1, 13, 1).reshape([4, 3])
    # no defined filename
    with pytest.raises(ValueError):
        buildings_var.to_nc()

    buildings_var.to_nc(file=to_file_variable)

    buildings_id = NCDFVariable(
        name="buildings_id",
        dimensions=(y_dim, x_dim),
        values=ma.MaskedArray(range(10, 22)).reshape([4, 3]),
        data_type="i",
        fill_value=9999,
        long_name="buildings id",
        units="",
        coordinates="E_UTM N_UTM lon lat",
        grid_mapping="crs",
        file=to_file_variable,
    )

    buildings_id.to_nc()

    # compare output files
    assert ncdf_equal("tests/01_netcdf_data/output/dimension.nc", to_file_dimension)
    assert ncdf_equal("tests/01_netcdf_data/output/variable.nc", to_file_variable)
