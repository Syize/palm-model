#!/usr/bin/env python3
# ------------------------------------------------------------------------------ #
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
# ------------------------------------------------------------------------------ #
#
# Description:
# ------------
# Create test data set to test installation.
#
# ------------------------------------------------------------------------------ #

"""Executable to create a test data set for installation test."""

import numpy as np
from netCDF4 import Dataset


def main():
    """Create a test data set for the installation test."""
    # Settings
    file_name = "test_data.nc"

    dimension_length = 10
    fill_value = -9999

    dimension_name_x = "x"
    dimension_name_y = "y"
    variable_name = "Band1"

    # Create file
    nc_file = Dataset(file_name, "w", format="NETCDF4")

    # Create dimensions
    nc_file.createDimension(dimension_name_x, dimension_length)
    nc_file.createDimension(dimension_name_y, dimension_length)

    x = nc_file.createVariable(dimension_name_x, "f4", (dimension_name_x,))
    y = nc_file.createVariable(dimension_name_y, "f4", (dimension_name_y,))

    x[:] = np.arange(0, dimension_length, 1)
    y[:] = np.arange(0, dimension_length, 1)

    # Create empty data array
    data = nc_file.createVariable(
        variable_name, "f4", (dimension_name_y, dimension_name_x), fill_value=fill_value
    )
    data[:, :] = data._FillValue
    data.grid_mapping = "crs"

    # Create grid-mapping variable
    crs = nc_file.createVariable("crs", "i1")
    crs.spatial_ref = ""
    crs.grid_mapping_name = "transverse_mercator"
    crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257222101
    crs.longitude_of_prime_meridian = 0.0
    crs.longitude_of_central_meridian = 15.0
    crs.scale_factor_at_central_meridian = 0.9996
    crs.latitude_of_projection_origin = 0.0
    crs.false_easting = 500000.0
    crs.false_northing = 0.0
    crs.epsg_code = "EPSG:25833"


if __name__ == "__main__":
    main()
