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

"""Module with objects that could be read or written to netCDF."""

import logging
import os
from pathlib import Path
from typing import ClassVar, List, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from netCDF4 import Dataset
from numpy import ma
from pydantic import BaseModel, field_validator

logger = logging.getLogger(__name__)


class NCDFDimension(BaseModel, arbitrary_types_allowed=True, validate_assignment=True):
    """A dimension that could written to netCDF with a corresponding dimension variable included.

    Its values can be stored in the attribute values or supplied when written.
    """

    name: str
    """Name."""
    data_type: Union[str, type]
    """Data type."""

    long_name: Optional[str] = None
    """Long name."""
    units: Optional[str] = None
    """Units."""
    standard_name: Optional[str] = None
    """Standard name."""
    values: Optional[npt.NDArray] = None
    """Values."""

    _metadata: ClassVar[List[str]] = ["long_name", "units", "standard_name"]
    """Metadata attributes."""

    @field_validator("data_type")
    @classmethod
    def _check_data_type(cls, v: Union[str, type]) -> Union[str, type]:
        """Check if the data type is a string or a type.

        Args:
            v: Data type.

        Raises:
            ValueError: Data type not valid.

        Returns:
            Validated data type.
        """
        if not (v is str or isinstance(v, str)):
            raise ValueError("Data type must be a string or str.")
        return v

    @property
    def size(self) -> int:
        """Get the size of the dimension.

        It is derived from the values attribute.

        Raises:
            ValueError: Values of dimension not defined.

        Returns:
            Size of the dimension.
        """
        if self.values is None:
            raise ValueError(f"Values of dimension {self.name} not defined.")
        return self.values.size

    def __len__(self) -> int:
        """Size of the dimension.

        Returns:
            Size of the dimension.
        """
        return self.size

    def to_dataset(self, nc_data: Dataset, values: Optional[npt.NDArray] = None) -> None:
        """Add dimension and its dimension variable to a netCDF Dataset if not already included.

        The values are either supplied in the function call or taken from the values
        attribute.

        Args:
            nc_data: NetCDF Dataset.
            values: Values to write. Defaults to None.

        Raises:
            ValueError: No data to write.
        """
        if self.name not in nc_data.dimensions:
            if values is not None:
                to_write = values
            elif self.values is not None:
                to_write = self.values
            else:
                raise ValueError(f"Values of dimension {self.name} not defined when writing.")

            logger.debug(f"Adding dimension {self.name}.")

            nc_data.createDimension(self.name, len(to_write))

            nc_var = nc_data.createVariable(self.name, self.data_type, self.name)
            nc_var[:] = to_write

            for attr in self._metadata:
                attr_value = getattr(self, attr)
                if attr_value is not None:
                    nc_var.setncattr(attr, attr_value)


class NCDFVariable(BaseModel, arbitrary_types_allowed=True, validate_assignment=True):
    """A variable that could written to netCDF.

    It includes its metadata and dimensions. Its values can be stored in the attribute values or
    supplied when written. A default file for input/output can be also stored.
    """

    name: str
    """Name."""
    dimensions: Tuple[NCDFDimension, ...]
    """Dimensions."""
    data_type: str
    """Data type."""
    fill_value: float
    """Fill value."""
    long_name: str
    """Long name."""
    units: str
    """Units."""

    values: Optional[ma.MaskedArray] = None
    """Values."""

    coordinates: Optional[str] = None
    """Coordinates."""
    grid_mapping: Optional[str] = None
    """Grid mapping."""
    lod: Optional[int] = None
    """Level of detail."""
    standard_name: Optional[str] = None
    """Standard name."""

    mandatory: bool = True
    """Whether the variable is mandatory in the static driver."""

    file: Optional[Path] = None
    """Default file for input/output."""

    _metadata: ClassVar[List[str]] = [
        "long_name",
        "units",
        "standard_name",
        "lod",
        "coordinates",
        "grid_mapping",
    ]
    """Metadata attributes."""

    def empty_array(self) -> ma.masked_array:
        """Create an empty array with the dimensions of the variable.

        Returns:
            Empty array.
        """
        return ma.masked_all([len(dim) for dim in self.dimensions])

    def to_nc(
        self,
        values: Optional[npt.ArrayLike] = None,
        file: Optional[Path] = None,
    ) -> None:
        """Write variable to a netCDF file.

        The netCDF file is openend and closed. The file is either specified in the function call or
        taken from the default file. If the variable is not mandatory and the values are all masked,
        the variable is not written. Otherwise, if the variable was not yet added, it is added with
        its dimensions, otherwise its values are overwritten. The values are either supplied in the
        function call or taken from the values attribute.

        Args:
            values: Values to write. Defaults to None.
            file: File to write to. Defaults to None.

        Raises:
            ValueError: Output file or values not defined.
            FileNotFoundError: Could not open file.
            NotImplementedError: Number of dimensions not implemented.
        """
        if values is not None:
            to_write = values
        elif self.values is not None:
            to_write = self.values
        else:
            raise ValueError(f"Values of variable {self.name} not defined.")

        if not self.mandatory:
            if isinstance(to_write, ma.MaskedArray) and to_write.mask.all():
                logger.debug(f"All values of variable {self.name} are masked. Skipping write.")
                return

        if file is not None:
            to_file = file
        elif self.file is not None:
            to_file = self.file
        else:
            raise ValueError(f"Output file for variable {self.name} not defined.")

        try:
            nc_data = Dataset(to_file, "a", format="NETCDF4")
        except FileNotFoundError:
            logger.critical(f"Could not open file {to_file}.")
            raise

        logger.debug(f"Writing array {self.name} to file {to_file}.")
        if self.name not in nc_data.variables:
            for nc_dim in self.dimensions:
                nc_dim.to_dataset(nc_data)

            nc_var = nc_data.createVariable(
                self.name,
                self.data_type,
                (o.name for o in self.dimensions),
                fill_value=self.fill_value,
            )

            for attr in self._metadata:
                attr_value = getattr(self, attr)
                if attr_value is not None:
                    nc_var.setncattr(attr, attr_value)

        else:
            nc_var = nc_data.variables[self.name]

        # When writing, the data type of to_write will be automatically adjusted to self.data_type
        # as defined above. For masked arrays, this includes also masked data and the fill value. If
        # these for the static driver irrelevant values are outside of the target data type, a
        # warning will be raised, e.g. "RuntimeWarning: invalid value encountered in cast" or
        # "RuntimeWarning: overflow encountered in cast". In order to avoid this, data that is
        # masked and the fill value are set to the fill value of the variable. TODO: Check if this
        # is still necessary.
        if isinstance(to_write, ma.MaskedArray):
            to_write = ma.MaskedArray(
                np.where(to_write.mask, self.fill_value, to_write.data),
                mask=to_write.mask,
                fill_value=self.fill_value,
            )

        if len(self.dimensions) == 1:
            nc_var[:] = to_write
        elif len(self.dimensions) == 2:
            nc_var[:, :] = to_write
        elif len(self.dimensions) == 3:
            nc_var[:, :, :] = to_write
        elif len(self.dimensions) == 4:
            nc_var[:, :, :, :] = to_write
        elif len(self.dimensions) == 5:
            nc_var[:, :, :, :, :] = to_write
        else:
            raise NotImplementedError

        nc_data.close()

    def from_nc(
        self, file: Optional[Path] = None, allow_nonexistent: bool = False
    ) -> ma.MaskedArray:
        """Get values from a netCDF file.

        The file is either specified in the function call or taken from the default file.

        Args:
            file: File to read from. Defaults to None.
            allow_nonexistent: If True, return empty array when variable does not exist in the file.
              Defaults to False.

        Raises:
            ValueError: Input file not defined.
            FileNotFoundError: Could not open file.

        Returns:
            Values from the netCDF file. If allow_nonexistent is True, an empty array is returned
            when the variable does not exist.
        """
        if file is not None:
            from_file = file
        elif self.file is not None:
            from_file = self.file
        else:
            raise ValueError(f"Input file for variable {self.name} not defined.")

        try:
            nc_data = Dataset(from_file, "r", format="NETCDF4")
        except FileNotFoundError:
            logger.critical(f"Could not open file {from_file}.")
            raise

        try:
            tmp_array = nc_data.variables[self.name][:]
        except KeyError:
            if allow_nonexistent:
                tmp_array = self.empty_array()
            else:
                logger.critical(f"Variable {self.name} not found in file {from_file}.")
                raise

        nc_data.close()

        return tmp_array


class NCDFCoordinateReferenceSystem(BaseModel, validate_assignment=True):
    """A coordinate reference system that can be written to a netCDF file."""

    long_name: str
    """Long name."""
    grid_mapping_name: str
    """Grid mapping name."""
    semi_major_axis: float
    """Semi-major axis."""
    inverse_flattening: float
    """Inverse flattening."""
    longitude_of_prime_meridian: float
    """Longitude of prime meridian."""
    longitude_of_central_meridian: float
    """Longitude of central meridian."""
    scale_factor_at_central_meridian: float
    """Scale factor at central meridian."""
    latitude_of_projection_origin: float
    """Latitude of projection origin."""
    false_easting: float
    """False easting."""
    false_northing: float
    """False northing."""
    spatial_ref: str
    """Spatial reference."""
    units: str
    """Units."""
    epsg_code: str
    """EPSG code."""

    file: Optional[Path] = None
    """Default file for input/output."""

    def to_nc(self, file: Optional[Path] = None) -> None:
        """Write CRS to a netCDF file.

        The file is openend and closed. The file is either specified in the function call or taken
        from the default `file`.

        Args:
            file: File to write to. Defaults to None.

        Raises:
            ValueError: Output file not defined.
            FileNotFoundError: Could not open file.
        """
        if file is not None:
            to_file = file
        elif self.file is not None:
            to_file = self.file
        else:
            raise ValueError("Output file for CRS not defined")

        try:
            nc_data = Dataset(to_file, "a", format="NETCDF4")
        except FileNotFoundError:
            logger.critical(f"Could not open file {to_file}.")
            raise

        logger.debug(f"Writing crs to file {to_file}.")

        nc_var = nc_data.createVariable("crs", "i")

        nc_var.long_name = self.long_name
        nc_var.grid_mapping_name = self.grid_mapping_name
        nc_var.semi_major_axis = self.semi_major_axis
        nc_var.inverse_flattening = self.inverse_flattening
        nc_var.longitude_of_prime_meridian = self.longitude_of_prime_meridian
        nc_var.longitude_of_central_meridian = self.longitude_of_central_meridian
        nc_var.scale_factor_at_central_meridian = self.scale_factor_at_central_meridian
        nc_var.latitude_of_projection_origin = self.latitude_of_projection_origin
        nc_var.false_easting = self.false_easting
        nc_var.false_northing = self.false_northing
        nc_var.spatial_ref = self.spatial_ref
        nc_var.units = self.units
        nc_var.epsg_code = self.epsg_code

        nc_data.close()


def remove_existing_file(file: Path) -> None:
    """Remove a file if it exists.

    Args:
        file: File to remove.
    """
    try:
        os.remove(file)
        logger.debug(f"Removing {file}.")
    except FileNotFoundError:
        pass
