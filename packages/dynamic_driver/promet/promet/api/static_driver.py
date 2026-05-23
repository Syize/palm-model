import os
import numpy
import netCDF4

from promet.api.logging import PrometException


class StaticDriver:
    """used for loading domain info from static driver file."""

    def __init__(
            self,
            filepath: str,
    ):
        self.filepath = os.path.normpath(os.path.expandvars(os.path.expanduser(filepath)))

    def __enter__(self):
        self.ncfile = netCDF4.Dataset(self.filepath, 'r')

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.ncfile.close()

    def load(self) -> bool:
        try:
            with netCDF4.Dataset(self.filepath, 'r') as nc_file:
                self.dimension_names = list(nc_file.dimensions.keys())
                self.variable_names = list(nc_file.variables.keys())
            return False
        except FileNotFoundError as e:
            raise PrometException(
                f'Missing "{self.filepath}"',
                id='INwPBw'
            )

    def get_global_attribute(
            self,
            attribute: str
    ) -> str:
        with netCDF4.Dataset(self.filepath, 'r') as nc_file:
            return nc_file.getncattr(attribute)

    def get_dimension(
            self,
            dimension: str
    ):
        with netCDF4.Dataset(self.filepath, 'r') as ncfile:
            result = dict(
                size=ncfile.dimensions[dimension].size,
                units=ncfile.variables[dimension].units,
                values=ncfile.variables[dimension][:],
            )
        return result

    def get_terrain_array(self) -> numpy.ndarray:
        with netCDF4.Dataset(self.filepath, 'r') as ncfile:
            try:
                array = ncfile.variables['zt'][:]
            except KeyError as e:
                raise PrometException(
                    f'Missing "zt" in {self.filepath}.',
                    id='1eHciA',
                )
        return array

    def get_metadata(self) -> dict:
        return dict(
            dimension_names=self.dimension_names,
            variable_names=self.variable_names,
        )
