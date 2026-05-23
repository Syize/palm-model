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
# Copyright 1997-2024  Leibniz Universitaet Hannover
# Copyright 2022-2024  Technische Universitaet Berlin

"""Variables and methods for domains."""

import inspect
import logging
from enum import Enum
from math import floor
from pathlib import Path
from typing import Callable, Optional, Tuple, Type, TypedDict, cast

import numpy as np
import numpy.typing as npt
import rasterio.warp as riowp
from netCDF4 import Dataset, Variable
from numpy import ma

from palm_csd import StatusLogger
from palm_csd.csd_config import (
    NBUILDING_SURFACE_LAYER,
    CSDConfig,
    CSDConfigAttributes,
    CSDConfigDomain,
    CSDConfigInput,
    IndexBuildingGeneralPars,
    IndexBuildingIndoorPars,
    IndexBuildingSurfaceLevel,
    IndexBuildingSurfaceType,
    ScalingMethods,
    defaults,
)
from palm_csd.geo_converter import GeoConverter
from palm_csd.lcz import LCZTypes
from palm_csd.netcdf_data import (
    NCDFCoordinateReferenceSystem,
    NCDFDimension,
    NCDFVariable,
    remove_existing_file,
)
from palm_csd.tools import DefaultMinMax

# Module logger. In __init__.py, it is ensured that the logger is a StatusLogger. For type checking,
# do explicit cast.
logger = cast(StatusLogger, logging.getLogger(__name__))


class CSDDomain:
    """A domain that stores all its configurations, output dimensions and variables."""

    name: str
    """Name of the domain."""

    config: CSDConfigDomain
    """Domain configuration."""
    input_config: CSDConfigInput
    """Input configuration."""
    attributes: CSDConfigAttributes
    """Attributes configuration"""

    # TODO: Python 3.11: Use Self to indicate same type as class
    parent: Optional["CSDDomain"]
    """Parent domain."""

    geo_converter: Optional[GeoConverter]
    """GeoConverter to handle geographic data and transformation."""

    file_output: Path
    """Output path."""

    rotation_angle: float
    """Rotation angle."""

    replace_invalid_input_values: bool
    """Replace invalid input values in input files by the respective default value."""

    downscaling_method: ScalingMethods
    """Methods for downscaling."""
    upscaling_method: ScalingMethods
    """Methods for upscaling."""

    origin_x: Optional[float]
    """x-coordinate of the left border of the lower-left grid point of the PALM domain in the
    custom CRS."""
    origin_y: Optional[float]
    """y-coordinate of the lower border of the lower-left grid point of the PALM domain in the
    custom CRS."""
    origin_lon: Optional[float]
    """Longitude of the left border of the lower-left grid point of the PALM domain."""
    origin_lat: Optional[float]
    """Latitude of the lower border of the lower-left grid point of the PALM domain."""
    origin_z: Optional[float]
    """Reference height in m above sea level after DHHN2016."""

    x0: Optional[int]
    """Lowest x-index used for reading input data."""
    y0: Optional[int]
    """Lowest y-index used for reading input data."""
    x1: Optional[int]
    """Highest x-index used for reading input data."""
    y1: Optional[int]
    """Highest y-index used for reading input data."""

    x: NCDFDimension
    """x dimension."""
    y: NCDFDimension
    """y dimension."""
    z: NCDFDimension
    """z dimension used for buildings."""
    zlad: NCDFDimension
    """z dimension used for resolved vegetation."""

    nsurface_fraction: NCDFDimension
    """Surface fraction dimension."""

    building_general_par: NCDFDimension
    """Parameter dimension of building general."""
    building_indoor_par: NCDFDimension
    """Parameter dimension of building indoor."""
    building_surface_layer: NCDFDimension
    """Dimension of building surface layer."""
    building_surface_level: NCDFDimension
    """Dimension of building surface level."""
    building_surface_type: NCDFDimension
    """Dimension of building surface type."""

    nvegetation_pars: NCDFDimension
    """Parameter dimension of vegetation_pars"""
    nwater_pars: NCDFDimension
    """Parameter dimension of water_pars"""

    lat: NCDFVariable
    """Latitude."""
    lon: NCDFVariable
    """Longitude."""

    x_global: NCDFVariable
    """Global x coordinates of all dimensions."""
    y_global: NCDFVariable
    """Global x coordinates of all dimensions."""

    E_UTM: NCDFVariable
    """East UTM coordinates."""
    N_UTM: NCDFVariable
    """North UTM coordinates."""

    zt: NCDFVariable
    """Terrain height relative to origin_z."""

    buildings_2d: NCDFVariable
    """Building height."""
    building_id: NCDFVariable
    """Building ID."""
    building_type: NCDFVariable
    """Building type."""
    buildings_3d: NCDFVariable
    """3D building representation with 0 and 1."""

    surface_fraction: NCDFVariable
    """Surface fraction."""

    vegetation_type: NCDFVariable
    """Vegetation type."""
    pavement_type: NCDFVariable
    """Pavement type."""
    water_type: NCDFVariable
    """Water type."""
    soil_type: NCDFVariable
    """Soil type."""
    street_type: NCDFVariable
    """Street type."""
    street_crossing: NCDFVariable
    """Street crossing."""

    building_albedo_type: NCDFVariable
    """Building albedo type."""
    building_emissivity: NCDFVariable
    """Building emissivity."""
    building_fraction: NCDFVariable
    """Building fraction."""
    building_general_pars: NCDFVariable
    """Building general parameters."""
    building_heat_capacity: NCDFVariable
    """Building heat capacity."""
    building_heat_conductivity: NCDFVariable
    """Building heat conductivity."""
    building_indoor_pars: NCDFVariable
    """Building indoor parameters."""
    building_lai: NCDFVariable
    """Building lai."""
    building_roughness_length: NCDFVariable
    """Building roughness length."""
    building_roughness_length_qh: NCDFVariable
    """Building roughness length for moisture and heat."""
    building_thickness: NCDFVariable
    """Building thickness."""
    building_transmissivity: NCDFVariable
    """Building window transmissivity."""

    vegetation_pars: NCDFVariable
    """Vegetation parameters."""
    water_pars: NCDFVariable
    """Water parameters."""

    nuc: NCDFDimension
    """Urban class."""
    streetdir: NCDFDimension
    """Street direction."""
    z_uhl: NCDFDimension
    """Height urban half level."""

    lad: NCDFVariable
    """Leaf area density."""
    bad: NCDFVariable
    """Basal area density."""
    tree_id: NCDFVariable
    """Tree ID."""
    tree_type: NCDFVariable
    """Tree type."""

    fr_urb: NCDFVariable
    """Fraction of urban area."""
    fr_urbcl: NCDFVariable
    """Fraction of urban classes."""
    fr_streetdir: NCDFVariable
    """Fraction of street directions."""
    street_width: NCDFVariable
    """Street width."""
    building_width: NCDFVariable
    """Building width."""
    building_height: NCDFVariable
    """Building height."""

    def __init__(
        self,
        name: str,
        config: CSDConfig,
        parent: Optional["CSDDomain"] = None,
        gis_debug_output: bool = False,
    ) -> None:
        """Initialize domain.

        Copy configurations, initialize geo converter if needed, set output file name and initialize
        dimensions and variables.

        Args:
            name: Name of the domain.
            config: palm_csd configuration.
            parent: Parent domain configuration. Defaults to None.
            gis_debug_output: Write out reprojected data for debugging. Defaults to False.

        Raises:
            ValueError: When using geo converter, geo converter of parent is None.
        """
        self.name = name

        # configurations
        self.config = config.domain_dict[name]
        self.input_config = config.input_of_domain(self.name)

        self.attributes = config.attributes

        self.parent = parent

        # converter for geo data
        if config.settings.epsg is None:
            self.geo_converter = None
        else:
            # Find root parent
            if self.parent is not None:
                parent_geoconverter = self.parent.geo_converter
                if parent_geoconverter is None:
                    raise ValueError("Parent domain has no geo converter")
                parent_tmp = self.parent
                while parent_tmp.parent is not None:
                    parent_tmp = parent_tmp.parent
                root_parent_geoconverter = parent_tmp.geo_converter
                if root_parent_geoconverter is None:
                    raise ValueError("Root parent domain has no geo converter")
            else:
                parent_geoconverter = None
                root_parent_geoconverter = None
            logger.info(f"Setting up of coordinate calculation for domain {self.name}.")
            self.geo_converter = GeoConverter(
                self.config,
                config.settings,
                config.output,
                parent_geoconverter,
                root_parent_geoconverter,
                self.name,
                debug_output=gis_debug_output,
            )

        self.replace_invalid_input_values = config.settings.replace_invalid_input_values

        # set output file name: file_out + domain name
        # TODO: use with_stem with Python 3.9
        self.file_output = config.output.file_out.with_name(
            f"{config.output.file_out.name}_{self.name}"
        )

        self.rotation_angle = config.settings.rotation_angle

        if (
            self.config.input_lower_left_x is not None
            and self.config.input_lower_left_y is not None
        ):
            self.x0 = int(floor(self.config.input_lower_left_x / self.config.pixel_size))
            self.y0 = int(floor(self.config.input_lower_left_y / self.config.pixel_size))
            self.x1 = self.x0 + self.config.nx
            self.y1 = self.y0 + self.config.ny
        else:
            if self.input_config.any_netcdf():
                logger.critical_raise(
                    "input_lower_left_x and input_lower_left_y must be set "
                    + f"in the domain section of domain {self.name} for netCDF input."
                )
            self.x0 = None
            self.y0 = None
            self.x1 = None
            self.y1 = None

        self.parent = parent

        self.downscaling_method = config.settings.downscaling_method
        self.upscaling_method = config.settings.upscaling_method

        self.check_consistency()

        self._initialize_dimensions()
        self._initialize_variables()

    def _initialize_dimensions(self) -> None:
        """Initialize dimensions."""
        self.x = NCDFDimension(
            name="x",
            datatype="f4",
            standard_name="projection_x_coordinate",
            long_name="x",
            units="m",
        )
        self.y = NCDFDimension(
            name="y",
            datatype="f4",
            standard_name="projection_y_coordinate",
            long_name="y",
            units="m",
        )
        self.z = NCDFDimension(name="z", datatype="f4", long_name="z", units="m")

        self.nsurface_fraction = NCDFDimension(name="nsurface_fraction", datatype="i")

        self.building_general_par = NCDFDimension(
            name="building_general_par",
            datatype=str,
            values=_enum_to_str_array(IndexBuildingGeneralPars),
        )
        self.building_indoor_par = NCDFDimension(
            name="building_indoor_par",
            datatype=str,
            values=_enum_to_str_array(IndexBuildingIndoorPars),
        )
        self.building_surface_layer = NCDFDimension(
            name="building_surface_layer",
            datatype="i",
            values=np.arange(1, NBUILDING_SURFACE_LAYER + 1),
        )
        self.building_surface_level = NCDFDimension(
            name="building_surface_level",
            datatype=str,
            values=_enum_to_str_array(IndexBuildingSurfaceLevel),
        )
        self.building_surface_type = NCDFDimension(
            name="building_surface_type",
            datatype=str,
            values=_enum_to_str_array(IndexBuildingSurfaceType),
        )

        self.nvegetation_pars = NCDFDimension(
            name="nvegetation_pars", datatype="i", values=ma.arange(0, 12)
        )
        self.nwater_pars = NCDFDimension(name="nwater_pars", datatype="i", values=np.arange(0, 7))

        self.zlad = NCDFDimension(name="zlad", datatype="f4")

        self.nuc = NCDFDimension(name="nuc", datatype="i")
        self.streetdir = NCDFDimension(name="streetdir", datatype="i")
        self.z_uhl = NCDFDimension(name="z_uhl", datatype="f4")

    def _initialize_variables(self):
        """Initialize variables."""
        dimensions_yx = (self.y, self.x)
        dimensions_zladyx = (self.zlad, self.y, self.x)
        dimensions_levelyx = (
            self.building_surface_level,
            self.y,
            self.x,
        )
        dimensions_typelayeryx = (
            self.building_surface_type,
            self.building_surface_layer,
            self.y,
            self.x,
        )

        # variables
        self.lat = self._variable_float(
            name="lat",
            dimensions=dimensions_yx,
            long_name="latitude",
            standard_name="latitude",
            units="degrees_north",
            add_spatial_metadata=False,
        )
        self.lon = self._variable_float(
            name="lon",
            dimensions=dimensions_yx,
            long_name="longitude",
            standard_name="longitude",
            units="degrees_east",
            add_spatial_metadata=False,
        )

        self.x_global = self._variable_float(
            name="x_UTM",
            dimensions=(self.x,),
            long_name="easting",
            standard_name="projection_x_coordinate",
            units="m",
            add_spatial_metadata=False,
        )

        self.y_global = self._variable_float(
            name="y_UTM",
            dimensions=(self.y,),
            long_name="northing",
            standard_name="projection_y_coordinate",
            units="m",
            add_spatial_metadata=False,
        )

        self.E_UTM = self._variable_float(
            name="E_UTM",
            dimensions=dimensions_yx,
            long_name="easting",
            standard_name="projection_x_coordinate",
            units="m",
            add_spatial_metadata=False,
        )

        self.N_UTM = self._variable_float(
            name="N_UTM",
            dimensions=dimensions_yx,
            long_name="northing",
            standard_name="projection_y_coordinate",
            units="m",
            add_spatial_metadata=False,
        )

        self.zt = self._variable_float(
            name="zt",
            dimensions=dimensions_yx,
            long_name="orography",
            units="m",
        )

        self.buildings_2d = self._variable_float(
            name="buildings_2d",
            dimensions=dimensions_yx,
            long_name="buildings",
            units="m",
            lod=1,
        )

        self.building_id = self._variable_int(
            name="building_id",
            dimensions=dimensions_yx,
            long_name="building id",
            units="",
        )

        self.building_type = self._variable_byte(
            name="building_type",
            dimensions=dimensions_yx,
            long_name="building type",
            units="",
        )

        self.buildings_3d = self._variable_byte(
            name="buildings_3d",
            dimensions=(self.z, self.y, self.x),
            long_name="buildings 3d",
            units="",
            lod=2,
        )

        self.surface_fraction = self._variable_float(
            name="surface_fraction",
            dimensions=(self.nsurface_fraction, self.y, self.x),
            long_name="surface fraction",
            units="1",
        )

        self.vegetation_type = self._variable_byte(
            name="vegetation_type",
            dimensions=dimensions_yx,
            long_name="vegetation type",
            units="",
        )

        self.pavement_type = self._variable_byte(
            name="pavement_type",
            dimensions=dimensions_yx,
            long_name="pavement type",
            units="",
        )

        self.water_type = self._variable_byte(
            name="water_type",
            dimensions=dimensions_yx,
            long_name="water type",
            units="",
        )

        self.soil_type = self._variable_byte(
            name="soil_type",
            dimensions=dimensions_yx,
            long_name="soil type",
            units="",
        )

        self.street_type = self._variable_byte(
            name="street_type",
            dimensions=dimensions_yx,
            long_name="street type",
            units="",
        )

        self.street_crossing = self._variable_byte(
            name="street_crossing",
            dimensions=dimensions_yx,
            long_name="street crossings",
            units="",
        )

        self.building_albedo_type = self._variable_int(
            name="building_albedo_type",
            dimensions=(self.building_surface_type, self.y, self.x),
            long_name="building surface albedo types",
            units="",
        )

        self.building_emissivity = self._variable_float(
            name="building_emissivity",
            dimensions=(self.building_surface_type, self.y, self.x),
            long_name="building surface emissivity",
            units="1",
        )

        self.building_fraction = self._variable_float(
            name="building_fraction",
            dimensions=(self.building_surface_type, self.y, self.x),
            long_name="building surface fractions",
            units="1",
        )

        self.building_general_pars = self._variable_float(
            name="building_general_pars",
            dimensions=(self.building_general_par, self.y, self.x),
            long_name="general building parameters",
            units="",
        )

        self.building_heat_capacity = self._variable_float(
            name="building_heat_capacity",
            dimensions=dimensions_typelayeryx,
            long_name="building surface layer heat capacities",
            units="J m-3 K-1",
        )

        self.building_heat_conductivity = self._variable_float(
            name="building_heat_conductivity",
            dimensions=dimensions_typelayeryx,
            long_name="building surface layer heat conductivities",
            units="W m-1 K-1",
        )

        self.building_indoor_pars = self._variable_float(
            name="building_indoor_pars",
            dimensions=(self.building_indoor_par, self.y, self.x),
            long_name="building indoor parameters",
            units="",
        )

        self.building_lai = self._variable_float(
            name="building_lai",
            dimensions=dimensions_levelyx,
            long_name="building surface lai",
            units="m2 m-2",
        )

        self.building_roughness_length = self._variable_float(
            name="building_roughness_length",
            dimensions=dimensions_levelyx,
            long_name="building surface roughness lengths",
            units="m",
        )

        self.building_roughness_length_qh = self._variable_float(
            name="building_roughness_length_qh",
            dimensions=dimensions_levelyx,
            long_name="building surface roughness lengths for moisture and heat",
            units="m",
        )

        self.building_thickness = self._variable_float(
            name="building_thickness",
            dimensions=dimensions_typelayeryx,
            long_name="building surface layer thicknesses",
            units="m",
        )

        self.building_transmissivity = self._variable_float(
            name="building_transmissivity",
            dimensions=dimensions_levelyx,
            long_name="building window transmissivities",
            units="1",
        )

        self.vegetation_pars = self._variable_float(
            name="vegetation_pars",
            dimensions=(self.nvegetation_pars, self.y, self.x),
            long_name="vegetation_pars",
            units="",
        )

        self.water_pars = self._variable_float(
            name="water_pars",
            dimensions=(self.nwater_pars, self.y, self.x),
            long_name="water_pars",
            units="",
        )

        self.lad = self._variable_float(
            name="lad",
            dimensions=dimensions_zladyx,
            long_name="leaf area density",
            units="m2 m-3",
        )

        self.bad = self._variable_float(
            name="bad",
            dimensions=dimensions_zladyx,
            long_name="basal area density",
            units="m2 m-3",
        )

        self.tree_id = self._variable_int(
            name="tree_id",
            dimensions=dimensions_zladyx,
            long_name="tree id",
            units="",
        )

        self.tree_type = self._variable_byte(
            name="tree_type",
            dimensions=dimensions_zladyx,
            long_name="tree type",
            units="",
        )

        self.fr_urb = self._variable_float(
            name="fr_urb",
            dimensions=dimensions_yx,
            long_name="fraction of urban area",
            units="1",
        )
        self.fr_urbcl = self._variable_float(
            name="fr_urbcl",
            dimensions=(self.nuc, self.y, self.x),
            long_name="fraction of urban classes",
            units="1",
        )
        self.fr_streetdir = self._variable_float(
            name="fr_streetdir",
            dimensions=(self.nuc, self.streetdir, self.y, self.x),
            long_name="fraction of street directions",
            units="1",
        )
        self.street_width = self._variable_float(
            name="street_width",
            dimensions=(self.nuc, self.streetdir, self.y, self.x),
            long_name="street width",
            units="m",
        )
        self.building_width = self._variable_float(
            name="building_width",
            dimensions=(self.nuc, self.streetdir, self.y, self.x),
            long_name="building width",
            units="m",
        )
        self.building_height = self._variable_float(
            name="building_height",
            dimensions=(self.nuc, self.streetdir, self.z_uhl, self.y, self.x),
            long_name="building height",
            units="m",
        )

    def check_consistency(self) -> None:
        """Check consistency of domain configuration.

        Raises:
            ValueError: Neither all coordinate inputs provided nor target CRS defined.
        """
        # Geographical coordinate input

        if (
            self.input_config.file_x_UTM is None
            or self.input_config.file_y_UTM is None
            or self.input_config.file_lon is None
            or self.input_config.file_lat is None
        ):
            if self.geo_converter is None:
                logger.critical_raise(
                    f"Not all coordinate inputs provided for domain {self.name}, "
                    + "but target coordinate reference system not defined by EPSG code "
                    + " to calculate coordinates.",
                )

    class DefaultMetadataVariable(TypedDict, total=False):
        """Helper dictionary for default metadata in variables."""

        datatype: str
        """Datatype."""
        fillvalue: float
        """Fill value."""
        file: Path
        """Input/Output file."""
        res_orig: float
        """Original resolution."""
        coordinates: str
        """Coordinates."""
        grid_mapping: str
        """Grid mapping."""

    def _variable_float(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        add_spatial_metadata: bool = True,
    ) -> NCDFVariable:
        """Helper function that returns a float variable with some predefined attributes.

        Args:
            name: Name.
            dimensions: Dimensions.
            long_name: Long name.
            units: Unit.
            standard_name: Standard name. Defaults to None.
            lod: Level of detail. Defaults to None.
            add_spatial_metadata: Add res_orig, coordinates, grid_mapping. Defaults to True.

        Returns:
            float variable with given and predefined attributes.
        """
        default_values: CSDDomain.DefaultMetadataVariable = {
            "datatype": "f4",
            "fillvalue": -9999.0,
            "file": self.file_output,
        }
        if add_spatial_metadata:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def _variable_int(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        add_spatial_metadata: bool = True,
    ) -> NCDFVariable:
        """Helper function that returns an int variables with some predefined attributes.

        Args:
            name: Name.
            dimensions: Dimensions.
            long_name: Long name.
            units: Unit.
            standard_name: Standard name. Defaults to None.
            lod: Level of detail. Defaults to None.
            add_spatial_metadata: Add res_orig, coordinates, grid_mapping. Defaults to True.

        Returns:
            int variable with given and predefined attributes.
        """
        default_values: CSDDomain.DefaultMetadataVariable = {
            "datatype": "i",
            "fillvalue": -9999,
            "file": self.file_output,
        }
        if add_spatial_metadata:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def _variable_byte(
        self,
        name: str,
        dimensions: Tuple[NCDFDimension, ...],
        long_name: str,
        units: str,
        standard_name: Optional[str] = None,
        lod: Optional[int] = None,
        add_spatial_metadata: bool = True,
    ) -> NCDFVariable:
        """Helper function that returns a byte variables with some predefined attributes.

        Args:
            name: Name.
            dimensions: Dimensions.
            long_name: Long name.
            units: Unit.
            standard_name: Standard name. Defaults to None.
            lod: Level of detail. Defaults to None.
            add_spatial_metadata: Add res_orig, coordinates, grid_mapping. Defaults to True.

        Returns:
            byte variable with given and predefined attributes.
        """
        default_values: CSDDomain.DefaultMetadataVariable = {
            "datatype": "b",
            "fillvalue": -127,
            "file": self.file_output,
        }
        if add_spatial_metadata:
            default_values.update(
                {
                    "res_orig": self.config.pixel_size,
                    "coordinates": "E_UTM N_UTM lon lat",
                    "grid_mapping": "crs",
                }
            )

        return NCDFVariable(
            name=name,
            dimensions=dimensions,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
            lod=lod,
            **default_values,
        )

    def remove_existing_output(self) -> None:
        """Remove configured output file if it exists."""
        remove_existing_file(self.file_output)

    def write_global_attributes(self) -> None:
        """Write global attributes to the netCDF.

        Attributes are written to self.file_output. None attributes are not added.
        """
        logger.debug(f"Writing global attributes to file {self.file_output}.")

        nc_data = Dataset(self.file_output, "a", format="NETCDF4")

        nc_data.setncattr("Conventions", "CF-1.7")

        all_attributes = vars(self.attributes)
        for attribute in all_attributes:
            if all_attributes[attribute] is not None:
                nc_data.setncattr(attribute, all_attributes[attribute])

        # add additional attributes
        for attribute in [
            "rotation_angle",
            "origin_x",
            "origin_y",
            "origin_lon",
            "origin_lat",
            "origin_z",
        ]:
            if getattr(self, attribute) is not None:
                nc_data.setncattr(attribute, getattr(self, attribute))
            else:
                raise Exception(f"Attribute {attribute} not set.")

        nc_data.close()

    def write_crs_to_file(self) -> None:
        """Write coordinate reference system information in CF convention to the netCDF.

        CRS data is written to self.file_output. Values are taken from geo_converter's dst_crs.
        """
        logger.debug(f"Writing crs to file {self.file_output}.")

        if self.geo_converter is None:
            raise ValueError("geoconverter must not be None.")
        crs_dict = self.geo_converter.dst_crs_to_cf()

        try:
            nc_data = Dataset(self.file_output, "a", format="NETCDF4")
        except FileNotFoundError:
            logger.critical(f"Could not open file {self.file_output}.")
            raise

        nc_var = nc_data.createVariable("crs", "i")

        # Add long_name
        nc_var.setncattr("long_name", "coordinate reference system")

        # Add crs information
        for key, value in crs_dict.items():
            nc_var.setncattr(key, value)

        nc_data.close()

    def read_nc_3d(
        self,
        file: Optional[Path],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
        y0: Optional[int] = None,
        y1: Optional[int] = None,
        z0: Optional[int] = None,
        z1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 3d raster data from a netCDF file.

        The file is openend and closed. If file is None, the values of the returned array are all
        masked. The default boundary coordinates are taken from the containing domain. If complete,
        the full variable is read.

        Args:
            file: Input file.
            varname: Variable name. Defaults to None.
            complete: If true, read complete data. Defaults to False.
            x0: Lowest x index. Defaults to None.
            x1: Highest x index. Defaults to None.
            y0: Lowest y index. Defaults to None.
            y1: Highest y index. Defaults to None.
            z0: Lowest z index. Defaults to None.
            z1: Highest z index. Defaults to None.

        Raises:
            NotImplementedError: z0 or z1 is None and complete is False.
            ValueError: file is None and complete is True.

        Returns:
            3d variable data.
        """
        if x0 is None:
            x0 = self.x0
        if x1 is None:
            x1 = self.x1
        if y0 is None:
            y0 = self.y0
        if y1 is None:
            y1 = self.y1
        if z0 is None:
            if not complete:
                raise NotImplementedError
        if z1 is None:
            if not complete:
                raise NotImplementedError

        if file is not None:
            try:
                nc_data = Dataset(file, "r", format="NETCDF4")
            except FileNotFoundError:
                logger.critical(f"Could not open file {file}.")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 3)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:, :, :]
            else:
                tmp_array = variable[z0 : (z1 + 1), y0 : (y1 + 1), x0 : (x1 + 1)]  # type: ignore
            nc_data.close()
        else:
            if complete:
                raise ValueError("file needs to be given when complete==True.")
            tmp_array = ma.masked_all((z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1))  # type: ignore

        return tmp_array

    def read_nc_2d(
        self,
        file: Optional[Path],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
        y0: Optional[int] = None,
        y1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 2d raster data from a netCDF file.

        The file is openend and closed. If file is None, the values of the returned array are all
        masked. The default boundary coordinates are taken from the containing domain. If complete,
        the full variable is read.

        Args:
            file: Input file.
            varname: Variable name. Defaults to None.
            complete: If true, read complete data. Defaults to False.
            x0: Lowest x index. Defaults to None.
            x1: Highest x index. Defaults to None.
            y0: Lowest y index. Defaults to None.
            y1: Highest y index. Defaults to None.

        Raises:
            ValueError: file is None and complete is True.

        Returns:
            2d variable data.
        """
        if x0 is None:
            x0 = self.x0
        if x1 is None:
            x1 = self.x1
        if y0 is None:
            y0 = self.y0
        if y1 is None:
            y1 = self.y1

        if file is not None:
            try:
                nc_data = Dataset(file, "r", format="NETCDF4")
            except FileNotFoundError:
                logger.critical(f"Could not open file {file}.")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 2)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:, :]
            else:
                if x0 is None or x1 is None or y0 is None or y1 is None:
                    raise ValueError("x0, x1, y0 and y1 must not be None.")
                tmp_array = variable[y0 : (y1 + 1), x0 : (x1 + 1)]
            nc_data.close()
        else:
            if complete:
                raise ValueError("file needs to be given when complete==True.")

            if y0 is None or y1 is None:
                shape_y = self.config.ny + 1
            else:
                shape_y = y1 - y0 + 1
            if x0 is None or x1 is None:
                shape_x = self.config.nx + 1
            else:
                shape_x = x1 - x0 + 1

            tmp_array = ma.masked_all((shape_y, shape_x))

        return tmp_array

    def read_nc_1d(
        self,
        file: Optional[Path],
        varname: Optional[str] = None,
        complete: bool = False,
        x0: Optional[int] = None,
        x1: Optional[int] = None,
    ) -> ma.MaskedArray:
        """Read a 1d raster data from a netCDF file.

        The file is openend and closed. If file is None, the values of the returned array are all
        masked. The default boundary coordinates are taken from the containing domain. If complete,
        the full variable is read.

        Args:
            file: Input file.
            varname: Variable name. Defaults to None.
            complete: If true, read complete data. Defaults to False.
            x0: Lowest x index. Defaults to None.
            x1: Highest x index. Defaults to None.

        Raises:
            ValueError: file is None and complete is True.

        Returns:
            1d variable data.
        """
        if x0 is None:
            x0 = self.x0
        if x1 is None:
            x1 = self.x1

        if file is not None:
            try:
                nc_data = Dataset(file, "r", format="NETCDF4")
            except FileNotFoundError:
                logger.critical(f"Could not open file {file}.")
                raise

            if varname is None:
                variable = _find_variable_name(nc_data, 2)
            else:
                variable = nc_data.variables[varname]

            if complete:
                tmp_array = variable[:]
            else:
                if x0 is None or x1 is None:
                    raise ValueError("x0 and x1 must not be None.")
                tmp_array = variable[x0 : (x1 + 1)]
            nc_data.close()
        else:
            if complete:
                raise ValueError("file needs to be given when complete==True")
            if x0 is None or x1 is None:
                raise ValueError("x0 and x1 must not be None.")
            tmp_array = ma.masked_all(x1 - x0 + 1)

        return tmp_array

    def read_nc_crs(
        self, file: Optional[Path] = None, varname: Optional[str] = None
    ) -> NCDFCoordinateReferenceSystem:
        """Read coordinate reference system from a netCDF file.

        The file is openend and closed. If file is None, self.input_config.file_x_UTM is used.

        Args:
            file: Input file. Defaults to None.
            varname: Variable name. Defaults to None.

        Raises:
            ValueError: Both file and input_config.file_x_UTM are None.

        Returns:
            Coordinate reference system from file.
        """
        if file is not None:
            from_file = file
        elif self.input_config.file_x_UTM is not None:
            from_file = self.input_config.file_x_UTM
        else:
            raise ValueError("file or input_config.file_x_UTM needs to be not None")

        try:
            nc_data = Dataset(from_file, "r", format="NETCDF4")
        except FileNotFoundError:
            logger.critical(f"Could not open file {from_file}.")
            raise

        if varname is None:
            variable = _find_variable_name(nc_data, 2)
        else:
            variable = nc_data.variables[varname]
        crs_from_file = nc_data.variables[variable.grid_mapping]

        # Get EPSG code from crs
        try:
            epsg_code = crs_from_file.epsg_code
        except AttributeError:
            epsg_code = "unknown"
            if crs_from_file.spatial_ref.find("ETRS89", 0, 100) and crs_from_file.spatial_ref.find(
                "UTM", 0, 100
            ):
                if crs_from_file.spatial_ref.find("28N", 0, 100) != -1:
                    epsg_code = "EPSG:25828"
                elif crs_from_file.spatial_ref.find("29N", 0, 100) != -1:
                    epsg_code = "EPSG:25829"
                elif crs_from_file.spatial_ref.find("30N", 0, 100) != -1:
                    epsg_code = "EPSG:25830"
                elif crs_from_file.spatial_ref.find("31N", 0, 100) != -1:
                    epsg_code = "EPSG:25831"
                elif crs_from_file.spatial_ref.find("32N", 0, 100) != -1:
                    epsg_code = "EPSG:25832"
                elif crs_from_file.spatial_ref.find("33N", 0, 100) != -1:
                    epsg_code = "EPSG:25833"
                elif crs_from_file.spatial_ref.find("34N", 0, 100) != -1:
                    epsg_code = "EPSG:25834"
                elif crs_from_file.spatial_ref.find("35N", 0, 100) != -1:
                    epsg_code = "EPSG:25835"
                elif crs_from_file.spatial_ref.find("36N", 0, 100) != -1:
                    epsg_code = "EPSG:25836"
                elif crs_from_file.spatial_ref.find("37N", 0, 100) != -1:
                    epsg_code = "EPSG:25837"

        crs_var = NCDFCoordinateReferenceSystem(
            long_name="coordinate reference system",
            grid_mapping_name=crs_from_file.grid_mapping_name,
            semi_major_axis=crs_from_file.semi_major_axis,
            inverse_flattening=crs_from_file.inverse_flattening,
            longitude_of_prime_meridian=crs_from_file.longitude_of_prime_meridian,
            longitude_of_central_meridian=crs_from_file.longitude_of_central_meridian,
            scale_factor_at_central_meridian=crs_from_file.scale_factor_at_central_meridian,
            latitude_of_projection_origin=crs_from_file.latitude_of_projection_origin,
            false_easting=crs_from_file.false_easting,
            false_northing=crs_from_file.false_northing,
            spatial_ref=crs_from_file.spatial_ref,
            units="m",
            epsg_code=epsg_code,
            file=self.file_output,
        )

        nc_data.close()

        return crs_var

    def read_transform_2d(
        self,
        name: str,
        file: Optional[Path] = None,
        default_min_max: Optional[DefaultMinMax] = None,
        all_or_none_missing: bool = False,
        initialize_default: bool = False,
        resampling_downscaling: riowp.Resampling = riowp.Resampling.nearest,
        resampling_upscaling: riowp.Resampling = riowp.Resampling.nearest,
        compatibility_resampling_downscaling: Optional[riowp.Resampling] = None,
        compatibility_resampling_upscaling: Optional[riowp.Resampling] = None,
        warning_point_data: bool = False,
        mod_func: Optional[Callable[..., ma.MaskedArray]] = None,
        **kwargs,
    ) -> ma.MaskedArray:
        """Read 2d raster data, and apply geographic transformation and data modification.

        If file is a .nc file, assume it is a 2d netCDF file and its values are returned. Otherwise,
        assume it is a general GIS raster file with defined but arbitrary projection. It is cut to
        the target grid if the grids align, otherwise it is reprojected to the output projection
        with the supplied resampling method. If mod_func is supplied, this function is applied to
        the raster with the **kwargs as input. If not, the raster's first band is used.

        Minimum and maximum values are check with the default_min_max or the values from defaults.
        If replace_invalid_input_values is True, out of range values are replaced by the default
        value. If False, an error is raised. If initialize_default is True, masked values are
        replaced by the default value.

        Args:
            name: Variable name.
            file: Input file. Defaults to None.
            default_min_max: Default, minimum and maximum values. Defaults to None.
            all_or_none_missing: If true, raise an error if both missing and defined values are
              found. Defaults to False.
            initialize_default: If true, replace masked values by the default. Defaults to False.
            resampling_downscaling: Resampling downscaling method. Defaults to
              riowp.Resampling.nearest.
            resampling_upscaling: Resampling upscaling method. Defaults to riowp.Resampling.nearest.
            compatibility_resampling_downscaling: Masked values of this resampling method should be
              applied to the output when downscaling. Defaults to None.
            compatibility_resampling_upscaling: Masked values of this resampling method should be
              applied to the output when upscaling. Defaults to None.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            mod_func: Function to modify the data. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func.

        Raises:
            ValueError: Non netCDF file and geo_converter not set.
            ValueError: Read data not within range and replace_invalid_input_values is False.
            ValueError: No default value defined when replacement is necessary.

        Returns:
            Read and modified 2d raster data.
        """
        if file is None:
            file = getattr(self.input_config, f"file_{name}", None)

        if default_min_max is None:
            default_min_max = defaults[name]

        # If input_file is a netcdf file, read it directly; this handles also None.
        if file is None or file.suffix == ".nc":
            raster_values = self.read_nc_2d(file)

        # Otherwise, assume other GIS raster formats.
        else:
            if self.geo_converter is None:
                raise ValueError("geo_converter not set.")
            raster_values = self.geo_converter.read_to_dst(
                file,
                resampling_downscaling=resampling_downscaling,
                resampling_upscaling=resampling_upscaling,
                compatibility_resampling_downscaling=compatibility_resampling_downscaling,
                compatibility_resampling_upscaling=compatibility_resampling_upscaling,
                warning_point_data=warning_point_data,
                name=name,
            )

            # Apply modification function if supplied.
            if mod_func is not None:
                raster_values = mod_func(raster_values, **kwargs)
            else:
                raster_values = raster_values[0, :, :]

            # Flip raster vertically to convert from GIS to netcdf convention.
            raster_values = ma.MaskedArray(np.flipud(raster_values))

        # Check values.
        below_minimum = raster_values < default_min_max.minimum
        above_maximum = raster_values > default_min_max.maximum
        replacement = default_min_max.default if default_min_max.default is not None else ma.masked
        if default_min_max.minimum is not None and ma.any(below_minimum):
            if not self.replace_invalid_input_values:
                logger.critical(
                    f"In {name}, {ma.sum(below_minimum)} values are smaller than "
                    + f"minimum value {default_min_max.minimum}.\n"
                    + "Enable replace_invalid_input_values for automatic replacement.\n"
                    + "Alternatively, adjust the minimum defined in "
                    + "palm_csd/data/value_defaults.csv."
                )
                raise ValueError("Invalid input values found in {name}.")
            logger.warning(
                f"In {name}, replacing {ma.sum(below_minimum)} values "
                + f"smaller than minimum value {default_min_max.minimum} "
                + f"by default {replacement}."
            )
            logger.debug_indent(
                "If this is unintended, disable replace_invalid_input_values or "
                + "adjust the minimum in palm_csd/data/value_defaults.csv."
            )
            raster_values = ma.where(below_minimum, replacement, raster_values)
        if default_min_max.maximum is not None and ma.any(above_maximum):
            if not self.replace_invalid_input_values:
                logger.critical(
                    f"In {name}, {ma.sum(above_maximum)} values are larger than "
                    + f"maximum value {default_min_max.maximum}.\n"
                    + "Enable replace_invalid_input_values for automatic replacement.\n"
                    + "Alternatively, adjust the maximum defined in "
                    + "palm_csd/data/value_defaults.csv."
                )
                raise ValueError("Invalid input values found in {name}.")
            logger.warning(
                f"In {name}, replacing {ma.sum(above_maximum)} values "
                + f"larger than maximum value {default_min_max.maximum} "
                + f"by default {replacement}."
            )
            logger.debug_indent(
                "If this is unintended, disable replace_invalid_input_values or "
                + "adjust the maximum in palm_csd/data/value_defaults.csv."
            )
            raster_values = ma.where(above_maximum, replacement, raster_values)

        if raster_values.mask.any():
            if all_or_none_missing and not raster_values.mask.all():
                logger.critical_raise(
                    "Found both missing and defined values in {name}.\n"
                    + f"Missing values in {name} are only allowed if all values are missing.",
                )
            if initialize_default:
                if default_min_max.default is None:
                    logger.critical_raise(
                        f"No default value defined for {name} to replace missing values.\n"
                        + "Set a default in palm_csd/data/value_defaults.csv.",
                    )
                logger.debug(
                    f"In {name}, replacing {ma.sum(ma.getmaskarray(raster_values))} missing values "
                    + f"by default {default_min_max.default}."
                )
                raster_values = ma.where(raster_values.mask, default_min_max.default, raster_values)

        return raster_values

    def _read_named_categorical(
        self,
        all_or_none_missing: bool = False,
        initialize_default: bool = False,
        warning_point_data: bool = False,
        mod_func: Optional[Callable[..., ma.MaskedArray]] = None,
        **kwargs,
    ) -> ma.MaskedArray:
        """Read categorical raster data.

        The resampling methods are meant for categorical data (e.g. building type). For both,
        downscaling and upscaling, nearest neighbour resampling is used.

        Args:
            all_or_none_missing: If true, raise an error if both missing and defined values are
              found. Defaults to False.
            initialize_default: If true, replace masked values by the default. Defaults to False.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            mod_func: Function to modify the data. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func.

        Returns:
            Read and modified 2d raster data.
        """
        variable = _calling_read_variable_name()

        return self.read_transform_2d(
            name=variable,
            resampling_downscaling=self.downscaling_method["categorical"],
            resampling_upscaling=self.upscaling_method["categorical"],
            compatibility_resampling_downscaling=riowp.Resampling.nearest,
            compatibility_resampling_upscaling=riowp.Resampling.nearest,
            all_or_none_missing=all_or_none_missing,
            initialize_default=initialize_default,
            warning_point_data=warning_point_data,
            mod_func=mod_func,
            **kwargs,
        )

    def _read_named_continuous(
        self,
        all_or_none_missing: bool = False,
        initialize_default: bool = False,
        warning_point_data: bool = False,
        mod_func: Optional[Callable[..., ma.MaskedArray]] = None,
        **kwargs,
    ) -> ma.MaskedArray:
        """Read continuous raster data.

        The resampling methods are meant for (ususally) continuous data (e.g. terrain height). For
        downscaling, the user defined resampling method is used; for upscaling, the average is used.

        Args:
            all_or_none_missing: If true, raise an error if both missing and defined values are
              found. Defaults to False.
            initialize_default: If true, replace masked values by the default. Defaults to False.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            mod_func: Function to modify the data. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func.

        Returns:
            Read and modified 2d raster data.
        """
        variable = _calling_read_variable_name()

        return self.read_transform_2d(
            name=variable,
            resampling_downscaling=self.downscaling_method["continuous"],
            resampling_upscaling=self.upscaling_method["continuous"],
            compatibility_resampling_downscaling=riowp.Resampling.nearest,
            compatibility_resampling_upscaling=riowp.Resampling.nearest,
            all_or_none_missing=all_or_none_missing,
            initialize_default=initialize_default,
            warning_point_data=warning_point_data,
            mod_func=mod_func,
            **kwargs,
        )

    def _read_named_discrete(
        self,
        all_or_none_missing: bool = False,
        initialize_default: bool = False,
        warning_point_data: bool = True,
        mod_func: Optional[Callable[..., ma.MaskedArray]] = None,
        **kwargs,
    ) -> ma.MaskedArray:
        """Read discrete raster data.

        The resampling methods are meant for (ususally) discrete data (e.g. single tree data). For
        downscaling, the user defined resampling method is used; for upscaling, the average is used.

        Args:
            all_or_none_missing: If true, raise an error if both missing and defined values are
              found. Defaults to False.
            initialize_default: If true, replace masked values by the default. Defaults to False.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            mod_func: Function to modify the data. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func.

        Returns:
            Read and modified 2d raster data.
        """
        variable = _calling_read_variable_name()

        return self.read_transform_2d(
            name=variable,
            resampling_downscaling=self.downscaling_method["discrete"],
            resampling_upscaling=self.upscaling_method["discrete"],
            compatibility_resampling_downscaling=riowp.Resampling.nearest,
            compatibility_resampling_upscaling=riowp.Resampling.nearest,
            all_or_none_missing=all_or_none_missing,
            initialize_default=initialize_default,
            warning_point_data=warning_point_data,
            mod_func=mod_func,
            **kwargs,
        )

    def _read_named_discontinuous(
        self,
        all_or_none_missing: bool = False,
        initialize_default: bool = False,
        warning_point_data: bool = False,
        mod_func: Optional[Callable[..., ma.MaskedArray]] = None,
        **kwargs,
    ) -> ma.MaskedArray:
        """Read discontinuous raster data.

        The resampling methods are meant for (ususally) discontinous data (e.g. building height).
        For downscaling, the user defined resampling method is used; for upscaling, the average is
        used.

        Args:
            all_or_none_missing: If true, raise an error if both missing and defined values are
              found. Defaults to False.
            initialize_default: If true, replace masked values by the default. Defaults to False.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            mod_func: Function to modify the data. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func.

        Returns:
            Read and modified 2d raster data.
        """
        variable = _calling_read_variable_name()

        return self.read_transform_2d(
            name=variable,
            resampling_downscaling=self.downscaling_method["discontinuous"],
            resampling_upscaling=self.upscaling_method["discontinuous"],
            compatibility_resampling_downscaling=riowp.Resampling.nearest,
            compatibility_resampling_upscaling=riowp.Resampling.nearest,
            all_or_none_missing=all_or_none_missing,
            initialize_default=initialize_default,
            warning_point_data=warning_point_data,
            mod_func=mod_func,
            **kwargs,
        )

    def read_buildings_2d(self) -> ma.MaskedArray:
        """Read building height.

        Returns:
            Building height raster data.
        """
        return self._read_named_discontinuous()

    def read_building_id(self) -> ma.MaskedArray:
        """Read building ID.

        Returns:
            Building ID raster data.
        """
        return self._read_named_categorical()

    def read_building_type(self) -> ma.MaskedArray:
        """Read building type.

        Default values are applied if necessary.

        Returns:
            Building type raster data.
        """
        return self._read_named_categorical(initialize_default=True)

    def read_bridges_2d(self) -> ma.MaskedArray:
        """Read bridge height.

        Returns:
            Bridge height raster data.
        """
        return self._read_named_discontinuous()

    def read_bridges_id(self) -> ma.MaskedArray:
        """Read bridge ID.

        Returns:
            Bridge ID raster data.
        """
        return self._read_named_categorical()

    def read_lai(self) -> ma.MaskedArray:
        """Read Leaf Area Index.

        Returns:
            Leaf Area Index raster data.
        """
        return self._read_named_discontinuous()

    def read_tree_crown_diameter(self) -> ma.MaskedArray:
        """Read tree crown diameter.

        Returns:
            Tree crown diameter raster data.
        """
        return self._read_named_discrete()

    def read_tree_height(self) -> ma.MaskedArray:
        """Read tree height.

        Returns:
            Tree height raster data.
        """
        return self._read_named_discrete()

    def read_tree_trunk_diameter(self) -> ma.MaskedArray:
        """Read tree trunk diameter.

        Returns:
            Tree trunk diameter raster data.
        """
        return self._read_named_discrete()

    def read_tree_type(self) -> ma.MaskedArray:
        """Read tree type.

        Returns:
            Tree type raster data.
        """
        return self._read_named_categorical(warning_point_data=True)

    def read_pavement_type(self) -> ma.MaskedArray:
        """Read pavement type.

        Returns:
            Pavement type raster data.
        """
        return self._read_named_categorical()

    def read_patch_height(self) -> ma.MaskedArray:
        """Read vegetation patch height.

        Returns:
            Vegetation patch height raster data.
        """
        return self._read_named_discontinuous()

    def read_patch_type(self) -> ma.MaskedArray:
        """Read vegetation patch type.

        Returns:
            Vegetation patch type raster data.
        """
        return self._read_named_categorical()

    def read_street_crossings(self) -> ma.MaskedArray:
        """Read street crossing.

        Returns:
            Street crossing raster data.
        """
        return self._read_named_categorical()

    def read_street_type(self) -> ma.MaskedArray:
        """Read street type.

        Returns:
            Street type raster data.
        """
        return self._read_named_categorical()

    def read_soil_type(self) -> ma.MaskedArray:
        """Read soil type.

        Returns:
            Soil type raster data.
        """
        return self._read_named_categorical(initialize_default=True)

    def read_water_temperature(self) -> ma.MaskedArray:
        """Read water temperature.

        Returns:
            Water temperature raster data.
        """
        return self._read_named_continuous()

    def read_water_type(self) -> ma.MaskedArray:
        """Read water type.

        Returns:
            Water type raster data.
        """
        return self._read_named_categorical()

    def read_vegetation_height(self) -> ma.MaskedArray:
        """Read vegetation height.

        Returns:
            Vegetation height raster data.
        """
        return self._read_named_discontinuous()

    def read_vegetation_on_roofs(self) -> ma.MaskedArray:
        """Read vegetation on roof.

        Returns:
            Vegetation on roof raster data.
        """
        return self._read_named_categorical()

    def read_vegetation_type(self) -> ma.MaskedArray:
        """Read vegetation type.

        Returns:
            Vegetation type raster data.
        """
        return self._read_named_categorical()

    def read_zt(self) -> ma.MaskedArray:
        """Read terrain height.

        Returns:
            Terrain height raster data.
        """
        return self._read_named_continuous(
            all_or_none_missing=True,
            initialize_default=True,
        )

    def read_lcz(self, lcz_types: LCZTypes) -> ma.MaskedArray:
        """Read Local Climate Zones data.

        If the file has 3 bands, assume it is a rgb raster and convert it to lcz index.

        Args:
            lcz_types: LCZ types.

        Returns:
            LCZ index raster data.
        """

        def lcz_raster_to_index(raster: ma.MaskedArray, lcz_types: LCZTypes) -> ma.MaskedArray:
            """If `raster` has 3 bands, assume it is a rgb raster and convert it to lcz index."""
            if len(raster.shape) == 3 and raster.shape[0] >= 3:
                # assume rgb values that need to be converted to lcz index
                return lcz_types.lcz_rgb_to_index(raster)

            return raster

        return self._read_named_categorical(
            mod_func=lcz_raster_to_index,
            lcz_types=lcz_types,
        )


def _calling_read_variable_name() -> str:
    """Get the name of to read variable from the calling function.

    Assumes the calling function is named read_<variable_name>.

    Returns:
        Variable name.
    """
    return inspect.stack()[2].function.replace("read_", "")


def _enum_to_str_array(enum: Type[Enum]) -> npt.NDArray[np.str_]:
    """Convert an enumeration to an array of lowercase string values.

    Args:
        enum: The enumeration type to convert.

    Returns:
        Array containing the lowercase string values of the enumeration.
    """
    return np.array([element.name.lower() for element in enum])


def _find_variable_name(nc_data: Dataset, ndim: int) -> Variable:
    """Find the variable of the input data set with the given number of dimensions.

    Exclude dimension variables. Assume the files includes just one suitable variable.

    Args:
        nc_data: netCDF data set.
        ndim: Variable with this number of dimensions to search for.

    Raises:
        ValueError: No suitable variable found.
        ValueError: More than one suitable variable found.

    Returns:
        Variable with the given number of dimensions.
    """
    dimension_names = list(nc_data.dimensions.keys())

    variables_correct_dim = []
    for name, variable in nc_data.variables.items():
        if len(variable.dimensions) == ndim and name not in dimension_names:
            variables_correct_dim.append(name)

    nfound = len(variables_correct_dim)
    if nfound == 0:
        raise ValueError(f"No suitable variable with {ndim} dimension found.")
    elif nfound > 1:
        raise ValueError(f"Found {nfound} suitable variables when expecting 1.")
    return nc_data.variables[variables_correct_dim[0]]
