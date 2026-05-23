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

"""Variables and methods for domains."""

import logging
from enum import Enum
from math import floor
from pathlib import Path
from typing import (
    Callable,
    Dict,
    List,
    Literal,
    Optional,
    Tuple,
    Type,
    TypedDict,
    TypeVar,
    Union,
    cast,
)

import geopandas as gpd
import numpy as np
import numpy.ma as ma
import numpy.ma.core as ma_core
import numpy.typing as npt
import pandas as pd
import pandas.api.extensions as pde
import pandas.api.typing as pdtypes
import rasterio.warp as riowp
from netCDF4 import Dataset, Variable

from palm_csd import StatusLogger
from palm_csd.constants import (
    INPUT_DATA_INFO,
    NBUILDING_SURFACE_LAYER,
    NGHOST_POINTS,
    ColumnGeneral,
    FillValue,
    IndexBuildingGeneralPars,
    IndexBuildingIndoorPars,
    IndexBuildingSurfaceLevel,
    IndexBuildingSurfaceType,
    IndexPavementType,
    IndexVegetationType,
    IndexWaterType,
    InputData,
    get_parent_input_data,
)
from palm_csd.csd_config import (
    CSDConfig,
    CSDConfigAttributes,
    CSDConfigDomain,
    CSDConfigInput,
    ScalingMethods,
    value_defaults,
)
from palm_csd.geo_converter import GeoConverter
from palm_csd.lcz import LCZTypes
from palm_csd.netcdf_data import (
    NCDFCoordinateReferenceSystem,
    NCDFDimension,
    NCDFVariable,
    remove_existing_file,
)
from palm_csd.tools import Node
from palm_csd.vegetation import CanopyGenerator, tree_defaults

# Module logger. In __init__.py, it is ensured that the logger is a StatusLogger. For type checking,
# do explicit cast.
logger = cast(StatusLogger, logging.getLogger(__name__))


class CSDDomain(Node["CSDDomain"]):
    """A domain that stores all its configurations, output dimensions and variables."""

    name: str
    """Name of the domain."""

    config: CSDConfigDomain
    """Domain configuration."""
    input_config: CSDConfigInput
    """Input configuration."""
    attributes: CSDConfigAttributes
    """Attributes configuration"""

    canopy_generator: CanopyGenerator
    """Canopy generator."""

    geo_converter: Optional[GeoConverter]
    """GeoConverter to handle geographic data and transformation."""

    file_output: Path
    """Output path."""

    input_surface_polygons: Optional[Dict[Path, gpd.GeoDataFrame]]
    """Input surface polygons."""

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

    buildings_2d_removed: Optional[npt.NDArray[np.bool_]] = None
    """2D buildings that were removed due to building free border."""

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
        # Initialize Node part of the domain including parent to allow tree structure.
        super().__init__(parent)

        self.name = name

        # configurations
        self.config = config.domain_dict[name]
        self.input_config = config.input_of_domain(self.name)

        self.attributes = config.attributes

        # Create CanopyGenerator with the LAD method and parameters from the configuration.
        self.canopy_generator = CanopyGenerator(
            method=config.settings.lad_method,
            alpha_Metal2003=config.settings.lad_alpha,
            beta_Metal2003=config.settings.lad_beta,
            z_max_rel_LM2004=config.settings.lad_z_max_rel,
            dz=self.config.dz,
            pixel_size=self.config.pixel_size,
            season=config.settings.season,
            height_rel_resolved_vegetation_lower_threshold=config.settings.height_rel_resolved_vegetation_lower_threshold,
            lai_tree_lower_threshold=config.settings.lai_tree_lower_threshold,
            remove_low_lai_tree=self.config.remove_low_lai_tree,
        )

        # converter for geo data
        if config.settings.epsg is None:
            self.geo_converter = None
        else:
            # Find root parent
            parent = self.get_parent()
            if parent is not None:
                parent_geoconverter = parent.geo_converter
                if parent_geoconverter is None:
                    raise ValueError("Parent domain has no geo converter")
                root_parent_geoconverter = self.find_root().geo_converter
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
        self.file_output = config.output.file_out.with_stem(
            f"{config.output.file_out.stem}_{self.name}"
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

        self.input_surface_polygons = self._read_vector_data("surfaces")

        self._initialize_dimensions()
        self._initialize_variables()

    def _read_vector_data(
        self, element_type: Literal["surfaces", "trees"]
    ) -> Optional[Dict[Path, gpd.GeoDataFrame]]:
        """Read vector data from input files and process them.

        The vector data is read according to the element type's input configuration. Only required
        columns are kept and renamed to the internal names. Surface type columns that include
        several types (e.g. vegetation and pavement) are expanded to multiple columns.

        Args:
            element_type: Type of the element to read, either "trees" or "surfaces".

        Returns:
            Dictionary with file names as keys and GeoDataFrames as values. If no input files are
            found, return None.
        """
        # Assign input vector files, shape type and coordinate columns based on element type.
        shape_type: Literal["Point", "Polygon"]
        if element_type == "trees":
            if "trees" not in self.input_config.files.keys():
                return None
            vector_files = self.input_config.files["trees"]
            shape_type = "Point"
            columns_coordinates = [
                ColumnGeneral.geometry,
                ColumnGeneral.x_index,
                ColumnGeneral.y_index,
            ]
        elif element_type == "surfaces":
            if "surfaces" not in self.input_config.files.keys():
                return None
            vector_files = self.input_config.files["surfaces"]
            shape_type = "Polygon"
            columns_coordinates = [ColumnGeneral.geometry]
        else:
            raise ValueError(f"Unknown element type {element_type} for reading vector data.")

        if self.geo_converter is None:
            raise ValueError("geo_converter not set.")

        # Get user-defined columns from input configuration.
        # All columns
        columns_input = [key for key in self.input_config.columns.keys()]
        # Only columns that should be expanded according to the value's dict
        columns_expand = [
            key for key, value in self.input_config.columns.items() if isinstance(value, Dict)
        ]

        # Create mapping from column names to types and IDs. This will be used to expand single
        # columns with potentially multiple values into multiple columns according to the user input
        # in columns_expand.
        # Target types
        column_to_surface_type = {
            InputData.pavement_type: IndexPavementType,
            InputData.vegetation_type: IndexVegetationType,
            InputData.water_type: IndexWaterType,
        }
        # Mapping from surface type member names to their respective type and ID, for example:
        # "bare_soil" -> "vegetation_type", 1
        member_to_type = {}
        member_to_id = {}
        for surface_type_column, surface_type_members in column_to_surface_type.items():
            member_to_type.update(
                {member.name: surface_type_column for member in surface_type_members}
            )
            member_to_id.update({member.name: member.value for member in surface_type_members})

        # Read vector files and process them.
        input_shapes = {}
        for vector_file in vector_files:
            if not isinstance(vector_file, Path):
                raise ValueError(f"{element_type} file must be a Path object.")
            shapes = self.geo_converter.read_shp_to_dst(
                vector_file, shape_type=shape_type, name=vector_file.stem
            )

            # Make all column names lowercase.
            shapes.columns = shapes.columns.str.lower()

            # Keep only columns that are keys in self.input_config.columns (case-insensitive)
            columns_to_keep = [col for col in shapes.columns if col in columns_input]
            columns_to_keep += columns_coordinates
            shapes = cast(gpd.GeoDataFrame, shapes[columns_to_keep])

            # Rename columns from what is in the input file to internal name according to the user
            # input in columns.
            columns_rename = {
                key: value
                for key, value in self.input_config.columns.items()
                if isinstance(value, str) and key in shapes.columns
            }
            shapes = cast(gpd.GeoDataFrame, shapes.rename(columns=columns_rename))

            # Expand columns to multiple surface types if needed.
            for column in columns_expand:
                if column not in shapes.columns:
                    continue
                # Convert content of column to expand to lowercase. This would only work on string
                # and potentially also on object columns. Catch AttributeError if failed. Content
                # should be int then, float might also work. TODO: Do we need to check this case?
                try:
                    shapes[column] = shapes[column].str.lower()
                except AttributeError:
                    pass
                # Apply user defined mapping to a specific surface type member (e.g. "bare_soil").
                shapes["type"] = shapes[column].map(self.input_config.columns[column])  # type: ignore
                # Get the surface type (e.g. "vegetation_type" for "bare_soil").
                shapes["surface_type"] = shapes["type"].map(member_to_type)
                # Get ID (e.g. 1 for "bare_soil").
                shapes["type_id"] = shapes["type"].map(member_to_id).astype(pd.Int8Dtype())

                # Create new columns for each surface type in column_to_surface_type.
                for surface_type in column_to_surface_type.keys():
                    surface_type_present = shapes["surface_type"] == surface_type
                    if surface_type_present.any():
                        shapes[surface_type] = shapes["type_id"].where(surface_type_present)

                shapes = shapes.drop(columns=["surface_type", "type", "type_id", column])

            input_shapes[vector_file] = shapes

        return input_shapes

    def _initialize_dimensions(self) -> None:
        """Initialize dimensions."""
        self.x = NCDFDimension(
            name="x",
            data_type="f4",
            standard_name="projection_x_coordinate",
            long_name="x",
            units="m",
        )
        self.y = NCDFDimension(
            name="y",
            data_type="f4",
            standard_name="projection_y_coordinate",
            long_name="y",
            units="m",
        )
        self.z = NCDFDimension(name="z", data_type="f4", long_name="z", units="m")

        self.nsurface_fraction = NCDFDimension(name="nsurface_fraction", data_type="i")

        self.building_general_par = NCDFDimension(
            name="building_general_par",
            data_type=str,
            values=_enum_to_str_array(IndexBuildingGeneralPars),
        )
        self.building_indoor_par = NCDFDimension(
            name="building_indoor_par",
            data_type=str,
            values=_enum_to_str_array(IndexBuildingIndoorPars),
        )
        self.building_surface_layer = NCDFDimension(
            name="building_surface_layer",
            data_type="i",
            values=np.arange(1, NBUILDING_SURFACE_LAYER + 1),
        )
        self.building_surface_level = NCDFDimension(
            name="building_surface_level",
            data_type=str,
            values=_enum_to_str_array(IndexBuildingSurfaceLevel),
        )
        self.building_surface_type = NCDFDimension(
            name="building_surface_type",
            data_type=str,
            values=_enum_to_str_array(IndexBuildingSurfaceType),
        )

        # TODO: Use IndexVegetationPars values when PALM accept that.
        self.nvegetation_pars = NCDFDimension(
            name="nvegetation_pars", data_type="i", values=ma.arange(0, 12)
        )
        self.nwater_pars = NCDFDimension(name="nwater_pars", data_type="i", values=np.arange(0, 7))

        self.zlad = NCDFDimension(
            name="zlad",
            data_type="f4",
            long_name="z coordinate for resolved vegetation",
            units="m",
        )

        self.nuc = NCDFDimension(name="nuc", data_type="i")
        self.streetdir = NCDFDimension(name="streetdir", data_type="i")
        self.z_uhl = NCDFDimension(name="z_uhl", data_type="f4")

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
            mandatory=False,
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

        if not (
            "x_utm" in self.input_config.files
            and "y_utm" in self.input_config.files
            and "lon" in self.input_config.files
            and "lat" in self.input_config.files
        ):
            if self.geo_converter is None:
                logger.critical_raise(
                    f"Not all coordinate inputs provided for domain {self.name}, "
                    + "but target coordinate reference system not defined by EPSG code "
                    + "to calculate coordinates.",
                )

    def overlaps(self, other: "CSDDomain") -> Optional[bool]:
        """Checks if this domains overlaps with another domain.

        Args:
            other: Other domain.

        Returns:
            True if the domains overlap, False if they do not overlap, None if geo_converters are
            not defined.
        """
        if self.geo_converter is None or other.geo_converter is None:
            return None

        # Boundary boxes considering ghost points
        self_bbox = self.geo_converter.boundary(nghost_points=NGHOST_POINTS)
        other_bbox = other.geo_converter.boundary(nghost_points=NGHOST_POINTS)

        # The two domains do not overlap when they are next to each other.
        return not (
            self_bbox[2] <= other_bbox[0]  # self_box right of other_box
            or self_bbox[0] >= other_bbox[2]  # self_box left of other_box
            or self_bbox[3] <= other_bbox[1]  # self_box below other_box
            or self_bbox[1] >= other_bbox[3]  # self_box above other_box
        )

    class DefaultMetadataVariable(TypedDict, total=False):
        """Helper dictionary for default metadata in variables."""

        data_type: str
        """Data type."""
        fill_value: float
        """Fill value."""
        file: Path
        """Input/Output file."""
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
        mandatory: bool = True,
    ) -> NCDFVariable:
        """Helper function that returns a float variable with some predefined attributes.

        Args:
            name: Name.
            dimensions: Dimensions.
            long_name: Long name.
            units: Unit.
            standard_name: Standard name. Defaults to None.
            lod: Level of detail. Defaults to None.
            add_spatial_metadata: Add coordinates, grid_mapping. Defaults to True.
            mandatory: Whether the variable is mandatory. Defaults to True.

        Returns:
            float variable with given and predefined attributes.
        """
        default_values: CSDDomain.DefaultMetadataVariable = {
            "data_type": "f4",
            "fill_value": FillValue.FLOAT,
            "file": self.file_output,
        }
        if add_spatial_metadata:
            default_values.update(
                {
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
            mandatory=mandatory,
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
        mandatory: bool = True,
    ) -> NCDFVariable:
        """Helper function that returns an int variables with some predefined attributes.

        Args:
            name: Name.
            dimensions: Dimensions.
            long_name: Long name.
            units: Unit.
            standard_name: Standard name. Defaults to None.
            lod: Level of detail. Defaults to None.
            add_spatial_metadata: Add coordinates, grid_mapping. Defaults to True.
            mandatory: Whether the variable is mandatory. Defaults to True.

        Returns:
            int variable with given and predefined attributes.
        """
        default_values: CSDDomain.DefaultMetadataVariable = {
            "data_type": "i",
            "fill_value": FillValue.INTEGER,
            "file": self.file_output,
        }
        if add_spatial_metadata:
            default_values.update(
                {
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
            mandatory=mandatory,
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
        mandatory: bool = True,
    ) -> NCDFVariable:
        """Helper function that returns a byte variables with some predefined attributes.

        Args:
            name: Name.
            dimensions: Dimensions.
            long_name: Long name.
            units: Unit.
            standard_name: Standard name. Defaults to None.
            lod: Level of detail. Defaults to None.
            add_spatial_metadata: Add coordinates, grid_mapping. Defaults to True.
            mandatory: Whether the variable is mandatory. Defaults to True.

        Returns:
            byte variable with given and predefined attributes.
        """
        default_values: CSDDomain.DefaultMetadataVariable = {
            "data_type": "b",
            "fill_value": FillValue.BYTE,
            "file": self.file_output,
        }
        if add_spatial_metadata:
            default_values.update(
                {
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
            mandatory=mandatory,
            **default_values,
        )

    def is_resolved_vegetation(self) -> Union[np.bool_, npt.NDArray[np.bool_]]:
        """Array indicating if the grid cell includes resolved vegetation.

        Both, LAD and BAD values are checked.

        Returns:
            Boolean 3d array with True if grid cell includes resolved vegetation.
        """
        is_rv: Union[np.bool_, npt.NDArray[np.bool_]] = np.False_
        if self.lad.values is not None:
            is_rv = is_rv | ~self.lad.values.mask
        if self.bad.values is not None:
            is_rv = is_rv | ~self.bad.values.mask

        return is_rv

    def is_resolved_vegetation2d(self) -> Union[np.bool_, npt.NDArray[np.bool_]]:
        """Array indicating if the column above each pixel includes resolved vegetation.

        Returns:
            Boolean 2d array with True if column above each pixel includes resolved vegetation.
        """
        return self.is_resolved_vegetation().any(axis=0)

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

        The file is openend and closed. If file is None, self.input_config.files["x_utm"] is used.

        Args:
            file: Input file. Defaults to None.
            varname: Variable name. Defaults to None.

        Raises:
            ValueError: Both file and input_config.files["x_utm"] are None.

        Returns:
            Coordinate reference system from file.
        """
        if file is not None:
            from_file = file
        elif "x_utm" in self.input_config.files:
            from_file = self.input_config.files["x_utm"][0]
        else:
            raise ValueError("file or input_config.files['x_utm'] needs to be not None")

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

    MaskedArrayOrSeries = TypeVar("MaskedArrayOrSeries", pd.Series, ma.MaskedArray)

    def _validate_and_process_values(
        self,
        name: str,
        data: MaskedArrayOrSeries,
        all_or_none_missing: bool,
        initialize_default: bool,
    ) -> MaskedArrayOrSeries:
        """Validate and process values for a given column or field.

        Args:
            name: Name of the column or field.
            data: Data to validate and process, either as a pandas Series or a masked array.
            all_or_none_missing: If True, all values must be either missing or defined.
            initialize_default: If True, missing values are replaced by the default value.

        Raises:
            ValueError: Minimum or maximum are not defined.
            ValueError: Found values below or above minimum or maximum and
              replace_invalid_input_values is False.
            ValueError: Found both missing and defined values in the data when
              all_or_none_missing is True.
            ValueError: No default value defined for missing values when
              initialize_default is True.

        Returns:
            Checked and processed data, with the same type as the input data.
        """
        # Default values and allowed minimum and maximum values
        matches = [vd for vd in value_defaults.keys() if name.startswith(vd)]
        if len(matches) == 0:
            raise ValueError(f"No default values defined for {name}.")
        if len(matches) > 1:
            raise ValueError(f"Multiple default values defined for {name}: {matches}.")
        default_min_max = value_defaults[matches[0]]
        if default_min_max.minimum is None:
            raise ValueError(f"Minimum for {name} is not set.")
        if default_min_max.maximum is None:
            raise ValueError(f"Maximum for {name} is not set.")

        # In a masked array, ma.masked is assigned to indicate a masked value. In the Series,
        # missing values should be np.nan for float columns and pd.NA for integer columns. They are
        # converted to the fitting type. NOTE: The type below will be checked when potentielly
        # rasterizing the polygons. Adjust there as well.
        na_value: Union[float, pdtypes.NAType, ma_core.MaskedConstant]
        astype: Union[pde.ExtensionDtype, npt.DTypeLike]
        if isinstance(data, ma.MaskedArray):
            na_value = ma.masked
        else:
            if isinstance(default_min_max.minimum, int) and isinstance(
                default_min_max.maximum, int
            ):
                na_value = pd.NA
                if default_min_max.minimum < -128 or default_min_max.maximum > 127:
                    astype = pd.Int32Dtype()
                else:
                    astype = pd.Int8Dtype()
            else:
                na_value = np.nan
                astype = np.float64
            data = pd.to_numeric(data).astype(astype)

        below_minimum = data < default_min_max.minimum
        above_maximum = data > default_min_max.maximum
        replacement = default_min_max.default if default_min_max.default is not None else na_value

        # Check for values below minimum and above maximum, and deal with them.
        if below_minimum.any():
            if not self.replace_invalid_input_values:
                logger.critical(
                    f"In {name}, {below_minimum.sum()} values are "
                    + f"smaller than minimum value {default_min_max.minimum}.\n"
                    + "Enable replace_invalid_input_values for automatic replacement.\n"
                    + "Alternatively, adjust the minimum defined in "
                    + "palm_csd/data/value_defaults.csv."
                )
                raise ValueError(f"Invalid input values found in {name}.")
            logger.warning(
                f"In {name}, replacing {below_minimum.sum()} values "
                + f"smaller than minimum value {default_min_max.minimum} "
                + f"by default {replacement}."
            )
            logger.debug_indent(
                "If this is unintended, disable replace_invalid_input_values or "
                + "adjust the minimum in palm_csd/data/value_defaults.csv."
            )
            data[below_minimum] = replacement
        if above_maximum.any():
            if not self.replace_invalid_input_values:
                logger.critical(
                    f"In {name}, {above_maximum.sum()} values are "
                    + f"larger than maximum value {default_min_max.maximum}.\n"
                    + "Enable replace_invalid_input_values for automatic replacement.\n"
                    + "Alternatively, adjust the maximum defined in "
                    + "palm_csd/data/value_defaults.csv."
                )
                raise ValueError(f"Invalid input values found in {name}.")
            logger.warning(
                f"In {name}, replacing {above_maximum.sum()} values "
                + f"larger than maximum value {default_min_max.maximum} "
                + f"by default {replacement}."
            )
            logger.debug_indent(
                "If this is unintended, disable replace_invalid_input_values or "
                + "adjust the maximum in palm_csd/data/value_defaults.csv."
            )
            data[above_maximum] = replacement

        # Check for missing values and deal with them.
        if isinstance(data, ma.MaskedArray):
            na_values = data.mask
        else:
            na_values = data.isna()
        if na_values.any():
            if all_or_none_missing and not na_values.all():
                logger.critical_raise(
                    f"Found both missing and defined values in {name}.\n"
                    + f"Missing values in {name} are only allowed if all values are missing.",
                )
            if initialize_default:
                if default_min_max.default is None:
                    logger.critical_raise(
                        f"No default value defined for {name} to replace missing values.\n"
                        + "Set a default in palm_csd/data/value_defaults.csv.",
                    )
                logger.debug(
                    f"In column {name}, replacing {na_values.sum()} "
                    + f"missing values by default {default_min_max.default}."
                )
                data[na_values] = default_min_max.default

        return data

    def _select_check_vector_data(
        self,
        name: str,
        column_types: List[str],
        shape_input: Optional[Dict[Path, gpd.GeoDataFrame]] = None,
        all_or_none_missing_list: Optional[List[bool]] = None,
        initialize_default_list: Optional[List[bool]] = None,
        keep_shape_without_data: bool = False,
        mod_func_before: Optional[Callable[..., gpd.GeoDataFrame]] = None,
        mod_func_after: Optional[Callable[..., gpd.GeoDataFrame]] = None,
        **kwargs,
    ) -> Optional[gpd.GeoDataFrame]:
        """Process read-in vector data: Select columns, check values, and apply modifications.

        Use the supplied shape_input or the input_polygons from the domain. For each input shape,
        use only the columns of interest. If keep_shape_without_data is True, shapes (rows) are kept
        even if all other columns do not have a defined value (e.g. trees are defined by the mere
        presence of a point element, even if other data in the column is missing). If False, shapes
        are dropped if all other columns do not have a defined value.

        Minimum and maximum values are checked with the values from defaults. If
        replace_invalid_input_values is True, out of range values are replaced by the default value.
        If False, an error is raised. If initialize_default is True, missing values are replaced by
        the default value.

        Before or after checking the values, a modification function can be applied to the
        GeoDataFrame.

        Args:
            name: Variable name.
            column_types: List of ColumnVectorData enum values to read. These are used to
              determine the columns to read from the input GeoDataFrame.
            shape_input: List of input polygons to read. If None, the polygons from the domain's
              input_config are used. Defaults to None.
            all_or_none_missing_list: List with elements if true, raise an error if both missing and
              defined values are found. Defaults to False for each column type.
            initialize_default_list: List with elements if true, replace masked values by the
              default. Defaults to False for each column type.
            keep_shape_without_data: If True, keep shapes (rows) in the input, even if all other
              columns do not have a defined value. Defaults to False.
            mod_func_before: Function to modify the data before checking. Defaults to None.
            mod_func_after: Function to modify the data before checking. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func_before and mod_func_after.

        Raises:
            ValueError: geo_converter not set.
            ValueError: Read data not within range and replace_invalid_input_values is False.
            ValueError: No default value defined when replacement is necessary.

        Returns:
            Read and modified 2d raster data.
        """
        # If shape_input not explicitely set, use the domain's input polygons.
        if shape_input is None:
            if self.input_surface_polygons is None:
                return None
            shape_input = self.input_surface_polygons

        # Generate defaults when None is given.
        if all_or_none_missing_list is None:
            all_or_none_missing_list = [False] * len(column_types)
        if initialize_default_list is None:
            initialize_default_list = [False] * len(column_types)

        # Process the input polygons. Keep only columns of interest and drop NaN rows. Concatenate
        # to one GeoDataFrame.
        inputs = []
        for input_path, input_polygon in shape_input.items():
            # Columns to get from the input polygon that are present in the input polygon, and
            # columns including coordinates and geometry.
            existing_columns = [col for col in column_types if col in input_polygon.columns]
            columns_coordinates = [
                col
                for col in [ColumnGeneral.x_index, ColumnGeneral.y_index, ColumnGeneral.geometry]
                if col in input_polygon.columns
            ]
            columns_to_keep = existing_columns + columns_coordinates

            # If keep_shape_without_data, keep shaped even when only coordinates are present.
            # Otherwise, at least one of the data columns must be present.
            if keep_shape_without_data:
                filtered_polygon = input_polygon.dropna(subset=columns_to_keep, how="all")
            else:
                filtered_polygon = input_polygon.dropna(subset=existing_columns, how="all")
            if not filtered_polygon.empty:
                self.input_config.add_used_file(input_path)
                inputs.append(filtered_polygon[columns_to_keep])
        if not inputs:
            return None
        input = pd.concat(inputs, ignore_index=True)

        if not isinstance(input, gpd.GeoDataFrame):
            raise ValueError(f"Input {name} is not a GeoDataFrame.")

        # Apply modification function if supplied.
        if mod_func_before is not None:
            input = mod_func_before(input, **kwargs)

        # Check remaining columns.
        for column_index in range(len(column_types)):
            column_type = column_types[column_index]

            # Do not process "name" columns.
            if column_type.endswith("name"):
                if column_type not in input.columns:
                    # Assign missing values to string column.
                    input[column_type] = pd.Series(dtype=pd.StringDtype())
                continue

            # Check column. If it is not present, generate a new column with NaN values for the
            # check function.
            input[column_type] = self._validate_and_process_values(
                column_type,
                input.get(column_type, default=pd.Series([np.nan] * len(input))),
                all_or_none_missing_list[column_index],
                initialize_default_list[column_index],
            )

        # Apply modification function if supplied.
        if mod_func_after is not None:
            input = mod_func_after(input, **kwargs)

        return input

    def _read_check_raster_data(
        self,
        name: str,
        all_or_none_missing: bool = False,
        initialize_default: bool = False,
        resampling_downscaling: riowp.Resampling = riowp.Resampling.nearest,
        resampling_upscaling: riowp.Resampling = riowp.Resampling.nearest,
        compatibility_resampling_downscaling: Optional[riowp.Resampling] = None,
        compatibility_resampling_upscaling: Optional[riowp.Resampling] = None,
        warning_point_data: bool = False,
        mod_func_before: Optional[Callable[..., ma.MaskedArray]] = None,
        mod_func_after: Optional[Callable[..., ma.MaskedArray]] = None,
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
            mod_func_before: Function to modify the data before checking. Defaults to None.
            mod_func_after: Function to modify the data after checking. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func_before and mod_func_after.

        Raises:
            ValueError: Non netCDF file and geo_converter not set.
            ValueError: Read data not within range and replace_invalid_input_values is False.
            ValueError: No default value defined when replacement is necessary.

        Returns:
            Read and modified 2d raster data.
        """
        # Get file from input_config and mark it as read.
        file_list = self.input_config.files.get(name)
        if file_list is not None:
            file = file_list[0]
            self.input_config.add_used_file(file)
        else:
            file = None

        # If input_file is a netcdf file, read it directly; this handles also None.
        if file is None or file.suffix == ".nc":
            raster_values = self.read_nc_2d(file)

        # Otherwise, assume other GIS raster formats.
        else:
            if self.geo_converter is None:
                raise ValueError("geo_converter not set.")
            raster_values = self.geo_converter.read_raster_to_dst(
                file,
                resampling_downscaling=resampling_downscaling,
                resampling_upscaling=resampling_upscaling,
                compatibility_resampling_downscaling=compatibility_resampling_downscaling,
                compatibility_resampling_upscaling=compatibility_resampling_upscaling,
                warning_point_data=warning_point_data,
                name=name,
            )

            # Apply modification function if supplied.
            if mod_func_before is not None:
                raster_values = mod_func_before(raster_values, **kwargs)
            else:
                raster_values = raster_values[0, :, :]

            # Flip raster vertically to convert from GIS to netcdf convention.
            raster_values = ma.MaskedArray(np.flipud(raster_values))

        # Check values.
        raster_values = self._validate_and_process_values(
            name,
            raster_values,
            all_or_none_missing,
            initialize_default,
        )

        if mod_func_after is not None:
            raster_values = mod_func_after(raster_values, **kwargs)

        return raster_values

    def _rasterize_columns(
        self,
        name: str,
        input: gpd.GeoDataFrame,
        column_types: List[str],
    ) -> List[ma.MaskedArray]:
        if self.geo_converter is None:
            raise ValueError("geo_converter not set.")

        rasterized: List[ma.MaskedArray] = []
        dtype_rasterization: npt.DTypeLike
        for column_type in column_types:
            if input[column_type].dtype == pd.Int32Dtype():
                dtype_rasterization = np.int32
            elif input[column_type].dtype == pd.Int8Dtype():
                dtype_rasterization = np.int8
            elif input[column_type].dtype == np.float64:
                dtype_rasterization = np.float64
            else:
                raise ValueError(
                    f"Unexpected data type {input[column_type].dtype} for column {column_type}."
                )

            rasterized.append(
                ma.MaskedArray(
                    np.flipud(
                        self.geo_converter.rasterize(
                            input, column_type, dtype=dtype_rasterization, name=name
                        )
                    )
                )
            )
        return rasterized

    def _select_check_rasterize_vector_data(
        self,
        name: str,
        column_types: List[str],
        initialize_default_list: Optional[List[bool]] = None,
        mod_func_before: Optional[Callable[..., gpd.GeoDataFrame]] = None,
        mod_func_after: Optional[Callable[..., gpd.GeoDataFrame]] = None,
        **kwargs,
    ) -> Optional[List[ma.MaskedArray]]:
        """Process read-in vector data and rasterize it.

        Args:
            name: Name of the variable to read.
            column_types: List of ColumnVectorData enum values to read. These are used to
              determine the columns to read from the input GeoDataFrame.
            initialize_default_list: List of booleans, if true, replace masked values by the
              default. Defaults to None.
            mod_func_before: Function to modify the data before checking. Defaults to None.
            mod_func_after: Function to modify the data after checking. Defaults to None.
            **kwargs: Additional keyword arguments for mod_func_before and mod_func_after.

        Returns:
            List of rasterized columns as masked arrays, or None if no data is found.
        """
        vector = self._select_check_vector_data(
            name=name,
            column_types=column_types,
            initialize_default_list=initialize_default_list,
            mod_func_before=mod_func_before,
            mod_func_after=mod_func_after,
            **kwargs,
        )
        if vector is None:
            return None
        return self._rasterize_columns(name, vector, column_types)

    @staticmethod
    def _create_ids(df: gpd.GeoDataFrame, id_column: str) -> gpd.GeoDataFrame:
        """Create unique IDs.

        Args:
            df: Vector data frame.
            id_column: Column name for IDs.

        Returns:
            Data frame with unique IDs.
        """
        max_id = df[id_column].max(skipna=True)
        if pd.isna(max_id):
            max_id = 0
        else:
            max_id = int(max_id)

        nan_indices = df[id_column].isna()
        df.loc[nan_indices, id_column] = np.arange(max_id + 1, max_id + 1 + nan_indices.sum())
        return df

    def read(self, input_data: str, **kwargs) -> ma.MaskedArray:
        """Read input data, process it and return it as a masked array.

        Args:
            input_data: The input data to read.
            **kwargs: Additional arguments for mod_func.

        Raises:
            ValueError: If the input data is not a member of InputData.
            ValueError: If the input data is not supported.

        Returns:
            A masked array containing the read data.
        """
        if input_data == InputData.tree_type_name:
            raise ValueError(f"{input_data} can only be read in tree vector data.")

        if input_data == InputData.lcz:

            def _lcz_raster_to_index(raster: ma.MaskedArray, lcz_types: LCZTypes) -> ma.MaskedArray:
                """If raster has 3 bands, assume it is a rgb raster and convert it to lcz index."""
                if len(raster.shape) == 3 and raster.shape[0] >= 3:
                    # assume rgb values that need to be converted to lcz index
                    return lcz_types.lcz_rgb_to_index(raster)
                return raster[0, :, :]

            mod_func_raster_before = _lcz_raster_to_index
        else:
            mod_func_raster_before = None

        # input_data could refer to a specific layer or level, e.g. heat_capacity_wall_gfl_3. Get
        # general reading properties from parent input data, e.g. heat_capacity.
        parent_input_data = get_parent_input_data(input_data)
        all_or_none_missing = INPUT_DATA_INFO[parent_input_data].all_or_none_missing
        interpolation_type = INPUT_DATA_INFO[parent_input_data].type
        initialize_default = INPUT_DATA_INFO[parent_input_data].initialize_default
        warning_point_data = INPUT_DATA_INFO[parent_input_data].warning_point_data

        vector_data = self._select_check_rasterize_vector_data(
            name=input_data,
            column_types=[input_data],
            all_or_none_missing_list=[all_or_none_missing],
            initialize_default_list=[initialize_default],
        )
        if vector_data is not None:
            return vector_data[0]
        else:
            return self._read_check_raster_data(
                name=input_data,
                resampling_downscaling=self.downscaling_method[interpolation_type],
                resampling_upscaling=self.upscaling_method[interpolation_type],
                compatibility_resampling_downscaling=riowp.Resampling.nearest,
                compatibility_resampling_upscaling=riowp.Resampling.nearest,
                all_or_none_missing=all_or_none_missing,
                initialize_default=initialize_default,
                warning_point_data=warning_point_data,
                mod_func_before=mod_func_raster_before,
                **kwargs,
            )

    def read_buildings(self) -> Optional[List[ma.MaskedArray]]:
        """Read building shape data.

        Returns:
            Building data.
        """
        return self._select_check_rasterize_vector_data(
            name="buildings",
            column_types=[
                InputData.buildings_2d,
                InputData.building_id,
                InputData.building_type,
            ],
            initialize_default_list=[False, False, True],
            mod_func_after=self._create_ids,
            id_column=InputData.building_id,
        )

    def read_bridges(self) -> Optional[List[ma.MaskedArray]]:
        """Read bridge shape data.

        Returns:
            Bridge data.
        """
        return self._select_check_rasterize_vector_data(
            name="bridges",
            column_types=[
                InputData.bridges_2d,
                InputData.bridges_id,
            ],
            initialize_default_list=[False, False],
            mod_func_after=self._create_ids,
            id_column=InputData.bridges_id,
        )

    def read_trees(self) -> Optional[gpd.GeoDataFrame]:
        """Read tree shape data.

        Returns:
            Tree data.
        """
        shape_trees = self._read_vector_data("trees")
        if shape_trees is None:
            return None

        def type_name_to_type(points: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
            """Convert tree type name to tree type.

            Args:
                points: Tree data.

            Returns:
                Tree data with tree type.
            """
            if InputData.tree_type_name not in points.columns:
                return points

            type_names = points[InputData.tree_type_name].values
            types = []

            for type_name in type_names:
                if pd.isna(type_name) or type_name == "":
                    types.append(0)
                    continue
                type_match = [
                    index
                    for index, tree in enumerate(tree_defaults)
                    if type_name.startswith(tree.species)
                ]
                if len(type_match) == 1:
                    types.append(type_match[0])
                elif len(type_match) > 1:
                    logger.debug(
                        f"Multiple tree types found for {type_name}.\n" + "Using the first one."
                    )
                    types.append(type_match[0])
                elif len(type_match) == 0:
                    logger.debug(f"No tree type found for {type_name}.\n" + "Using the default.")
                    types.append(0)

            points[InputData.tree_type] = np.where(
                points[InputData.tree_type].isna(), types, points[InputData.tree_type]
            )
            return points

        return self._select_check_vector_data(
            name="trees",
            shape_input=shape_trees,
            column_types=[
                InputData.tree_crown_diameter,
                InputData.tree_height,
                InputData.tree_lai,
                InputData.tree_shape,
                InputData.tree_trunk_diameter,
                InputData.tree_type,
                InputData.tree_type_name,
            ],
            keep_shape_without_data=True,
            mod_func_after=type_name_to_type,
        )


def _enum_to_str_array(enum: Type[Enum]) -> npt.NDArray[np.str_]:
    """Convert an enumeration to an array of string values.

    Args:
        enum: The enumeration type to convert.

    Returns:
        Array containing the lowercase string values of the enumeration.
    """
    return np.array([element.name for element in enum])


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
