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

"""GIS tools.

Note that PALM as well the rest of palm_csd considers the (0, 0) point to be the lower left corner,
while the GIS convention is to have the (0, 0) point in the upper left corner of the domain. This
module takes care of this difference.
"""

import logging
from math import floor
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Literal, Optional, Tuple, cast

import affine
import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pyproj
import rasterio as rio
import rasterio.features as riof
import rasterio.transform as riotf
import rasterio.warp as riowp
from numpy import ma
from shapely.geometry import Polygon

from palm_csd import StatusLogger
from palm_csd.constants import ColumnGeneral, FillValue
from palm_csd.csd_config import CSDConfigDomain, CSDConfigOutput, CSDConfigSettings

# module logger, StatusLogger already set in __init__.py so cast only for type checking
logger = cast(StatusLogger, logging.getLogger(__name__))

# Suppress Pyogrio logging messages.
pyogrio_logger = logging.getLogger("pyogrio")
pyogrio_logger.setLevel(logging.WARNING)


class GeoConverter:
    """Class for GIS data processing.

    It includes methods to convert between different coordinate systems, select target areas, do
    interpolation and create CRS information.

    This class stores all required information for GIS data procession to the destination domain. In
    particular, it stores the destination's Coordinate Reference System (CRS) and the transform
    matrix. The transform matrix is an affine transformation matrix that maps pixel locations in
    (row, col) coordinates to (x, y) spatial positions. The product of this matrix and (0, 0), the
    row and column coordinates of the upper left corner of the dataset, is the spatial position of
    the upper left corner. It also supports rotation.
    """

    domain_name: Optional[str]
    """Domain name."""

    dst_crs: rio.CRS
    """Destination Coordinate Reference System."""
    dst_transform: rio.Affine
    """Destination transform matrix."""
    crs_wgs84: ClassVar[rio.CRS] = rio.CRS.from_epsg(4326)
    """WGS84 Coordinate Reference System."""

    pixel_size: float
    """Pixel size in meters."""
    rotation_angle: float
    """Rotation angle in degrees."""

    dst_width: int
    """Width of the destination domain."""
    dst_height: int
    """Height of the destination domain."""

    lower_left_x: float
    """Distance (m) on x-axis between the lower-left corner of the nested domain and the lower-left
    corner of the root parent domain."""
    lower_left_y: float
    """Distance (m) on y-axis between the lower-left corner of the nested domain and the lower-left
    corner of the root parent domain."""

    debug: bool
    """Write out reprojected data for debugging."""
    debug_file_prefix: Path
    """Prefix for debug output files."""

    origin_x: float
    """x-coordinate of the left border of the lower-left grid point in the destination CRS."""
    origin_y: float
    """y-coordinate of the lower border of the lower-left grid point in the destination CRS."""
    origin_lon: float
    """Longitude of the left border of the lower-left grid point in WGS84."""
    origin_lat: float
    """Latitude of the lower border of the lower-left grid point in WGS84."""

    corner_ll: Tuple[float, float]
    """Coordinates of the lower-left corner of the lower-left grid point of the domain in the
    destination CRS."""
    corner_lr: Tuple[float, float]
    """Coordinates of the lower-right corner of the lower-right grid point of the domain in the
    destination CRS."""
    corner_ul: Tuple[float, float]
    """Coordinates of the upper-left corner of the upper-left grid point of the domain in the
    destination CRS."""
    corner_ur: Tuple[float, float]
    """Coordinates of the upper-right corner of the upper-right grid point of the domain in the
    destination CRS."""

    corner_ll_cr: Tuple[float, float]
    """Coordinates of the centre of the lower-left grid point of the domain in the destination
    CRS."""
    corner_lr_cr: Tuple[float, float]
    """Coordinates of the centre of the lower-right grid point of the domain in the destination
    CRS."""
    corner_ul_cr: Tuple[float, float]
    """Coordinates of the centre of the upper-left grid point of the domain in the destination
    CRS."""
    corner_ur_cr: Tuple[float, float]
    """Coordinates of the centre of the upper-right grid point of the domain in the destination
    CRS."""

    ignore_input_georeferencing: bool
    """Ignore input file's georeferencing."""

    input_x0: Optional[float]
    """x index of the lower-left corner in the input data."""
    input_y0: Optional[float]
    """y coordinate of the lower-left corner in the input data."""

    def __init__(
        self,
        domain_config: CSDConfigDomain,
        settings: CSDConfigSettings,
        output: CSDConfigOutput,
        parent: Optional["GeoConverter"] = None,
        root_parent: Optional["GeoConverter"] = None,
        domain_name: Optional[str] = None,
        debug_output: bool = False,
    ):
        """Initialize the GeoConverter instance.

        Args:
            domain_config: Configuration of the domain.
            settings: General setting.
            output: General output settings.
            parent: GeoConverter of the parent. Defaults to None.
            root_parent: GeoConverter of the root domain. Defaults to None.
            domain_name: Name of the domain. Defaults to None.
            debug_output: Write out reprojected data for debugging. Defaults to False.

        Raises:
            ValueError: EPSG code not defined in settings.
            ValueError: If ignore_input_georeferencing, input_lower_left_x or input_lower_left_y not
              set.
            ValueError: Parent and child are not compatible.
            ValueError: origin_x/y or origin_lon/lat not set when required.
        """
        self.domain_name = domain_name

        if settings.epsg is None:
            raise ValueError("No EPSG code given.")
        self.dst_crs = rio.CRS.from_epsg(settings.epsg)

        self.pixel_size = domain_config.pixel_size
        self.rotation_angle = settings.rotation_angle

        self.dst_width = domain_config.nx + 1
        self.dst_height = domain_config.ny + 1

        # Ignore georeferencing of input data? If yes, input_lower_left_x/y in m have to be given
        # and we calculate the respective coordinates here. Otherwise, we don't need them.
        self.ignore_input_georeferencing = settings.ignore_input_georeferencing
        if self.ignore_input_georeferencing:
            if domain_config.input_lower_left_x is None or domain_config.input_lower_left_y is None:
                raise ValueError("input_lower_left_x and input_lower_left_y have to be given.")
            self.input_x0 = int(floor(domain_config.input_lower_left_x / domain_config.pixel_size))
            self.input_y0 = int(floor(domain_config.input_lower_left_y / domain_config.pixel_size))
        else:
            self.input_x0 = None
            self.input_y0 = None

        if parent is None:
            self.lower_left_x = 0.0
            self.lower_left_y = 0.0
            if domain_config.origin_x is not None and domain_config.origin_y is not None:
                self.origin_x = domain_config.origin_x
                self.origin_y = domain_config.origin_y
        else:
            if root_parent is None:
                raise ValueError("Root parent not given.")
            self._set_lower_left_with_parent(domain_config, parent, root_parent)

        # origin_x/y set above? Calculate consistent origin_lon/lat.
        if hasattr(self, "origin_x") and hasattr(self, "origin_y"):
            tmp_lon, tmp_lat = self.transform_points_to_wgs84([self.origin_x], [self.origin_y])
            self.origin_lon = tmp_lon[0]
            self.origin_lat = tmp_lat[0]
        # Calculate origin_x/y from origin_lon/lat.
        elif domain_config.origin_lon is not None and domain_config.origin_lat is not None:
            self.origin_lon = domain_config.origin_lon
            self.origin_lat = domain_config.origin_lat

            tmp_x, tmp_y = self.transform_points_from_wgs84(
                [domain_config.origin_lon], [domain_config.origin_lat]
            )
            self.origin_x = tmp_x[0]
            self.origin_y = tmp_y[0]
        else:
            logger.critical_indent_raise(
                "A complete set of origin_x/y or origin_lon/lat not given "
                + f"for domain {self.domain_name}.",
            )

        # Top left corner (border of the pixel) of the domain is needed for the transform.
        top_left_x = self.origin_x
        top_left_y = self.origin_y + self.dst_height * self.pixel_size

        # Unrotated transform.
        self.dst_transform = riotf.from_origin(
            top_left_x, top_left_y, self.pixel_size, self.pixel_size
        )

        # Rotation with rotation_angle (clockwise) around (x=0, y=self.dst_height).
        # These coordinates are relative to top left corner.
        rotation = affine.Affine.rotation(self.rotation_angle, (0, self.dst_height))
        self.dst_transform = self.dst_transform * rotation

        # Destination CRS coordinates of the corners of the corner pixel of the domain
        self.corner_ul = riotf.xy(self.dst_transform, 0, 0, offset="ul")
        self.corner_ur = riotf.xy(self.dst_transform, 0, self.dst_width - 1, offset="ur")
        self.corner_ll = riotf.xy(self.dst_transform, self.dst_height - 1, 0, offset="ll")
        self.corner_lr = riotf.xy(
            self.dst_transform, self.dst_height - 1, self.dst_width - 1, offset="lr"
        )
        # Destination CRS coordinates of the centres of the corner pixel of the domain
        self.corner_ul_cr = riotf.xy(self.dst_transform, 0, 0)
        self.corner_ur_cr = riotf.xy(self.dst_transform, 0, self.dst_width - 1)
        self.corner_ll_cr = riotf.xy(self.dst_transform, self.dst_height - 1, 0)
        self.corner_lr_cr = riotf.xy(self.dst_transform, self.dst_height - 1, self.dst_width - 1)

        # Debug settings for writing out the reprojected data.
        self.debug = debug_output
        self.debug_file_prefix = (
            output.file_out.parent / "verbose_gis_output" / output.file_out.stem
        )

    # TODO: Use Self with Python 3.11.
    def _set_lower_left_with_parent(
        self, domain_config: CSDConfigDomain, parent: "GeoConverter", root_parent: "GeoConverter"
    ) -> None:
        """Set lower_left_x/y and origin_x/y compatible with parent.

        Args:
            domain_config: Configuration of the domain.
            parent: GeoConverter of the parent. Defaults to None.
            root_parent: GeoConverter of the root domain. Defaults to None.

        Raises:
            ValueError: Parent and child are not compatible.
        """
        # Check compatibility of parent and child.
        if parent.dst_crs != self.dst_crs:
            raise ValueError("Parent and child have different CRS.")

        if parent.rotation_angle != self.rotation_angle:
            raise ValueError("Parent and child have different rotation angle.")

        if not parent.pixel_size % self.pixel_size == 0.0:
            logger.critical_indent_raise(
                f"Pixel size of parent domain {parent.domain_name} is not an integer multiple "
                + f"of the pixel size of child domain {self.domain_name}.",
            )

        if not (self.pixel_size * self.dst_width) % parent.pixel_size == 0.0:
            logger.critical_indent_raise(
                f"x size of child domain {self.domain_name} does not completely fill "
                + f"pixels of parent {parent.domain_name}.",
            )

        if not (self.pixel_size * self.dst_height) % parent.pixel_size == 0.0:
            logger.critical_indent_raise(
                f"y size of child domain {self.domain_name} does not completely fill "
                + f"pixels of parent {parent.domain_name}.",
            )

        # Child coordinates relative to root parent. If not given, calculate lower_left_x/y from
        # origin_x/y or origin_lon/lat. Further checks below.
        if domain_config.lower_left_x is not None and domain_config.lower_left_y is not None:
            self.lower_left_x = domain_config.lower_left_x
            self.lower_left_y = domain_config.lower_left_y
        else:
            # Calculate lower_left_x/y from origin_*.
            # Preliminay origin_x/y.
            if domain_config.origin_x is not None and domain_config.origin_y is not None:
                origin_x_prelim = domain_config.origin_x
                origin_y_prelim = domain_config.origin_y
            elif domain_config.origin_lon is not None and domain_config.origin_lat is not None:
                tmp_x, tmp_y = self.transform_points_from_wgs84(
                    [domain_config.origin_lon], [domain_config.origin_lat]
                )
                origin_x_prelim = tmp_x[0]
                origin_y_prelim = tmp_y[0]
            else:
                raise ValueError("Not all required input needed for geo conversion given.")

            # Preliminary lower left corner of child domain, possibly not compatible with its
            # parent.
            # origin_x/y -> lower_left_x/y: rotate around root_parent.origin_x/y
            #                               with -rotation_angle
            lower_left_x_prelim, lower_left_y_prelim = self._rotate(
                origin_x_prelim - root_parent.origin_x,
                origin_y_prelim - root_parent.origin_y,
                -self.rotation_angle,
            )

            # Calculate lower left corner of child domain compatible with its parent.
            self.lower_left_x = (
                np.round((lower_left_x_prelim - parent.lower_left_x) / parent.pixel_size)
                * parent.pixel_size
                + parent.lower_left_x
            )
            self.lower_left_y = (
                np.round((lower_left_y_prelim - parent.lower_left_y) / parent.pixel_size)
                * parent.pixel_size
                + parent.lower_left_y
            )

            logger.info_indent("Position relative the root parent domain:", hierarchy=1)
            logger.info_indent(f"lower_left_x: {self.lower_left_x}", hierarchy=2)
            logger.info_indent(f"lower_left_y: {self.lower_left_y}", hierarchy=2)

            if self.lower_left_x != lower_left_x_prelim or self.lower_left_y != lower_left_y_prelim:
                logger.warning_indent(
                    "Adjusted origin_x/y or origin_lon/lat to be compatible with parent.",
                    hierarchy=1,
                )

        # Check position relative to parent.
        self._check_clearance(self.lower_left_x, parent.lower_left_x, parent, "left x")
        self._check_clearance(self.lower_left_y, parent.lower_left_y, parent, "bottom y")
        self._check_clearance(
            parent.lower_left_x + parent.pixel_size * parent.dst_width,
            self.lower_left_x + self.pixel_size * self.dst_width,
            parent,
            "right x",
        )
        self._check_clearance(
            parent.lower_left_y + parent.pixel_size * parent.dst_height,
            self.lower_left_y + self.pixel_size * self.dst_height,
            parent,
            "top y",
        )

        # Calculate origin_x/y from lower_left_x/y.
        # lower_left_x/y -> origin_x/y: rotate with rotation_angle and
        #                               add root_parent.origin_x/y
        lower_left_x_rot, lower_left_y_rot = self._rotate(
            self.lower_left_x, self.lower_left_y, self.rotation_angle
        )
        self.origin_x = root_parent.origin_x + lower_left_x_rot
        self.origin_y = root_parent.origin_y + lower_left_y_rot

    def _check_clearance(
        self,
        lower_left_larger: float,
        lower_left_smaller: float,
        parent: "GeoConverter",
        side: str,
    ) -> None:
        clearence_message = (
            "The clearance must be at least 4 parent-grid cells if using the "
            + "Wicker-Skamarock advection scheme and at least 2 parent-grid cells if using the "
            + "Piaseck-Williams advection scheme."
        )

        # Check position relative to parent.
        diff_lower_left = lower_left_larger - lower_left_smaller
        diff_lower_left_ratio = int(diff_lower_left / parent.pixel_size)
        # Check if child domain is aligned with parent domain.
        if diff_lower_left % parent.pixel_size != 0.0:
            logger.critical_indent_raise(
                f"{side.capitalize()} border does not align with "
                + "grid of parent {parent.domain_name}.",
            )
        # Check if there is enough space at the border to the parent domain.
        if diff_lower_left_ratio < 2:
            logger.critical_indent_raise(
                f"Not enough space at the {side} border to parent domain {parent.domain_name} "
                + f"with {diff_lower_left_ratio} parent-grid "
                + f"{'cell' if diff_lower_left_ratio == 1 else 'cells'}.\n"
                + clearence_message,
            )
        if diff_lower_left_ratio < 4:
            logger.warning_indent(
                f"Not enough space at the {side} border to parent domain {parent.domain_name} "
                + "for Wicker-Skamarock advection scheme "
                + f"with {diff_lower_left_ratio} parent-grid cells.\n"
                + clearence_message,
            )

    @staticmethod
    def _transform_points(
        src_crs: rio.CRS,
        dst_crs: rio.CRS,
        x: npt.ArrayLike,
        y: npt.ArrayLike,
    ) -> Tuple[List[float], List[float]]:
        """Transform points from one CRS to another.

        Args:
            src_crs: Source Coordinate Reference System.
            dst_crs: Destination Coordinate Reference System.
            x: Source x-coordinates.
            y: Source y-coordinates.

        Returns:
            Destination coordinates.
        """
        x_transform, y_transform, *_ = riowp.transform(src_crs, dst_crs, x, y)
        return x_transform, y_transform

    @staticmethod
    def _rotate(x: float, y: float, angle: float) -> Tuple[float, float]:
        """Rotate point anti-clockwise.

        Args:
            x: Source x-coordinate.
            y: Source y-coordinate.
            angle: Angle in degrees.

        Returns:
            Rotated coordinates.
        """
        return affine.Affine.rotation(angle) * (x, y)

    def global_palm_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate x and y coordinates of the domain in the root parent domain.

        Returns:
            Global x and y coordinates in the root parent domain.
        """
        x_global = np.arange(0.5, self.dst_width + 0.5) * self.pixel_size + self.lower_left_x
        y_global = np.arange(0.5, self.dst_height + 0.5) * self.pixel_size + self.lower_left_y
        return x_global, y_global

    def boundary(self, nghost_points: int = 0) -> Tuple[float, float, float, float]:
        """Calculate the boundary of the domain in the root parent coordinate system.

        Args:
            nghost_points: Number of ghost points to consider. Defaults to 0.

        Returns:
            Lower left x, lower left y, upper right x, upper right y.
        """
        lower_left_gp_x = self.lower_left_x - nghost_points * self.pixel_size
        lower_left_gp_y = self.lower_left_y - nghost_points * self.pixel_size
        upper_right_gp_x = self.lower_left_x + (nghost_points + self.dst_width) * self.pixel_size
        upper_right_gp_y = self.lower_left_y + (nghost_points + self.dst_height) * self.pixel_size

        return lower_left_gp_x, lower_left_gp_y, upper_right_gp_x, upper_right_gp_y

    def geographic_coordinates(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Calculate x and y coordinates of the domain in the destination CRS and in WGS84.

        Returns:
            x-corrdinates in destination CRS, y-coordinates in destination CRS, longitudes in WGS84,
            latitudes in WGS84.
        """
        # Create meshgrid of coordinates and flip rows to get correct order.
        cols, rows = np.meshgrid(np.arange(self.dst_width), np.arange(self.dst_height))
        rows = np.flipud(rows)

        # x and y coordinates in destination CRS.
        xs, ys = riotf.xy(self.dst_transform, rows, cols)
        x_coord = np.reshape(xs, cols.shape)
        y_coord = np.reshape(ys, rows.shape)

        # Transform to wgs84 for longitudes and latitudes.
        lon, lat = self.transform_points_to_wgs84(x_coord.flatten(), y_coord.flatten())
        lon_coord = np.array(lon).reshape(self.dst_height, self.dst_width)
        lat_coord = np.array(lat).reshape(self.dst_height, self.dst_width)

        return x_coord, y_coord, lon_coord, lat_coord

    def transform_points_from_wgs84(
        self, x: npt.ArrayLike, y: npt.ArrayLike
    ) -> Tuple[List[float], List[float]]:
        """Transform points from WGS84 to destination CRS.

        Args:
            x: x-coordinates in WGS84.
            y: y-coordinates in WGS84.

        Returns:
            x-coordinates in destination CRS, y-coordinates in destination CRS.
        """
        x_transform, y_transform = self._transform_points(
            self.crs_wgs84,
            self.dst_crs,
            x,
            y,
        )
        return x_transform, y_transform

    def transform_points_to_wgs84(
        self, x: npt.ArrayLike, y: npt.ArrayLike
    ) -> Tuple[List[float], List[float]]:
        """Transform points from destination CRS to WGS84.

        Args:
            x: x-coordinates in destination CRS.
            y: y-coordinates in destination CRS.

        Returns:
            x-coordinates in WGS84, y-coordinates in WGS84.
        """
        x_transform, y_transform, *_ = self._transform_points(self.dst_crs, self.crs_wgs84, x, y)
        return x_transform, y_transform

    def dst_crs_to_cf(self) -> Dict[str, Any]:
        """CF compliant CRS information.

        This method uses a pyproj routine. epsg_code and units are added if not already present.

        Raises:
            ValueError: Unknown linear unit.

        Returns:
            CF compliant CRS information.
        """
        crs = pyproj.CRS(self.dst_crs).to_cf()
        # PALM wants the epsg_code and units, add it if not already present.
        if "epsg_code" not in crs:
            crs["epsg_code"] = f"EPSG:{self.dst_crs.to_epsg()}"
        if "units" not in crs:
            unit = self.dst_crs.linear_units
            if unit in ["metre", "meter", "m"]:
                crs["units"] = "m"
            else:
                raise ValueError(f"Unknown linear unit {unit}")
        return crs

    def read_shp_to_dst(
        self, file: Path, shape_type: Literal["Point", "Polygon"], name: Optional[str] = None
    ) -> gpd.GeoDataFrame:
        """Read a raster file and cut or reproject it to the destination domain.

        Open the input file and check its grid. If required, reproject the data to the destination
        grid. Cut the data to the destination grid.

        Args:
            file: File with raster data.
            shape_type: Type of the shape to read. Either "Point" or "Polygon".
            resampling_downscaling: Resampling method if downscaling is needed. Defaults to
              riowp.Resampling.nearest.
            resampling_upscaling: Resampling method if upscaling is needed. Defaults to
              riowp.Resampling.nearest.
            compatibility_resampling_downscaling: Masked values of this resampling method should be
              applied to the output when downscaling. Defaults to None.
            compatibility_resampling_upscaling: Masked values of this resampling method should be
              applied to the output when upscaling. Defaults to None.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            name: Variable name for debug output. Defaults to None.

        Returns:
            Raster data in the destination domain.
        """
        logger.debug(f"Reading {name} shape input.")
        # TODO: Cutting to the destination grid is only done later to avoid issues when
        # reprojecting. Doing it here might save some memory. Check if it is worth it.
        # NOTE: integer columns with NULL values are read in as floats.
        vector_input = gpd.read_file(file)

        # Convert multi-part geometries to multiple single geometries.
        vector_input = vector_input.explode()

        # Remove wrong geometries.
        is_shape_type = vector_input.geometry.type == shape_type
        if not is_shape_type.all():
            logger.warning(f"Removing non-{shape_type} geometries from {name}.")
            shapes = vector_input[is_shape_type]
        else:
            shapes = vector_input

        if not isinstance(shapes, gpd.GeoDataFrame):
            logger.critical_indent_raise(f"Input file {name} is not a GeoDataFrame.")

        output_str = "cut"
        if shapes.crs != self.dst_crs:
            logger.info_indent(f"Reprojecting {name} shape input due to different CRS.")
            shapes = shapes.to_crs(self.dst_crs)
            if shapes is None:
                logger.critical_indent_raise(f"Reprojection of {name} failed.")
            output_str = "reprojected"

        # Select points in the destination domain.
        if self.rotation_angle != 0.0:
            # Create a polygon from the domain corners
            domain_polygon = Polygon(
                [self.corner_ll, self.corner_lr, self.corner_ur, self.corner_ul]
            )

            # Select points that are contained within or cross the domain polygon
            shapes = shapes[shapes.geometry.intersects(domain_polygon)]
        else:
            # For non-rotated domains, keep using .cx for efficiency
            # See https://github.com/python/mypy/issues/2410 for the mypy issue
            shapes = shapes.cx[
                self.corner_ll[0] : self.corner_ur[0], self.corner_ll[1] : self.corner_ur[1]  # type: ignore[misc]
            ]

        # After processing, it should still be a GeoDataFrame. mypy wants this check.
        if not isinstance(shapes, gpd.GeoDataFrame):
            logger.critical_indent_raise(f"Input file {name} is not a GeoDataFrame anymore.")

        if shape_type == "Point":
            # Add coordinate indices to the points.
            # Warn if columns already exist.
            for coord_column in [ColumnGeneral.x_index, ColumnGeneral.y_index]:
                if coord_column in shapes.columns:
                    logger.warning(
                        f"{coord_column} is already a column in {name} shape input.\n"
                        + "It will be overwritten."
                    )
            if len(shapes) == 0:
                shapes[ColumnGeneral.x_index] = []
                shapes[ColumnGeneral.y_index] = []
            else:
                # Calculate the indices of the points in the destination domain. apply returns a
                # series of tuples, which are unpackage with *. zip aggregates the tuples into two
                # lists.
                y_indices, x_indices = np.array(
                    riotf.rowcol(self.dst_transform, shapes.geometry.x, shapes.geometry.y)
                )
                # Points directly on the right or bottom border are set to the last index.
                x_indices[x_indices == self.dst_width] = self.dst_width - 1
                y_indices[y_indices == self.dst_height] = self.dst_height - 1

                # All other points should be inside the domain because of the cx selection above.
                if np.any(x_indices < 0) or np.any(x_indices >= self.dst_width):
                    logger.critical_indent_raise(
                        f"Points in {name} shape input are outside of the destination domain in x."
                    )
                if np.any(y_indices < 0) or np.any(y_indices >= self.dst_height):
                    logger.critical_indent_raise(
                        f"Points in {name} shape input are outside of the destination domain in y."
                    )

                # Convert from top left origin to bottom left origin.
                y_indices = self.dst_height - 1 - y_indices

                # Update DataFrame.
                shapes[ColumnGeneral.x_index] = x_indices
                shapes[ColumnGeneral.y_index] = y_indices

        # Write out the reprojected data if debug is enabled.
        if self.debug:
            if name is None:
                raise ValueError("name has to be given.")

            if self.domain_name is None:
                raise ValueError("domain_name is not set.")

            self.debug_file_prefix.parent.mkdir(exist_ok=True)
            output = self.debug_file_prefix.with_stem(
                f"{self.debug_file_prefix.stem}_{name}-{output_str}_{self.domain_name}"
            ).with_suffix(".shp")
            shapes.to_file(output)
            logger.debug_indent(f"Result written to {output}.", hierarchy=1)

        return shapes

    def rasterize(
        self,
        vector_data: gpd.GeoDataFrame,
        column: str,
        dtype: npt.DTypeLike,
        name: Optional[str] = None,
    ) -> ma.MaskedArray:
        """Rasterize vector data.

        Args:
            vector_data: Vector data to rasterize.
            column: Column in the vector data to rasterize.
            dtype: Data type of the rasterized data.
            name: Variable name for debug output. Defaults to None.

        Returns:
            Rasterized data.
        """
        logger.debug(f"Rasterizing {column} shape input.")
        # Rasterize the vector data.
        shapes = zip(vector_data.geometry, vector_data[column])
        # TODO: generate masked array directly.
        rasterized = riof.rasterize(
            shapes,
            out_shape=(self.dst_height, self.dst_width),
            transform=self.dst_transform,
            fill=FillValue.from_dtype(dtype),
            all_touched=False,
            dtype=dtype,
        )
        rasterized_ma: ma.MaskedArray = ma.MaskedArray(
            rasterized,
            mask=(rasterized == FillValue.from_dtype(dtype)),
            fill_value=FillValue.from_dtype(dtype),
        )

        # Write out the reprojected data if debug is enabled.
        if self.debug:
            if name is None:
                raise ValueError("name has to be given.")

            if self.domain_name is None:
                raise ValueError("domain_name is not set.")

            self.debug_file_prefix.parent.mkdir(exist_ok=True)
            output_reprojected = self.debug_file_prefix.with_stem(
                f"{self.debug_file_prefix.stem}_{name}_{column}-rasterized_{self.domain_name}"
            ).with_suffix(".tif")
            self.write_dst_geotiff(output_reprojected, rasterized_ma)
            logger.debug_indent(f"Result written to {output_reprojected}.", hierarchy=1)

        return rasterized_ma

    def read_raster_to_dst(
        self,
        file: Path,
        resampling_downscaling: riowp.Resampling = riowp.Resampling.nearest,
        resampling_upscaling: riowp.Resampling = riowp.Resampling.nearest,
        compatibility_resampling_downscaling: Optional[riowp.Resampling] = None,
        compatibility_resampling_upscaling: Optional[riowp.Resampling] = None,
        warning_point_data: bool = False,
        name: Optional[str] = None,
    ) -> ma.MaskedArray:
        """Read a raster file and cut or reproject it to the destination domain.

        Open the input file and check its grid. If its CRS and grid align with the destination grid,
        cut the data to the destination grid. Otherwise reproject the data to the destination grid.

        Args:
            file: File with raster data.
            resampling_downscaling: Resampling method if downscaling is needed. Defaults to
              riowp.Resampling.nearest.
            resampling_upscaling: Resampling method if upscaling is needed. Defaults to
              riowp.Resampling.nearest.
            compatibility_resampling_downscaling: Masked values of this resampling method should be
              applied to the output when downscaling. Defaults to None.
            compatibility_resampling_upscaling: Masked values of this resampling method should be
              applied to the output when upscaling. Defaults to None.
            warning_point_data: Warn if single point data is reprojected. Defaults to False.
            name: Variable name for debug output. Defaults to None.

        Returns:
            Raster data in the destination domain.
        """
        src = rio.open(file)

        # Ignore georeferencing if requested and just cut the input to the destination grid.
        if self.ignore_input_georeferencing:
            logger.info(f"Cut {name} input to target grid ignoring georeferencing.")
            return self.cut_aligned_raster_to_dst(src, name)

        # Determine whether to use upscaling or downscaling resampling.
        resolution_src = self.resolution_in_dst(src)
        resolution_ratio = resolution_src / self.pixel_size
        if resolution_ratio < 0.9:
            resampling = resampling_upscaling
            compatibility_resampling = compatibility_resampling_upscaling
            resampling_str = "upscaling"
        else:
            resampling = resampling_downscaling
            compatibility_resampling = compatibility_resampling_downscaling
            resampling_str = "downscaling"
        debug_message = (
            f"Average input resolution is {resolution_src:.2f} m.\n"
            + f"Ratio of input/output resolutions is {resolution_ratio:.1f} "
            + f"so {resampling_str} resampling is applied."
        )

        # Reproject if CRSs are different. This already includes the correct spatial extensions.
        if src.crs != self.dst_crs:
            logger.info(
                f"Reprojecting {name} input due to different CRS "
                + f"using {resampling.name} resampling."
            )
            logger.debug_indent(debug_message)
            if warning_point_data:
                logger.warning_indent(
                    "Reprojection of point data may lead to unexpected results.\n"
                    + f"Consider supplying {name} in the output domain's CRS and grid."
                )
            return self.reproject_raster_to_dst(
                src,
                resampling=resampling,
                compatibility_resampling=compatibility_resampling,
                name=name,
            )

        # Reproject if resolution or rotation are different. This already includes the correct
        # spatial extensions.
        if not (
            np.isclose(src.transform.a, self.dst_transform.a)
            and np.isclose(src.transform.b, self.dst_transform.b)
            and np.isclose(src.transform.d, self.dst_transform.d)
            and np.isclose(src.transform.e, self.dst_transform.e)
        ):
            logger.info(
                f"Changing grid of {name} input due to different resolution or rotation "
                + f"using {resampling.name} resampling."
            )
            logger.debug_indent(debug_message)
            if warning_point_data:
                logger.warning_indent(
                    "Changing grid of point data may lead to unexpected results.\n"
                    + f"Consider supplying {name} in the output domain's CRS and grid."
                )
            return self.reproject_raster_to_dst(
                src,
                resampling=resampling,
                compatibility_resampling=compatibility_resampling,
                name=name,
            )

        # Reproject if domains are not aligned. This already includes the correct spatial
        # extensions.
        # Check if exact coordinate values of self.origin_x/y within src are integer
        origin_y_coord, origin_x_coord = riotf.rowcol(
            src.transform, self.origin_x, self.origin_y, op=lambda _: _
        )
        if not (
            abs(origin_x_coord - np.round(origin_x_coord)) < 0.2
            and abs(origin_y_coord - np.round(origin_y_coord)) < 0.2
        ):
            logger.info(
                f"Changing grid of {name} input due to non-aligned grids "
                + f"using {resampling.name} resampling."
            )
            logger.debug_indent(debug_message)
            if warning_point_data:
                logger.warning_indent(
                    "Changing grid of point data may lead to unexpected results.\n"
                    + f"Consider supplying {name} in the output domain's CRS and grid."
                )
            return self.reproject_raster_to_dst(src, resampling=resampling, name=name)

        # Cut if everything is aligned.
        logger.info(f"Cut {name} input to target grid.")
        return self.cut_aligned_raster_to_dst(src, name)

    def resolution_in_dst(self, src: rio.DatasetReader) -> float:
        """Calculate the resolution of the source data in units of the destination domain.

        If the CRSs of the source and the destination domain are the same, the resolution is
        calculated directly.

        If the CRSs are different, the resolution in the destination CRS is calculate by determining
        the approximate part of the source domain that covers the destination domain and calculating
        the average resolution in this part. First, the corners of the destination domain in the
        source CRS are calculated. The corners are likely not aligned with the source grid so an
        enclosing rectangle is calculated in the source CRS. The x/y resolution in the destination
        CRS is then the distance of the corners in the destination CRS divided by the respective
        number of pixels. The average of the two resolutions is calculated using the geometric mean
        to avoid bias towards the coarser resolutions.

        Args:
            src: Source raster.

        Returns:
            Geometric mean of the resolution in x and y direction.
        """
        # If the CRSs are the same, resolution is given by size divided by number of pixels.
        if src.crs == self.dst_crs:
            resolution_x = np.abs(src.bounds.right - src.bounds.left) / src.width
            resolution_y = np.abs(src.bounds.top - src.bounds.bottom) / src.height
            # Geometric mean of the resolutions in x and y direction
            return np.sqrt(resolution_x * resolution_y)

        # If the CRSs are different, calculate the resolution in the destination CRS.

        # The destination (dst) domain with its corners ul, ur, ll, lr in the source (src) dataset:
        # ---------------------------
        # |         1----ur---3=max |
        # | src     |   / \   |     |
        # |         |  /   \  |     |
        # |         | /     \ |     |
        # |         ul       \|     |
        # |         | \  dst  lr    |
        # |         |  \     /|     |
        # |         |   \   / |     |
        # |         |    \ /  |     |
        # |     min=0-----ll--2     |
        # ---------------------------

        # Calculate the values of the corners of the destination (dst) domain in the source (src)
        # CRS.
        x_corners_dst = [self.corner_ll[0], self.corner_ul[0], self.corner_lr[0], self.corner_ur[0]]
        y_corners_dst = [self.corner_ll[1], self.corner_ul[1], self.corner_lr[1], self.corner_ur[1]]
        x_corners_src, y_corners_src = self._transform_points(
            src_crs=self.dst_crs, dst_crs=src.crs, x=x_corners_dst, y=y_corners_dst
        )

        # Components of enclosing rectangle 0-3 of the destination domain corners in source CRS
        corner_min_src_x = np.min(x_corners_src)
        corner_min_src_y = np.min(y_corners_src)
        corner_max_src_x = np.max(x_corners_src)
        corner_max_src_y = np.max(y_corners_src)

        # Corresponding integer coordinates in the source raster (argument op defaults to floor).
        # Due to floor, the coordinates likely do not exactly correspond to the corner_min/max_src.
        corner_min_src_coord_y, corner_min_src_coord_x = riotf.rowcol(
            src.transform,
            corner_min_src_x,
            corner_min_src_y,
        )
        corner_max_src_coord_y, corner_max_src_coord_x = riotf.rowcol(
            src.transform,
            corner_max_src_x,
            corner_max_src_y,
        )
        # Check for ints to please static type checking.
        if not (
            isinstance(corner_min_src_coord_x, (int, np.integer))
            and isinstance(corner_min_src_coord_y, (int, np.integer))
            and isinstance(corner_max_src_coord_x, (int, np.integer))
            and isinstance(corner_max_src_coord_y, (int, np.integer))
        ):
            raise ValueError("Corner coordinates are not int.")

        # If the min/max coordinates are the same, add 1 to the max coordinate to get a non-zero
        # distance.
        if corner_max_src_coord_x == corner_min_src_coord_x:
            corner_max_src_coord_x += 1
        if corner_max_src_coord_y == corner_min_src_coord_y:
            corner_max_src_coord_y += 1

        # Update exact enclosing rectangle 0-3 in source CRS.
        corner_min_src_x, corner_min_src_y = riotf.xy(
            src.transform, corner_min_src_coord_y, corner_min_src_coord_x
        )
        corner_max_src_x, corner_max_src_y = riotf.xy(
            src.transform, corner_max_src_coord_y, corner_max_src_coord_x
        )

        # Calculate corners of the enclosing rectangle 0-3 in dst CRS.
        x_corners_src = [corner_min_src_x, corner_min_src_x, corner_max_src_x, corner_max_src_x]
        y_corners_src = [corner_min_src_y, corner_max_src_y, corner_min_src_y, corner_max_src_y]
        x_corners_dst, y_corners_dst = self._transform_points(
            src_crs=src.crs, dst_crs=self.dst_crs, x=x_corners_src, y=y_corners_src
        )

        # Resolution of src in dst is given by distance of the corners in dst CRS divided by number
        # of pixels. For the former, an average of the respective two sides of the enclosing
        # rectangle is used.
        resolution_x = np.mean(
            [
                np.sqrt(
                    (x_corners_dst[2] - x_corners_dst[0]) ** 2
                    + (y_corners_dst[2] - y_corners_dst[0]) ** 2
                ),
                np.sqrt(
                    (x_corners_dst[3] - x_corners_dst[1]) ** 2
                    + (y_corners_dst[3] - y_corners_dst[1]) ** 2
                ),
            ]
        ) / np.abs(corner_max_src_coord_x - corner_min_src_coord_x)
        resolution_y = np.mean(
            [
                np.sqrt(
                    (x_corners_dst[1] - x_corners_dst[0]) ** 2
                    + (y_corners_dst[1] - y_corners_dst[0]) ** 2
                ),
                np.sqrt(
                    (x_corners_dst[3] - x_corners_dst[2]) ** 2
                    + (y_corners_dst[3] - y_corners_dst[2]) ** 2
                ),
            ]
        ) / np.abs(corner_max_src_coord_y - corner_min_src_coord_y)

        # Geometric mean of the resolutions in x and y direction
        return np.sqrt(resolution_x * resolution_y)

    def reproject_raster_to_dst(
        self,
        src: rio.DatasetReader,
        resampling: riowp.Resampling = riowp.Resampling.nearest,
        compatibility_resampling: Optional[riowp.Resampling] = None,
        name: Optional[str] = None,
    ) -> ma.MaskedArray:
        """Reproject the input data to the destination domain.

        Read the data and reproject to the destination CRS with the given resampling method. If
        required, the result is masked where nearest neighbour interpolation would be masked. This
        might be necessary because the chosen resampling method might behave differently than
        nearest neighbour resampling and consistent behaviour might be required. If the debug option
        was enabled while initializing the class, the reprojected data will be written to a file in
        the debug folder.

        Args:
            src: Input raster.
            resampling: Resampling method. Defaults to riowp.Resampling.nearest.
            compatibility_resampling: Masked values of this resampling method should be applied to
              the output.
            name: Variable name for debug output. Defaults to None.

        Raises:
            NotImplementedError: Input data has wrong number of dimensions.

        Returns:
            Reprojected raster data in the destination domain.
        """
        src_values = src.read(masked=True)
        src_values = _correct_fill_value(src_values)

        # Create empty reprojection target.
        ndim = len(src_values.shape)
        if ndim == 2:
            dst_filled = np.empty((self.dst_height, self.dst_width), dtype=src_values.dtype)
        elif ndim == 3:
            dst_filled = np.empty(
                (src_values.shape[0], self.dst_height, self.dst_width), dtype=src_values.dtype
            )
        else:
            logger.critical_indent_raise(
                f"Cannot handle input file for {name} with {ndim} dimensions.",
                exception_type=NotImplementedError,
            )

        # Handle masked entries manually until https://github.com/rasterio/rasterio/issues/2575 is
        # fixed.
        src_filled = src_values.filled(fill_value=src_values.fill_value)
        riowp.reproject(
            source=src_filled,
            destination=dst_filled,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src_values.fill_value,
            dst_transform=self.dst_transform,
            dst_crs=self.dst_crs,
            dst_nodata=src_values.fill_value,
            resampling=resampling,
        )
        dst: ma.MaskedArray = ma.masked_array(
            dst_filled,
            mask=(dst_filled == src_values.fill_value),
            fill_value=src_values.fill_value,
            dtype=src_values.dtype,
        )

        # Mask dst where compatibility resampling would be masked for consistency.
        if compatibility_resampling is not None and resampling != compatibility_resampling:
            # Calculate results as if compatibility resampling would be used.
            dst_filled_comp = np.empty_like(dst_filled)
            riowp.reproject(
                source=src_filled,
                destination=dst_filled_comp,
                src_transform=src.transform,
                src_crs=src.crs,
                src_nodata=src_values.fill_value,
                dst_transform=self.dst_transform,
                dst_crs=self.dst_crs,
                dst_nodata=src_values.fill_value,
                resampling=compatibility_resampling,
            )
            dst_comp: ma.MaskedArray = ma.masked_array(
                dst_filled_comp,
                mask=(dst_filled_comp == src_values.fill_value),
                fill_value=src_values.fill_value,
                dtype=src_values.dtype,
            )

            # Mask like nearest neighbour.
            ma.masked_where(dst_comp.mask, dst, copy=False)

        # Write out the reprojected data if debug is enabled.
        if self.debug:
            if name is None:
                raise ValueError("name has to be given.")

            if self.domain_name is None:
                raise ValueError("domain_name is not set.")

            self.debug_file_prefix.parent.mkdir(exist_ok=True)
            output_reprojected = self.debug_file_prefix.with_stem(
                f"{self.debug_file_prefix.stem}_{name}-reprojected_{self.domain_name}"
            ).with_suffix(".tif")
            self.write_dst_geotiff(output_reprojected, dst)
            logger.debug_indent(f"Result written to {output_reprojected}.", hierarchy=1)

        return dst

    def cut_aligned_raster_to_dst(
        self,
        src: rio.DatasetReader,
        name: Optional[str] = None,
    ) -> ma.MaskedArray:
        """Cut input data to the destination domain assuming the same CRS and alignment.

        If the debug option was enabled while initializing the class, the reprojected data will be
        written to a file in the debug folder.

        Args:
            src: Input raster.
            name: Variable name for debug output. Defaults to None.

        Returns:
            Cut raster data in the destination domain.
        """
        if self.ignore_input_georeferencing:
            # Coordinates from top left corner of the domain relative to top left corner of the
            # input data.
            input_y1 = src.shape[0] - self.input_y0 - self.dst_height
            # Window covering the target domain.
            read_window = rio.windows.Window(
                self.input_x0,
                input_y1,
                self.dst_width,
                self.dst_height,
            )
        else:
            # Coordinates from top left corner of the domain relative to top left corner of the
            # input data. Use center here to avoid rounding errors and minor shifts between the
            # grids.
            corner_index_ul_y, corner_index_ul_x = riotf.rowcol(
                src.transform, self.corner_ul_cr[0], self.corner_ul_cr[1]
            )

            # Window covering the target domain.
            read_window = rio.windows.Window(
                corner_index_ul_x,
                corner_index_ul_y,
                self.dst_width,
                self.dst_height,
            )

        # Read only data within window.
        dst = src.read(masked=True, window=read_window)
        dst = _correct_fill_value(dst)

        if self.debug:
            if name is None:
                raise ValueError("name has to be given.")

            if self.domain_name is None:
                raise ValueError("domain_name is not set.")

            self.debug_file_prefix.parent.mkdir(exist_ok=True)
            output_cut = self.debug_file_prefix.with_stem(
                f"{self.debug_file_prefix.stem}_{name}-cut_{self.domain_name}"
            ).with_suffix(".tif")
            self.write_dst_geotiff(output_cut, dst)
            logger.debug_indent(f"Result written to {output_cut}.", hierarchy=1)

        return dst

    def write_dst_geotiff(self, file: Path, dst: ma.MaskedArray) -> None:
        """Write raster data of the size of the domain to a GeoTIFF file.

        Args:
            file: Output file.
            dst: Output data.

        Raises:
            ValueError: Output data has wrong shape or size.
        """
        if len(dst.shape) == 2:
            if self.dst_height != dst.shape[0] or self.dst_width != dst.shape[1]:
                raise ValueError("dst has wrong shape.")
            count = 1
        elif len(dst.shape) == 3:
            if self.dst_height != dst.shape[1] or self.dst_width != dst.shape[2]:
                raise ValueError("dst has wrong shape.")
            count = dst.shape[0]
        else:
            raise ValueError("dst has wrong shape.")

        with rio.open(
            file,
            "w",
            driver="GTiff",
            height=self.dst_height,
            width=self.dst_width,
            count=count,
            dtype=dst.dtype,
            crs=self.dst_crs,
            transform=self.dst_transform,
            nodata=dst.fill_value,
            compress="DEFLATE",
        ) as output_file:
            if len(dst.shape) == 2:
                output_file.write(dst, count)
            else:
                output_file.write(dst)


def _correct_fill_value(raster_values: ma.MaskedArray) -> ma.MaskedArray:
    """Correct the fill value of the raster data.

    If the input raster is of integer dtype, check if the fill value is within the dtype boundaries
    and extend dtype if necessary. This is necessary e.g. when reading a uint8 geotiff without
    fill value, the fill value is set to 999999.

    Args:
        raster_values: Input data.

    Returns:
        Raster data with corrected fill value.
    """
    if np.issubdtype(raster_values.dtype, np.integer):
        if (
            raster_values.fill_value > np.iinfo(raster_values.dtype).max
            or raster_values.fill_value < np.iinfo(raster_values.dtype).min
        ):
            # change dtype, ma.masked_array is needed for the correct type
            raster_values = ma.masked_array(
                raster_values.astype(np.result_type(raster_values.dtype, raster_values.fill_value))
            )

    return raster_values
