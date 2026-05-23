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

"""Main module of palm_csd to create a PALM static driver.

create_driver is the main routine, which calls the necessary functions to create a static driver in
this module.
"""

import logging
import math
from itertools import combinations
from os import PathLike
from typing import Dict, List, Optional, Union, cast

import numpy as np
import numpy.ma as ma
import pandas as pd
import yaml

from palm_csd import (
    StatusLogger,
    csd_domain,
    geo_converter,
    netcdf_data,
    statistics,
    tools,
    vegetation,
)
from palm_csd.constants import (
    INPUT_DATA_EXPANDED,
    NBUILDING_SURFACE_LAYER,
    VT_HIGH_VEGETATION,
    VT_NO_PLANTS,
    ColumnGeneral,
    IndexBuildingGeneralPars,
    IndexBuildingSurfaceLevel,
    IndexBuildingSurfaceType,
    IndexBuildingType,
    IndexVegetationPars,
    IndexVegetationType,
    IndexWaterPars,
    InputData,
)
from palm_csd.csd_config import (
    CSDConfig,
    CSDConfigInput,
    CSDConfigSettings,
    value_defaults,
)
from palm_csd.csd_domain import (
    CSDDomain,
)
from palm_csd.lcz import LCZTypes
from palm_csd.statistics import static_driver_statistics
from palm_csd.tools import (
    blend_array_2d,
    check_consistency_3,
    check_consistency_4,
    height_to_z_grid,
    interpolate_2d,
    is_missing,
    ma_isin,
)
from palm_csd.vegetation import DomainTree

# Module logger. In __init__.py, it is ensured that the logger is a StatusLogger. For type checking,
# do explicit cast.
logger = cast(StatusLogger, logging.getLogger(__name__))


def create_driver(
    input_configuration_file: Union[str, PathLike],
    verbose: Optional[Dict[str, bool]] = None,
    show_plot: bool = False,
    pdf: bool = False,
    png: bool = False,
) -> None:
    """Main routine for creating the static driver.

    Args:
        input_configuration_file: Input configuration YAML file.
        verbose: Dictionary of debug flags and if they are enabled. Defaults to None.
        show_plot: Show a plot of the result. Defaults to False.
        pdf: Save the plot of the static driver as PDF. Defaults to False.
        png: Save the plot of the static driver as PNG. Defaults to False.
    """
    # If verbose is None, set to empty dictionary to simplify further processing.
    if verbose is None:
        verbose = {}

    logger.status("Reading configuration.")

    # Load yml configuration file.
    try:
        with open(input_configuration_file, "r", encoding="utf-8") as file:
            input_configuration_dict = yaml.safe_load(file)
    except FileNotFoundError:
        logger.critical(f"Configuration file {input_configuration_file} not found.")
        raise

    # Read configuration file and set parameters accordingly.
    config = CSDConfig(input_configuration_dict)

    logger.status("Initializing domains.")

    def add_domain_and_parents(name: str, domains: Dict[str, CSDDomain]) -> None:
        """Recursively create domain and its parents as CSDDomain objects.

        They are added to the domains with parents first.

        Raises:
            ValueError: Parent domain of name set but not found.

        Args:
            name: Name of the domain to create and add.
            domains: Dictionary of domain names and CSDDomain objects.
        """
        if name not in domains:
            parent = None
            parent_name = config.domain_dict[name].domain_parent
            if parent_name is not None:
                if parent_name not in config.domain_dict:
                    logger.critical_raise(
                        f"Parent domain {parent_name} of domain {name} not found."
                    )
                add_domain_and_parents(parent_name, domains)
                parent = domains[parent_name]
            domains[name] = CSDDomain(
                name, config, parent, gis_debug_output=verbose.get("gis", False)
            )

    # Create domains and add them to the dictionary. Parents are added first.
    domains: Dict[str, CSDDomain] = {}
    for name in config.domain_dict:
        add_domain_and_parents(name, domains)

    # The root parent should be the first domain in the dictionary. Just to be sure, apply
    # find_root.
    domains_root = next(iter(domains.values())).find_root()

    # Check for overlap of domains. This is done recursively for all children of a domain starting
    # from the root parent.
    check_overlap(domains_root)

    # Initialize LCZ types and update the values from the configuration.
    lcz_types = LCZTypes(config.settings.season, config.lcz.height_geometric_mean)
    lcz_types.update_defaults(config.lcz)

    # Set debug output for the modules depending on the configuration.
    if verbose.get("gis", False):
        geo_converter.logger.setLevel(logging.DEBUG)
    if verbose.get("io", False):
        netcdf_data.logger.setLevel(logging.DEBUG)
        csd_domain.logger.setLevel(logging.DEBUG)
    if verbose.get("misc", False):
        logger.setLevel(logging.DEBUG)
        tools.logger.setLevel(logging.DEBUG)
        statistics.logger.setLevel(logging.DEBUG)
    if verbose.get("vegetation", False):
        vegetation.logger.setLevel(logging.DEBUG)

    # Find the minium of all terrain heights. This value will be subtracted from all domain's
    # terrain heights.
    zt_min = minimum_terrain_height(domains)

    # Loop over domains, domains are independent of each other except terrain height. Potentially,
    # the terrain height of the parent domain is needed. By traversing the domains tree and starting
    # from the root parent, it is ensured that parents are always dealt with before the children.
    for domain in domains_root.traverse():
        if domain.parent is not None:
            log_str_parent = f" WITH PARENT DOMAIN {domain.parent.name}"
        else:
            log_str_parent = ""
        logger.status(f"WORKING ON DOMAIN {domain.name}" + log_str_parent + ".")

        domain.remove_existing_output()

        if domain.config.lcz_input == "full":
            process_coordinates(domain, zt_min)

            process_lcz(domain, lcz_types)

            if (
                domain.config.water_temperature is not None
                or domain.input_config.file_water_temperature is not None
            ):
                process_water_temperature(domain)

            domain.write_global_attributes()

        else:  # standard case
            process_coordinates(domain, zt_min)

            process_buildings_bridges(domain, config.settings)

            process_types(domain)
            process_street_type_crossing(domain)

            process_resolved_vegetation(domain, config.settings)

            if (
                domain.config.water_temperature is not None
                or domain.input_config.file_water_temperature is not None
            ):
                process_water_temperature(domain)

            consistency_check_update_surface_fraction(domain)
            domain.write_global_attributes()

    print_unused_input_files(config.input_dict)

    # Calculate statistics from the output netcdf files.
    for domain in domains.values():
        logger.status(f"STATISTICS OF DOMAIN {domain.name}.")
        # Generate plot file if requested.
        if pdf:
            plot_file = domain.file_output.with_suffix(".pdf")
        elif png:
            plot_file = domain.file_output.with_suffix(".png")
        else:
            plot_file = None

        static_driver_statistics(
            domain.file_output,
            show_plot=show_plot,
            plot_file=plot_file,
            plot_title=f"Domain: {domain.name}",
        )


def check_overlap(domain: CSDDomain) -> None:
    """Check if children of a domain overlap and check recursively their children.

    It is assumed that a child is fully covered by its parent. Then, only children of the same
    parent could possibly overlap when the parent does not overlap with one of its sibling. This
    function is called recursively for all children.

    Args:
        domain: Parent domain to start the check.
    """
    children = domain.get_children()
    for child1, child2 in combinations(children, 2):
        if child1.overlaps(child2):
            logger.warning(
                f"Domains {child1.name} and {child2.name} overlap. "
                + "Only one-way nesting is allowed."
            )

    for child in children:
        check_overlap(child)


def minimum_terrain_height(domains: Dict[str, CSDDomain]) -> float:
    """Calculate minimum terrain height of given domains.

    Args:
        domains: All domains to consider.

    Returns:
        Minimum terrain height of all domains.
    """
    logger.status("Calculating minimum terrain height of all domains.")

    zt_min = math.inf
    for domain in domains.values():
        zt = domain.read(InputData.zt)
        zt_min = min(zt_min, min(zt.flatten()))

    logger.info(f"Shifting down all domains by minimum terrain height of {zt_min:0.2f} m.")
    return zt_min


def process_coordinates(domain: CSDDomain, zt_min: float) -> None:
    """Process coordinates and terrain height of a domain.

    When the geo_converter is defined in the domain, it is used to calculate the coordinates of the
    domain. Otherwise, the coordinates are read from input data. The coordinates are written to the
    result file. z_min is subtracted from the terrain height. If interpolate_terrain set to True,
    the terrain height is adopted to the parent's terrain height. The terrain height is kept in
    memory and written to the result file.

    Args:
        domain: Domain to process.
        zt_min: Height to subtract from the terrain height.

    Raises:
        ValueError: x0 and y0 are not defined when geo_converter is not defined.
        ValueError: If interpolate_terrain set to True, but the parent domain is not defined.
        ValueError: If interpolate_terrain set to True, but the parent domain's x and y coordinates
          are not set.
        ValueError: If interpolate_terrain set to True, but the parent domain's terrain height is
          not set.
        ValueError: If interpolate_terrain set to True, but the parent does not fully cover the
          child domain.
    """
    logger.status("Processing coordinates.")

    # Use origin_x and origin_y to calculate UTM and lon/lat coordinates
    if domain.geo_converter is not None:
        domain.origin_x = domain.geo_converter.origin_x
        domain.origin_y = domain.geo_converter.origin_y
        domain.origin_lon = domain.geo_converter.origin_lon
        domain.origin_lat = domain.geo_converter.origin_lat

        # Global x and y coordinates (cell centre) relative to root parent domain
        # Used only for zt interpolation below
        x_global, y_global = domain.geo_converter.global_palm_coordinates()
        domain.x_global.values = ma.MaskedArray(x_global)
        domain.y_global.values = ma.MaskedArray(y_global)

        # Coordinates
        e_utm, n_utm, lon, lat = domain.geo_converter.geographic_coordinates()

        # Write CRS
        domain.write_crs_to_file()

    else:
        # Get coordinates near origin
        if domain.x0 is None or domain.y0 is None:
            raise ValueError(f"Domain {domain.name} has no x0 or y0 defined")

        x_utm_origin = domain.read_nc_2d(
            domain.input_config.files["x_utm"][0],
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        domain.input_config.add_used_file(domain.input_config.files["x_utm"][0])

        y_utm_origin = domain.read_nc_2d(
            domain.input_config.files["y_utm"][0],
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        domain.input_config.add_used_file(domain.input_config.files["y_utm"][0])

        lat_origin = domain.read_nc_2d(
            domain.input_config.files["lat"][0],
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        domain.input_config.add_used_file(domain.input_config.files["lat"][0])

        lon_origin = domain.read_nc_2d(
            domain.input_config.files["lon"][0],
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        domain.input_config.add_used_file(domain.input_config.files["lon"][0])

        # Calculate position of origin. Added as global attributes later
        domain.origin_x = float(x_utm_origin[0, 0]) - 0.5 * (
            float(x_utm_origin[0, 1]) - float(x_utm_origin[0, 0])
        )
        domain.origin_y = float(y_utm_origin[0, 0]) - 0.5 * (
            float(y_utm_origin[1, 0]) - float(y_utm_origin[0, 0])
        )
        domain.origin_lon = float(lon_origin[0, 0]) - 0.5 * (
            float(lon_origin[0, 1]) - float(lon_origin[0, 0])
        )
        domain.origin_lat = float(lat_origin[0, 0]) - 0.5 * (
            float(lat_origin[1, 0]) - float(lat_origin[0, 0])
        )

        # Read x and y values
        domain.x_global.values = domain.read_nc_1d(domain.input_config.files["x_utm"][0], "x")
        domain.y_global.values = domain.read_nc_1d(
            domain.input_config.files["y_utm"][0], "y", x0=domain.y0, x1=domain.y1
        )

        # Read and write lon, lat and UTM coordinates
        lat = domain.read_nc_2d(domain.input_config.files["lat"][0])
        lon = domain.read_nc_2d(domain.input_config.files["lon"][0])

        e_utm = domain.read_nc_2d(domain.input_config.files["x_utm"][0])
        n_utm = domain.read_nc_2d(domain.input_config.files["y_utm"][0])

        # Write CRS
        crs = domain.read_nc_crs()
        crs.to_nc()

    # Shift x and y coordinates for x and y local cell centre coordinates of domain
    # Used as output dimensions
    domain.x.values = (
        domain.x_global.values
        - min(domain.x_global.values.flatten())
        + domain.config.pixel_size / 2.0
    )
    domain.y.values = (
        domain.y_global.values
        - min(domain.y_global.values.flatten())
        + domain.config.pixel_size / 2.0
    )

    domain.lat.to_nc(lat)
    domain.lon.to_nc(lon)

    domain.E_UTM.to_nc(e_utm)
    domain.N_UTM.to_nc(n_utm)

    # Read and process terrain height (zt). Its values are stored in the domain object to be
    # available for potential child domain.
    domain.zt.values = domain.read(InputData.zt)
    domain.zt.values = domain.zt.values - zt_min
    domain.origin_z = float(zt_min)

    # If necessary, interpolate parent domain terrain height on child domain grid and blend
    # the two.
    if domain.config.interpolate_terrain:
        if domain.parent is None:
            raise ValueError("Interpolation of terrain height requires a parent domain")
        if domain.parent.x_global.values is None or domain.parent.y_global.values is None:
            raise ValueError(f"x_utm or y_utm of parent {domain.parent.name} not calculated")
        if domain.parent.zt.values is None:
            raise ValueError(f"zt of parent {domain.parent.name} not calculated")

        tmp_x0 = np.searchsorted(domain.parent.x_global.values, domain.x_global.values[0]) - 1
        tmp_y0 = np.searchsorted(domain.parent.y_global.values, domain.y_global.values[0]) - 1
        tmp_x1 = np.searchsorted(domain.parent.x_global.values, domain.x_global.values[-1]) + 1
        tmp_y1 = np.searchsorted(domain.parent.y_global.values, domain.y_global.values[-1]) + 1

        if tmp_x0 < 0:
            raise ValueError(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the left border"
            )
        if tmp_y0 < 0:
            raise ValueError(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the bottom border"
            )
        if tmp_x1 > domain.parent.x_global.values.shape[0]:
            raise ValueError(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the right border"
            )
        if tmp_y1 > domain.parent.y_global.values.shape[0]:
            raise ValueError(
                f"Parent {domain.parent.name} not fully covering "
                + f"child {domain.name} on the top border"
            )

        tmp_x = domain.parent.x_global.values[tmp_x0:tmp_x1]
        tmp_y = domain.parent.y_global.values[tmp_y0:tmp_y1]

        zt_parent = domain.parent.zt.values[tmp_y0:tmp_y1, tmp_x0:tmp_x1]

        # Interpolate array and bring to PALM grid of child domain.
        zt_ip = interpolate_2d(
            zt_parent, tmp_x, tmp_y, domain.x_global.values, domain.y_global.values
        )
        zt_ip = height_to_z_grid(zt_ip, domain.parent.config.dz)

        # Shift the child terrain height according to the parent mean terrain height.
        # mypy wants us to check again if zt.values is None, not sure why. Let's do it.
        if domain.zt.values is None:
            raise ValueError(f"Domain {domain.name} has undefined zt values")
        z_mean = np.mean(domain.zt.values)
        z_mean_parent = np.mean(zt_ip)
        logger.debug(f"Average domain height: {z_mean:0.2f} m.")
        logger.debug(f"Avergage covered parent domain height: {z_mean_parent:0.2f} m.")
        dz_mean = z_mean - z_mean_parent
        logger.info(
            f"Shifting down terrain height by {dz_mean:0.2f} m to adjust for parent domain height."
        )
        domain.zt.values = domain.zt.values - dz_mean
        if domain.zt.values is None:
            raise ValueError(f"Domain {domain.name} has undefined zt values")

        # Blend over the parent and child terrain height within a radius of 50 px (or less if
        # domain is smaller than 50 px).
        domain.zt.values = ma.MaskedArray(
            blend_array_2d(domain.zt.values, zt_ip, min(50, min(domain.zt.values.shape) * 0.5))
        )

    # If necessary, bring terrain height to PALM's vertical grid. This is either forced by
    # the user or implicitly by using interpolation for a child domain.
    if domain.zt.values is None:
        raise ValueError(f"Domain {domain.name} has undefined zt values")
    if domain.config.use_palm_z_axis:
        domain.zt.values = ma.MaskedArray(height_to_z_grid(domain.zt.values, domain.config.dz))

    domain.zt.to_nc()


def process_buildings_bridges(domain: CSDDomain, settings: CSDConfigSettings) -> None:
    """Process buildings and bridges of a domain.

    The building height, id and type is read from the input data and checked for consistency.
    Optionally, the 3d building field is calculated. All data is written to the result file. Read
    bridge height and id from the input data and check for consistency. If bridges are present and
    thick enough to be represented by the z grid, they are added to the 3d building field;
    buildings2d is not adjusted. The building type is set to the Bridge type. If both, buildings and
    bridges are defined at a pixel, the building information is chosen. All data is written to
    the result file. The building parameters set in the configuration are written to the respective
    variables, where the building pixels are set.

    Args:
        domain: Domain to process.
        settings: General settings from the configuration.

    Raises:
        ValueError: Building IDs are missing for some building pixels.
        ValueError: Bridge IDs are missing for some bridge pixels.
        ValueError: Length of building parameter data does not match number of building surface
          layers.
        ValueError: Invalid value for a key in the building parameter data.
        ValueError: Invalid dimension for the building parameter output variable.
    """
    logger.status("Processing buildings and bridges.")

    buildings = domain.read_buildings()

    if buildings is not None:
        buildings_2d, building_id, building_type = buildings
    else:
        buildings_2d = domain.read(InputData.buildings_2d)
        building_id = domain.read(InputData.building_id)
        building_type = domain.read(InputData.building_type)

    # Remove buildings in border area if requested and store where buildings were removed.
    if domain.config.building_free_border_width > 0.0:
        buildings_2d_original_mask = ma.getmaskarray(buildings_2d).copy()
        n_border = int(np.ceil(domain.config.building_free_border_width / domain.config.pixel_size))
        logger.info(
            f"Applying building free border of {domain.config.building_free_border_width:0.2f} m "
            + f"({n_border} pixels)."
        )
        buildings_2d[:n_border, :] = ma.masked
        buildings_2d[-n_border:, :] = ma.masked
        buildings_2d[:, :n_border] = ma.masked
        buildings_2d[:, -n_border:] = ma.masked
        domain.buildings_2d_removed = ma.getmaskarray(buildings_2d) != buildings_2d_original_mask

    # Check if there is a building_id (no default value applied) for all buildings_2d pixels.
    building_without_id = ma.getmaskarray(building_id)[~ma.getmaskarray(buildings_2d)]
    logger.critical_argwhere_raise(
        "Building ID missing for",
        building_without_id,
        "building pixels defined by buildings_2d.",
    )

    if buildings_2d.mask.all():
        logger.info("No buildings in domain.")

    buildings_2d_small = buildings_2d < 0.5 * domain.config.dz
    logger.warning_argwhere(
        "Found",
        buildings_2d_small,
        "building pixels with height < 1/2 dz.\n"
        + "They will be treated by PALM as a flat surface.",
    )

    # Apply building mask to building_id and building_type.
    building_id.mask = buildings_2d.mask.copy()
    building_type.mask = buildings_2d.mask.copy()

    # Express bridge depth in terms of grid height to minimize discretization error.
    bridge_depth_grid = round(domain.config.bridge_depth / domain.config.dz) * domain.config.dz
    if bridge_depth_grid == 0:
        bridges_2d = ma.masked_all_like(buildings_2d)
        logger.warning("Bridge depth < 1/2 dz. Bridges will not be added.")
    else:
        bridges = domain.read_bridges()
        if bridges is not None:
            bridges_2d, bridges_id = bridges
        else:
            bridges_2d = domain.read(InputData.bridges_2d)
            bridges_id = domain.read(InputData.bridges_id)

        # Check if there is a bridges_id (no default value applied) for all bridges_2d pixels.
        bridge_without_id = ma.getmaskarray(bridges_id)[~ma.getmaskarray(bridges_2d)]
        logger.critical_argwhere_raise(
            "Bridge ID missing for", bridge_without_id, "bridge pixels defined by bridges_2d."
        )

        if bridges_2d.mask.all():
            logger.info("No bridges in domain.")

        bridges_2d_small = bridges_2d < 0.5 * domain.config.dz
        logger.warning_argwhere(
            "Found",
            bridges_2d_small,
            "bridge pixels with height < 1/2 dz.\n"
            + "They will be treated by PALM as a flat surface.",
        )

        buildings_bridges_overlap = ~ma.getmaskarray(buildings_2d) & ~ma.getmaskarray(bridges_2d)
        logger.warning_argwhere(
            "Buildings and bridges are overlapping at",
            buildings_bridges_overlap,
            "pixels.\n" + "Prefering building information at these pixels.",
        )

        bridges_id.mask = bridges_2d.mask.copy()
        building_id = ma.where(buildings_2d.mask & ~bridges_2d.mask, bridges_id, building_id)
        building_type = ma.where(
            buildings_2d.mask & ~bridges_2d.mask, IndexBuildingType.bridges, building_type
        )

    domain.buildings_2d.to_nc(buildings_2d)
    domain.building_id.to_nc(building_id)
    domain.building_type.to_nc(building_type)

    # Create 3d buildings if necessary. Add bridge pixels to building layer.
    if domain.config.generate_buildings_3d or (bridge_depth_grid > 0 and not bridges_2d.mask.all()):
        if not domain.config.generate_buildings_3d:
            logger.info("Creating 3D buildings due to the presence of bridges.")

        # Calculate maximum height of buildings and bridges, 0 if no buildings and bridges present.
        # Fill masked values with 0 before taking the maximum to avoid UserWarning about converting
        # a masked element to nan.
        z_max = np.max((ma.filled(buildings_2d.max(), 0.0), ma.filled(bridges_2d.max(), 0.0)))

        # z array for 3D buildings, z[1:] are at the centre of the grid cells
        z = np.arange(0, ma.ceil(z_max / domain.config.dz) + 1) * domain.config.dz
        z[1:] = z[1:] - domain.config.dz * 0.5

        # 3D buildings from buildings_2d
        # cell centre heights -  -  -
        # cell border heights -------
        # buildings_2d assigned to cell around z_k ////////
        #
        # -  -  -  z_k+1
        # ////////
        # -------- zw_k
        # ////////
        # -  -  -  z_k
        #
        # -------- zw_k-1
        #
        # discretization error on average 0:
        # z_k <= buildings_2d <= zw_k: overestimation of building height by up to 0.5 dz
        # zw_k <= buildings_2d < z_k+1: underestimation of building height by up to 0.5 dz
        #
        # Check mask to avoid masked values as result of ma.where.
        buildings_3d = ma.where(
            ~ma.getmaskarray(buildings_2d)[np.newaxis, :, :]
            & (z[:, np.newaxis, np.newaxis] <= buildings_2d.data[np.newaxis, :, :]),
            1,
            0,
        )

        # Add bridges to building layer.
        # Check mask to avoid masked values as result of ma.where.
        buildings_3d = ma.where(
            ~ma.getmaskarray(bridges_2d)[np.newaxis, :, :]
            & (z[:, np.newaxis, np.newaxis] > bridges_2d.data[np.newaxis, :, :] - bridge_depth_grid)
            & (z[:, np.newaxis, np.newaxis] <= bridges_2d.data[np.newaxis, :, :]),
            1,
            buildings_3d,
        )

        domain.z.values = z
        domain.buildings_3d.to_nc(buildings_3d)

    def normalize_building_fractions(fractions: np.ma.MaskedArray) -> None:
        """Normalize building surface fractions for wall, window and green of each surface level."""
        for surface_level in ["gfl", "agfl", "roof"]:
            mask_wall = ma.getmaskarray(fractions)[
                IndexBuildingSurfaceType[f"wall_{surface_level}"], :, :
            ]
            mask_windows = ma.getmaskarray(fractions)[
                IndexBuildingSurfaceType[f"window_{surface_level}"], :, :
            ]
            mask_green = ma.getmaskarray(fractions)[
                IndexBuildingSurfaceType[f"green_{surface_level}"], :, :
            ]
            n_undefined = mask_wall.astype(int) + mask_windows.astype(int) + mask_green.astype(int)
            to_zero = (n_undefined == 1) | (n_undefined == 2)
            # Warn when setting undefined fractions to 0.
            logger.warning_argwhere(
                f"Setting undefined {surface_level} building fractions to 0 at",
                to_zero,
                "pixels.",
            )
            fractions[IndexBuildingSurfaceType[f"wall_{surface_level}"], :, :] = ma.where(
                to_zero & mask_wall,
                0.0,
                fractions[IndexBuildingSurfaceType[f"wall_{surface_level}"], :, :],
            )
            fractions[IndexBuildingSurfaceType[f"window_{surface_level}"], :, :] = ma.where(
                to_zero & mask_windows,
                0.0,
                fractions[IndexBuildingSurfaceType[f"window_{surface_level}"], :, :],
            )
            fractions[IndexBuildingSurfaceType[f"green_{surface_level}"], :, :] = ma.where(
                to_zero & mask_green,
                0.0,
                fractions[IndexBuildingSurfaceType[f"green_{surface_level}"], :, :],
            )

            norm = (
                fractions[IndexBuildingSurfaceType[f"wall_{surface_level}"], :, :]
                + fractions[IndexBuildingSurfaceType[f"window_{surface_level}"], :, :]
                + fractions[IndexBuildingSurfaceType[f"green_{surface_level}"], :, :]
            )
            # Warn if norm is zero and set to masked.
            norm_zero = norm == 0.0
            logger.warning_argwhere(
                f"Removing {surface_level} building fraction at",
                norm_zero,
                "pixels with sum of fractions being zero.",
            )
            ma.masked_where(
                norm_zero,
                fractions[IndexBuildingSurfaceType[f"wall_{surface_level}"], :, :],
                copy=False,
            )
            ma.masked_where(
                norm_zero,
                fractions[IndexBuildingSurfaceType[f"window_{surface_level}"], :, :],
                copy=False,
            )
            ma.masked_where(
                norm_zero,
                fractions[IndexBuildingSurfaceType[f"green_{surface_level}"], :, :],
                copy=False,
            )
            ma.masked_where(norm_zero, norm, copy=False)

            # Warn if deviation from 1 is larger than 0.01. Do normalization anyway.
            to_normalize = ma.abs(norm - 1.0) > 0.01
            logger.warning_argwhere(
                f"Normalizing {surface_level} building fractions at",
                to_normalize,
                "pixels.",
            )
            fractions[IndexBuildingSurfaceType[f"wall_{surface_level}"], :, :] /= norm
            fractions[IndexBuildingSurfaceType[f"window_{surface_level}"], :, :] /= norm
            fractions[IndexBuildingSurfaceType[f"green_{surface_level}"], :, :] /= norm

    def add_building_parslike(variable_name: InputData) -> None:
        """Add global and local input data to building pars like where buildings are present.

        Global and local input data as given by raster and vector input are is written to the result
        file. Values are only applied to buildings. The variable_name gives the general input group,
        all corresponding input variables with specific levels and layers are processed.

        Args:
            variable_name: Variable group.

        Raises:
            ValueError: Length of input data does not match number of building surface layers.
            ValueError: Invalid value for a key in the input data.
            ValueError: Invalid dimension for the output variable.
        """
        # Global input data, which will be applied to all building pixels
        input_global = getattr(domain.config, variable_name)
        # Variables to read for local input data
        input_local_variable = [
            var for var in domain.input_config.files.keys() if var.startswith(variable_name)
        ] + [
            var
            for var in domain.input_config.columns.values()
            if isinstance(var, str) and var.startswith(variable_name)
        ]

        # Output variable and field to fill
        output_variable = getattr(domain, variable_name)
        output_field = output_variable.empty_array()
        if output_field.ndim not in (3, 4):
            raise ValueError(f"Invalid dimension for {output_variable.name}")

        # First, apply global input data everywhere.
        if input_global is not None:
            # with surface layers
            if output_field.ndim == 4:
                for key, value in input_global.items():
                    values: Union[List[float], List[int]]
                    if isinstance(value, int):
                        values = [value] * NBUILDING_SURFACE_LAYER
                    elif isinstance(value, float):
                        values = [value] * NBUILDING_SURFACE_LAYER
                    elif isinstance(value, List):
                        if len(value) != NBUILDING_SURFACE_LAYER:
                            raise ValueError(
                                f"Length of input data for {key} does not match "
                                + f"number of layers ({NBUILDING_SURFACE_LAYER})"
                            )
                        values = value
                    else:
                        raise ValueError(f"Invalid value for {key}")
                    for i, v in enumerate(values):
                        output_field[key, i, :, :] = v
            # without surface layers
            else:
                for key, value in input_global.items():
                    if not isinstance(value, (int, float)):
                        raise ValueError(f"Invalid value for {key}")
                    output_field[key, :, :] = value

        # Second, apply local input data.
        if variable_name == InputData.building_lai and settings.use_lai_for_roofs:
            # Use general lai for roofs.
            lai = domain.read(InputData.lai)
            output_field[IndexBuildingSurfaceLevel.roof, :, :] = ma.where(
                ma.getmaskarray(output_field)[IndexBuildingSurfaceLevel.roof, :, :],
                lai,
                output_field[IndexBuildingSurfaceLevel.roof, :, :],
            )
        for iv in input_local_variable:
            input_local = domain.read(iv)
            if input_local.mask.all():
                continue
            # Find all variables that input_local corresponds to.
            rows = INPUT_DATA_EXPANDED[INPUT_DATA_EXPANDED["name"].str.startswith(iv)]
            if rows.empty:
                raise ValueError(f"Input variable {iv} not found in expanded variable table")
            # Apply input_local to all corresponding output variables.
            for row in rows.itertuples(index=False):
                if output_field.ndim == 4:
                    if not isinstance(row.layer, (int, np.integer)):
                        raise ValueError(f"Layer value for {iv} is not integer-like: {row.layer!r}")
                    output_field[row.level, row.layer - 1, :, :] = ma.where(
                        input_local.mask,
                        output_field[row.level, row.layer - 1, :, :],
                        input_local,
                    )
                else:
                    output_field[row.level, :, :] = ma.where(
                        input_local.mask, output_field[row.level, :, :], input_local
                    )

        if output_field.mask.all():
            return

        # Finally, mask output field where no buildings are present and write to file.
        if output_field.ndim == 4:
            ma.masked_where(
                np.tile(
                    ma.getmaskarray(buildings_2d),
                    (output_field.shape[0], output_field.shape[1], 1, 1),
                ),
                output_field,
                copy=False,
            )
        else:
            ma.masked_where(
                np.tile(ma.getmaskarray(buildings_2d), (output_field.shape[0], 1, 1)),
                output_field,
                copy=False,
            )
            # When one fraction at a surface level is defined, set the missing of that triplet
            # (wall, window, green) fractions to 0. Normalize all building fractions to 1.
            if variable_name == InputData.building_fraction:
                normalize_building_fractions(output_field)
            # Mask LAI on roofs where green fraction on roofs is masked or zero.
            if variable_name == InputData.building_lai:
                fractions = domain.building_fraction.from_nc(allow_nonexistent=True)
                ma.masked_where(
                    fractions[IndexBuildingSurfaceType.green_roof, :, :].filled(0.0) == 0.0,
                    output_field[IndexBuildingSurfaceLevel.roof, :, :],
                    copy=False,
                )
            # Mask green type of roofs where green fraction on roofs is masked or zero.
            if variable_name == InputData.building_general_pars:
                fractions = domain.building_fraction.from_nc(allow_nonexistent=True)
                ma.masked_where(
                    fractions[IndexBuildingSurfaceType.green_roof, :, :].filled(0.0) == 0.0,
                    output_field[IndexBuildingGeneralPars.green_type_roof, :, :],
                    copy=False,
                )

        if output_field.mask.all():
            return

        output_variable.to_nc(output_field)

    add_building_parslike(InputData.building_albedo_type)
    add_building_parslike(InputData.building_emissivity)
    add_building_parslike(InputData.building_fraction)
    add_building_parslike(InputData.building_general_pars)  # building_fraction before
    add_building_parslike(InputData.building_heat_capacity)
    add_building_parslike(InputData.building_heat_conductivity)
    add_building_parslike(InputData.building_indoor_pars)
    add_building_parslike(InputData.building_lai)  # building_fraction before
    add_building_parslike(InputData.building_roughness_length)
    add_building_parslike(InputData.building_roughness_length_qh)
    add_building_parslike(InputData.building_thickness)
    add_building_parslike(InputData.building_transmissivity)


def process_lcz(
    domain: CSDDomain,
    lcz_types: LCZTypes,
) -> None:
    """Process LCZ data of a domain.

    Read LCZ data and derive surface_fraction, vegetation_type, LAI, pavement_type and water_type
    from it. soil_type is also read. If DCEP fields should be also calculated, fr_urb, fr_urbcl,
    fr_streetdir, street_width, building_width and building_height are also derived. All data is
    written to the result file.

    Args:
        domain: Domain to process.
        lcz_types: LCZ types.
    """
    lcz_type = domain.read(InputData.lcz, lcz_types=lcz_types)
    # LAI data input
    lai = domain.read(InputData.lai)
    # LAI from LCZ table
    lai_lcz = lcz_types.lai_from_lcz_map(lcz_type)
    # Use LAI from LCZ table if LAI data is missing
    lai = ma.where(lai.mask, lai_lcz, lai)

    soil_type = domain.read(InputData.soil_type)

    water_type = lcz_types.water_type_from_lcz_map(lcz_type)
    vegetation_type = lcz_types.vegetation_type_from_lcz_map(lcz_type)
    pavement_type = ma.masked_all_like(water_type)

    lai.mask = ma.mask_or(vegetation_type.mask, lai.mask)

    vegetation_pars = domain.vegetation_pars.empty_array()
    vegetation_pars[IndexVegetationPars.lai, :, :] = lai

    # Create surface_fraction array.
    domain.nsurface_fraction.values = ma.arange(0, 3)
    surface_fraction = ma.ones((domain.nsurface_fraction.size, domain.y.size, domain.x.size))

    # Remove soil_type for pixels with no vegetation_type and no pavement_type.
    soil_type = ma.where(vegetation_type.mask & pavement_type.mask, ma.masked, soil_type)

    surface_fraction[0, :, :] = ma.where(vegetation_type.mask, 0.0, 1.0)
    surface_fraction[1, :, :] = ma.where(pavement_type.mask, 0.0, 1.0)
    surface_fraction[2, :, :] = ma.where(water_type.mask, 0.0, 1.0)

    domain.surface_fraction.to_nc(surface_fraction)
    domain.vegetation_type.to_nc(vegetation_type)
    domain.vegetation_pars.to_nc(vegetation_pars)
    domain.pavement_type.to_nc(pavement_type)
    domain.water_type.to_nc(water_type)
    domain.soil_type.to_nc(soil_type)

    if domain.config.dcep:
        urban_fraction = lcz_types.urban_fraction_from_lcz_map(lcz_type)
        urban_class = lcz_types.urban_class_fraction_from_lcz_map(lcz_type)
        street_direction_fraction = lcz_types.street_direction_fraction_from_lcz_map(
            lcz_type, domain.config.udir
        )
        street_width = lcz_types.street_width_from_lcz_map(lcz_type, domain.config.udir)
        building_width = lcz_types.building_width_from_lcz_map(lcz_type, domain.config.udir)
        building_height = lcz_types.building_height_from_lcz_map(
            lcz_type, domain.config.z_uhl, domain.config.udir
        )

        domain.nuc.values = ma.arange(0, 1)
        domain.streetdir.values = ma.masked_array(domain.config.udir)
        domain.z_uhl.values = ma.masked_array(domain.config.z_uhl)

        domain.fr_urb.to_nc(urban_fraction)
        domain.fr_urbcl.to_nc(urban_class)
        domain.fr_streetdir.to_nc(street_direction_fraction)
        domain.street_width.to_nc(street_width)
        domain.building_width.to_nc(building_width)
        domain.building_height.to_nc(building_height)


def process_types(domain: CSDDomain) -> None:
    """Process vegetation, water, pavement and soil types of a domain.

    Read vegetation type, water_type, pavement_type and soil_type and make fields consistent. All
    data is written to the result file.

    Args:
        domain: Domain to process.

    Raises:
        ValueError: Several surface types defined for one pixel.
        ValueError: No surface types defined for some pixels and not replace_invalid_input_values
          set to True.
        ValueError: No default vegetation type defined when replacing invalid input values.
    """
    logger.status("Processing surface types.")

    vegetation_type = domain.read(InputData.vegetation_type)
    pavement_type = domain.read(InputData.pavement_type)
    water_type = domain.read(InputData.water_type)
    soil_type = domain.read(InputData.soil_type)
    # Use buildings_2d because it does not include bridges unlike building type
    building_height = domain.buildings_2d.from_nc()

    # Make arrays consistent
    # Set vegetation type to masked for pixel where a pavement type is set.
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~pavement_type.mask)

    # Set vegetation type to masked for pixel where a building type is set or originally buildings
    # set.
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~building_height.mask)
    if domain.buildings_2d_removed is not None:
        vegetation_type.mask = ma.mask_or(vegetation_type.mask, domain.buildings_2d_removed)

    # Set vegetation type to masked for pixel where a water type is set.
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~water_type.mask)

    # Remove pavement for pixels with buildings.
    pavement_type.mask = ma.mask_or(pavement_type.mask, ~building_height.mask)

    # Set pavement where buildings were removed in border area.
    if domain.buildings_2d_removed is not None:
        pavement_type = ma.where(
            domain.buildings_2d_removed,
            domain.config.building_free_border_pavement_type,
            pavement_type,
        )

    # Remove pavement for pixels with water.
    pavement_type.mask = ma.mask_or(pavement_type.mask, ~water_type.mask)

    # Remove water for pixels with buildings and originally buildings.
    water_type.mask = ma.mask_or(water_type.mask, ~building_height.mask)
    if domain.buildings_2d_removed is not None:
        water_type.mask = ma.mask_or(water_type.mask, domain.buildings_2d_removed)

    # Check for consistency and fill empty fields with default vegetation type.
    # number of not masked types per pixel
    n_type = np.count_nonzero(
        [
            ~ma.getmaskarray(vegetation_type),
            ~ma.getmaskarray(building_height),
            ~ma.getmaskarray(pavement_type),
            ~ma.getmaskarray(water_type),
        ],
        axis=0,
    )
    if (n_type == 1).all():
        logger.debug("Surface types are consistent for all pixels.")
    else:
        n_type_overdefined = n_type > 1
        if n_type_overdefined.any():
            type_overdefined = list(map(tuple, np.argwhere(n_type_overdefined)))
            logger.critical("Multiple surface types defined for some pixels.")
            logger.debug("Coordinates of overdefined pixels:")
            logger.debug(", ".join(map(str, type_overdefined)))
            raise ValueError("Inconsistent surface types.")

        # Only pixels with n_type==0 left, define vegetation here.
        n_type_undefined = n_type == 0
        type_undefined = list(map(tuple, np.argwhere(n_type_undefined)))
        if not domain.replace_invalid_input_values:
            logger.critical_raise(
                "No surface types defined for some pixels. "
                + "Enable replace_invalid_input_values for automatic replacement.",
            )
        if value_defaults["vegetation_type"].default is None:
            raise ValueError("No default vegetation type defined.")
        logger.warning(
            f"Setting default vegetation_type {value_defaults['vegetation_type'].default} "
            + f"for {np.sum(n_type_undefined)} pixels without surface type."
        )
        logger.debug_indent("Coordinates of undefined surface types:")
        logger.debug_indent(", ".join(map(str, type_undefined)))
        vegetation_type = ma.where(
            n_type_undefined, value_defaults["vegetation_type"].default, vegetation_type
        )

    # Remove soil_type for pixels without vegetation_type and pavement_type.
    soil_type = ma.where(vegetation_type.mask & pavement_type.mask, ma.masked, soil_type)

    # Create surface_fraction array.
    domain.nsurface_fraction.values = ma.arange(0, 3)
    surface_fraction = ma.ones((domain.nsurface_fraction.size, domain.y.size, domain.x.size))

    surface_fraction[0, :, :] = ma.where(vegetation_type.mask, 0.0, 1.0)
    surface_fraction[1, :, :] = ma.where(pavement_type.mask, 0.0, 1.0)
    surface_fraction[2, :, :] = ma.where(water_type.mask, 0.0, 1.0)

    domain.surface_fraction.to_nc(surface_fraction)
    domain.vegetation_type.to_nc(vegetation_type)
    domain.pavement_type.to_nc(pavement_type)
    domain.water_type.to_nc(water_type)
    domain.soil_type.to_nc(soil_type)


def process_street_type_crossing(domain: CSDDomain) -> None:
    """Process street type and street crossings of a given domain.

    Read street type and street crossings and make fields consistent with pavement type. All data is
    written to the result file.

    Args:
        domain: Domain to process.
    """
    logger.status("Processing street types and street crossings.")

    street_type = domain.read(InputData.street_type)
    pavement_type = domain.pavement_type.from_nc()
    street_type.mask = ma.mask_or(pavement_type.mask, street_type.mask)

    domain.street_type.to_nc(street_type)

    street_crossings = domain.read(InputData.street_crossings)
    domain.street_crossing.to_nc(street_crossings)


def process_resolved_vegetation(
    domain: CSDDomain,
    settings: CSDConfigSettings,
) -> None:
    """Process resolved vegetation of a given domain.

    Call function to process single trees and vegetation patches. Adopt vegetation_type and LAI
    accordingly. All data is written to the result file.

    Args:
        domain: Domain to process.
        settings: General settings.
        canopy_generator: Canopy generator to generate LAD and BAD fields.
    """
    if domain.config.generate_single_trees:
        process_single_trees(domain, settings)
    if domain.config.generate_vegetation_patches:
        process_vegetation_patches(domain, settings)

    lai = domain.read(InputData.lai)
    vegetation_height = domain.read(InputData.vegetation_height)
    vegetation_type = domain.vegetation_type.from_nc()

    # Replace high vegetation types that have LAD/BAD above it by user specified vegetation type.
    vegetation_type = ma.where(
        domain.is_resolved_vegetation2d() & ma_isin(vegetation_type, VT_HIGH_VEGETATION),
        settings.vegetation_type_below_trees,
        vegetation_type,
    )

    # Treat remaining high vegetation pixels.
    if domain.config.replace_high_vegetation_types:
        # Replace all remaining high vegetation pixels by short grass. High vegetation pixels could
        # still remain when generate_vegetation_patches=False or for pixels with low patch height,
        # which is not considered when creating vegetation patch LAD.
        high_vegetation_to_short_grass = ma_isin(vegetation_type, VT_HIGH_VEGETATION)
        if domain.config.generate_vegetation_patches:
            logger.warning_argwhere(
                "After processing vegetation patches,",
                high_vegetation_to_short_grass,
                "high vegetation pixels left. Replacing with short grass.",
            )
        else:
            logger.warning_argwhere(
                "Replacing",
                high_vegetation_to_short_grass,
                "pixels by short grass.\n"
                + "Consider generate_vegetation_patches: True or "
                + "replace_high_vegetation_types: False.",
            )
    else:
        # Remove vegetation_type when vegetation height indicates low vegetation
        high_vegetation_to_short_grass = (
            vegetation_height < settings.height_high_vegetation_lower_threshold
        ).filled(False) & ma_isin(vegetation_type, VT_HIGH_VEGETATION)
        logger.warning_argwhere(
            "Replacing",
            high_vegetation_to_short_grass,
            "high vegetation pixels with a vegetation height "
            + f"smaller than {settings.height_high_vegetation_lower_threshold}m by short grass.",
        )
    vegetation_type = ma.where(
        high_vegetation_to_short_grass,
        IndexVegetationType.short_grass,
        vegetation_type,
    )

    # Derive LAI used by LSM. This includes only low vegetation when replace_high_vegetation_types,
    # both low and high vegetation if not. LSM LAI has vegetation type specific defaults so we do
    # not have to define a value here.
    # Set default LAI for pixels with resolved vegetation.
    if settings.lai_low_vegetation_default is None:
        lai_lsm = ma.where(domain.is_resolved_vegetation2d(), ma.masked, lai)
    else:
        lai_lsm = ma.where(
            domain.is_resolved_vegetation2d(), settings.lai_low_vegetation_default, lai
        )

    # Fill remaining high vegetation pixels without LAI or with LAI = 0 with default value for high
    # vegetation.
    if settings.lai_high_vegetation_default is not None:
        lai_lsm = ma.where(
            (lai_lsm.mask | (lai_lsm == 0.0).filled(False))
            & ~domain.is_resolved_vegetation2d()
            & ma_isin(vegetation_type, VT_HIGH_VEGETATION),
            settings.lai_high_vegetation_default,
            lai_lsm,
        )

    # Fill remaining (low vegetation) pixels without LAI or with LAI = 0 with default value for low
    # vegetation.
    if settings.lai_low_vegetation_default is not None:
        lai_lsm = ma.where(
            lai_lsm.mask
            | (lai_lsm == 0.0).filled(False) & ~ma_isin(vegetation_type, VT_HIGH_VEGETATION),
            settings.lai_low_vegetation_default,
            lai_lsm,
        )

    # Remove lai for pixels that have no vegetation_type or no plants in general.
    ma.masked_where(
        vegetation_type.mask | ma_isin(vegetation_type, VT_NO_PLANTS).filled(False),
        lai_lsm,
        copy=False,
    )

    vegetation_pars = domain.vegetation_pars.empty_array()
    vegetation_pars[IndexVegetationPars.lai, :, :] = lai_lsm

    domain.vegetation_pars.to_nc(vegetation_pars)
    domain.vegetation_type.to_nc(vegetation_type)

    # Write results to file and remove from memory.
    if domain.zlad.values is not None:
        domain.lad.to_nc()
        domain.bad.to_nc()
        domain.tree_id.to_nc()
        domain.tree_type.to_nc()

        domain.lad.values = None
        domain.bad.values = None
        domain.tree_id.values = None
        domain.tree_type.values = None


def process_single_trees(
    domain: CSDDomain,
    settings: CSDConfigSettings,
) -> None:
    """Process single trees of a given domain.

    Read the single tree data and identify single trees. Create DomainTree objects for each single
    tree. Create empty domain global LAD and BAD fields. Add LAD and BAD of each single tree
    directly into the respective global one. Keep the global LAD and BAD fields in memory for
    further processing.

    Args:
        domain: Domain to process.
        settings: General settings.
        canopy_generator: Canopy generator to generate LAD and BAD fields.
    """
    logger.status("Processing single trees.")

    # Domain-wide fields for if respective tree value is not defined
    lai = domain.read(InputData.lai)
    vegetation_height = domain.read(InputData.vegetation_height)

    trees: List[DomainTree] = []

    # Tree data from shape file
    tree_points = domain.read_trees()
    if tree_points is not None:
        # Use shape data.
        number_of_trees = len(tree_points)
        if number_of_trees == 0:
            logger.info("No single trees found.")
            return
        logger.info(f"Found {number_of_trees} trees.")

        # Create a DomainTree for each single tree
        for tree_point in tree_points.itertuples(index=False):
            # Use vegetation_height if tree height is not defined. vegetation_height does not define
            # a single tree.
            tree_height = getattr(tree_point, InputData.tree_height)
            if settings.use_vegetation_height_for_trees and pd.isna(tree_height):
                tree_height = vegetation_height[
                    getattr(tree_point, ColumnGeneral.y_index),
                    getattr(tree_point, ColumnGeneral.x_index),
                ]

            # Use LAI if tree lai is not defined. lai does not define a single tree.
            tree_lai = getattr(tree_point, InputData.tree_lai)
            if settings.use_lai_for_trees and pd.isna(tree_lai):
                tree_lai = lai[
                    getattr(tree_point, ColumnGeneral.y_index),
                    getattr(tree_point, ColumnGeneral.x_index),
                ]
                # If value is missing in the input LAI field, estimate tree_lai depending on
                # tree_height. Use lai_per_vegetation_height * tree_height if
                # estimate_lai_from_height, otherwise use the corresponding
                # lai_low/high_vegetation_default.
                if is_missing(tree_lai) and not is_missing(tree_height):
                    if domain.config.estimate_lai_from_vegetation_height:
                        tree_lai = settings.lai_per_vegetation_height * tree_height
                    else:
                        if tree_height < settings.height_high_vegetation_lower_threshold:
                            if settings.lai_low_vegetation_default is not None:
                                tree_lai = settings.lai_low_vegetation_default
                        else:
                            if settings.lai_high_vegetation_default is not None:
                                tree_lai = settings.lai_high_vegetation_default

            tree = domain.canopy_generator.generate_tree(
                i=getattr(tree_point, ColumnGeneral.x_index),
                j=getattr(tree_point, ColumnGeneral.y_index),
                type=getattr(tree_point, InputData.tree_type),
                shape=getattr(tree_point, InputData.tree_shape),
                height=tree_height,
                lai=tree_lai,
                crown_diameter=getattr(tree_point, InputData.tree_crown_diameter),
                trunk_diameter=getattr(tree_point, InputData.tree_trunk_diameter),
            )
            if tree is not None:
                trees.append(tree)
    else:
        # Use raster data.
        # Read all tree parameters from file. They are defined at the centre of the tree.
        # Data correction and modification is done in generate_tree below.
        tree_crown_diameter_centre = domain.read(InputData.tree_crown_diameter)
        tree_height_centre = domain.read(InputData.tree_height)
        tree_lai_centre = domain.read(InputData.tree_lai)
        tree_shape_centre = domain.read(InputData.tree_shape)
        tree_trunk_diameter_centre = domain.read(InputData.tree_trunk_diameter)
        tree_type_centre = domain.read(InputData.tree_type)

        # Centre of a tree?
        tree_pixels = np.where(
            ~ma.getmaskarray(tree_crown_diameter_centre)
            | ~ma.getmaskarray(tree_height_centre)
            | ~ma.getmaskarray(tree_lai_centre)
            | ~ma.getmaskarray(tree_shape_centre)
            | ~ma.getmaskarray(tree_trunk_diameter_centre)
            | ~ma.getmaskarray(tree_type_centre),
            True,
            False,
        )

        number_of_trees = np.sum(tree_pixels)
        if number_of_trees == 0:
            logger.info("No single trees found.")
            return
        logger.info(f"Found {number_of_trees} trees.")

        # Create a DomainTree for each single tree
        for j, i in np.argwhere(tree_pixels):
            # Use vegetation_height if tree_height_centre is not defined. vegetation_height does not
            # define a single tree.
            tree_height = tree_height_centre[j, i]
            if settings.use_vegetation_height_for_trees and ma.is_masked(tree_height):
                tree_height = vegetation_height[j, i]

            # Use LAI if tree_lai_centre is not defined. lai does not define a single tree.
            tree_lai = tree_lai_centre[j, i]
            if settings.use_lai_for_trees and ma.is_masked(tree_lai):
                tree_lai = lai[j, i]

            tree = domain.canopy_generator.generate_tree(
                i=i,
                j=j,
                type=tree_type_centre[j, i],
                shape=tree_shape_centre[j, i],
                height=tree_height,
                lai=tree_lai,
                crown_diameter=tree_crown_diameter_centre[j, i],
                trunk_diameter=tree_trunk_diameter_centre[j, i],
            )
            if tree is not None:
                trees.append(tree)

    domain.canopy_generator.check_tree_counters()

    if not trees:
        logger.warning("No valid trees left.")
        return

    max_tree_height = max(tree.height for tree in trees)

    # Create array for vegetation canopy heights, might be extended later for vegetation patches
    zlad = ma.arange(
        0,
        math.floor(max_tree_height / domain.config.dz) * domain.config.dz + 2 * domain.config.dz,
        domain.config.dz,
    )
    zlad[1:] = zlad[1:] - 0.5 * domain.config.dz
    domain.zlad.values = zlad

    # Create common arrays for LAD and BAD as well as arrays for tree IDs and types
    # use NCDFVariable values to save for next routine
    domain.lad.values = domain.lad.empty_array()
    domain.bad.values = domain.bad.empty_array()
    domain.tree_id.values = domain.tree_id.empty_array()
    domain.tree_type.values = domain.tree_type.empty_array()

    for tree in trees:
        domain.canopy_generator.add_tree_to_3d_fields(
            tree,
            domain.lad.values,
            domain.bad.values,
            domain.tree_id.values,
            domain.tree_type.values,
        )

    # Remove LAD volumes that are inside buildings
    if not domain.config.overhanging_trees:
        buildings_2d = domain.buildings_2d.from_nc()
        building_col_3d = np.repeat(
            ~ma.getmaskarray(buildings_2d)[np.newaxis, :, :], domain.lad.values.shape[0], axis=0
        )

        ma.masked_where(building_col_3d, domain.lad.values, copy=False)
        ma.masked_where(building_col_3d, domain.bad.values, copy=False)
        ma.masked_where(building_col_3d, domain.tree_id.values, copy=False)
        ma.masked_where(building_col_3d, domain.tree_type.values, copy=False)


def process_vegetation_patches(
    domain: CSDDomain,
    settings: CSDConfigSettings,
) -> None:
    """Process vegetation patches of a given domain.

    Read patch type, vegetation height (interpreted as patch height), vegetation type and LAI.
    Vegetation patches are identified by three criteria:

    1. patch height >= height_rel_resolved_vegetation_lower_threshold * dz
    2. if replace_high_vegetation_types, high vegetation vegetation type with either undefined
      patch height or patch height >= height_rel_resolved_vegetation_lower_threshold * dz
    3. defined patch type with either undefined patch height or
      patch height >= height_rel_resolved_vegetation_lower_threshold * dz

    If overhanging_trees is False, all identified patch pixels without vegetation type are
    disregarded. Pixels with already defined LAD or BAD values are ignored. For the found vegetation
    patch pixels, calculate LAD fields and add to the existing global LAD field or create new global
    LAD and BAD fields. Keep the global LAD and BAD fields in memory for further processing.

    Args:
        domain: Domain to process.
        settings: General settings.
        canopy_generator: Canopy generator to generate LAD and BAD fields.
    """
    logger.status("Processing vegetation patches.")

    # Load vegetation patch related data.
    patch_height = domain.read(InputData.vegetation_height)
    vegetation_type = domain.vegetation_type.from_nc()
    lai = domain.read(InputData.lai)
    patch_type = domain.read(InputData.patch_type)

    # Identify vegetation patches.

    # Option 1: patch_height >= height_rel_resolved_vegetation_lower_threshold * dz.
    # In order to produce a warning message about the ignored pixels with patch height <
    # height_rel_resolved_vegetation_lower_threshold * dz, find all pixels with patch height first
    # and compare with the filtered pixels.
    patch_prelim = ~ma.getmaskarray(patch_height)
    patch = (
        patch_height >= settings.height_rel_resolved_vegetation_lower_threshold * domain.config.dz
    ).filled(False)
    # Low height pixels with patch height < height_rel_resolved_vegetation_lower_threshold * dz
    low_height = patch_prelim & ~patch
    logger.warning_argwhere(
        "Found",
        low_height,
        "pixels with patch height < "
        + f"{settings.height_rel_resolved_vegetation_lower_threshold} dz.\n"
        + "They are not treated as vegetation patches.",
    )
    # Remember when pixels are removed from the vegetation patches
    filtered = low_height.any()

    # Option 2: high vegetation vegetation type with either undefined patch height or patch_height
    # >= height_rel_resolved_vegetation_lower_threshold * dz.
    # In order to produce a warning message about the ignored pixels with patch height <
    # height_rel_resolved_vegetation_lower_threshold * dz, find all pixels with high vegetation
    # first and compare with the filtered pixels.
    if domain.config.replace_high_vegetation_types:
        high_vegetation_prelim = ma_isin(vegetation_type, VT_HIGH_VEGETATION).filled(False)
        high_vegetation = high_vegetation_prelim & ~low_height
        logger.warning_argwhere(
            "Of the the pixels with patch height < "
            + f"{settings.height_rel_resolved_vegetation_lower_threshold} * dz,",
            high_vegetation_prelim & ~high_vegetation,
            "pixels are of a high vegetation type.\n"
            + "They are not treated as vegetation patches.",
        )
        # Add high vegetation pixels.
        patch = patch | high_vegetation

    # Option 3: defined patch type with either undefined patch height or patch height >=
    # RESOLVED_VEGETATION_HEIGHT_MIN_REL * dz.
    # In order to produce a warning message about the ignored pixels with patch height <
    # RESOLVED_VEGETATION_HEIGHT_MIN_REL dz, find all pixels with patch type first and compare with
    # the filtered pixels.
    type_defined_prelim = ~ma.getmaskarray(patch_type)
    type_defined = type_defined_prelim & ~low_height
    logger.warning_argwhere(
        "Of the the pixels with patch height < "
        + f"{settings.height_rel_resolved_vegetation_lower_threshold} * dz,",
        type_defined_prelim & ~type_defined,
        "pixels have a defined patch type.\n" + "They are not treated as vegetation patches.",
    )
    # Add defined patch type pixels.
    patch = patch | type_defined

    # If overhanging trees are not allowed, remove pixels with undefined vegetation type.
    if not domain.config.overhanging_trees:
        vegetation_type_defined = ~ma.getmaskarray(vegetation_type)
        logger.warning_argwhere(
            f"Of the {'remaining ' if filtered else ''}identified vegetation patch pixels,",
            patch & ~vegetation_type_defined,
            "have undefined vegetation type.\n"
            + "They are not treated as vegetation patches "
            + "due to overhanging_trees: False.",
        )
        patch = patch & vegetation_type_defined
        # Remember when pixels are removed from the vegetation patches
        filtered = filtered or (patch & ~vegetation_type_defined).any()

    # Remove pixels with already set trees.
    single_trees = domain.is_resolved_vegetation2d()
    logger.warning_argwhere(
        f"Of the {'remaining ' if filtered else ''}identified vegetation patch pixels,",
        patch & single_trees,
        "are already covered by single tree pixels.\n"
        + "They are not treated as vegetation patches.",
    )
    patch = patch & ~single_trees
    # Remember when pixels are removed from the vegetation patches
    filtered = filtered or (patch & single_trees).any()

    # Check if any vegetation patches are left.
    if not patch.any():
        logger.info("No vegetation patches found.")
        return

    # Define missing patch type, patch height and LAI after vegetation patch identification.

    # For missing patch type, use vegetation type in case of high vegetation vegetation type. Use
    # negative value to differentiate from tree type values. This gives at least some information on
    # the patch type.
    patch_type = ma.where(
        patch_type.mask & ma_isin(vegetation_type, VT_HIGH_VEGETATION).filled(False),
        -vegetation_type,
        patch_type,
    )
    # For still missing patch type, set default value.
    patch_type = ma.where(patch_type.mask, value_defaults["patch_type"].default, patch_type)

    # Set missing patch heights to default.
    patch_height = ma.where(patch_height.mask, settings.patch_height_default, patch_height)

    # For missing LAI, use lai_per_vegetation_height * patch_height or, if defined, default values
    # for low and high vegetation.
    if domain.config.estimate_lai_from_vegetation_height:
        lai = ma.where(lai.mask, settings.lai_per_vegetation_height * patch_height, lai)
    else:
        if settings.lai_low_vegetation_default is not None:
            lai = ma.where(
                lai.mask
                & (patch_height < settings.height_high_vegetation_lower_threshold).filled(False),
                settings.lai_low_vegetation_default,
                lai,
            )
        if settings.lai_high_vegetation_default is not None:
            lai = ma.where(
                lai.mask
                & (patch_height >= settings.height_high_vegetation_lower_threshold).filled(False),
                settings.lai_high_vegetation_default,
                lai,
            )

    patch_no_lai = patch & lai.mask
    logger.warning_argwhere(
        f"Of the {'remaining ' if filtered else ''}identified vegetation patch pixels,",
        patch_no_lai,
        "have no LAI values, neither from lai input nor derived from vegetation_height.\n"
        + "They are not treated as vegetation patches. Alternatively, define "
        + "lai_low_vegetation_default and/or lai_high_vegetation_default.",
    )
    patch = patch & ~patch_no_lai
    # Remember when pixels are removed from the vegetation patches.
    filtered = filtered or (patch & single_trees).any()

    # Only keep pixels with vegetation patches.
    ma.masked_where(~patch, patch_type, copy=False)
    ma.masked_where(~patch, patch_height, copy=False)
    ma.masked_where(~patch, lai, copy=False)

    # Calculate result fields for vegetation patches.
    lad_patch, patch_id, patch_type_3d = domain.canopy_generator.process_patch(
        patch_height,
        patch_type,
        lai,
    )

    # Update global resolved vegetation fields. Check if resolved vegetation is already present.
    # If so, merge current and former data.
    if (
        domain.zlad.values is not None
        and domain.tree_id.values is not None
        and domain.tree_type.values is not None
        and domain.lad.values is not None
        and domain.bad.values is not None
    ):
        # Need to merge data.
        # Check zlad size and adjust data if necessary.
        nz_diff = lad_patch.shape[0] - domain.zlad.values.size
        if nz_diff < 0:
            # If former resolved vegetation fields are larger than current ones, extend current
            # ones.
            fillup = ma.masked_all((-nz_diff, domain.y.size, domain.x.size))
            lad_patch = ma.concatenate((lad_patch, fillup), axis=0)
            patch_id = ma.concatenate((patch_id, fillup), axis=0)
            patch_type_3d = ma.concatenate((patch_type_3d, fillup), axis=0)
        elif nz_diff > 0:
            # If current resolved vegetation fields are larger than former ones, extend former
            # ones.
            zlad = ma.arange(lad_patch.shape[0]) * domain.config.dz
            zlad[1:] = zlad[1:] - 0.5 * domain.config.dz
            domain.zlad.values = zlad

            fillup = ma.masked_all((nz_diff, domain.y.size, domain.x.size))
            domain.lad.values = ma.MaskedArray(ma.concatenate((domain.lad.values, fillup), axis=0))
            domain.bad.values = ma.MaskedArray(ma.concatenate((domain.bad.values, fillup), axis=0))
            domain.tree_id.values = ma.MaskedArray(
                ma.concatenate((domain.tree_id.values, fillup), axis=0)
            )
            domain.tree_type.values = ma.MaskedArray(
                ma.concatenate((domain.tree_type.values, fillup), axis=0)
            )

        # Add current fields to former fields.
        # Use negative patch_id to distinguish from tree_id.
        domain.tree_id.values = ma.where(~lad_patch.mask, -1.0 * patch_id, domain.tree_id.values)
        domain.tree_type.values = ma.where(~lad_patch.mask, patch_type_3d, domain.tree_type.values)
        domain.lad.values = ma.where(~lad_patch.mask, lad_patch, domain.lad.values)
    else:
        # Create global resolved vegetation fields.

        zlad = ma.arange(lad_patch.shape[0]) * domain.config.dz
        zlad[1:] = zlad[1:] - 0.5 * domain.config.dz
        domain.zlad.values = zlad

        # Use negative patch_id to distinguish from tree_id.
        domain.tree_id.values = -1.0 * patch_id
        domain.tree_type.values = patch_type_3d
        domain.lad.values = lad_patch
        domain.bad.values = ma.masked_all_like(lad_patch)


def process_water_temperature(domain: CSDDomain) -> None:
    """Process water temperatures of a given domain.

    Read water type and water temperature. Use config values and input water temperature to set
    output water temperature. Write water_pars to the result file.

    Args:
        domain: Domain to process.
    """
    logger.status("Processing water temperatures.")

    # Read water type from output file and create water_pars.
    water_type = domain.water_type.from_nc()
    water_pars = domain.water_pars.empty_array()

    # Set specific water temperature per type as assigned in config.
    if domain.config.water_temperature is not None:
        for (
            water_type_index,
            water_temperature_from_config,
        ) in domain.config.water_temperature.items():
            water_pars[IndexWaterPars.water_temperature, :, :] = ma.where(
                water_type == water_type_index,
                water_temperature_from_config,
                water_pars[IndexWaterPars.water_temperature, :, :],
            )

    # Set water temperature based on input file.
    water_temperature_from_file = domain.read(InputData.water_temperature)
    if not water_temperature_from_file.mask.all():
        water_temperature_from_file.mask = ma.mask_or(
            water_temperature_from_file.mask, water_type.mask
        )
        water_pars[IndexWaterPars.water_temperature, :, :] = ma.where(
            ~water_temperature_from_file.mask,
            water_temperature_from_file,
            water_pars[IndexWaterPars.water_temperature, :, :],
        )

    domain.water_pars.to_nc(water_pars)


def consistency_check_update_surface_fraction(domain: CSDDomain) -> None:
    """Do consistency check and update surface fractions for a given domain.

    Args:
        domain: Domain to process.
    """
    vegetation_type = domain.vegetation_type.from_nc()
    pavement_type = domain.pavement_type.from_nc()
    building_type = domain.building_type.from_nc()
    water_type = domain.water_type.from_nc()
    soil_type = domain.soil_type.from_nc()

    # Check for consistency and fill empty fields with default vegetation type.
    consistency_array, test = check_consistency_4(
        vegetation_type, building_type, pavement_type, water_type
    )

    # Check for consistency and fill empty fields with default vegetation type.
    consistency_array, test = check_consistency_3(vegetation_type, pavement_type, soil_type)

    surface_fraction = domain.surface_fraction.from_nc()
    surface_fraction[0, :, :] = ma.where(vegetation_type.mask, 0.0, 1.0)
    surface_fraction[1, :, :] = ma.where(pavement_type.mask, 0.0, 1.0)
    surface_fraction[2, :, :] = ma.where(water_type.mask, 0.0, 1.0)
    domain.surface_fraction.to_nc(surface_fraction)


def print_unused_input_files(inputs: Dict[str, CSDConfigInput]) -> None:
    """Print unread files.

    Args:
        inputs: Input files.
    """
    for name, input in inputs.items():
        unused = input.unused_file()
        unused_len = len(unused)
        if unused_len > 0:
            unused_items = "\n".join(f"  {item[0]}: {item[1]}" for item in unused)
            logger.warning(
                f"In input_{name}, the following "
                + f"{'input was' if unused_len == 1 else 'inputs were'} not used:\n"
                + unused_items
            )
