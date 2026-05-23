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

"""Main module of palm_csd to create a PALM static driver.

create_driver is the main routine, which calls the necessary functions to create a static driver in
this module.
"""

import logging
import math
from os import PathLike
from typing import Dict, List, Optional, Union, cast

import numpy as np
import numpy.ma as ma
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
from palm_csd.csd_config import (
    NBUILDING_SURFACE_LAYER,
    VT_HIGH_VEGETATION,
    CSDConfig,
    CSDConfigSettings,
    IndexBuildingType,
    IndexVegetationType,
    IndexWaterPars,
    defaults,
)
from palm_csd.csd_domain import (
    CSDDomain,
    IndexBuildingSurfaceLevel,
    IndexBuildingSurfaceType,
)
from palm_csd.lcz import LCZTypes
from palm_csd.statistics import static_driver_statistics
from palm_csd.tools import (
    blend_array_2d,
    check_consistency_3,
    check_consistency_4,
    height_to_z_grid,
    interpolate_2d,
    ma_isin,
)
from palm_csd.vegetation import CanopyGenerator, DomainTree

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
    config.update_defaults()

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

    # Initialize domain tree's defaults from the tree database.
    DomainTree.populate_defaults()

    # Create CanopyGenerator with the LAD method and parameters from the configuration.
    canopy_generator = CanopyGenerator(
        method=config.settings.lad_method,
        alpha_Metal2003=config.settings.lad_alpha,
        beta_Metal2003=config.settings.lad_beta,
        z_max_rel_LM2004=config.settings.lad_z_max_rel,
    )

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
    # the terrain height of the parent domain is needed. domains is constructed in such a way that
    # parents are always dealt with before the children.
    for domain in domains.values():
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

            process_buildings_bridges(domain)

            process_types(domain)
            process_street_type_crossing(domain)

            if domain.config.vegetation_on_roofs:
                process_vegetation_roof(domain, config.settings)
            process_resolved_vegetation(domain, config.settings, canopy_generator)

            if (
                domain.config.water_temperature is not None
                or domain.input_config.file_water_temperature is not None
            ):
                process_water_temperature(domain)

            consistency_check_update_surface_fraction(domain)
            domain.write_global_attributes()

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
        zt = domain.read_zt()
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
        e_UTM, n_UTM, lon, lat = domain.geo_converter.geographic_coordinates()

        # Write CRS
        domain.write_crs_to_file()

    else:
        # Get coordinates near origin
        if domain.x0 is None or domain.y0 is None:
            raise ValueError(f"Domain {domain.name} has no x0 or y0 defined")
        x_UTM_origin = domain.read_nc_2d(
            domain.input_config.file_x_UTM,
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        y_UTM_origin = domain.read_nc_2d(
            domain.input_config.file_y_UTM,
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        lat_origin = domain.read_nc_2d(
            domain.input_config.file_lat,
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )
        lon_origin = domain.read_nc_2d(
            domain.input_config.file_lon,
            x0=domain.x0,
            x1=domain.x0 + 1,
            y0=domain.y0,
            y1=domain.y0 + 1,
        )

        # Calculate position of origin. Added as global attributes later
        domain.origin_x = float(x_UTM_origin[0, 0]) - 0.5 * (
            float(x_UTM_origin[0, 1]) - float(x_UTM_origin[0, 0])
        )
        domain.origin_y = float(y_UTM_origin[0, 0]) - 0.5 * (
            float(y_UTM_origin[1, 0]) - float(y_UTM_origin[0, 0])
        )
        domain.origin_lon = float(lon_origin[0, 0]) - 0.5 * (
            float(lon_origin[0, 1]) - float(lon_origin[0, 0])
        )
        domain.origin_lat = float(lat_origin[0, 0]) - 0.5 * (
            float(lat_origin[1, 0]) - float(lat_origin[0, 0])
        )

        # Read x and y values
        domain.x_global.values = domain.read_nc_1d(domain.input_config.file_x_UTM, "x")
        domain.y_global.values = domain.read_nc_1d(
            domain.input_config.file_y_UTM, "y", x0=domain.y0, x1=domain.y1
        )

        # Read and write lon, lat and UTM coordinates
        lat = domain.read_nc_2d(domain.input_config.file_lat)
        lon = domain.read_nc_2d(domain.input_config.file_lon)

        e_UTM = domain.read_nc_2d(domain.input_config.file_x_UTM)
        n_UTM = domain.read_nc_2d(domain.input_config.file_y_UTM)

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

    domain.E_UTM.to_nc(e_UTM)
    domain.N_UTM.to_nc(n_UTM)

    # Read and process terrain height (zt). Its values are stored in the domain object to be
    # available for potential child domain.
    domain.zt.values = domain.read_zt()
    domain.zt.values = domain.zt.values - zt_min
    domain.origin_z = float(zt_min)

    # If necessary, interpolate parent domain terrain height on child domain grid and blend
    # the two.
    if domain.config.interpolate_terrain:
        if domain.parent is None:
            raise ValueError("Interpolation of terrain height requires a parent domain")
        if domain.parent.x_global.values is None or domain.parent.y_global.values is None:
            raise ValueError(f"x_UTM or y_UTM of parent {domain.parent.name} not calculated")
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


def process_buildings_bridges(domain: CSDDomain) -> None:
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

    Raises:
        ValueError: Building IDs are missing for some building pixels.
        ValueError: Bridge IDs are missing for some bridge pixels.
        ValueError: Length of building parameter data does not match number of building surface
          layers.
        ValueError: Invalid value for a key in the building parameter data.
        ValueError: Invalid dimension for the building parameter output variable.
    """
    logger.status("Processing buildings and bridges.")

    buildings_2d = domain.read_buildings_2d()
    building_id = domain.read_building_id()
    building_type = domain.read_building_type()

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
        bridges_2d = domain.read_bridges_2d()
        bridges_id = domain.read_bridges_id()

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

        buildings_bridges_overlap = ~buildings_2d.mask & ~bridges_2d.mask
        logger.warning_argwhere(
            "Buildings and bridges are overlapping at",
            buildings_bridges_overlap,
            "pixels.\n" + "Prefering building information at these pixels.",
        )

        bridges_id.mask = bridges_2d.mask.copy()
        building_id = ma.where(buildings_2d.mask & ~bridges_2d.mask, bridges_id, building_id)
        building_type = ma.where(
            buildings_2d.mask & ~bridges_2d.mask, IndexBuildingType.BRIDGES, building_type
        )

    domain.buildings_2d.to_nc(buildings_2d)
    domain.building_id.to_nc(building_id)
    domain.building_type.to_nc(building_type)

    # Create 3d buildings if necessary. Add bridge pixels to building layer.
    if domain.config.buildings_3d or (bridge_depth_grid > 0 and not bridges_2d.mask.all()):
        if not domain.config.buildings_3d:
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

    def add_where_buildings(
        input: Optional[
            Union[
                Dict[int, float],
                Dict[int, int],
                Dict[int, List[float]],
                Dict[int, List[int]],
            ]
        ],
        output_variable: netcdf_data.NCDFVariable,
    ):
        """Add input data to output_variable where buildings are present.

        The input data is a dictionary with keys corresponding to the output layer name. The values
        represent the building surface layers if they are a list. The data is written to the result
        file.

        Args:
            input: Input data. key: output layer index, value: singe value or building surface
              layer values.
            output_variable: Output variable to write the data to.

        Raises:
            ValueError: Length of input data does not match number of building surface layers.
            ValueError: Invalid value for a key in the input data.
            ValueError: Invalid dimension for the output variable.
        """
        if input is None:
            return

        output_field = output_variable.empty_array()
        # with surface layers
        if output_field.ndim == 4:
            for key, value in input.items():
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
                    output_field[key, i, :, :] = ma.where(~buildings_2d.mask, v, ma.masked)
        # without surface layers
        elif output_field.ndim == 3:
            for key, value in input.items():
                if not isinstance(value, (int, float)):
                    raise ValueError(f"Invalid value for {key}")
                output_field[key, :, :] = ma.where(~buildings_2d.mask, value, ma.masked)
        else:
            raise ValueError(f"Invalid dimension for {output_variable.name}")

        output_variable.to_nc(output_field)

    add_where_buildings(domain.config.building_albedo_type, domain.building_albedo_type)
    add_where_buildings(domain.config.building_emissivity, domain.building_emissivity)
    add_where_buildings(domain.config.building_fraction, domain.building_fraction)
    add_where_buildings(domain.config.building_general_pars, domain.building_general_pars)
    add_where_buildings(domain.config.building_heat_capacity, domain.building_heat_capacity)
    add_where_buildings(domain.config.building_heat_conductivity, domain.building_heat_conductivity)
    add_where_buildings(domain.config.building_indoor_pars, domain.building_indoor_pars)
    add_where_buildings(domain.config.building_lai, domain.building_lai)
    add_where_buildings(domain.config.building_roughness_length, domain.building_roughness_length)
    add_where_buildings(
        domain.config.building_roughness_length_qh, domain.building_roughness_length_qh
    )
    add_where_buildings(domain.config.building_thickness, domain.building_thickness)
    add_where_buildings(domain.config.building_transmissivity, domain.building_transmissivity)


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
    lcz_type = domain.read_lcz(lcz_types)
    # LAI data input
    lai = domain.read_lai()
    # LAI from LCZ table
    lai_lcz = lcz_types.lai_from_lcz_map(lcz_type)
    # Use LAI from LCZ table if LAI data is missing
    lai = ma.where(lai.mask, lai_lcz, lai)

    soil_type = domain.read_soil_type()

    water_type = lcz_types.water_type_from_lcz_map(lcz_type)
    vegetation_type = lcz_types.vegetation_type_from_lcz_map(lcz_type)
    pavement_type = ma.masked_all_like(water_type)

    lai.mask = ma.mask_or(vegetation_type.mask, lai.mask)

    domain.nvegetation_pars.values = ma.arange(0, 12)
    vegetation_pars = ma.masked_all((domain.nvegetation_pars.size, domain.y.size, domain.x.size))
    vegetation_pars[1, :, :] = lai

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

    vegetation_type = domain.read_vegetation_type()
    pavement_type = domain.read_pavement_type()
    water_type = domain.read_water_type()
    soil_type = domain.read_soil_type()
    # Use buildings_2d because it does not include bridges unlike building type
    building_height = domain.buildings_2d.from_nc()

    # Make arrays consistent
    # #1 Set vegetation type to masked for pixel where a pavement type is set
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~pavement_type.mask)

    # #2 Set vegetation type to masked for pixel where a building type is set
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~building_height.mask)

    # #3 Set vegetation type to masked for pixel where a water type is set
    vegetation_type.mask = ma.mask_or(vegetation_type.mask, ~water_type.mask)

    # #4 Remove pavement for pixels with buildings
    pavement_type.mask = ma.mask_or(pavement_type.mask, ~building_height.mask)

    # #5 Remove pavement for pixels with water.
    pavement_type.mask = ma.mask_or(pavement_type.mask, ~water_type.mask)

    # #6 Remove water for pixels with buildings
    water_type.mask = ma.mask_or(water_type.mask, ~building_height.mask)

    # Correct vegetation_type when a vegetation height is available and is indicative of low
    # vegetation
    vegetation_height = domain.read_vegetation_height()

    # Correct vegetation_type depending on vegetation_height.
    # ma.where gives ma.masked when its first argument is ma.masked. We don't want this for
    # vegetation_height here so do an extra check and use .data.
    vegetation_type = ma.where(
        (~vegetation_height.mask)
        & (vegetation_height.data == 0.0)
        & ma_isin(vegetation_type, VT_HIGH_VEGETATION),
        3,
        vegetation_type,
    )
    ma.masked_where(
        (vegetation_height == 0.0) & ma_isin(vegetation_type, VT_HIGH_VEGETATION),
        vegetation_height,
        copy=False,
    )

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
        if defaults["vegetation_type"].default is None:
            raise ValueError("No default vegetation type defined.")
        logger.warning(
            f"Setting default vegetation_type {defaults['vegetation_type'].default} "
            + f"for {np.sum(n_type_undefined)} pixels without surface type."
        )
        logger.debug_indent("Coordinates of undefined surface types:")
        logger.debug_indent(", ".join(map(str, type_undefined)))
        vegetation_type = ma.where(
            n_type_undefined, defaults["vegetation_type"].default, vegetation_type
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

    street_type = domain.read_street_type()
    pavement_type = domain.pavement_type.from_nc()
    street_type.mask = ma.mask_or(pavement_type.mask, street_type.mask)

    domain.street_type.to_nc(street_type)

    street_crossings = domain.read_street_crossings()
    domain.street_crossing.to_nc(street_crossings)


def process_vegetation_roof(domain: CSDDomain, settings: CSDConfigSettings) -> None:
    """Process vegetation on roofs of a given domain.

    Read vegetation on roofs and make fields consistent with building fraction and LAI. All data is
    written to the result file.

    Args:
        domain: Domain to process.
        settings: General settings.
    """
    logger.status("Processing vegetation on roofs.")

    green_roofs = domain.read_vegetation_on_roofs()
    buildings_2d = domain.buildings_2d.from_nc()

    # Adjust building fraction possibly defined above.
    building_fraction = domain.building_fraction.from_nc(allow_nonexistent=True)
    building_lai = domain.building_lai.from_nc(allow_nonexistent=True)

    # Assign green fraction on roofs.
    building_fraction[IndexBuildingSurfaceType.GREEN_ROOF, :, :] = ma.where(
        (~buildings_2d.mask) & ~green_roofs.mask & (green_roofs.data != 0.0),
        1.0,
        building_fraction[IndexBuildingSurfaceType.GREEN_ROOF, :, :],
    )

    # Set wall fraction to 0 where green fraction defined.
    building_fraction[IndexBuildingSurfaceType.WALL_ROOF, :, :] = ma.where(
        ~(building_fraction.mask[IndexBuildingSurfaceType.GREEN_ROOF, :, :]),
        0.0,
        building_fraction[IndexBuildingSurfaceType.WALL_ROOF, :, :],
    )

    # Set window fraction to 0 where green fraction defined.
    building_fraction[IndexBuildingSurfaceType.WINDOW_ROOF, :, :] = ma.where(
        ~(building_fraction.mask[IndexBuildingSurfaceType.GREEN_ROOF, :, :]),
        0.0,
        building_fraction[IndexBuildingSurfaceType.WINDOW_ROOF, :, :],
    )

    # Assign leaf area index for vegetation on roofs.
    building_lai[IndexBuildingSurfaceLevel.ROOF, :, :] = ma.where(
        (~(building_fraction.mask[IndexBuildingSurfaceType.GREEN_ROOF, :, :]))
        & (green_roofs >= 0.5),
        settings.lai_roof_intensive,
        building_lai[IndexBuildingSurfaceLevel.ROOF, :, :],
    )
    building_lai[IndexBuildingSurfaceLevel.ROOF, :, :] = ma.where(
        (~(building_fraction.mask[IndexBuildingSurfaceType.GREEN_ROOF, :, :]))
        & (green_roofs < 0.5),
        settings.lai_roof_extensive,
        building_lai[IndexBuildingSurfaceLevel.ROOF, :, :],
    )

    domain.building_fraction.to_nc(building_fraction)
    domain.building_lai.to_nc(building_lai)


def process_resolved_vegetation(
    domain: CSDDomain,
    settings: CSDConfigSettings,
    canopy_generator: CanopyGenerator,
) -> None:
    """Process resolved vegetation of a given domain.

    Call function to process single trees and vegetation patches. Adopt vegetation_type and LAI
    accordingly. All data is written to the result file.

    Args:
        domain: Domain to process.
        settings: General settings.
        canopy_generator: Canopy generator to generate LAD and BAD fields.
    """
    if domain.config.street_trees:
        process_single_trees(domain, settings, canopy_generator)
    if domain.config.generate_vegetation_patches:
        process_vegetation_patches(domain, settings, canopy_generator)

    lai = domain.read_lai()
    vegetation_type = domain.vegetation_type.from_nc()

    # Remove high vegetation wherever it is replaced by a leaf area density. This should effectively
    # remove all high vegetation pixels.
    if domain.lad.values is not None:
        vegetation_type = ma.where(
            ~ma.getmaskarray(domain.lad.values)[0, :, :] & ~vegetation_type.mask,
            settings.vegetation_type_below_trees,
            vegetation_type,
        )

    # If desired, remove all high vegetation.
    if not domain.config.allow_high_vegetation:
        vegetation_type = ma.where(
            ~vegetation_type.mask & ma_isin(vegetation_type, VT_HIGH_VEGETATION),
            3,
            vegetation_type,
        )

    if domain.lad.values is not None:
        # Set default low LAI for pixels with an LAD (short grass below trees).
        lai_low = ma.where(
            ma.getmaskarray(domain.lad.values)[0, :, :], lai, settings.lai_low_vegetation_default
        )
    else:
        lai_low = lai

    # Fill low vegetation pixels without LAI set or with LAI = 0 with default value.
    lai_low = ma.where(
        lai_low.mask | (~lai_low.mask & (lai_low.data == 0.0)),
        settings.lai_low_vegetation_default,
        lai_low,
    )

    # Remove lai for pixels that have no vegetation_type.
    ma.masked_where(
        vegetation_type.mask | (vegetation_type == IndexVegetationType.BARE_SOIL),
        lai_low,
        copy=False,
    )

    vegetation_pars = domain.vegetation_pars.empty_array()
    vegetation_pars[1, :, :] = lai_low

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
    canopy_generator: CanopyGenerator,
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

    lai = domain.read_lai()

    # Read all tree parameters from file. They are defined at the centre of the tree.
    # Data correction and modification is done in generate_tree below
    tree_height_centre = domain.read_tree_height()
    tree_crown_diameter_centre = domain.read_tree_crown_diameter()
    tree_trunk_diameter_centre = domain.read_tree_trunk_diameter()
    tree_type_centre = domain.read_tree_type()

    # Centre of a tree?
    tree_pixels = np.where(
        ~ma.getmaskarray(tree_height_centre)
        | ~ma.getmaskarray(tree_type_centre)
        | ~ma.getmaskarray(tree_crown_diameter_centre)
        | ~ma.getmaskarray(tree_trunk_diameter_centre),
        True,
        False,
    )
    number_of_trees = np.sum(tree_pixels)

    if number_of_trees == 0:
        logger.info("No street trees found.")
        return

    # Create a DomainTree for each single tree
    logger.info(f"Found {number_of_trees} trees.")
    trees: List[DomainTree] = []
    # counter for tree IDs and adjusted trees
    DomainTree.reset_counter()
    for i in range(0, len(domain.x)):
        for j in range(0, len(domain.y)):
            if tree_pixels[j, i]:
                tree = DomainTree.generate_tree(
                    i=i,
                    j=j,
                    type=tree_type_centre[j, i],
                    shape=ma.masked,
                    height=tree_height_centre[j, i],
                    lai=lai[j, i],
                    crown_diameter=tree_crown_diameter_centre[j, i],
                    trunk_diameter=tree_trunk_diameter_centre[j, i],
                    config=domain.config,
                    settings=settings,
                )
                if tree is not None:
                    trees.append(tree)

    DomainTree.check_counter(domain.config, settings)

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
        canopy_generator.add_tree_to_3d_fields(
            tree,
            domain.lad.values,
            domain.bad.values,
            domain.tree_id.values,
            domain.tree_type.values,
            domain.config,
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
    canopy_generator: CanopyGenerator,
) -> None:
    """Process vegetation patches of a given domain.

    Read patch type, height and LAI and identify vegetation patches. Calculate LAD fields and add to
    the existing global LAD field or create new global LAD and BAD fields. Keep the global LAD and
    BAD fields in memory for further processing.

    Args:
        domain: Domain to process.
        settings: General settings.
        canopy_generator: Canopy generator to generate LAD and BAD fields.
    """
    logger.status("Processing vegetation patches.")

    patch_height = domain.read_patch_height()

    vegetation_type = domain.vegetation_type.from_nc()
    lai = domain.read_lai()
    patch_type_2d = domain.read_patch_type()

    # patch_type_2d: use high vegetation vegetation_type if patch_type is missing
    # Note: vegetation_type corrected above.
    patch_type_2d = ma.where(
        patch_type_2d.mask & ma_isin(vegetation_type, VT_HIGH_VEGETATION),
        -vegetation_type,
        patch_type_2d,
    )
    # patch_type_2d: use default value for the rest of the pixels.
    patch_type_2d = ma.where(patch_type_2d.mask, defaults["patch_type"].default, patch_type_2d)

    # Initialize lai_high with lai or default value for high vegetation, but masked.
    lai_high: ma.MaskedArray = ma.MaskedArray(
        np.where(lai.mask, settings.lai_high_vegetation_default, lai), mask=True
    )

    # Unmask all high vegetation pixels.
    lai_high.mask = np.where(
        ma_isin(vegetation_type, VT_HIGH_VEGETATION)
        & (patch_height.mask | (patch_height.data >= domain.config.dz)),
        False,
        lai_high.mask,
    )

    # Mask all pixels where street trees were already set.
    if domain.lad.values is not None:
        lai_high.mask = np.where(
            ma.getmaskarray(domain.lad.values)[0, :, :],
            lai_high.mask,
            True,
        )

    # Unmask all pixels where short grass is defined, but where a patch_height >= dz is found, as
    # high vegetation (often the case in backyards).
    lai_high.mask = np.where(
        ~patch_height.mask
        & (patch_height.data >= domain.config.dz)
        & ~vegetation_type.mask
        & (vegetation_type.data == 3),
        False,
        lai_high.mask,
    )

    # If overhanging trees are allowed, unmask pixels with patch_height > dz that are not included
    # in vegetation_type.
    if domain.config.overhanging_trees:
        lai_high.mask = np.where(
            ~patch_height.mask & (patch_height.data >= domain.config.dz),
            False,
            lai_high.mask,
        )

    # Define a patch height wherever it is missing.
    patch_height_high = ma.where(patch_height.mask, settings.patch_height_default, patch_height)

    # Remove pixels where street trees were already set.
    if domain.lad.values is not None:
        ma.masked_where(~ma.getmaskarray(domain.lad.values)[0, :, :], patch_height_high, copy=False)

    # Remove patch heights that have no lai_high value.
    ma.masked_where(lai_high.mask, patch_height_high, copy=False)

    # For missing LAI values, set either the high vegetation default or the low vegetation default.
    lai_high = ma.where(
        lai_high.mask & ~patch_height.mask & (patch_height.data > 2.0),
        settings.lai_high_vegetation_default,
        lai_high,
    )
    lai_high = ma.where(
        lai_high.mask & ~patch_height.mask & (patch_height.data <= 2.0),
        settings.lai_low_vegetation_default,
        lai_high,
    )

    if ma.max(patch_height_high) >= (2.0 * domain.config.dz):
        # Result fields for vegetation patches
        lad_patch, patch_id, patch_types = canopy_generator.process_patch(
            domain.config.dz,
            patch_height_high,
            patch_type_2d,
            lai_high,
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
                patch_types = ma.concatenate((patch_types, fillup), axis=0)
            elif nz_diff > 0:
                # If current resolved vegetation fields are larger than former ones, extend former
                # ones.
                zlad = ma.arange(lad_patch.shape[0]) * domain.config.dz
                zlad[1:] = zlad[1:] - 0.5 * domain.config.dz
                domain.zlad.values = zlad

                fillup = ma.masked_all((nz_diff, domain.y.size, domain.x.size))
                domain.lad.values = ma.MaskedArray(
                    ma.concatenate((domain.lad.values, fillup), axis=0)
                )
                domain.bad.values = ma.MaskedArray(
                    ma.concatenate((domain.bad.values, fillup), axis=0)
                )
                domain.tree_id.values = ma.MaskedArray(
                    ma.concatenate((domain.tree_id.values, fillup), axis=0)
                )
                domain.tree_type.values = ma.MaskedArray(
                    ma.concatenate((domain.tree_type.values, fillup), axis=0)
                )

            # Add current fields to former fields.
            # Use negative patch_id to distinguish from tree_id.
            domain.tree_id.values = ma.where(
                domain.lad.values.mask, -1.0 * patch_id, domain.tree_id.values
            )
            domain.tree_type.values = ma.where(
                domain.lad.values.mask, patch_types, domain.tree_type.values
            )
            domain.lad.values = ma.where(domain.lad.values.mask, lad_patch, domain.lad.values)
        else:
            # Create global resolved vegetation fields.

            zlad = ma.arange(lad_patch.shape[0]) * domain.config.dz
            zlad[1:] = zlad[1:] - 0.5 * domain.config.dz
            domain.zlad.values = zlad

            # Use negative patch_id to distinguish from tree_id.
            domain.tree_id.values = -1.0 * patch_id
            domain.tree_type.values = patch_types
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
            water_pars[IndexWaterPars.WATER_TEMPERATURE, :, :] = ma.where(
                water_type == water_type_index,
                water_temperature_from_config,
                water_pars[IndexWaterPars.WATER_TEMPERATURE, :, :],
            )

    # Set water temperature based on input file.
    if domain.input_config.file_water_temperature is not None:
        water_temperature_from_file = domain.read_water_temperature()
        water_temperature_from_file.mask = ma.mask_or(
            water_temperature_from_file.mask, water_type.mask
        )
        water_pars[IndexWaterPars.WATER_TEMPERATURE, :, :] = ma.where(
            ~water_temperature_from_file.mask,
            water_temperature_from_file,
            water_pars[IndexWaterPars.WATER_TEMPERATURE, :, :],
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
