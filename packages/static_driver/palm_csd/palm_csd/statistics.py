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

"""Show statistics and produce plots of a static driver netCDF file."""

import logging
import re
from enum import IntEnum
from pathlib import Path
from typing import Dict, List, Optional, Tuple, TypedDict, Union, cast

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib import colors
from netCDF4 import Dataset

from palm_csd import StatusLogger
from palm_csd.constants import (
    COLORS,
    VT_NO_PLANTS,
    IndexBuildingType,
    IndexPavementType,
    IndexVegetationType,
    IndexWaterType,
    ValuePlot,
)
from palm_csd.tools import ma_isin

logger = cast(StatusLogger, logging.getLogger(__name__))


def plot_static(
    nc_static: Dataset,
    show: bool = False,
    detailed: bool = False,
    output: Optional[Path] = None,
    title: Optional[str] = None,
    width: Optional[float] = None,
    height: Optional[float] = None,
    geo_referenced: bool = False,
) -> None:
    """Create a basic summary plot of a static driver showing surface types and buildings.

    If the width and height are both not specified, width = aspect_ratio * height + 2 is assumed,
    where aspect_ratio is the ratio of the grid points in x and y direction.

    Args:
        nc_static: Static driver netCDF file.
        show: Show the plot on screen. Defaults to False.
        detailed: If True, show detailed surface types instead of basic categories.
        output: Plot output path. Defaults to None.
        title: Plot title. Defaults to None.
        width: Plot width. Defaults to None.
        height: Plot height. Defaults to None.
        geo_referenced: If True, use georeferenced UTM coordinates.
    """
    # Get dimensions from the NetCDF file.
    n_x = len(nc_static.dimensions["x"])
    n_y = len(nc_static.dimensions["y"])
    aspect_ratio = n_x / n_y

    # Assumption: width = aspect_ratio * height + 2.
    if height is None:
        if width is None:
            height = 6
        else:
            height = (width - 2.0) / aspect_ratio
    if width is None:
        width = aspect_ratio * height + 2.0

    if geo_referenced:
        if "N_UTM" not in nc_static.variables or "E_UTM" not in nc_static.variables:
            logger.critical_raise(
                "Georeferenced UTM coordinates not found. Use pixel coordinates instead."
            )
        e_utm = nc_static.variables["E_UTM"][:]
        n_utm = nc_static.variables["N_UTM"][:]

        # Assume e_utm is constant along the y-axis and n_utm is constant along the x-axis.
        x_dim = e_utm[0, :]
        y_dim = n_utm[:, 0]
        x_label = "Easting (m)"
        y_label = "Northing (m)"
        if not np.all(e_utm == x_dim[np.newaxis, :]) or not np.all(n_utm == y_dim[:, np.newaxis]):
            # Try rotated coordinates.
            x_dim = n_utm[0, :]
            y_dim = e_utm[:, 0]
            x_label = "Northing (m)"
            y_label = "Easting (m)"
            if not np.all(n_utm == x_dim[np.newaxis, :]) or not np.all(
                e_utm == y_dim[:, np.newaxis]
            ):
                logger.critical_raise(
                    "UTM coordinates are not consistent along one dimension. "
                    + "Use pixel coordinates instead."
                )
    else:
        # Extract coordinate variables to compute grid spacing.
        x_dim = nc_static.variables["x"][:] if "x" in nc_static.variables else None
        y_dim = nc_static.variables["y"][:] if "y" in nc_static.variables else None
        x_label = "x (m)"
        y_label = "y (m)"

    if x_dim is None or y_dim is None:
        logger.warning(
            "Cannot calculate grid spacing, x or y coordinates are missing. "
            + "Using pixels as plot units."
        )
        extent: Optional[Tuple[float, float, float, float]] = None
        x_label = "x"
        y_label = "y"
    else:
        # Calculate horizontal grid spacing if the coordinates exist.
        dx = np.diff(x_dim).mean()
        dy = np.diff(y_dim).mean()
        extent = (
            x_dim[0] - dx / 2,
            x_dim[-1] + dx / 2,
            y_dim[0] - dy / 2,
            y_dim[-1] + dy / 2,
        )

    # Mapping from indices to variable type
    map_type = {
        IndexWaterType: "water_type",
        IndexVegetationType: "vegetation_type",
        IndexPavementType: "pavement_type",
        IndexBuildingType: "building_type",
        ValuePlot: "none",
    }

    # Helper functions to map types to unique ranges
    def map_variable(data, name: Optional[str] = None, subtypes: bool = False):
        """Map variable to unique range.

        This function maps the given data to a unique range based on the type of the variable. For
        example, a IndexWaterType.lake with value 1 would be mapped to 101 if detailed is True and
        to 100 if False. This ensures that all classes have unique values if detailed is True or
        have the same base value if False. Only the vegetation types that are not considered
        vegetation are always separated.

        Args:
            data: Values to map.
            name: Variable type name. Defaults to None.
            subtypes: Whether to include subtypes in the mapping. Defaults to False.

        Raises:
            ValueError: If name is not given for mapping variable.

        Returns:
            Mapped values.
        """
        # Try to get the variable name from the map_type when supplying an IntEnum
        # (e.g. IndexWaterType.lake).
        if isinstance(data, IntEnum):
            name = map_type[type(data)]

        if name is None:
            raise ValueError("name not given for mapping variable.")

        if subtypes:
            # Each surface type is mapped to a unique value.
            return getattr(ValuePlot, name).value + data
        else:
            # Each surface type is mapped to the same base value except for the non-vegetation
            # vegetation types.
            if name == "vegetation_type":
                if isinstance(data, ma.MaskedArray):
                    return ma.where(
                        ma_isin(data, VT_NO_PLANTS),
                        getattr(ValuePlot, name).value + data,
                        getattr(ValuePlot, name).value,
                    )
                elif data in VT_NO_PLANTS:
                    return getattr(ValuePlot, name).value + data
            return getattr(ValuePlot, name).value

    # Colors and labels dictionary. The keys are the mapped values, the values are dictionary of
    # color and label.
    dict_plot = {}
    for surface_colors in COLORS.values():
        for surface_type, color in surface_colors.items():
            mapped_value = map_variable(surface_type, subtypes=True)

            # Generate nice labels.
            label = surface_type.name.capitalize()
            # Replace underscores between numbers (e.g., 1951_2000 -> 1951-2000)
            label = re.sub(r"(?<=\d)_(?=\d)", "-", label)
            # Replace remaining underscores with spaces
            label = label.replace("_", " ")

            dict_plot[mapped_value] = {
                "color": color,
                "label": label,
            }

    # Update the base labels, which are not covered by the above rule.
    dict_plot[ValuePlot.vegetation_resolved]["label"] = "Vegetation (resolved)"
    dict_plot[ValuePlot.building_type]["label"] = "Buildings"
    dict_plot[ValuePlot.pavement_type]["label"] = "Pavement"
    dict_plot[ValuePlot.vegetation_type]["label"] = "Vegetation (flat)"
    dict_plot[ValuePlot.water_type]["label"] = "Water"
    if detailed:
        dict_plot[ValuePlot.vegetation_resolved]["color"] = "#21522F"

    # Create raster of surface types with their mapped values.
    raster = np.full((n_y, n_x), ValuePlot.none, dtype=np.uint16)
    for variable in ["building_type", "pavement_type", "vegetation_type", "water_type"]:
        if variable in nc_static.variables:
            data = nc_static.variables[variable][:]
            if data is not None:
                raster = np.where(
                    ma.getmaskarray(data), raster, map_variable(data, variable, subtypes=detailed)
                )
        else:
            logger.debug(f"{variable} does not exist in data.variables.")

    # Overwrite the raster where resolved vegetation is.
    resolved_vegetation_present = np.full_like(raster, False)
    for variable_tree in ["lad", "bad"]:
        if variable_tree in nc_static.variables:
            tree_field = nc_static.variables[variable_tree][:, :]
            resolved_vegetation_present = np.logical_or(
                (tree_field[1:, :, :] > 0.0).any(axis=0).filled(False),
                resolved_vegetation_present,
            )
        else:
            logger.debug(f"{variable_tree} does not exist in data.variables.")
    raster = np.where(resolved_vegetation_present, ValuePlot.vegetation_resolved, raster)

    # Create a colormap and labels based on the combined dictionary.
    unique_values = np.unique(raster)
    existing_colors = [dict_plot[val]["color"] for val in unique_values]
    existing_labels = [dict_plot[val]["label"] for val in unique_values]
    cmap = colors.ListedColormap(existing_colors)

    # Define the boundaries, tick positions and labels dynamically. Note that unique_values might
    # not be continuous.
    if len(unique_values) == 1:
        boundaries = np.array([unique_values[0] - 0.5, unique_values[0] + 0.5])
        tick_positions = unique_values
    else:
        boundaries = np.append(unique_values - 0.5, unique_values[-1] + 0.5)
        tick_positions = 0.5 * (boundaries[1:] + boundaries[:-1])
    norm = colors.BoundaryNorm(boundaries, cmap.N)
    tick_labels = existing_labels

    # Plot the raster with dynamic settings.
    plt.figure(figsize=(width, height))
    if extent is None:
        plt.imshow(raster, cmap=cmap, norm=norm, origin="lower", interpolation="nearest")
    else:
        plt.imshow(
            raster, cmap=cmap, norm=norm, origin="lower", extent=extent, interpolation="nearest"
        )

    # Ticks need to be set for the range of the colors but are visually removed.
    cbar = plt.colorbar(ticks=tick_positions, label="Surface types", boundaries=boundaries)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.tick_params(axis="y", which="both", length=0)

    if title is not None:
        plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to {output}.")
    if show:
        plt.show()
    plt.close()


def static_driver_statistics(
    nc_file: Union[str, Path],
    show_plot: bool = False,
    plot_file: Optional[Union[str, Path]] = None,
    plot_title: Optional[str] = None,
    plot_width: Optional[float] = None,
    plot_height: Optional[float] = None,
    detailed_plot: bool = False,
    geo_referenced: bool = False,
) -> None:
    """Calculate and print statistics of a static driver netCDF file.

    Args:
        nc_file: Path to the static driver netCDF file.
        show_plot: Show a plot. Defaults to False.
        plot_file: Path to plot. Defaults to None.
        plot_title: Plot title. Defaults to None.
        plot_width: Plot width. Defaults to None.
        plot_height: Plot height. Defaults to None.
        detailed_plot: If True, show detailed surface types instead of basic categories.
        geo_referenced: If True, use georeferenced UTM coordinates.
    """
    nc_static = Dataset(nc_file, mode="r")

    # Extract dimensions.
    dims = nc_static.dimensions
    n_x = dims["x"].size if "x" in dims else None
    n_y = dims["y"].size if "y" in dims else None
    n_zbuild = dims["z"].size if "z" in dims else None
    n_zlad = dims["zlad"].size if "zlad" in dims else None

    # Extract coordinate variables to compute grid spacing.
    x_dim = nc_static.variables["x"][:] if "x" in nc_static.variables else None
    y_dim = nc_static.variables["y"][:] if "y" in nc_static.variables else None
    zbuild_dim = nc_static.variables["z"][:] if "z" in nc_static.variables else None
    zlad_dim = nc_static.variables["zlad"][:] if "zlad" in nc_static.variables else None

    # Calculate horizontal grid spacing if the coordinates exist.
    dx = np.diff(x_dim).mean() if x_dim is not None else None
    dy = np.diff(y_dim).mean() if y_dim is not None else None

    # Calculate number of grid cells in domain.
    if n_x is not None and n_y is not None:
        n_domain = n_x * n_y
    else:
        n_domain = None

    if zbuild_dim is not None and len(zbuild_dim) >= 2:
        if not np.isclose(zbuild_dim[0], 0.0):
            logger.warning(
                f"z grid starts at the non-zero value {zbuild_dim[0]}. "
                + "Skipping dz calculation."
            )
            dzbuild = None
        else:
            dzbuild = 2 * zbuild_dim[1]
    else:
        dzbuild = None

    if zlad_dim is not None and len(zlad_dim) >= 2:
        if not np.isclose(zlad_dim[0], 0.0):
            logger.warning(
                f"z LAD grid starts at the non-zero value {zlad_dim[0]}. "
                + "Skipping dz calculation for LAD."
            )
            dzlad = None
        else:
            dzlad = 2 * zlad_dim[1]
    else:
        dzlad = None

    if dzbuild is not None and dzlad is not None and not np.isclose(dzbuild, dzlad):
        logger.warning("z gridspacing differ in building and lad fields. Using building dz.")

    if dx is not None and dy is not None:
        area_grid_cell = dx * dy
    else:
        area_grid_cell = None

    vol_grid_cell = None
    if area_grid_cell is not None:
        if dzbuild is not None:
            vol_grid_cell = area_grid_cell * dzbuild
        elif dzlad is not None:
            vol_grid_cell = area_grid_cell * dzlad

    # Calculate surface type fractions.
    static_vars = ["vegetation_type", "pavement_type", "building_type", "water_type", "tree_type"]
    result_types = {}
    max_count = 0
    for var in static_vars:
        # Check if variable is in dataset, skip if not.
        if var not in nc_static.variables:
            continue

        if var == "tree_type":
            variable_data = nc_static.variables[var][0, :, :]
        else:
            variable_data = nc_static.variables[var][:]

        n_type = (~ma.getmaskarray(variable_data)).sum()
        max_count = max(max_count, n_type)

        # Calculate the fraction for each surface type
        fraction = n_type / n_domain

        # Store both the count and fraction in the results dictionary
        result_types[var] = {"count": n_type, "fraction": fraction}

    # Calculate building statistics.
    if "buildings_2d" in nc_static.variables:
        buildings_2d = nc_static.variables["buildings_2d"][:]
        mean_buildings = buildings_2d.mean()
        max_buildings = buildings_2d.max()
        if mean_buildings is ma.masked or max_buildings is ma.masked:
            mean_buildings = None
            max_buildings = None
    else:
        mean_buildings = None
        max_buildings = None

    # Calculate resolved vegetation statistics.
    TreeStats = TypedDict(
        "TreeStats",
        {
            "non_na_sum": ma.MaskedArray,
            "fraction": ma.MaskedArray,
            "sum": ma.MaskedArray,
            "total_sum": ma.MaskedArray,
            "mean": ma.MaskedArray,
        },
    )
    result_tree: Dict[str, Optional[TreeStats]] = {}
    max_total_sum = 0.0
    if (
        n_x is not None
        and n_y is not None
        and n_zlad is not None
        and vol_grid_cell is not None
        and "zlad" in nc_static.variables
        and ("lad" in nc_static.variables or "bad" in nc_static.variables)
    ):
        zlad = nc_static.variables["zlad"][:]

        # Array to store if vegetation is present. Will be True if either LAD or BAD > 0.0.
        resolved_vegetation_present = np.full((n_zlad, n_y, n_x), False)

        for variable_tree in ["lad", "bad"]:
            if variable_tree in nc_static.variables:
                tree_field = nc_static.variables[variable_tree][:, :]

                if tree_field.shape[0] > 1:
                    # Check if vegetation is present and replace ma.masked with False because
                    # np.logical_or(ma.masked, True) is ma.masked.
                    resolved_vegetation_present = np.logical_or(
                        (tree_field > 0.0).filled(False),
                        resolved_vegetation_present,
                    )

                non_na_sum = (~ma.getmaskarray(tree_field)).sum(axis=(1, 2))
                sum = (tree_field * vol_grid_cell).sum(axis=(1, 2)).filled(0)
                total_sum = sum.sum()
                result_tree[variable_tree] = {
                    "non_na_sum": non_na_sum,
                    "fraction": non_na_sum / n_domain,
                    "sum": sum,
                    "total_sum": total_sum,
                    "mean": tree_field.mean(axis=(1, 2)),
                }
                max_total_sum = max(max_total_sum, total_sum)

            else:
                result_tree[variable_tree] = None

        if "tree_type" in nc_static.variables:
            tree_type = nc_static.variables["tree_type"][:]
            fraction_tree_type = (
                ~ma.getmaskarray(tree_type)[resolved_vegetation_present]
            ).sum() / resolved_vegetation_present.sum()
        else:
            fraction_tree_type = 0.0

        if "tree_id" in nc_static.variables:
            tree_id = nc_static.variables["tree_id"][:]
            fraction_tree_id = (
                ~ma.getmaskarray(tree_id)[resolved_vegetation_present]
            ).sum() / resolved_vegetation_present.sum()
        else:
            fraction_tree_id = 0.0

        # Calculate 2d fraction of LAD/BAD > 0.0.
        n_resolved_vegetation = resolved_vegetation_present.any(axis=0).sum()
        max_count = max(max_count, n_resolved_vegetation)

        # Calculate the fraction for each surface type
        fraction = n_resolved_vegetation / n_domain

        # Store both the count and fraction in the results dictionary
        result_types["resolved_vegetation"] = {"count": n_resolved_vegetation, "fraction": fraction}

    else:
        zlad = None
        fraction_tree_type = None
        fraction_tree_id = None

    # Print results.
    # Grid parameters
    if n_zbuild is not None:
        n_zbuild_str = f", nz Buildings: {n_zbuild}"
    else:
        n_zbuild_str = ""

    if n_zlad is not None:
        n_zlad_str = f", nz LAD: {n_zlad}"
    else:
        n_zlad_str = ""

    if dzbuild is not None:
        dz_str = f", dz: {dzbuild} m"
    elif dzlad is not None:
        dz_str = f", dz: {dzlad} m"
    else:
        dz_str = ""

    if vol_grid_cell is not None:
        vol_str = f", cell volume: {vol_grid_cell} m³"
    else:
        vol_str = ""

    logger.info("Dimensions and grid spacing:")
    logger.info_indent(
        f"(nx+1): {n_x}, (ny+1): {n_y}"
        + n_zbuild_str
        + n_zlad_str
        + "\n"
        + f"dx: {dx} m, dy: {dy} m"
        + dz_str
        + vol_str,
    )

    if mean_buildings is not None and max_buildings is not None:
        logger.info("Buildings:")
        logger.info_indent(
            f"Mean building height: {mean_buildings:.1f} m\n"
            + f"Max. building height: {max_buildings:.1f} m",
        )
    else:
        logger.info("No buildings found.")

    # Fractions
    def fraction_count_str(variable: str) -> str:
        """Helper function to format the fraction and count of a surface type."""
        var = result_types.get(variable, {})
        return (
            f"{var.get('fraction', 0):5.1%} "
            + f"(count: {var.get('count', 0):{len(str(max_count))}})"
        )

    logger.info("Fraction surface types:")
    logger.info_indent(
        f"Vegetation:   {fraction_count_str('vegetation_type')}\n"
        + f"Pavement:     {fraction_count_str('pavement_type')}\n"
        + f"Buildings:    {fraction_count_str('building_type')}\n"
        + f"Water:        {fraction_count_str('water_type')}",
    )

    # LAD and BAD fields
    logger.info("Resolved vegetation:")
    logger.info_indent(f"2D coverage:  {fraction_count_str('resolved_vegetation')}")

    if zlad is not None:
        if fraction_tree_type is not None:
            logger.info_indent(f"3D cells with defined type: {fraction_tree_type:5.1%}")
        if fraction_tree_id is not None:
            logger.info_indent(f"3D cells with defined ID:   {fraction_tree_id:5.1%}")

        len_max_total_sum = len(f"{max_total_sum:.0f}")
        if result_tree["lad"] is not None:
            logger.info_indent(
                f"Total leaf area:  {result_tree['lad']['sum'].sum():{len_max_total_sum}.0f} m²",
            )
        if result_tree["bad"] is not None:
            logger.info_indent(
                f"Total basal area: {result_tree['bad']['sum'].sum():{len_max_total_sum}.0f} m²",
            )

        def tree_variable_table_str(dict_tree: TreeStats) -> List[List[str]]:
            """Helper function to format the table for a tree variable."""
            len_count = len(f"{ma.max(dict_tree['non_na_sum'][1:])}")

            # Initialize the table with headers.
            table = [["coverage (count)"], ["total area"], ["av. density"]]

            # Populate the table columns.
            for i in range(1, len(zlad)):
                table[0].append(
                    f"{dict_tree['fraction'][i]:5.1%} ({dict_tree['non_na_sum'][i]:{len_count}})"
                )
                table[1].append(f"{dict_tree['sum'][i]:.1f} m²")
                if dict_tree["mean"][i] is ma.masked:
                    table[2].append("    --    ")
                else:
                    table[2].append(f"{dict_tree['mean'][i]:.2f} m²/m³")

            return table

        # Initialize the table with the zlad column.
        table = [["      "] + [f"{zlad[i]} m" for i in range(1, len(zlad))]]
        # Add the tree variables to the table.
        if result_tree["lad"] is not None:
            table += tree_variable_table_str(result_tree["lad"])
        if result_tree["bad"] is not None:
            table += tree_variable_table_str(result_tree["bad"])

        # Maximum column width for alignment
        max_lengths = [max(len(item) for item in sublist) for sublist in table]

        # Table header
        header = f"{' Height':^{max_lengths[0]}}"
        i = 1
        if result_tree["lad"] is not None:
            len_columns = max_lengths[i] + max_lengths[i + 1] + max_lengths[i + 2] + 6
            header += f" | {'Leaf area density (LAD)':^{len_columns}}"
            i += 3
        if result_tree["lad"] is not None:
            len_columns = max_lengths[i] + max_lengths[i + 1] + max_lengths[i + 2] + 6
            header += f" | {'Basal area density (BAD)':^{len_columns}}"

        # Calculate the number of rows, assuming all columns have the same number of rows.
        num_rows = len(table[0])

        logger.info_indent("\n" + header)

        for i in range(num_rows):
            # Row left-aligned based on max_lengths.
            row_values = [f"{table[col][i]:>{max_lengths[col]}}" for col in range(len(table))]
            # Join the column values with a separator and print/log the row.
            row_str = " | ".join(row_values)
            logger.info("   " + row_str)
            # Header separator
            if i == 0:
                logger.info_indent("-" * len(row_str))

    else:
        logger.info_indent("No BAD or LAD fields in model domain.")

    if show_plot or plot_file is not None:
        plot_static(
            nc_static,
            show=show_plot,
            detailed=detailed_plot,
            output=Path(plot_file) if plot_file else None,
            title=plot_title,
            width=plot_width,
            height=plot_height,
            geo_referenced=geo_referenced,
        )

    # Close the NetCDF file.
    nc_static.close()
