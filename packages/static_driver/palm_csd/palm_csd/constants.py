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

"""Module of constants and Enums."""

from enum import Enum, IntEnum, auto
from typing import Dict, List, Literal, NamedTuple, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import pandas.api.typing as pdtyping

NGHOST_POINTS = 3
"""Number of ghost points in PALM."""

# number of wall layers fixed in PALM
NBUILDING_SURFACE_LAYER = 4
"""Number of layers in building surfaces."""


# TODO: Use native StrEnum when using Python 3.11
class StrEnum(str, Enum):
    """StrEnum where enum.auto() returns the lower-cased member name."""

    @staticmethod
    def _generate_next_value_(name: str, start: int, count: int, last_values: list) -> str:
        return name.lower()

    def __str__(self) -> str:
        """Return the value of the enum."""
        return self.value


class ColumnGeneral(StrEnum):
    """General column names."""

    geometry = auto()
    """Coordinates in a geopandas GeoDataFrame"""
    x_index = auto()
    """x-index of the tree"""
    y_index = auto()
    """y-index of the tree"""


class InputData(StrEnum):
    """Names of input data, either single raster files or columns in InputDataVector."""

    bridges_2d = auto()
    """height of the bridge"""
    bridges_id = auto()
    """bridge id"""
    building_albedo_type = auto()
    """building albedo type"""
    building_emissivity = auto()
    """building emissivity"""
    building_fraction = auto()
    """building fraction"""
    building_general_pars = auto()
    """building general parameters"""
    building_heat_capacity = auto()
    """building heat capacity"""
    building_heat_conductivity = auto()
    """building heat conductivity"""
    building_id = auto()
    """building id"""
    building_indoor_pars = auto()
    """building indoor parameters"""
    building_lai = auto()
    """building leaf area index"""
    building_roughness_length = auto()
    """building roughness length"""
    building_roughness_length_qh = auto()
    """building roughness length (qh)"""
    building_thickness = auto()
    """building thickness"""
    building_transmissivity = auto()
    """building transmissivity"""
    building_type = auto()
    """building type"""
    buildings_2d = auto()
    """height of the building"""
    lat = auto()
    """latitude"""
    lai = auto()
    """leaf area index"""
    lcz = auto()
    """local climate zone"""
    lon = auto()
    """longitude"""
    patch_type = auto()
    """patch type"""
    pavement_type = auto()
    """type of the pavement"""
    soil_type = auto()
    """soil type"""
    street_crossings = auto()
    """street crossings"""
    street_type = auto()
    """street type"""
    tree_crown_diameter = auto()
    """crown diameter of the tree"""
    tree_height = auto()
    """height of the tree"""
    tree_lai = auto()
    """leaf area index of the tree"""
    tree_shape = auto()
    """shape of the tree"""
    tree_trunk_diameter = auto()
    """trunk diameter of the tree"""
    tree_type = auto()
    """type of the tree"""
    tree_type_name = auto()
    """type name of the tree"""
    vegetation_height = auto()
    """height of the vegetation"""
    vegetation_type = auto()
    """vegetation type"""
    water_temperature = auto()
    """water temperature"""
    water_type = auto()
    """water type"""
    x_utm = auto()
    """x-coordinate in UTM"""
    y_utm = auto()
    """y-coordinate in UTM"""
    zt = auto()
    """height of the surface above the ground"""


def get_parent_input_data(variable: str) -> InputData:
    """Get the corresponding InputData enum member for a given variable name.

    For example, "building_heat_capacity_wall_gfl_1" will return InputData.building_heat_capacity.

    Args:
        variable: The variable name to match against InputData enum members.

    Raises:
        ValueError: If the variable does not match any InputData member.
        ValueError: If the variable is ambiguous and can be expanded to multiple keys.

    Returns:
        The corresponding InputData enum member.
    """
    keys_starting_with_key = [k for k in InputData if variable.startswith(k)]
    if not keys_starting_with_key:
        raise ValueError(f"Key {variable} cannot be expanded to any key in {list(InputData)}.")
    if len(keys_starting_with_key) > 1:
        raise ValueError(
            f"Key {variable} is ambiguous and can be expanded to {keys_starting_with_key}."
        )
    return keys_starting_with_key[0]


class InputDataVector(StrEnum):
    """Names of input data."""

    surfaces = auto()
    """surface vector input"""
    trees = auto()
    """trees vector input"""


class InputDataInfo(NamedTuple):
    """Information how to process input data."""

    type: Literal["continuous", "discontinuous", "categorical", "discrete"]
    initialize_default: bool = False
    all_or_none_missing: bool = False
    warning_point_data: bool = False


INPUT_DATA_INFO = {
    InputData.bridges_2d: InputDataInfo("discontinuous"),
    InputData.bridges_id: InputDataInfo("categorical"),
    InputData.building_albedo_type: InputDataInfo("categorical"),
    InputData.building_emissivity: InputDataInfo("discontinuous"),
    InputData.building_fraction: InputDataInfo("discontinuous"),
    InputData.building_general_pars: InputDataInfo("discontinuous"),
    InputData.building_heat_capacity: InputDataInfo("discontinuous"),
    InputData.building_heat_conductivity: InputDataInfo("discontinuous"),
    InputData.building_id: InputDataInfo("categorical"),
    InputData.building_indoor_pars: InputDataInfo("discontinuous"),
    InputData.building_lai: InputDataInfo("discontinuous"),
    InputData.building_roughness_length: InputDataInfo("discontinuous"),
    InputData.building_roughness_length_qh: InputDataInfo("discontinuous"),
    InputData.building_thickness: InputDataInfo("discontinuous"),
    InputData.building_transmissivity: InputDataInfo("discontinuous"),
    InputData.building_type: InputDataInfo("categorical", initialize_default=True),
    InputData.buildings_2d: InputDataInfo("discontinuous"),
    InputData.lai: InputDataInfo("discontinuous"),
    InputData.lat: InputDataInfo("continuous"),
    InputData.lcz: InputDataInfo("categorical", initialize_default=True),
    InputData.lon: InputDataInfo("continuous"),
    InputData.patch_type: InputDataInfo("categorical"),
    InputData.pavement_type: InputDataInfo("categorical"),
    InputData.soil_type: InputDataInfo("categorical", initialize_default=True),
    InputData.street_crossings: InputDataInfo("categorical"),
    InputData.street_type: InputDataInfo("categorical"),
    InputData.tree_crown_diameter: InputDataInfo("discrete", warning_point_data=True),
    InputData.tree_height: InputDataInfo("discrete", warning_point_data=True),
    InputData.tree_lai: InputDataInfo("discrete", warning_point_data=True),
    InputData.tree_shape: InputDataInfo("discrete", warning_point_data=True),
    InputData.tree_trunk_diameter: InputDataInfo("discrete", warning_point_data=True),
    InputData.tree_type: InputDataInfo("categorical", warning_point_data=True),
    InputData.vegetation_height: InputDataInfo("discontinuous"),
    InputData.vegetation_type: InputDataInfo("categorical"),
    InputData.water_temperature: InputDataInfo("continuous"),
    InputData.water_type: InputDataInfo("categorical"),
    InputData.x_utm: InputDataInfo("continuous"),
    InputData.y_utm: InputDataInfo("continuous"),
    InputData.zt: InputDataInfo("continuous", all_or_none_missing=True, initialize_default=True),
}
"""Input data information."""


class IndexBuildingSurfaceLevel(IntEnum):
    """Index for building surface levels."""

    gfl = 0
    """ground floor level"""
    agfl = 1
    """above ground floor level"""
    roof = 2
    """roof"""


class IndexBuildingGeneralPars(IntEnum):
    """Index for building general parameters."""

    height_gfl = 0
    """ground floor level height"""
    green_type_roof = 1
    """type of green roof"""


class IndexBuildingSurfaceType(IntEnum):
    """Index for building surface types."""

    wall_gfl = 0
    """wall ground floor level"""
    wall_agfl = 1
    """wall above ground floor level"""
    wall_roof = 2
    """wall roof"""
    window_gfl = 3
    """window ground floor level"""
    window_agfl = 4
    """window above ground floor level"""
    window_roof = 5
    """window roof"""
    green_gfl = 6
    """green on wall ground floor level"""
    green_agfl = 7
    """green on wall above ground floor level"""
    green_roof = 8
    """green on roof ground floor level"""


class IndexBuildingIndoorPars(IntEnum):
    """Index for building indoor parameters."""

    indoor_temperature_summer = 0
    """indoor target summer temperature"""
    indoor_temperature_winter = 1
    """indoor target winter temperature"""
    shading_window = 2
    """shading factor"""
    g_window = 3
    """g-value windows"""
    u_window = 4
    """u-value windows"""
    airflow_unoccupied = 5
    """basic airflow without occupancy of the room"""
    airflow_occupied = 6
    """additional airflow dependent on occupancy of the room"""
    heat_recovery_efficiency = 7
    """heat recovery efficiency"""
    effective_surface = 8
    """dynamic parameter specific effective surface"""
    inner_heat_storage = 9
    """dynamic parameter innner heat storage"""
    ratio_surface_floor = 10
    """ratio internal surface/floor area"""
    heating_capacity_max = 11
    """maximal heating capacity"""
    cooling_capacity_max = 12
    """maximal cooling capacity"""
    heat_gain_high = 13
    """additional internal heat gains dependent on occupancy of the room"""
    heat_gain_low = 14
    """basic internal heat gains without occupancy of the room"""
    height_storey = 15
    """storey height"""
    height_ceiling_construction = 16
    """ceiling construction height"""
    heating_factor = 17
    """anthropogenic heat output factor for heating"""
    cooling_factor = 18
    """anthropogenic heat output factor for cooling"""


class IndexPavementType(IntEnum):
    """Pavement type index."""

    asphalt_concrete_mix = 1
    """asphalt concrete mix"""
    asphalt = 2
    """asphalt"""
    concrete = 3
    """concrete"""
    sett = 4
    """sett"""
    paving_stones = 5
    """paving stones"""
    cobblestone = 6
    """cobblestone"""
    metal = 7
    """metal"""
    wood = 8
    """wood"""
    gravel = 9
    """gravel"""
    fine_gravel = 10
    """fine gravel"""
    pebblestone = 11
    """pebblestone"""
    woodchips = 12
    """woodchips"""
    tartan = 13
    """tartan"""
    artifical_turf = 14
    """artifical turf"""
    clay = 15
    """clay"""


class IndexStreetType(IntEnum):
    """Index for street types."""

    unclassified = 1
    """unclassified"""
    cycleway = 2
    """cycleway"""
    footway_pedestrian = 3
    """footway / pedestrian"""
    path = 4
    """path"""
    track = 5
    """track"""
    living_street = 6
    """living street"""
    service = 7
    """service"""
    residential = 8
    """residential"""
    tertiary = 9
    """tertiary"""
    tertiary_link = 10
    """tertiary link"""
    secondary = 11
    """secondary"""
    secondary_link = 12
    """secondary link"""
    primary = 13
    """primary"""
    primary_link = 14
    """primary link"""
    trunk = 15
    """trunk"""
    trunk_link = 16
    """trunk link"""
    motorway = 17
    """motorway"""
    motorway_link = 18
    """motorway link"""
    raceway = 19
    """raceway"""


class IndexVegetationType(IntEnum):
    """Vegetation type index."""

    user_defined = 0
    """user defined"""
    bare_soil = 1
    """bare soil"""
    crops_mixed_farming = 2
    """crops, mixed farming"""
    short_grass = 3
    """short grass"""
    evergreen_needleleaf_trees = 4
    """evergreen needleleaf trees"""
    deciduous_needleleaf_trees = 5
    """deciduous needleleaf trees"""
    evergreen_broadleaf_trees = 6
    """evergreen broadleaf trees"""
    deciduous_broadleaf_trees = 7
    """deciduous broadleaf trees"""
    tall_grass = 8
    """tall grass"""
    desert = 9
    """desert"""
    tundra = 10
    """tundra"""
    irrigated_crops = 11
    """irrigated crops"""
    semidesert = 12
    """semidesert"""
    ice_caps_glaciers = 13
    """ice caps and glaciers"""
    bogs_marshes = 14
    """bogs and marshes"""
    evergreen_shrubs = 15
    """evergreen shrubs"""
    deciduous_shrubs = 16
    """deciduous shrubs"""
    mixed_forest_woodland = 17
    """mixed forest/woodland"""
    interrupted_forest = 18
    """interrupted forest"""


VT_HIGH_VEGETATION = [
    IndexVegetationType.evergreen_needleleaf_trees,
    IndexVegetationType.deciduous_needleleaf_trees,
    IndexVegetationType.evergreen_broadleaf_trees,
    IndexVegetationType.deciduous_broadleaf_trees,
    IndexVegetationType.mixed_forest_woodland,
    IndexVegetationType.interrupted_forest,
]
"""Vegetation types that are considered high vegetation"""

VT_NO_PLANTS = [
    IndexVegetationType.bare_soil,
    IndexVegetationType.desert,
    IndexVegetationType.ice_caps_glaciers,
]
"""Vegetation types that are assumed to be without plants."""


class IndexVegetationPars(IntEnum):
    """Vegetation parameters index."""

    canopy_resistance_min = 0
    """Minimum canopy resistance"""
    lai = 1
    """Leaf area index"""
    vegetation_coverage = 2
    """Vegetation coverage"""
    canopy_resistance_coefficient = 3
    """Canopy resistance coefficient"""
    roughness_length = 4
    """Roughness length for momentum"""
    roughness_length_qh = 5
    """Roughness length for heat and moisture"""
    heat_conductivity_stable = 6
    """Skin layer heat conductivity (stable conditions)"""
    heat_conductivity_unstable = 7
    """Skin layer heat conductivity (unstable conditions)"""
    fraction_shortwave_soil = 8
    """Fraction of incoming shortwave radiation transmitted directly to the soil (not implemented
    yet)"""
    heat_capacity = 9
    """Heat capacity of the surface"""
    albedo_type = 10
    """Albedo type"""
    emissivity = 11
    """Surface emissivity """


class IndexBuildingType(IntEnum):
    """Index for building types."""

    residential_1950 = 1
    """Residential, before 1950"""
    residential_1951_2000 = 2
    """Residential, 1951 -- 2000"""
    residential_2001 = 3
    """Residential, after 2001"""
    office_1950 = 4
    """Office, before 1950"""
    office_1951_2000 = 5
    """Office, 1951 -- 2000"""
    office_2001 = 6
    """Office, after 2001"""
    bridges = 7
    """Bridges"""


class IndexSoilType(IntEnum):
    """Index for soil types."""

    coarse = 1
    """Coarse"""
    medium = 2
    """Medium"""
    medium_fine = 3
    """Medium-fine"""
    fine = 4
    """Fine"""
    very_fine = 5
    """Very fine"""
    organic = 6
    """Organic"""


class IndexWaterPars(IntEnum):
    """Index for water parameters."""

    water_temperature = 0
    """water temperature"""
    roughness_length = 1
    """roughness length for momentum"""
    roughness_length_qh = 2
    """roughness length for heat"""
    heat_conductivity_stable = 3
    """heat conductivity between skin layer and water (stable conditions)"""
    heat_conductivity_unstable = 4
    """heat conductivity between skin layer and water (unstable conditions)"""
    albedo_type = 5
    """albedo type"""
    emissivity = 6
    """surface emissivity"""


class IndexWaterType(IntEnum):
    """Index for water types."""

    lake = 1
    """Lake"""
    river = 2
    """River"""
    ocean = 3
    """Ocean"""
    pond = 4
    """Pond"""
    fountain = 5
    """Fountain"""


def _expand_input_data(variable: str) -> pd.DataFrame:
    """Expand input data variable into a DataFrame with names and indices."""
    records: List[Dict[str, Union[str, int, pdtyping.NAType]]] = []
    if variable in ["building_heat_capacity", "building_heat_conductivity", "building_thickness"]:
        for surface_type in IndexBuildingSurfaceType:
            for layer in range(1, NBUILDING_SURFACE_LAYER + 1):
                records.append(
                    {
                        "name": f"{variable}_{surface_type.name}_{layer}",
                        "level": surface_type.value,
                        "layer": layer,
                    }
                )
    elif variable in [
        "building_albedo_type",
        "building_emissivity",
        "building_fraction",
    ]:
        for surface_type in IndexBuildingSurfaceType:
            records.append(
                {
                    "name": f"{variable}_{surface_type.name}",
                    "level": surface_type.value,
                    "layer": pd.NA,
                }
            )
    elif variable in [
        "building_lai",
        "building_roughness_length",
        "building_roughness_length_qh",
        "building_transmissivity",
    ]:
        for surface_level in IndexBuildingSurfaceLevel:
            records.append(
                {
                    "name": f"{variable}_{surface_level.name}",
                    "level": surface_level.value,
                    "layer": pd.NA,
                }
            )
    elif variable == "building_general_pars":
        for general_par in IndexBuildingGeneralPars:
            records.append(
                {
                    "name": f"{variable}_{general_par.name}",
                    "level": general_par.value,
                    "layer": pd.NA,
                }
            )
    elif variable == "building_indoor_pars":
        for indoor_par in IndexBuildingIndoorPars:
            records.append(
                {
                    "name": f"{variable}_{indoor_par.name}",
                    "level": indoor_par.value,
                    "layer": pd.NA,
                }
            )
    else:
        # Default: just return the variable as a single-row DataFrame
        records = [{"name": variable, "level": pd.NA, "layer": pd.NA}]
    return pd.DataFrame(records).astype({"name": "string", "level": "Int8", "layer": "Int8"})


INPUT_DATA_EXPANDED = pd.concat([_expand_input_data(name) for name in InputData], ignore_index=True)


class ValuePlot(IntEnum):
    """Base values for plotting."""

    none = 0

    water_type = 100
    vegetation_resolved = 200
    vegetation_type = 300
    pavement_type = 400
    building_type = 500


# Colors separated by categories as hashing is based on the int of IntEnum.
COLORS: Dict[str, Dict] = {
    "base": {
        ValuePlot.none: "#000000",
        ValuePlot.vegetation_resolved: "#44a03e",
        ValuePlot.building_type: "#d26666",
        ValuePlot.pavement_type: "#bfbfbf",
        ValuePlot.vegetation_type: "#b2df8a",
        ValuePlot.water_type: "#418bbd",
    },
    # Water types - variations of "#418bbd" (blue)
    "water": {
        IndexWaterType.lake: "#5DA8E5",
        IndexWaterType.river: "#1A5E8E",
        IndexWaterType.ocean: "#083C68",
        IndexWaterType.pond: "#90CFFF",
        IndexWaterType.fountain: "#98E3FF",
    },
    "vegetation": {
        # Vegetation types - variations of "#b2df8a" (green)
        IndexVegetationType.user_defined: "#7AB890",
        IndexVegetationType.short_grass: "#9DC183",
        IndexVegetationType.tall_grass: "#7AB055",
        IndexVegetationType.evergreen_needleleaf_trees: "#2D5C3E",
        IndexVegetationType.deciduous_needleleaf_trees: "#4B7F41",
        IndexVegetationType.evergreen_broadleaf_trees: "#1E3F29",
        IndexVegetationType.deciduous_broadleaf_trees: "#5D8C44",
        IndexVegetationType.mixed_forest_woodland: "#3D7A5F",
        IndexVegetationType.interrupted_forest: "#598C69",
        IndexVegetationType.evergreen_shrubs: "#8ABE6D",
        IndexVegetationType.deciduous_shrubs: "#A1D082",
        # Soil/desert-like vegetation types - earth tones
        IndexVegetationType.bare_soil: "#C4A484",
        IndexVegetationType.desert: "#D2B48C",
        IndexVegetationType.tundra: "#A9B277",
        IndexVegetationType.semidesert: "#BFB06C",
        IndexVegetationType.crops_mixed_farming: "#C2D275",
        IndexVegetationType.irrigated_crops: "#A4C639",
        IndexVegetationType.bogs_marshes: "#5F9970",
        IndexVegetationType.ice_caps_glaciers: "#E0E8F0",
    },
    # Pavement types
    "pavement": {
        IndexPavementType.asphalt_concrete_mix: "#505050",
        IndexPavementType.asphalt: "#666666",
        IndexPavementType.concrete: "#DDDDDD",
        IndexPavementType.sett: "#9A9A9A",
        IndexPavementType.paving_stones: "#ABABAB",
        IndexPavementType.cobblestone: "#808080",
        IndexPavementType.metal: "#A0A5C3",
        IndexPavementType.wood: "#BC9B7A",
        IndexPavementType.gravel: "#777777",
        IndexPavementType.fine_gravel: "#B5B5B5",
        IndexPavementType.pebblestone: "#707070",
        IndexPavementType.woodchips: "#C4A682",
        IndexPavementType.tartan: "#B67F7F",
        IndexPavementType.artifical_turf: "#8AAB8A",
        IndexPavementType.clay: "#C49878",
    },
    # Building types - variations of "#d26666" (red)
    "building": {
        IndexBuildingType.residential_1950: "#8A2727",
        IndexBuildingType.residential_1951_2000: "#C92D2D",
        IndexBuildingType.residential_2001: "#F09090",
        IndexBuildingType.office_1950: "#6A2A2A",
        IndexBuildingType.office_1951_2000: "#9A5050",
        IndexBuildingType.office_2001: "#C07878",
        IndexBuildingType.bridges: "#B33E80",
    },
}
"""Colors for different building and surface types."""


class FillValue:
    """Fill values for different data types."""

    BYTE = -127
    FLOAT = -9999.0
    INTEGER = -9999

    @classmethod
    def from_dtype(cls, dtype: npt.DTypeLike) -> float:
        """Return the fill value for a given data type."""
        if dtype == np.int32:
            return cls.INTEGER
        if dtype == np.float64:
            return cls.FLOAT
        if dtype == np.int8:
            return cls.BYTE
        raise ValueError(f"Unknown dtype {dtype}")


def plot_color_overview(figsize=(12, 10)):
    """Plot an overview of all color schemes defined in COLORS.

    Args:
        figsize: Figure size as (width, height) in inches. Defaults to (12, 10).

    Returns:
        The figure object containing the color overview
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # Create figure with subplots - one for each category
    fig = plt.figure(figsize=figsize)
    categories = list(COLORS.keys())
    n_categories = len(categories)

    # Determine layout based on number of categories
    n_rows = (n_categories + 1) // 2

    # Create subplots
    for i, category in enumerate(categories):
        ax = plt.subplot(n_rows, 2, i + 1)
        ax.set_title(f"{category.capitalize()} Colors", fontweight="bold")
        ax.axis("off")

        colors = COLORS[category]

        # Calculate how many rows we need for this category
        items_per_row = 1 if category == "base" else 3

        col_spacing = 0.33

        for j, (key, color) in enumerate(colors.items()):
            row = j // items_per_row
            col = j % items_per_row

            # Position and size of color rectangle
            rect = Rectangle(
                (col * col_spacing, 1 - (row + 1) * 0.15),
                0.25,
                0.1,
                facecolor=color,
                edgecolor="black",
            )
            ax.add_patch(rect)

            # Add key name and color code
            key_name = key.name if isinstance(key, Enum) else f"ValuePlot.{key.name}"
            ax.text(
                col * col_spacing + 0.27,
                1 - (row + 1) * 0.15 + 0.05,
                f"{key_name}\n{color}",
                va="center",
            )

    plt.tight_layout()
    plt.show()

    return fig


if __name__ == "__main__":
    plot_color_overview()
