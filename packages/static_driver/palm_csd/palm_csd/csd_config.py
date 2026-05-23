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

"""Module to handle the palm_csd configuration.

The reading of the input yaml file and storing of the several objects is done by CSDConfig. Several
classes CSDConfig* are defined to handle a specific section of the palm_csd configuration. The
validity of the input values and the existance of the input files is checked by pydantic. The
default values are read from a csv file.
"""

import csv
import logging
import os
from enum import Enum
from importlib.resources import files

# remove Dict here and replace by dict below once Python >=3.9 could be used
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Optional, Set, Tuple, Type, TypeVar, Union, cast

import geopandas as gpd
import rasterio.warp as riowp
from pydantic import (
    AfterValidator,
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    StringConstraints,
    ValidationError,
    ValidationInfo,
    field_validator,
    model_validator,
)
from pydantic.types import PathType
from typing_extensions import Annotated, TypedDict

from palm_csd import StatusLogger
from palm_csd.constants import (
    INPUT_DATA_EXPANDED,
    NBUILDING_SURFACE_LAYER,
    IndexBuildingGeneralPars,
    IndexBuildingIndoorPars,
    IndexBuildingSurfaceLevel,
    IndexBuildingSurfaceType,
    IndexPavementType,
    IndexWaterType,
    InputDataVector,
)
from palm_csd.tools import DefaultMinMax

# module logger, StatusLogger already set in __init__.py so cast only for type checking
logger = cast(StatusLogger, logging.getLogger(__name__))


def _populate_defaults() -> Dict[str, DefaultMinMax]:
    """Read default, minimum and maximum values from value_defaults.csv.

    The csv file data/value_defaults.csv is read and the values are returned as a dictionary.

    Returns:
        Dictionary with the variable names as keys and the DefaultMinMax objects as values.
    """
    defaults = {}
    # Read csv from palm_csd.data.
    with files("palm_csd.data").joinpath("value_defaults.csv").open() as default_min_max_csv:
        reader = csv.DictReader(default_min_max_csv)
        for row in reader:
            # Process each row. Each row is a dict of strings.
            # remove whitespaces
            row = {key.strip(): value.strip() for key, value in row.items()}
            # just the default, minimum and maximum values
            row_data = {key: row[key] for key in ["default", "minimum", "maximum"]}
            # The values in row_data are strings. These are handled by the before validator in
            # the DefaultMinMax class. This breaks static type checking. Explicitly pass it to
            # the model_validate method to get the correct types.
            defaults[row["variable"]] = DefaultMinMax.model_validate(row_data)

    return defaults


value_defaults = _populate_defaults()
"""Default, minimum and maximum value for the variables in the keys."""
# TODO: implement some kind of unpacking for defaults, current issues:
# - mypy does not like it and does not recognize the default value
# - harder to differentiate between no default and None default


class CSDConfigElement(BaseModel):
    """Basic class for configuration elements.

    This class defines a counter to check how often the class has been used.
    """

    # configuration of this and inherited classes
    # frozen=True: no new fields can be modified after creation
    # extra="forbid": extra fields are not allowed
    # validate_default=True: validate default values
    model_config = ConfigDict(frozen=True, extra="forbid", validate_default=True)

    _type: ClassVar[str]
    """Type of the configuration class."""
    counter: ClassVar[int] = 0
    """Number of initialized objects."""
    _unique: ClassVar[bool] = True
    """Indicate whether only one instance of this class is allowed."""

    @model_validator(mode="after")
    def _validate_unique_counter(self: "CSDConfigElement") -> "CSDConfigElement":
        """Handle and validate the number of instances counter.

        This method checks if there already is an instance of the class, when there should only be
        one. Otherwise, the instance counter is increased.

        Args:
            self: The instance.

        Raises:
            ValueError: More than one instance when only one is allowed.

        Returns:
            Unmodified instance.
        """
        if type(self)._unique and type(self).counter == 1:
            raise ValueError(f"More than 1 configuration section of type {type(self)._type} found")

        # increase number of processed configs
        type(self).counter += 1

        return self

    @classmethod
    def _reset_counter(cls) -> None:
        """Reset the counter for the class."""
        cls.counter = 0


class CSDConfigAttributes(CSDConfigElement):
    """Class for global attributes in the static driver."""

    _type = "attributes"

    author: Optional[str] = None
    """Author of the static driver."""
    contact_person: Optional[str] = None
    """Contact person for the static driver."""
    institution: Optional[str] = None
    """Institution where the static driver was made."""
    acronym: Optional[str] = None
    """Institutional acronym."""
    campaign: Optional[str] = None
    """Campaign the static driver belongs to."""
    location: Optional[str] = None
    """Geo-location of the static driver content."""
    site: Optional[str] = None
    """Site description of the static driver content."""
    comment: Optional[str] = None
    """General comment."""
    data_content: Optional[str] = None
    """Data content."""
    dependencies: Optional[str] = None
    """Dependencies."""
    keywords: Optional[str] = None
    """Keywords."""
    references: Optional[str] = None
    """References."""
    source: Optional[str] = None
    """List of data sources used to generate the driver."""
    palm_version: Optional[float] = None
    """REMOVED: PALM version."""
    origin_time: Optional[str] = None
    """Reference point in time."""

    @field_validator("palm_version")
    @classmethod
    def _message_removed_palm_version(
        cls, value: Optional[float], info: ValidationInfo
    ) -> Optional[float]:
        return _validate_removed_conf(
            value,
            info,
        )


def _default_not_none(name: str) -> float:
    """Check if defaults[name].default is not None and return it.

    Args:
        name: The variable name.

    Raises:
        ValueError: The default value of name is None.

    Returns:
        The default value of the variable name.
    """
    value = value_defaults[name].default
    if value is None:
        raise ValueError(f"Default value of {name} must not be None")
    return value


ScalingMethods = TypedDict(
    "ScalingMethods",
    {
        "categorical": riowp.Resampling,
        "continuous": riowp.Resampling,
        "discontinuous": riowp.Resampling,
        "discrete": riowp.Resampling,
    },
)
"""Typed dictionary for the resampling methods."""

_downscaling_method_default: ScalingMethods = {
    "categorical": riowp.Resampling.nearest,
    "continuous": riowp.Resampling.bilinear,
    "discontinuous": riowp.Resampling.nearest,
    "discrete": riowp.Resampling.nearest,
}
"""Default resampling methods for downscaling."""
_upscaling_method_default: ScalingMethods = {
    "categorical": riowp.Resampling.mode,
    "continuous": riowp.Resampling.average,
    "discontinuous": riowp.Resampling.average,
    "discrete": riowp.Resampling.average,
}
"""Default resampling methods for upscaling."""


def _expand_scaling(
    input_value: Optional[
        Union[
            Dict[str, str],
            Dict[str, riowp.Resampling],
            Dict[str, Union[str, riowp.Resampling]],
            str,
            riowp.Resampling,
        ]
    ],
    default: ScalingMethods,
) -> ScalingMethods:
    """Expand scaling input.

    The input_value is expanded to the typed dictionary ScalingMethods. The input_value can be a
    single value or a dictionary. A string key of the dictionary is expanded to the full dictionary
    with a Resampling enum value.

    Args:
        input_value: Input value to expand.
        default: Default values to expand to.

    Raises:
        ValueError: Unknown resampling method or key in input_value.

    Returns:
        Expanded dictionary.
    """
    output_value = default.copy()

    if isinstance(input_value, str):
        for key in output_value.keys():
            try:
                output_value[key] = riowp.Resampling[input_value]  # type: ignore
            # https://github.com/python/mypy/issues/7178
            except KeyError:
                raise ValueError(f"Unknown resampling method {input_value} for {key}.")
    elif isinstance(input_value, riowp.Resampling):
        for key in output_value.keys():
            output_value[key] = input_value  # type: ignore
    elif isinstance(input_value, dict):
        for key, value in input_value.items():
            if key not in output_value.keys():
                raise ValueError(f"Unknown key {key} for resampling method.")
            if isinstance(value, str):
                try:
                    output_value[key] = riowp.Resampling[value]  # type: ignore
                except KeyError:
                    raise ValueError(f"Unknown resampling method {value} for {key}.")
            elif isinstance(value, riowp.Resampling):
                output_value[key] = value  # type: ignore
            else:
                raise ValueError(f"Unknown type {type(value)} for resampling method.")
    else:
        raise ValueError(f"Unknown type {type(input_value)} for resampling method.")

    return output_value


def _validate_scaling(input: ScalingMethods) -> ScalingMethods:
    """Validate scaling input.

    The input is checked if the keys are valid and the values are of type Resampling. For the
    categorical type, only nearest and mode are accepted.

    Args:
        input: Input value to validate.
        enum: Enum class to check against.
        info: ValidationInfo with attribute field_name.

    Raises:
        ValueError: Unknown key or values.

    Returns:
        Validated input value.
    """
    for key, value in input.items():
        if key not in ScalingMethods.__annotations__.keys():
            raise ValueError(f"Unknown key {key} for resampling method.")
        if not isinstance(value, riowp.Resampling):
            raise ValueError(f"Value for key {key} must be a string.")
        if key == "categorical":
            _check_string(value.name, ["nearest", "mode"])
        else:
            _check_string(value.name, [method.name for method in riowp.Resampling])
    return input


class CSDConfigSettings(CSDConfigElement):
    """Class for settings in the static driver."""

    _type = "settings"

    bridge_width: Optional[float] = None
    """REMOVED: Depth of the bridge."""
    ignore_input_georeferencing: bool = False
    """Ignore input file's georeferencing."""

    downscaling_method: Annotated[
        ScalingMethods,
        BeforeValidator(lambda x: _expand_scaling(x, _downscaling_method_default)),
        AfterValidator(_validate_scaling),
    ] = _downscaling_method_default
    """Resampling method for downscaling."""

    upscaling_method: Annotated[
        ScalingMethods,
        BeforeValidator(lambda x: _expand_scaling(x, _upscaling_method_default)),
        AfterValidator(_validate_scaling),
    ] = _upscaling_method_default
    """Resampling method for upscaling."""

    height_high_vegetation_lower_threshold: float = Field(
        default=_default_not_none("height_high_vegetation_lower_threshold"),
        ge=value_defaults["height_high_vegetation_lower_threshold"].minimum,
        le=value_defaults["height_high_vegetation_lower_threshold"].maximum,
    )
    """Lower threshold for high vegetation height."""

    height_rel_resolved_vegetation_lower_threshold: float = Field(
        default=_default_not_none("height_rel_resolved_vegetation_lower_threshold"),
        ge=value_defaults["height_rel_resolved_vegetation_lower_threshold"].minimum,
        le=value_defaults["height_rel_resolved_vegetation_lower_threshold"].maximum,
    )
    """Minimum vegetation height for resolved vegetation relative to dz."""

    lai_roof_extensive: Optional[float] = None
    """REMOVED: LAI of extensive roof greening."""

    lai_roof_intensive: Optional[float] = None
    """REMOVED: LAI of intensive roof greening."""

    lai_tree_lower_threshold: float = Field(
        default=_default_not_none("lai_tree_lower_threshold"),
        ge=value_defaults["lai_tree_lower_threshold"].minimum,
        le=value_defaults["lai_tree_lower_threshold"].maximum,
    )
    """Lower threshold for tree LAI."""

    lai_per_vegetation_height: float = Field(
        default=_default_not_none("lai_per_vegetation_height"),
        ge=value_defaults["lai_per_vegetation_height"].minimum,
        le=value_defaults["lai_per_vegetation_height"].maximum,
    )

    lai_low_vegetation_default: Optional[float] = Field(
        default=value_defaults["lai_low_vegetation_default"].default,
        ge=value_defaults["lai_low_vegetation_default"].minimum,
        le=value_defaults["lai_low_vegetation_default"].maximum,
    )
    """Default LAI for low vegetation."""

    lai_high_vegetation_default: Optional[float] = Field(
        default=value_defaults["lai_high_vegetation_default"].default,
        ge=value_defaults["lai_high_vegetation_default"].minimum,
        le=value_defaults["lai_high_vegetation_default"].maximum,
    )
    """Default LAI for high vegetation."""

    lad_alpha: float = Field(
        default=_default_not_none("lad_alpha"),
        ge=value_defaults["lad_alpha"].minimum,
        le=value_defaults["lad_alpha"].maximum,
    )
    """alpha parameter for LAD profile of Markkanen et al. (2003)."""

    lad_beta: float = Field(
        default=_default_not_none("lad_beta"),
        ge=value_defaults["lad_beta"].minimum,
        le=value_defaults["lad_beta"].maximum,
    )
    """beta parameter for LAD profile of Markkanen et al. (2003)."""

    lad_method: str = "Metal2003"
    """Method for LAD profile calculation."""

    lad_z_max_rel: float = Field(
        default=_default_not_none("lad_z_max_rel"),
        ge=value_defaults["lad_z_max_rel"].minimum,
        le=value_defaults["lad_z_max_rel"].maximum,
    )
    """Relative height of the maximum LAD for LAD profile of Lalic and Mihailovic (2004)."""

    patch_height_default: float = Field(
        default=_default_not_none("patch_height_default"),
        ge=value_defaults["patch_height_default"].minimum,
        le=value_defaults["patch_height_default"].maximum,
    )
    """Default vegetation patch height."""

    replace_invalid_input_values: bool = True
    """Replace invalid input values in input files by the respective default value."""

    rotation_angle: float = Field(
        default=_default_not_none("rotation_angle"),
        ge=value_defaults["rotation_angle"].minimum,
        le=value_defaults["rotation_angle"].maximum,
    )
    """Rotation angle of the domains."""

    season: str = "summer"
    """Season for LAD profile calculation."""

    soil_type_default: int = Field(
        default=int(_default_not_none("soil_type")),
        ge=value_defaults["soil_type"].minimum,
        le=value_defaults["soil_type"].maximum,
    )
    """Default soil type."""

    use_lai_for_roofs: bool = True
    """Use lai for vegetation on roofs."""

    use_lai_for_trees: bool = True
    """Use lai for trees."""

    use_vegetation_height_for_trees: bool = True
    """Use vegetation height for trees."""

    vegetation_type_below_trees: int = Field(
        default=int(_default_not_none("vegetation_type_below_trees")),
        ge=value_defaults["vegetation_type_below_trees"].minimum,
        le=value_defaults["vegetation_type_below_trees"].maximum,
    )
    """Vegetation type below trees."""

    debug: Optional[Union[bool, List[str]]] = None
    """REMOVED: Debug switches."""

    epsg: Optional[int] = None
    """EPSG code for georeferencing."""

    @field_validator("bridge_width")
    @classmethod
    def _message_old_bridge_width(
        cls, value: Optional[float], info: ValidationInfo
    ) -> Optional[float]:
        return _validate_removed_conf(
            value, info, "It is now called bridge_depth and is set in the domain section."
        )

    @field_validator("debug")
    @classmethod
    def _message_former_debug(
        cls, value: Optional[Union[bool, List[str]]], info: ValidationInfo
    ) -> Optional[Union[bool, List[str]]]:
        return _validate_removed_conf(value, info, "Use command line switch -v/--verbose instead.")

    @field_validator("lai_roof_extensive")
    @classmethod
    def _message_old_lai_roof_extensive(
        cls, value: Optional[float], info: ValidationInfo
    ) -> Optional[float]:
        return _validate_removed_conf(
            value, info, "Define LAI with building_lai for surface level roof instead."
        )

    @field_validator("lai_roof_intensive")
    @classmethod
    def _message_old_lai_roof_intensive(
        cls, value: Optional[float], info: ValidationInfo
    ) -> Optional[float]:
        return _validate_removed_conf(
            value, info, "Define LAI with building_lai for surface level roof instead."
        )

    @field_validator("season")
    @classmethod
    def _check_season(cls, value: str) -> str:
        """Check if season is summer or winter.

        Args:
            value: Input season value.

        Returns:
            The valid season value.
        """
        _check_string(value, ["summer", "winter"])
        return value

    @field_validator("lad_method")
    @classmethod
    def _check_lad_method(cls, value: str) -> str:
        """Check if lad_method is Metal2003 or LM2004.

        Args:
            value: Input lad_method value.

        Returns:
            The valid lad_method value.
        """
        _check_string(value, ["Metal2003", "LM2004"])
        return value


def _expand_user_path(path: Optional[Path]) -> Optional[Path]:
    """Expand user path.

    ~ in path is expanded to the user's home directory.

    Args:
        path: Input path.

    Returns:
        Output path with expanded ~.
    """
    return path.expanduser() if path is not None else None


FileLike = TypeVar("FileLike", Path, str)


def _prepend_path_to_file(file: Optional[FileLike], info: ValidationInfo) -> Optional[FileLike]:
    """Prepend path from ValidationInfo to file.

    The path value of info.data is prepended to the file if it exists.

    Args:
        file: Input file.
        info: ValidationInfo with already validated values in data attribute.

    Returns:
        file with prepended path.
    """
    if file is None:
        return None
    if "path" not in info.data:
        return file
    if isinstance(file, str):
        return os.path.join(info.data["path"], file) if info.data["path"] else file
    elif isinstance(file, Path):
        return Path(info.data["path"]) / file if info.data["path"] else file
    else:
        raise ValueError(f"Unknown type {type(file)} for file.")


# ~ expansion in path according to https://github.com/pydantic/pydantic/issues/7990
DirPathExpanded = Annotated[Path, AfterValidator(_expand_user_path), PathType("dir")]
"""Existing directory Path with expanded ~."""

FilePathPrependedExpanded = Annotated[
    Path,
    BeforeValidator(_prepend_path_to_file),
    AfterValidator(_expand_user_path),
    PathType("file"),
]
"""Existing file Path with prepended path attribute and expanded ~."""

PathPrependedExpanded = Annotated[
    Path, BeforeValidator(_prepend_path_to_file), AfterValidator(_expand_user_path)
]
"""Path with prepended path attribute and expanded ~. Can exist or not."""


class CSDConfigOutput(CSDConfigElement):
    """Class for output configuration of the static driver.

    path is prepended to file_out during validation.
    """

    _type = "output"

    path: Optional[DirPathExpanded] = None
    """Path for output files."""
    file_out: PathPrependedExpanded
    """Output file name."""
    version: Optional[int] = None
    """Version of the output file."""

    @field_validator("file_out")
    @classmethod
    def _check_file_out(cls, path: Path) -> Path:
        """Check if path is a file and parent directory exists.

        Args:
            path: file_out Path.

        Raises:
            ValueError: path is a directory.
            ValueError: path's parent directory does not exist.

        Returns:
            Validated path.
        """
        # path must be a file
        if path.is_dir():
            raise ValueError(f"Path {path} is a directory.")
        # parent directory must exist
        if not path.parent.exists():
            raise ValueError(f"Parent directory {path.parent} does not exist.")
        return path


StrLowercaseStrip = Annotated[str, StringConstraints(to_lower=True, strip_whitespace=True)]
IntStrLowercaseStrip = Union[int, StrLowercaseStrip]


class CSDConfigInput(CSDConfigElement):
    """Class for input data configuration for the static driver.

    path is prepended to all file_* Paths during validation. All input files are checked for
    existance..
    """

    _type = "input"
    _unique = False

    pixel_size: Optional[float] = None
    """REMOVED: (Intendend) Pixel size of the input data."""

    path: Optional[DirPathExpanded] = None
    """Path for input files."""

    files: Dict[str, List[FilePathPrependedExpanded]]

    file_buildings_2d: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for building height."""
    file_bridges_2d: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for bridge height."""
    file_building_id: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for building ID."""
    file_bridges_id: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for bridge ID."""
    file_building_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for building type."""
    file_vegetation_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for vegetation type."""
    file_vegetation_height: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for vegetation height."""
    file_pavement_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for pavement type."""
    file_water_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for water type."""
    file_street_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for street type."""
    file_street_crossings: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for street crossings."""
    file_vegetation_on_roofs: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for vegetation on roofs."""
    file_patch_height: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for vegetation patch height."""
    file_x_UTM: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for x UTM coordinates."""
    file_y_UTM: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for y UTM coordinates."""
    file_lat: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for latitude coordinates."""
    file_lon: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for longitude coordinates."""
    file_lai: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for Leaf Area Index."""
    file_lcz: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for Local Climate Zones."""
    file_patch_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for vegetation patch type."""
    file_soil_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for soil type."""
    file_tree_crown_diameter: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for tree crown diameter."""
    file_tree_height: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for tree height."""
    file_tree_lai: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for tree LAI."""
    file_tree_shape: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for tree shape."""
    file_tree_trunk_diameter: Optional[FilePathPrependedExpanded] = None
    """IREMOVED: nput file for tree trunk diameter."""
    file_tree_type: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for tree type."""
    file_water_temperature: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for water temperature."""
    file_zt: Optional[FilePathPrependedExpanded] = None
    """REMOVED: Input file for terrain height."""

    columns: Dict[
        StrLowercaseStrip, Union[StrLowercaseStrip, Dict[IntStrLowercaseStrip, StrLowercaseStrip]]
    ] = {}
    """Assignment of input column names to column types or subtypes."""

    _used_file: Set[Path] = set()
    """List of used input files."""

    @field_validator("file_patch_height")
    @classmethod
    def _message_removed_file_patch_height(
        cls, value: Optional[FilePathPrependedExpanded], info: ValidationInfo
    ) -> Optional[FilePathPrependedExpanded]:
        return _validate_removed_conf(
            value,
            info,
            "Use vegetation_height in files to set (among other things) "
            + "the vegetation patch height.",
        )

    @field_validator(
        "file_bridges_2d",
        "file_bridges_id",
        "file_building_id",
        "file_building_type",
        "file_buildings_2d",
        "file_lat",
        "file_lai",
        "file_lcz",
        "file_lon",
        "file_patch_type",
        "file_pavement_type",
        "file_soil_type",
        "file_street_crossings",
        "file_street_type",
        "file_tree_crown_diameter",
        "file_tree_height",
        "file_tree_lai",
        "file_tree_shape",
        "file_tree_trunk_diameter",
        "file_tree_type",
        "file_vegetation_height",
        "file_vegetation_on_roofs",
        "file_vegetation_type",
        "file_water_temperature",
        "file_water_type",
        "file_x_UTM",
        "file_y_UTM",
        "file_zt",
    )
    @classmethod
    def _message_removed_separate_file_entries(
        cls, value: Optional[FilePathPrependedExpanded], info: ValidationInfo
    ) -> Optional[FilePathPrependedExpanded]:
        return _validate_removed_conf(
            value,
            info,
        )

    @field_validator("pixel_size")
    @classmethod
    def _message_removed_pixel_size(
        cls, value: Optional[float], info: ValidationInfo
    ) -> Optional[float]:
        return _validate_removed_conf(
            value, info, "Use input in each domain section to define its input section."
        )

    @field_validator("files", mode="before")
    @classmethod
    def _file_elements_to_list(cls, value: Any) -> Optional[Dict]:
        """Wrap single value in a list if it is not already a list.

        Args:
            value: Value to check.

        Returns:
            Dictionary with list values if value is not None, otherwise None.
        """
        if value is None:
            return None
        if not isinstance(value, dict):
            raise ValueError(f"Expected a dictionary for files, got {type(value)} instead.")
        for key, val in value.items():
            if not isinstance(val, list):
                value[key] = [val]
        return value

    @field_validator("columns", mode="after")
    @classmethod
    def _columns_keys_valid(
        cls,
        value: Dict[
            StrLowercaseStrip,
            Union[StrLowercaseStrip, Dict[IntStrLowercaseStrip, StrLowercaseStrip]],
        ],
    ) -> Dict[
        StrLowercaseStrip, Union[StrLowercaseStrip, Dict[IntStrLowercaseStrip, StrLowercaseStrip]]
    ]:
        valid_variables = set(INPUT_DATA_EXPANDED["name"]).union({col for col in InputDataVector})
        invalid = []
        for key in value.values():
            if isinstance(key, str):
                matches = {k for k in valid_variables if k.startswith(key)}
                if len(matches) == 0:
                    invalid.append(key)
        if len(invalid) > 0:
            raise ValueError(
                f"Invalid input column{'s' if len(invalid) > 1 else ''} found: {invalid}"
            )
        return value

    @field_validator("files", mode="after")
    @classmethod
    def _files_keys_valid(
        cls, value: Dict[str, List[FilePathPrependedExpanded]]
    ) -> Dict[str, List[FilePathPrependedExpanded]]:
        valid_variables = set(INPUT_DATA_EXPANDED["name"]).union({col for col in InputDataVector})
        invalid = []
        for key in value.keys():
            matches = {k for k in valid_variables if k.startswith(key)}
            if len(matches) == 0:
                invalid.append(key)
        if len(invalid) > 0:
            raise ValueError(
                f"Invalid input file{'s' if len(invalid) > 1 else ''} found: {invalid}"
            )
        # Check if value lists have only one element for all entries except those in InputDataVector
        for key, file_list in value.items():
            if key not in {col.value for col in InputDataVector} and len(file_list) > 1:
                raise ValueError(
                    f"Input file '{key}' should be a single file, but got {len(file_list)} given."
                )
        return value

    @model_validator(mode="after")
    def _check_columns_exist(self) -> "CSDConfigInput":
        """Check if all columns exists in vector data."""
        if self.columns == {}:
            return self
        input_columns: Set[str] = set()
        for input in InputDataVector:
            if input.value in self.files:
                for vinput in self.files[input]:
                    gdf = gpd.read_file(vinput, rows=1)
                    input_columns = input_columns.union({col.lower() for col in gdf.columns})

        for column in self.columns.keys():
            if column not in input_columns:
                raise ValueError(
                    f"Column '{column}' not found in input vector data. "
                    + "Please check the input files and the columns configuration."
                )
        return self

    def any_netcdf(self) -> bool:
        """Check if any input file is a netCDF file.

        Returns:
            True if any input file is a netCDF file, False otherwise.
        """
        for paths in self.files.values():
            for p in paths:
                if p.suffix == ".nc":
                    return True
        return False

    def add_used_file(self, file: Path) -> None:
        """Add file to list of used files.

        Args:
            file: File to add.
        """
        self._used_file.add(file)

    def unused_file(self) -> List[Tuple[str, Path]]:
        """Return list of unread files.

        Returns:
            List of unread files.
        """
        unused_file_path = []
        for file, paths in self.files.items():
            for p in paths:
                if p not in self._used_file:
                    unused_file_path.append((file, p))

        return unused_file_path


FloatInt = Union[float, int]


def _expand_parslike(
    input_value: Optional[
        Union[
            Dict[str, Union[FloatInt, List[FloatInt]]],
            Dict[int, Union[FloatInt, List[FloatInt]]],
            Dict[Union[str, int], Union[FloatInt, List[FloatInt]]],
            Union[FloatInt, List[FloatInt]],
        ]
    ],
    enum: Type[Enum],
    nlayer: int = 1,
) -> Optional[Dict[int, Union[FloatInt, List[FloatInt]]]]:
    """Expand parslike input with attribute indices according to enum and nlayer.

    The input_value is expanded to a dictionary with indices according to the enum and the number of
    layers. The input_value can be a single value, a list of values or a dictionary. A string key of
    the dictionary is expanded to the full attribute name and then translated to the integer index.
    The value is expanded to the number of layers, if necessary.

    Example:
        _expand_parslike(1.8, IndexBuildingSurfaceType, nlayer=4) returns
        {
            0: [1.8, 1.8, 1.8, 1.8],
            1: [1.8, 1.8, 1.8, 1.8],
            2: [1.8, 1.8, 1.8, 1.8],
            3: [1.8, 1.8, 1.8, 1.8],
            4: [1.8, 1.8, 1.8, 1.8],
            5: [1.8, 1.8, 1.8, 1.8],
            6: [1.8, 1.8, 1.8, 1.8],
            7: [1.8, 1.8, 1.8, 1.8],
            8: [1.8, 1.8, 1.8, 1.8],
        }

        _expand_parslike({"wall": 1.8}, IndexBuildingSurfaceType, nlayer=4) returns
        {
            0: [1.8, 1.8, 1.8, 1.8],
            1: [1.8, 1.8, 1.8, 1.8],
            2: [1.8, 1.8, 1.8, 1.8],
        }

    Args:
        input_value: Input value to expand.
        enum: Enum class to expand to.
        nlayer: Number of layers to expand to. Defaults to 1.

    Raises:
        ValueError: nlayer must be at least 1.
        ValueError: Length of value list in input_value does not match number of layer.
        ValueError: Key in input_value cannot be expanded to any key in enum.
        ValueError: Key in input_value is not a valid key in enum.
        ValueError: Unknown type in input_value.

    Returns:
        Expanded dictionary.
    """
    if input_value is None:
        return None

    if nlayer < 1:
        raise ValueError("Layer must be at least 1.")

    enum_key_strings = [element.name for element in enum]
    enum_key_values = [element.value for element in enum]

    output_value: Union[FloatInt, List[FloatInt]]
    if isinstance(input_value, (int, float, list)):
        if isinstance(input_value, list):
            if len(input_value) == nlayer:
                output_value = input_value
            elif len(input_value) == 1:
                output_value = input_value * nlayer
            else:
                raise ValueError(
                    f"Length of value list {len(input_value)} does not match "
                    + f"number of layer {nlayer}."
                )
        else:
            if nlayer == 1:
                output_value = input_value
            else:
                output_value = [input_value] * nlayer
        return {key: output_value for key in enum_key_values}

    if not isinstance(input_value, dict):
        raise ValueError(f"Unknown type {type(input_value)} for parslike value.")

    dict_value = {}
    for key, value in input_value.items():
        if isinstance(value, (int, float)):
            if nlayer == 1:
                output_value = value
            else:
                output_value = [value] * nlayer
        elif isinstance(value, list):
            if len(value) == nlayer:
                output_value = value
            elif len(value) == 1:
                output_value = value * nlayer
            else:
                raise ValueError(
                    f"Length of value list {len(value)} does not match number of layer {nlayer}."
                )
        else:
            raise ValueError(f"Unknown type {type(value)} for parslike value.")

        if isinstance(key, str):
            keys_starting_with_key = [k for k in enum_key_strings if k.startswith(key)]
            if not keys_starting_with_key:
                raise ValueError(f"Key {key} cannot be expanded to any key in {enum_key_strings}.")
            for k in keys_starting_with_key:
                dict_value[enum[k].value] = output_value
        elif isinstance(key, int):
            if key not in enum_key_values:
                raise ValueError(f"Key {key} is not included in the valid keys {enum_key_values}.")
            dict_value[key] = output_value
        else:
            raise ValueError(f"Unknown type {type(key)} for parslike key.")

    return dict_value


def _validate_parslike(
    input: Optional[Dict[int, Union[FloatInt, List[FloatInt]]]],
    enum: Type[Enum],
    info: ValidationInfo,
) -> Optional[Dict[int, Union[FloatInt, List[FloatInt]]]]:
    """Validate parslike input named info.field_name with data from defaults.

    Args:
        input: Input value to validate.
        enum: Enum class to check against.
        info: ValidationInfo with attribute field_name.

    Raises:
        ValueError: Field name in info is None.
        ValueError: Key in input does not match any key in defaults.

    Returns:
        Validated input value.
    """
    if input is None:
        return input

    if info.field_name is None:
        raise ValueError("Field name is None. Should not happen.")

    for key, value in input.items():
        full_key = f"{info.field_name}_{enum(key).name}"
        keys_starting_with_key = [k for k in value_defaults.keys() if full_key.startswith(k)]
        if not keys_starting_with_key:
            raise ValueError(
                f"Key {full_key} cannot be expanded to any key in {list(value_defaults.keys())}."
            )
        if len(keys_starting_with_key) > 1:
            # take the longest match as a list
            keys_starting_with_key = sorted(keys_starting_with_key, key=len, reverse=True)[0:1]
        default = value_defaults[keys_starting_with_key[0]]
        if isinstance(value, list):
            for v in value:
                _check_within_range(v, default)
        else:
            _check_within_range(value, default)
    return input


DictBuildingGeneralFloat = Annotated[
    Dict[int, float],
    BeforeValidator(lambda x: _expand_parslike(x, IndexBuildingGeneralPars)),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexBuildingGeneralPars, info)),
]
"""Dictionary with building general parameters as float values."""

DictBuildingIndoorFloat = Annotated[
    Dict[int, float],
    BeforeValidator(lambda x: _expand_parslike(x, IndexBuildingIndoorPars)),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexBuildingIndoorPars, info)),
]
"""Dictionary with building indoor parameters as float values."""

DictBuildingSurfaceLevelFloat = Annotated[
    Dict[int, float],
    BeforeValidator(lambda x: _expand_parslike(x, IndexBuildingSurfaceLevel)),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexBuildingSurfaceLevel, info)),
]
"""Dictionary with building surface level parameters as float values."""

DictBuildingSurfaceTypeInt = Annotated[
    Dict[int, int],
    BeforeValidator(lambda x: _expand_parslike(x, IndexBuildingSurfaceType)),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexBuildingSurfaceType, info)),
]
"""Dictionary with building surface type parameters as integer values."""

DictBuildingSurfaceTypeFloat = Annotated[
    Dict[int, float],
    BeforeValidator(lambda x: _expand_parslike(x, IndexBuildingSurfaceType)),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexBuildingSurfaceType, info)),
]
"""Dictionary with building surface type parameters as float values."""

DictBuildingSurfaceTypeFloatLayer = Annotated[
    Dict[int, List[float]],
    BeforeValidator(
        lambda x: _expand_parslike(x, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER)
    ),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexBuildingSurfaceType, info)),
]
"""Dictionary with building surface type parameters as list of float values, one for each layer."""

DictWaterTypeFloat = Annotated[
    Dict[int, float],
    BeforeValidator(lambda x: _expand_parslike(x, IndexWaterType)),
    AfterValidator(lambda x, info: _validate_parslike(x, IndexWaterType, info)),
]
"""Dictionary with water type parameters as float values."""


class CSDConfigDomain(CSDConfigElement):
    """Class for domain configuration of the static driver.

    All parslike values are expanded to the full attribute names according to the Index* enums.
    """

    _type = "domain"
    _unique = False

    domain_parent: Optional[str] = None
    """Name of the parent domain."""
    input: Optional[str] = None
    """Name of the input section for the domain."""

    pixel_size: float = Field(
        ge=value_defaults["pixel_size"].minimum,
        le=value_defaults["pixel_size"].maximum,
    )
    """Grid spacing in x and y direction."""

    nx: int = Field(
        ge=value_defaults["nx"].minimum,
        le=value_defaults["nx"].maximum,
    )
    """Number of grid points -1 in x direction."""

    ny: int = Field(
        ge=value_defaults["ny"].minimum,
        le=value_defaults["ny"].maximum,
    )
    """Number of grid points -1 in y direction."""

    dz: float = Field(
        ge=value_defaults["dz"].minimum,
        le=value_defaults["dz"].maximum,
    )
    """Grid spacing in z direction."""

    input_lower_left_x: Optional[float] = Field(
        default=value_defaults["input_lower_left_x"].default,
        ge=value_defaults["input_lower_left_x"].minimum,
        le=value_defaults["input_lower_left_x"].maximum,
    )
    """Distance (m) on x-axis between the lower-left domain corner and the lower-left input data
    corner."""

    input_lower_left_y: Optional[float] = Field(
        default=value_defaults["input_lower_left_y"].default,
        ge=value_defaults["input_lower_left_y"].minimum,
        le=value_defaults["input_lower_left_y"].maximum,
    )
    """Distance (m) on y-axis between the lower-left domain corner and the lower-left input data
    corner."""

    lower_left_x: Optional[float] = Field(
        default=value_defaults["lower_left_x"].default,
        ge=value_defaults["lower_left_x"].minimum,
        le=value_defaults["lower_left_x"].maximum,
    )
    """Distance (m) on x-axis between the lower-left corner of the nested domain and the lower-left
    corner of the root parent domain."""

    lower_left_y: Optional[float] = Field(
        default=value_defaults["lower_left_y"].default,
        ge=value_defaults["lower_left_y"].minimum,
        le=value_defaults["lower_left_y"].maximum,
    )
    """Distance (m) on y-axis between the lower-left corner of the nested domain and the lower-left
    corner of the root parent domain."""

    origin_lon: Optional[float] = Field(
        default=value_defaults["origin_lon"].default,
        ge=value_defaults["origin_lon"].minimum,
        le=value_defaults["origin_lon"].maximum,
    )
    """Longitude of the left border of the lower-left grid point of the PALM domain in WGS84."""

    origin_lat: Optional[float] = Field(
        default=value_defaults["origin_lat"].default,
        ge=value_defaults["origin_lat"].minimum,
        le=value_defaults["origin_lat"].maximum,
    )
    """Latitude of the lower border of the lower-left grid point of the PALM domain in WGS84."""

    origin_x: Optional[float] = None
    """x-coordinate of the left border of the lower-left grid point of the PALM domain in the
    custom CRS."""
    origin_y: Optional[float] = None
    """y-coordinate of the lower border of the lower-left grid point of the PALM domain in the
    custom CRS."""

    bridge_depth: float = Field(
        default=_default_not_none("bridge_depth"),
        ge=value_defaults["bridge_depth"].minimum,
        le=value_defaults["bridge_depth"].maximum,
    )
    """Structural depth of bridge pixels."""

    buildings_3d: Optional[bool] = None
    """REMOVED: Use generate 3D buildings field."""

    building_free_border_width: float = Field(
        default=_default_not_none("building_free_border_width"),
        ge=value_defaults["building_free_border_width"].minimum,
        le=value_defaults["building_free_border_width"].maximum,
    )
    """Width of the border (in meters) around the domain that is free of buildings."""

    building_free_border_pavement_type: int = Field(
        default=int(_default_not_none("building_free_border_pavement_type")),
        ge=value_defaults["building_free_border_pavement_type"].minimum,
        le=value_defaults["building_free_border_pavement_type"].maximum,
    )
    """Substitute pavement type where buildings are removed."""

    generate_buildings_3d: bool = False
    """Generate 3D buildings field."""

    allow_high_vegetation: Optional[bool] = None
    """REMOVED: Allow high vegetation in vegetation_type."""
    estimate_lai_from_vegetation_height: bool = True
    """Calculate resolved LAI from vegetation height."""
    generate_single_trees: bool = True
    """Generate LAD and BAD fields of single trees."""
    generate_vegetation_on_roofs: Optional[bool] = None
    """REMOVED: Generate vegetation on roofs."""
    generate_vegetation_patches: bool = True
    """Generate vegetation patches."""
    overhanging_trees: bool = True
    """Allow overhanging trees."""
    remove_low_lai_tree: bool = False
    """Remove trees with low Lead Area Index."""
    replace_high_vegetation_types: bool = True
    """Replace high vegetation in vegetation_type."""
    street_trees: Optional[bool] = None
    """REMOVED: Generate LAD and BAD fields of single trees."""
    vegetation_on_roofs: Optional[bool] = None
    """REMOVED: Generate vegetation on roofs."""

    interpolate_terrain: bool = False
    """Interpolate and blend terrain height between parent and child domains."""
    use_palm_z_axis: bool = False
    """Raster terrain height on the z-grid of PALM."""

    lcz_input: Optional[str] = None
    """Use Local Climate Zone input."""
    dcep: bool = False
    """Generate DCEP fields."""
    z_uhl: List[float] = Field(default=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0])
    """Height of the DCEP urban half layers."""
    udir: List[float] = Field(default=[0.0, 90.0])
    """Street orientations for DCEP."""

    building_albedo_type: Optional[DictBuildingSurfaceTypeInt] = None
    """Albedo type of building surfaces."""
    building_emissivity: Optional[DictBuildingSurfaceTypeFloat] = None
    """Emissivity of building surfaces."""
    building_fraction: Optional[DictBuildingSurfaceTypeFloat] = None
    """Fraction of building surfaces."""
    building_general_pars: Optional[DictBuildingGeneralFloat] = None
    """General parameters of buildings."""
    building_heat_capacity: Optional[DictBuildingSurfaceTypeFloatLayer] = None
    """Heat capacity of building surfaces."""
    building_heat_conductivity: Optional[DictBuildingSurfaceTypeFloatLayer] = None
    """Heat conductivity of building surfaces."""
    building_indoor_pars: Optional[DictBuildingIndoorFloat] = None
    """Indoor parameters of buildings."""
    building_lai: Optional[DictBuildingSurfaceLevelFloat] = None
    """Leaf Area Index of building surfaces."""
    building_roughness_length: Optional[DictBuildingSurfaceLevelFloat] = None
    """Roughness length of building surfaces."""
    building_roughness_length_qh: Optional[DictBuildingSurfaceLevelFloat] = None
    """Roughness length for moisture and heat of building surfaces."""
    building_thickness: Optional[DictBuildingSurfaceTypeFloatLayer] = None
    """Thickness of layers of building surfaces."""
    building_transmissivity: Optional[DictBuildingSurfaceLevelFloat] = None
    """Transmissivity of building windows."""

    water_temperature: Optional[DictWaterTypeFloat] = None
    """Water temperature for water types."""
    water_temperature_per_water_type: Optional[DictWaterTypeFloat] = None
    """REMOVED: Water temperature for water types."""

    @field_validator("building_free_border_pavement_type", mode="before")
    @classmethod
    def _expand_building_free_border_pavement_type(cls, values: Any) -> Any:
        if isinstance(values, str):
            enum_key_strings = [element.name for element in IndexPavementType]
            matching_keys = [k for k in enum_key_strings if k == values]
            if not matching_keys:
                raise ValueError(
                    f"Value {values} cannot be expanded to any key in {enum_key_strings}."
                )
            if len(matching_keys) > 1:
                raise ValueError(
                    f"Value {values} is ambiguous and matches multiple keys in "
                    + f"{enum_key_strings}."
                )
            return IndexPavementType[matching_keys[0]].value
        return values

    @field_validator("z_uhl")
    @classmethod
    def _check_z_uhl(cls, values: List[float]) -> List[float]:
        """Check if z_uhl is monotonously increasing and within range.

        Args:
            values: Input z_uhl values.

        Raises:
            ValueError: values[0] is not 0.0.
            ValueError: values are not monotonously increasing.

        Returns:
            Validated z_uhl values.
        """
        if values[0] != 0.0:
            raise ValueError("z_uhl[0] must be 0.0")
        for i in range(1, len(values)):
            _check_within_range(values[i], value_defaults["z_uhl"])
            if (values[i] - values[i - 1]) <= 0.0:
                raise ValueError(
                    f"z_uhl not monotonously increasing from {values[i - 1]} " + f"to {values[i]}."
                )
        return values

    @field_validator("udir")
    @classmethod
    def _check_udir(cls, values: List[float]) -> List[float]:
        """Check if udir values are within range.

        Args:
            values: Input udir values.

        Returns:
            Validated udir values.
        """
        for value in values:
            _check_within_range(value, value_defaults["udir"])
        return values

    @field_validator("allow_high_vegetation")
    @classmethod
    def _message_old_allow_high_vegetation(
        cls, value: Optional[bool], info: ValidationInfo
    ) -> Optional[bool]:
        return _validate_removed_conf(
            value,
            info,
            "It is replaced by replace_high_vegetation_types with in inverted meaning.",
        )

    @field_validator("generate_vegetation_on_roofs")
    @classmethod
    def _message_old_generate_vegetation_on_roofs(
        cls, value: Optional[bool], info: ValidationInfo
    ) -> Optional[bool]:
        return _validate_removed_conf(
            value,
            info,
            "Define building_fraction for surface type wall_roof, window_roof and green_roof, "
            + "or building_lai for surface level roof instead.",
        )

    @field_validator("buildings_3d")
    @classmethod
    def _message_old_buildings_3d(
        cls, value: Optional[bool], info: ValidationInfo
    ) -> Optional[bool]:
        return _validate_removed_conf(
            value,
            info,
            "It is replaced by generate_buildings_3d.",
        )

    @field_validator("street_trees")
    @classmethod
    def _message_old_street_trees(
        cls, value: Optional[bool], info: ValidationInfo
    ) -> Optional[bool]:
        return _validate_removed_conf(
            value,
            info,
            "It is now called generate_single_trees.",
        )

    @field_validator("vegetation_on_roofs")
    @classmethod
    def _message_old_vegetation_on_roofs(
        cls, value: Optional[bool], info: ValidationInfo
    ) -> Optional[bool]:
        return _validate_removed_conf(
            value,
            info,
            "Define building_fraction for surface type wall_roof, window_roof and green_roof, "
            + "or building_lai for surface level roof instead.",
        )

    @field_validator("water_temperature_per_water_type")
    @classmethod
    def _message_old_water_temperature_per_water_type(
        cls, values: Optional[Dict[int, float]], info: ValidationInfo
    ) -> Optional[Dict[int, float]]:
        return _validate_removed_conf(
            values,
            info,
            "It is now called water_temperature and "
            + "supports named water types and single floats for all water types.",
        )

    @model_validator(mode="after")
    def _validate_interpolation_palm_z_axis(self: "CSDConfigDomain") -> "CSDConfigDomain":
        """Check if interpolate_terrain and use_palm_z_axis are consistent.

        Args:
            self: Instance of CSDConfigDomain.

        Raises:
            ValueError: interpolate_terrain set to True and use_palm_z_axis set to False.

        Returns:
            Validated instance of CSDConfigDomain.
        """
        if self.interpolate_terrain and not self.use_palm_z_axis:
            raise ValueError("use_palm_z_axis needs to be true for interpolate_terrain")
        return self


class CSDConfigLCZ(CSDConfigElement):
    """Class for domain configuration of the static driver."""

    _type = "lcz"

    height_geometric_mean: bool = True
    """Use geometric mean of building heights for LCZ classification."""

    compact_highrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Compact highrise LCZ class."""
    compact_midrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Compact midrise LCZ class."""
    compact_lowrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Compact lowrise LCZ class."""
    open_highrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Open highrise LCZ class."""
    open_midrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Open midrise LCZ class."""
    open_lowrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Open lowrise LCZ class."""
    lightweight_lowrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Lightweight lowrise LCZ class."""
    large_lowrise: Optional[Dict[str, float]] = None
    """Updated parameters of the Large lowrise LCZ class."""
    sparsely_built: Optional[Dict[str, float]] = None
    """Updated parameters of the Sparsely built LCZ class."""
    heavy_industry: Optional[Dict[str, float]] = None
    """Updated parameters of the Heavy industry LCZ class."""
    dense_trees: Optional[Dict[str, float]] = None
    """Updated parameters of the Dense trees LCZ class."""
    scattered_trees: Optional[Dict[str, float]] = None
    """Updated parameters of the Scattered trees LCZ class."""
    bush_scrub: Optional[Dict[str, float]] = None
    """Updated parameters of the Bush scrub LCZ class."""
    low_plants: Optional[Dict[str, float]] = None
    """Updated parameters of the Low plants LCZ class."""
    bare_rock_or_paved: Optional[Dict[str, float]] = None
    """Updated parameters of the Bare rock or paved LCZ class."""
    bare_soil_or_sand: Optional[Dict[str, float]] = None
    """Updated parameters of the Bare soil or sand LCZ class."""
    water: Optional[Dict[str, float]] = None
    """Updated parameters of the Water LCZ class."""


class CSDConfig:
    """Class to collect all palm_csd configuration sections."""

    attributes: CSDConfigAttributes
    """attributes section of the configuration."""
    settings: CSDConfigSettings
    """settings section of the configuration."""
    output: CSDConfigOutput
    """output section of the configuration."""
    input_dict: Dict[str, CSDConfigInput]
    """input sections of the configuration."""
    domain_dict: Dict[str, CSDConfigDomain]
    """domain sections of the configuration."""
    lcz: CSDConfigLCZ
    """lcz section of the configuration."""

    def __init__(self, configuration_dict: Dict) -> None:
        """Create the configuration sections from nested dictionary data.

        Args:
            configuration_dict: Nested dictionary data.

        Raises:
            ValueError: The configuration section names include too many separators.
            ValueError: The configuration section name is unknown.
            ValidationError: Errors in the configuration section.
        """
        self.domain_dict = {}
        self.input_dict = {}
        counter_input = 0
        for section_key, section_value in configuration_dict.items():
            key_splitted = section_key.split("_")
            if key_splitted[0] == "input" or key_splitted[0] == "domain":
                if len(key_splitted) > 2:
                    logger.critical_raise(
                        f"The section {section_key} in the configuration file "
                        + "includes too many separators '_'.",
                    )
            else:
                if len(key_splitted) > 1:
                    logger.critical_raise(
                        f"The section {section_key} in the configuration file "
                        + "includes separators '_'.\n"
                        + "Only the sections 'input' and 'domain' are allowed to have separators.",
                    )
            try:
                # use match/case for the following once Python >=3.10 could be used
                if key_splitted[0] == "attributes":
                    self.attributes = CSDConfigAttributes(**section_value)
                elif key_splitted[0] == "settings":
                    self.settings = CSDConfigSettings(**section_value)
                elif key_splitted[0] == "output":
                    self.output = CSDConfigOutput(**section_value)
                elif key_splitted[0] == "input":
                    if len(key_splitted) == 1:
                        name = f"input_{counter_input}"
                        counter_input += 1
                        logger.debug(f"Assigning name {name} to unnamed input.")
                    else:
                        name = key_splitted[1]
                    if name in self.input_dict:
                        logger.critical_raise(f"Input {name} is defined more than once.")
                    self.input_dict[name] = CSDConfigInput(**section_value)
                elif key_splitted[0] == "domain":
                    if len(key_splitted) == 1:
                        name = "root"
                        logger.info(f"Assigning name {name} to unnamed domain.")
                    else:
                        name = key_splitted[1]
                    if name in self.domain_dict:
                        logger.critical_raise(f"Domain {name} is defined more than once.")
                    self.domain_dict[name] = CSDConfigDomain(**section_value)
                elif key_splitted[0] == "lcz":
                    self.lcz = CSDConfigLCZ(**section_value)
                else:
                    logger.critical_raise(
                        f"Unknown section {key_splitted[0]} in the configuration file."
                    )
            except ValidationError as errors:
                n_error = len(errors.errors())
                logger.critical(
                    f"In configuration section {section_key}, "
                    + f"{n_error} {'error' if n_error == 1 else 'errors'} occurred."
                )
                for error in errors.errors():
                    if error["type"] == "extra_forbidden":
                        logger.critical(f"Unknown key: {error['loc'][0]}")
                    elif error["type"] == "greater_than_equal":
                        logger.critical(f"{error['loc'][0]}: {error['msg']}.")
                        logger.critical(
                            f"If the value {error['input']} should be used, "
                            + "adjust the minimum in palm_csd/data/value_defaults.csv."
                        )
                    elif error["type"] == "less_than_equal":
                        logger.critical(f"{error['loc'][0]}: {error['msg']}.")
                        logger.critical(
                            f"If the value {error['input']} should be used, "
                            + "adjust the maximum in palm_csd/data/value_defaults.csv."
                        )
                    else:
                        message = f"{error['loc'][0]}: " if error["loc"] else ""
                        if "ctx" in error:
                            logger.critical(f"{message}{error['ctx']['error']}")
                        else:
                            logger.critical(f"{message}{error['msg']}")
                raise

        # if not LCZ section is given, use default values
        if CSDConfigLCZ.counter == 0:
            self.lcz = CSDConfigLCZ()

        self._update_value_defaults()

    def _update_value_defaults(self) -> None:
        """Update value defaults with input settings."""
        value_defaults["soil_type"].default = self.settings.soil_type_default

    def input_of_domain(self, domain_name: str) -> CSDConfigInput:
        """Find fitting input configuration for a domain.

        The domain_name and the input option are used to identify the input. If there is only one
        input section, this is used.

        Args:
            domain_name: Domain name

        Raises:
            ValueError: No fitting input configuration section was found.

        Returns:
            Fitting Input configuration.
        """
        domain_config = self.domain_dict[domain_name]

        # Check if either domain name or a given input in domain configuration fits with input
        # section names.
        if domain_config.input is not None:
            input_name = domain_config.input
        else:
            input_name = domain_name
        try:
            return self.input_dict[input_name]
        except KeyError:
            if domain_config.input is not None:
                logger.critical_raise(
                    f"Input section {domain_config.input} of domain {domain_name} not found.",
                )

        # If there is only one input section, use it.
        if len(self.input_dict) == 1:
            return list(self.input_dict.values())[0]

        # If no fitting input section was found, raise an error.
        logger.critical_raise(
            f"For the domain {domain_name} with pixel_size {domain_config.pixel_size}, "
            + "no fitting input configuration section was found.",
        )

    def input_of_parent_domain(self, domain_config: CSDConfigDomain) -> Optional[CSDConfigInput]:
        """Find fitting input configuration for the parent of a domain.

        Args:
            domain_config: Domain configuration of which to use the parent.

        Returns:
            Input configuration of the parent domain configuration if a parent exists. If not,
            return None.
        """
        domain_parent = domain_config.domain_parent
        if domain_parent is None:
            return None
        return self.input_of_domain(domain_parent)


def reset_all_config_counters() -> None:
    """Reset counters of all config classes."""
    CSDConfigAttributes._reset_counter()
    CSDConfigDomain._reset_counter()
    CSDConfigInput._reset_counter()
    CSDConfigLCZ._reset_counter()
    CSDConfigOutput._reset_counter()
    CSDConfigSettings._reset_counter()


T = TypeVar("T")


def _validate_removed_conf(
    value: T, validation_info: ValidationInfo, user_info: Optional[str] = None
) -> T:
    """Validate that value is None.

    Args:
        value: Value to validate.
        validation_info: ValidationInfo with attribute field_name.
        user_info: Additional information for the user. Defaults to None.

    Raises:
        ValueError: value is not None.

    Returns:
        Input value.
    """
    if value is not None:
        message = f"{validation_info.field_name} is not used anymore."
        if user_info is not None:
            message += "\n" + user_info
        raise ValueError(message)
    return value


def _validate_deprecated_conf(
    value: T, validation_info: ValidationInfo, user_info: Optional[str] = None
) -> T:
    """Issue warning that the configuration option is deprecated.

    Args:
        value: Value to validate.
        validation_info: ValidationInfo with attribute field_name.
        user_info: Additional information for the user. Defaults to None.

    Raises:
        ValueError: value is not None.

    Returns:
        Input value.
    """
    if value is not None:
        message = f"{validation_info.field_name} is deprecated."
        if user_info is not None:
            message += "\n" + user_info
        logger.warning(message)
    return value


def _check_string(value: str, valid: List[str]) -> None:
    """Check if value is in valid. Raise ValueError if not.

    Args:
        value: Value to check.
        valid: List of valid values.

    Raises:
        ValueError: value not in valid list.
    """
    if value not in valid:
        raise ValueError(f"Unknown value {value}. Valid values are {' '.join(valid)}.")


def _check_within_range(value: float, range: DefaultMinMax) -> None:
    """Check if value is within the minimum and maximum of range. Raise ValueError if not.

    Args:
        value: Value to check.
        range: Reference minimum und maximum.

    Raises:
        ValueError: value is smaller than minimum or larger than maximum.
    """
    if range.minimum is not None and value < range.minimum:
        raise ValueError(
            f"Value {value} is smaller than minimum {range.minimum}. type=greater_than_equal"
        )
    if range.maximum is not None and value > range.maximum:
        raise ValueError(
            f"Value {value} is larger than maximum {range.maximum}. type=less_than_equal"
        )
