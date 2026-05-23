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

"""Tools to handle Local Climate Zones (LCZs)."""

from importlib.resources import files
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import numpy.ma as ma
import scipy.stats
from pydantic import BaseModel, Field, model_validator

from palm_csd.csd_config import CSDConfigLCZ, value_defaults
from palm_csd.tools import DefaultMinMax


class LCZ(BaseModel, arbitrary_types_allowed=True, frozen=True):
    """Class to store the properties of a local climate zone (LCZ).

    This class all three fractions (building, impervious, pervious) although their sum is 1. Stewart
    and Oke (2012) define minimum and maximum values for each fraction independently. In order to
    fully use these conditions, separate fractions are implemented.
    """

    aspect_ratio: DefaultMinMax
    """Aspect ratio."""
    building_plan_area_fraction: DefaultMinMax
    """Building plan area fraction."""
    impervious_plan_area_fraction: DefaultMinMax
    """Impervious plan area fraction."""
    pervious_plan_area_fraction: DefaultMinMax
    """Pervious plan area fraction."""
    height_roughness_elements: DefaultMinMax
    """Height of roughness elements."""

    index: int = Field(
        ge=value_defaults["lcz"].minimum,
        le=value_defaults["lcz"].maximum,
    )
    """Index."""

    r: int = Field(
        ge=value_defaults["color_r"].minimum,
        le=value_defaults["color_r"].maximum,
    )
    """Red value in color coding."""
    g: int = Field(
        ge=value_defaults["color_g"].minimum,
        le=value_defaults["color_g"].maximum,
    )
    """Green value in color coding."""
    b: int = Field(
        ge=value_defaults["color_b"].minimum,
        le=value_defaults["color_b"].maximum,
    )
    """Blue value in color coding."""

    vegetation_type: ma.MaskedArray
    """Assigned vegetation type."""
    water_type: ma.MaskedArray
    """Assigned water type."""
    lai: ma.MaskedArray
    """Assigned Leaf Area Index."""

    @model_validator(mode="after")
    def _check_consistency(self: "LCZ") -> "LCZ":
        """Check if the fractions sum up to 1.

        Args:
            self: Instance of the LCZ class.

        Raises:
            ValueError: Fractions are undefined.
            ValueError: Fractions do not sum up to 1.

        Returns:
            Validated instance of the LCZ class.
        """
        if (
            self.building_plan_area_fraction.default is None
            or self.impervious_plan_area_fraction.default is None
            or self.pervious_plan_area_fraction.default is None
        ):
            raise ValueError("Fractions undefined.")
        fractions_sum = (
            self.building_plan_area_fraction.default
            + self.impervious_plan_area_fraction.default
            + self.pervious_plan_area_fraction.default
        )
        if abs(fractions_sum - 1.0) > 1.0e-13:
            raise ValueError(f"Fractions do not sum up to 1 but to {fractions_sum}.")
        return self

    def set_fractions(
        self,
        building_plan_area_fraction: Optional[float] = None,
        impervious_plan_area_fraction: Optional[float] = None,
        pervious_plan_area_fraction: Optional[float] = None,
    ) -> None:
        """Set the fractions ensuring that they sum up to 1.

        Args:
            building_plan_area_fraction: Building plan area fraction. Defaults to None.
            impervious_plan_area_fraction: Impervious plan area fraction. Defaults to None.
            pervious_plan_area_fraction: Pervious plan area fraction. Defaults to None.

        Raises:
            ValueError: Not enough values given to set the fractions.
            ValueError: Fractions do not sum up to 1.
        """
        if building_plan_area_fraction is not None:
            self.building_plan_area_fraction.default = building_plan_area_fraction
            residual = 1.0 - building_plan_area_fraction
            if impervious_plan_area_fraction is not None:
                self.impervious_plan_area_fraction.default = impervious_plan_area_fraction
                residual -= impervious_plan_area_fraction
                if pervious_plan_area_fraction is not None:
                    self.pervious_plan_area_fraction.default = pervious_plan_area_fraction
                    self._check_consistency
                else:
                    self.pervious_plan_area_fraction.default = residual
            else:
                if pervious_plan_area_fraction is not None:
                    self.pervious_plan_area_fraction.default = pervious_plan_area_fraction
                    residual -= pervious_plan_area_fraction
                    self.impervious_plan_area_fraction.default = residual
                else:
                    if self.impervious_plan_area_fraction.default is None:
                        raise ValueError("imperious_plan_area_fraction is None.")
                    self.pervious_plan_area_fraction.default = (
                        residual - self.impervious_plan_area_fraction.default
                    )
        else:
            if impervious_plan_area_fraction is not None:
                self.impervious_plan_area_fraction.default = impervious_plan_area_fraction
                residual = 1.0 - impervious_plan_area_fraction
                if pervious_plan_area_fraction is not None:
                    self.pervious_plan_area_fraction.default = pervious_plan_area_fraction
                    residual -= pervious_plan_area_fraction
                    self.building_plan_area_fraction.default = residual
                else:
                    if self.building_plan_area_fraction.default is None:
                        raise ValueError("building_plan_area_fraction is None.")
                    self.pervious_plan_area_fraction.default = (
                        residual - self.building_plan_area_fraction.default
                    )
            else:
                if pervious_plan_area_fraction is not None:
                    raise ValueError(
                        "Only modifying pervious_plan_area_fraction."
                        + "Do not know how to adjust the other."
                    )
                # else case: all three arguments are None so do nothing


class LCZTypes:
    """Class to store the different LCZ types and converter functions.

    The user has the choice to use the arithmetic or the geometric mean when calculating the
    building height distribution. Here, we use the fact that the geometric mean is the exponential
    of the arithmetic mean of the logarithm of the values.
    """

    compact_highrise: LCZ
    """Compact high-rise LCZ."""
    compact_midrise: LCZ
    """Compact mid-rise LCZ."""
    compact_lowrise: LCZ
    """Compact low-rise LCZ."""
    open_highrise: LCZ
    """Open high-rise LCZ."""
    open_midrise: LCZ
    """Open mid-rise LCZ."""
    open_lowrise: LCZ
    """Open low-rise LCZ."""
    lightweight_lowrise: LCZ
    """Lightweight low-rise LCZ."""
    large_lowrise: LCZ
    """Large low-rise LCZ."""
    sparsely_built: LCZ
    """Sparsely built LCZ."""
    heavy_industry: LCZ
    """Heavy industry LCZ."""
    dense_trees: LCZ
    """Dense trees LCZ."""
    scattered_trees: LCZ
    """Scattered trees LCZ."""
    bush_scrub: LCZ
    """Bush scrub LCZ."""
    low_plants: LCZ
    """Low plants LCZ."""
    bare_rock_or_paved: LCZ
    """Bare rock or paved LCZ."""
    bare_soil_or_sand: LCZ
    """Bare soil or sand LCZ."""
    water: LCZ
    """Water LCZ."""

    vegetation_like: List[LCZ]
    """Vegetation-like LCZs."""
    urban_like: List[LCZ]
    """Urban-like LCZs."""

    index: Dict[int, LCZ]  # List would be faster but starts at 0
    """Index to LCZ mapping."""

    _mean_f1: Callable
    """Helper function applied before mean calculation."""
    _mean_f2: Callable
    """Helper function applied after mean calculation."""
    _mean_height: Callable
    """Helper function to calculate the mean height."""

    def __init__(self, season: str, height_geometric_mean: bool) -> None:
        """Initialize the LCZTypes class.

        The LCZs are defined in the `lcz_definitions.csv` file and the mapping of LCZs to vegetation
        and water types is set in the `lcz_mappings.csv` file.

        Args:
            season: Season for which the LAI is calculated.
            height_geometric_mean: If True, the geometric mean is used for the building height. If
              False, the arithmetic mean is used.

        Raises:
            ValueError: If the season is not 'summer' or 'winter'.
            ValueError: LCZ definitions are inconsistent.
        """
        filling_values_temp = -9999
        if season == "summer":
            lai_label = "lai_summer"
        elif season == "winter":
            lai_label = "lai_winter"
        else:
            raise ValueError(f"Season must either be 'summer' or 'winter' instead of {season}.")

        if height_geometric_mean:
            # The geometric mean is the exponential of the arithmetic mean of the logarithm of the
            # values. Steward and Oke (2012) use geometric mean.
            self._mean_f1 = np.log
            self._mean_f2 = np.exp
        else:
            # Arithmetic mean.
            self._mean_f1 = lambda x: x
            self._mean_f2 = lambda x: x
        self._mean_height = lambda x: self._mean_f2(np.mean(self._mean_f1(x)))

        # Read the LCZ definition:
        #   min and max values based on Steward and Oke (2012)
        #   RGB values based on WUDAPT convention
        #   default values based on W2W convention
        with files("palm_csd.data").joinpath("lcz_definitions.csv").open() as lcz_csv:
            lcz_definitions = np.genfromtxt(
                lcz_csv,
                delimiter=",",
                missing_values="None",
                dtype=None,
                # names=False,
                skip_header=1,
                encoding="utf-8",
            )
        # Mapping of LCZ to vegetation and water properties.
        with files("palm_csd.data").joinpath("lcz_mappings.csv").open() as lcz_csv:
            lcz_mappings = np.genfromtxt(
                lcz_csv,
                delimiter=",",
                missing_values="None",
                filling_values=filling_values_temp,
                dtype=None,
                names=True,
                skip_header=0,
                encoding="utf-8",
            )

        self.index = {}
        for lcz_definition, lcz_mapping in zip(lcz_definitions, lcz_mappings):
            lcz_name = lcz_definition[0]
            # Default height is the mean of the min and max height.
            height_minimum = lcz_definition[14]
            height_default = lcz_definition[15]
            height_maximum = lcz_definition[16]
            if np.isnan(height_default):
                if np.isnan(height_minimum) or np.isnan(height_maximum):
                    raise ValueError(
                        "lcz_definitions.csv: if height_roughness_elements default is"
                        + " None, the respective min and max values need to be given."
                    )
                height_default = self._mean_height([height_minimum, height_maximum])
            lcz = LCZ(
                index=lcz_definition[1],
                aspect_ratio=DefaultMinMax(
                    minimum=lcz_definition[2], default=lcz_definition[3], maximum=lcz_definition[4]
                ),
                building_plan_area_fraction=DefaultMinMax(
                    minimum=lcz_definition[5], default=lcz_definition[6], maximum=lcz_definition[7]
                ),
                impervious_plan_area_fraction=DefaultMinMax(
                    minimum=lcz_definition[8], default=lcz_definition[9], maximum=lcz_definition[10]
                ),
                pervious_plan_area_fraction=DefaultMinMax(
                    minimum=lcz_definition[11],
                    default=lcz_definition[12],
                    maximum=lcz_definition[13],
                ),
                height_roughness_elements=DefaultMinMax(
                    minimum=height_minimum,
                    default=height_default,
                    maximum=height_maximum,
                ),
                r=lcz_definition[17],
                g=lcz_definition[18],
                b=lcz_definition[19],
                vegetation_type=ma.masked_equal(
                    lcz_mapping["vegetation_type"], filling_values_temp
                ),
                water_type=ma.masked_equal(lcz_mapping["water_type"], filling_values_temp),
                lai=ma.masked_equal(lcz_mapping[lai_label], filling_values_temp),
            )
            setattr(self, lcz_name, lcz)
            self.index[lcz.index] = lcz

        self.vegetation_like = [
            self.dense_trees,
            self.scattered_trees,
            self.bush_scrub,
            self.low_plants,
            self.bare_rock_or_paved,
            self.bare_soil_or_sand,
        ]

        self.urban_like = [
            self.compact_highrise,
            self.compact_midrise,
            self.compact_lowrise,
            self.open_highrise,
            self.open_midrise,
            self.open_lowrise,
            self.lightweight_lowrise,
            self.large_lowrise,
            self.sparsely_built,
            self.heavy_industry,
        ]

    def update_defaults(self, lcz_config: CSDConfigLCZ) -> None:
        """Update the default values of the LCZs from the configuration.

        Args:
            lcz_config: User configuration of the LCZs.
        """
        for lcz_name, lcz_properties in vars(lcz_config).items():
            if lcz_properties is not None and isinstance(lcz_properties, dict):
                # get all entries including "fraction" and set all fractions at the same time
                fraction_properties = {k: v for (k, v) in lcz_properties.items() if "fraction" in k}
                lcz = getattr(self, lcz_name)
                lcz.set_fractions(**fraction_properties)
                # apply the other properties directly
                other_properties = {
                    k: v for (k, v) in lcz_properties.items() if "fraction" not in k
                }
                for prop, value in other_properties.items():
                    prop_object = getattr(lcz, prop)
                    if isinstance(prop_object, DefaultMinMax):
                        prop_object.default = value
                    else:
                        setattr(lcz, prop, value)

    def lcz_rgb_to_index(self, rgb: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the LCZ index from an LCZ RGB array.

        Args:
            rgb: RGB-array with shape (3, ...).

        Raises:
            ValueError: RGB array does not have the correct shape.

        Returns:
            LCZ index array.
        """
        if rgb.shape[0] != 3:
            raise ValueError("rgb array must have shape (3, ...).")
        r, g, b = rgb
        lcz_index = ma.masked_all(r.shape, dtype="uint8")
        for lcz in self.index.values():
            if isinstance(lcz, LCZ):
                lcz_index = ma.where(
                    (r == lcz.r) & (g == lcz.g) & (b == lcz.b), lcz.index, lcz_index
                )
        return lcz_index

    def lcz_index_to_rgb(self, lcz_index: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the LCZ RGB array from an LCZ index array.

        Args:
            lcz_index: LCZ index array.

        Returns:
            LCZ RGB-array with shape (3, ...).
        """
        rgb = ma.masked_all((3,) + lcz_index.shape, dtype="uint8")
        u, inv = np.unique(lcz_index, return_inverse=True)
        for i, band in enumerate(["r", "g", "b"]):
            rgb[i, ...] = np.array([getattr(self.index[x], band) for x in u])[inv].reshape(
                lcz_index.shape
            )
        return rgb

    def value_from_lcz_map(
        self,
        lcz_map: ma.MaskedArray,
        lcz_func: Callable[..., ma.MaskedArray],
        **kwargs: Any,
    ) -> ma.MaskedArray:
        """General function to calculate values from a LCZ array.

        The input function is applied to each pixel with **kwargs as additional arguments.

        Args:
            lcz_map: LCZ index array.
            lcz_func: Function to calculate the value from a LCZ.
            **kwargs: Additional arguments for the function.

        Returns:
            Value array with the same shape as the LCZ index array.
        """
        # Run function to get dimension and dtype.
        test_value = lcz_func(self.compact_highrise, **kwargs)

        # Calculate the value for each unique LCZ only once:
        # Get unique LCZs from lcz_map.
        unique_lcz, inverse = np.unique(lcz_map, return_inverse=True)
        # Results of each unique LCZ.
        unique_results = ma.masked_all(test_value.shape + unique_lcz.shape, dtype=test_value.dtype)
        for i, u in enumerate(unique_lcz):
            unique_results[..., i] = lcz_func(self.index[u], **kwargs)

        # Get results for each pixel from unique results.
        value = unique_results[..., inverse].reshape(test_value.shape + lcz_map.shape)
        return value

    def water_type_from_lcz_map(self, lcz_map: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the water type from an LCZ index array.

        Args:
            lcz_map: LCZ index array.

        Returns:
            Water type array with the same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.water_type_from_lcz)

    def water_type_from_lcz(self, lcz: LCZ) -> ma.MaskedArray:
        """Calculate the water type for a Local Climate Zone.

        Args:
            lcz: Local Climate Zone.

        Returns:
            Water type.
        """
        return lcz.water_type

    def vegetation_type_from_lcz_map(self, lcz_map: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the vegetation type from an LCZ index array.

        Args:
            lcz_map: LCZ index array.

        Returns:
            Vegetation type array with the same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.vegetation_type_from_lcz)

    def vegetation_type_from_lcz(self, lcz: LCZ) -> ma.MaskedArray:
        """Calculate the vegetation type for a Local Climate Zone.

        Args:
            lcz: Local Climate Zone.

        Returns:
            Vegetation type.
        """
        return lcz.vegetation_type

    def urban_fraction_from_lcz_map(self, lcz_map: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the urban fraction from an LCZ index array.

        Args:
            lcz_map: LCZ index array.

        Returns:
            Urban fraction array with the same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.urban_fraction_from_lcz)

    def urban_fraction_from_lcz(self, lcz: LCZ) -> ma.MaskedArray:
        """Calculate the urban fraction for a Local Climate Zone.

        Args:
            lcz: Local Climate Zone.

        Returns:
            Urban fraction.
        """
        if (
            lcz.building_plan_area_fraction.default is None
            or lcz.impervious_plan_area_fraction.default is None
        ):
            raise ValueError("Building or impervious fraction is None.")
        urban_fraction: ma.MaskedArray = ma.MaskedArray(
            lcz.building_plan_area_fraction.default + lcz.impervious_plan_area_fraction.default
        )
        return urban_fraction

    def urban_class_fraction_from_lcz_map(self, lcz_map: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the urban class from an LCZ index array.

        Args:
            lcz_map: LCZ index array.

        Returns:
            Urban class with the same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.urban_class_fraction_from_lcz)

    def urban_class_fraction_from_lcz(self, lcz: LCZ) -> ma.MaskedArray:
        """Calculate the urban class for a Local Climate Zone.

        We assume that there is only one urban class.

        Args:
            lcz: Local Climate Zone.

        Returns:
            Urban fraction.
        """
        urban_class_fraction = ma.masked_all(1)
        if lcz in self.urban_like:
            urban_class_fraction[0] = 1.0
        return urban_class_fraction

    def street_direction_fraction_from_lcz_map(
        self, lcz_map: ma.MaskedArray, udir: List[float]
    ) -> ma.MaskedArray:
        """Calculate the street direction fraction from an LCZ array.

        We assume that there is only one urban class.

        Args:
            lcz_map: LCZ index array.
            udir: List of street directions.

        Returns:
            Street direction fraction array for each urban class and street direction, otherwise the
            same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.street_direction_fraction_from_lcz, udir=udir)

    def street_direction_fraction_from_lcz(self, lcz: LCZ, udir: List[float]) -> ma.MaskedArray:
        """Calculate the street direction fraction for a Local Climate Zone.

        We assume that there is only one urban class.

        Args:
            lcz: Local Climate Zone.
            udir: List of street directions.

        Returns:
            Street direction fraction array for each urban class and street direction.
        """
        street_direction_fraction = ma.masked_all((1, len(udir)))
        if lcz in self.urban_like:
            street_direction_fraction[0, :] = ma.repeat(1.0 / len(udir), len(udir))
        return street_direction_fraction

    def street_width_from_lcz_map(
        self, lcz_map: ma.MaskedArray, udir: List[float]
    ) -> ma.MaskedArray:
        """Calculate the street width from an LCZ array.

        We assume that there is only one urban class.

        Args:
            lcz_map: LCZ index array.
            udir: List of street directions.

        Returns:
            Street width array for each urban class and street direction, otherwise the
            same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.street_width_from_lcz, udir=udir)

    def street_width_from_lcz(self, lcz: LCZ, udir: List[float]) -> ma.MaskedArray:
        """Calculate the street width for a Local Climate Zone.

        We assume that there is only one urban class.

        Args:
            lcz: Local Climate Zone.
            udir: List of street directions.

        Returns:
            Street width array for each urban class and street direction.
        """
        if lcz.height_roughness_elements.default is None or lcz.aspect_ratio.default is None:
            raise ValueError("Height or aspect ratio is None.")
        # shape: one urban class, udir
        street_width = ma.masked_all((1, len(udir)))
        if lcz in self.urban_like:
            street_width[0, :] = ma.repeat(
                lcz.height_roughness_elements.default / lcz.aspect_ratio.default,
                len(udir),
            )
        return street_width

    def building_width_from_lcz_map(
        self, lcz_map: ma.MaskedArray, udir: List[float]
    ) -> ma.MaskedArray:
        """Calculate the building width from an LCZ array.

        We assume that there is only one urban class.

        Args:
            lcz_map: LCZ index array.
            udir: List of street directions.

        Returns:
            Building width array for each urban class and street direction, otherwise the
            same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.building_width_from_lcz, udir=udir)

    def building_width_from_lcz(self, lcz: LCZ, udir: List[float]) -> ma.MaskedArray:
        """Calculate the building width for a Local Climate Zone.

        We assume that there is only one urban class.

        Args:
            lcz: Local Climate Zone.
            udir: List of street directions.

        Returns:
            Building width array for each urban class and street direction.
        """
        if (
            lcz.building_plan_area_fraction.default is None
            or lcz.impervious_plan_area_fraction.default is None
        ):
            raise ValueError("Building or impervious fraction is None.")
        # shape: one urban class, udir
        building_width = ma.masked_all((1, len(udir)))
        if lcz in self.urban_like:
            building_width[0, :] = ma.repeat(
                lcz.building_plan_area_fraction.default / lcz.impervious_plan_area_fraction.default,
                len(udir),
            )
            building_width *= self.street_width_from_lcz(lcz, udir)
        return building_width

    def building_height_from_lcz_map(
        self, lcz_map: ma.MaskedArray, z_uhl: List[float], udir: List[float]
    ) -> ma.MaskedArray:
        """Calculate the building height fraction from an LCZ array.

        We assume that there is only one urban class.

        Args:
            lcz_map: LCZ index array.
            z_uhl: List of urban half-level heights.
            udir: List of street directions.

        Returns:
            Building height fraction array for each urban class, street direction, and urban
            half-level, otherwise the same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(
            lcz_map, self.building_height_from_lcz, z_uhl=z_uhl, udir=udir
        )

    def building_height_from_lcz(
        self, lcz: LCZ, z_uhl: List[float], udir: List[float]
    ) -> ma.MaskedArray:
        """Calculate the building height fractions for a Local Climate Zone.

        A truncated normal distribution of the building height is assumend:
          mean: lcz.height_roughness_elements.default
          standard deviation: 0.25 * (lcz.height_roughness_elements.maximum -
                                      lcz.height_roughness_elements.minimum)
          values are truncated at +- 2 standard deviations from mean
            (if `default` not set, truncated at `maximum` and `minimum`)

        We assume that there is only one urban class.

        Args:
            lcz: Local Climate Zone.
            z_uhl: List of urban half-level heights.
            udir: List of street directions.

        Returns:
            Building height fraction array for each urban class, street direction, and urban
            half-level.
        """
        building_height = ma.masked_all((1, len(udir), len(z_uhl)))
        if lcz in self.urban_like:
            if (
                lcz.height_roughness_elements.maximum is not None
                and lcz.height_roughness_elements.minimum is not None
            ):
                sd = 0.25 * (
                    self._mean_f1(lcz.height_roughness_elements.maximum)
                    - self._mean_f1(lcz.height_roughness_elements.minimum)
                )
            else:
                raise ValueError("No height range specified for LCZ.")

            # Height distribution. It is a truncated normal distribution with
            # mean at default height and sd as standard deviation
            # truncated at +- 2 standard deviations from mean (including the borders)
            # default behaviour: default = 0.5 * (maximum + minimum) so
            #                    truncated at maximum and minimum
            height_distribution = scipy.stats.truncnorm(
                -2.0000000000001,
                +2.0000000000001,
                loc=self._mean_f1(lcz.height_roughness_elements.default),
                scale=sd,
            )
            z_uhl_array = np.array(z_uhl)
            # Half of layer thickness.
            dz_uhl = (z_uhl_array[1:] - z_uhl_array[:-1]) / 2.0
            # Assume same layer thickness above last height
            dz_uhl = np.append(dz_uhl, dz_uhl[-1])
            # Urban full layers.
            z_ufl = z_uhl_array + dz_uhl
            # Add 0 as lower border.
            z_ufl = np.insert(z_ufl, 0, 0.0)
            # height_prob = np.zeros_like(z_uhl_array)
            # Height probability for each layer:
            #   integral over pdf between layer borders = difference of cdf at borders
            upper_cdf = height_distribution.cdf(self._mean_f1(z_ufl[1:]))
            # z_ufl[0] = 0.0
            #    with geometric mean: log(z_ufl[0]) -> -inf so cdf(log(z_ufl[0])) -> 0.0
            if self._mean_f1 == np.log:
                lower_cdf = height_distribution.cdf(self._mean_f1(z_ufl[1:-1]))
                lower_cdf = np.insert(lower_cdf, 0, 0.0)
            else:
                lower_cdf = height_distribution.cdf(self._mean_f1(z_ufl[:-1]))
            height_prob = upper_cdf - lower_cdf
            # Normalize because pdf not completely covered
            height_prob /= height_prob.sum()
            building_height[..., :] = height_prob
        return building_height

    def lai_from_lcz_map(self, lcz_map: ma.MaskedArray) -> ma.MaskedArray:
        """Calculate the Leaf Area Index from an LCZ index array.

        Args:
            lcz_map: LCZ index array.

        Returns:
            LAI array with the same shape as the LCZ index array.
        """
        return self.value_from_lcz_map(lcz_map, self.lai_from_lcz)

    def lai_from_lcz(self, lcz: LCZ) -> ma.MaskedArray:
        """Calculate the Leaf Area Index for a Local Climate Zone.

        Args:
            lcz: Local Climate Zone.

        Returns:
            LAI.
        """
        return lcz.lai
