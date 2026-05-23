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

"""Module for resolved vegetation."""

import logging
from importlib.resources import files
from math import pi
from typing import Callable, Final, List, Optional, Tuple, Union, cast

import numpy as np
import numpy.ma as ma
import numpy.ma.core as ma_core
import pandas.api.typing as pdtypes
import scipy
import scipy.integrate as integrate
from pydantic import BaseModel, Field, field_validator

from palm_csd import StatusLogger
from palm_csd.csd_config import (
    _check_string,
    _default_not_none,
    value_defaults,
)
from palm_csd.tools import is_missing

# TODO Python 3.11 supports Self in typing
# from typing import Self

# Module logger. In __init__.py, it is ensured that the logger is a StatusLogger. For type checking,
# do explicit cast.
logger = cast(StatusLogger, logging.getLogger(__name__))


class Tree(BaseModel, validate_assignment=True):
    """A tree object with size parameters."""

    shape: int = Field(
        ge=value_defaults["tree_shape"].minimum,
        le=value_defaults["tree_shape"].maximum,
    )
    """General shape of the tree:
    1 sphere or ellipsoid
    2 cylinder
    3 cone
    4 inverted cone
    5 paraboloid (rounded cone)
    6 inverted paraboloid (inverted rounded cone)
    """
    crown_ratio: float = Field(
        ge=value_defaults["tree_crown_ratio"].minimum,
        le=value_defaults["tree_crown_ratio"].maximum,
    )
    """Ratio of maximum crown height to maximum crown diameter."""
    crown_diameter: float = Field(
        ge=value_defaults["tree_crown_diameter"].minimum,
        le=value_defaults["tree_crown_diameter"].maximum,
    )
    """Crown diameter (m)."""
    height: float = Field(
        ge=value_defaults["tree_height"].minimum,
        le=value_defaults["tree_height"].maximum,
    )
    """Total height of the tree including trunk (m)."""
    bad_scale: float = Field(
        ge=value_defaults["bad_scale"].minimum,
        le=value_defaults["bad_scale"].maximum,
    )
    """Ratio of basal area in the crown area to the leaf area."""
    trunk_diameter: float = Field(
        ge=value_defaults["tree_trunk_diameter"].minimum,
        le=value_defaults["tree_trunk_diameter"].maximum,
    )
    """Trunk diameter at breast height (1.4 m) (m)."""

    z_max_rel: float = Field(
        ge=value_defaults["lad_z_max_rel"].minimum,
        le=value_defaults["lad_z_max_rel"].maximum,
    )
    """Height where the leaf area density is maximum relative to total tree height, only used with
    Lalic and Mihailovic (2004) method."""
    alpha: float = Field(
        ge=value_defaults["lad_alpha"].minimum,
        le=value_defaults["lad_alpha"].maximum,
    )
    """LAD profile parameter alpha only used with Markkanen et al. (2003) method."""
    beta: float = Field(
        ge=value_defaults["lad_beta"].minimum,
        le=value_defaults["lad_beta"].maximum,
    )
    """LAD profile parameter alpha only used with Markkanen et al. (2003) method."""


class ReferenceTree(Tree):
    """A tree object with additional reference parameters."""

    species: str
    """Species."""
    lai_summer: float = Field(
        ge=value_defaults["lai"].minimum,
        le=value_defaults["lai"].maximum,
    )
    """Default leaf area index fully leafed."""
    lai_winter: float = Field(
        ge=value_defaults["lai"].minimum,
        le=value_defaults["lai"].maximum,
    )
    """Default winter-time leaf area index."""


class DomainTree(Tree):
    """A tree object with all the necessary parameters to calculate the LAD and BAD fields."""

    lai: float = Field(
        ge=value_defaults["lai"].minimum,
        le=value_defaults["lai"].maximum,
    )
    """Actual leaf area index."""
    i: int = Field(
        ge=0,
        le=value_defaults["nx"].maximum,
    )
    """x coordinate of the tree."""
    j: int = Field(
        ge=0,
        le=value_defaults["ny"].maximum,
    )
    """y coordinate of the tree."""
    id: int = Field(
        ge=value_defaults["tree_id"].minimum,
        le=value_defaults["tree_id"].maximum,
    )
    """ID of the tree."""
    type: int = Field(
        ge=value_defaults["tree_type"].minimum,
        le=value_defaults["tree_type"].maximum,
    )
    """Type of the tree."""


def _populate_defaults() -> List[ReferenceTree]:
    """Read default tree species data from file."""
    with files("palm_csd.data").joinpath("tree_defaults.csv").open() as tree_csv:
        tree_data = np.genfromtxt(
            tree_csv,
            delimiter=",",
            dtype=None,
            names=True,
            skip_header=12,
            encoding="utf-8",
        )
    defaults = []
    for tree in tree_data:
        defaults.append(
            ReferenceTree(
                species=tree["species"],
                shape=tree["shape"],
                crown_ratio=tree["crown_ratio"],
                crown_diameter=tree["crown_diameter"],
                height=tree["height"],
                lai_summer=tree["lai_summer"],
                lai_winter=tree["lai_winter"],
                z_max_rel=tree["z_max_rel"],
                alpha=tree["alpha"],
                beta=tree["beta"],
                bad_scale=tree["bad_scale"],
                trunk_diameter=tree["trunk_diameter"],
            )
        )
    return defaults


tree_defaults = _populate_defaults()


class CanopyGenerator(BaseModel):
    """Generate a canopy for vegatation patches.

    Two methods are possible:
    - LM2004: Lalic and Mihailovic (2004)
    - Metal2003: Markkanen et al. (2003)
    """

    N_BELOW_Z_MAX_LM2004: Final = 6.0
    """Exponent below z_max in Lalic and Mihailovic (2004)."""
    N_ABOVE_Z_MAX_LM2004: Final = 0.5
    """Exponent above z_max in Lalic and Mihailovic (2004)."""

    dz: float
    """Grid spacing in z direction."""
    pixel_size: float
    """Pixel size in m."""

    method: str = "Metal2003"
    """Method to calculate the LAD profile."""
    alpha_Metal2003: float = Field(
        default=_default_not_none("lad_alpha"),
        ge=value_defaults["lad_alpha"].minimum,
        le=value_defaults["lad_alpha"].maximum,
    )
    """LAD profile parameter alpha for Markkanen et al. (2003)."""
    beta_Metal2003: float = Field(
        default=_default_not_none("lad_beta"),
        ge=value_defaults["lad_beta"].minimum,
        le=value_defaults["lad_beta"].maximum,
    )
    """LAD profile parameter beta for Markkanen et al. (2003)."""
    z_max_rel_LM2004: float = Field(
        default=_default_not_none("lad_z_max_rel"),
        ge=value_defaults["lad_z_max_rel"].minimum,
        le=value_defaults["lad_z_max_rel"].maximum,
    )
    """Height of maximum LAD divided by patch height for Lalic and Mihailovic (2004)."""

    # Private underscore methods not treated as fields by Pydantic.
    _lad_norm_fun: Callable[..., ma.MaskedArray]
    """Function to calculate normalized 3D LAD (LAD_max / LAI * h)."""
    _lad_max_norm_fun: Callable[..., float]
    """Function to calculate maximum normalized LAD (LAD_max / LAI * h)."""
    _z_max_rel_fun: Callable[..., float]
    """Function to calculate height of maximum LAD divided by patch height."""

    season: str = "summer"
    """Season for which the LAI is used."""

    height_rel_resolved_vegetation_lower_threshold: float = Field(
        default=_default_not_none("height_rel_resolved_vegetation_lower_threshold"),
        ge=value_defaults["height_rel_resolved_vegetation_lower_threshold"].minimum,
        le=value_defaults["height_rel_resolved_vegetation_lower_threshold"].maximum,
    )
    """Lower threshold for relative tree height below which trees are ignored."""
    lai_tree_lower_threshold: float = Field(
        default=_default_not_none("lai_tree_lower_threshold"),
        ge=value_defaults["lai_tree_lower_threshold"].minimum,
        le=value_defaults["lai_tree_lower_threshold"].maximum,
    )
    """Lower threshold for tree LAI below which trees are ignored."""
    remove_low_lai_tree: bool = False
    """Remove trees with low LAI."""

    shallow_tree_count: int = 0
    """Counter of trees with height < 1/2 dz."""
    low_lai_count: int = 0
    """Counter of trees with low LAI."""
    mod_count: int = 0
    """Counter of modified trees."""
    id_count: int = 0
    """Tree ID counter."""

    tree_sphere_extinction: float = 0.6
    """Sphere extinction coefficient of the LAD/BAD generator; experimental."""
    tree_cone_extinction: float = 0.2
    """Cone extinction coefficient of the LAD/BAD generator; experimental."""

    @field_validator("method")
    @classmethod
    def _check_method(cls, value: str) -> str:
        """Check if the method is valid.

        Args:
            value: Value to check.

        Returns:
            Validated value.
        """
        _check_string(value, ["Metal2003", "LM2004"])
        return value

    @field_validator("season")
    @classmethod
    def _check_season(cls, value: str) -> str:
        """Check if the season is valid.

        Args:
            value: Value to check.

        Returns:
            Validated value.
        """
        _check_string(value, ["summer", "winter"])
        return value

    def __init__(self, **kwargs) -> None:
        """Use standard __init__ and, depending on `self.method`, set the appropriate functions.

        Raises:
            ValueError: Unknown patch LAD method.
        """
        super().__init__(**kwargs)
        if self.method == "LM2004":
            self._lad_norm_fun = self._lad_norm_LM2004
            self._lad_max_norm_fun = self._lad_max_norm_LM2004
            self._z_max_rel_fun = self._z_max_rel_LM2004
        elif self.method == "Metal2003":
            self._lad_norm_fun = self._lad_norm_Metal2003
            self._lad_max_norm_fun = self._lad_max_norm_Metal2003
            self._z_max_rel_fun = self._z_max_rel_Metal2003
        else:
            raise ValueError(f"Unknown patch LAD method {self.method}.")

    def _integral_LM2004(
        self, z_rel: Union[ma.MaskedArray, float], n: float, **kwargs: float
    ) -> Union[ma.MaskedArray, float]:
        """Indefinite integral per patch height h of the LAD profile in Lalic and Mihailovic (2004).

        According to Wolfram Alpha, the indefinite integral I of (1) Lalic and Mihailovic (2004) is:
        I = -e^n n^(1 - n) (h - zm) Γ(n - 1, ((h - zm) n)/(h - z)) + constant
        where Γ(s, x) is the non-regularized upper incomplete gamma function.
        Here, we use zm/h and z/h instead of zm and z, respectively, thus, we get I/h.

        The integral is evaluated at input height with values up to a constant. This constant
        cancels out when calculating the difference of the integrals at different heights.

        Args:
            z_rel: Height z divided by patch height.
            n: Exponent in the LAD profile.
            **kwargs: Optionally set z_max_rel, otherwise the instance attribute is used.

        Returns:
            Indefinite integral at z_rel of the LAD profile.
        """
        # Use z_max_rel from input or instance default.
        z_max_rel = kwargs.get("z_max_rel", self.z_max_rel_LM2004)

        # 1 - zm/h
        dphzm = 1.0 - z_max_rel
        # 1 - z/h
        dphz = 1.0 - z_rel

        # (h - zm)/(h - z) = (1 - zm/h)/(1 - z/h)
        # Avoid division by zero with the where in ma.divide.
        x = ma.where(dphz == 0.0, np.inf, ma.divide(dphzm * n, dphz, where=(dphz != 0.0)))

        # Non-regularized upper incomplete gamma, in scipy not defined for negative first argument.
        if n < 1.0:
            # For negative first argument, use recurrence relation
            # Γ(s+1, x) = sΓ(s, x) + x^s e^-x  .
            # Checked with mpmath.gammainc.
            # Multiply regularized upper incomplete gammaincc with gamma to get non-regularized.
            gamma_inc = (
                scipy.special.gammaincc(n, x) * scipy.special.gamma(n) - x ** (n - 1.0) * np.exp(-x)
            ) / (n - 1.0)
        else:
            # Multiply regularized upper incomplete gammaincc with gamma to get non-regularized.
            gamma_inc = scipy.special.gammaincc(n - 1.0, x) * scipy.special.gamma(n - 1.0)

        return -np.exp(n) * n ** (1.0 - n) * dphzm * gamma_inc

    def _lad_norm_LM2004(self, z_rel: ma.MaskedArray, **kwargs: float) -> ma.MaskedArray:
        """Calculate normalized 3D LAD for vegetation patches following Lalic and Mihailovic (2004).

        The Lalic and Mihailovic (2004) LAD profile is assumed from height 0 to z_rel. The LAD of
        the highest, likely not completely filled grid cell is distributed over the whole grid cell.
        For example, for a grid spacing of 3m and a patch height of 10m, the LAD of the height 9m
        (z_rel=0.9) to 10m (z_rel=0.9) is distributed over the whole 9m (z_rel=0.9) to 12m
        (z_rel=1.2) grid cell.

        The calculated values are normalized by LAI and patch height (LAD / LAI * h). The average
        LAD of a grid cell is calculated by integrating the LAD profile over the height of this cell
        and dividing by the cell height. Thus, the LAD is proportional to the integral of (1) of
        Lalic and Mihailovic (2004) at the upper boundary minus the integral at the lower boundary
        of the cell. The integral is calculated in _integral_LM2004. Note that the integral depends
        on the exponent n in the LAD profile, which is different above and below the maximum LAD
        height. In particular, the grid cell that includes zm needs to be split into two parts.

        Args:
            z_rel: Height z divided by patch height.
            **kwargs: Optionally set z_max_rel, otherwise the instance attribute is used.

        Returns:
            Normalized 3D LAD.
        """
        # Use z_max_rel from input or instance default.
        z_max_rel = kwargs.get("z_max_rel", self.z_max_rel_LM2004)

        # Size of vertical grid cells. Here, the full grid cell height is calculated. This is
        # relevant for the highest grid point with z_rel>=1.0, i.e. the upper grid cell level is
        # potentially above the patch height.
        dz_rel = z_rel[1:, ...] - z_rel[:-1, ...]
        # Limit z_rel to 1.0.
        # 1.0 should only be reached for the top grid points. This ensures that the profile is
        # calculated correctly. The normalization by dz_rel below ensures that the LAD is
        # distributed correctly over the whole grid cell.
        z_rel = ma.where(z_rel > 1.0, 1.0, z_rel)

        # Prepare pdf integral array.
        pdf_int = ma.masked_all_like(z_rel)

        # Different factor n in LAD PDF above and below zm.
        below_z_max_rel = z_rel < z_max_rel
        above_z_max_rel = ~below_z_max_rel & ~z_rel.mask
        below_z_max_rel = below_z_max_rel & ~z_rel.mask  # type: ignore
        pdf_int[below_z_max_rel] = self._integral_LM2004(
            z_rel[below_z_max_rel], self.N_BELOW_Z_MAX_LM2004
        )
        pdf_int[above_z_max_rel] = self._integral_LM2004(
            z_rel[above_z_max_rel], self.N_ABOVE_Z_MAX_LM2004
        )

        # Grid cells that include z_max_rel need to be split into two parts.
        # Find grid cells that include z_max_rel.
        i_z_max_rel = ma.apply_along_axis(np.searchsorted, 0, z_rel, z_max_rel)
        # Calculate integral for zm coming from below and from above.
        int_zm_low = self._integral_LM2004(z_max_rel, self.N_BELOW_Z_MAX_LM2004)
        int_zm_up = self._integral_LM2004(z_max_rel, self.N_ABOVE_Z_MAX_LM2004)

        # Prepare output array.
        lad_norm = ma.masked_all_like(z_rel)

        # LAD proportional to integral upper - integral lower boundary.
        lad_norm[1:, ...] = pdf_int[1:, ...] - pdf_int[:-1, ...]

        # Cells that include zm need to be corrected:
        # LAD proportional to
        #      integral upper - integral zm upper + integral zm lower - integral lower .
        lad_norm[i_z_max_rel, ...] = lad_norm[i_z_max_rel, ...] - int_zm_up + int_zm_low

        # Normalize LAD.
        lad_norm[1:, ...] = lad_norm[1:, ...] / dz_rel[:, ...] / np.sum(lad_norm[1:, ...], axis=0)

        # Set LAD to zero at lowest level for columns that include LAD values.
        lad_norm[0, ...] = ma.where(np.any(~lad_norm.mask, axis=0), 0.0, ma.masked)

        return lad_norm

    def _lad_norm_Metal2003(self, z_rel: ma.MaskedArray, **kwargs: float) -> ma.MaskedArray:
        """Calculate normalized 3D LAD for vegetation patches following Markkanen et al. (2003).

        The Markkanen et al. (2003) LAD profile is assumed from height 0 to z_rel. The LAD of the
        highest, likely not completely filled grid cell is distributed over the whole grid cell. For
        example, for a grid spacing of 3m and a patch height of 10m, the LAD of the height 9m
        (z_rel=0.9) to 10m (z_rel=0.9) is distributed over the whole 9m (z_rel=0.9) to 12m
        (z_rel=1.2) grid cell.

        The calculated values are normalized by LAI and patch height (LAD / LAI * h). The average
        LAD of a grid cell is calculated by integrating the LAD profile over the height of this cell
        and dividing by the cell height. Thus, the LAD is proportional to the integral of (1) of
        Markkanen et al. (2003) at the upper boundary minus the integral at the lower boundary of
        the cell.

        Args:
            z_rel: Height z divided by patch height.
            **kwargs: Optionally set alpha and beta, otherwise the instance attribute is used.

        Returns:
            Normalized 3D LAD.
        """
        # Use alpha and beta from input or instance default.
        alpha = kwargs.get("alpha", self.alpha_Metal2003)
        beta = kwargs.get("beta", self.beta_Metal2003)

        # Size of vertical grid cells. Here, the full grid cell height is calculated. This is
        # relevant for the highest grid point with z_rel>=1.0, i.e. the upper grid cell level is
        # potentially above the patch height.
        dz_rel = z_rel[1:, ...] - z_rel[:-1, ...]
        # Limit z_rel to 1.0.
        # 1.0 should only be reached for the top grid points. This ensures that the profile is
        # calculated correctly. The normalization by dz_rel below ensures that the LAD is
        # distributed correctly over the whole grid cell.
        z_rel = ma.where(z_rel > 1.0, 1.0, z_rel)

        # Prepare output array.
        lad_norm = ma.masked_all_like(z_rel)

        # Regularized (normalized) incomplete beta function corresponds to integral 0 to z_rel dz/h
        # of Markannen et al. (2003), Eq. (1) and, thus, the fraction of leafs between 0 and z_rel.
        # We would need the integral dz instead dz/h, thus we miss a factor h here, which cancels
        # out below.
        pdf_int = scipy.special.betainc(alpha, beta, z_rel)

        # Fraction of leaves between height levels, already normalized.
        lad_norm[1:, ...] = (pdf_int[1:, ...] - pdf_int[:-1, ...]) / dz_rel[:, ...]

        # Set LAD to zero at lowest level for columns that include LAD values.
        lad_norm[0, ...] = ma.where(np.any(~lad_norm.mask, axis=0), 0.0, ma.masked)

        return lad_norm

    def _z_max_rel_Metal2003(self, **kwargs: float) -> float:
        """Calculate relative height of maximum LAD following Markkanen et al. (2003).

        The height of the maximum LAD is divided by patch height.

        Args:
            **kwargs: Optionally set alpha and beta, otherwise the instance attribute is used.

        Returns:
            Relative height of maximum LAD.
        """
        # Use alpha and beta from input or instance default.
        alpha = kwargs.get("alpha", self.alpha_Metal2003)
        beta = kwargs.get("beta", self.beta_Metal2003)

        return (alpha - 1.0) / (alpha + beta - 2.0)

    def _z_max_rel_LM2004(self, **kwargs: float) -> float:
        """Return relative height of maximum LAD following Lalic and Mihailovic (2004).

        The height of the maximum LAD is divided by patch height. This number is the input parameter
        of this model.

        Returns:
            Relative height of maximum LAD.
        """
        return kwargs.get("z_max_rel", self.z_max_rel_LM2004)

    def _lad_max_norm_Metal2003(self, **kwargs: float) -> float:
        """Calculate maximum normalized LAD following Markkanen et al. (2003).

        The calculated values are normalized by LAI and patch height (LAD_max / LAI * h).

        Returns:
            Maximum normalized LAD.
        """
        # Use alpha and beta from input or instance default.
        alpha = kwargs.get("alpha", self.alpha_Metal2003)
        beta = kwargs.get("beta", self.beta_Metal2003)

        # Insert z_max_rel into (1) of Markkanen et al. (2003).
        z_max_rel = self._z_max_rel_Metal2003(alpha=alpha, beta=beta)
        return (
            z_max_rel ** (alpha - 1.0)
            * (1.0 - z_max_rel) ** (beta - 1.0)
            / scipy.special.beta(alpha, beta)
        )

    def _lad_max_norm_LM2004(self, **kwargs: float) -> float:
        """Calculate maximum normalized LAD following Lalic and Mihailovic (2004).

        The calculated values are normalized by LAI and patch height (LAD_max / LAI * h).

        Returns:
            Maximum normalized LAD.
        """
        # Use z_max_rel from input or instance default.
        z_max_rel = kwargs.get("z_max_rel", self.z_max_rel_LM2004)

        # LAD maximum derived from the total integral in (2) of Lalic and Mihailovic (2004).
        pdf_int_0 = self._integral_LM2004(0.0, self.N_BELOW_Z_MAX_LM2004)
        pdf_int_zmdph1 = self._integral_LM2004(z_max_rel, self.N_BELOW_Z_MAX_LM2004)
        pdf_int_zmdph2 = self._integral_LM2004(z_max_rel, self.N_ABOVE_Z_MAX_LM2004)
        pdf_int_1 = self._integral_LM2004(1.0, self.N_ABOVE_Z_MAX_LM2004)

        pdf_int = pdf_int_zmdph1 - pdf_int_0 + pdf_int_1 - pdf_int_zmdph2

        return 1.0 / pdf_int

    def process_patch(
        self,
        patch_height: ma.MaskedArray,
        patch_type_2d: ma.MaskedArray,
        patch_lai: ma.MaskedArray,
    ) -> Tuple[ma.MaskedArray, ma.MaskedArray, ma.MaskedArray]:
        """Calculate 3D leaf area density for vegetation patches according to the chosen method.

        The LAD profiles are assumed from height 0m to patch height. The LAD of the highest, likely
        not completely filled grid cell is distributed over the whole grid cell. For example, for a
        grid spacing of 3m and a patch height of 10m, the LAD of the height 9m to 10m is distributed
        over the whole 9m to 12m grid cell.

        Args:
            patch_height: Vegetation patch height.
            patch_type_2d: Vegetation patch type.
            patch_lai: Vegetation patch leaf area index.

        Returns:
            3d leaf area density, patch id, and patch type.
        """
        # Maximum canopy height in patch and outside, and corresponding vertical grid indices
        canopy_height_max = ma.max(patch_height)
        iz_1d = np.arange(0, np.ceil(canopy_height_max / self.dz) + 1, dtype="int")

        # Array of number of vertical grid points in vegetation patch
        iz_max = ma.ceil(patch_height / self.dz).astype("int")

        # Broadcast iz_1d to 3D array and mask values larger than iz_max.
        iz = np.broadcast_to(iz_1d[:, np.newaxis, np.newaxis], iz_1d.shape + patch_height.shape)
        iz = ma.masked_where(iz[:, :, :] > iz_max[np.newaxis, :, :], iz, copy=False)

        # Grid point height divided by patch height.
        z_rel = iz * self.dz / patch_height[np.newaxis, :, :]

        # Calculate normalized LAD with chosen method.
        lad_3d = self._lad_norm_fun(z_rel)

        # Average LAD
        lad_av = patch_lai / patch_height

        # Calculate LAD: Fraction of LAI between height levels, divided by dz. The (not included)
        # patch height in the integral cancels out.
        lad_3d[:, :, :] = lad_av[np.newaxis, :, :] * lad_3d[:, :, :]

        # Set patch_id and patch_type where LAD is defined.
        patch_id_3d = ma.where(~lad_3d.mask, 1, ma.masked)
        patch_type_3d = ma.where(~lad_3d.mask, patch_type_2d[np.newaxis, :, :], ma.masked)

        return lad_3d, patch_id_3d, patch_type_3d

    def z_max_rel(self) -> float:
        """Height of maximum LAD divided by patch height.

        Returns:
            Height of maximum LAD divided by patch height.
        """
        return self._z_max_rel_fun()

    def lad_max_norm(self) -> float:
        """Maximum normalized LAD (LAD_max / LAI * h).

        Returns:
            Maximum normalized LAD (LAD_max / LAI * h).
        """
        return self._lad_max_norm_fun()

    def generate_tree(
        self,
        i,
        j,
        type: Union[float, ma_core.MaskedConstant, pdtypes.NAType],
        shape: Union[float, ma_core.MaskedConstant, pdtypes.NAType],
        height: Union[float, ma_core.MaskedConstant, pdtypes.NAType],
        lai: Union[float, ma_core.MaskedConstant, pdtypes.NAType],
        crown_diameter: Union[float, ma_core.MaskedConstant, pdtypes.NAType],
        trunk_diameter: Union[float, ma_core.MaskedConstant, pdtypes.NAType],
    ) -> Optional[DomainTree]:
        """Generate a tree object.

        Input values are checked and default values are used if necessary. Depending on the set-up,
        trees with too low tree_height or tree_lai are removed.
        """
        # Check for missing data in the input and set default values if needed.
        if is_missing(type):
            if value_defaults["tree_type"].default is None:
                raise ValueError("Tree type default must be defined.")
            type_checked = int(value_defaults["tree_type"].default)
        else:
            type_checked = int(type)  # type: ignore

        if is_missing(shape):
            shape_checked = tree_defaults[type_checked].shape
        else:
            shape_checked = int(shape)  # type: ignore

        if is_missing(height):
            height_checked = tree_defaults[type_checked].height
        else:
            height_checked = float(height)  # type: ignore

        if is_missing(lai):
            lai_checked = getattr(tree_defaults[type_checked], f"lai_{self.season}")
        else:
            lai_checked = float(lai)  # type: ignore

        if is_missing(crown_diameter):
            crown_diameter_checked = tree_defaults[type_checked].crown_diameter
        else:
            crown_diameter_checked = float(crown_diameter)  # type: ignore

        if is_missing(trunk_diameter):
            trunk_diameter_checked = tree_defaults[type_checked].trunk_diameter
        else:
            trunk_diameter_checked = float(trunk_diameter)  # type: ignore

        # Very small trees are ignored.
        if height_checked <= (self.height_rel_resolved_vegetation_lower_threshold * self.dz):
            self.shallow_tree_count += 1
            logger.debug_indent(
                f"Removed low tree with height = {height_checked:0.1f} at ({i}, {j})."
            )
            return None

        # Check tree_lai.
        # Tree LAI lower than threshold?
        if lai_checked < self.lai_tree_lower_threshold:
            # Deal with low lai tree
            self.mod_count += 1
            if self.remove_low_lai_tree:
                # Skip this tree
                logger.debug_indent(f"Removed tree with LAI = {lai_checked:0.3f} at ({i}, {j}).")
                return None
            else:
                # Use type specific default
                lai_checked = getattr(tree_defaults[type_checked], f"lai_{self.season}")
                logger.debug_indent(f"Adjusted tree to LAI = {lai_checked:0.3f} at ({i}, {j}).")

        # Warn about a tree with lower LAI than we would expect in winter.
        if lai_checked < tree_defaults[type_checked].lai_winter:
            self.low_lai_count += 1
            logger.debug_indent(
                f"Found tree with LAI = {lai_checked:0.3f} (tree-type specific default winter LAI "
                + f"of {tree_defaults[type_checked].lai_winter:0.2}) at ({i}, {j})."
            )

        # Assign values that are not defined as user input from lookup table.
        crown_ratio_checked = tree_defaults[type_checked].crown_ratio
        z_max_rel_checked = tree_defaults[type_checked].z_max_rel
        alpha_checked = tree_defaults[type_checked].alpha
        beta_checked = tree_defaults[type_checked].beta
        bad_scale_checked = tree_defaults[type_checked].bad_scale

        # Increase the tree ID counter.
        self.id_count += 1

        return DomainTree(
            type=type_checked,
            shape=shape_checked,
            crown_ratio=crown_ratio_checked,
            crown_diameter=crown_diameter_checked,
            height=height_checked,
            z_max_rel=z_max_rel_checked,
            alpha=alpha_checked,
            beta=beta_checked,
            bad_scale=bad_scale_checked,
            trunk_diameter=trunk_diameter_checked,
            lai=lai_checked,
            i=i,
            j=j,
            id=self.id_count,
        )

    def check_tree_counters(self) -> None:
        """Print a summary of the tree and tree LAI adjustments.

        Args:
            config: Domain configuration
            settings: Settings configuration
        """
        if self.shallow_tree_count > 0:
            logger.warning(f"Removed {self.shallow_tree_count} low trees with height < 1/2 dz.")
        if self.mod_count > 0:
            if self.remove_low_lai_tree:
                logger.warning(f"Removed {self.mod_count} trees due to low LAI.")
            else:
                logger.warning(
                    f"Adjusted LAI of {self.mod_count} trees below lai_tree_lower_threshold "
                    + f"using tree-type specific default {self.season} LAI."
                )
        if self.low_lai_count > 0:
            logger.warning(
                f"Found {self.low_lai_count} trees with LAI lower then the "
                + "tree-type specific default winter LAI."
            )
            logger.info_indent(
                "Consider adjusting lai_tree_lower_threshold and remove_low_lai_tree."
            )

    def add_tree_to_3d_fields(
        self,
        tree: DomainTree,
        lad_global: ma.MaskedArray,
        bad_global: ma.MaskedArray,
        id_global: ma.MaskedArray,
        type_global: ma.MaskedArray,
    ) -> None:
        """Generate and store the LAD and BAD profile for a single tree.

        The values are stored in the input arrays. Also store tree id and type with the same extent
        as the LAD and BAD arrays in the respective input.

        Args:
            tree: Single tree object.
            lad_global: Global 3D leaf area density.
            bad_global: Global 3D basal area density.
            id_global: Global 3D tree id.
            type_global: Global 3D tree type.

        Raises:
            ValueError: Unknown tree shape.
        """
        # Calculate crown height and height of the crown center.
        crown_height = tree.crown_ratio * tree.crown_diameter
        if crown_height > tree.height:
            crown_height = tree.height

        crown_center = tree.height - crown_height * 0.5

        # Calculate height of maximum LAD.
        z_lad_max = tree.z_max_rel * tree.height

        # Calculate the maximum LAD after Lalic and Mihailovic (2004).
        lad_max_part_1 = integrate.quad(
            lambda z: ((tree.height - z_lad_max) / (tree.height - z)) ** self.N_BELOW_Z_MAX_LM2004
            * np.exp(
                self.N_BELOW_Z_MAX_LM2004 * (1.0 - (tree.height - z_lad_max) / (tree.height - z))
            ),
            0.0,
            z_lad_max,
        )
        lad_max_part_2 = integrate.quad(
            lambda z: ((tree.height - z_lad_max) / (tree.height - z)) ** self.N_ABOVE_Z_MAX_LM2004
            * np.exp(
                self.N_ABOVE_Z_MAX_LM2004 * (1.0 - (tree.height - z_lad_max) / (tree.height - z))
            ),
            z_lad_max,
            tree.height,
        )

        lad_max = tree.lai / (lad_max_part_1[0] + lad_max_part_2[0])

        # Define position of tree and its output domain.
        nx = int(tree.crown_diameter / self.pixel_size) + 2
        nz = int(tree.height / self.dz) + 2

        # Add one grid point if diameter is an odd value.
        if (tree.crown_diameter % 2.0) != 0.0:
            nx = nx + 1

        # Create local domain of the tree's LAD.
        x = np.arange(0, nx * self.pixel_size, self.pixel_size)
        x[:] = x[:] - 0.5 * self.pixel_size
        y = x

        z = np.arange(0, nz * self.dz, self.dz)
        z[1:] = z[1:] - 0.5 * self.dz

        # Define center of the tree position inside the local LAD domain.
        location_x = x[int(nx / 2)]
        location_y = y[int(nx / 2)]

        # Calculate LAD profile after Lalic and Mihailovic (2004).
        # Will be later used for normalization.
        lad_profile = np.arange(0, nz, 1.0)
        lad_profile[:] = 0.0

        for k in range(1, nz - 1):
            if (z[k] > 0.0) & (z[k] < z_lad_max):
                n = self.N_BELOW_Z_MAX_LM2004
            else:
                n = self.N_ABOVE_Z_MAX_LM2004

            lad_profile[k] = (
                lad_max
                * ((tree.height - z_lad_max) / (tree.height - z[k])) ** n
                * np.exp(n * (1.0 - (tree.height - z_lad_max) / (tree.height - z[k])))
            )

        # Create lad array and populate according to the specific tree shape.
        # NOTE This is still experimental
        lad_local = ma.masked_all((nz, nx, nx))
        bad_local = ma.copy(lad_local)

        # Branch for spheres and ellipsoids.
        # A symmetric LAD sphere is created assuming an LAD extinction towards the center of the
        # tree, representing the effect of sunlight extinction and thus less leaf mass inside the
        # tree crown.
        # NOTE Extinction coefficients are experimental.
        if tree.shape == 1:
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(0, nz):
                        r_test = np.sqrt(
                            (x[i] - location_x) ** 2 / (tree.crown_diameter * 0.5) ** 2
                            + (y[j] - location_y) ** 2 / (tree.crown_diameter * 0.5) ** 2
                            + (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** (2)
                        )
                        if r_test <= 1.0:
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.tree_sphere_extinction * (1.0 - r_test)
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for cylinder shapes
        elif tree.shape == 2:
            k_min = int((crown_center - crown_height * 0.5) / self.dz)
            k_max = int((crown_center + crown_height * 0.5) / self.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        r_test = np.sqrt(
                            (x[i] - location_x) ** 2 / (tree.crown_diameter * 0.5) ** 2
                            + (y[j] - location_y) ** 2 / (tree.crown_diameter * 0.5) ** 2
                        )
                        if r_test <= 1.0:
                            r_test3 = np.sqrt(
                                (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** 2
                            )
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.tree_sphere_extinction * (1.0 - max(r_test, r_test3))
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for cone shapes
        elif tree.shape == 3:
            k_min = int((crown_center - crown_height * 0.5) / self.dz)
            k_max = int((crown_center + crown_height * 0.5) / self.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k - k_min
                        r_test = (
                            (x[i] - location_x) ** 2
                            + (y[j] - location_y) ** 2
                            - ((tree.crown_diameter * 0.5) ** 2 / crown_height**2)
                            * (z[k_rel] - crown_height) ** 2
                        )
                        if r_test <= 0.0:
                            r_test2 = np.sqrt(
                                (x[i] - location_x) ** 2 / (tree.crown_diameter * 0.5) ** 2
                                + (y[j] - location_y) ** 2 / (tree.crown_diameter * 0.5) ** 2
                            )
                            r_test3 = np.sqrt(
                                (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** 2
                            )
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.tree_cone_extinction
                                * (1.0 - max((r_test + 1.0), r_test2, r_test3))
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for inverted cone shapes.
        # TODO: what is r_test2 and r_test3 used for? Debugging needed!
        elif tree.shape == 4:
            k_min = int((crown_center - crown_height * 0.5) / self.dz)
            k_max = int((crown_center + crown_height * 0.5) / self.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k_max - k
                        r_test = (
                            (x[i] - location_x) ** 2
                            + (y[j] - location_y) ** 2
                            - ((tree.crown_diameter * 0.5) ** 2 / crown_height**2)
                            * (z[k_rel] - crown_height) ** 2
                        )
                        if r_test <= 0.0:
                            r_test2 = np.sqrt(
                                (x[i] - location_x) ** 2 / (tree.crown_diameter * 0.5) ** 2
                                + (y[j] - location_y) ** 2 / (tree.crown_diameter * 0.5) ** 2
                            )
                            r_test3 = np.sqrt(
                                (z[k] - crown_center) ** 2 / (crown_height * 0.5) ** 2
                            )
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.tree_cone_extinction * (-r_test)
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for paraboloid shapes
        elif tree.shape == 5:
            k_min = int((crown_center - crown_height * 0.5) / self.dz)
            k_max = int((crown_center + crown_height * 0.5) / self.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k - k_min
                        r_test = (
                            (x[i] - location_x) ** 2 + (y[j] - location_y) ** (2)
                        ) * crown_height / (tree.crown_diameter * 0.5) ** 2 - z[k_rel]
                        if r_test <= 0.0:
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.tree_cone_extinction * (-r_test)
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        # Branch for inverted paraboloid shapes
        elif tree.shape == 6:
            k_min = int((crown_center - crown_height * 0.5) / self.dz)
            k_max = int((crown_center + crown_height * 0.5) / self.dz)
            for i in range(0, nx):
                for j in range(0, nx):
                    for k in range(k_min, k_max):
                        k_rel = k_max - k
                        r_test = (
                            (x[i] - location_x) ** 2 + (y[j] - location_y) ** (2)
                        ) * crown_height / (tree.crown_diameter * 0.5) ** 2 - z[k_rel]
                        if r_test <= 0.0:
                            lad_local[k, j, i] = lad_max * np.exp(
                                -self.tree_cone_extinction * (-r_test)
                            )
                        else:
                            lad_local[k, j, i] = ma.masked

        else:
            raise ValueError("Unknown tree shape.")

        # Leave if no LAD was generated.
        if ma.all(lad_local.mask):
            return

        # Indicate a defined LAD in a column by setting lowest value to 0.
        for i in range(0, nx):
            for j in range(0, nx):
                if ma.any(~ma.getmaskarray(lad_local)[:, j, i]):
                    lad_local[0, j, i] = 0.0

        # Normalize according to LAI input. Calculate the average LAI, which should be tree.lai.
        # Correct accordingly.
        lai_avg = (lad_local.sum(axis=0) * self.dz).mean()
        if lai_avg == 0.0:
            return
        lad_local = lad_local * tree.lai / lai_avg

        # Create BAD array and populate.
        # TODO: revise as low LAD inside the foliage does not result in low BAD values.
        bad_local = (1.0 - (lad_local / (ma.max(lad_local) + 0.01))) * 0.1

        # Overwrite grid cells that are occupied by the tree trunk
        radius = tree.trunk_diameter * 0.5
        for i in range(0, nx):
            for j in range(0, nx):
                for k in range(0, nz):
                    if z[k] <= crown_center:
                        r_test = np.sqrt((x[i] - location_x) ** 2 + (y[j] - location_y) ** 2)
                        if r_test == 0.0:
                            if tree.trunk_diameter <= self.pixel_size:
                                bad_local[k, j, i] = radius**2 * pi
                            else:
                                # WORKAROUND: divide remaining circle area over the 8 surrounding
                                # valid_pixels
                                bad_local[k, j - 1 : j + 2, i - 1 : i + 2] = radius**2 * pi / 8.0
                                # for the central pixel fill the pixel
                                bad_local[k, j, i] = self.pixel_size**2
                        # elif ( r_test <= radius ):
                        # TODO: calculate circle segment of grid points cut by the grid

        # Calculate the position of the local 3d tree array within the full
        # domain in order to achieve correct mapping and cutting off at the edges
        # of the full domain
        lad_loc_nx = int(len(x) / 2)
        lad_loc_ny = int(len(y) / 2)
        lad_loc_nz = int(len(z))

        odd_x = int(len(x) % 2)
        odd_y = int(len(y) % 2)

        ind_l_x = max(0, (tree.i - lad_loc_nx))
        ind_l_y = max(0, (tree.j - lad_loc_ny))
        ind_r_x = min(lad_global.shape[2] - 1, tree.i + lad_loc_nx - 1 + odd_x)
        ind_r_y = min(lad_global.shape[1] - 1, tree.j + lad_loc_ny - 1 + odd_y)

        out_l_x = ind_l_x - (tree.i - lad_loc_nx)
        out_l_y = ind_l_y - (tree.j - lad_loc_ny)
        out_r_x = len(x) - 1 + ind_r_x - (tree.i + lad_loc_nx - 1 + odd_x)
        out_r_y = len(y) - 1 + ind_r_y - (tree.j + lad_loc_ny - 1 + odd_y)

        is_tree = ~ma.getmaskarray(lad_local) | ~ma.getmaskarray(bad_local)

        lad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            ~ma.getmaskarray(lad_local)[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            lad_local[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            lad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
        bad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            ~ma.getmaskarray(bad_local)[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            bad_local[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            bad_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
        id_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            is_tree[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            tree.id,
            id_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
        type_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1] = ma.where(
            is_tree[0:lad_loc_nz, out_l_y : out_r_y + 1, out_l_x : out_r_x + 1],
            tree.type,
            type_global[0:lad_loc_nz, ind_l_y : ind_r_y + 1, ind_l_x : ind_r_x + 1],
        )
