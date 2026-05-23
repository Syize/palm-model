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

"""Support routines for palm_csd."""

import logging
from typing import Any, Generator, Generic, List, Optional, Tuple, Type, TypeVar, Union, cast

import numpy as np
import numpy.ma as ma
import numpy.ma.core as ma_core
import numpy.typing as npt
import pandas as pd
import pandas.api.typing as pdtypes
from pydantic import BaseModel, model_validator
from scipy.interpolate import RectBivariateSpline

# module logger
logger = logging.getLogger(__name__)


def blend_array_2d(
    array1: Union[npt.NDArray, ma.MaskedArray],
    array2: Union[npt.NDArray, ma.MaskedArray],
    radius: float,
) -> Union[npt.NDArray, ma.MaskedArray]:
    """Blend over the parent and child terrain height within a given radius.

    Args:
        array1: Array 1.
        array2: Array 2.
        radius: Radius of the blending.

    Returns:
        Blended array.
    """
    gradient_matrix = np.ones(array1.shape)
    radius = int(radius)

    for j in range(0, radius):
        gradient_matrix[:, j] = float(j) / float(radius)
        gradient_matrix[:, -j - 1] = float(j) / float(radius)
        gradient_matrix[j, :] = float(j) / float(radius)
        gradient_matrix[-j - 1, :] = float(j) / float(radius)

    for j in range(0, radius):
        for i in range(0, radius):
            gradient_matrix[j, i] = max(
                1 - np.sqrt((i - (0.0 + radius)) ** 2 + (j - (0.0 + radius)) ** 2) / radius, 0
            )
            gradient_matrix[-j - 1, i] = max(
                1 - np.sqrt((i - (0.0 + radius)) ** 2 + (j - (0.0 + radius)) ** 2) / radius, 0
            )
            gradient_matrix[j, -i - 1] = max(
                1 - np.sqrt((i - (0.0 + radius)) ** 2 + (j - (0.0 + radius)) ** 2) / radius, 0
            )
            gradient_matrix[-j - 1, -i - 1] = max(
                1 - np.sqrt((i - (0.0 + radius)) ** 2 + (j - (0.0 + radius)) ** 2) / radius, 0
            )

    array_blended = array1 * gradient_matrix + (1.0 - gradient_matrix) * array2

    return array_blended


def interpolate_2d(
    array: npt.NDArray, x1: npt.NDArray, y1: npt.NDArray, x2: npt.NDArray, y2: npt.NDArray
):
    """Linearly interpolate array(x1, y1) to array(x2, y2).

    Args:
        array: Array with original values defined for x1 and y1.
        x1: x-coordinates of input array.
        y1: y-coordinates of input array.
        x2: Target x-coordinates.
        y2: Target y-coordinates.

    Returns:
        Array with interpolated values defined for x2 and y2.
    """
    # Create interpolation object that approximates f(x,y)=array(x1, y1). Based on the deprecated
    # interp2d.
    tmp_int2d = RectBivariateSpline(x1, y1, array.T, kx=1, ky=1)
    # Apply f to x2 and y2.
    array_ip = tmp_int2d(x2.astype(float), y2.astype(float)).T

    # Round values to avoid numerical issues resulting from 7.499999999999999 vs. 7.5.
    return np.round(array_ip, 14)


def height_to_z_grid(array: npt.NDArray, dz: float) -> npt.NDArray:
    """Discretize height array to z grid defined by dz.

    Args:
        array: Input array.
        dz: Height of grid cells.

    Raises:
        ValueError: Negative values in array.

    Returns:
        Array with discretized heights.
    """
    if np.any(array < 0):
        raise ValueError("All array values need to be larger or equal 0")

    # z grid starting at zero.
    k_tmp = np.arange(0, max(array.flatten()) + dz * 2, dz)
    k_tmp[1:] = k_tmp[1:] - dz * 0.5

    # Index of k_tmp with smaller value than array.
    index_smaller = np.searchsorted(k_tmp, array, side="right") - 1  # "right" to for height 0

    # Return height smaller than array plus half grid cell.
    return k_tmp[index_smaller] + dz * 0.5


def check_consistency_3(
    array1: ma.MaskedArray, array2: ma.MaskedArray, array3: ma.MaskedArray
) -> Tuple[npt.NDArray, np.bool_]:
    """Check consistency.

    Args:
        array1: Array 1.
        array2: Array 2.
        array3: Array 3.

    Returns:
        Array with the number of unmasked arrays at each point and a boolean indicating if the
        arrays are consistent.
    """
    # Todo: is -1 for array3 correct?
    tmp_array = (
        np.where(ma.getmaskarray(array1), 0, 1)
        + np.where(ma.getmaskarray(array2), 0, 1)
        + np.where(ma.getmaskarray(array3), 0, -1)
    )

    test = np.any(tmp_array != 0)
    if test:
        logger.warning("soil_type array is not consistent!")
        logger.warning(
            "max: " + str(max(tmp_array.flatten())) + ", min: " + str(min(tmp_array.flatten()))
        )
    else:
        logger.debug("soil_type array is consistent!")
    return tmp_array, test


# Check if at every point only one of the arrays is not masked
def check_consistency_4(
    array1: ma.MaskedArray, array2: ma.MaskedArray, array3: ma.MaskedArray, array4: ma.MaskedArray
) -> Tuple[npt.NDArray, np.bool_]:
    """Check if at every point only one of the arrays is not masked.

    Args:
        array1: Input array 1.
        array2: Input array 2.
        array3: Input array 3.
        array4: Input array 4.

    Returns:
        Array with the number of unmasked arrays at each point and a boolean indicating if the
        arrays are consistent.
    """
    tmp_array = (
        np.where(ma.getmaskarray(array1), 0, 1)
        + np.where(ma.getmaskarray(array2), 0, 1)
        + np.where(ma.getmaskarray(array3), 0, 1)
        + np.where(ma.getmaskarray(array4), 0, 1)
    )

    test = np.any(tmp_array != 1)
    if test:
        logger.warning("*_type arrays are not consistent!")
        logger.warning(f"max: {max(tmp_array.flatten())}, min: {min(tmp_array.flatten())}")
    else:
        logger.debug("*_type arrays are consistent!")
    return tmp_array, test


def is_missing(value: Union[float, pdtypes.NAType, ma_core.MaskedConstant]) -> bool:
    """Check if a value is missing.

    Args:
        value: Value to check.

    Returns:
        True if the value is missing, False otherwise.
    """
    # pd.isna returns bool for scalar value so type should be fine.
    return ma.is_masked(value) or pd.isna(value)  # type: ignore


# TODO: use just ma.isin() once the bug in https://github.com/numpy/numpy/issues/19877 or
#  https://stackoverflow.com/questions/69160969/ is fixed
def ma_isin(array: ma.MaskedArray, comparison: npt.ArrayLike) -> ma.MaskedArray:
    """Check for each element of the input array if it is in comparison.

    If a value in the input array is masked, the result is masked as well.

    Args:
        array: Input array to check.
        comparison: Comparison values.

    Returns:
        Boolean masked array with the same shape as the input array.
    """
    return ma.MaskedArray(data=np.isin(array, comparison), mask=array.mask.copy())


class DefaultMinMax(BaseModel, validate_assignment=True):
    """Class to store a default value together with minimum and maximum values.

    It is ensured that the default values is between the minimum and maximum.
    """

    minimum: Optional[Union[float, int]]
    """Minimum value and lower boundary for the default value."""
    maximum: Optional[Union[float, int]]
    """Maximum value and upper boundary for the default value."""
    default: Optional[Union[float, int]]
    """Default value."""

    @model_validator(mode="before")
    @classmethod
    def _adapt_number_type(cls, data: Any) -> Any:
        """Convert input to int or float if possible. Prefer int.

        Args:
            data: Input values. Assumed to be a dictionary.

        Returns:
            Dictionary with values converted to int or float if possible.
        """

        def is_int_like(s: Optional[Any]) -> bool:
            """Check if input is None or representable by an integer.

            Args:
                s: Input value.

            Raises:
                ValueError: Input value not None and not a number.

            Returns:
                True if input is None, "None", a string that can be convertetd to an integer or an
                integer, False otherwise.
            """
            if s is None:
                return True
            if isinstance(s, str):
                if s.strip() == "None":
                    return True
                try:
                    float(s)
                except ValueError:
                    raise ValueError(f"{s} is not None and is no number.")
                return "." not in s

            return isinstance(s, int)

        def input_to_type(
            number: Optional[Any],
            number_type: Union[Type[int], Type[float]] = float,
        ) -> Optional[float]:
            """Convert input to number_type or None.

            Args:
                number: Input value.
                number_type: int or float. Defaults to float.

            Raises:
                ValueError: Input value not convertible to number_type.

            Returns:
                Converted input value according to number_type or None.
            """
            if number is None or number == "None":
                return None
            try:
                result = number_type(number)
            except TypeError:
                raise ValueError(f"Input must be {number_type} or convertible to {number_type}")
            if np.isnan(result):
                return None
            return result

        if not isinstance(data, dict):
            raise ValueError("Expected a dictionary as input.")

        # If all values are int like, use int to store the values, otherwise float
        number_type: Union[Type[int], Type[float]]
        if all(is_int_like(value) for value in data.values()):
            number_type = int
        else:
            number_type = float
        for key, value in data.items():
            data[key] = input_to_type(value, number_type)

        return data

    @model_validator(mode="after")
    def _check_min_max(self: "DefaultMinMax") -> "DefaultMinMax":
        """Check if minimum <= default <= maximum.

        Args:
            self: Instance of DefaultMinMax.

        Raises:
            ValueError: Not minimum <= default <= maximum.

        Returns:
            Validated instance of DefaultMinMax.
        """
        if self.minimum is not None:
            if self.default is not None:
                if not self.minimum <= self.default:
                    raise ValueError("minimum is greater than default.")
            if self.maximum is not None:
                if not self.minimum <= self.maximum:
                    raise ValueError("minimum is larger than maximum.")
        if self.maximum is not None:
            if self.default is not None:
                if not self.default <= self.maximum:
                    raise ValueError("maximum is smaller than default.")
        return self


T = TypeVar("T", bound="Node")


class Node(Generic[T]):
    """Node of a tree structure."""

    def __init__(self, parent: Optional[T] = None):
        """Create a new node.

        Args:
            parent: Parent node. Defaults to None.
        """
        self._parent = parent
        self._children: List[T] = []

        if parent is not None:
            parent.add_child(self)

    def add_child(self, child: T) -> None:
        """Add a child to the node.

        Args:
            child: Child node.
        """
        self._children.append(child)
        child._parent = self

    def remove_child(self, child: T) -> None:
        """Remove a child from the node.

        Args:
            child: Child node.
        """
        self._children.remove(child)
        child._parent = None

    def get_children(self) -> List[T]:
        """Get the children of the node.

        Returns:
            List of children.
        """
        return self._children

    def get_parent(self) -> Optional[T]:
        """Get the parent of the node.

        Returns:
            Parent node or None if the node has no parent.
        """
        return self._parent

    def find_root(self) -> T:
        """Find the root node of the tree.

        Returns:
            Root node of the tree.
        """
        node = cast(T, self)
        while node is not None and node._parent is not None:
            node = node._parent
        return node

    def traverse(self) -> Generator[T, None, None]:
        """Traverse the tree starting from the current node.

        Yields:
            Nodes in the tree starting from the current node.
        """
        yield cast(T, self)
        for child in self._children:
            yield from child.traverse()
