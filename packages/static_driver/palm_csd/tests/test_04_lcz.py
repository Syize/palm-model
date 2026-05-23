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

"""Local Climate Zone (LCZ) related tests."""

from math import sqrt
from typing import Tuple

import numpy as np
import pytest

from palm_csd.lcz import LCZTypes


@pytest.fixture(scope="module")
def lcz_types_geometric_arithmetic(request) -> Tuple[LCZTypes, bool]:
    """Create LCZTypes object for summer and geometric or arithmetic averaging."""
    lcz_types = LCZTypes("summer", request.param)
    return lcz_types, request.param


@pytest.mark.parametrize("lcz_types_geometric_arithmetic", [True, False], indirect=True)
def test_lcz_types(lcz_types_geometric_arithmetic: Tuple[LCZTypes, bool]):
    """Test that LCZTypes object is created correctly."""
    lcz_types = lcz_types_geometric_arithmetic[0]
    geometric_mean = lcz_types_geometric_arithmetic[1]
    if geometric_mean:
        # averaging function is geometric mean
        assert lcz_types._mean_height([2, 3]) == sqrt(2 * 3)
    else:
        # averaging function is arithmetic mean
        assert lcz_types._mean_height([2, 3]) == 2.5

    # all urban LCZ should have default height as average of min and max
    for lcz in lcz_types.urban_like:
        assert lcz.height_roughness_elements.default == lcz_types._mean_height(
            [lcz.height_roughness_elements.minimum, lcz.height_roughness_elements.maximum]
        )


@pytest.mark.parametrize("lcz_types_geometric_arithmetic", [True, False], indirect=True)
def test_building_heights(lcz_types_geometric_arithmetic: Tuple[LCZTypes, bool]):
    """Check height distribution of buildings of all LCZ classes."""
    lcz_types = lcz_types_geometric_arithmetic[0]

    u_dir = [0.0, 90.0]
    z_uhl = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0]
    z_uhl_array = np.array(z_uhl)

    # Half of layer thickness
    dz_uhl = (z_uhl_array[1:] - z_uhl_array[:-1]) / 2.0
    # Assume same layer thickness above last height
    dz_uhl = np.append(dz_uhl, dz_uhl[-1])
    # Urban full layers
    z_ufl = z_uhl + dz_uhl

    for lcz in lcz_types.index.values():
        height_distribution = lcz_types.building_height_from_lcz(lcz, z_uhl, u_dir)
        if lcz in lcz_types.urban_like:
            assert not height_distribution.mask.all()
            assert height_distribution.shape == (1, len(u_dir), len(z_uhl))
            for d in range(len(u_dir)):
                assert height_distribution[0, d, :].sum() == pytest.approx(1.0)
            # default behaviour: values between min and max, 0.0 otherwise
            #                    largest value at default height
            assert lcz.height_roughness_elements.minimum is not None
            assert lcz.height_roughness_elements.maximum is not None
            assert lcz.height_roughness_elements.default is not None
            index_min = np.searchsorted(z_ufl, lcz.height_roughness_elements.minimum)
            index_max = np.searchsorted(z_ufl, lcz.height_roughness_elements.maximum)
            index_mean = np.searchsorted(z_ufl, lcz.height_roughness_elements.default)
            assert np.all(height_distribution[0, :, 0:index_min] == 0.0)
            assert np.all(height_distribution[0, :, (index_max + 1) :] == 0.0)
            assert np.all(height_distribution[0, :, index_min : (index_max + 1)] > 0.0)
            for d in range(len(u_dir)):
                assert np.all(height_distribution[0, d, index_mean] >= height_distribution[0, d, :])
        else:
            assert height_distribution.mask.all()
