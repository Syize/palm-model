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

"""Run tests for vegetation patches."""

import itertools
from typing import Tuple

import numpy as np
import numpy.ma as ma
import numpy.ma.core as mac
import pytest

from palm_csd.vegetation import CanopyGenerator


@pytest.fixture(scope="module")
def patch_input_fields() -> Tuple[ma.MaskedArray, ma.MaskedArray, ma.MaskedArray]:
    """Vegetation patch input fields.

    Returns:
        Vegetation patch height, type, and LAI.
    """
    rng = np.random.default_rng(12345)
    nx = 4
    ny = 3

    patch_height = ma.masked_array(np.arange(0, nx * ny * 2.0, 2.0) + rng.random(nx * ny)).reshape(
        (ny, nx)
    )
    # ma.ones_like works but is currently not included in stub file. Use mac.ones_like for now to
    # avoid mypy error.
    patch_type = mac.ones_like(patch_height) * 14
    patch_lai = ma.masked_array(np.arange(0, nx * ny * 1.0, 1.0) + rng.random(nx * ny)).reshape(
        (ny, nx)
    )

    return patch_height, patch_type, patch_lai


@pytest.mark.parametrize("z_max_rel", np.round(np.arange(1, 20) * 0.05, 2))
def test_process_patch_LM2004(patch_input_fields, z_max_rel):
    """Test the process_patch method with the LM2004 method.

    Args:
        patch_input_fields: Trigger the fixture to provide the input fields.
        z_max_rel: Relative maximum height of the canopy.
    """
    dz = 2.0
    canopy_generator = CanopyGenerator(
        method="LM2004", z_max_rel_LM2004=z_max_rel, dz=dz, pixel_size=1.0
    )

    patch_height = patch_input_fields[0]
    patch_type = patch_input_fields[1]
    patch_lai = patch_input_fields[2]

    lad_3d, patch_id_3d, patch_type_3d = canopy_generator.process_patch(
        patch_height, patch_type, patch_lai
    )

    np.testing.assert_allclose(ma.sum(lad_3d * dz, axis=0), patch_lai)
    assert ma.all(lad_3d >= 0.0)
    assert ma.all(lad_3d.mask == patch_id_3d.mask)
    assert ma.all(lad_3d.mask == patch_type_3d.mask)


alphas = np.arange(1.0, 20.0, 5.0)
betas = np.arange(1.0, 10.0, 5.0)
combinations = list(itertools.product(alphas, betas))
names = [str(x) + " " + str(y) for x, y in combinations]


@pytest.mark.parametrize("alpha,beta", combinations, ids=names)
def test_process_patch_Metal2002(patch_input_fields, alpha, beta):
    """Test the process_patch method with the Metal2002 method.

    Args:
        patch_input_fields: Trigger the fixture to provide the input fields.
        alpha: Alpha parameter.
        beta: Beta parameter.
    """
    dz = 2.0
    canopy_generator = CanopyGenerator(
        method="Metal2003", alpha_Metal2003=alpha, beta_Metal2003=beta, dz=dz, pixel_size=1.0
    )

    patch_height = patch_input_fields[0]
    patch_type = patch_input_fields[1]
    patch_lai = patch_input_fields[2]

    lad_3d, patch_id_3d, patch_type_3d = canopy_generator.process_patch(
        patch_height, patch_type, patch_lai
    )

    np.testing.assert_allclose(ma.sum(lad_3d * dz, axis=0), patch_lai)
    assert ma.all(lad_3d >= 0.0)
    assert ma.all(lad_3d.mask == patch_id_3d.mask)
    assert ma.all(lad_3d.mask == patch_type_3d.mask)
