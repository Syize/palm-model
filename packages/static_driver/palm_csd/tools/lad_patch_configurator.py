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

"""Visualize the height-dependence of the Leaf Area Density (LAD) for different formulations."""

import sys
from pathlib import Path

import hvplot.xarray
import numpy as np
import numpy.ma as ma
import xarray as xr

sys.path.append(str(Path(__file__).parent.parent))
from palm_csd.vegetation import CanopyGenerator  # noqa: E402

# input values, most of them are not important for the plot
dz = 0.2
pixel_size = 1.0
patch_height = ma.ones((1, 1)) * 20.0
patch_lai = ma.ones((1, 1)) * 8.0

patch_type = ma.ones((1, 1)) * 14

z_rel = np.linspace(0, 1.0, int(patch_height[0, 0] / dz) + 1)
z_rel_centre = 0.5 * (z_rel[1:] + z_rel[:-1])


lad_Metal2002 = xr.DataArray(
    dims=["z_rel", "alpha", "beta"],
    # data=np.ones((100, 100, 100)),
    coords={
        "z_rel": z_rel_centre,
        "alpha": np.arange(1.0, 20.0, 0.5),
        "beta": np.arange(1.0, 20.0, 0.5),
    },
)

# prepare data for different alpha and beta values
for alpha in lad_Metal2002.alpha.values:
    for beta in lad_Metal2002.beta.values:
        canopy_generator = CanopyGenerator(
            method="Metal2003",
            alpha_Metal2003=alpha,
            beta_Metal2003=beta,
            dz=dz,
            pixel_size=pixel_size,
        )
        lad_tmp, *_ = canopy_generator.process_patch(patch_height, patch_type, patch_lai)

        lad_Metal2002.loc[dict(alpha=alpha, beta=beta)] = (
            lad_tmp[1:, 0, 0] * patch_height[0, 0] / patch_lai[0, 0]
        )


lad_LM2004 = xr.DataArray(
    dims=["z_rel", "z_max_rel"],
    # data=np.ones((100, 100, 100)),
    coords={
        "z_rel": z_rel_centre,
        "z_max_rel": np.linspace(0.05, 0.95, 19),
    },
)

# prepare data for different alpha and beta values
for z_max_rel in lad_LM2004.z_max_rel.values:
    canopy_generator = CanopyGenerator(
        method="LM2004",
        z_max_rel_LM2004=z_max_rel,
        dz=dz,
        pixel_size=pixel_size,
    )
    lad_tmp, *_ = canopy_generator.process_patch(patch_height, patch_type, patch_lai)

    lad_LM2004.loc[dict(z_max_rel=z_max_rel)] = (
        lad_tmp[1:, 0, 0] * patch_height[0, 0] / patch_lai[0, 0]
    )


# combine data into one dataset
lads = xr.Dataset({"Metal2002": lad_Metal2002, "LM2004": lad_LM2004})

# plot
plot = lads.hvplot(
    x="z_rel",
    groupby=["alpha", "beta", "z_max_rel"],
    fields={"alpha": {"default": 5}, "beta": {"default": 3}, "z_max_rel": {"default": 0.7}},
    xlim=(0, 1.0),
    xlabel="z / h",
    ylabel="LAD / LAI * h",
    title="Height-dependence of the Leaf Area Density (LAD)",
    group_label="Formulation",
)
hvplot.show(plot)
