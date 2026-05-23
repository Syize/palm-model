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

"""Test the GeoConverter class."""

import itertools
import os
from pathlib import Path
from typing import Generator, Tuple

import pytest
import rasterio as rio
import rasterio.warp as riowp

from palm_csd.csd_config import CSDConfigDomain, CSDConfigOutput, CSDConfigSettings
from palm_csd.geo_converter import GeoConverter
from tests.tools import geotiff_equal

test_folder = Path("tests/05_geo_converter/")
test_folder_input = test_folder / "input/"
test_folder_output = test_folder / "output/"


@pytest.fixture(scope="module")
def configs_converter(
    request,
) -> Generator[
    Tuple[GeoConverter, CSDConfigDomain, CSDConfigSettings, CSDConfigOutput], None, None
]:
    """Create configuration and GeoConverter objects.

    Args:
        request: Test mode and rotation angle.

    Raises:
        ValueError: Unknown test mode.

    Yields:
        GeoConverter, domain configuration, settings configuration, output configuration.
    """
    mode = request.param[0]  # relation of input data to target grid
    rotation_angle = request.param[1]

    if mode in ["aligned", "wgs84", "downscaling", "upscaling"]:
        # aligned: aligned input and target grid
        # wgs84: different projection of input grid
        # downscaling: input grid with coarser resolution
        # upscaling: input grid with finer resolution
        origin_x = 389891.5
        origin_y = 5819849.5
    elif mode == "shifted":
        # shifted: input grid shifted relative to target grid
        origin_x = 389892.5
        origin_y = 5819850.5
    else:
        raise ValueError(f"Unknown mode {mode}")

    if mode == "upscaling":
        pixel_size = 10
        nx = 19
        ny = 29
        dz = 10
    else:
        pixel_size = 3
        nx = 69
        ny = 94
        dz = 3

    settings_config = CSDConfigSettings(epsg=25833, rotation_angle=rotation_angle)
    domain_config = CSDConfigDomain(
        pixel_size=pixel_size,
        nx=nx,
        ny=ny,
        dz=dz,
        origin_x=origin_x,
        origin_y=origin_y,
    )
    output_config = CSDConfigOutput(
        path=test_folder,
        file_out=Path("static_driver"),
    )
    gc = GeoConverter(
        domain_config,
        settings_config,
        output_config,
        debug_output=True,
        domain_name=f"{mode}-{rotation_angle}",
    )
    yield gc, domain_config, settings_config, output_config
    CSDConfigDomain._reset_counter()
    CSDConfigSettings._reset_counter()
    CSDConfigOutput._reset_counter()


mode = ["aligned"]  # see configs_converter for meaning
angles = [0, 30, 165, 200, 320]  # rotation angle
combinations = list(itertools.product(mode, angles))
names = [x + " " + str(y) for x, y in combinations]


@pytest.mark.parametrize(
    "configs_converter",
    combinations,
    ids=names,
    indirect=True,
)
def test_geo_converter_attributes(
    configs_converter: Tuple[GeoConverter, CSDConfigDomain, CSDConfigSettings, CSDConfigOutput],
):
    """Check the attributes of the GeoConverter object.

    Args:
        configs_converter: GeoConverter, domain configuration, settings configuration, output
          configuration.
    """
    gc = configs_converter[0]
    domain_config = configs_converter[1]
    settings_config = configs_converter[2]

    # check simple attributes
    assert gc.pixel_size == domain_config.pixel_size
    assert gc.rotation_angle == settings_config.rotation_angle

    assert gc.dst_width == domain_config.nx + 1
    assert gc.dst_height == domain_config.ny + 1

    assert gc.lower_left_x == 0
    assert gc.lower_left_y == 0

    assert gc.origin_x == domain_config.origin_x
    assert gc.origin_y == domain_config.origin_y

    # coordinate lat lon
    assert gc.origin_lon == pytest.approx(13.377263228261443, abs=1e-14)
    assert gc.origin_lat == pytest.approx(52.51762092679638, abs=1e-14)

    # origin_x, origin_y should not be rotated
    x, y = rio.transform.xy(gc.dst_transform, gc.dst_height - 1, 0, offset="ll")
    assert x == gc.origin_x and y == gc.origin_y


def _read_project_check(
    gc: GeoConverter,
    variable: str,
    resampling_downscaling: riowp.Resampling,
    resampling_upscaling: riowp.Resampling,
    compatibility_resampling_downscaling: riowp.Resampling,
    compatibility_resampling_upscaling: riowp.Resampling,
):
    """Read, project and check a GeoTIFF file.

    Args:
        gc: GeoConverter object.
        variable: Variable name.
        resampling_downscaling: Resampling method for downscaling.
        resampling_upscaling: Resampling method for upscaling.
        compatibility_resampling_downscaling: Masked values of this resampling method should be
            applied to the output when downscaling.
        compatibility_resampling_upscaling: Masked values of this resampling method should be
            applied to the output when upscaling.

    Raises:
        ValueError: Undefined domain name.
    """
    if gc.domain_name is None:
        raise ValueError("Domain name is None.")
    if gc.domain_name.startswith("wgs84"):
        file_input = test_folder_input / f"Berlin_{variable}_3m_DLR_WGS84.tif"
        file_output = f"static_driver_{variable}-reprojected_{gc.domain_name}.tif"
    elif gc.domain_name.startswith("downscaling"):
        file_input = test_folder_input / f"Berlin_{variable}_15m_DLR.tif"
        file_output = f"static_driver_{variable}-reprojected_{gc.domain_name}.tif"
    elif gc.domain_name.startswith("upscaling") or gc.domain_name.startswith("shifted"):
        file_input = test_folder_input / f"Berlin_{variable}_3m_DLR.tif"
        file_output = f"static_driver_{variable}-reprojected_{gc.domain_name}.tif"
    else:
        file_input = test_folder_input / f"Berlin_{variable}_3m_DLR.tif"
        if gc.rotation_angle == 0:
            file_output = f"static_driver_{variable}-cut_{gc.domain_name}.tif"
        else:
            file_output = f"static_driver_{variable}-reprojected_{gc.domain_name}.tif"

    gc.read_to_dst(
        file_input,
        name=variable,
        resampling_downscaling=resampling_downscaling,
        resampling_upscaling=resampling_upscaling,
        compatibility_resampling_downscaling=compatibility_resampling_downscaling,
        compatibility_resampling_upscaling=compatibility_resampling_upscaling,
    )
    assert geotiff_equal(test_folder_output / file_output, test_folder / file_output)
    os.remove(test_folder / file_output)


mode = [
    "aligned",
    "shifted",
    "wgs84",
    "downscaling",
    "upscaling",
]  # see configs_converter for meaning
angles = [0, 165]  # rotation angle
combinations = list(itertools.product(mode, angles))
names = [x + " " + str(y) for x, y in combinations]


@pytest.mark.parametrize(
    "configs_converter",
    combinations,
    ids=names,
    indirect=True,
)
def test_geo_converter_transform(
    configs_converter: Tuple[GeoConverter, CSDConfigDomain, CSDConfigSettings, CSDConfigOutput],
):
    """Check the reading and transformation of GeoTIFF files.

    Args:
        configs_converter: GeoConverter, domain configuration, settings configuration, output
          configuration.
    """
    gc = configs_converter[0]

    _read_project_check(
        gc,
        "terrain_height",
        resampling_downscaling=riowp.Resampling.bilinear,
        resampling_upscaling=riowp.Resampling.average,
        compatibility_resampling_downscaling=riowp.Resampling.nearest,
        compatibility_resampling_upscaling=riowp.Resampling.nearest,
    )

    _read_project_check(
        gc,
        "building_height",
        resampling_downscaling=riowp.Resampling.nearest,
        resampling_upscaling=riowp.Resampling.average,
        compatibility_resampling_downscaling=riowp.Resampling.nearest,
        compatibility_resampling_upscaling=riowp.Resampling.nearest,
    )

    _read_project_check(
        gc,
        "building_type",
        resampling_downscaling=riowp.Resampling.nearest,
        resampling_upscaling=riowp.Resampling.mode,
        compatibility_resampling_downscaling=riowp.Resampling.nearest,
        compatibility_resampling_upscaling=riowp.Resampling.nearest,
    )

    _read_project_check(
        gc,
        "tree_height",
        resampling_downscaling=riowp.Resampling.nearest,
        resampling_upscaling=riowp.Resampling.nearest,
        compatibility_resampling_downscaling=riowp.Resampling.nearest,
        compatibility_resampling_upscaling=riowp.Resampling.nearest,
    )
