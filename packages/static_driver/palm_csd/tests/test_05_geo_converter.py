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

"""Test the GeoConverter class."""

import itertools
from pathlib import Path
from typing import Generator, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
import rasterio as rio
import rasterio.warp as riowp
from shapely.geometry import MultiPoint, Polygon
from shapely.geometry import Polygon as ShapelyPolygon

from palm_csd.csd_config import CSDConfigDomain, CSDConfigOutput, CSDConfigSettings
from palm_csd.geo_converter import GeoConverter
from tests.tools import geotiff_equal

test_folder = Path("tests/05_geo_converter/")
test_folder_input = test_folder / "input/"
test_folder_output = test_folder / "output/"

tree_shp = test_folder_input / "trees.shp"
tree_shp_wgs84 = test_folder_input / "trees_wgs84.shp"
tree_shp_z = test_folder_input / "trees_z.shp"
tree_shp_multi = test_folder_input / "trees_multi.shp"
tree_shp_multiz = test_folder_input / "trees_multiz.shp"


@pytest.fixture(scope="function")
def configs_converter(
    request: pytest.FixtureRequest,
    tmp_path: Path,
) -> Generator[
    Tuple[GeoConverter, CSDConfigDomain, CSDConfigSettings, CSDConfigOutput], None, None
]:
    """Create configuration and GeoConverter objects.

    Args:
        request: Test mode and rotation angle.
        tmp_path: Temporary path for storing the configuration file.

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
        path=tmp_path,
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

    assert gc.corner_ll == (gc.origin_x, gc.origin_y)

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

    gc.read_raster_to_dst(
        file_input,
        name=variable,
        resampling_downscaling=resampling_downscaling,
        resampling_upscaling=resampling_upscaling,
        compatibility_resampling_downscaling=compatibility_resampling_downscaling,
        compatibility_resampling_upscaling=compatibility_resampling_upscaling,
    )
    assert geotiff_equal(
        gc.debug_file_prefix.parent / file_output, test_folder_output / file_output
    )


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


# Tree test data
pixel_size = 1.0
epsg = 25833
epsg_wgs84 = 4326
nx = 54
ny = 19
origin_x = 400000.0
origin_y = 5800000.0

trees = pd.DataFrame(
    {
        "Kronend": [5.5, np.nan, 7.2, 8.1, 5.9, 6.3],
        "Hoehe": [9.2, 13.5, 14.8, np.nan, 15.9, 16.3],
        "Form": [pd.NA, 2, 3, 4, 5, 6],
        "Typ": [10, 20, pd.NA, pd.NA, pd.NA, 60],
        "Spezies": [pd.NA, "Cladrastis", "Elaeagnus whatever", "Juglans etc", "", ""],
        "x_index": [10, 20, 30, 40, 50, 60],
        "y_index": [10, 10, 10, 10, 10, 10],
        "z": [0, 2, 4, 6, 8, 10],
    }
)
trees["Form"] = trees["Form"].astype("Int32")
trees["Typ"] = trees["Typ"].astype("Int32")
trees["Spezies"] = trees["Spezies"].astype("string")

trees_multi = pd.DataFrame(
    {
        "Kronend": [5.5],
        "Hoehe": [9.2],
        "Form": [3],
        "Typ": [10],
        "Spezies": ["Cladrastis"],
    }
)
trees_multi["Spezies"] = trees_multi["Spezies"].astype("string")


@pytest.fixture(scope="function")
def configs_converter_trees() -> Generator[
    Tuple[GeoConverter, CSDConfigDomain, CSDConfigSettings, CSDConfigOutput], None, None
]:
    """Create configuration and GeoConverter objects.

    Yields:
        GeoConverter, domain configuration, settings configuration, output configuration.
    """
    settings_config = CSDConfigSettings(epsg=epsg)
    domain_config = CSDConfigDomain(
        pixel_size=pixel_size,
        nx=nx,
        ny=ny,
        dz=1.0,
        origin_x=origin_x,
        origin_y=origin_y,
    )
    output_config = CSDConfigOutput(
        path=Path(test_folder),
        file_out=Path("static_driver"),
    )
    gc = GeoConverter(
        domain_config,
        settings_config,
        output_config,
        debug_output=False,
        domain_name="domain",
    )
    yield gc, domain_config, settings_config, output_config
    CSDConfigDomain._reset_counter()
    CSDConfigSettings._reset_counter()
    CSDConfigOutput._reset_counter()


def create_tree_shape() -> None:
    """Create a shapefile with trees."""
    # Create shapefile with trees with x, y coordinates
    trees_geo = trees.copy()
    trees_geo["geometry"] = gpd.points_from_xy(
        (trees_geo["x_index"] + 0.5) * pixel_size + origin_x,
        (trees_geo["y_index"] + 0.5) * pixel_size + origin_y,
    )
    trees_geo = trees_geo.drop(columns=["x_index", "y_index", "z"])
    trees_geo = gpd.GeoDataFrame(data=trees_geo, crs=f"EPSG:{epsg}")
    trees_geo.to_file(tree_shp)

    # Create shapefile with trees in WGS84
    trees_geo.to_crs(epsg_wgs84).to_file(tree_shp_wgs84)

    # Create shapefile with trees with z coordinate
    trees_geo = trees.copy()
    trees_geo["geometry"] = gpd.points_from_xy(
        (trees_geo["x_index"] + 0.5) * pixel_size + origin_x,
        (trees_geo["y_index"] + 0.5) * pixel_size + origin_y,
        trees_geo["z"],
    )
    trees_geo = trees_geo.drop(columns=["x_index", "y_index", "z"])
    trees_geo = gpd.GeoDataFrame(data=trees_geo, crs=f"EPSG:{epsg}")
    trees_geo.to_file(tree_shp_z)

    # Create shapefile with trees with multipoint geometry
    trees_points = [
        (
            (row["x_index"] + 0.5) * pixel_size + origin_x,
            (row["y_index"] + 0.5) * pixel_size + origin_y,
        )
        for _, row in trees.iterrows()
    ]
    multipoint = MultiPoint(trees_points)
    # Create a GeoDataFrame with a single row containing the MultiPoint
    trees_multi_geo = gpd.GeoDataFrame(data=trees_multi, geometry=[multipoint], crs=f"EPSG:{epsg}")
    trees_multi_geo.to_file(tree_shp_multi)

    # Create shapefile with trees with multipoint geometry and z coordinate
    trees_points_z = [
        (
            (row["x_index"] + 0.5) * pixel_size + origin_x,
            (row["y_index"] + 0.5) * pixel_size + origin_y,
            10.0,
        )
        for _, row in trees.iterrows()
    ]
    multipoint_z = MultiPoint(trees_points_z)
    # Create a GeoDataFrame with a single row containing the MultiPoint
    trees_multiz_geo = gpd.GeoDataFrame(
        data=trees_multi, geometry=[multipoint_z], crs=f"EPSG:{epsg}"
    )
    trees_multiz_geo.to_file(tree_shp_multiz)


# Uncomment to create test trees data.
# create_tree_shape()


@pytest.mark.parametrize(
    "tree_shp",
    [tree_shp, tree_shp_wgs84, tree_shp_z],
    ids=["shp", "shp_wgs84", "shp_z"],
)
def test_reading_trees(configs_converter_trees, tree_shp):
    """Test reading of tree data.

    Args:
        configs_converter_trees: GeoConverter, domain configuration, settings configuration, output
        tree_shp: Path to tree shapefile.
    """
    gc = configs_converter_trees[0]
    # filter trees outside of domain
    trees_filtered = trees[trees["x_index"] <= nx].drop(columns="z")

    trees_read = gc.read_shp_to_dst(tree_shp, shape_type="Point").drop(columns="geometry")
    # For comparision, convert to correct dtype
    trees_read["Form"] = trees_read["Form"].astype("Int32")
    trees_read["Typ"] = trees_read["Typ"].astype("Int32")
    trees_read["Spezies"] = trees_read["Spezies"].astype("string")
    # Original "" is converted to NaN
    assert pd.isna(trees_read["Spezies"][4])
    trees_read.loc[4, "Spezies"] = ""
    pdt.assert_frame_equal(
        trees_read.sort_index(axis="columns"),
        trees_filtered.sort_index(axis="columns"),
        check_dtype=False,
    )


@pytest.mark.parametrize(
    "tree_shp", [tree_shp_multi, tree_shp_multiz], ids=["shp_multi", "shp_multiz"]
)
def test_reading_trees_multi(configs_converter_trees, tree_shp):
    """Test reading of tree data with multipoints.

    Args:
        configs_converter_trees: GeoConverter, domain configuration, settings configuration, output
        tree_shp: Path to tree shapefile.
    """
    gc = configs_converter_trees[0]
    # filter trees outside of domain
    trees_filtered = trees[trees["x_index"] <= nx].drop(columns="z")

    trees_read = gc.read_shp_to_dst(tree_shp, "Point").drop(columns="geometry")
    # For comparision, convert to correct dtype
    trees_read["Spezies"] = trees_read["Spezies"].astype("string")

    pdt.assert_series_equal(
        trees_read.iloc[0, trees_read.columns.get_indexer(trees_multi.columns)], trees_multi.iloc[0]
    )

    pdt.assert_series_equal(
        trees_read["x_index"], trees_filtered["x_index"], check_index=False, check_dtype=False
    )
    pdt.assert_series_equal(
        trees_read["y_index"], trees_filtered["y_index"], check_index=False, check_dtype=False
    )


@pytest.fixture(scope="function")
def polygon_fixture(request, tmp_path):
    """Parametrized fixture for polygon shapefiles in UTM or WGS84."""

    def make_polygon_gdf(crs):
        poly = Polygon([(400010, 5800010), (400020, 5800010), (400020, 5800020), (400010, 5800020)])
        data = {
            "id": [1],
            "name": ["TestPolygon"],
            "value": [42.0],
            "geometry": [poly],
        }
        return gpd.GeoDataFrame(data, crs=crs)

    kind = request.param
    gdf = make_polygon_gdf(f"EPSG:{epsg}")
    if kind == "polygon_shp":
        shp_path = tmp_path / "polygon.shp"
        gdf.to_file(shp_path)
        return shp_path, gdf
    elif kind == "polygon_shp_wgs84":
        gdf_wgs84 = gdf.to_crs(epsg=epsg_wgs84)
        shp_path = tmp_path / "polygon_wgs84.shp"
        gdf_wgs84.to_file(shp_path)
        return shp_path, gdf_wgs84
    else:
        raise ValueError(f"Unknown polygon_fixture param: {kind}")


@pytest.mark.parametrize(
    "polygon_fixture",
    ["polygon_shp", "polygon_shp_wgs84"],
    ids=["utm", "wgs84"],
    indirect=True,
)
def test_reading_polygon(configs_converter_trees, polygon_fixture):
    """Test reading of polygon data, including reprojection from WGS84."""
    gc = configs_converter_trees[0]
    shp_path, original_df = polygon_fixture

    # Read polygon shapefile
    polygons_read = gc.read_shp_to_dst(shp_path, shape_type="Polygon")

    # If WGS84, reproject original to target CRS for comparison
    if str(original_df.crs) != f"EPSG:{epsg}":
        original_df = original_df.to_crs(epsg=epsg)

    # Compare attribute columns (ignore geometry order)
    pd.testing.assert_frame_equal(
        polygons_read.drop(columns="geometry").reset_index(drop=True),
        original_df.drop(columns="geometry").reset_index(drop=True),
        check_dtype=False,
    )
    # Compare geometry ignoring vertex order
    for g1, g2 in zip(polygons_read.geometry, original_df.geometry):
        assert isinstance(g1, ShapelyPolygon) and isinstance(g2, ShapelyPolygon)
        assert g1.equals_exact(g2, tolerance=1e-6) or g1.equals(g2)


@pytest.fixture(scope="function")
def polygon_partial_outside_fixture(request, tmp_path):
    """Parametrized fixture for polygons partly or fully outside the domain."""
    kind = request.param

    def make_polygon_gdf(coords, crs):
        poly = Polygon(coords)
        data = {
            "id": [1],
            "name": ["TestPolygon"],
            "value": [42.0],
            "geometry": [poly],
        }
        return gpd.GeoDataFrame(data, crs=crs)

    # Domain: origin_x=400000, origin_y=5800000, nx=54, ny=19, pixel_size=1.0
    # So domain covers x: 400000-400054, y: 5800000-5800019

    # Polygon partly inside domain
    if kind == "partial":
        coords = [(400050, 5800015), (400060, 5800015), (400060, 5800025), (400050, 5800025)]
        gdf = make_polygon_gdf(coords, f"EPSG:{epsg}")
        shp_path = tmp_path / "polygon_partial.shp"
        gdf.to_file(shp_path)
        return shp_path, gdf

    # Polygon completely outside domain
    elif kind == "outside":
        coords = [(400100, 5800100), (400110, 5800100), (400110, 5800110), (400100, 5800110)]
        gdf = make_polygon_gdf(coords, f"EPSG:{epsg}")
        shp_path = tmp_path / "polygon_outside.shp"
        gdf.to_file(shp_path)
        return shp_path, gdf

    else:
        raise ValueError(f"Unknown polygon_partial_outside_fixture param: {kind}")


@pytest.mark.parametrize(
    "polygon_partial_outside_fixture",
    ["partial", "outside"],
    ids=["partial_in_domain", "outside_domain"],
    indirect=True,
)
def test_reading_polygon_partial_and_outside(
    configs_converter_trees, polygon_partial_outside_fixture
):
    """Test reading of polygons that are partly or fully outside the domain."""
    gc = configs_converter_trees[0]
    shp_path, original_df = polygon_partial_outside_fixture

    polygons_read = gc.read_shp_to_dst(shp_path, shape_type="Polygon")

    if "partial" in shp_path.name:
        # Should be included (geometry may be clipped, but at least one feature present)
        assert len(polygons_read) == 1
    elif "outside" in shp_path.name:
        # Should not be included
        assert len(polygons_read) == 0
