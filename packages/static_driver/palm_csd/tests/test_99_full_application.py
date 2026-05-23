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

"""Run the Berlin test cases."""

import itertools
import os
import shutil
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pytest
import rasterio as rio
import rasterio.transform as riotf
from netCDF4 import Dataset

from palm_csd.create_driver import create_driver
from palm_csd.geo_converter import GeoConverter
from tests.tools import (
    add_root_n02,
    add_to_stem,
    modify_configuration,
    modify_configuration_output,
    ncdf_equal,
)

test_folder = Path("tests/99_full_application/")
test_folder_ref = test_folder / "output/"

epsg = 25833
# coordinates of the root domain
origin_x_root = 386891.5
origin_y_root = 5818569.0
origin_lon_root = 13.333505900572572
origin_lat_root = 52.50549956009953
# coordinates of the child domain
# aligned coordinates fit to the root domain
# invalid coordinates are slightly shifted and do not align
# when using the GeoConverter, the invalid coordinates should be corrected to the aligned ones
origin_x_nest_aligned = 389456.5
origin_y_nest_aligned = 5819499.0
origin_lon_nest_aligned = 13.3709716985343
origin_lat_nest_aligned = 52.5143830934953
origin_x_nest_invalid = 389462.5
origin_y_nest_invalid = 5819496.0
origin_lon_nest_invalid = 13.371061076093232
origin_lat_nest_invalid = 52.51435735065383
lower_left_x_nest = 2565.0
lower_left_y_nest = 930.0

# parameters of the input files needed for geotiff conversion
fields_input_data = [
    "bridges_height",
    "bridges_id",
    "bridges_pavement_type",
    "bridges_street_type",
    "building_height",
    "building_id",
    "building_type",
    "leaf_area_index",
    "pavement_type",
    "soil_type",
    "trees_trunk_clean",
    "street_crossings",
    "street_type",
    "terrain_height",
    "trees_age",
    "trees_crown_clean",
    "trees_height_clean",
    "trees_species",
    "trees_trunk_clean",
    "trees_type",
    "vegetation_on_roofs",
    "vegetation_patch_age",
    "vegetation_patch_dbh",
    "vegetation_patch_height",
    "vegetation_patch_type",
    "vegetation_type",
    "water_id",
    "water_type",
]
top_left_x_input_data = {3: 389443 - 1.5, 15: 386749 - 7.5}
top_left_y_input_data = {3: 5820226 + 1.5, 15: 5820676 + 7.5}

fields_debug_output = [
    "bridges_2d",
    "bridges_id",
    "building_id",
    "buildings_2d",
    "building_type",
    "lai",
    "pavement_type",
    "soil_type",
    "street_crossings",
    "street_type",
    "tree_crown_diameter",
    "tree_height",
    "tree_trunk_diameter",
    "tree_type",
    "vegetation_height",
    "vegetation_on_roofs",
    "vegetation_type",
    "water_type",
    "zt",
]

TO_DELETE = List[List[str]]
TO_SET = List[Tuple[List[str], Union[str, float]]]


@pytest.fixture
def configuration(tmp_path: Path) -> Tuple[Path, Path, Path]:
    """Generate a configuration file.

    Args:
        tmp_path: Temporary path for the configuration file and output files.

    Returns:
        Configuration file, output file and reference file.
    """
    config_in = test_folder / "berlin_tiergarten.yml"
    config_out = tmp_path / "berlin_tiergarten.yml"
    file_out = "berlin_tiergarten"
    file_ref = test_folder_ref / "berlin_tiergarten"

    modify_configuration_output(config_in, config_out, tmp_path)

    return config_out, tmp_path / file_out, file_ref


@pytest.fixture
def configuration_geo_referenced(tmp_path: Path) -> Tuple[Path, Path, Path]:
    """Generate a configuration file for geo-referenced data.

    Args:
        tmp_path: Temporary path for the configuration file and output files.

    Returns:
        Configuration file, output file and reference file.
    """
    config_in = test_folder / "berlin_tiergarten_geo_referenced.yml"
    config_out = tmp_path / "berlin_tiergarten_geo_referenced.yml"
    file_out = "berlin_tiergarten_geo_referenced"
    file_ref = test_folder_ref / "berlin_tiergarten_geo_referenced"

    modify_configuration_output(config_in, config_out, tmp_path)

    return config_out, tmp_path / file_out, file_ref


@pytest.fixture
def configuration_trees_input(tmp_path: Path) -> Tuple[Path, Path, Path, Path]:
    """Generate a configuration file with trees input for the nest.

    Args:
        tmp_path: Temporary path for the configuration file and output files.

    Returns:
        Configuration file, output file and reference file.
    """
    config_in = test_folder / "berlin_tiergarten.yml"
    config_out = tmp_path / "berlin_tiergarten_trees.yml"
    file_out = "berlin_tiergarten_trees"
    file_ref = test_folder_ref / "berlin_tiergarten"
    file_ref_diff_nest = test_folder_ref / "diff_berlin_tiergarten_trees_N02"

    # remove coordinate inputs from config
    to_delete = []
    for file_input in ["input_15m", "input_3m"]:
        for file in ["x_utm", "y_utm", "lon", "lat"]:
            to_delete.append((file_input, "files", file))

    # set new values depending on the run parameter
    to_set: TO_SET = [(["output", "file_out"], file_out)]
    to_set.extend(
        [
            (["settings", "epsg"], epsg),
            (["input_3m", "files", "trees"], "Berlin_trees_N02.shp"),
            (["input_3m", "columns", "kronedurch"], "tree_crown_diameter"),
            (["input_3m", "columns", "baumhoehe"], "tree_height"),
            (["input_3m", "columns", "stammdurch"], "tree_trunk_diameter"),
            (["input_3m", "columns", "art_bot"], "tree_type_name"),
            (["domain_root", "origin_x"], origin_x_root),
            (["domain_root", "origin_y"], origin_y_root),
            (["domain_N02", "origin_x"], origin_x_nest_aligned),
            (["domain_N02", "origin_y"], origin_y_nest_aligned),
        ]
    )
    modify_configuration(config_in, config_out, to_delete, to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    return config_out, tmp_path / file_out, file_ref, file_ref_diff_nest


@pytest.fixture
def configuration_no_coordinates_input(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Tuple[Path, Path, Path]:
    """Generate a configuration file with geograpophic parameters but without coordinate files.

    Depending on the request, origin_x/y, origin_lon/lat or epsg_code are added.

    Args:
        request: param attribute with the run parameter.
        tmp_path: Temporary path for the configuration file and output files.

    Raises:
        ValueError: Unknown run parameter.

    Returns:
        Configuration file, output file and reference file.
    """
    run = request.param

    config_in = test_folder / "berlin_tiergarten.yml"
    config_out = tmp_path / f"berlin_tiergarten_no_coordinates_{run}.yml"
    file_out = f"berlin_tiergarten_no_coordinates_{run}"
    file_ref = test_folder_ref / "berlin_tiergarten"

    # remove coordinate inputs from config
    to_delete = []
    for file_input in ["input_15m", "input_3m"]:
        for file in ["x_utm", "y_utm", "lon", "lat"]:
            to_delete.append((file_input, "files", file))

    # set new values depending on the run parameter
    to_set: TO_SET = [(["output", "file_out"], file_out)]
    # to_set.extend([[["settings", "rotation_angle"], rotation_angle]])
    if run == "full":
        to_set.extend(
            [
                (["settings", "epsg"], epsg),
                (["domain_root", "origin_x"], origin_x_root),
                (["domain_root", "origin_y"], origin_y_root),
                (["domain_root", "origin_lon"], origin_lon_root),
                (["domain_root", "origin_lat"], origin_lat_root),
                (["domain_N02", "origin_x"], origin_x_nest_aligned),
                (["domain_N02", "origin_y"], origin_y_nest_aligned),
                (["domain_N02", "origin_lon"], origin_lon_nest_aligned),
                (["domain_N02", "origin_lat"], origin_lat_nest_aligned),
            ]
        )
    elif run == "origin_xy":
        to_set.extend(
            [
                (["settings", "epsg"], epsg),
                (["domain_root", "origin_x"], origin_x_root),
                (["domain_root", "origin_y"], origin_y_root),
                (["domain_N02", "origin_x"], origin_x_nest_aligned),
                (["domain_N02", "origin_y"], origin_y_nest_aligned),
            ]
        )
    elif run == "origin_lonlat":
        to_set.extend(
            [
                (["settings", "epsg"], epsg),
                (["domain_root", "origin_lon"], origin_lon_root),
                (["domain_root", "origin_lat"], origin_lat_root),
                (["domain_N02", "origin_lon"], origin_lon_nest_aligned),
                (["domain_N02", "origin_lat"], origin_lat_nest_aligned),
            ]
        )
    elif run == "no_epsg":
        to_set.extend(
            [
                (["domain_root", "origin_x"], origin_x_root),
                (["domain_root", "origin_y"], origin_y_root),
                (["domain_N02", "origin_x"], origin_x_nest_aligned),
                (["domain_N02", "origin_y"], origin_y_nest_aligned),
            ]
        )
    elif run == "no_origin":
        to_set.extend(
            [
                (["settings", "epsg"], epsg),
                (["domain_root", "origin_x"], origin_x_root),
                (["domain_root", "origin_lat"], origin_lat_root),
                (["domain_N02", "origin_y"], origin_y_nest_aligned),
                (["domain_N02", "origin_lon"], origin_lon_nest_aligned),
            ]
        )
    else:
        raise ValueError("Unknown parameter")

    modify_configuration(config_in, config_out, to_delete, to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    return config_out, tmp_path / file_out, file_ref


@pytest.fixture
def configuration_rotation(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Tuple[Path, Path, Path]:
    """Generate a configuration file with rotation and an adjusted result file.

    The result file is a combination from the default one and a diff file.

    Args:
        request: param attribute with the run and rotation_angle parameter.
        tmp_path: Temporary path for the configuration file and output files.

    Raises:
        ValueError: Unknown run parameter.

    Returns:
        Configuration file, output file and reference file.
    """
    run = request.param[0]
    rotation_angle = request.param[1]

    config_in = test_folder / "berlin_tiergarten.yml"
    config_out = tmp_path / f"berlin_tiergarten_no_coordinates_{rotation_angle}_{run}.yml"
    file_out = f"berlin_tiergarten_no_coordinates_{rotation_angle}_{run}"
    # reference file with be generated from the `orig` non-rotated case and a diff
    file_ref_orig = test_folder_ref / "berlin_tiergarten"
    file_ref_diff = test_folder_ref / f"diff_berlin_tiergarten_no_coordinates_{rotation_angle}"
    ref_dir = tmp_path / "output"
    file_ref = ref_dir / f"berlin_tiergarten_no_coordinates_{rotation_angle}_{run}"

    ref_dir.mkdir()

    # need to move the nested domain to be included in the parent domain
    # with the following, the nested domain is at the same relative position within the
    # parent domain as in the non-rotated case
    # all results except coordinates should be equal to non-rotated case
    diff_origin_x_invalid_rot, diff_origin_y_invalid_rot = GeoConverter._rotate(
        origin_x_nest_invalid - origin_x_root, origin_y_nest_invalid - origin_y_root, rotation_angle
    )
    origin_x_nest_rot_invalid = origin_x_root + diff_origin_x_invalid_rot
    origin_y_nest_rot_invalid = origin_y_root + diff_origin_y_invalid_rot

    # repeat for corrected coordinates
    diff_origin_x_aligned_rot, diff_origin_y_aligned_rot = GeoConverter._rotate(
        origin_x_nest_aligned - origin_x_root, origin_y_nest_aligned - origin_y_root, rotation_angle
    )
    origin_x_nest_rot_aligned = origin_x_root + diff_origin_x_aligned_rot
    origin_y_nest_rot_aligned = origin_y_root + diff_origin_y_aligned_rot

    lon, lat = GeoConverter._transform_points(
        rio.CRS.from_epsg(epsg),
        GeoConverter.crs_wgs84,
        [origin_x_nest_rot_aligned],
        [origin_y_nest_rot_aligned],
    )
    origin_lon_nest_rot_aligned = lon[0]
    origin_lat_nest_rot_aligned = lat[0]

    # remove coordinate inputs from config
    to_delete = []
    for file_input in ["input_15m", "input_3m"]:
        for file in ["x_utm", "y_utm", "lon", "lat"]:
            to_delete.append((file_input, "files", file))

    # set new values, use invalid origin_?_nest here to test the correction
    to_set: TO_SET = [(["output", "file_out"], file_out)]
    to_set.extend(
        [
            (["settings", "rotation_angle"], rotation_angle),
            (["settings", "epsg"], epsg),
            (["domain_root", "origin_x"], origin_x_root),
            (["domain_root", "origin_y"], origin_y_root),
        ]
    )
    if run == "lower_left":
        to_set.extend(
            [
                (["domain_N02", "lower_left_x"], lower_left_x_nest),
                (["domain_N02", "lower_left_y"], lower_left_y_nest),
            ]
        )
    elif run == "origin_xy":
        to_set.extend(
            [
                (["domain_N02", "origin_x"], origin_x_nest_rot_invalid),
                (["domain_N02", "origin_y"], origin_y_nest_rot_invalid),
            ]
        )
    else:
        raise ValueError("Unknown parameter")

    modify_configuration(config_in, config_out, to_delete, to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    # generate result file from non-rotated case and diff file
    # these diff files are generated after the runs with the following command:
    # for i in berlin_tiergarten_no_coordinates_*; do
    #     ncks -O -L9 -v E_UTM,N_UTM,lat,lon $i output/diff_$i
    # done
    for nest in ["_root", "_N02"]:
        shutil.copyfile(add_to_stem(file_ref_orig, nest), add_to_stem(file_ref, nest))
        ds = Dataset(add_to_stem(file_ref, nest), "a")
        ds_diff = Dataset(add_to_stem(file_ref_diff, nest), "r")
        for variable in ["E_UTM", "N_UTM", "lat", "lon"]:
            ds[variable][:] = ds_diff[variable][:]
        ds.rotation_angle = ds_diff.rotation_angle
        if nest == "_N02":
            ds.origin_x = origin_x_nest_rot_aligned
            ds.origin_y = origin_y_nest_rot_aligned
            ds.origin_lon = origin_lon_nest_rot_aligned
            ds.origin_lat = origin_lat_nest_rot_aligned
        ds.close()
        ds_diff.close()

    return config_out, tmp_path / file_out, file_ref


@pytest.fixture
def configuration_geotiff_input(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Tuple[Path, Path, Path]:
    """Generate a configuration file for geotiff input and geographic parameters.

    Args:
        request: Includes param attribute with the ignore_georeferencing parameter.
        tmp_path: Temporary directory for the test.

    Returns:
        Configuration file, the output file and the reference file.
    """
    ignore_georeferencing = request.param

    config_in = test_folder / "berlin_tiergarten.yml"
    config_out = tmp_path / f"berlin_tiergarten_geotiff_input_{ignore_georeferencing}.yml"
    file_out = f"berlin_tiergarten_geotiff_input_{ignore_georeferencing}"
    file_ref = test_folder_ref / "berlin_tiergarten"

    folder_geotiff = tmp_path / "input"

    # remove coordinate inputs from config
    to_delete: List[Tuple] = []
    for file_input in ["input_15m", "input_3m"]:
        for file in ["x_utm", "y_utm", "lon", "lat"]:
            to_delete.append((file_input, "files", file))
    if not ignore_georeferencing:
        for domain in ["domain_root", "domain_N02"]:
            for lower_left in ["input_lower_left_x", "input_lower_left_y"]:
                to_delete.append((domain, lower_left))

    # set new values depending
    to_set: TO_SET = [(["output", "file_out"], file_out)]
    to_set.extend(
        [
            (["settings", "epsg"], epsg),
            (["domain_root", "origin_x"], origin_x_root),
            (["domain_root", "origin_y"], origin_y_root),
            (["domain_N02", "origin_x"], origin_x_nest_aligned),
            (["domain_N02", "origin_y"], origin_y_nest_aligned),
        ]
    )
    to_set.extend(
        [
            (["settings", "ignore_input_georeferencing"], ignore_georeferencing),
        ]
    )

    to_replace = []
    for section in ["input_15m", "input_3m"]:
        to_set.append(([section, "path"], str(folder_geotiff)))
        to_replace.append(([section], ".nc", ".tif"))

    modify_configuration(config_in, config_out, to_delete, to_set, to_replace)
    modify_configuration_output(config_out, config_out, tmp_path)

    # Create the geotiff files.
    os.makedirs(folder_geotiff, exist_ok=True)

    crs = rio.CRS.from_epsg(epsg)
    for resolution in [3, 15]:
        top_left_x = top_left_x_input_data[resolution]
        top_left_y = top_left_y_input_data[resolution]

        for field in fields_input_data:
            file_nc = test_folder / f"input/Berlin_{field}_{resolution}m_DLR.nc"
            data_nc = Dataset(file_nc, "r")
            variable = data_nc.variables["Band1"]

            transform = riotf.from_origin(top_left_x, top_left_y, resolution, resolution)

            file_geotiff = folder_geotiff / f"Berlin_{field}_{resolution}m_DLR.tif"
            with rio.open(
                file_geotiff,
                "w",
                driver="GTiff",
                height=variable.shape[0],
                width=variable.shape[1],
                count=1,
                dtype=variable.dtype,
                crs=crs,
                transform=transform,
                nodata=variable._FillValue,
            ) as output_file:
                output_file.write(np.flipud(variable[:, :]), 1)
            data_nc.close()

    return config_out, tmp_path / file_out, file_ref


@pytest.fixture
def configuration_no_high_vegetation_elements(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Tuple[Path, Path, Path, Path]:
    """Generate a configuration file without certain high vegetation elements.

    No new references files are created because the diff might different dimension sizes. Instead
    yield path to the diff file.

    Args:
        request: param attribute with the run parameter.
        tmp_path: Temporary path for the configuration file and output files.

    Raises:
        ValueError: Unknown run parameter.

    Returns:
        Configuration file, output file, reference file and diff file.
    """
    run = request.param

    config_in = test_folder / "berlin_tiergarten.yml"
    config_out = tmp_path / f"berlin_tiergarten_{run}.yml"
    file_out = f"berlin_tiergarten_{run}"
    # Reference file and diff
    file_ref_orig = test_folder_ref / "berlin_tiergarten"
    file_ref_diff = test_folder_ref / f"diff_berlin_tiergarten_{run}"

    to_delete = []
    to_set: TO_SET = [(["output", "file_out"], file_out)]
    if run == "no_trees":
        for file_input in ["input_15m", "input_3m"]:
            for file in [
                "tree_height",
                "tree_crown_diameter",
                "tree_trunk_diameter",
                "tree_type",
            ]:
                to_delete.append((file_input, "files", file))
        to_set.extend([(["domain_root", "generate_single_trees"], False)])
        to_set.extend([(["domain_N02", "generate_single_trees"], False)])
    elif run == "no_patches":
        to_set.extend([(["domain_root", "generate_vegetation_patches"], False)])
        to_set.extend([(["domain_root", "replace_high_vegetation_types"], False)])
        to_set.extend([(["domain_N02", "generate_vegetation_patches"], False)])
        to_set.extend([(["domain_N02", "replace_high_vegetation_types"], False)])
    else:
        raise ValueError("Unknown parameter")

    modify_configuration(config_in, config_out, to_delete, to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    return config_out, tmp_path / file_out, file_ref_orig, file_ref_diff


@pytest.mark.usefixtures("config_counters")
def test_complete_run(configuration: Tuple[Path, Path, Path]):
    """Run the Berlin test case and compare with correct output."""
    create_driver(configuration[0], pdf=True, verbose={"gis": True})

    output_root = add_to_stem(configuration[1], "_root")
    output_nest = add_to_stem(configuration[1], "_N02")

    output_root_ref = add_to_stem(configuration[2], "_root")
    output_nest_ref = add_to_stem(configuration[2], "_N02")

    assert ncdf_equal(
        output_root_ref,
        output_root,
    ), "Root driver does not comply with reference"
    assert os.path.exists(output_root.parent / "berlin_tiergarten_root.pdf")

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
    ), "Nest driver does not comply with reference"
    assert os.path.exists(output_nest.parent / "berlin_tiergarten_N02.pdf")


@pytest.mark.usefixtures("config_counters")
def test_complete_run_geo_referenced(configuration_geo_referenced: Tuple[Path, Path, Path]):
    """Run the Berlin test case and compare with correct output."""
    create_driver(configuration_geo_referenced[0], pdf=True, verbose={"gis": True})

    output_root = add_to_stem(configuration_geo_referenced[1], "_root")
    output_nest = add_to_stem(configuration_geo_referenced[1], "_N02")

    output_root_ref = add_to_stem(configuration_geo_referenced[2], "_root")
    output_nest_ref = add_to_stem(configuration_geo_referenced[2], "_N02")

    assert ncdf_equal(
        output_root_ref,
        output_root,
    ), "Root driver does not comply with reference"
    assert os.path.exists(output_root.parent / "berlin_tiergarten_geo_referenced_root.pdf")

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
    ), "Nest driver does not comply with reference"
    assert os.path.exists(output_nest.parent / "berlin_tiergarten_geo_referenced_N02.pdf")


@pytest.mark.usefixtures("config_counters")
def test_trees_input(configuration_trees_input: Tuple[Path, Path, Path, Path]):
    """Run the Berlin test case with trees input for the nest."""
    create_driver(configuration_trees_input[0], verbose={"gis": True})

    output_root, output_nest = add_root_n02(configuration_trees_input[1])
    output_root_ref, output_nest_ref = add_root_n02(configuration_trees_input[2])

    output_nest_ref_diff = configuration_trees_input[3]

    # exclude crs from comparison because it is generate with pyproj and not read from file
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
        fields_exclude=[
            "lad",
            "bad",
            "tree_id",
            "tree_type",
            "vegetation_pars",
            "vegetation_type",
            "zlad",
        ],
    ), "Nest driver (except vegetation fields) does not comply with reference"

    # check only vegetation fields
    assert ncdf_equal(
        output_nest_ref_diff,
        output_nest,
        check_metadata=False,
        fields_only=[
            "lad",
            "bad",
            "tree_id",
            "tree_type",
            "vegetation_pars",
            "vegetation_type",
            "zlad",
        ],
    ), "Nest driver vegetation fields does not comply with reference"

    # check crs manually
    with Dataset(output_root) as nc_data:
        crs_root = nc_data.variables["crs"].__dict__
    with Dataset(output_nest) as nc_data:
        crs_nest = nc_data.variables["crs"].__dict__

    crs_ref = {
        "long_name": "coordinate reference system",
        "crs_wkt": 'PROJCRS["ETRS89 / UTM zone 33N",BASEGEOGCRS["ETRS89",'
        'DATUM["European Terrestrial Reference System 1989",'
        'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
        'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4258]],'
        'CONVERSION["UTM zone 33N",METHOD["Transverse Mercator",ID["EPSG",9807]],'
        'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],'
        'ID["EPSG",8801]],PARAMETER["Longitude of natural origin",15,'
        'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
        'PARAMETER["Scale factor at natural origin",0.9996,SCALEUNIT["unity",1],'
        'ID["EPSG",8805]],PARAMETER["False easting",500000,LENGTHUNIT["metre",1],'
        'ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],'
        'ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],'
        'AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",25833]]',
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.314140356,
        "inverse_flattening": 298.257222101,
        "reference_ellipsoid_name": "GRS 1980",
        "longitude_of_prime_meridian": 0.0,
        "prime_meridian_name": "Greenwich",
        "geographic_crs_name": "ETRS89",
        "horizontal_datum_name": "European Terrestrial Reference System 1989",
        "projected_crs_name": "ETRS89 / UTM zone 33N",
        "grid_mapping_name": "transverse_mercator",
        "latitude_of_projection_origin": 0.0,
        "longitude_of_central_meridian": 15.0,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "epsg_code": "EPSG:25833",
        "units": "m",
    }

    assert crs_root == crs_ref, "Root crs does not comply with reference"
    assert crs_nest == crs_ref, "Nest crs does not comply with reference"


@pytest.mark.parametrize(
    "configuration_no_coordinates_input",
    ["full", "origin_xy", "origin_lonlat"],
    indirect=True,
)
@pytest.mark.usefixtures("config_counters")
def test_no_coordinates_successful(configuration_no_coordinates_input: Tuple[Path, Path, Path]):
    """Run the Berlin test case without coordinate input but information to calculate it."""
    create_driver(configuration_no_coordinates_input[0], verbose={"gis": True})

    output_root, output_nest = add_root_n02(configuration_no_coordinates_input[1])
    output_root_ref, output_nest_ref = add_root_n02(configuration_no_coordinates_input[2])

    # exclude crs from comparison because it is generate with pyproj and not read from file
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Nest driver does not comply with reference"

    # check crs manually
    with Dataset(output_root) as nc_data:
        crs_root = nc_data.variables["crs"].__dict__
    with Dataset(output_nest) as nc_data:
        crs_nest = nc_data.variables["crs"].__dict__

    crs_ref = {
        "long_name": "coordinate reference system",
        "crs_wkt": 'PROJCRS["ETRS89 / UTM zone 33N",BASEGEOGCRS["ETRS89",'
        'DATUM["European Terrestrial Reference System 1989",'
        'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
        'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4258]],'
        'CONVERSION["UTM zone 33N",METHOD["Transverse Mercator",ID["EPSG",9807]],'
        'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],'
        'ID["EPSG",8801]],PARAMETER["Longitude of natural origin",15,'
        'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
        'PARAMETER["Scale factor at natural origin",0.9996,SCALEUNIT["unity",1],'
        'ID["EPSG",8805]],PARAMETER["False easting",500000,LENGTHUNIT["metre",1],'
        'ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],'
        'ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],'
        'AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",25833]]',
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.314140356,
        "inverse_flattening": 298.257222101,
        "reference_ellipsoid_name": "GRS 1980",
        "longitude_of_prime_meridian": 0.0,
        "prime_meridian_name": "Greenwich",
        "geographic_crs_name": "ETRS89",
        "horizontal_datum_name": "European Terrestrial Reference System 1989",
        "projected_crs_name": "ETRS89 / UTM zone 33N",
        "grid_mapping_name": "transverse_mercator",
        "latitude_of_projection_origin": 0.0,
        "longitude_of_central_meridian": 15.0,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "epsg_code": "EPSG:25833",
        "units": "m",
    }

    assert crs_root == crs_ref, "Root crs does not comply with reference"
    assert crs_nest == crs_ref, "Nest crs does not comply with reference"


@pytest.mark.parametrize(
    "configuration_no_coordinates_input",
    ["no_epsg", "no_origin"],
    indirect=True,
)
@pytest.mark.usefixtures("config_counters")
def test_no_coordinates_failing(configuration_no_coordinates_input):
    """Run the Berlin test case without coordinate input and not enough information to calculate it.

    This should raise an error.
    """
    with pytest.raises(ValueError):
        create_driver(configuration_no_coordinates_input[0])


run = ["lower_left", "origin_xy"]
angles = [30, 165, 200, 320]
combinations = list(itertools.product(run, angles))
names = [x + " " + str(y) for x, y in combinations]


@pytest.mark.parametrize(
    "configuration_rotation",
    combinations,
    ids=names,
    indirect=True,
)
@pytest.mark.usefixtures("config_counters")
def test_rotation(configuration_rotation: Tuple[Path, Path, Path]):
    """Run the Berlin test case with rotation."""
    create_driver(configuration_rotation[0])

    output_root, output_nest = add_root_n02(configuration_rotation[1])
    output_root_ref, output_nest_ref = add_root_n02(configuration_rotation[2])

    # exclude crs from comparison because it is generate with pyproj and not read from file
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Nest driver does not comply with reference"

    # check crs manually
    with Dataset(output_root) as nc_data:
        crs_root = nc_data.variables["crs"].__dict__
    with Dataset(output_nest) as nc_data:
        crs_nest = nc_data.variables["crs"].__dict__

    crs_ref = {
        "long_name": "coordinate reference system",
        "crs_wkt": 'PROJCRS["ETRS89 / UTM zone 33N",BASEGEOGCRS["ETRS89",'
        'DATUM["European Terrestrial Reference System 1989",'
        'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
        'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4258]],'
        'CONVERSION["UTM zone 33N",METHOD["Transverse Mercator",ID["EPSG",9807]],'
        'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],'
        'ID["EPSG",8801]],PARAMETER["Longitude of natural origin",15,'
        'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
        'PARAMETER["Scale factor at natural origin",0.9996,SCALEUNIT["unity",1],'
        'ID["EPSG",8805]],PARAMETER["False easting",500000,LENGTHUNIT["metre",1],'
        'ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],'
        'ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],'
        'AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",25833]]',
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.314140356,
        "inverse_flattening": 298.257222101,
        "reference_ellipsoid_name": "GRS 1980",
        "longitude_of_prime_meridian": 0.0,
        "prime_meridian_name": "Greenwich",
        "geographic_crs_name": "ETRS89",
        "horizontal_datum_name": "European Terrestrial Reference System 1989",
        "projected_crs_name": "ETRS89 / UTM zone 33N",
        "grid_mapping_name": "transverse_mercator",
        "latitude_of_projection_origin": 0.0,
        "longitude_of_central_meridian": 15.0,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "epsg_code": "EPSG:25833",
        "units": "m",
    }

    assert crs_root == crs_ref, "Root crs does not comply with reference"
    assert crs_nest == crs_ref, "Nest crs does not comply with reference"


@pytest.mark.parametrize(
    "configuration_geotiff_input",
    [True, False],
    indirect=True,
)
@pytest.mark.usefixtures("config_counters")
def test_geotiff_input(configuration_geotiff_input: Tuple[Path, Path, Path]):
    """Run the Berlin test case with geotiff input."""
    create_driver(configuration_geotiff_input[0], verbose={"gis": True})

    output_root, output_nest = add_root_n02(configuration_geotiff_input[1])
    output_root_ref, output_nest_ref = add_root_n02(configuration_geotiff_input[2])

    # exclude crs from comparison because it is generate with pyproj and not read from file
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        metadata_exclude_regex_paths=["crs"],
    ), "Nest driver does not comply with reference"

    # check crs manually
    with Dataset(output_root) as nc_data:
        crs_root = nc_data.variables["crs"].__dict__
    with Dataset(output_nest) as nc_data:
        crs_nest = nc_data.variables["crs"].__dict__

    crs_ref = {
        "long_name": "coordinate reference system",
        "crs_wkt": 'PROJCRS["ETRS89 / UTM zone 33N",BASEGEOGCRS["ETRS89",'
        'DATUM["European Terrestrial Reference System 1989",'
        'ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],'
        'PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4258]],'
        'CONVERSION["UTM zone 33N",METHOD["Transverse Mercator",ID["EPSG",9807]],'
        'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],'
        'ID["EPSG",8801]],PARAMETER["Longitude of natural origin",15,'
        'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],'
        'PARAMETER["Scale factor at natural origin",0.9996,SCALEUNIT["unity",1],'
        'ID["EPSG",8805]],PARAMETER["False easting",500000,LENGTHUNIT["metre",1],'
        'ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],'
        'ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],'
        'AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],ID["EPSG",25833]]',
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.314140356,
        "inverse_flattening": 298.257222101,
        "reference_ellipsoid_name": "GRS 1980",
        "longitude_of_prime_meridian": 0.0,
        "prime_meridian_name": "Greenwich",
        "geographic_crs_name": "ETRS89",
        "horizontal_datum_name": "European Terrestrial Reference System 1989",
        "projected_crs_name": "ETRS89 / UTM zone 33N",
        "grid_mapping_name": "transverse_mercator",
        "latitude_of_projection_origin": 0.0,
        "longitude_of_central_meridian": 15.0,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "epsg_code": "EPSG:25833",
        "units": "m",
    }

    assert crs_root == crs_ref, "Root crs does not comply with reference"
    assert crs_nest == crs_ref, "Nest crs does not comply with reference"


@pytest.mark.parametrize(
    "configuration_no_high_vegetation_elements",
    ["no_trees", "no_patches"],
    indirect=True,
)
@pytest.mark.usefixtures("config_counters")
def test_no_high_vegetation_elements(
    configuration_no_high_vegetation_elements: Tuple[Path, Path, Path, Path],
):
    """Run the Berlin test case without different high vegetation elements."""
    create_driver(configuration_no_high_vegetation_elements[0])

    output_root, output_nest = add_root_n02(configuration_no_high_vegetation_elements[1])

    # Treat the vegetation fields separately because replacing these fields only
    # in the original reference file is difficult due to possible difference in zlad

    # reference files (except vegetation fields)
    output_root_ref, output_nest_ref = add_root_n02(configuration_no_high_vegetation_elements[2])

    # reference files vegetation fields only
    output_root_ref_diff, output_nest_ref_diff = add_root_n02(
        configuration_no_high_vegetation_elements[3]
    )

    # exclude vegetation fields from comparison
    assert ncdf_equal(
        output_root_ref,
        output_root,
        metadata_significant_digits=4,
        fields_exclude=[
            "lad",
            "bad",
            "tree_id",
            "tree_type",
            "vegetation_pars",
            "vegetation_type",
        ],
    ), "Root driver (except vegetation fields) does not comply with reference"

    # check only vegetation fields
    assert ncdf_equal(
        output_root_ref_diff,
        output_root,
        check_metadata=False,
        fields_only=[
            "lad",
            "bad",
            "tree_id",
            "tree_type",
            "vegetation_pars",
            "vegetation_type",
        ],
    ), "Root driver vegetation fields does not comply with reference"

    # exclude vegetation fields from comparison
    assert ncdf_equal(
        output_nest_ref,
        output_nest,
        metadata_significant_digits=4,
        fields_exclude=[
            "lad",
            "bad",
            "tree_id",
            "tree_type",
            "vegetation_pars",
            "vegetation_type",
            "zlad",
        ],
    ), "Nest driver (except vegetation fields) does not comply with reference"

    # check only vegetation fields
    assert ncdf_equal(
        output_nest_ref_diff,
        output_nest,
        check_metadata=False,
        fields_only=[
            "lad",
            "bad",
            "tree_id",
            "tree_type",
            "vegetation_pars",
            "vegetation_type",
            "zlad",
        ],
    ), "Nest driver vegetation fields does not comply with reference"
