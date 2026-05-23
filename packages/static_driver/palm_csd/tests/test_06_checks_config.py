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

"""Test the CSDConfig classes."""

import copy
import os
import tempfile
from pathlib import Path
from typing import Dict, Generator, List, Tuple, Type, Union

import pytest
import rasterio.warp as riowp
import yaml
from pydantic import ValidationError

from palm_csd import csd_config
from palm_csd.constants import (
    NBUILDING_SURFACE_LAYER,
    IndexBuildingSurfaceLevel,
    IndexBuildingSurfaceType,
    InputData,
)
from palm_csd.csd_config import (
    CSDConfig,
    CSDConfigAttributes,
    CSDConfigDomain,
    CSDConfigElement,
    CSDConfigInput,
    CSDConfigOutput,
    CSDConfigSettings,
    _expand_parslike,
    _expand_scaling,
    _upscaling_method_default,
    _validate_scaling,
    reset_all_config_counters,
)
from palm_csd.csd_domain import CSDDomain
from tests.tools import modify_configuration

test_folder = Path("tests/99_full_application/")
test_configuration = test_folder / "berlin_tiergarten.yml"


@pytest.fixture
def valid_configuration() -> dict:
    """Return the valid tiergarten configuration from test_99."""
    with open(test_configuration, "r", encoding="utf-8") as file:
        conf = yaml.safe_load(file)
    # add input and domain section
    conf["input"] = conf["input_15m"]
    conf["domain"] = conf["domain_root"]
    return conf


TO_DELETE = List[List[str]]
TO_SET = List[Tuple[List[str], Union[str, float]]]


@pytest.fixture
def domain_wrong_range_tree_trunk_diameter(
    request: pytest.FixtureRequest,
) -> Generator[Tuple[CSDDomain, bool], None, None]:
    """Generate a CSDDomain with the tree trunk diameter wrong range input file."""
    config_in = test_configuration

    to_set: TO_SET = [
        (
            ["input_15m", "files", "tree_trunk_diameter"],
            "Berlin_trees_trunk_clean_15m_DLR_wrong_range.nc",
        ),
        (
            ["input_3m", "files", "tree_trunk_diameter"],
            "Berlin_trees_trunk_clean_3m_DLR_wrong_range.nc",
        ),
        (["settings", "replace_invalid_input_values"], request.param),
    ]

    config_dic = modify_configuration(config_in, to_set=to_set)
    config = CSDConfig(config_dic)

    yield CSDDomain("root", config), request.param

    reset_all_config_counters()


@pytest.mark.parametrize(
    "csdconfig_class",
    [CSDConfigOutput, CSDConfigDomain],
)
def test_missing_settings(csdconfig_class: Type[CSDConfigElement]):
    """Test that a configuration without settings raises an error."""
    with pytest.raises(ValidationError, match="type=missing"):
        csdconfig_class()


@pytest.mark.parametrize(
    "csdconfig_class",
    [
        CSDConfigAttributes,
        CSDConfigSettings,
        CSDConfigOutput,
        CSDConfigInput,
        CSDConfigDomain,
    ],
)
def test_invalid_settings(valid_configuration: dict, csdconfig_class: Type[CSDConfigElement]):
    """Test that an invalid configuration raises an error."""
    configuration_valid = valid_configuration[csdconfig_class._type]
    instance_valid = csdconfig_class(**configuration_valid)
    csdconfig_class._reset_counter()

    for variable in instance_valid.__dict__.keys():
        # bridge_depth in CSDSettings is not used anymore
        # TODO: remove once warning is removed
        if variable == "bridge_depth" and csdconfig_class == CSDConfigSettings:
            continue
        # pixel_size in CSDConfigInput is not used anymore
        # TODO: remove once warning is removed
        if variable == "pixel_size" and csdconfig_class == CSDConfigInput:
            continue
        # lai_roof_extensive in CSDConfigInput is not used anymore
        # TODO: remove once warning is removed
        if variable == "lai_roof_extensive" and csdconfig_class == CSDConfigSettings:
            continue
        # lai_roof_intensive in CSDConfigInput is not used anymore
        # TODO: remove once warning is removed
        if variable == "lai_roof_intensive" and csdconfig_class == CSDConfigSettings:
            continue
        if variable in csd_config.value_defaults:
            minimum = csd_config.value_defaults[variable].minimum
            maximum = csd_config.value_defaults[variable].maximum
            if minimum is not None:
                configuration_invalid = copy.deepcopy(configuration_valid)
                # for list variables, the ranges for the list elements are checked
                if isinstance(getattr(instance_valid, variable), list):
                    configuration_invalid[variable] = copy.deepcopy(
                        getattr(instance_valid, variable)
                    )
                    configuration_invalid[variable][0] = minimum - 1
                else:
                    configuration_invalid[variable] = minimum - 1
                with pytest.raises(
                    ValidationError, match="(type=greater_than_equal)|(type=value_error)"
                ):
                    csdconfig_class(**configuration_invalid)
                csdconfig_class._reset_counter()
            if maximum is not None:
                configuration_invalid = copy.deepcopy(configuration_valid)
                # for list variables, the ranges for the list elements are checked
                if isinstance(getattr(instance_valid, variable), list):
                    configuration_invalid[variable] = copy.deepcopy(
                        getattr(instance_valid, variable)
                    )
                    configuration_invalid[variable][-1] = maximum + 1
                else:
                    configuration_invalid[variable] = maximum + 1
                with pytest.raises(ValidationError, match="type=less_than_equal"):
                    csdconfig_class(**configuration_invalid)
                csdconfig_class._reset_counter()


@pytest.mark.parametrize("domain_wrong_range_tree_trunk_diameter", [True, False], indirect=True)
def test_invalid_input_values(domain_wrong_range_tree_trunk_diameter: Tuple[CSDDomain, bool]):
    """Test that an invalid input value raises an error."""
    domain, replace_invalid = domain_wrong_range_tree_trunk_diameter

    if replace_invalid:
        # no default tree trunk value so invalid values masked
        tree_trunk = domain.read(InputData.tree_trunk_diameter)
        assert tree_trunk.mask.all() >= 0
    else:
        # error when invalid values are found and not to be replaced
        with pytest.raises(ValueError):
            domain.read(InputData.tree_trunk_diameter)


def test_paths() -> None:
    """Test the path handling of the configuration classes."""
    # We need to suppress the warning that string literals are not allowed as Path
    # arguments in the pydantic models. str is necessary for the tests.

    # CSDConfigOutput successful

    file_out = "test_output"
    # create existing temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        # path and file
        CSDConfigOutput(path=tmpdirname, file_out=file_out)  # type: ignore
        CSDConfigOutput._reset_counter()
        # only file
        CSDConfigOutput(file_out=tmpdirname + "/" + file_out)  # type: ignore
        CSDConfigOutput._reset_counter()
        # path with ~ and file
        tmpdirname_relative_home = "~/" + os.path.relpath(tmpdirname, os.path.expanduser("~"))
        CSDConfigOutput(path=tmpdirname_relative_home, file_out=file_out)  # type: ignore
        CSDConfigOutput._reset_counter()
        # only file with ~
        CSDConfigOutput(file_out=tmpdirname_relative_home + "/" + file_out)  # type: ignore
        CSDConfigOutput._reset_counter()

    # CSDConfigOutput failure

    # folder of tmpdirname does not exist now
    with pytest.raises(ValidationError, match="type=path_not_directory"):
        # path does not exist
        CSDConfigOutput(path=tmpdirname, file_out=file_out)  # type: ignore
    with pytest.raises(ValidationError, match="type=value_error"):
        # parent of file_out does not exist
        CSDConfigOutput(file_out=tmpdirname + "/" + file_out)  # type: ignore

    # CSDConfigInput successful

    file_zt = "Berlin_terrain_height_15m_DLR.nc"
    path = "tests/99_full_application/input"
    # path and file
    CSDConfigInput(path=path, files={"zt": file_zt})  # type: ignore
    CSDConfigInput._reset_counter()
    # only file
    CSDConfigInput(files={"zt": path + "/" + file_zt})  # type: ignore
    CSDConfigInput._reset_counter()
    # path with ~ and file
    path_relative_home = "~/" + os.path.relpath(path, os.path.expanduser("~"))
    CSDConfigInput(path=path_relative_home, files={"zt": file_zt})  # type: ignore
    CSDConfigInput._reset_counter()
    # only file with ~
    CSDConfigInput(files={"zt": path_relative_home + "/" + file_zt})  # type: ignore
    CSDConfigInput._reset_counter()

    # CSDConfigInput failure

    with pytest.raises(ValidationError, match="type=path_not_directory"):
        # path does not exist
        CSDConfigInput(path=tmpdirname, files={"zt": file_out})  # type: ignore
    with pytest.raises(ValidationError, match="type=path_not_file"):
        # file does not exist
        CSDConfigInput(files={"zt": tmpdirname + "/" + file_zt})  # type: ignore


def test_expand_validate_scaling() -> None:
    """Test the _expand_scaling and _validate_scaling functions."""
    # single value for all values
    scaling_str = "nearest"
    scaling_resampling = riowp.Resampling.nearest
    scaling_expanded = {
        "categorical": riowp.Resampling.nearest,
        "continuous": riowp.Resampling.nearest,
        "discontinuous": riowp.Resampling.nearest,
        "discrete": riowp.Resampling.nearest,
    }
    result = _expand_scaling(scaling_str, _upscaling_method_default)
    assert result == scaling_expanded
    result = _expand_scaling(scaling_resampling, _upscaling_method_default)
    assert result == scaling_expanded
    assert _validate_scaling(result) == result

    scaling_str = "average"
    scaling_resampling = riowp.Resampling.average
    scaling_expanded = {
        "categorical": riowp.Resampling.average,
        "continuous": riowp.Resampling.average,
        "discontinuous": riowp.Resampling.average,
        "discrete": riowp.Resampling.average,
    }
    result = _expand_scaling(scaling_str, _upscaling_method_default)
    assert result == scaling_expanded
    result = _expand_scaling(scaling_resampling, _upscaling_method_default)
    assert result == scaling_expanded
    # categorical should not be average
    with pytest.raises(ValueError):
        _validate_scaling(result)

    scaling_dict = {"categorical": "nearest", "discrete": "bilinear"}
    scaling_resampling = {
        "categorical": riowp.Resampling.nearest,
        "discrete": riowp.Resampling.bilinear,
    }
    scaling_expanded = {
        "categorical": riowp.Resampling.nearest,
        "continuous": riowp.Resampling.average,
        "discontinuous": riowp.Resampling.average,
        "discrete": riowp.Resampling.bilinear,
    }
    result = _expand_scaling(scaling_dict, _upscaling_method_default)
    assert result == scaling_expanded
    result = _expand_scaling(scaling_resampling, _upscaling_method_default)
    assert result == scaling_expanded
    assert _validate_scaling(result) == result


def test_expand_parslike() -> None:
    """Test the `_expand_parslike` function."""
    # single value for all values
    building_heat_conductivity_float: Union[float, List[float]] = 1.8
    building_heat_conductivity_expanded = {
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
    assert (
        _expand_parslike(
            building_heat_conductivity_float,
            IndexBuildingSurfaceType,
            nlayer=NBUILDING_SURFACE_LAYER,
        )
        == building_heat_conductivity_expanded
    )

    # single list for all values
    building_heat_conductivity_float = [1.8, 1.8, 1.8, 1.9]
    building_heat_conductivity_expanded = {
        0: [1.8, 1.8, 1.8, 1.9],
        1: [1.8, 1.8, 1.8, 1.9],
        2: [1.8, 1.8, 1.8, 1.9],
        3: [1.8, 1.8, 1.8, 1.9],
        4: [1.8, 1.8, 1.8, 1.9],
        5: [1.8, 1.8, 1.8, 1.9],
        6: [1.8, 1.8, 1.8, 1.9],
        7: [1.8, 1.8, 1.8, 1.9],
        8: [1.8, 1.8, 1.8, 1.9],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity_float,
            IndexBuildingSurfaceType,
            nlayer=NBUILDING_SURFACE_LAYER,
        )
        == building_heat_conductivity_expanded
    )

    # single value for all wall
    building_heat_conductivity: Dict[str, Union[float, List[float]]] = {"wall": 1.8}
    building_heat_conductivity_expanded = {
        0: [1.8, 1.8, 1.8, 1.8],
        1: [1.8, 1.8, 1.8, 1.8],
        2: [1.8, 1.8, 1.8, 1.8],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )

    # single value for wall_agfl
    building_heat_conductivity = {"wall_agfl": 1.8}
    building_heat_conductivity2: Dict[int, Union[float, List[float]]] = {1: 1.8}
    building_heat_conductivity_expanded = {
        1: [1.8, 1.8, 1.8, 1.8],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )
    assert (
        _expand_parslike(
            building_heat_conductivity2, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )

    # combination
    building_heat_conductivity = {"wall_agfl": 1.8, "green": [1.5, 1.6, 1.7, 1.8]}
    building_heat_conductivity2 = {
        1: 1.8,
        6: [1.5, 1.6, 1.7, 1.8],
        7: [1.5, 1.6, 1.7, 1.8],
        8: [1.5, 1.6, 1.7, 1.8],
    }
    building_heat_conductivity_expanded = {
        1: [1.8, 1.8, 1.8, 1.8],
        6: [1.5, 1.6, 1.7, 1.8],
        7: [1.5, 1.6, 1.7, 1.8],
        8: [1.5, 1.6, 1.7, 1.8],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )
    assert (
        _expand_parslike(
            building_heat_conductivity2, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )

    # mixed single value for all wall
    building_heat_conductivity3: Dict[Union[int, str], Union[float, List[float]]] = {
        "wall": 1.8,
        6: 1.9,
    }
    building_heat_conductivity_expanded = {
        0: [1.8, 1.8, 1.8, 1.8],
        1: [1.8, 1.8, 1.8, 1.8],
        2: [1.8, 1.8, 1.8, 1.8],
        6: [1.9, 1.9, 1.9, 1.9],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity3, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )

    # mixed single value for wall_agfl
    building_heat_conductivity3 = {"wall_agfl": 1.8, 6: 1.9}
    building_heat_conductivity_expanded = {
        1: [1.8, 1.8, 1.8, 1.8],
        6: [1.9, 1.9, 1.9, 1.9],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity3, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )

    # mixed single list for wall_agfl
    building_heat_conductivity3 = {"wall_agfl": [1.8, 1.8, 1.9, 1.10], 6: 1.9}
    building_heat_conductivity_expanded = {
        1: [1.8, 1.8, 1.9, 1.10],
        6: [1.9, 1.9, 1.9, 1.9],
    }
    assert (
        _expand_parslike(
            building_heat_conductivity3, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
        == building_heat_conductivity_expanded
    )

    # too long list for wall_agfl
    building_heat_conductivity = {"wall_agfl": [1.8, 1.8, 1.9, 1.10, 1.11]}
    building_heat_conductivity2 = {1: [1.8, 1.8, 1.9, 1.10, 1.11]}
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity2, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )

    # too short list for wall_agfl
    building_heat_conductivity = {
        "wall_agfl": [
            1.8,
            1.8,
        ]
    }
    building_heat_conductivity2 = {
        1: [
            1.8,
            1.8,
        ]
    }
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity2, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )

    # too long list for wall
    building_heat_conductivity = {"wall": [1.8, 1.8, 1.9, 1.10, 1.11]}
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )

    building_heat_conductivity = {
        "wall": [
            1.8,
            1.8,
        ]
    }
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )

    # too long list for all values
    building_heat_conductivity_float = [1.8, 1.8, 1.9, 1.10, 1.11]
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity_float,
            IndexBuildingSurfaceType,
            nlayer=NBUILDING_SURFACE_LAYER,
        )

    # too short list for all values
    building_heat_conductivity_float = [
        1.8,
        1.8,
    ]
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )

    # wrong entry
    building_heat_conductivity = {"blub": 1.8}
    building_heat_conductivity2 = {10: 1.8}
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )
    with pytest.raises(ValueError):
        _expand_parslike(
            building_heat_conductivity2, IndexBuildingSurfaceType, nlayer=NBUILDING_SURFACE_LAYER
        )

    # no layers, single value
    building_lai_float: Union[float, List[float]] = 3.0
    building_lai_expanded = {0: 3.0, 1: 3.0, 2: 3.0}
    assert _expand_parslike(building_lai_float, IndexBuildingSurfaceLevel) == building_lai_expanded

    # no layers
    building_lai: Dict[str, Union[float, List[float]]] = {"agfl": 3.0}
    building_lai2: Dict[int, Union[float, List[float]]] = {1: 3.0}
    building_lai_expanded = {1: 3.0}
    assert _expand_parslike(building_lai, IndexBuildingSurfaceLevel) == building_lai_expanded
    assert _expand_parslike(building_lai2, IndexBuildingSurfaceLevel) == building_lai_expanded

    # no layers, too long
    building_lai_float = [3.0, 4.0]
    with pytest.raises(ValueError):
        _expand_parslike(building_lai_float, IndexBuildingSurfaceLevel)

    # no layers, too long
    building_lai = {"agfl": [3.0, 4.0]}
    building_lai2 = {1: [3.0, 4.0]}
    with pytest.raises(ValueError):
        _expand_parslike(building_lai, IndexBuildingSurfaceLevel)
    with pytest.raises(ValueError):
        _expand_parslike(building_lai2, IndexBuildingSurfaceLevel)


@pytest.fixture
def configuration_dict_one_input() -> Generator[Dict, None, None]:
    """Generate a configuration with only one input section.

    Returns:
        Configuration dictionary.
    """
    to_delete = [
        [
            "input_3m",
        ],
        ["domain_root", "input"],
        ["domain_N02", "input"],
    ]

    yield modify_configuration(config_in=test_configuration, to_delete=to_delete)

    reset_all_config_counters()


def test_one_input(configuration_dict_one_input: Dict):
    """Test the configuration with only one input section."""
    config = CSDConfig(configuration_dict_one_input)
    for domain_name in ["root", "N02"]:
        input = config.input_of_domain(domain_name)
        assert input == config.input_dict["15m"]


@pytest.fixture
def configuration_dict_named_input() -> Generator[Dict, None, None]:
    """Generate a configuration dictionary with named input in domain sections.

    Returns:
        Configuration dictionary.
    """
    with open(test_configuration, "r", encoding="utf-8") as file:
        yield yaml.safe_load(file)

    reset_all_config_counters()


def test_named_input(configuration_dict_named_input: Dict):
    """Test the configuration with named input section."""
    config = CSDConfig(configuration_dict_named_input)
    input_name = {"root": "15m", "N02": "3m"}
    for domain_name in ["root", "N02"]:
        input = config.input_of_domain(domain_name)
        assert input == config.input_dict[input_name[domain_name]]
