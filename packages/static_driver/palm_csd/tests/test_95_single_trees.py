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

"""Run the vegetation patch application with combinations of input data."""

from pathlib import Path
from typing import Tuple

import pytest

from palm_csd.create_driver import create_driver
from tests.tools import modify_configuration, ncdf_equal

test_folder = Path("tests/95_single_trees/")
test_folder_ref = Path("tests/95_single_trees/output/")


@pytest.fixture
def configuration(request: pytest.FixtureRequest, tmp_path: Path) -> Tuple[Path, Path, Path]:
    """Generate a configuration file for the single trees test."""
    lai = request.param

    config_in = test_folder / "single_trees.yml"
    config_out = tmp_path / "single_trees.yml"

    file_out = f"single_trees_{lai}"
    file_out_root = file_out + "_root"

    to_set = [
        (["output", "path"], str(tmp_path)),
        (["output", "file_out"], str(file_out)),
        (["settings", "use_lai_for_trees"], lai),
    ]

    modify_configuration(config_in, config_out, to_set=to_set)

    return config_out, tmp_path / file_out_root, test_folder_ref / file_out_root


@pytest.mark.usefixtures("config_counters")
@pytest.mark.parametrize("configuration", [True, False], indirect=True)
def test_single_trees(configuration: Tuple[Path, Path, Path]):
    """Run the Berlin test case and compare with correct output."""
    create_driver(configuration[0])

    assert ncdf_equal(
        configuration[2],
        configuration[1],
    )
