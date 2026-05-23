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

"""Check of overlap and clearance of domains."""

import logging
from pathlib import Path
from typing import Literal, Tuple

import pytest
from pytest import LogCaptureFixture

from palm_csd.create_driver import create_driver
from tests.tools import modify_configuration, modify_configuration_output

test_folder = Path("tests/96_domains/")


@pytest.fixture
def configuration_low_clearance(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Tuple[Literal["left", "right", "bottom", "top"], Path]:
    """Generate configuration files for low clearance tests.

    Args:
        request: Which run.
        tmp_path: Temporary path for storing the configuration file.

    Raises:
        ValueError: If the run is unknown.

    Returns:
        run, config_file
    """
    run = request.param

    config_in = test_folder / "low_clearance.yml"
    config_out = tmp_path / f"low_clearance_{run}.yml"

    if run == "left":
        to_set = [
            (["domain_N02", "lower_left_x"], 2),
            (["domain_N02", "lower_left_y"], 50),
        ]
    elif run == "right":
        to_set = [
            (["domain_N02", "lower_left_x"], 198 - 10),
            (["domain_N02", "lower_left_y"], 50),
        ]
    elif run == "bottom":
        to_set = [
            (["domain_N02", "lower_left_x"], 50),
            (["domain_N02", "lower_left_y"], 2),
        ]
    elif run == "top":
        to_set = [
            (["domain_N02", "lower_left_x"], 50),
            (["domain_N02", "lower_left_y"], 398 - 10),
        ]
    else:
        raise ValueError("Unknown parameter")

    modify_configuration(config_in, config_out, to_set=to_set)
    modify_configuration_output(config_out, config_out, tmp_path)

    return run, config_out


@pytest.mark.usefixtures("config_counters")
@pytest.mark.parametrize(
    "configuration_low_clearance",
    ["left", "right", "bottom", "top"],
    indirect=True,
)
def test_low_clearance(
    configuration_low_clearance: Tuple[Literal["left", "right", "bottom", "top"], Path],
):
    """Check if an error for low clearance is raised.

    Args:
        configuration_low_clearance: Which run to test and the configuration file.
    """
    run = configuration_low_clearance[0]
    config_file = configuration_low_clearance[1]
    with pytest.raises(ValueError, match=f"Not enough space at the {run}"):
        create_driver(config_file)


@pytest.mark.usefixtures("config_counters")
def test_overlap(caplog: LogCaptureFixture):
    """Check if an error for overlapping domains is raised."""
    with caplog.at_level(logging.WARNING):
        create_driver(test_folder / "overlap.yml")
    assert any("overlap" in record.message for record in caplog.records)
