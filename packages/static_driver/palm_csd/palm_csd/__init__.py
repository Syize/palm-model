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

"""Package of palm_csd, a tool to create static-driver files for the PALM model.

palm_csd is a tool to create static-driver files for the PALM model. It reads input data and
calculates the necessary fields for PALM from it.

This module initializes the logging system for palm_csd with the custom Logger StatusLogger and the
custom Formatter ColorFormatter.
"""

import logging

import pandas as pd

from palm_csd.logger import STATUS, ColorFormatter, StatusLogger

debug_logging = False
"""Debug logging mode with additional output not required by the user."""

# Use new approach of pandas.
pd.options.mode.copy_on_write = True

# Ensure that getLogger returns a StatusLogger.
logging.setLoggerClass(StatusLogger)

# Set-up of logging. By modifying the root logger, all loggers will use the same settings unless set
# otherwise.
logger_root = logging.getLogger()

# Set the logger to output to console by adding a StreamHandler.
handler = logging.StreamHandler()
logger_root.addHandler(handler)

# Set the logger to have color and indention by using ColorFormatter.
formatter = ColorFormatter(debug_logging)
handler.setFormatter(formatter)

# Set log-level default: show STATUS, INFO, WARNING, ERROR and CRITICAL. DEBUG is not shown by
# default.
logger_root.setLevel(STATUS)
