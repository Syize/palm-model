#!/usr/bin/env python3
# ------------------------------------------------------------------------------ #
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
# ------------------------------------------------------------------------------ #
#
# Description:
# ------------
# Create static-driver files for the PALM model system.
#
# ------------------------------------------------------------------------------ #

"""Main executable for palm_csd, a tool to create static-driver files for the PALM model."""

import argparse

from palm_csd.create_driver import create_driver

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a static driver from data input",
        epilog="See docs/palm_csd.md for more information.\n"
        + "Example usage:\n"
        + "  %(prog)s config.yml               # standard output, no plots\n"
        + "  %(prog)s -v all config.yml        # verbose all parts\n"
        + "  %(prog)s config.yml -v            # verbose all parts\n"
        + "  %(prog)s -v -s config.yml         # verbose all parts, show plot\n"
        + "  %(prog)s -v gis -v io config.yml  # verbose GIS and IO part\n"
        + "  %(prog)s -s --pdf config.yml      # show plot and store in PDF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # either --png or --pdf can be used to store a plot
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--png",
        action="store_true",
        help='store a plot as PNG by adding ".png" to output name',
    )
    group.add_argument(
        "--pdf",
        action="store_true",
        help='store a plot as PDF by adding ".pdf" to output name',
    )

    parser.add_argument(
        "-s", "--show", action="store_true", help="show a plot after static driver generation"
    )

    verbose_choices = ["all", "gis", "io", "misc", "vegetation"]
    parser.add_argument(
        "-v",
        "--verbose",
        action="append",
        const="all",
        nargs="?",
        choices=verbose_choices,
        help="produce additional output in selected processing parts",
    )
    parser.add_argument("ymlconfig", help="configuration file in yaml format")

    args = parser.parse_args()

    # Create a dictionary with verbose settings.
    #
    # If -v is not set, args.verbose is None and all verbose settings are set to False. If -v is set
    # without arguments, args.verbose is ["all"] and all verbose settings are set to True. If -v is
    # set with arguments, args.verbose is a list of strings and only the verbose settings that are
    # in the list are set to True.
    verbose = {}
    for v in verbose_choices:
        verbose[v] = args.verbose is not None and ("all" in args.verbose or v in args.verbose)

    create_driver(args.ymlconfig, verbose=verbose, show_plot=args.show, pdf=args.pdf, png=args.png)
