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
# Print statistics of a PALM static driver and plot it.
#
# ------------------------------------------------------------------------------ #

"""Executable for a general static driver statistics tool."""

import argparse
import logging
import os

from palm_csd import ColorFormatter
from palm_csd.statistics import static_driver_statistics

if __name__ == "__main__":
    # unset the indention used in palm_csd log messages
    if isinstance(logging.root.handlers[0].formatter, ColorFormatter):
        logging.root.handlers[0].formatter.default_indent_non_status = 0

    parser = argparse.ArgumentParser(
        description="Print statistics and plot a static driver",
        epilog="See docs/palm_csd.md for more information.\n"
        + "Example usage:\n"
        + "  %(prog)s file.nc                          # print only statistics\n"
        + "  %(prog)s -s file.nc                       # statistics and shown plot\n"
        + '  %(prog)s --pdf -t "My Domain" config.yml  # show plot and store in PDF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # either --png or --pdf can be used to store a plot
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--png",
        action="store_true",
        help='store a plot as PNG by adding ".png" to ncfile name',
    )
    group.add_argument(
        "--pdf",
        action="store_true",
        help='store a plot as PDF by adding ".pdf" to ncfile name',
    )

    parser.add_argument("-s", "--show", action="store_true", help="show a plot")
    parser.add_argument(
        "-d",
        "--detailed",
        action="store_true",
        help="detailed surface types instead of basic categories",
    )
    parser.add_argument(
        "-g", "--geo-referenced", action="store_true", help="use georeferenced UTM coordinates"
    )
    parser.add_argument(
        "-H",
        "--height",
        action="store",
        help="plot width in inches",
    )
    parser.add_argument(
        "-W",
        "--width",
        action="store",
        help="plot width in inches",
    )
    parser.add_argument(
        "-t",
        "--title",
        action="store",
        help="plot title",
    )
    parser.add_argument("ncfile", help="static driver netCDF file")
    args = parser.parse_args()

    # generate plot file name if requested
    base, _ = os.path.splitext(args.ncfile)
    if args.png:
        plot_file = base + ".png"
    elif args.pdf:
        plot_file = base + ".pdf"
    else:
        plot_file = None

    width = float(args.width) if args.width is not None else None
    height = float(args.height) if args.height is not None else None

    static_driver_statistics(
        args.ncfile,
        show_plot=args.show,
        plot_file=plot_file,
        plot_title=args.title,
        plot_width=width,
        plot_height=height,
        detailed_plot=args.detailed,
        geo_referenced=args.geo_referenced,
    )
