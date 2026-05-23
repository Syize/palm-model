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

"""Tools for logging.

StatusLogger has a new log level STATUS between INFO and DEBUG.
ColorFormatter adds colors to the log messages and supports indention to show different levels of
messages.
"""

import logging
from typing import Callable, Final, NoReturn, Optional, Type

import numpy as np
import numpy.typing as npt

# Add new log level STATUS.
STATUS: Final = logging.INFO - 5
"""STATUS log level between INFO and DEBUG."""
logging.addLevelName(STATUS, "STATUS")


class StatusLogger(
    logging.Logger
):  # assume logger class was not modified, mypy does not like logging.getLoggerClass()
    """Standard Logger with added level STATUS and optional hierarchy information."""

    def debug_indent(self, msg, hierarchy: int = 1, *args, **kwargs) -> None:
        """Log message with level DEBUG and indention of the defined hierarchy.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            *args: Variable length argument list for the debug method.
            **kwargs: Arbitrary keyword arguments for the debug method.
        """
        super().debug(msg, extra={"hierarchy": hierarchy}, stacklevel=2, *args, **kwargs)

    def status(self, msg, *args, **kwargs) -> None:
        """Log message with level STATUS.

        Args:
            msg: Message to log.
            *args: Variable length argument list for the _log method.
            **kwargs: Arbitrary keyword arguments for the _log method.
        """
        if self.isEnabledFor(STATUS):
            self._log(STATUS, msg, args, stacklevel=2, **kwargs)

    def status_indent(self, msg, hierarchy: int = 1, *args, **kwargs) -> None:
        """Log message with level STATUS and indention of the defined hierarchy.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            *args: Variable length argument list for the _log method.
            **kwargs: Arbitrary keyword arguments for the _log method.
        """
        if self.isEnabledFor(STATUS):
            self._log(STATUS, msg, args, extra={"hierarchy": hierarchy}, stacklevel=2, **kwargs)

    def info_indent(self, msg, hierarchy: int = 1, *args, **kwargs) -> None:
        """Log message with level INFO and indention of the defined hierarchy.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            *args: Variable length argument list for the info method.
            **kwargs: Arbitrary keyword arguments for the info method.
        """
        super().info(msg, extra={"hierarchy": hierarchy}, stacklevel=2, *args, **kwargs)

    def warning_indent(self, msg, hierarchy: int = 1, *args, **kwargs) -> None:
        """Log message with level WARNING and indention of the defined hierarchy.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            *args: Variable length argument list for the warning method.
            **kwargs: Arbitrary keyword arguments for the warning method.
        """
        super().warning(msg, extra={"hierarchy": hierarchy}, stacklevel=2, *args, **kwargs)

    def error_indent(self, msg, hierarchy: int = 1, *args, **kwargs) -> None:
        """Log message with level ERROR and indention of the defined hierarchy.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            *args: Variable length argument list for the error method.
            **kwargs: Arbitrary keyword arguments for the error method.
        """
        super().error(msg, extra={"hierarchy": hierarchy}, stacklevel=2, *args, **kwargs)

    def critical_indent(self, msg, hierarchy: int = 1, *args, **kwargs) -> None:
        """Log message with level CRITICAL and indention of the defined hierarchy.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            *args: Variable length argument list for the critical method.
            **kwargs: Arbitrary keyword arguments for the critical method.
        """
        super().critical(msg, extra={"hierarchy": hierarchy}, stacklevel=2, *args, **kwargs)

    def critical_raise(
        self, msg, exception_type: Type[Exception] = ValueError, *args, **kwargs
    ) -> NoReturn:
        """Log message with level CRITICAL. Raise exception.

        Args:
            msg: Message to log.
            exception_type: Exception type to raise. Defaults to ValueError.
            *args: Variable length argument list for the critical method.
            **kwargs: Arbitrary keyword arguments for the critical method.
        """
        super().critical(msg, stacklevel=2, *args, **kwargs)
        raise exception_type(msg)

    def critical_indent_raise(
        self, msg, hierarchy: int = 1, exception_type: Type[Exception] = ValueError, *args, **kwargs
    ) -> NoReturn:
        """Log message with level CRITICAL and indention of the defined hierarchy. Raise exception.

        Args:
            msg: Message to log.
            hierarchy: Indention hierarchy. Defaults to 1.
            exception_type: Exception type to raise. Defaults to ValueError.
            *args: Variable length argument list for the critical method.
            **kwargs: Arbitrary keyword arguments for the critical method.
        """
        super().critical(msg, extra={"hierarchy": hierarchy}, stacklevel=2, *args, **kwargs)
        raise exception_type(msg)

    def critical_argwhere(
        self,
        message: str,
        array: npt.NDArray[np.bool_],
        message2: Optional[str] = None,
        indent: int = 0,
    ) -> None:
        """Log with level CRITICAL/INFO details of the True elements in array if present.

        First, the number of True elements in the array is logged with level CRITICAL, surrounded by
        message and, optionally, message2. The coordinates of the True elements of the array are
        logged with level INFO.

        Args:
            message: First part of the message.
            array: Array to count and print the coordinates of the True elements.
            message2: Second part of the message. Defaults to None.
            indent: Indention hierarchy. Defaults to 0.
        """
        if not array.any():
            return
        self.critical_indent(f"{message} {array.sum()} {message2}", indent)
        self.info_indent("Coordinates:", indent + 1)
        self.info_indent(nonzero_element_coords_string(array), indent + 1)

    def critical_argwhere_raise(
        self,
        message: str,
        array: npt.NDArray[np.bool_],
        message2: Optional[str] = None,
        indent: int = 0,
        exception_type: Type[Exception] = ValueError,
    ) -> None:
        """Log and raise with level CRITICAL/INFO details of the True elements in array if present.

        First, the number of True elements in the array is logged with level CRITICAL, surrounded by
        message and, optionally, message2. The coordinates of the True elements of the array are
        logged with level INFO.

        Args:
            message: First part of the message.
            array: Array to count and print the coordinates of the True elements.
            message2: Second part of the message. Defaults to None.
            indent: Indention hierarchy. Defaults to 0.
            exception_type: Exception type to raise. Defaults to ValueError.
        """
        if not array.any():
            return
        self.critical_argwhere(message, array, message2, indent)
        raise exception_type(f"{message} {array.sum()} {message2}")

    def warning_argwhere(
        self,
        message: str,
        array: npt.NDArray[np.bool_],
        message2: Optional[str] = None,
        indent: int = 0,
    ) -> None:
        """Log with level WARNING/DEBUG details of the True elements in array if present.

        First, the number of True elements in the array is logged with level WARNING, surrounded by
        message and, optionally, message2. The coordinates of the True elements of the array are
        logged with level DEBUG.

        Args:
            message: First part of the message.
            array: Array to count and print the coordinates of the True elements.
            message2: Second part of the message. Defaults to None.
            indent: Indention hierarchy. Defaults to 0.
        """
        if not array.any():
            return
        self.warning_indent(f"{message} {array.sum()} {message2}", indent)
        self.debug_indent("Coordinates:", indent + 1)
        self.debug_indent(nonzero_element_coords_string(array), indent + 1)


def nonzero_element_coords_string(array: npt.ArrayLike) -> str:
    """Return a string with the coordinates of non-zero elements of the array.

    In particular, for boolean arrays, the coordinates of the True elements are produced.

    Args:
        array: Input array.

    Returns:
        String of the coordinates of non-zero elements of the array.
    """
    coords = list(map(tuple, np.argwhere(array)))
    if not coords:
        return ""

    # Calculate the maximum width for each column.
    max_widths = [max(len(str(coord[i])) for coord in coords) for i in range(len(coords[0]))]

    # Format each coordinate with the calculated widths.
    formatted_coords = [
        "(" + ", ".join(f"{coord[i]:>{max_widths[i]}}" for i in range(len(coord))) + ")"
        for coord in coords
    ]

    # Add a newline after every 5 tuples.
    lines = [", ".join(formatted_coords[i : i + 5]) for i in range(0, len(formatted_coords), 5)]
    return "\n".join(lines)


class ColorFormatter(logging.Formatter):
    """Logging Formatter to add colors, hierarchy spacing and, optionally, more detailed output.

    The _formats method is used to define the formatting for each log level. This includes the
    color, the hierarchy spacing according to the log level and, optionally, the debug information.
    """

    GREY: Final = "\x1b[38;20m"
    """Grey color terminal code."""
    YELLOW: Final = "\x1b[33;20m"
    """Yellow color terminal code."""
    RED: Final = "\x1b[31;20m"
    """RED color terminal code."""
    BOLD_RED: Final = "\x1b[31;1m"
    """Bold red color terminal code."""
    RESET: Final = "\x1b[0m"
    """Reset color terminal code."""

    format_debug = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    """Debug format string."""
    format_levelname = "%(levelname)s: %(message)s"
    """Format string with level name."""
    format_no_levelname = "%(message)s"
    """Format string without level name."""

    space_hierarchy: str = "  "
    """Indentation space for each hierarchy level."""
    default_indent_non_status: int = 1
    """Default indentation level for all log levels except STATUS."""

    _formats: Callable[[int, int], str]
    """Function to generate format string depending on loglevel and hierarchy."""

    def __init__(self, debug: bool = False) -> None:
        """Initialize the ColorFormatter. Set the format string function according to debug.

        Args:
            debug: If True, use debug formatting. Defaults to False.
        """
        super().__init__()
        if debug:
            self._formats = self.formats_debug
        else:
            self._formats = self.formats

    def _msg_mod(self, record: logging.LogRecord, hierarchy: int = 0) -> None:
        """Add indentation to each newline in the message for consistent hierarchy spacing.

        Args:
            record: Log record.
            hierarchy: Indention hierarchy. Defaults to 0.
        """
        if record.levelno == STATUS:
            space = self.space_hierarchy * (hierarchy)
        elif record.levelno == logging.WARNING:
            space = self.space_hierarchy * (hierarchy + self.default_indent_non_status) + " " * len(
                "WARNING: "
            )
        elif record.levelno == logging.ERROR:
            space = self.space_hierarchy * (hierarchy + self.default_indent_non_status) + " " * len(
                "ERROR: "
            )
        elif record.levelno == logging.CRITICAL:
            space = self.space_hierarchy * (hierarchy + self.default_indent_non_status) + " " * len(
                "CRITICAL: "
            )
        else:
            space = self.space_hierarchy * (hierarchy + self.default_indent_non_status)
        record.msg = record.getMessage().replace("\n", "\n" + space)

    def formats(self, levelno: int, hierarchy: int = 0) -> str:
        """Format string with color and hierarchy spacing.

        Args:
            levelno: Log level number.
            hierarchy: Indention hierarchy. Defaults to 0.

        Returns:
            Format string.
        """
        if levelno == logging.DEBUG:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_no_levelname
                + self.RESET
            )
        elif levelno == STATUS:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy)
                + self.format_no_levelname
                + self.RESET
            )
        elif levelno == logging.INFO:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_no_levelname
                + self.RESET
            )
        elif levelno == logging.WARNING:
            return (
                self.YELLOW
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_levelname
                + self.RESET
            )
        elif levelno == logging.ERROR:
            return (
                self.RED
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_levelname
                + self.RESET
            )
        elif levelno == logging.CRITICAL:
            return (
                self.BOLD_RED
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_levelname
                + self.RESET
            )
        else:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_no_levelname
                + self.RESET
            )

    def formats_debug(self, levelno: int, hierarchy: int = 0) -> str:
        """Format string with color, hierarchy spacing and debug information.

        Args:
            levelno: Log level number.
            hierarchy: Indention hierarchy. Defaults to 0.

        Returns:
            Format string.
        """
        if levelno == logging.DEBUG:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_debug
                + self.RESET
            )
        elif levelno == STATUS:
            return self.GREY + self.space_hierarchy * (hierarchy) + self.format_debug + self.RESET
        elif levelno == logging.INFO:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_debug
                + self.RESET
            )
        elif levelno == logging.WARNING:
            return (
                self.YELLOW
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_debug
                + self.RESET
            )
        elif levelno == logging.ERROR:
            return (
                self.RED
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_debug
                + self.RESET
            )
        elif levelno == logging.CRITICAL:
            return (
                self.BOLD_RED
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_debug
                + self.RESET
            )
        else:
            return (
                self.GREY
                + self.space_hierarchy * (hierarchy + self.default_indent_non_status)
                + self.format_debug
                + self.RESET
            )

    def format(self, record: logging.LogRecord) -> str:
        """Format the log message as text with color and hierarchy spacing.

        The log message is modified to add indentation for consistent hierarchy spacing and the
        formatting is adjusted to the log level.

        Args:
            record: Log record.

        Returns:
            Formatted log message.
        """
        self._msg_mod(record, getattr(record, "hierarchy", 0))
        self._style._fmt = self._formats(record.levelno, getattr(record, "hierarchy", 0))
        return super().format(record)
