---
title: Overview
---
# LES-Model

This is the user guide for the LES-Model.

---

!!! warning
    This site is Work in Progress.

    ToDo:

    - [ ] What is it? What is its purpose?
    - [ ] How to install?
    - [ ] How to run?
    - [ ] Example
    - [ ] Physical/technical background

## Currently available pages

- [Static Driver](static.md)<br>
  Note that all variable arrays defined in the static driver require the `_FillValue` attribute. This is a standard netCDF attribute. It is often need to write a value to represent undefined or missing values, e.g. for grid points where no buildings are defined. The fill value provides an appropriate value for this purpose. It should normally be outside the valid range of values of the respective variable and therefore treated as missing when read by generic applications (e.g. graphic software that reads the netCDF file) or PALM. Depending on the data type of the variable, recommended default values are <br>

  | Type     | suggested default value |
  | -------- | ----------------------- |
  | NC_FLOAT | -9999.f                 |
  | NC_INT   | -9999                   |
  | NC_BYTE  | -127b                   |
  
- [Dynamic Driver](dynamic.md)<br>
