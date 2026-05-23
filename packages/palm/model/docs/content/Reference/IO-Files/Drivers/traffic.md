---
title: Traffic driver
---

---

# Traffic Driver `PIDS_TRAFFIC`
The Traffic driver is used by the PALM's Traffic module. It is denoted by the suffix `_traffic`, e.g. `jobname_traffic`. The driver file, if found from the input data directory, will be automatically copied into job work directory by `palmrun`. The structure of the file is similar to VSRC chemistry driver. Timesteps can be of arbitrary variable length, the values are internally interpolated. The grid cells with active traffic are indexed linearly by dimension s. The number of car types can be selected according available data. Each category of the cars is represented by one average properties.     

## Global attributes
Setting global attributes is optional.

| Attribute      | Type     | Description                                                                                                  |
| -------------- | -------- | ------------------------------------------------------------------------------------------------------------ |
| author         | NC_CHAR  | First name, last name, email address.                                                                        |
| campaign       | NC_CHAR  | User-defined text, max. 12 characters.                                                                       |
| comment        | NC_CHAR  | User-defined text.                                                                                           |
| contact_person | NC_CHAR  | First name, last name, email address.                                                                        |
| creation_time  | NC_CHAR  | File creation date (UTC), format: YYYY-MM-DD hh:mm:ss +00.                                                   |
| title          | NC_CHAR  | Short description, e.g., "PALM-Traffic input file for scenario 1b".                                          |
| version        | NC_INT   | E.g. 1.                                                                                                      |

## Dimensions
The dimensions of the Traffic model driver file are defined as follows:

| Dimension name | Type     | Attributes                                                            | Description                                                                                                                                                                        |
|----------------| -------- |-----------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `time`         | NC_INT   | `long_name="number of timesteps"`                                     | Number of timesteps.                                                                                                                                                               |
| `s`            | NC_INT   | `long_name="number of traffic grid cells"`                            | Number of traffic grid cells.                                                                                                                                                      |
| `car_type`     | NC_INT   | `long_name="number of car types"`                                     | Number of car types.                                                                                                                                                               |
| `field_length` | NC_INT   | `long_name="length of the text fields timestamp and car_type"`        | Length of the text fields timestamp and car_type.                                                                                                                                  |



## Variables

| Dimension name                     | Type     | Attributes                                                                            | Description                                             |
|------------------------------------|----------|---------------------------------------------------------------------------------------|---------------------------------------------------------|
| `timestamp(time, field_length)`    | NC_CHAR  | `long_name="Timestamp"`                                                               | Timestamp (format the same as PALM VSR emission).       |
| `car_type(car_type, field_length)` | NC_CHAR  | `long_name="car type description"`                                                    | Car type description.                                   |
| `car_cd(car_type)`                 | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="car drag coefficient"`, `units="1"`              | Car drag coefficient (unitless).                        |
| `car_length(car_type)`             | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="car length"`, `units="m"`                        | Car length (m).                                         |
| `car_width(car_type)`              | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="car width"`, `units="m"`                         | Car width (m).                                          |
| `car_height(car_type)`             | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="car height"`, `units="m"`                        | Car height (m).                                         |
| `i(s)`                             | NC_INT   | `_FillValue=-999`(*), `long_name="i coordinate"`, `units="1"`                         | i coordinate of the grid cell.                          |
| `j(s)`                             | NC_INT   | `_FillValue=-999`(*), `long_name="j coordinate"`, `units="1"`                         | j coordinate of the grid cell.                          |
| `k(s)`                             | NC_INT   | `_FillValue=-999`(*), `long_name="k coordinate"`, `units="1"`                         | k coordinate of the grid cell.                          |
| `width(s)`                         | NC_INT   | `_FillValue=-999`(*), `long_name="width of the street"`, `units="m"`                  | Width of the traffic stream (m).                        |
| `slope(s)`                         | NC_INT   | `_FillValue=-999`(*), `long_name="slope of the street"`, `units="%"`                  | Slope of the traffic stream (%).                        |
| `dirx(s)`                          | NC_INT   | `_FillValue=-999`(*), `long_name="car movement direction x"`, `units="1"`             | x coordinate of car movement direction unit vector (1). |
| `diry(s)`                          | NC_INT   | `_FillValue=-999`(*), `long_name="car movement direction y"`, `units="1"`             | y coordinate of car movement direction unit vector (1). |
| `frac(s)`                          | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="fraction coefficient"`, `units="1"`              | Fraction of the traffic line in the grid cell.          |
| `intensity(time, s, car_type)`     | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="traffic intensity"`, `units="1/hour"`            | Traffic intensity (cars per hour).                      |
| `heat(time, s, car_type)`          | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="heat produced by car"`, `units="W"`              | Heat produced by car (W per car).                       |
| `speed(time, s, car_type)`         | NC_FLOAT | `_FillValue=-9999.f`(*), `long_name="speed of the car"`, `units="m/s"`                | Speed of the car (m/s).                                 |
 


## References
