# PALM Create Static Driver (`palm_csd`)

`palm_csd` is a Python 3 tool to generate the static driver for PALM. It offers two modes of operation for different application purposes of PALM:

1) building-resolving simulations and
2) non-building-resolving simulations with parameterized buildings.

For (1), rasterized input data for the urban, vegetation and soil properties are required while for (2) most of the properties can be derived from a Local Climate Zones map. The Local Climate Zone (LCZ) classification (Stewart and Oke 2012)[^stewart2012] consists of the 17 classes with 10 built and 7 land cover types.

## Execution

In order to create a static driver, follow the installation process described in the `README.md` and prepare a `palm_csd` configuration file in the **YAML format** (see below for details). After that, execute `palm_csd` as follows:

```bash
palm_csd <path/to/csd-config.yml>
```

`palm_csd` supports the following command line switches:

* `-h`, `--help`: Show the help message and exit.
* `--png`: Store a plot as PNG by adding ".png" to the output name.
* `--pdf`: Store a plot as PDF by adding ".pdf" to the output name.
* `-s,` `--show`: Show a plot after the generation of the static driver.
* <code>-v *part*</code>, <code>--verbose *part*</code>: Produce additional output in selected processing parts. Note that this switch always tries to consume one argument, thus, it cannot be used without argument directly in front of the input YAML configuration. The following <code>*part*</code>s are available:
  * `all`: All of processing parts described below. This is equivalent to `-v`/`--verbose` without argument.
  * `gis`: GIS-related operations. In particular, the results of GeoTIFF reprojection, resampling and cutting are saved for inspection in the output directory.
  * `io`: Input/Output operations.
  * `misc`: Miscellaneous operations not falling in the other categories.
  * `vegetation`: Generation of resolved LAD/BAD fields.

Here are some examples of how to use the command line switches:

```bash
palm_csd <path/to/csd-config.yml>               # standard output, no plots
palm_csd -v all <path/to/csd-config.yml>        # verbose all parts
palm_csd <path/to/csd-config.yml> -v            # verbose all parts
palm_csd -v -s <path/to/csd-config.yml>         # verbose all parts, show plot
palm_csd -v gis -v io <path/to/csd-config.yml>  # verbose GIS and IO part
palm_csd -s --pdf <path/to/csd-config.yml>      # show plot and store in PDF
```

The static driver will be written to the directory specified in the configuration file. During compilation of the driver, `palm_csd` will print some information to screen.

## Configuration and technical description

The configuration and the technical description of `palm_csd` is described separately for the two modes of application:

* [Detailed input for building-resolving simulations](detailed_input.md)
* [Local Climate Zone based input for non-building-resolving simulations](lcz_dcep.md)

## Static driver statistics and visualization tool

This package comes also with `static_driver_stats`, a tool calculate statistics of and visualize a static driver file in netCDF format. It is called as follows:

```bash
static_driver_stats <path/to/static_driver.nc>
```

It supports the following command line switches:

* `-h`, `--help`: Show the help message and exit.
* `--png`: Store a plot as PNG by adding ".png" to the static driver name.
* `--pdf`: Store a plot as PDF by adding ".pdf" to the static driver name.
* `-s,` `--show`: Show a plot of the static driver.
* <code>-H *HEIGHT*</code>, <code>--height *HEIGHT*</code>: Plot height in inches.
* <code>-W *WIDTH*</code>, <code>--width *WIDTH*</code>: Plot width in inches.
* <code>-t *TITLE*</code>, <code>--title *TITLE*</code>: Plot title.

Without explicit setting of the plot's width and height, a relationship of
$$
\textrm{width} = \textrm{aspect ratio} \cdot \textrm{height} + 2\,\textrm{inch}
$$
is assumed, where the aspect ratio is determined by the ratio of the number pixels in x and y direction in the static driver.

Here are some examples of how to use the command line switches:

```bash
static_driver_stats file.nc                          # print only statistics
static_driver_stats -s file.nc                       # statistics and shown plot
static_driver_stats --pdf -t "My Domain" config.yml  # show plot and store in PDF
```

## Developer guide

A description of the coding style and the automated testing is provided in the [developer guide](developer_guide.md).

[^stewart2012]: Stewart, I. D., and T. R. Oke. 2012. ‘Local Climate Zones for Urban Temperature Studies’. Bulletin of the American Meteorological Society 93 (12): 1879–1900. [doi: 10.1175/BAMS-D-11-00019.1](https://doi.org/10.1175/BAMS-D-11-00019.1).
