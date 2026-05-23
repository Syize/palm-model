# YAML configuration

Reference of the YAML configuration parameters

---

The YAML configuration file of `palm_csd` consists of the following sections, each serving a specific purpose:

- [`attributes`](#attributes-section): Optional global attributes that will be passed to the static driver file.
- [`settings`](#settings-section): Settings that are applied to all domains.
- [`output`](#output-section): Output settings for the static driver file(s).
- [`input`](#input-sections): One or several input data sets for one or several domains.
- [`domain`](#domain-sections): One or several domains to generate by `palm_csd`, each with its own settings.
- [`lcz`](#lcz-section): Local Climate Zone (LCZ) parameter adjustments.

## `attributes` section

This section includes optional global attributes that will be passed to the static driver file. Note that these global attributes have no effect on the PALM simulations.

### `acronym`

Entry type
: string

Default value
: `None`

Institutional acronym

### `author`

Entry type
: string

Default value
: `None`

Author of the static driver. Use the format: name, email

### `campaign`

Entry type
: string

Default value
: `None`

Information on measurement capaign (if applicable)

### `comment`

Entry type
: string

Default value
: `None`

Arbitrary text

### `contact_person`

Entry type
: string

Default value
: `None`

Contact person, format as for [`author`](#author)

### `data_content`

Entry type
: string

Default value
: `None`

Arbitrary text

### `dependencies`

Entry type
: string

Default value
: `None`

Arbitrary text

### `institution`

Entry type
: string

Default value
: `None`

Institution of the driver creator

### `keywords`

Entry type
: string

Default value
: `None`

Arbitrary keywords

### `location`

Entry type
: string

Default value
: `None`

Geo-location of the static driver content (if applicable)

### `origin_time`

Entry type
: string

Default value
: `None`

Reference point in time, format: `YYYY­-MM­-DD hh:mm:ss ZZZ`, e.g. `2000-01-01 11:00:00 +01` (1st January 2000, 11 am Central European Time)

### `references`

Entry type
: string

Default value
: `None`

Arbitrary text

### `site`

Entry type
: string

Default value
: `None`

Site description of the static driver content (if applicable)

### `source`

Entry type
: string

Default value
: `None`

List of data sources used to generate the driver

## `settings` section

This section includes global parameters used to create the static driver(s).

### `downscaling_method`

Entry type
: string or dictionary

Default value
: *see description*

Resampling algorithms for downscaling GeoTIFF input. Can be set for typically `categorical`, `continuous`, `discontinuous` and `discrete` variables to the string of an [algorithm supported by rasterio](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling).

### `epsg`

Entry type
: integer

Default value
: `None`

EPSG code of the coordinate reference system (CRS) of the output and the PALM simulation. Currently, only UTM CRSs were tested. If `None`, all netCDF coordinate input files in the `input` section have to be provided.

### `upscaling_method`

Entry type
: string or dictionary

Default value
: *see description*

Resampling algorithms for upscaling GeoTIFF input. Can be set for typically `categorical`, `continuous`, `discontinuous` and `discrete` variables to the string of an [algorithm supported by rasterio](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling).

### `ignore_input_georeferencing`

Entry type
: logical

Default value
: `False`

When reading GeoTIFF input, ignore its coordinate reference system (CRS), resulting in a similar behaviour as with netCDF input. In particular, both, `input_lower_left_x` and `input_lower_left_y` need to be set for each domain.

### `height_high_vegetation_lower_threshold`

Entry type
: float

Default value
: `2.0`

Unit
: `m`

Lower threshold of vegetation height for high vegetation. Vegetation with height < `height_high_vegetation_lower_threshold` is considered low vegetation.

### `height_rel_resolved_vegetation_lower_threshold`

Entry type
: float

Default value
: `0.5`

Lower threshold of vegetation height relative to the vertical grid spacing for resolved vegetation. Only vegetation (single trees and vegetation patches) with height > `height_rel_resolved_vegetation_lower_threshold` * `dz` is considered to be converted to resolved vegetation.

### `lai_high_vegetation_default`

Entry type
: float

Default value
: None

Data unit
: m<sup>2</sup> m<sup>-2</sup>

Default leaf area index for high vegetation according to [`height_high_vegetation_lower_threshold`](#height_high_vegetation_lower_threshold) used to generate the 3D leaf area density field. This value is used for all pixels for which no other leaf area density is available.

### `lai_low_vegetation_default`

Entry type
: float

Default value
: None

Data unit
: m<sup>2</sup> m<sup>-2</sup>

Default leaf area index for low vegetation according to [`height_high_vegetation_lower_threshold`](#height_high_vegetation_lower_threshold) used to generate the 3D leaf area density field. This parameter is also used for the `vegetation_type`-related LAI value in the `vegetation_pars` field. It is only used for all pixels for which no other leaf area density is available.

### `lai_tree_lower_threshold`

Entry type
: float

Default value
: `0.0`

Data unit
: m<sup>2</sup> m<sup>-2</sup>

Lower threshold of LAI for trees. Trees with LAI < `lai_tree_lower_threshold` are either removed or considered to have LAI = `lai_tree_lower_threshold`, depending on the setting `remove_low_lai_tree`.

### `lad_method`

Entry type
: string

Default value
: `Metal2003`

Approach to reconstruct the vertical LAD profiles for vegetation patches (parks, forests), where the canopy can be considered to be pseudo-1D and for which usually no information on individual trees is available. A tool to visualize the approaches is described below. Currently, `Metal2003` for Markkanen et al. (2003) and `LM2004` for Lalic and Mihailovic (2004) are supported. For the former, the parameters `lad_alpha` and `lad_beta` are considered, and for the latter `lad_z_max_rel`.

### `lad_alpha`

Entry type
: float

Default value
: `5.0`

Parameter for reconstruction of vertical LAD profiles based on tree shape parameters (alpha, beta) and the integral leaf area index after Markkanen et al. (2003). A tool to visualize the effect of this parameter is described below.

### `lad_beta`

Entry type
: float

Default value
: `3.0`

Parameter for reconstruction of vertical LAD profiles based on tree shape parameters (alpha, beta) and the integral leaf area index after Markkanen et al. (2003). A tool to visualize the effect of this parameter is described below.

### `lad_z_max_rel`

Entry type
: float

Default value
: `0.7`

Parameter for reconstruction of vertical LAD profiles after Lalic and Mihailovic (2004) togtether with the integral leaf area index. It represents the height of the maximum LAD relative to the patch height (zm/h). A tool to visualize the effect of this parameter is described below.

### `patch_height_default`

Entry type
: float

Default value
: `10.0`

Unit
: `m`

Default patch height, which is used in the canopy generator to process canopy patches (parks, forests) for which data for individual trees is usually lacking. This parameter comes into affect for data gaps where no other vegetation height is available.

### `replace_invalid_input_values`

Entry type
: logical

Default value
: `True`

If `True`, replace invalid input values. Currently, this includes replacing non-missing values that are outside of the valid value range as defined by `palm_csd/data/value_defaults.csv`, and using the default vegetation type when no other surface type or building is set for a pixel.

### `season`

Entry type
: string

Default value
: `summer`

As palm_csd can work with different sets of input data regarding leaf area index, this switch parameter can be set to either `summer` or `winter` to select the most suitable leaf area index input file to account for differences in leaf amount. Data for summer is usually from August (fully leaved), while data for winter is usually from April.

### `soil_type_default`

Entry type
: integer

Default value
: `3`

Default soil type used to fill data gaps in the soil type distribution.

### `rotation_angle`

Entry type
: float

Default value
: `0.0`

Unit
: `degrees`

Rotation angle of the model's North direction relative to geographical North (clockwise rotation). Only values between 0 and 360 are valid. This value overwrites the namelist parameter of the PALM run.

### `use_lai_for_roofs`

Entry type
: logical

Default value
: `True`

If set to `True`, the general LAI input is used for vegetation on roofs if the specific [`building_lai`](#building_lai) input for roof surfaces is not set.

### `use_lai_for_trees`

Entry type
: logical

Default value
: `True`

If set to `True`, the general LAI input is used as single tree LAI for each tree with undefined tree LAI. Only if both, the tree LAI and the general LAI, are undefined, the species-dependent default value is used. If set to `False`, only the tree LAI and the default values define the LAI of a single tree.

### `use_vegetation_height_for_trees`

Entry type
: logical

Default value
: `True`

If set to `True`, the vegetation height is used as single tree LAI for each tree with undefined tree height. Only if both, the tree height and the vegetation height, are undefined, the species-dependent default value is used. If set to `False`, only the tree height and the default values define the height of a single tree.

### `vegetation_type_below_trees`

Entry type
: integer

Default value
: `3`

If resolved vegetation pixels are added above high vegetation type pixels, the vegetation type pixels are changed to this value.

## `output` section

This section describes the location for the static driver output.

### `file_out`

Entry type
: string

Output file name. The final output will be stored under `path`/`file_out`_`domain`, where `domain` will be "root" for the parent (root) domain, and "N01", "N02", etc., for child domains N01, N02, etc., respectively. This parameter is mandatory.

### `path`

Entry type
: string

Default value
: `None`

Directory where the output file shall be stored. Note that the static driver can - depending on model domain size - be quite large (in the order of several GB).

### `version`

Entry type
: integer

Default value
: `None`

User-specific setting to track updates of a static driver. This value will be added as global attribute to the static driver.

## `input` section(s)

The configuration YAML can include several sets of input data for different domains. For each set of input data, an individual section must be provided and named accordingly (i.e. `input_root`, `input_N02`, etc.). If there is only one input data set, a name can also be omitted (i.e. `input`). One `input` section can be used by several domains.

### `columns`

Entry type
: dictionary

Mapping of the data columns in the [`surfaces`](#surfaces) or [`trees`](#trees) input files to the model variables. The entries can have two forms: 1) A direct mapping from an input column to a [`palm_csd` input variable](types.md#input-data). In this case, the whole input column includes only values of that specific input variable. Or 2), an input column followed by a mapping from a input column value to a specific palm_csd input variable value.

### `files`

Entry type
: dictionary

This dictionary contains the file names of the 2d input files. The keys are the names of the input files, and the values are the file paths relative to the [`path`](yaml.md#path-1) parameter. Both, [`surfaces`](#surfaces) and [`trees`](#trees) include one or a list of vector files (Shapefiles). All the other files are raster files (georeferenced such as GeoTIFF or netCDF).

#### `bridges_2d`

Entry type
: string

Data unit
: m

Bridge height. This is the height of the upper surface of the bridge-like structure. Its thickness is set in [`bridge_depth`](yaml.md#bridge_depth).

#### `bridges_id`

Entry type
: string

Bridge ids. Used to identify connected bridge-like structures in the model domain.

#### `building_id`

Entry type
: string

Building ids. Used to identify connected buildings in the model domain.

#### `building_type`

Entry type
: string

[Building type values](types.md#building-type).

#### `buildings_2d`

Entry type
: string

Data unit
: m

Building height.

#### `lai`

Entry type
: string

Data unit
: m<sup>2</sup> m<sup>-2</sup>

Leaf area index.

#### `lat`

Entry type
: string

Data unit
: degrees N

Latitude. DEPRECATED, let `palm_csd` calculate it using [`epsg`](#epsg) and [`origin_lon`](#origin_lon)/[`origin_lat`](#origin_lat) or [`origin_x`](#origin_x)/[`origin_y`](#origin_y).

#### `lcz`

Entry type
: string

Local Climate Zone classification. It could either consist of one Band indicating the LCZ class with values 1-17 or three bands with values indicating RGB colours with values 0-255. The RGB values can be customized in the [`lcz`](#lcz-section) section.

#### `lon`

Entry type
: string

Data unit
: degrees E

Longitude. DEPRECATED, let `palm_csd` calculate it using [`epsg`](#epsg) and [`origin_lon`](#origin_lon)/[`origin_lat`](#origin_lat) or [`origin_x`](#origin_x)/[`origin_y`](#origin_y).

#### `patch_type`

Entry type
: string

Tree type value of vegetation patches.

#### `pavement_type`

Entry type
: string

[Pavement type value](types.md#pavement-type).

#### `soil_type`

Entry type
: string

[Soil type value](types.md#soil-type).

#### `surfaces`

Entry type
: string *or* list of strings

Polygon shapefile surface data. One or several files can be given. Only polygons with a defined value in a data column with the appropriate mapping in `columns` will be considered.

#### `street_crossings`

Entry type
: string

Positions of street crossings that can be used by pedestrians (used for multi-agent model). `1` for crossings, a missing value otherwise.

#### `street_type`

Entry type
: string

[Street type value](types.md#street-type) (used for parameterized chemistry emissions and multi-agent model).

#### `tree_crown_diameter`

Entry type
: string

Data unit
: m

Tree crown diameter for single trees. For each tree, only one value can be given at the center of the tree location.

#### `tree_height`

Entry type
: string

Data unit
: m

Tree height for single trees. For each tree only one value can be given at the center of the tree location.

#### `tree_lai`

Entry type
: string

Data unit
: m<sup>2</sup> m<sup>-2</sup>

LAI specifically for single trees. For each tree only one value can be given at the center of the tree location.

#### `tree_shape`

Entry type
: string

Shape of single trees. For each tree only one value can be given at the center of the tree location.

#### `tree_trunk_diameter`

Entry type
: string

Data unit
: m

Trunk diameter at breast height for single trees. For each tree only one value can be given at the center of the tree location.

#### `tree_type`

Entry type
: string

Tree type according to the canopy generator tree inventory for single trees. For each tree only one value can be given at the center of the tree location.

#### `tree_type_name`

Entry type
: string

Tree type name according to the canopy generator tree inventory for single trees.

#### `trees`

Entry type
: string *or* list of strings

Point shapefile tree data. One or several files can be given. Every point represents a single tree. The columns are defined in the `columns`.

#### `vegetation_height`

Entry type
: string

Data unit
: m

Vegetation height.

#### `vegetation_type`

Entry type
: string

[Vegetation type values](types.md#vegetation-type).

#### `water_temperature`

Entry type
: string

Data unit
: K

Water temperature.

#### `water_type`

Entry type
: string

[Water type value](types.md#water-type).

#### `x_utm`

Entry type
: string

Data unit
: m

UTM x-coordinates. DEPRECATED, let `palm_csd` calculate it using [`epsg`](#epsg) and [`origin_lon`](#origin_lon)/[`origin_lat`](#origin_lat) or [`origin_x`](#origin_x)/[`origin_y`](#origin_y).

#### `y_utm`

Entry type
: string

Data unit
: `m`

UTM y-coordinates. DEPRECATED, let `palm_csd` calculate it using [`epsg`](#epsg) and [`origin_lon`](#origin_lon)/[`origin_lat`](#origin_lat) or [`origin_x`](#origin_x)/[`origin_y`](#origin_y).

#### `zt`

Entry type
: string

Data unit
: m

Terrain height.

### `path`

Entry type
: string

Directory where the input files reside.

## `domain` section(s)

This section contains settings for each model domain for the PALM run. If the name is omitted, the name `root` is assumed. In case of a nested run, the sections for the non-root domains must be named individually, e.g. `domain_N01`, `domain_N02`, etc. as it is done in the PALM parameter file.

The corresponding input data set for each domain can be defined in the `input` parameter. If not set, the name of both `input` and `domain` sections is used to find the matching input data set (e.g. `input_root` for `domain_root`). If there is only one `input` section, this is used.

If geographical coordinates of the output should be calculated, i.e. if they are not supplied in the input data with `file_x_UTM`, `file_y_UTM` etc., it is sufficient to either set `origin_x`/`origin_y` or `origin_lon`/`origin_lat`. Note that also `epsg` must be set in the `settings` section.

### `bridge_depth`

Entry type
: float

Default value
: `3.0`

Unit
: m

Vertical depth or thickness of all bridge elements in the domain. Bridges are treated as building grid cells in `buildings_3d`. The values in [`bridges_2d`](#bridges_2d) define for each pixel the maximum height of these grid cells and `bridge_depth` defines how far these grid cells extend downwards.

### `building_albedo_type`

Entry type
: dictionary

Albedo type for all buildings in the domain. Can be set for any [`building_surface_type`](types.md#building_surface_type).

### `building_emissivity`

Entry type
: dictionary

Entry unit
: 1

Emissivity for all buildings in the domain. Can be set for any [`building_surface_type`](types.md#building_surface_type).

### `building_fraction`

Entry type
: dictionary

Entry unit
: 1

Building surface fractions for all buildings in the domain. Can be set for any [`building_surface_type`](types.md#building_surface_type).

### `building_free_border_width`

Entry type
: float

Default value
: `0.0`

Entry unit
: m

Width of the clearance zone at the border of a the domain that is free of buildings. Buildings are replaced by pavement of type [`building_free_border_pavement_type`](#building_free_border_pavement_type).

### `building_free_border_pavement_type`

Entry type
: string or integer

Default value:
: `asphalt_concrete_mix`

[Pavement type](types.md#pavement_type) used as a replacement for buildings in the clearance zone defined by [`building_free_border_width`](#building_free_border_width).

### `building_general_pars`

Entry type
: dictionary

[General parameters](types.md#building_general_pars) for all the buildings in the domain.

### `building_heat_capacity`

Entry type
: dictionary

Entry unit
: J m<sup>-3</sup> K<sup>-1</sup>

Heat capacity of the urban surface layers for all buildings in the domain. Can be set for any [`building_surface_type`](types.md#building_surface_type) and any layer.

### `building_heat_conductivity`

Entry type
: dictionary

Entry unit
: W m<sup>-1</sup> K<sup>-1</sup>

Heat conductivity of the urban surface layers for all buildings in the domain. Can be set for any [`building_surface_type`](types.md#building_surface_type) and any layer.

### `building_indoor_pars`

Entry type
: dictionary

[Indoor parameters](types.md#building_indoor_pars) for all buildings in the domain.

### `building_lai`

Entry type
: dictionary

Entry unit
: m<sup>2</sup> m<sup>-2</sup>

LAI of the green fraction of the urban surfaces for all buildings in the domain. Can be set for any [`building_surface_level`](types.md#building_surface_level).

### `building_roughness_length`

Entry type
: dictionary

Entry unit
: m

Roughness length of the urban surfaces for all buildings. Can be set for any [`building_surface_level`](types.md#building_surface_level).

### `building_roughness_length_qh`

Entry type
: dictionary

Entry unit
: m

Roughness length for heat and moisture for all urban surfaces. Can be set for any [`building_surface_level`](types.md#building_surface_level).

### `building_thickness`

Entry type
: dictionary

Entry unit
: m

Layer thickness of the urban surfaces for all buildings in the domain. Can be set for any [`building_surface_type`](types.md#building_surface_type) and any level.

### `building_transmissivity`

Entry type
: dictionary

Entry unit
: 1

Window transmissivity of the urban surface for all buildings in the domain. Can be set for any [`building_surface_level`](types.md#building_surface_level).

### `buildings_3d`

Entry type
: logical

Default value
: `False`

Use 3D buildings via the `buildings_3d` array instead of `buildings_2d`. If bridges are present in the simulation domain, `buildings_3d` is generated in any case.

### `dcep`

Entry type  
: logical

Default value  
: `True`

Generate urban parameters for DCEP.

### `domain_parent`

Entry type
: string

Default value
: `None`

Name of the parent domain of the current domain. If the current domain is the root domain, do not set this parameter.

### `dz`

Entry type
: float

Unit
: `m`

Vertical grid spacing in PALM. This parameter is needed when `buildings_3d`, `generate_single_trees`, `generate_vegetation_patches`, `interpolate_terrain`, or `use_palm_z_axis` is used. This parameter is mandatory.

### `estimate_lai_from_vegetation_height`

Entry type
: logical

Default value
: `True`

If set to `True`, the LAI is estimated as [`lai_per_vegetation_height`](#lai_per_vegetation_height) * `vegetation_height`. This is only applied for pixels without LAI information.

### `generate_single_trees`

Entry type
: logical

Default value
: `True`

If set to `True`, information on individual (single) trees will be used to generate a 3D leaf area density and basal area density distribution for each tree. In contrast to vegetation patches, where a closed canopy is assumed and information is only distributed vertically for each pixel, single trees have a 3D shape that is mapped on the simulation domain.

### `generate_vegetation_patches`

Entry type
: logical

Default value
: `True`

If set to `True`, the embedded canopy generator will convert all surface pixels that contain high vegetation into a 3D leaf area density distribution. This applies to pixels where the patch height or the patch type are defined. If `replace_high_vegetation_types` is `True`, this applies also where `vegetation_type` is set to a high vegetation type. Note that only pixels with either undefined patch heights or patch heights `> 0.5*dz` are converted, while all other pixels will be parameterized via the `vegetation_type` field. (The height threshold can be adjusted with `height_rel_resolved_vegetation_lower_threshold`.)

### `input`

Entry type
: string

Default value
: `None`

Name of the `input` section to be used for this domain. This parameter is used to match the input data set with the domain. If not set, the name of both `input` and `domain` sections or the `pixel_size` parameter (deprecated) is used to find the matching input data set. If there is only one `input` section, this is used.

### `input_lower_left_x`

Entry type
: float

Default value
: `None`

Unit
: `m`

Distance along x-direction between the lower-left corner of the model domain and the lower-left corner of the input data. This parameter is used to shift the model domain with respect to the provided input data. Only needed for netCDF input data.

### `input_lower_left_y`

Entry type
: float

Default value
: `None`

Unit
: `m`

Distance along y-direction between the lower-left corner of the model domain and the lower-left corner of the input data. This parameter is used to shift the model domain with respect to the provided input data. Only needed for netCDF input data.

### `interpolate_terrain`

Entry type
: logical

Default value
: `False`

If set to `True`, the terrain height is interpolated and blended over between parent and child domains in order to avoid severe steps in terrain height due to different grid spacings between parent and child.

### `lai_per_vegetation_height`

Entry type
: float

Default value
: 0.2

Data unit
: m<sup>2</sup> m<sup>-3</sup>

Leaf area index per vegetation height used to estimate LAI from vegetation height for pixels without LAI information when [`estimate_lai_from_vegetation_height`](#estimate_lai_from_vegetation_height) is set to `True`.

### `lcz_input`

Entry type  
: string

Default value  
: `None`

Output parameters derived from LCZ data. Currently, only `None` and `full` are supported. `full` means that all parameters are derived from LCZ data except orography height.

### `lower_left_x`

Entry type
: float

Default value
: `None`

Unit
: `m`

Only for nested domains: Distance along x-direction between the lower-left corner of the nested domain and the lower-left corner of the root parent domain. This parameter is used to define the coordinates of origin of the nested domain. This parameter is not required if the origin is defined via `origin_x`/`origin_y` or `origin_lon`/`origin_lat`.

### `lower_left_y`

Entry type
: float

Default value
: `None`

Unit
: `m`

Only for nested domains: Distance along y-direction between the lower-left corner of the nested domain and the lower-left corner of the root parent domain. This parameter is used to define the coordinates of origin of the nested domain. This parameter is not required if the origin is defined via `origin_x`/`origin_y` or `origin_lon`/`origin_lat`.

### `nx`

Entry type
: integer

Number of grid points in x-direction. It equals the `nx` setting in the PALM parameter file so the actual number of grid points is `nx+1`. This parameter is mandatory.

### `ny`

Entry type
: integer

Number of grid points in y-direction. It equals the `ny` setting in the PALM parameter file so the actual number of grid points is `ny+1`. This parameter is mandatory.

### `origin_lat`

Entry type
: float

Default value
: `None`

Unit
: `degrees N`

Latitude of the lower border of the lower-left grid point of the PALM domain in WGS84.

### `origin_lon`

Entry type
: float

Default value
: `None`

Unit
: `degrees E`

Longitude of the left border of the lower-left grid point of the PALM domain in WGS84.

### `origin_x`

Entry type
: float

Default value
: `None`

Unit
: `m`

x-coordinate of the left border of the lower-left grid point of the PALM domain in the CRS defined by `epsg` in the `settings` section.

### `origin_y`

Entry type
: float

Default value
: `None`

Unit
: `m`

y-coordinate of the lower border of the lower-left grid point of the PALM domain in the CRS defined by `epsg` in the `settings` section.

### `overhanging_trees`

Entry type
: logical

Default value
: `True`

If set to `False`, no LAD volumes of trees are generated above surfaces without a vegetation type.

### `pixel_size`

Entry type
: float

Unit
: `m`

Size of a single pixel in x/y direction (equal to grid spacing in x and y). This parameter is mandatory.

### `remove_low_lai_tree`

Entry type
: logical

Default value
: `False`

If set to `True`, all trees with an LAI < `lai_tree_lower_threshold` are removed from the dataset. If set to `False`, those trees are considered with LAI = `lai_tree_lower_threshold`.

### `replace_high_vegetation_types`

Entry type
: logical

Default value
: `True`

If set to `True` and if `generate_vegetation_patches` is `True`, pixels where a high vegetation type was prescribed will be converted into a 3D leaf area density canopy using the canopy generator for undefined or defined patch heights `> 0.5*dz`. (This threshold can be adjusted with `height_rel_resolved_vegetation_lower_threshold`.) For lower patch heights, the respective pixels will be converted to short grass. If `generate_vegetation_patches` is `False`, all high vegetation pixels are converted to short grass. If set to `False`, it is allowed to have unresolved high vegetation classes according in the `vegetation_type` distribution. Still, after optional vegetation patch generation, high vegetation pixels with a vegetation height < 2.0 m are replaced by short grass. (This threshold can be adjusted with `height_high_vegetation_lower_threshold`.) Note that high vegetation classes can involve very large roughness lengths > 0.5 m. If the vertical grid spacing is close to or smaller than this threshold the PALM run will crash and/or does not provide meaningful results. It is generally recommended to set this parameter to `True` whenever the grid spacing is small enough to resolve canopy patches by 2 or more vertical grid levels.

### `udir`

Entry type  
: float

Default value  
: `[0.0, 90.0]`

Urban street directions (degrees).

### `use_palm_z_axis`

Entry type
: logical

Default value
: `False`

If set to `True`, the static driver will raster the input data on the z-grid of PALM for output. Note that PALM will convert continuous static driver data itself on its grid and apply additional filtering procedures. It is thus recommended to set this parameter to `False` unless `interpolate_terrain: True` in nested set-ups.

### `water_temperature`

Entry type
: float or dictionary

Default value
: `None`

Unit
: `K`

Water temperature in K for one or several water types as indicated by their name or their index 0 to 5. Also allows one value, which is applied to all water types.

### `z_uhl`

Entry type  
: float

Default value  
: `[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0]`

Unit
: `m`

Height of the urban DCEP layers.

## `lcz` section

This section contains settings to overwrite the default values assigned to each LCZ. The values consists of the parameters assigned by the LCZ classification scheme and the PALM parameters. If it is a LCZ-related value, it is checked if it is within the defined valid range. Furthermore, the kind of average used when interpreting the average building height can be set.

### `height_geometric_mean`

Entry type
: logical

Default value
: `True`

Use the geometrical mean (`True`) or the arithmetic mean (`False`) when calculating the building height distribution for DCEP.

### *class*

Entry type
: dictionary

Default value
: `None`

Setting of the parameters of the LCZ *class* in the form of `type: value`. Use the names of classes and parameters as used in the tables in [the technical description of the LCZ approach](lcz_dcep.md#definition-of-parameter-values).
