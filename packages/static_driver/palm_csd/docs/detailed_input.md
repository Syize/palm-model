# Detailed input for building resolving simulations

In the following, it is described how `palm_csd` can be used to create a static driver for building resolving simulations. Detailed, rasterized input data for the urban, vegetation and soil properties are required.

## Configuration file

This section describes how to set-up a configuration file for creating a static driver for PALM based on pre-processed input data. This data can be supplied as netCDF or any GIS raster format (GeoTIFF recommended and tested). When using GIS data input, `palm_csd` is able to reproject or align the data with the target grid. If the GIS input is already aligned with the target grid, it is read and cut without reprojection.

In the following, we will use exemplary data for Berlin, which is available via Open Access. The configuration file uses the YAML format. All variables with a default value can be omitted in the configuration file. Note that the `None` value of Python, which represents a non-defined value, is represented in the YAML file by `null`.

The configuration file consists of the following sections:

### `attributes` section

A set of global attributes can be defined that will be passed to the static driver file. The following attributes can be set:

| Variable         | Data type | Default value | Description|
|------------------|-----------|---------------|------------|
| `author`         | string | `None` | Author of the static driver. Use the format: name, email |
| `contact_person` | string | `None` | Contact person, format as for `author` |
| `acronym`        | string | `None` | Institutional acronym |
| `comment`        | string | `None` | Arbitrary text |
| `data_content`   | string | `None` | Arbitrary text |
| `dependencies`   | string | `None` | Arbitrary text |
| `keywords`       | string | `None` | Arbitrary keywords |
| `source`         | string | `None` | List of data sources used to generate the driver |
| `campaign`       | string | `None` | Information on measurement capaign (if applicable) |
| `location`       | string | `None` | Geo-location of the static driver content (if applicable) |
| `site`           | string | `None` | Site description of the static driver content (if applicable) |
| `institution`    | string | `None` | Institution of the driver creator |
| `references`     | string | `None` | Arbitrary text |
| `palm_version`   | float  | `None` | PALM version for which the driver was generated (for compatibility checks) |
| `origin_time`    | string | `None` | Reference point in time, format: `YYYY­-MM­-DD hh:mm:ss ZZZ`, e.g. `2000-01-01 11:00:00 +01` (1st January 2000, 11 am Central European Time) |

Note that these global attributes have no effect on the PALM simulations. Consequently, all attributes that are not explicitly set in the configuration file are omitted in the static driver.

Example:

```yml
attributes:
  author: Bjoern Maronga, maronga@muk.uni-hannover.de
  contact_person: Bjoern Maronga, maronga@muk.uni-hannover.de
  acronym: LUHimuk
  comment: created with palm_csd
  location: B
  site: Berlin Mitte
  institution: Leibniz University Hannover, Institute of Meterology and Climatology
  palm_version: 6.0
```

### `settings` section

This section describes global parameters used to create the static driver.

| Variable | Data type | Default value | Description |
|----------|-----------|---------------|-------------|
| `downscaling_method`      | string or dictionary  | *see description* | Resampling algorithms for downscaling GeoTIFF input. Can be set for typically `categorical`, `continuous`, `discontinuous` and `discrete` variables to the string of an [algorithm supported by rasterio](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling). See below for a detailed description.  |
| `upscaling_method`      | string or dictionary  | *see description* | Resampling algorithms for upscaling GeoTIFF input. Can be set for typically `categorical`, `continuous`, `discontinuous` and `discrete` variables to the string of an [algorithm supported by rasterio](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling). See below for a detailed description.  |
| `epsg`                        | integer | `None` | EPSG code of the coordinate reference system (CRS) of the output and the PALM simulation. Currently, only UTM CRSs were tested. If `None`, all netCDF coordinate input files in the `input` section have to be provided. |
| `ignore_input_georeferencing` | logical | `False` | When reading GeoTIFF input, ignore its coordinate reference system (CRS), resulting in a similar behaviour as with netCDF input. In particular, both, `input_lower_left_x` and `input_lower_left_y` need to be set for each domain. |
| `lai_roof_extensive`          | float   |  `0.8` | Leaf are index for green roofs with extensive vegetation, defined by setting the appropriate `building_pars` field. The value is assigned to all extensive green roofs in the model domain. |
| `lai_roof_intensive`          | float   |  `2.5` | Leaf are index for green roofs with intensive vegetation, defined by setting the appropriate `building_pars` field. The value is assigned to all intensive green roofs in the model domain. |
| `lai_high_vegetation_default` | float   |  `6.0` | Default leaf area index for (high) vegetation used to generate the 3D leaf area density field. This value is used for all pixels for which no other leaf area density is available (i.e. to fill missing data). |
| `lai_low_vegetation_default`  | float   |  `1.0` | Default leaf area index for (low) vegetation used to fill data gaps in the leaf area index distribution. This parameter only will a LOD2 leaf area index for parameterized vegetation via vegetation_type, i.e. through the `vegetation_pars` field. |
| `lai_tree_lower_threshold`    | float   |  `0.0` | Lower threshold of LAI for trees. Trees with LAI < `lai_tree_lower_threshold` are either removed or considered to have LAI = `lai_tree_lower_threshold`, depending on the setting `remove_low_lai_tree`. |
| `lad_method`                  | string  | `Metal2003` | Approach to reconstruct the vertical LAD profiles for vegetation patches (parks, forests), where the canopy can be considered to be pseudo-1D and for which usually no information on individual trees is available. A tool to visualize the approaches is described below. Currently, `Metal2003` for Markkanen et al. (2003) and `LM2004` for Lalic and Mihailovic (2004) are supported. For the former, the parameters `lad_alpha` and `lad_beta` are considered, and for the latter `lad_z_max_rel`. |
| `lad_alpha`                   | float   |  `5.0` | Parameter for reconstruction of vertical LAD profiles based on tree shape parameters (alpha, beta) and the integral leaf area index after Markkanen et al. (2003). A tool to visualize the effect of this parameter is described below. |
| `lad_beta`                    | float   |  `3.0` | Parameter for reconstruction of vertical LAD profiles based on tree shape parameters (alpha, beta) and the integral leaf area index after Markkanen et al. (2003). A tool to visualize the effect of this parameter is described below.|
| `lad_z_max_rel`               | float   |  `0.7` | Parameter for reconstruction of vertical LAD profiles after Lalic and Mihailovic (2004) togtether with the integral leaf area index. It represents the height of the maximum LAD relative to the patch height (zm/h). A tool to visualize the effect of this parameter is described below. |
| `patch_height_default`        | float   | `10.0` | Default patch height (in m), which is used in the canopy generator to process canopy patches (parks, forests) for which data for individual trees is usually lacking. This parameter comes into affect for data gaps where no other vegetation height is available. |
|`replace_invalid_input_values` | logical | `True` | If `True`, replace invalid input values. Currently, this includes replacing non-missing values that are outside of the valid value range as defined by `palm_csd/data/value_defaults.csv`, and using the default vegetation type when no other surface type or building is set for a pixel.  |
| `season`                      | string  | `summer` | As palm_csd can work with different sets of input data regarding leaf area index, this switch parameter can be set to either `summer` or `winter` to select the most suitable leaf area index input file to account for differences in leaf amount. Data for summer is usually from August (fully leaved), while data for winter is usually from April. |
| `rotation_angle`              | float   | `0.0`  | Rotation angle of the model's North direction relative to geographical North (clockwise rotation). This value overwrites the namelist parameter of the PALM run. |
| `vegetation_type_below_trees` | integer |  `3`   | If trees are added to the static driver, the vegetation type below the tree volumes is changed to this value. |

Example:

```yml
settings:
   downscaling_method:
     categorical: nearest
     continuous: cubic
   lai_roof_extensive: 3.0
   lai_roof_intensive: 1.5
   lai_high_vegetation_default: 5.0
   lai_low_vegetation_default: 1.0
   lai_alpha: 5.0
   lai_beta: 3.0
   patch_height_default: 10.0
   rotation_angle: 0.0
   season: summer
```

#### Choice of scaling methods

With `downscaling_method` and `upscaling_method`, the resampling algorithm can be chosen for the downscaling and upscaling of the input data when reprojecting or changing the grid.

| Type | Examples | Downscaling default | Upscaling default |
| ---- | -------- | ------------------- | ----------------- |
|`categorical` | building type, pavement type, soil type, street type, vegetation type, water type | `nearest` | `mode` |
|`continuous` | terrain height, water temperature | `bilinear` | `average` |
|`discontinuous` | building height, bridge height, leaf area index, vegetation height | `nearest` | `average` |
|`discrete` | single tree properties except tree type | `nearest` | `average` |

Note that the scaling of `discrete` single tree data does not guarantee to preserve single point quantities. Thus, it is recommended to supply this input data on the target grid. In order to preserve the values in the `categorical` data, only `nearest` and `mode` is allowed. For all other data types, all [algorithms supported by rasterio](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling) can be chosen.

The methods can be either set as a dictionary of the form `type: method` or as a single string. If a string is provided, the method is used for all types. If a dictionary is provided, the method is used for the respective type. If a method is not provided for a type, the default method is used.

Note that the different methods handle missing values differently. While `nearest` produces a missing value when the centre of the target pixel is closest to a missing value in the source data, the other methods calculate values as soon as a part of the target pixel is covered by a non-missing pixel in the source data. In order to ensure consistency between the different data types, the missing values of `nearest` are applied to all data types.

#### Visualization of the LAD profile approaches and their parameters

A simple tool to visualize the different approaches to generate the LAD distributions and their parameters is provided in `tools/plot_lad_profiles.py`. In addition to the Python packages installed for `palm_csd`, it requires the `xarray` and `hvplot` packages. The tool can be executed from the command line with

```bash
python3 tools/lad_patch_configurator.py
```

### `output` section

This section describes the location for the static driver output.

| Variable | Data type | Default value | Description |
|----------|-----------|---------------|-------------|
| `path`      | string  | `None` |Directory where the output file shall be stored. Note that the static driver can - depending on model domain size - be quite large (in the order of several GB). |
| `file_out`\*| string  |    |Output file name. The final output will be stored under `path`/`file_out`_`domain`, where `domain` will be "root" for the parent (root) domain, and "N01", "N02", etc., for child domains N01, N02, etc., respectively. |
| `version`   | integer | `None` |User-specific setting to track updates of a static driver. This value will be added as global attribute to the static driver. |

(*) This parameter is mandatory

Example:

```yml
output:
   path: /ldata2/MOSAIK/
   file_out: winter_iop1_test
   version: 1
```

### `input_ABC`, ..., `input_XYZ` sections

The configuration file can include several sets of input data for different domains. For each set of input data, an individual section must be provided and named accordingly (i.e. `input_root`, `input_N02`, etc.). If there is only one input data set, a name can also be omitted (i.e. `input`). The input files must be in netCDF or any GIS raster format (GeoTIFF recommended and tested). When using GIS files, the target coordinate reference system (CRS) must be defined in the `settings` section via `epsg` as well as the lower left corner of the target grid via `origin_x`/`origin_y` or `origin_lon`/`origin_lat` in the `domain` section. Note that `input_lower_left_x` and `input_lower_left_y` in the domain section are not needed if only GIS files are used.

The coordinate inputs `file_x_UTM`/`file_y_UTM` and `file_lon`/`file_lat` are optional. If they are not provided, `epsg` in the `settings` section and `origin_x`/`origin_y` or `origin_lon`/`origin_lat` in the `domain` section must be provided.

| Variable | Data type | Default value | Description |
|----------|-----------|---------------|-------------|
| `path`                      | string | `None` | Directory where the netCDF input files reside. |
| `pixel_size`                | float  | `None` | DEPRECATED: Horizontal target grid spacing (m) of a surface pixel. Only used to find the matching `input` section for each `domain`. Alternatively, name both, `input` and `domain` the same or use the `input` option in the `domain` section.  |
| `file_x_UTM`                | string | `None` | UTM x-coordinates for the simulation domain (m). |
| `file_y_UTM`                | string | `None` | UTM y-coordinates for the simulation domain (m). |
| `file_lat`                  | string | `None` | Latitude (degrees N) for the simulation domain. |
| `file_lon`                  | string | `None` | Longitude (degrees E) for the simulation domain. |
| `file_buildings_2d`         | string | `None` | 2D building height (m). |
| `file_building_id`          | string | `None` | Building ids. |
| `file_building_type`        | string | `None` | Building type distribution. |
| `file_bridges_2d`           | string | `None` | 2D map of bridge height (m). |
| `file_bridges_id`           | string | `None` | Bridge ids. |
| `file_lai`                  | string | `None` | Leaf area index. |
| `file_patch_height`         | string | `None` | 2D distribution of the vegetation canopy height. |
| `file_patch_type`           | string | `None` | 2D distribution of the vegetation type of vegetation patches. |
| `file_pavement_type`        | string | `None` | Pavement type distribution. |
| `file_soil_type`            | string | `None` | Soil type distribution. |
| `file_street_type`          | string | `None` | Street type distribution (used for parameterized chemistry emissions and multi-agent model). |
| `file_street_crossings`     | string | `None` | Street crossings (used for multi-agent model). |
| `file_tree_height`          | string | `None` | Tree height (m) for street trees. For each tree only one value can be given at the center of the tree location. |
| `file_tree_crown_diameter`  | string | `None` | Tree crown diameter (m). For each tree only one value can be given at the center of the tree location. |
| `file_tree_trunk_diameter`  | string | `None` | Trunk diameter at breast height (m). For each tree only one value can be given at the center of the tree location. |
| `file_tree_type`            | string | `None` | Tree type according to the canopy generator tree inventory. For each tree only one value can be given at the center of the tree location. |
| `file_vegetation_type`      | string | `None` | Vegetation type distribution. |
| `file_vegetation_height`    | string | `None` | Vegetation height (m). |
| `file_vegetation_on_roofs`  | string | `None` | 2D distribution of green roofs. Values can range between 0.0 - 1.0. Intensive vegetation is considered for values >= 0.5, while extensive vegetation is assumed for values < 0.5. |
| `file_water_temperature`    | string | `None` | 2D distribution of the vegetation canopy height. |
| `file_water_type`           | string | `None` | Water type distribution. |
| `file_zt`                   | string | `None` | Terrain height (m). |

(*) This parameter is mandatory

For a given target pixel size (i.e. horizontal grid spacing), only one set of input files can be provided. All input data must be two-dimensional (y,x). While the name of the input variable in each file is not prescribed, ensure that only one such variable is included in each file.

Example:

```yml
input_root:
   path: /ldata2/MOSAIK/Berlin_static_driver_data
   file_x_UTM: Berlin_CoordinatesUTM_x_15m_DLR.nc
   file_y_UTM: Berlin_CoordinatesUTM_y_15m_DLR.nc
   file_lat: Berlin_CoordinatesLatLon_y_15m_DLR.nc
   file_lon: Berlin_CoordinatesLatLon_x_15m_DLR.nc
   file_zt: Berlin_terrain_height_15m_DLR.nc
   file_buildings_2d: Berlin_building_height_15m_DLR.nc
   file_building_id: Berlin_building_id_15m_DLR.nc
   file_building_type: Berlin_building_type_15m_DLR.nc
   file_bridges_2d: Berlin_bridges_height_15m_DLR.nc
   file_bridges_id: Berlin_bridges_id_15m_DLR.nc
   file_lai:  Berlin_leaf_area_index_15m_DLR_WANG_summer.nc
   file_vegetation_type: Berlin_vegetation_type_15m_DLR.nc
   file_vegetation_height: Berlin_vegetation_patch_height_15m_DLR.nc
   file_pavement_type: Berlin_pavement_type_15m_DLR.nc
   file_water_type: Berlin_water_type_15m_DLR.nc
   file_soil_type: Berlin_soil_type_15m_DLR.nc
   file_street_type: Berlin_street_type_15m_DLR.nc
   file_street_crossings: Berlin_street_crossings_15m_DLR.nc
   file_tree_height: Berlin_trees_height_clean_15m.nc
   file_tree_crown_diameter: Berlin_tree_crown_15m_DLR.nc
   file_tree_trunk_diameter: Berlin_trees_trunk_clean_15m.nc
   file_tree_type: Berlin_trees_type_15m_DLR.nc
   file_patch_height: Berlin_vegetation_patch_height_15m_DLR.nc
   file_vegetation_on_roofs: Berlin_vegetation_on_roofs_15m_DLR.nc
```

### `domain_ABC`, ..., `domain_XYZ` sections

This section contains settings for each model domain for the PALM run. If the name is omitted, the name `root` is assumed. In case of a nested run, the sections for the non-root domains must be named individually, e.g. `domain_N01`, `domain_N02`, etc. as it is done in the PALM parameter file.

The corresponding input data set for each domain can be defined in the `input` parameter. If not set, the name of both `input` and `domain` sections is used to find the matching input data set (e.g. `input_root` for `domain_root`). If there is only one `input` section, this is used.

If geographical coordinates of the output should be calculated, i.e. if they are not supplied in the input data with `file_x_UTM`, `file_y_UTM` etc., it is sufficient to either set `origin_x`/`origin_y` or `origin_lon`/`origin_lat`. Note that also `epsg` must be set in the `settings` section.

| Variable | Data type | Default value | Description |
|----------|-----------|---------------|-------------|
| `pixel_size`\*               | float   |  | Size (in m) of a single pixel in x/y direction (equal to grid spacing in x and y). |
| `input_lower_left_x`              | float   | `None` | Distance (in m) along x-direction between the lower-left corner of the model domain and the lower-left corner of the input data. This parameter is used to shift the model domain with respect to the provided input data. Only needed for netCDF input data. |
| `input_lower_left_y`              | float   | `None` | Distance (in m) along y-direction between the lower-left corner of the model domain and the lower-left corner of the input data. This parameter is used to shift the model domain with respect to the provided input data. Only needed for netCDF input data. |
| `lower_left_x`               | float   | `None` | Only for nested domains: Distance (in m) along x-direction between the lower-left corner of the nested domain and the lower-left corner of the root parent domain. This parameter is used to define the coordinates of origin of the nested domain. This parameter is not required if the origin is defined via `origin_x`/`origin_y` or `origin_lon`/`origin_lat`. |
| `lower_left_y`               | float   | `None` | Only for nested domains: Distance (in m) along y-direction between the lower-left corner of the nested domain and the lower-left corner of the root parent domain. This parameter is used to define the coordinates of origin of the nested domain. This parameter is not required if the origin is defined via `origin_x`/`origin_y` or `origin_lon`/`origin_lat`. |
| `origin_x`                   | float   | `None` | x-coordinate of the left border of the lower-left grid point of the PALM domain in the CRS defined by `epsg` in the `settings` section. |
| `origin_y`                   | float   | `None` | y-coordinate of the lower border of the lower-left grid point of the PALM domain in the CRS defined by `epsg` in the `settings` section. |
| `origin_lon`                 | float   | `None` | Longitude of the left border of the lower-left grid point of the PALM domain in WGS84. |
| `origin_lat`                 | float   | `None` | Latitude of the lower border of the lower-left grid point of the PALM domain in WGS84. |
| `nx`\*                       | integer |  | Number of grid points in x-direction. It equals the `nx` setting in the PALM parameter file so the actual number of grid points is `nx+1`. |
| `ny`\*                       | integer |  | Number of grid points in y-direction. It equals the `ny` setting in the PALM parameter file so the actual number of grid points is `ny+1`.  |
| `dz`\*                       | float   |  | Vertical grid spacing in PALM (m). This parameter is needed when `buildings_3d`, `street_trees`, `canopy_patches`, `interpolate_terrain`, or `use_palm_z_axis` is used. |
| `input`                      | string  | `None` | Name of the `input` section to be used for this domain. This parameter is used to match the input data set with the domain. If not set, the name of both `input` and `domain` sections or the `pixel_size` parameter (deprecated) is used to find the matching input data set. If there is only one `input` section, this is used. |
| `bridge_depth`               | float   |  `3.0`  | Vertical depth or thickness (m) of all bridge elements in the domain. Bridges are treated as building grid cells in `buildings_3d`. The values in `file_bridges_2d` define for each pixel the maximum height of these grid cells and `bridge_depth` defines how far these grid cells extend downwards.  |
| `buildings_3d`               | logical | `False` | Use 3D buildings via the `buildings_3d` array instead of `buildings_2d`. If bridges are present in the simulation domain, `buildings_3d` is generated in any case. |
| `building_albedo_type`         | dictionary | `None` | Albedo type for all buildings. Possible settings are `wall_gfl`, `wall_agfl`, `wall_roof`, `window_gfl`, `window_agfl`, `window_roof`, `green_gfl`, `green_agfl`, `green_roof`. See below for details. |
| `building_emissivity`          | dictionary | `None` | Emissivity for all buildings. Possible settings are `wall_gfl`, `wall_agfl`, `wall_roof`, `window_gfl`, `window_agfl`, `window_roof`, `green_gfl`, `green_agfl`, `green_roof`. See below for details. |
| `building_fraction`            | dictionary | `None` | Building surface fractions for all buildings. Possible settings are `wall_gfl`, `wall_agfl`, `wall_roof`, `window_gfl`, `window_agfl`, `window_roof`, `green_gfl`, `green_agfl`, `green_roof`. See below for details. |
| `building_general_pars`        | dictionary | `None` | General parameters for all buildings. Possible settings are `height_gfl`, `green_type_roof`. See below for details. |
| `building_heat_capacity`       | dictionary | `None` | Heat capacity of the urban surface layers for all buildings. Possible settings are `wall_gfl`, `wall_agfl`, `wall_roof`, `window_gfl`, `window_agfl`, `window_roof`, `green_gfl`, `green_agfl`, `green_roof`. See below for details. |
| `building_heat_conductivity`   | dictionary | `None` | Heat conductivity of the urban surface layers for all buildings. Possible settings are `wall_gfl`, `wall_agfl`, `wall_roof`, `window_gfl`, `window_agfl`, `window_roof`, `green_gfl`, `green_agfl`, `green_roof`. See below for details. |
| `building_indoor_pars`         | dictionary | `None` | Indoor parameters for all buildings. Possible settings are `indoor_temperature_summer`, `indoor_temperature_winter`, `shading_window`, `g_window`, `u_window`, `airflow_unoccupied`, `airflow_occupied`, `heat_recovery_efficiency`, `effective_surface`, `inner_heat_storage`, `ratio_surface_floor`, `heating_capacity_max`, `cooling_capacity_max`, `heat_gain_high`, `heat_gain_low`, `height_storey`, `height_ceiling_construction`, `heating_factor`, `cooling_factor`. See below for details. |
| `building_lai`                 | dictionary | `None` | LAI of urban surfaces for all buildings. Possible settings are `gfl`, `agfl`, `roof`. See below for details. |
| `building_roughness_length`    | dictionary | `None` | Roughness length of the urban surfaces for all buildings. Possible settings are `gfl`, `agfl`, `roof`. See below for details. |
| `building_roughness_length_qh` | dictionary | `None` | Roughness length for heat and moisture for all urban surfaces. Possible settings are `gfl`, `agfl`, `roof`. See below for details. |
| `building_thickness`           | dictionary | `None` | Layer thickness of the urban surfaces for all buildings. Possible settings are `wall_gfl`, `wall_agfl`, `wall_roof`, `window_gfl`, `window_agfl`, `window_roof`, `green_gfl`, `green_agfl`, `green_roof`. See below for details. |
| `building_transmissivity`      | dictionary | `None` | Window transmissivity of the urban surface for all buildings. Possible settings are `gfl`, `agfl`, `roof`. See below for details. |
| `allow_high_vegetation`      | logical | `False` | If set to `True`, it is allowed to have unresolved high vegetation classes according in the `vegetation_type` distribution. Note that this can involve very large roughness lengths > 0.5 m. If the vertical grid spacing is close to or smaller than this threshold the PALM run will crash and/or does not provide meaningful results. It is generally recommended to set this parameter to `False` whenever the grid spacing is small enough to resolve canopy patches by 2 or more vertical grid levels. If set to `False` pixels where a high vegetation type was prescribed will be converted into a 3D leaf area density canopy using the canopy generator. |
| `generate_vegetation_patches`| logical | `True` | If set to `True`, the embedded canopy generator will convert all surface pixels that contain high vegetation into a 3D leaf area density distribution. This applies to pixels where `vegetation_type` is set to a high vegetation type, or where the vegetation height field suggests high vegetation. Note that only pixels with heights `> 2*dz` are converted, while all other pixels will be parameterized via the `vegetation_type` field. |
| `use_palm_z_axis`            | logical | `False` | If set to `True`, the static driver will raster the input data on the z-grid of PALM for output. Note that PALM will convert continuous static driver data itself on its grid and apply additional filtering procedures. It is thus recommended to set this parameter to `False` unless `interpolate_terrain: True` in nested set-ups.|
| `interpolate_terrain`        | logical | `False` | If set to `True`, the terrain height is interpolated and blended over between parent and child domains in order to avoid severe steps in terrain height due to different grid spacings between parent and child. |
| `domain_parent`              | string  | `None` | Name of the parent domain of the current domain. If the current domain is the root domain, do not set this parameter. |
| `vegetation_on_roofs`        | logical | `True` | If set to `True`, allow green roofs. |
| `street_trees`               | logical | `True` | If set to `True`, information on individual street trees will be used to generate a 3D leaf area density and basal area density distribution for each tree. In contrast to vegetation patches, where a closed canopy is assumed and information is only distributed vertically for each pixel, street trees have a 3D shape that is mapped on the simulation domain. |
| `overhanging_trees`          | logical | `True` | If set to `False`, no LAD volumes of trees are generated above surfaces without a vegetation type. |
| `remove_low_lai_tree`        | logical | `False` | If set to `True`, all trees with an LAI < `lai_tree_lower_threshold` are removed from the dataset. If set to `False`, those trees are considered with LAI = `lai_tree_lower_threshold`. |
| `water_temperature` | float or dictionary | `None` | Water temperature in K for one or several water types as indicated by their name or their index 0 to 5. Also allows one value, which is applied to all water types. |

(*) This parameter is mandatory

Example:

```yml
domain_root:
   pixel_size: 15.0
   origin_x: 19605
   origin_y: 20895
   nx: 199
   ny: 199
   bridge_depth: 3.0
   buildings_3d: False
   dz: 15.0
   allow_high_vegetation: True
   generate_vegetation_patches: True
   use_palm_z_axis: False
   interpolate_terrain: False
   vegetation_on_roofs: True
   street_trees: True
   water_temperature: 
      lake: 285
      fountain: 290
   building_albedo_type:
     window: 10
   building_emissivity:
     green_agfl: 0.9
   building_fraction:
     wall_gfl: 0.3
     green_gfl: 0.6
     window_gfl: 0.1
   building_general_pars:
     height_gfl: 3
   building_heat_capacity:
     wall_agfl: 1520000.
     wall_roof:  709000.
   building_heat_conductivity:
     wall_agfl: 2.1
     wall_roof: 0.7
   building_indoor_pars:
     heating_capacity_max: 3
   building_lai:
     agfl: 4
   building_roughness_length:
     agfl: 0.04
   building_roughness_length_qh:
     gfl: 0.03
   building_thickness:
     green_agfl: [0.2,0.3,0.1,0.05]
   building_transmissivity:
     roof: 0.2
```

#### Setting of building parameters

This section explains how to set the building parameters for *all* buildings of a domain in the YAML configuration file. While the respective values can be set for each domain separately, building parameters of single buildings or pixels cannot be set individually at the moment. In particular, raster input for building parameters is not yet supported.

The following parameters can *optionally* be set: `building_albedo_type`, `building_emissivity`, `building_fraction`, `building_general_pars`, `building_heat_capacity`, `building_heat_conductivity`, `building_indoor_pars`, `building_lai`, `building_roughness_length`, `building_roughness_length_qh`, `building_thickness`. The input is given using a `setting: input value` structure. The possible `setting`s are listed above for each parameter; please see the respective part in the PALM documentation for details. Note that also parts of these defined strings are allowed. For example, for `building_fraction`, a `window` setting will apply to `window_gfl`, `window_agfl` and `window_roof`.

The `input value`s can be either a single value or, in the case of parameters that describe the layers of an urban surface, a list of four values. These four values represent the four urban surface layers in PALM, beginning from the outermost layer. Each layer is characterized by a value from `building_heat_capacity`, from  `building_heat_conductivity` and from `building_thickness`. If only a single value is supplied for parameters of urban surface layers, this single values will be applied to all urban surfaces. In addition, it is also possible to set the value directly to the parameter, which will be expanded to all applicable settings.

Please refer to the PALM manual for a detailled description of the parameters and their respective unit.

*Here are some examples*:

One value set for a parameter like

```yml
   building_heat_conductivity: 1.8
```

will be expanded to

```yml
   building_heat_conductivity: 
      wall_gfl: [1.8, 1.8, 1.8, 1.8]
      wall_agfl: [1.8, 1.8, 1.8, 1.8]
      wall_roof: [1.8, 1.8, 1.8, 1.8]
      window_gfl: [1.8, 1.8, 1.8, 1.8]
      window_agfl: [1.8, 1.8, 1.8, 1.8]
      window_roof: [1.8, 1.8, 1.8, 1.8]
      green_gfl: [1.8, 1.8, 1.8, 1.8]
      green_agfl: [1.8, 1.8, 1.8, 1.8]
      green_roof: [1.8, 1.8, 1.8, 1.8]
```

In the case of a parameter that sets urban surface layer properties, also a list of four values like

```yml
   building_heat_conductivity: [1.6, 1.7, 1.8, 1.9]
```

will be expanded to

```yml
   building_heat_conductivity: 
      wall_gfl: [1.6, 1.7, 1.8, 1.9]
      wall_agfl: [1.6, 1.7, 1.8, 1.9]
      wall_roof: [1.6, 1.7, 1.8, 1.9]
      window_gfl: [1.6, 1.7, 1.8, 1.9]
      window_agfl: [1.6, 1.7, 1.8, 1.9]
      window_roof: [1.6, 1.7, 1.8, 1.9]
      green_gfl: [1.6, 1.7, 1.8, 1.9]
      green_agfl: [1.6, 1.7, 1.8, 1.9]
      green_roof: [1.6, 1.7, 1.8, 1.9]
```

A partial setting like

```yml
   building_heat_conductivity:
      wall: 1.8
```

will be expanded to

```yml
   building_heat_conductivity:
      wall_gfl: [1.8, 1.8, 1.8, 1.8]
      wall_agfl: [1.8, 1.8, 1.8, 1.8]
      wall_roof: [1.8, 1.8, 1.8, 1.8]
```

## Technical documentation

### Tree database

Default values for trees if individual parameters are not provided. Default data is derived as mean values from the tree database for Berlin, Germany.

| Index | Species | Shape | Crown height/width ratio(\*) | Crown diameter (m) | Height (m) | LAI summer(\*) | LAI winter(\*) | Height of maximum LAD (m) | LAD/BAD ratio(\*) | DBH (m) |
|-------|---------|-------|------------------------------|--------------------|------------|----------------|----------------|---------------------------|-------------------|---------|
|  0 | Default|         1.0 | 1.0 | 4.0| 12.0| 3.0| 0.8| 0.6| 0.025| 0.35|
|  1 | Abies|           3.0 | 1.0 | 4.0| 12.0| 3.0| 0.8| 0.6| 0.025| 0.80|
|  2 | Acer|            1.0 | 1.0 | 7.0| 12.0| 3.0| 0.8| 0.6| 0.025| 0.80|
|  3 | Aesculus|        1.0 | 1.0 | 7.0| 12.0| 3.0| 0.8| 0.6| 0.025| 1.00|
|  4 | Ailanthus|       1.0 | 1.0 | 8.5| 13.5| 3.0| 0.8| 0.6| 0.025| 1.30|
|  5 | Alnus|           3.0 | 1.0 | 6.0| 16.0| 3.0| 0.8| 0.6| 0.025| 1.20|
|  6 | Amelanchier|     1.0 | 1.0 | 3.0|  4.0| 3.0| 0.8| 0.6| 0.025| 1.20|
|  7 | Betula|          1.0 | 1.0 | 6.0| 14.0| 3.0| 0.8| 0.6| 0.025| 0.30|
|  8 | Buxus|           1.0 | 1.0 | 4.0|  4.0| 3.0| 0.8| 0.6| 0.025| 0.90|
|  9 | Calocedrus|      3.0 | 1.0 | 5.0| 10.0| 3.0| 0.8| 0.6| 0.025| 0.50|
| 10 | Caragana|        1.0 | 1.0 | 3.5|  6.0| 3.0| 0.8| 0.6| 0.025| 0.90|
| 11 | Carpinus|        1.0 | 1.0 | 6.0| 10.0| 3.0| 0.8| 0.6| 0.025| 0.70|
| 12 | Carya|           1.0 | 1.0 | 5.0| 17.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 13 | Castanea|        1.0 | 1.0 | 4.5|  7.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 14 | Catalpa|         1.0 | 1.0 | 5.5|  6.5| 3.0| 0.8| 0.6| 0.025| 0.70|
| 15 | Cedrus|          1.0 | 1.0 | 8.0| 13.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 16 | Celtis|          1.0 | 1.0 | 6.0|  9.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 17 | Cercidiphyllum|  1.0 | 1.0 | 3.0|  6.5| 3.0| 0.8| 0.6| 0.025| 0.80|
| 18 | Cercis|          1.0 | 1.0 | 2.5|  7.5| 3.0| 0.8| 0.6| 0.025| 0.90|
| 19 | Chamaecyparis|   5.0 | 1.0 | 3.5|  9.0| 3.0| 0.8| 0.6| 0.025| 0.70|
| 20 | Cladrastis|      1.0 | 1.0 | 5.0| 10.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 21 | Cornus|          1.0 | 1.0 | 4.5|  6.5| 3.0| 0.8| 0.6| 0.025| 1.20|
| 22 | Corylus|         1.0 | 1.0 | 5.0|  9.0| 3.0| 0.8| 0.6| 0.025| 0.40|
| 23 | Cotinus|         1.0 | 1.0 | 4.0|  4.0| 3.0| 0.8| 0.6| 0.025| 0.70|
| 24 | Crataegus|       3.0 | 1.0 | 3.5|  6.0| 3.0| 0.8| 0.6| 0.025| 1.40|
| 25 | Cryptomeria|     3.0 | 1.0 | 5.0| 10.0| 3.0| 0.8| 0.6| 0.025| 0.50|
| 26 | Cupressocyparis| 3.0 | 1.0 | 3.0|  8.0| 3.0| 0.8| 0.6| 0.025| 0.40|
| 27 | Cupressus|       3.0 | 1.0 | 5.0|  7.0| 3.0| 0.8| 0.6| 0.025| 0.40|
| 28 | Cydonia|         1.0 | 1.0 | 2.0|  3.0| 3.0| 0.8| 0.6| 0.025| 0.90|
| 29 | Davidia|         1.0 | 1.0 | 10.0| 14.0| 3.0| 0.8| 0.6| 0.025| 0.40|
| 30 | Elaeagnus|       1.0 | 1.0 | 6.5|  6.0| 3.0| 0.8| 0.6| 0.025| 1.20|
| 31 | Euodia|          1.0 | 1.0 | 4.5|  6.0| 3.0| 0.8| 0.6| 0.025| 0.90|
| 32 | Euonymus|        1.0 | 1.0 | 4.5|  6.0| 3.0| 0.8| 0.6| 0.025| 0.60|
| 33 | Fagus|           1.0 | 1.0 | 10.0| 12.5| 3.0| 0.8| 0.6| 0.025| 0.50|
| 34 | Fraxinus|        1.0 | 1.0 | 5.5| 10.5| 3.0| 0.8| 0.6| 0.025| 1.60|
| 35 | Ginkgo|          3.0 | 1.0 | 4.0|  8.5| 3.0| 0.8| 0.6| 0.025| 0.80|
| 36 | Gleditsia|       1.0 | 1.0 | 6.5| 10.5| 3.0| 0.8| 0.6| 0.025| 0.60|
| 37 | Gymnocladus|     1.0 | 1.0 | 5.5| 10.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 38 | Hippophae|       1.0 | 1.0 | 9.5|  8.5| 3.0| 0.8| 0.6| 0.025| 0.80|
| 39 | Ilex|            1.0 | 1.0 | 4.0|  7.5| 3.0| 0.8| 0.6| 0.025| 0.80|
| 40 | Juglans|         1.0 | 1.0 | 7.0|  9.0| 3.0| 0.8| 0.6| 0.025| 0.50|
| 41 | Juniperus|       5.0 | 1.0 | 3.0|  7.0| 3.0| 0.8| 0.6| 0.025| 0.90|
| 42 | Koelreuteria|    1.0 | 1.0 | 3.5|  5.5| 3.0| 0.8| 0.6| 0.025| 0.50|
| 43 | Laburnum|        1.0 | 1.0 | 3.0|  6.0| 3.0| 0.8| 0.6| 0.025| 0.60|
| 44 | Larix|           3.0 | 1.0 | 7.0| 16.5| 3.0| 0.8| 0.6| 0.025| 0.60|
| 45 | Ligustrum|       1.0 | 1.0 | 3.0|  6.0| 3.0| 0.8| 0.6| 0.025| 1.10|
| 46 | Liquidambar|     3.0 | 1.0 | 3.0|  7.0| 3.0| 0.8| 0.6| 0.025| 0.30|
| 47 | Liriodendron|    3.0 | 1.0 | 4.5|  9.5| 3.0| 0.8| 0.6| 0.025| 0.50|
| 48 | Lonicera|        1.0 | 1.0 | 7.0|  9.0| 3.0| 0.8| 0.6| 0.025| 0.70|
| 49 | Magnolia|        1.0 | 1.0 | 3.0|  5.0| 3.0| 0.8| 0.6| 0.025| 0.60|
| 50 | Malus|           1.0 | 1.0 | 4.5|  5.0| 3.0| 0.8| 0.6| 0.025| 0.30|
| 51 | Metasequoia|     5.0 | 1.0 | 4.5| 12.0| 3.0| 0.8| 0.6| 0.025| 0.50|
| 52 | Morus|           1.0 | 1.0 | 7.5| 11.5| 3.0| 0.8| 0.6| 0.025| 1.00|
| 53 | Ostrya|          1.0 | 1.0 | 2.0|  6.0| 3.0| 0.8| 0.6| 0.025| 1.00|
| 54 | Parrotia|        1.0 | 1.0 | 7.0|  7.0| 3.0| 0.8| 0.6| 0.025| 0.30|
| 55 | Paulownia|       1.0 | 1.0 | 4.0|  8.0| 3.0| 0.8| 0.6| 0.025| 0.40|
| 56 | Phellodendron|   1.0 | 1.0 | 13.5| 13.5| 3.0| 0.8| 0.6| 0.025| 0.50|
| 57 | Picea|           3.0 | 1.0 | 3.0| 13.0| 3.0| 0.8| 0.6| 0.025| 0.90|
| 58 | Pinus|           3.0 | 1.0 | 6.0| 16.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 59 | Platanus|        1.0 | 1.0 | 10.0| 14.5| 3.0| 0.8| 0.6| 0.025| 1.10|
| 60 | Populus|         1.0 | 1.0 | 9.0| 20.0| 3.0| 0.8| 0.6| 0.025| 1.40|
| 61 | Prunus|          1.0 | 1.0 | 5.0|  7.0| 3.0| 0.8| 0.6| 0.025| 1.60|
| 62 | Pseudotsuga|     3.0 | 1.0 | 6.0| 17.5| 3.0| 0.8| 0.6| 0.025| 0.70|
| 63 | Ptelea|          1.0 | 1.0 | 5.0|  4.0| 3.0| 0.8| 0.6| 0.025| 1.10|
| 64 | Pterocaria|      1.0 | 1.0 | 10.0| 12.0| 3.0| 0.8| 0.6| 0.025| 0.50|
| 65 | Pterocarya|      1.0 | 1.0 | 11.5| 14.5| 3.0| 0.8| 0.6| 0.025| 1.60|
| 66 | Pyrus|           3.0 | 1.0 | 3.0|  6.0| 3.0| 0.8| 0.6| 0.025| 1.80|
| 67 | Quercus|         1.0 | 1.0 | 8.0| 14.0| 3.1| 0.1| 0.6| 0.025| 0.40|
| 68 | Rhamnus|         1.0 | 1.0 | 4.5|  4.5| 3.0| 0.8| 0.6| 0.025| 1.30|
| 69 | Rhus|            1.0 | 1.0 | 7.0|  5.5| 3.0| 0.8| 0.6| 0.025| 0.50|
| 70 | Robinia|         1.0 | 1.0 | 4.5| 13.5| 3.0| 0.8| 0.6| 0.025| 0.50|
| 71 | Salix|           1.0 | 1.0 | 7.0| 14.0| 3.0| 0.8| 0.6| 0.025| 1.10|
| 72 | Sambucus|        1.0 | 1.0 | 8.0|  6.0| 3.0| 0.8| 0.6| 0.025| 1.40|
| 73 | Sasa|            1.0 | 1.0 | 10.0| 25.0| 3.0| 0.8| 0.6| 0.025| 0.60|
| 74 | Sequoiadendron|  5.0 | 1.0 | 5.5| 10.5| 3.0| 0.8| 0.6| 0.025| 1.60|
| 75 | Sophora|         1.0 | 1.0 | 7.5| 10.0| 3.0| 0.8| 0.6| 0.025| 1.40|
| 76 | Sorbus|          1.0 | 1.0 | 4.0|  7.0| 3.0| 0.8| 0.6| 0.025| 1.10|
| 77 | Syringa|         1.0 | 1.0 | 4.5|  5.0| 3.0| 0.8| 0.6| 0.025| 0.60|
| 78 | Tamarix|         1.0 | 1.0 | 6.0|  7.0| 3.0| 0.8| 0.6| 0.025| 0.50|
| 79 | Taxodium|        5.0 | 1.0 | 6.0| 16.5| 3.0| 0.8| 0.6| 0.025| 0.60|
| 80 | Taxus|           2.0 | 1.0 | 5.0|  7.5| 3.0| 0.8| 0.6| 0.025| 1.50|
| 81 | Thuja|           3.0 | 1.0 | 3.5|  9.0| 3.0| 0.8| 0.6| 0.025| 0.70|
| 82 | Tilia|           3.0 | 1.0 | 7.0| 12.5| 3.0| 0.8| 0.6| 0.025| 0.70|
| 83 | Tsuga|           3.0 | 1.0 | 6.0| 10.5| 3.0| 0.8| 0.6| 0.025| 1.10|
| 84 | Ulmus|           1.0 | 1.0 | 7.5| 14.0| 3.0| 0.8| 0.6| 0.025| 0.80|
| 85 | Zelkova|         1.0 | 1.0 | 4.0|  5.5| 3.0| 0.8| 0.6| 0.025| 1.20|
| 86 | Zenobia|         1.0 | 1.0 | 5.0|  5.0| 3.0| 0.8| 0.6| 0.025| 0.40|

(*) Preliminary parameter.

## Best practices

The following example is a best practice setting for a high-resolution (e.g. 1 m grid spacing) non-nested run in which most of the vegetation can be resolved via a 3D leaf area density distribution:

```yml
domain_root:
   pixel_size: 1.0
   origin_x: ...
   origin_y: ...
   nx: ...
   ny: ...
   dz: ...
   allow_high_vegetation: False
   buildings_3d: True
   generate_vegetation_patches: True
   use_palm_z_axis: False
   interpolate_terrain: False
   vegetation_on_roofs: True
   street_trees: True
```

For a nested run, the following settings should work nicely to avoid terrain height issues:

```yml
domain_root:
   pixel_size: 15.0
   dz: 15.0
   origin_x: ...
   origin_y: ...
   nx: ...
   ny: ...
   buildings_3d: False
   allow_high_vegetation: True
   generate_vegetation_patches: True
   use_palm_z_axis: False
   interpolate_terrain: False
   vegetation_on_roofs: False
   street_trees: True

domain_N02:
   domain_parent: root
   pixel_size: 1.0
   dz: 1.0
   origin_x: ...
   origin_y: ...
   nx: ...
   ny: ...
   buildings_3d: True
   allow_high_vegetation: False
   generate_vegetation_patches: True
   use_palm_z_axis: True
   interpolate_terrain: True
   vegetation_on_roofs: True
   street_trees: True
```

## Literature

* Heldens, W., Burmeister, C., Kanani-Sühring, F., Maronga, B., Pavlik, D., Sühring, M., Zeidler, J. and Esch, T. (2020): Geospatial input data for the PALM model system 6.0: model requirements, data sources and processing, Geosci. Model Dev., 13, 5833–5873, [doi: 10.5194/gmd-13-5833-2020](https://doi.org/10.5194/gmd-13-5833-2020).
* Lalic, B. and Mihailovic, D. T. (2004): An Empirical Relation Describing Leaf-Area Density inside the Forest for Environmental Modeling. Journal of Applied Meteorology, vol. 43, no. 4, pp. 641–645, [doi: 10.1175/1520-0450(2004)043<0641:AERDLD>2.0.CO;2](https://doi.org/10.1175/1520-0450(2004)043<0641:AERDLD>2.0.CO;2).
* Markkanen, T., Rannik, Ü., Marcolla, B., Cescatti, A. and Vesala, T. (2003): Footprints and Fetches for Fluxes over Forest Canopies with Varying Structure and Density. Boundary-Layer Meteorology 106, 437–459, [doi: 10.1023/A:1021261606719](https://doi.org/10.1023/A:1021261606719).
