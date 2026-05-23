# Overview of dimensions, parameters and input data

An overview of the names used by `palm_csd` for dimensions and other parameters, mainly for reference of the attributes

---

## Building parameter dimensions

### `building_surface_level`

| Attribute         | Value | Description              |
|-------------------|-------|--------------------------|
| `gfl`             | 0     | ground floor level       |
| `agfl`            | 1     | above ground floor level |
| `roof`            | 2     | roof                     |

### `building_surface_type`

| Attribute         | Value | Description                             |
|-------------------|-------|-----------------------------------------|
| `wall_gfl`        | 0     | wall ground floor level                 |
| `wall_agfl`       | 1     | wall above ground floor level           |
| `wall_roof`       | 2     | wall roof                               |
| `window_gfl`      | 3     | window ground floor level               |
| `window_agfl`     | 4     | window above ground floor level         |
| `window_roof`     | 5     | window roof                             |
| `green_gfl`       | 6     | green on wall ground floor level        |
| `green_agfl`      | 7     | green on wall above ground floor level  |
| `green_roof`      | 8     | green on roof ground floor level        |

## Building parameters

### `building_general_pars`

| Attribute         | Value | Description                  |
|-------------------|-------|------------------------------|
| `height_gfl`      | 0     | ground floor level height    |
| `green_type_roof` | 1     | type of green roof           |

### `building_indoor_pars`

| Attribute                     | Value | Description                                           |
|-------------------------------|-------|-------------------------------------------------------|
| `indoor_temperature_summer`   | 0     | indoor target summer temperature                      |
| `indoor_temperature_winter`   | 1     | indoor target winter temperature                      |
| `shading_window`              | 2     | shading factor                                        |
| `g_window`                    | 3     | g-value windows                                       |
| `u_window`                    | 4     | u-value windows                                       |
| `airflow_unoccupied`          | 5     | basic airflow without occupancy of the room           |
| `airflow_occupied`            | 6     | additional airflow dependent on occupancy of the room |
| `heat_recovery_efficiency`    | 7     | heat recovery efficiency                              |
| `effective_surface`           | 8     | dynamic parameter specific effective surface          |
| `inner_heat_storage`          | 9     | dynamic parameter inner heat storage                  |
| `ratio_surface_floor`         | 10    | ratio internal surface/floor area                     |
| `heating_capacity_max`        | 11    | maximal heating capacity                              |
| `cooling_capacity_max`        | 12    | maximal cooling capacity                              |
| `heat_gain_high`              | 13    | additional internal heat gains dependent on occupancy |
| `heat_gain_low`               | 14    | basic internal heat gains without occupancy           |
| `height_storey`               | 15    | storey height                                         |
| `height_ceiling_construction` | 16    | ceiling construction height                           |
| `heating_factor`              | 17    | anthropogenic heat output factor for heating          |
| `cooling_factor`              | 18    | anthropogenic heat output factor for cooling          |

## Building and surface types, and their parameters

### `building_type`

| Attribute               | Value | Description                |
|-------------------------|-------|----------------------------|
| `residential_1950`      | 1     | Residential, before 1950   |
| `residential_1951_2000` | 2     | Residential, 1951 -- 2000  |
| `residential_2001`      | 3     | Residential, after 2001    |
| `office_1950`           | 4     | Office, before 1950        |
| `office_1951_2000`      | 5     | Office, 1951 -- 2000       |
| `office_2001`           | 6     | Office, after 2001         |
| `bridges`               | 7     | Bridges                    |

### `pavement_type`

| Attribute               | Value | Description          |
|-------------------------|-------|----------------------|
| `asphalt_concrete_mix`  | 1     | asphalt concrete mix |
| `asphalt`               | 2     | asphalt              |
| `concrete`              | 3     | concrete             |
| `sett`                  | 4     | sett                 |
| `paving_stones`         | 5     | paving stones        |
| `cobblestone`           | 6     | cobblestone          |
| `metal`                 | 7     | metal                |
| `wood`                  | 8     | wood                 |
| `gravel`                | 9     | gravel               |
| `fine_gravel`           | 10    | fine gravel          |
| `pebblestone`           | 11    | pebblestone          |
| `woodchips`             | 12    | woodchips            |
| `tartan`                | 13    | tartan               |
| `artificial_turf`       | 14    | artificial turf      |
| `clay`                  | 15    | clay                 |

### `soil_type`

| Attribute     | Value | Description |
|---------------|-------|-------------|
| `coarse`      | 1     | Coarse      |
| `medium`      | 2     | Medium      |
| `medium_fine` | 3     | Medium-fine |
| `fine`        | 4     | Fine        |
| `very_fine`   | 5     | Very fine   |
| `organic`     | 6     | Organic     |

### `street_type`

| Attribute              | Value | Description                   |
|------------------------|-------|-------------------------------|
| `unclassified`         | 1     | Unclassified                  |
| `cycleway`             | 2     | Cycleway                      |
| `footway_pedestrian`   | 3     | Footway or pedestrian area    |
| `path`                 | 4     | Path                          |
| `track`                | 5     | Track                         |
| `living_street`        | 6     | Living street                 |
| `service`              | 7     | Service                       |
| `residential`          | 8     | Residential                   |
| `tertiary`             | 9     | Tertiary                      |
| `tertiary_link`        | 10    | Tertiary link                 |
| `secondary`            | 11    | Secondary                     |
| `secondary_link`       | 12    | Secondary link                |
| `primary`              | 13    | Primary                       |
| `primary_link`         | 14    | Primary link                  |
| `trunk`                | 15    | Trunk                         |
| `trunk_link`           | 16    | Trunk link                    |
| `motorway`             | 17    | Motorway                      |
| `motorway_link`        | 18    | Motorway link                 |
| `raceway`              | 19    | Raceway                       |

### `vegetation_pars`

| Attribute                        | Value | Description                                                                 |
|----------------------------------|-------|-------------------------------------------------------------------------|
| `canopy_resistance_min`          | 0     | Minimum canopy resistance                                               |
| `lai`                            | 1     | Leaf area index                                                         |
| `vegetation_coverage`            | 2     | Vegetation coverage                                                     |
| `canopy_resistance_coefficient`  | 3     | Canopy resistance coefficient                                           |
| `roughness_length`               | 4     | Roughness length for momentum                                           |
| `roughness_length_qh`            | 5     | Roughness length for heat and moisture                                  |
| `heat_conductivity_stable`       | 6     | Skin layer heat conductivity (stable conditions)                        |
| `heat_conductivity_unstable`     | 7     | Skin layer heat conductivity (unstable conditions)                      |
| `fraction_shortwave_soil`        | 8     | Fraction of incoming shortwave radiation transmitted directly to the soil (not implemented yet) |
| `heat_capacity`                  | 9     | Heat capacity of the surface                                            |
| `albedo_type`                    | 10    | Albedo type                                                             |
| `emissivity`                     | 11    | Surface emissivity                                                      |

### `vegetation_type`

| Attribute                     | Value | Description                    |
|-------------------------------|-------|--------------------------------|
| `user_defined`                | 0     | user defined                   |
| `bare_soil`                   | 1     | bare soil                      |
| `crops_mixed_farming`         | 2     | crops, mixed farming           |
| `short_grass`                 | 3     | short grass                    |
| `evergreen_needleleaf_trees`  | 4     | evergreen needleleaf trees     |
| `deciduous_needleleaf_trees`  | 5     | deciduous needleleaf trees     |
| `evergreen_broadleaf_trees`   | 6     | evergreen broadleaf trees      |
| `deciduous_broadleaf_trees`   | 7     | deciduous broadleaf trees      |
| `tall_grass`                  | 8     | tall grass                     |
| `desert`                      | 9     | desert                         |
| `tundra`                      | 10    | tundra                         |
| `irrigated_crops`             | 11    | irrigated crops                |
| `semidesert`                  | 12    | semidesert                     |
| `ice_caps_glaciers`           | 13    | ice caps and glaciers          |
| `bogs_marshes`                | 14    | bogs and marshes               |
| `evergreen_shrubs`            | 15    | evergreen shrubs               |
| `deciduous_shrubs`            | 16    | deciduous shrubs               |
| `mixed_forest_woodland`       | 17    | mixed forest/woodland          |
| `interrupted_forest`          | 18    | interrupted forest             |

#### High vegetation

The following vegetation types are considered as high (grown) vegetation with large roughness lengths.

- `evergreen_needleleaf_trees`
- `deciduous_needleleaf_trees`
- `evergreen_broadleaf_trees`
- `deciduous_broadleaf_trees`
- `mixed_forest_woodland`
- `interrupted_forest`

### `water_pars`

| Attribute                    | Value | Description                                                    |
|------------------------------|-------|----------------------------------------------------------------|
| `water_temperature`          | 0     | water temperature                                              |
| `roughness_length`           | 1     | roughness length for momentum                                  |
| `roughness_length_qh`        | 2     | roughness length for heat                                      |
| `heat_conductivity_stable`   | 3     | heat conductivity between skin layer and water (stable conditions)   |
| `heat_conductivity_unstable` | 4     | heat conductivity between skin layer and water (unstable conditions) |
| `albedo_type`                | 5     | albedo type                                                    |
| `emissivity`                 | 6     | surface emissivity                                             |

### `water_type`

| Attribute   | Value | Description |
|-------------|-------|-------------|
| `lake`      | 1     | Lake        |
| `river`     | 2     | River       |
| `ocean`     | 3     | Ocean       |
| `pond`      | 4     | Pond        |
| `fountain`  | 5     | Fountain    |

## Input data

| Attribute             | Description                            |
|-----------------------|----------------------------------------|
| `bridges_2d`          | Height of the bridge                   |
| `bridges_id`          | Bridge id                              |
| `building_id`         | Building id                            |
| `building_type`       | Building type                          |
| `buildings_2d`        | Height of the building                 |
| `lat`                 | Latitude                               |
| `lai`                 | Leaf area index                        |
| `lcz`                 | Local climate zone                     |
| `lon`                 | Longitude                              |
| `patch_type`          | Patch type                             |
| `pavement_type`       | Type of the pavement                   |
| `soil_type`           | Soil type                              |
| `street_crossings`    | Street crossings                       |
| `street_type`         | Street type                            |
| `tree_crown_diameter` | Crown diameter of the tree             |
| `tree_height`         | Height of the tree                     |
| `tree_lai`            | Leaf area index of the tree            |
| `tree_shape`          | Shape of the tree                      |
| `tree_trunk_diameter` | Trunk diameter of the tree             |
| `tree_type`           | Type of the tree                       |
| `tree_type_name`      | Type name of the tree                  |
| `vegetation_height`   | Height of the vegetation               |
| `vegetation_on_roofs` | Vegetation on roofs                    |
| `vegetation_type`     | Vegetation type                        |
| `water_temperature`   | Water temperature                      |
| `water_type`          | Water type                             |
| `x_utm`               | X-coordinate in UTM                    |
| `y_utm`               | Y-coordinate in UTM                    |
| `zt`                  | Height of the surface above the ground |
