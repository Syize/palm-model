In the following, the defined structure of the NetCDF input file `[run_id]_uv` will be specified.
An example input file can be found under *.../model/tests/cases/urban_environment/INPUT/urban_environment_uv* .

### Dimensions

Following dimensions are defined in the UV input file. For the LOD1 approach, `sun_zenith` and `wavelength` are mandatory, while for the LOD2 approach `sun_zenith`, `wavelength`, `view_azimuth` and `view_zenith` are mandatory. The dimension size of `view_azimuth` and `view_zenith` must correspond to the discrete number of spherical angles along azimuth and zenith in the RTM (see [raytrace_discrete_azims](/Reference/LES_Model/Namelists/#radiation_parameters--raytrace_discrete_azims) and [raytrace_discrete_elevs](/Reference/LES_Model/Namelists/#radiation_parameters--raytrace_discrete_elevs) in namelist [radiation_parameters](/Reference/LES_Model/Namelists/#radiation-parameters)).

| Dimension name | Value      | Description                                      |
|----------------|------------|--------------------------------------------------|
| `sun_zenith`   | at least 2 | number of sun-zenith angles. Recommended: >=90   |
| `wavelength`   | 64         | number of wavelength. Recommended: >= 120        |
| `view_azimuth` | at least 1 | number of view azimuth angles. Recommended: >=60 |
| `view_zenith`  | at least 1 | number of view zenith angles. Recommended: >=90  |


### Variables

Following variables are defined in the UV input file. For the LOD1 approach, `uv_dir_irradiance` and `uv_diff_irradiance` are mandatory, while for the LOD2 approach `uv_radiance` is mandatory. All variables are defined as `NC_FLOAT`.

| Variable name        | Unit             | Dimension(s)                                                 | Description                 |
|----------------------|------------------|--------------------------------------------------------------|-----------------------------|
| `uv_dir_irradiance`  | mW m-2 nm-1      | (`sun_zenith`, `wavelength`)                                 | direct spectral irradiance  |
| `uv_diff_irradiance` | mW m-2 nm-1      | (`sun_zenith`, `wavelength`)                                 | diffuse spectral irradiance |
| `uv_radiance`        | mW m-2 nm-1 sr-1 | (`sun_zenith`, `wavelength`, `view_zenith`, `view_azimuth` ) | spectral radiance           |
