---
title: Overview
---
# UV Radiation Model

---

!!! warning
    This site is  Work in Progress.

## Overview

The UV-radiation model enables the simulation of (erythemally-weighted) UV-irradiances near the surface, in order to quantify the UV-exposure at a given position in a complex environment and to evaluate possible mitigation strategies, e.g. by additional sun protection such as awnings. Therefore, shading of the direct and diffuse portion of UV-radiation by buildings and trees, transmission of UV-radiation by trees, as well as multiple Lambertian reflections are considered based on PALMs radiative transfer model (RTM, [Krc et al. (2021)](https://gmd.copernicus.org/articles/14/3095/2021/)). To model UV-irradiances near the surface, the UV-radiation model takes the undisturbed externally provided and spectrally resolved UV-radiation at the top of the urban layer (the uppermost rooftop level in the model domain) as boundary condition and models the impact of buildings, trees, 3D obstacles and reflections. The external UV-radiation is prescribed by a NetCDF file named `[run_id]_uv` in the respective `INPUT` directory of the simulation run, where `run_id` is the run-identifier of the simulation. The PALM model system offers the pre-processor tool `uv2palm` to create that file based on the libRadtran tool `uvspec` ([Mayer and Kylling, 2005](https://acp.copernicus.org/articles/5/1855/2005/), [Emde et al. 2016](https://gmd.copernicus.org/articles/9/1647/2016/)).

Two modelling approaches have been implemented to model the UV-irradiation:
i) An LOD1 approach where the shading of direct radiation is considered via a geometric approach, while the incoming diffuse radiation is reduced by an effective sky-view factor according to [Krc et al. (2021)](https://gmd.copernicus.org/articles/14/3095/2021/). This approach is based on an isotropic distribution of the diffuse radiation portion over the sky and requires input of externally provided, spectrally resolved portions of diffuse and direct irradiation in the UV spectral range (280 - 400 nm) at the top of the urban layer.

ii) An LOD2 approach, where the shading of direct and diffuse radiation is modelled based on a spherical-angle dependent method, taking into account a possibly anisotropic distribution of the diffuse radiation portion. This approach requires input of externally provided, spectrally resolved radiances in the UV spectral range (280 - 400 nm) at the top of the urban layer.

Since the incoming UV-radiation strongly depends on the sun-zenith angle. The actual sun position in the model is determined and the corresponding data is read from the UV input file. If the model sun-zenith angle is in between two sun-zenith angles in the external data set, the external UV-radiation data is linearly interpolated in between.

The UV-radiation model assumes UV-specific values of the albedo for the multiple Lambertian reflections. UV-specific albedo values are taken from [Turner and Parisi (2018)](https://www.mdpi.com/1660-4601/15/7/1507).


## Structure of NetCDF Input File

See [the format description in the reference section](/Reference/LES_Model/Iofiles/PIDS_UV).

## Creation of NetCDF Input File

The pre-processor `uv2palm` is developed to create external UV-radiation scenarios based on the libRadtran package `uvspec`. `uv2palm` requires a libRadtran [installation](http://www.libradtran.org/doku.php?id=download). `uv2palm` is designed to be user friendly, so that there is no need to setup and carry-out own `uvspec` simulations. Instead, `uv2palm` assumes a pre-configured `uvspec` input file in the background and the user can only specify single parameters on top of this. `uv2palm` and the considered UV-scenario can be configured by a `.yml` file.
In the following, an example of an `uv2palm` input file is given:

```
data_files_path: /home/libRadtran-2.0.4/data

atmosphere_file: /home/libRadtran-2.0.4/data/atmmod/afglms.dat
source_solar: /home/libRadtran-2.0.4/data/solar_flux/atlas_plus_modtran

albedo: 0.05
altitude: 0.5
mol_modify: O3 350. DU
pressure: 1023
day_of_year: 139
sza:                  # zenith angle of the sun
  lower_limit: 0
  upper_limit: 90
  step_size: 1        # 1 degree steps in SZA
phi0: 0.0             # sun is in the South

umu:
  lower_limit: 0
  upper_limit: 90
  step_size: 1         # 1 degree steps in view zenith

phi:
  lower_limit: 0
  upper_limit: 360
  step_size: 6         # 6 degree steps in view azimuth

wavelength:
  lower_limit: 280
  upper_limit: 400
  step_size: 1.0       # 1 nm stepsize

# customized output of lambda, spectral irradiance direct and diffuse,
# global irradiance, as well as spectral radiances
output_user: [lambda, edir, edn, eglo, uu]

aerosol_default: true
aerosol_haze: 1          # rural type aerosols in the lowest 2km
aerosol_visibility: 100  # visibility in km
```

### Installation of uv2palm

A description of `uv2palm` installation will follow soon.

### Usage of uv2palm

`uv2palm` can be used in the following way:

`python3 -m uv2palm --uv-setup tests/example_configurations/example-01_fast.yml --output-driver ./tmp/driver.nc`


## Usage

The UV-radiation model requires the building-surface model (see namelist [urban_surface_parameters](/Reference/LES_Model/Namelists/#urban-surface-parameters)) and the land-surface model (see namelist [land_surface_parameters](/Reference/LES_Model/Namelists/#land-surface-parameters)) switched on. Also the plant-canopy model (see namelist [plant_canopy_parameters](/Reference/LES_Model/Namelists/#plant-canopy-parameters)) and the radiation model (see namelist [radiation_parameters](/Reference/LES_Model/Namelists/#radiation-parameters)) with parameter [radiation_interactions_on](/Reference/LES_Model/Namelists/#radiation_parameters--radiation_interactions_on) = *.T.* need to be switched on.
The UV-model itself is switched on via the namelist [uv_radiation_parameters](/Reference/LES_Model/Namelists/#uv-radiation-parameters).
The UV-model can be used in two ways:

- Users that are **only** interested in UV-radiation, can use the model spinup to compute the UV-irradiation. This requires a non-zero [spinup_time](/Reference/LES_Model/Namelists/#initialization_parameters--spinup_time) and [data_output_during_spinup](/Reference/LES_Model/Namelists/#initialization_parameters--data_output_during_spinup) = *.T.*. With this, the UV-model is invoked at each data-output step (see [dt_do2d_xy](/Reference/LES_Model/Namelists/#runtime_parameters--dt_do2d_xy)) during the spinup and outputs data. This way, no three-dimensional flow simulation needs to be carried-out (set [end_time](/Reference/LES_Model/Namelists/#runtime_parameters--end_time) to a small non-zero value), saving computational resources.

- Users that are interested in UV-radiation and other quantities from a three-dimensional flow simulation need to carry-out a standard simulation. In this case, the UV-model is also invoked during the flow simulation at each 2D x-y data-output step (see [dt_do2d_xy](/Reference/LES_Model/Namelists/#runtime_parameters--dt_do2d_xy)).

As mentioned before, the UV-radiation model can be switched-on and controlled with the namelist [uv_radiation_parameters](/Reference/LES_Model/Namelists/#uv-radiation-parameters). Therein, users can choose the approach how the wavelength-integrated UV-irradiance is computed, i.e. from the LOD1 approach where the incoming diffuse radiation is assumed to be isotropic, or from the spherical-angle dependent LOD2 approach where the incoming diffuse radiation can be also distributed anisotropically. The LOD1 approach requires [uv_integration_method](/Reference/LES_Model/Namelists/#uv_radiation_parameters--uv_integration_method) = *'from_irradiance'*, while the LOD2 approach requires [uv_integration_method](/Reference/LES_Model/Namelists/#uv_radiation_parameters--uv_integration_method) = *'from_radiance'*. Both approaches can be combined in a simulation, which requires [uv_integration_method](/Reference/LES_Model/Namelists/#uv_radiation_parameters--uv_integration_method) = *'from_irradiance from_radiance'*.
Further on, the number of Lambertian reflection steps at mutually visible horizontal and vertical surfaces can be specificed via the parameter [num_reflections](/Reference/LES_Model/Namelists/#uv_radiation_parameters--num_reflections). Sensitivity tests revealed that after 3 to 5 reflections the UV-irradiances does not significantly change anymore.

The UV-radiation model in PALM can be started for example using the following minimal example setup:

```fortran
&initialization_parameters
    dx         = 0.25,
    dy         = 0.25,
    dz         = 0.25,
    nx         = 199,
    ny         = 199,
    nz         = 120,

    latitude         = 50.0,
    longitude        = 0.0,
    origin_date_time = '2023-06-21 00:00:00 +02',
    rotation_angle   = 0.0,

    initializing_actions = 'set_constant_profiles',

    spinup_pt_amplitude = 1.0,
    spinup_pt_mean      = 278.15,
    spinup_time         = 86400.0,
    data_output_during_spinup = .T.,

    topography = 'read_from_file',

    allow_roughness_limitation = .T.,
/

&runtime_parameters
    end_time   = 0.000001, ! Simulation time of 3D simulation

    dt_do2d_xy         = 100.0,

    section_xy         = 1,
    skip_time_data_output = 0.0,
    netcdf_data_format    = 5,

    data_output = 'uv_ewir1*_xy',   'uv_ewir2*_xy',
                  'uv_ir1*_xy',     'uv_ir2*_xy',

/

&radiation_parameters
    radiation_scheme = 'clear-sky',
    dt_radiation = 100.0,

    raytrace_discrete_azims = 60,
    raytrace_discrete_elevs = 90,

    surface_reflections = .T.,
    localized_raytracing = .T.,
 /

 &land_surface_parameters
    soil_temperature       = 293.5, 293.6, 293.1 293.1, 293.1, 293.1 293.1, 293.1,
    soil_moisture          = 0.2,   0.2,   0.2,  0.2,   0.2,   0.2,   0.2,   0.2,
 /

 &urban_surface_parameters
 /


 &plant_canopy_parameters
    canopy_mode       = 'read_from_file',
 /

 &uv_radiation_parameters
   uv_integration_method = 'from_radiance from_irradiance',
   num_reflections = 3,
 /
```

## Output quantities

You may output 2d-horizontal arrays of the wavelength integrated and/or erythemally weighted UV irradiance that are calculated based on the LOD1 and LOD2 integration methods given via [uv_integration_method](/Reference/LES_Model/Namelists/#uv-radiation-parameters--uv_integration_method). See [table of UV output quantities](/Reference/LES_Model/Output_quantities/#uv-quantities) for allowed quantity names.

**Note**, time-averaged output of UV-related quantities is not possible.

