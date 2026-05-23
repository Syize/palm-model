---
title: Overview
---
# I/O Files Format Reference

---

The model works with a set of files, which are all located in the temporary working directory and which are either expected at the beginning of the run (the so-called input files) or which are produced during the run and/or at the end of the simulation (output files). All default files used by the model are listed below.

- `File Suffix` gives the suffix of the file in the respective user folder (e.g. `INPUT`, `OUTPUT`, or `MONITORING`). <br>In case of nested setups, each nest requires or generates a separate file that is labeled by an additional suffix `_N02` (for the first child), `_N03` (for the second child), etc. The root domain files have no additional suffix. <br>In case of coupled atmosphere-ocean runs (see [coupled atmosphere-ocean](../Modules/Ocean/#atmosphere-ocean-coupling)), the ocean files (need to) have the additional suffix `_O`. The atmosphere files do not have an additional suffix.

- `Local Name` gives the respective file name in the temporary working folder (this name is used in the respective `OPEN` statement of the Fortran code). By default, file names are always written in capital letters. If the name is followed by a slash (`/`), the local name may refer to a folder instead of a file, in which one file per core is expected. File names in the folder consist of a `_`, followed by the six digit core ID, e.g. `_000000`, `_000001`, `_000002` when running a job on three cores. These files contain data of the respective subdomains. Core dependent data in a folder are only expected / generated in case of [restart_data_format](../Namelists/#initialization_parameters--restart_data_format) = *'fortran_binary'*.

- `I/O` indicates whether it is an input or output file.

- `Format` indicates the file format.

For internal use, the model may open a set of further files, which do not contain any usable information. Some of them are contained in this list. No file suffix is given for these files.

<br>
<br>


## Namelist Parameters
| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_p3d`    | `PARIN`  | I | ASCII / Fortran NAMELIST |

Contains the setup parameters for [model steering](../Scripts/palmrun/palmrun_script.md). See [example_cbl_p3d](https://gitlab.palm-model.org/palm/model-/blob/master/tests/cases/example_cbl/INPUT/example_cbl_p3d) for a typical parameter set for the convective boundary layer.

<br>
<br>


## Static Driver

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_static`    | `PIDS_STATIC`  |    I      | netCDF |

The static input file contains all static information, such as topography, geographical latitude and longitude, surface properties, and vegetation information. More detailed information on individual input variables is provided [here](Drivers/static.md).

<br>
<br>


## Dynamic Driver

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :-------: | :------- |
| `_dynamic`    | `PIDS_DYNAMIC`  |    I      | netCDF |

The dynamic input file contains dynamic information on the initial state of the atmosphere or time-dependent boundary conditions. More detailed information on individual input variables is provided [here](Drivers/dynamic.md).

<br>
<br>


## Run Control

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_rc`    | `RUN_CONTROL`  | O | ASCII |

This file contains the so-called run control output of the model. At a given temporal interval, determined via the runtime parameter [dt_run_control](../Namelists/#runtime_parameters--dt_run_control), one record (line) with values of certain control parameters at the specific time step is written to this file. Additionally, a new record is always written, whenever the time step of the model has changed. All data and quantities given in this file always refer to the entire model domain.

If the 1D-model is switched on for the initialization of the 3D-model via [initializing_actions](../Namelists/#initialization_parameters--initializing_actions) = *'set_1d-model_profiles'*, records are likewise written into this file at temporal intervals given via [dt_run_control_1d](../Namelists/#runtime_parameters--dt_run_control_1d).

In case of initial runs, beside the time step information, file [`RUN_CONTROL`](#run-control) also contains information about the selected model setup (parameter values, etc.) at the beginning of the file. This information is written at the beginnning of the run. It corresponds to the content of file [`HEADER`](#header), except data concerning the consumed CPU time (because that is only known at the end of a run). Via runtime parameter [force_print_header](../Namelists/#runtime_parameters--force_print_header) = *.T.* setup information is also written for restart runs.

The meaning of the individual columns of the run control time step output is given in the table below. The first column of this table shows the respective heading of the appropriate column in the [`RUN_CONTROL`](#run-control) file: 

|             |   |
|-------------|---|
|  RUN        |  Serial-number of the job in the job chain. The initial run has the number 0, restart runs accordingly have incremented numbers.  |
|  ITER.      |  Number of time steps accomplished since the start of the initial run (`t`=*0*).  |
|  HH:MM:SS   |  Time (in hours: minutes: seconds) since the start of the initial run (`t`=*0*).  |
|  DT (E)     |  Time step (in s). The appended character indicates the reason of timestep limitation due to:<br> (A) advection (CFL criterion)<br> (D) diffusion criterion<br> (S) diffusion criterion in surface energy balance models (see [land surface](../Namelists/#land-surface-parameters) and [urban surface](../Namelists/#urban-surface-parameters) parameters)<br> (P) precipitaton (see [bulk cloud model parameters](../Namelists/#bulk-cloud-parameters))<br> (X) setting of [dt_max](../Namelists/#runtime_parameters--dt_max)<br> (F) fixing the time step via setting [dt](../Namelists/#runtime_parameters--dt).  |
|  UMAX       |  Maximum absolute velocity (u-component) (in m/s). If a random disturbance has been added to this velocity component (see [create_disturbances](../Namelists/#runtime_parameters--create_disturbances)), the character D is appended to the velocity value.  |
|  VMAX       |  Maximum absolute velocity (v-component) (in m/s). If a random disturbance has been added to this velocity component (see [create_disturbances](../Namelists/#runtime_parameters--create_disturbances)), the character D is appended to the velocity value.  |
|  WMAX       |  Maximum absolute velocity (w-component) (in m/s).  |
|  U*         |  Horizontal average of friction velocity in the constant flux layer (in m/s).  |
|  W*         |  Convective velocity scale (in m/s). The assumed boundary layer height is determined via the criterion of minimum heat flux.  |
|  THETA*     |  Horizontal average of characteristic temperature in the constant flux layer (in K).  |
|  Z_I        |  Height of the convective boundary layer (in m), determined via the heat flux minimum criterion.  |
|  ENERG.     |  Total energy of the resolved-scale velocity, i.e. mean flow plus resolved-scale turbulence (in m^2^/s^2^), normalized with the total number of grid points.  |
|  DISTENERG  |  Energy of the resolved-scale turbulence (in m^2^/s^2^), normalized with the total number of grid points. **Please note**: For small number of horizontal grid points it is possible that DISTENERG > ENERG, which can be attributed to insufficient statistics.  |
|  DIVOLD     |  Divergence of the velocity field (sum of absolute values at all grid points) (in 1/s) before the pressure solver has been called, normalized with the total number of grid points.  |
|  DIVNEW     |  Divergence of the velocity field (sum of absolute values at all grid points) (in 1/s) after the pressure solver has been called, normalized with the total number of grid points.  |
|  UMAX (KJI) |  Indices of the grid point with the maximum absolute u-component of the velocity (sequence: k, j, i).  |
|  VMAX (KJI) |  Indices of the grid point with the maximum absolute v-component of the velocity (sequence: k, j, i).  |
|  WMAX (KJI) |  Indices of the grid point with the maximum absolute w-component of the velocity (sequence: k, j, i).  |
|  ADVECX     |  Distance (in km) the coordinate system has been moved along x-direction with Galilei-Transformation switched on (see [galilei_tranformation](../Namelists/#initialization_parameters--galilei_transformation)).  |
|  ADVECY     |  Distance (in km) the coordinate system has been moved along y-direction with Galilei-Transformation switched on (see [galilei_tranformation](../Namelists/#initialization_parameters--galilei_transformation)).  |
|  MGCYC      |  Number of multigrid cycles used if [psolver](../Namelists/#initialization_parameters--psolver) = *'multigrid'* has been set (see also [mg_cycles](../Namelists/#initialization_parameters--mg_cycles)). |

<br>
The meaning of the individual columns of the 1D run cuntrol time step output is given in the table below. The first column of this table shows the respective heading of the appropriate column: 

|           |                                                                                                     |
|-----------|-----------------------------------------------------------------------------------------------------|
|  ITER.    |  Number of time steps accomplished so far.                                                          |
|  HH:MM:SS |  Time (in hours: minutes: seconds).                                                                 |
|  DT       |  Time step (in s).                                                                                  |
|  UMAX     |  Maximum absolute velocity (u-component) (in m/s).                                             |
|  VMAX     |  Maximum absolute velocity (v-component) (in m/s).                                             |
|  U*       |  Friction velocity (in m/s).                                                                        |
|  ALPHA    |  Angle of the wind vector (to the x-axis) at the top of the constant flux layer (k=1) (in degrees). |
|  ENERG.   |  Kinetic energy (in m^2^/s^2^ ), averaged over all grid points.                  |

<br>
<br>


## Header File

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_header`    | `HEADER`  |    O      | ASCII |

Information about the selected model parameters (physical and numerical values) as well as general information about the run.

This file contains the values of all important steering parameters (numerical procedures, computational grid and model dimensions, boundary conditions, physical dimensions, turbulence quantities, actions during the simulation, 1D-model-parameters), as well as data concerning the selected output quantities. The headlines of the file list the program version used, date and time of the beginning of the run, the name of the executing host, the run identifier, and the number of the run (number of the restart run). With parallel runs the number of cores as well as the assigned virtual processor grid are displayed, too. After these headlines run time and time step information appear (point of starting time, related to `t = 0` of the initial run, end-time, time actually reached, CPU time, etc.). If a model run is incorrectly terminated (e.g. run time error or excess of the permitted CPU time), information over the time reached and the necessary CPU time is missing (to understand: the file HEADER is written twice by the model; once briefly after beginning of the simulation (naturally here the information over the time reached is missing etc.) and then short before the normal end of the simulation. The second, now complete output overwrites the first output.).

The information for the core model is followed by information about the modules that have been activated (including the user-interface module). If a module is not displayed here, it has not been activated.

<br>
<br>


## CPU Measures

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_cpu`    | `CPU_MEASURES`  |    O      | ASCII |


Time measurements are accomplished through the subroutine `cpu_log.f90`. The file CPU_MEASURES contains a header with some basic information of the respective run (e.g. model version, run identifier, executing host, date, etc.), followed  by two tables with data of CPU times required by certain model parts. The parts analyzed in the first table do not overlap and the CPU times should therefore approximately sum up to the “total” value given in the first line of this table. In the second table (heading: special measures) overlaps, in particular with the parts listed in the first table, are possible.

For each model part, the columns show how much CPU time was needed (absolutely and relative), and how often the part was called. With runs on several cores, the average values over all cores are indicated. In the case of parallel runs, information for those cores with the largest and smallest CPU time follow, as well as the standard deviation of all cores from the average value. Assuming that the code parallelizes well, the CPU times on the individual cores should vary only little and the standard deviation should be small.

Below the two tables additional information may be given about I/O transfer speed and file sizes, both  for output 3d netCDF and restart data.

<br>
<br>


## Restart (Input)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_d3d`    | `BININ/`  | I | Binary |

Binary data, which are read by the model at the beginning of a [restart run](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/runs). The appropriate file must have been written by the preceding job of the job chain (see [`BINOUT`](#restart-output)) or by a precursor run. This file contains the [initial parameters](../Namelists/#initialization-parameters) of the job chain, arrays of the prognostic and diagnostic variables as well as those parameters and variables for plots of horizontally averaged vertical profiles (see [data_output_pr](../Output_quantities/#vertical-profile-quantities)), which have been determined by the job chain so far.

Depending on the [restart_data_format](../Namelists/#initialization_parameters--restart_data_format) a single file or a folder with core-dependent files is expected.

<br>
<br>


## Restart (Output)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_d3d`    | `BINOUT/`  | O | Binary |

Binary data, which are written at the end of a run and required for initializing the next [restart runs](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/runs). This output file is then read in as file [`BININ`](#restart-input). It contains the [initial parameters](../Namelists/#initialization-parameters) of the model run, arrays of the prognostic and diagnostic variables as well as those parameters determined so far during a job chain and variables for plots of horizontally averaged vertical profiles (see [data_output_pr](../Output_quantities/#vertical-profile-quantities)).

File [`BINOUT`](#restart-output) is only written if the activation string `restart` has been set via [palmrun](../Scripts/palmrun/palmrun_script.md) option `-a`. Internally, this will set the `ENVPAR` namelist variable `write_binary` to *.T.*.

With large number of grid points, the file size of [`BINOUT`](#restart-output) (or the sized of files residing in folder [`BINOUT/`](#restart-output)) may become very large and should be stored (if available) on file-systems usually provided on cluster-systems to store such files (i.e. not under `$HOME`).

<br>
<br>


## Vertical Profiles (ASCII)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_list_pr`    | `LIST_PROFIL`  |    O      | ASCII |

**Attention: This file is deprecated and will not be available in future releases.**

File contains horizontally (and possibly temporally) averaged vertical profiles of some model variables. The quantities saved are internally set and can not be determined by the user. At each output time (see [dt_dopr_listing](../Namelists/#runtime_parameters--dt_dopr_listing)) two tables are written to the file consecutively. The first contains variables which are defined on the scalar / u-v-grid-levels, the second table contains variables which are defined on the w-grid-levels. If subdomains have been defined via initialization parameter [statistic_regions](../Namelists/#initialization_parameters--statistic_regions)), the tables are output for each individual subdomain. The name of the respective subdomain (see [region](../Namelists/#user_parameters--region)) appears in the header of the respective table.

In each case the tables consist of a header, followed by the profiles in separate columns. The header contains some basic information about the respective run (model version, run identifier, number of the job in the job chain, executing host, date, time of the beginning of the run, name of the subdomain, output time, kind of averaging). On the u-v-level the following columns are written: 

|      |                                                                                                                   |
|------|-------------------------------------------------------------------------------------------------------------------|
|  k   |  Vertical grid point index.                                                                                       |
|  zu  |  Height of the grid point level (in m).                                                                           |
|  u   |  u-component of the velocity (in m/s).                                                                       |
|  du  |  Deviation of the u-component from the initial profile at the time `t = 0` (initialization profile) (in m/s).     |
|  v   |  v-component of the velocity (in m/s).                                                                       |
|  dv  |  Deviation of the v-component from the initial profile at the time `t = 0` (initialization profile) (in m/s).     |
|  pt  |  Potential temperature (in K).                                                                                    |
|  dpt |  Deviation of potential temperature from the initial profile at the time `t = 0` (initialization profile) (in K). |
|  e   |  Turbulent kinetic energy (subgrid-scale) (in m^2^/s^2^).                                                         |
|  Km  |  Eddy diffusivity for momentum (in m^2^/s).                                                        |
|  Kh  |  Eddy diffusivity for heat (in m^2^/s).                                                            |
|  l   |  Mixing length (in m).                                                                                            |

On the w-level the following columns are written: 

|        |                                                                                               |
|--------|-----------------------------------------------------------------------------------------------|
|  k     |  Vertical grid point index.                                                                   |
|  zw    |  Height of the grid point level (in m).                                                       |
|  w'pt' |  Vertical subgrid-scale kinematic heat flux (in K m/s).                                       |
|  wpt   |  Vertical total ( subgrid-scale + resolved) kinematic heat flux (in K m/s).                   |
|  w'u'  |  u-component of the vertical subgrid-scale momentum flux (in m^2^/s^2^).                      |
|  wu    |  u-component of the total vertical momentum flux ( subgrid-scale + resolved) (in m^2^/s^2^).  |
|  w'v'  |  v-component of the vertical subgrid-scale momentum flux (in m^2^/s^2^).                      |
|  wv    |  v-component of the total vertical momentum flux ( subgrid-scale + resolved) (in m^2^/s^2^).  |

<br>
<br>


## Vertical Profiles (ASCII, 1D Model)

**Attention: This file is deprecated and will not be available in future releases.**

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| -    | `LIST_PROFIL_1D`  |    O      | ASCII |

This file contains the vertical profiles calculated by the 1D-model at the beginning of initial runs. The given quantities are internally set and can not be determined by the user. At the respective output times (see [dt_pr_1d](../Namelists/#initialization_parameters--dt_pr_1d)) a table with the following information is written to the file: The table header contains some basic information of the respective run (model version, run identifier, number of the job in the job chain (always *00*, because the 1D-model is only switched on for initial runs), executing host, date, time of the beginning of the run, output time). Afterwards, the following columns appear: 

|      |                                              |
|------|----------------------------------------------|
|  k   |  Vertical grid point index.                  |
|  zu  |  Height of the grid point level (in m).      |
|  u   |  u-component of the velocity (in m/s).       |
|  v   |  v-component of the velocity (in m/s).       |
|  pt  |  Potential temperature (in K).               |
|  e   |  Turbulent kinetic energy (in m^2^/s^2^).    |
|  rif |  Flux Richardson number.                     |
|  Km  |  Eddy diffusivity for momentum (in m^2^/s).  |
|  Kh  |  Eddy diffusivity for heat (in m^2^/s).      |
|  l   |  Mixing length (in m).                       |

<br>
<br>


## Vertical Profiles (netCDF)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_pr`    | `DATA_1D_PR_NETCDF`  |    O      | netCDF |

This file contains data horizontally averaged vertical profiles (see [data_output_pr](../Namelists/#runtime_parameters--data_output_pr)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Vertical Profiles (ASCII, for STG)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_stg`    | `STG_PROFILES`  |    I      | ASCII |

This file contains vertical profiles to be used as input for the [synthetic turbulence generator](../Modules/Synthetic_turbulence/). The first line consists of a header which describes the meaning of the columns in which the data is given. This header line is ignored by PALM. The quantities to be given in the colummns are (from left to right: 

|         |                                              |
|---------|----------------------------------------------|
| k       |  vertical grid point index                   |
| L(u,y)  |  length scale of u along y-direction         |
| L(u,z)  |  length scale of u along z-direction         |
| T(u)    |  time scale of u                             |
| L(v,y)  |  length scale of v along y-direction         |
| L(v,z)  |  length scale of v along z-direction         |
| T(v)    |  time scale of v                             |
| L(w,y)  |  length scale of w along y-direction         |
| L(w,z)  |  length scale of w along z-direction         |
| T(w)    |  time scale of w                             |
| uu      |  Reynolds stress component R11               |
| uv      |  Reynolds stress component R12               |
| vv      |  Reynolds stress component R22               |
| uw      |  Reynolds stress component R13               |
| vw      |  Reynolds stress component R23               |
| ww      |  Reynolds stress component R33               |
| U       |  mean wind speed along x-direction           |
| V       |  mean wind speed along y-direction           |
| W       |  mean wind speed along z-direction           |
| TKE     |  mean subgrid-scale turbulent kinetic energy |

See [STG_PROFILES](https://palm.muk.uni-hannover.de/trac/attachment/wiki/doc/app/iofiles/STG_PROFILES) for an example input file.

<br>
<br>


## Time Series

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_ts`    | `DATA_1D_TS_NETCDF`  |    O      | netCDF |

This file contains time series data (see [dt_dots](../Namelists/#runtime_parameters--dt_dots)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Spectra

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_sp`    | `DATA_1D_SP_NETCDF`  |    O      | netCDF |

This file contains data of horizontal spectra (see [data_output_sp](../Namelists/#spectra_parameters--data_output_sp)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## XY Sections (Binary)

**Attention: This file is deprecated and will not be available in future releases.**

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| -    | `PLOT2D_XY_<core-id>`  |    O      | Binary |

This local file contains xy cross section data of the chosen output quantities (see [data_output](../Namelists/#runtime_parameters--data_output)) in Fortran binary format. It is internally read by the post-processing tool `combine_plot_fields.x`, that is automatically called by [palmrun](../Scripts/palmrun/palmrun_script.md) after the model has finished.

With parallel runs and setting of [data_output_2d_on_each_pe](../Namelists/#runtime_parameters--data_output_2d_on_each_pe) = *.T.* and [netcdf_data_format](../Namelists/#runtime_parameters--netcdf_data_format) <= *4*, each core writes the data of its subdomain to a separate file with name `PLOT2D_XY_<core-id>`, where `<core-id>` is a six digit number (e.g. `PLOT2D_XY_000000`). These individual files are read and sampled into one final by `combine_plot_fields.x`. The tool writes informative messages about performed actions to the job protocol (also in case that no files have been found). `combine_plot_fields.x` also treats files of other cross sections (xz and/or yz) and 3d-output ([`PLOT3D_DATA`](#3d-data-binary)), if existing.

<br>
<br>


## XY Sections (netCDF)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_xy`    | `DATA_2D_XY_NETCDF`  |    O      | netCDF |

This file contains horizontal cross section data (see [data_output](../Namelists/#runtime_parameters--data_output)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## XY Sections (netCDF, averaged)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_av_xy`    | `DATA_2D_XY_AV_NETCDF`  |    O      | netCDF |

This file contains temporally averaged horizontal cross section data (see [data_output](../Namelists/#runtime_parameters--data_output)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## XZ Sections (Binary)

**Attention: This file is deprecated and will not be available in future releases.**

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| -    | `PLOT2D_XZ_<core-id>`  |    O      | Binary |

This local file contains xz cross section data of the chosen output quantities (see [data_output](../Namelists/#runtime_parameters--data_output)) in Fortran binary format.

For more information, see description of local file [`PLOT2D_XY`](#xy-sections-binary).

<br>
<br>


## XZ Sections (netCDF)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_xz`    | `DATA_2D_XZ_NETCDF`  |    O      | netCDF |

This file contains vertical (xz) cross section data (see [data_output](../Namelists/#runtime_parameters--data_output)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## XZ Sections (netCDF, averaged)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_av_xz`    | `DATA_2D_XZ_AV_NETCDF`  |    O      | netCDF |

This file contains vertical (xz) cross section data (see [data_output](../Namelists/#runtime_parameters--data_output)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## YZ Sections (Binary)

**Attention: This file is deprecated and will not be available in future releases.**

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| -    | `PLOT2D_YZ_<core-id>`  |    O      | Binary |

This local file contains yz cross section data of the chosen output quantities (see [data_output](../Namelists/#runtime_parameters--data_output)) in Fortran binary format.

For more information, see description of local file [`PLOT2D_XY`](#xy-sections-binary).

<br>
<br>


## YZ Sections (netCDF)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_yz`    | `DATA_2D_YZ_NETCDF`  |    O      | netCDF |

This file contains vertical (yz) cross section data (see [data_output](../Namelists/#runtime_parameters--data_output)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## YZ Sections (netCDF, averaged)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_av_yz`    | `DATA_2D_YZ_AV_NETCDF`  |    O      | netCDF |

This file contains vertical (yz) cross section data (see [data_output](../Namelists/#runtime_parameters--data_output)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## 3D Data (Binary)

**Attention: This file is deprecated and will not be available in future releases.**

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| -    | `PLOT3D_DATA_<core-id>`  |   O     | Binary |

This local file contains 3d-data in Fortran binary format of quantities selected via [data_output](../Namelists/#runtime_parameters--data_output)).

For more information, see description of local file [`PLOT2D_XY`](#xy-sections-binary).

<br>
<br>


## 3D Data (netCF)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_3d`    | `DATA_3D_NETCDF`  |    O      | netCDF |

This file contains 3d-data in netCDF format of quantities selected via [data_output](../Namelists/#runtime_parameters--data_output).

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## 3D Data (netCDF, averaged)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_av_3d`    | `DATA_3D_AV_NETCDF`  |    O      | netCDF |

This file contains time averaged 3d-data in netCDF format of quantities selected via (see [data_output](../Namelists/#runtime_parameters--data_output)).

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Masked Data

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_masked_XX`  | `DATA_MASK_XX_NETCDF`  |    O      | netCDF |

This file contains masked data in netCDF format of quantities set via [data_output_masks](../Namelists/#runtime_parameters--data_output_masks). `XX` is the respective mask number (`01`, `02`, etc.). One file is generated for each mask.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Masked Data (averaged)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_av_masked_XX`  | `DATA_MASK_XX_AV_NETCDF`  |    O      | netCDF |

This file contains temporally averaged masked data in netCDF format of quantities set via [data_output_masks](../Namelists/#runtime_parameters--data_output_masks). `XX` is the respective mask number (`01`, `02`, etc.). One file is generated for each mask.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Surface Data (for Paraview)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_surf_bin`    | `SURFACE_DATA_BIN/`  |   O     | Binary |

These files contains unstructured surface data at time intervals given by [dt_dosurf](../Namelists/#surface_data_output_parameters--dt_dosurf), and are generated if [to_vtk](../Namelists/#surface_data_output_parameters--to_vtk) = *.T.* has been set. Data needs to be further processed to convert it to Paraview-readable ASCII files (see postprocessing routine [surface_output_to_vtk](https://palm.muk.uni-hannover.de/trac/browser/palm/trunk/UTIL/surface_output_processing/surface_output_to_vtk.f90)).

<br>
<br>


## Surface Data (for Paraview, Averaged)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_av_surf_bin`    | `SURFACE_DATA_AV_BIN/`  |   O     | Binary |

These files contains unstructured surface data at time intervals given by [dt_dosurf_av](../Namelists/#surface_data_output_parameters--dt_dosurf_av), and are generated if [to_vtk](../Namelists/#surface_data_output_parameters--to_vtk) = *.T.* has been set. Data needs to be further processed to convert it to Paraview-readable ASCII files (see postprocessing routine [surface_output_to_vtk](https://palm.muk.uni-hannover.de/trac/browser/palm/trunk/UTIL/surface_output_processing/surface_output_to_vtk.f90)).

<br>
<br>


## Surface Data (for Spinup, Input)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_spinup`    | `SPINUPIN`  |    I      | Binary (MPI) |

Binary data for storing spinup surface data.

<br>
<br>


## Surface Data (for Spinup, Output)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_spinup`    | `SPINUPOUT`  |    O      | Binary (MPI) |

Binary data for storing spinup surface data.

<br>
<br>


## Surface View Factor (Input)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_svf`    | `SVFIN/`  |    I      | Binary |

Binary data for storing surface view factors needed by the RTM.

<br>
<br>


## Surface View Factor (Output)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_svf`    | `SVFOUT/`  |    O      | Binary |

Binary data for storing surface view factors generated by the RTM.

<br>
<br>


## Particle Statistics for Debugging

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_prt_info`    | `PARTICLE_INFOS/`  |    O      | ASCII |

This file is for debugging purposes and is created in case the Lagrangian particle model is active and [write_particle_statistics](../Namelists/#particle_parameters--write_particle_statistics) = *.T.* has been set. It contains statistical information about the number of particles that have been generated, and about the number of particles that have been exchanged between the respective neighbouring PEs. This information is output after every timestep.

<br>
<br>


## Particle Data (Binary)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_prt_bin`    | `PARTICLE_DATA/`  |    O     | Binary |

This file is generated in case the Lagrangian particle model is active and [dt_write_particle_data](../Namelists/#particle_parameters--dt_write_particle_data) has been set. One file per core (subdomain) is created.

The first record of a file contains an identification string (PALM release, run identifier, etc., 80 characters long). The second record also contains a version string (80 characters long) followed by two records with particle group information and a record containing the index bounds of the 3d subdomain arrays (see source code `file check_open.f90`). Then, for each output time, two records follow, where the first one gives the output time and the second one the number of particles that are output for the subdomain. After that, particle data is output, one record per grid box (if particles exist in that box) of the respective subdomain. See source code file `lagrangian_particle_model_mod.f90` for the respective output statements. A Fortran TYPE structure is used for storing the particle attributes. See source code file `mod_particle_attributes.f90` for the detailed TYPE structure.

To read data from this file requires a Fortran program with `READ` statements that exactly match the `WRITE` statements that have been used to generate the records described above.

<br>
<br>


## Particle Restart (Input)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_rprt`    | `PARTICLE_RESTART_DATA_IN/`  |    I      | Binary |

Binary data read at the beginning of a [restart run](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/runs). The respective file must have been generated by the precedent job of a job chain (see [`PARTICLE_RESTART_DATA_OUT`](#particle-restart-output). This file is required for restart runs if the Lagrangian particle model is switched on. It contains all particle information (particle positions, velocities, etc.) from the end of the preceding run.

In case of [restart_data_format](../Namelists/#initialization_parameters--restart_data_format) = *'fortran_binary'*, the number of cores used in a restart run must be identical to the number used in the precedent run.

<br>
<br>


## Particle Restart (Output)

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_rprt`    | `PARTICLE_RESTART_DATA_OUT/`  |    O      | Binary |

Binary data written at the end of a run, to be used by a restart run, where it is read as file [`PARTICLE_RESTART_DATA_IN`](#particle-restart-input). The file is only written if the Lagrangian particle model is active.

<br>
<br>


## Large Scale Forcing

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_lsf`    | `LSF_DATA`  |    I      | ASCII |

The file contains column-separated surface variables `shf`, `qsws`, `pt_surface`, `q_surface`, `surface_pressure`, and profiles of `ug`, `vg`, `w_subs`, `td_lsa_lpt`, `td_lsa_q` (large scale advection tendencies), `td_sub_lpt`, `td_sub_q` (large scale subsidence tendencies), at different times. The data is usually taken from / provided by measurements or larger scale models.

In case of [large_scale_forcing](../Namelists/#initialization_parameters--large_scale_forcing) = *.T.*, data is read in free floating point format. The hash character (`#`) is used as a special steering character. The data in the file has to be organized in the following way: 

- The file content has to start with **exactly** three lines, beginning with a `#`, where for example header information can be given.
- The first line not beginning with a `#` contains surface data which must be column-separated in the order: `time`, `shf`, `qsws`, `pt_surface`, `q_surface`, `surface_pressure`.
- At least one blank line after the surface data must follow to indicate the end of the surface data.
- A line beginning with `# *time*` indicates the beginning of a profile data set at the given time.

Profile data following `# *time*` has to be column-separated in the order: `zu` (height), `ug`, `vg`, `w_subs`, `td_lsa_lpt`, `td_lsa_q`, `td_sub_lpt`, `td_sub_q`.

The heights given in column `zu` do not have to match PALM's grid. Data are linearly interpolated to the `zu`-grid of PALM. Also, linear interpolation in time of the surface and profile data is done, if required.

See below an example for the general file layout:
```
# Data obtained from ~/hdcp2/COSMO/042013_p0.nc, t=561-567h
#          time        shf           qsws  pt_surface  q_surface  surface_pressure
#           (s)    (K m/s)    (m/s kg/kg)         (K)    (kg/kg)             (hPa)
           0.00    0.03998    0.000014616    291.1135   0.007065         1011.1792
       10800.00    0.18522    0.000037693    296.3799   0.007050         1010.5699
       21600.00    0.16072    0.000043956    295.7625   0.007415         1009.1753
            .         .           .             .           .                .
            .         .           .             .           .                .
            .         .           .             .           .                .


#     zu (m)    ug (m/s)    vg (m/s)    w_subs (m/s)  td_lsa_lpt (K/s)  td_lsa_q (kg/kgs)  td_sub_lpt (K/s)  td_sub_q (kg/kgs)
# 0.00
      9.8874    2.317242    2.532310   -0.012626      0.0000485535      0.0000000021       0.0000000         0.0000000000
     35.2900    2.317242    2.532310   -0.015071      0.0000458420      0.0000000020      -0.0001100        -0.0000000241
     72.1618    2.739493    2.652984   -0.017274      0.0000426586      0.0000000015      -0.0000482        -0.0000000126
    120.8500    3.179652    2.856771   -0.017966      0.0000397725      0.0000000011      -0.0000271        -0.0000000087
    181.7068    3.856441    3.111020   -0.017831      0.0000374494      0.0000000010      -0.0000176        -0.0000000073
    255.0584    4.699064    3.438433   -0.017074      0.0000356241      0.0000000013      -0.0000120        -0.0000000070
       .           .            .           .              .                 .                 .                  .
       .           .            .           .              .                 .                 .                  .
   6195.9199   12.441278   -2.188596   -0.013782     -0.0001402120     -0.0000000303       0.0000536        -0.0000000109
   6673.4595   12.857915   -1.855247   -0.016264      0.0000284057     -0.0000000670       0.0000785        -0.0000000096
       .           .            .           .              .                 .                 .                  .
       .           .            .           .              .                 .                 .                  .


# 10800.00
      9.8882    2.932734    0.453175   -0.014016      0.0000159810      0.0000000443       0.0000000         0.0000000000
     35.2749    2.932734    0.453175   -0.017078      0.0000123939      0.0000000508      -0.0001734        -0.0000000402
     72.1310    3.293148    0.409822   -0.020240      0.0000090166      0.0000000533      -0.0000818        -0.0000000196
    120.8031    3.548535    0.466816   -0.022120      0.0000057484      0.0000000544      -0.0000493        -0.0000000128
    181.6363    3.946456    0.530473   -0.023343      0.0000028124      0.0000000547      -0.0000348        -0.0000000101
    254.9711    4.402065    0.498268   -0.024055      0.0000001112      0.0000000540      -0.0000267        -0.0000000088
       .           .            .           .              .                 .                 .                  .
   6195.4304   10.566041   -2.318027    0.014784      0.0000591095     -0.0000000120      -0.0000498         0.0000000040
   6672.8375   12.129893   -2.050756   -0.001662     -0.0001437306      0.0000000326       0.0000084        -0.0000000010
       .           .            .           .              .                 .                 .                  .
       .           .            .           .              .                 .                 .                  .
```

An complete example parameter file (p3d) with large scale forcing can be found [here](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/examples/lsf).

<br>
<br>


## Nudging

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_nudge`    | `NUDGING_DATA`  |    I      | ASCII |

The file contains column-separated profiles of the prognostic variables `u`, `v`, `w`, `lpt`, `q`, taken from measurements or larger scale models, to be used for nudging. Additionally, height information and the nudging time scale `tnudge` have to be provided.

In case of [nudging](../Namelists/#initialization_parameters--nudging) = *.T.* the data are read in free floating point format. The hash character `#` is used as a special steering character. Data in the file has to be organized as follows:

- The file content can start with an arbitrary number of lines beginning with a `#`, where for example header information can be stored.
- A line beginning with `# *time*` indicates the beginning of a profile data set at the given time.
- Profile data following `# *time*` has to be column separated in the order: `zu` (height), `tau` (nudging time scale), `u`, `v`, `w`, `lpt`, `q`.

The heights given in column `zu` do not have to match PALM's grid. Data are linearly interpolated to the `zu`-grid of PALM. Also, linear interpolation in time of the surface and profile data is done, if required

If quantities are set to *-999999.9* for all heights and all time levels, no nudging is applied to them.

See below an example for the general file layout:
```
# Data obtained from ~/hdcp2/COSMO/042013_p0.nc, t=561-567h
#    zu (m)   tau (s)     u (m/s)     v (m/s)     w (m/s)      lpt (K)     q (kg/kg)
# 0.00
     9.8874   21600.0    3.623120    1.888156    0.000000   287.692993   0.006222516
    35.2900   21600.0    4.259634    2.231079    0.000000   287.383728   0.006158236
    72.1618   21600.0    4.606365    2.421291    0.000000   287.238007   0.006122927
   120.8500   21600.0    4.876022    2.568060    0.000000   287.144562   0.006095806
   181.7068   21600.0    5.129545    2.700137    0.000000   287.072601   0.006069813
   255.0584   21600.0    5.402528    2.828126    0.000000   287.012238   0.006040839
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
  6195.9199   21600.0   11.302397   -2.185184    0.000000   314.965851   0.000958589
  6673.4595   21600.0   10.486557   -2.779184    0.000000   317.305450   0.000541379
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .


# 10800.00
     9.8882   21600.0    4.122205    2.008313    0.000000   291.437469   0.005792995
    35.2749   21600.0    4.786541    2.331582    0.000000   291.013519   0.005695062
    72.1310   21600.0    5.142901    2.507198    0.000000   290.804626   0.005646179
   120.8031   21600.0    5.404265    2.636267    0.000000   290.667511   0.005611916
   181.6363   21600.0    5.628885    2.745806    0.000000   290.560150   0.005582496
   254.9711   21600.0    5.841899    2.846067    0.000000   290.467194   0.005553886
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
  6195.4304   21600.0   11.760574    0.330661    0.000000   314.826996   0.001249985
  6672.8375   21600.0   11.387577   -0.290878    0.000000   316.955536   0.000936071
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
      .           .            .         .       .          .          .
```
An example parameter file (`_p3d`) with nudging can be found [here](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/examples/lsf).

<br>
<br>


## Topography

**Attention: This file is deprecated and will not be available in future releases.**

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_topo`    | `TOPOGRAPHY_DATA`  |    I      | ASCII |

Two-dimensional rastered topography height information (in m above ground).

In case of [topography](../Namelists/#initialization_parameters--topography) = *'read_from_file'* topography height information (in m above ground) is read for each grid point in a free floating point format from this file, if no static driver file [`PIDS_STATIC`](#static-driver) is provided. The ASCII file format is [ESRI grid](http://en.wikipedia.org/wiki/ESRI_grid#ASCII) - also known as [ARC/INFO ASCII GRID](http://en.wikipedia.org/wiki/ESRI_grid#ASCII) - **without the header**. The data on file is laid out naturally, i.e. in W-E orientation horizontally and in S-N orientation vertically, and must thus be organized as follows: 

- each line contains height information in m from i = 0, ..., nx,
- the top line contains height information in m for j = ny (North), the bottom line for j = 0 (South),
- individual data must be separated by at least one blank.

Layout sketch:

```
                        N
    (0,ny)   (1,ny)   (2,ny)   ... (nx,ny)       top of file
    (0,ny-1) (1,ny-1) (2,ny-1) ... (nx,ny-1)  
W   (0,ny-2) (1,ny-2) (2,ny-2) ... (nx,ny-2)   E 
                        : 
                        : 
    (0,0)    (1,0)    (2,0)    ... (nx,0)        bottom of file 
                        S
```

 Example for a 50m tall building surrounded by a 12.5m tall podium on flat ground: 

```
   0    0    0    0    0    0    0    0    0    0    0    0    0
   0 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5    0    0
   0 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5    0    0
   0 12.5 12.5   50   50   50   50   50   50 12.5 12.5    0    0
   0 12.5 12.5   50   50   50   50   50   50 12.5 12.5    0    0
   0 12.5 12.5   50   50   50   50   50   50 12.5 12.5    0    0
   0 12.5 12.5   50   50   50   50   50   50 12.5 12.5    0    0
   0 12.5 12.5   50   50   50   50   50   50 12.5 12.5    0    0
   0 12.5 12.5   50   50   50   50   50   50 12.5 12.5    0    0
   0 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5    0    0
   0 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5    0    0
   0    0    0    0    0    0    0    0    0    0    0    0    0
```

These data must exactly match the horizontal grid. Due to the staggered grid, the topography may be displaced by -0.5 [dx](../Namelists/#initialization_parameters--dx) in x-direction and -0.5 [dy](../Namelists/#initialization_parameters--dy) in y-direction depending on the parameter [topography_grid_convention](../Namelists/#initialization_parameters--topography_grid_convention).

Alternatively, one may add code to the user interface subroutine `user_init_grid`.

Please note, that the recommended way to provide topography information is via the static driver file [`PIDS_STATIC`](#static-driver-input-file) in netCDF format.

<br>
<br>


## UV Radiation

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_uv`    | `PIDS_UV`  |    I      | netCDF |

The file contains static information on UV radiation properties such as sun-zenith angles or wavelengths. More detailed information on individual input variables is provided [here](PIDS_UV.md).

<br>
<br>


## Chemistry

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_chemistry`    | `PIDS_CHEM`  |    I      | netCDF |

The file contains static and dynamic information on chemical species and emissions. More detailed information on individual input variables is provided [here](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chememi).

<br>
<br>

## Traffic driver for Traffic module

| File Suffix | Local Name     | I/O  | Format      |
|:------------|:---------------| :------- |:--------|
| `_traffic`  | `PIDS_TRAFFIC` |    I     | netCDF  |

File contains information about traffic. The detailed information on individual input variables is provided [here](Drivers/traffic.md).

<br>
<br>


## Virtual Measurement Setup

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_vmeas`    | `PIDS_VM`  |    I      | netCDF |

File to setup virtual measurement locations. The input file contains the coordiate information of a measurement, its type, the sampled variables, as well as further attributes. The radiation input file contains static and dynamic information on chemical species and emissions. More detailed information on individual input variables is provided [here](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/vm).

<br>
<br>


## Wind Turbine Data

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :-------: | :------- |
| `_wtm`    | `WTM_DATA`  |    I      | ASCII |

This file contains data for PALM's built-in wind turbine model. The tabulated data is from the publicly available NREL 5 MW reference turbine (see Jonkman et al. 2009: Technical Report NREL/TP-500-38060, doi:10.2172/947422).

<br>
<br>


## Particle Time Series

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_pts`    | `DATA_1D_PTS_NETCDF`  |    O      | netCDF |

This file contains time series of particle quantities (see [dt_prel](../Namelists/#particle_parameters--dt_prel)) in netCDF format.

In case of using more than one particle group (see [number_of_particle_groups](../Namelists/#particle_parameters--number_of_particle_groups)), separate time series are output for each of the groups. The long names of the variables in the netCDF file containing the respective time series all end with the string `PG ##`, where `##` is the number of the respective particle group (01, 02, etc.).

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Agent Data

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_agt`    | `DATA_AGT_NETCDF `  |    O      | netCDF |

This file contains all output data from the [Multi Agent System (MAS)](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/mas) in netCDF format.

Each variable other than *time* (1D, contains what the name suggests) has two dimensions (*time* and *agent_number*) so that at each output time step defined via [dt_write_agent_data](../Namelists/#agent_parameters--dt_write_agent_data) all agent variables are output for each agent.

The variables are: 

- unique agent ID (ag_id)
- agent position (ag_x, ag_y)
- wind speed at agent position (ag_wind)
- temperature at agent position (ag_temp)
- agent group (ag_group, see [number_of_agent_groups](../Namelists/#agent_parameters--number_of_agent_groups))

**Very important**: Due to agent data structure, agents may not always be sorted in the same succession. E.g. `ag_id( time = i, agent_number = j)` is not guaranteed to be the same as `ag_id(time = i + 1, agent_number = j)`. Thus, if individual agent data is needed, sorting by `ag_id` in postprocessing is required.

For further information, have a look at parameters [dim_size_agtnum_manual](../Namelists/#agent_parameters--dim_size_agtnum_manual), [dim_size_factor_agtnum](../Namelists/#agent_parameters--dim_size_factor_agtnum) and [dt_write_agent_data](../Namelists/#agent_parameters--dt_write_agent_data).

<br>
<br>


## Agent Navigation Data Input File

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_nav`    | `NAVIGATION_DATA `  |    I      | Binary |

This file contains the navigation data from the [Agent Preprocessing Tool](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/mas/agent_preprocessing) in Fortran binary format.

<br>
<br>


## Virtual Flight Measurements

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :------- | :------- |
| `_vf`    | `DATA_1D_FL_NETCDF`  |    O      | netCDF |

This file contains data of the space-time series obtained from virtual flight measurements (see [virtual_flight_parameters](../Namelists/#virtual-flight-parameters)) in netCDF format.

See the description of [PALM-netCDF-output](../../../Guide/LES_Model/NetCDF_data_output/steering.md) for more detailed information.

<br>
<br>


## Wind Turbine Time Series

| File Suffix | Local Name  | I/O  | Format   |
| :-------- | :------- | :-------: | :------- |
| `_wtm`    | `DATA_1D_TS_WTM_NETCDF`  |    O      | netCDF |

This file contains time series in netCDF format of wind turbine quantities for each defined turbine at time intervals set via [dt_data_output_wtm](../Namelists/#wind_turbine_parameters--dt_data_output_wtm). Data of following quantities are output.

|             |                                           |
|-------------|-------------------------------------------|
| Time        |  simulated time in s                      |
| RSpeed      |  rotor speed in rad/s                     |
| GSpeed      |  generator speed in rad/s                 |
| GenTorque   |  generator torque in Nm                   |
| AeroTorque  |  aerodynamic torque in Nm                 |
| Pitch       |  pitch angle of the rotor blades in °     |
| Power(Gen)  |  electrical generator power in W          |
| Power(Rot)  |  mechanical power (rotor) in W            |
| RotThrust   |  rotor thrust in N                        |
| WDirection  |  wind direction (0° = west)               |
| YawOrient   |  yaw angle (0° = perpendicular to x axis) |

<br>
<br>


## Parallel I/O
On multicore systems, many of the files are read and/or written by one core only (processing element 0, PE0). These files have processor-independent content (and therefore they can be read or written by other PEs as well). However, certain files may have processor-dependent content. For the binary I/O of restart data (e.g. local file [BINOUT](#BINOUT)), each PE reads/writes only the data of its subdomain. So each processing element writes into its own file with its own file name. These files are located in a subdirectory of the temporary working directory. The respective file names are built from the underline ("_") and the six digit processor ID. The data written for restart runs would be e.g. written to the files BINOUT/_000000 (PE0), BINOUT/_000001 (PE1), BINOUT/_000002 (PE2) etc. Such files with processor-dependent content on parallel machines are marked in the above list by a slash character  ("/") at the end of the local file name. If appropriate output files are to be copied through [palmrun](../Scripts/palmrun/palmrun_script.md) to permanent files, and/or files with processor-dependent content are supposed to be copied as input files into the temporary working directory, you have to indicate a special file attribute in the appropriate file connection statement (see `arpe`, `flpe` in the [palmrun](../Scripts/palmrun/palmrun_script.md) description). Then the permanent file name will be interpreted as a directory name, in which the input files are expected and/or to which output files are written. The file names in these directories are always named _0000, _0001, _0002 etc.

In case of [restart_data_format](../Namelists/#initialization_parameters--restart_data_format) = *'fortran_binary'*, depending on the configuration of the underlying file-system (e.g. Lustre) and the capacity of the I/O-hardware, simultaneous output to a larger number of files (i.e. if a larger number of PEs >1000 is used) may lead to severe problems as job aborts or even to a complete crash of the whole system. In order to circumvent this problem, the maximum number of parallel I/O streams (i.e. the number of files which are simultaneously written), can be limited by using the [palmrun](../Scripts/palmrun/palmrun_script.md)-option `-w <max # of streams>`, where `<max # of streams>` should be smaller or equal e.g. *256*. The (parallel) output is then done in a sequential order for blocks of *256* PEs.


## Other available file format descriptions

- [Drivers](Drivers/index.md)
- [UV-radiation](PIDS_UV.md)
- [SLUrb driver](PIDS_SLURB.md)

