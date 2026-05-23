---
title: Masked data output
---

# Masked data output
---

This feature allows to output quantities at different mask locations, i.e. arbitrary subsets of the total model domain. Subsets can be 3d volumes, 2d cross sections, or even 0d or 1d data at any position and of any amount.
NetCDF4 parallel I/O is not available for masked data output.
Default quantities (instantaneous and time averaged) are declared with the parameter [data_output_masks](../../../../Reference/LES_Model/Namelists/#runtime_parameters--data_output_masks), user-defined quantities can be output with [data_output_masks_user](../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_masks_user) (see also [user-defined output quantities](../../LES_Model/Modules/User_Interface/output.md#user-defined-output-quantities)).

Terrain-following masked output allows to output masked quantities at a specified height above the surface (see [mask_k_over_surface](../../../../Reference/LES_Model/Namelists/#user_parameters--mask_k_over_surface)). 

## Output Steering

The [runtime parameters](../../../../Reference/LES_Model/Namelists/#runtime-parameters) that are listed below steer the output of those quantities that have been set via [data_output_masks](../../../../Reference/LES_Model/Namelists/#runtime_parameters--data_output_masks) and/or [data_output_masks_user](../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_masks_user):

 Following parameters define the masks. Each mask can be steered with a separate set of x-, y- and z-parameters. By default all gridpoints along the respective direction are output.

| | |
|---|---|
| [mask_x](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_x)               | all x-coordinates of mask locations (in multiples of mask scale)                                        |
| [mask_y](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_y)               | all y-coordinates of mask locations (in multiples of mask scale)                                        |
| [mask_z](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_z)               | all z-coordinates of mask locations (in multiples of mask scale)                                        |
| [mask_x_loop](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_x_loop)          | loop begin, end and stride for x-coordinates of mask locations for masks (in multiples of mask scale)   |
| [mask_y_loop](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_y_loop)          | loop begin, end and stride for y-coordinates of mask locations for masks (in multiples of mask scale)   |
| [mask_z_loop](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_z_loop)          | loop begin, end and stride for z-coordinates of mask locations for masks (in multiples of mask scale)   |
| [mask_k_over_surface](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_k_over_surface)  | vertical grid index above a surface to use for terrain-following masked output                          |

<br>
The scaling lengths can be used to scale the parameters which defines the masks (those listed above). The scaling lengths apply for all masks.

| | |
|---|---|
| [mask_scale_x](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_x) | scaling length for masked data output in x-direction |
| [mask_scale_y](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_y) | scaling length for masked data output in y-direction |
| [mask_scale_z](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_z) | scaling length for masked data output in z-direction |

<br>
The time intervals of the output times for each mask are determined with [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask). An individual time interval for output of temporally averaged data can be assigned via parameter [dt_data_output_av](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_data_output_av). The length of the averaging interval is controlled by [averaging_interval](../../../../Reference/LES_Model/Namelists/#runtime_parameters--averaging_interval). Output of masked data can be switched off until the interval defined via [skip_time_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--skip_time_domask) has passed.

No particular parameters are available for steering the time averaged output of each of the masks separately.


By default, a maximum number of *300* different masks can be defined. Each mask is output to one separate local file: 
    
   - Instantaneous data of mask 1 is output to file `DATA_MASK_NETCDF_M001` in the local temporary working directory (permanent file suffix `_masked_M001`), time averaged data to file `DATA_MASK_AV_NETCDF_M001` (`_av_masked_M001`).
   - Instantaneous data of mask 2 is output to file `DATA_MASK_NETCDF_M002` in the local temporary working directory (permanent file suffix `_masked_M002`), time averaged data to file `DATA_MASK_AV_NETCDF_M002` (`_av_masked_M002`).
   - ...

If more than *300* masks shall be defined, the Fortran parameter `max_masks` that is defined in source code file [modules.f90](https://gitlab.palm-model.org/palm/model/-/blob/master/src/modules.f90) needs to be increased manually and the model has to be re-compiled via `palmbuild`. This way, `max_masks` can be increased to a maximum value of *999*.

## Examples

The following examples are given based on the example setup for a convective boundary layer [example_cbl_p3d](https://gitlab.palm-model.org/palm/model/-/blob/master/tests/cases/example_cbl/INPUT/example_cbl_p3d) with a model domain size of *2000* m x *2000* m x *2000* m and a spatial resolution of *50.0* m. 

1. [Output of one mask](#1.-output-of-one-mask)
2. [Output of two different masks](#2.-output-of-two-different-masks)
3. [Output of three different masks](#3.-output-of-three-different-masks)
4. [Output of three different masks with scaling length](#4.-output-of-three-different-masks-with-scaling-length)
5. [Output of four different masks with user-defined quantities](#5.-output-of-four-different-masks-with-user-defined-quantities)
6. [Output of terrain-following mask](#6.-output-of-terrain-following-mask)

### 1. Output of one mask

**Mask 1**: 3d volume data with an extension of *500* m x *200* m x *2000* m from x = *0* m to x = *500* m, from y = *800* m to y = *1000* m, both at every *50* m, and from the bottom to the top of the model domain. Output is for instantaneous data of the three wind components.

```Fortran
&initialization_parameters

    nx = 39, ny = 39, nz = 40,
    dx = 50.0, dy = 50.0, dz = 50.0,
    .../

&runtime_parameters

    ...
    data_output_masks(1,:) = 'u','v','w',
   
    mask_x_loop(1,:) = 0., 500. ,50. ,
    mask_y_loop(1,:) = 800., 1000., 50. ,/
```
If [mask_x](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_x) , [mask_y](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_y) , [mask_z](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_z)  or [mask_x_loop](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_x_loop), [mask_y_loop](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_y_loop), [mask_z_loop](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_z_loop) are not set, all gridpoints along the corresponding direction are output.
Output is done at time intervals given by [dt_data_output](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_data_output), as [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask) has not been set in this example.

### 2. Output of two different masks

**Mask 1**: xy cross section from x = *500* m to x = *1000* m and from y = *1000* m to *2000* m, both at every *50* m, and at heights *0* m, *50* m, *100* m, *500* m, *1000* m, and *1200* m. Output of instantaneous and time averaged data of the three wind components.

**Mask 2**: Every second gridpoint in all directions. Output of instantaneous data of the potential temperature.

```Fortran
&initialization_parameters

    nx = 39, ny = 39, nz = 40,
    dx = 50.0, dy = 50.0, dz = 50.0,
    .../

&runtime_parameters

    ...
    data_output_masks(1,:) = 'u','v','w','u_av','v_av','w_av',
    data_output_masks(2,:) = 'theta',
   
    mask_x_loop(1,:) = 500.,1000.,50. ,
    mask_y_loop(1,:) = 1000.,2000.,50. ,
    mask_z(1,:) = 0.,50.,100.,500.,1000.,1200.,

    mask_x_loop(2,:) = 0.,2000.,100.,
    mask_y_loop(2,:) = 0.,2000.,100.,
    mask_z_loop(2,:) = 0.,2000.,100.,

    dt_domask = 600.,1800.,
    skip_time_domask = 3600.,3600.,/
```
Output starts after 1h simulated time (see [skip_time_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--skip_time_domask)) for both masks and is done every 10 minutes for **mask 1** and every 30 minutes for **mask 2** (see [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask)). 

### 3. Output of three different masks

**Mask 1**: 1d data along x-direction at y = *50* m, *200* m, *500* m and *1000* m and every fifth grid point in z-dircetion. Output of instantaneous data of the w-velocity component.  
**Mask 2**: Output of the potential temperature at position x = *0* m, y = *500* m, z = *200* m, *300* m and *400* m.  
**Mask 3**: xz cross section at y = *400* m, *450* m, *600* m. Output of time averaged data of the three wind components.

```Fortran
&initialization_parameters

    nx = 39, ny = 39, nz = 40,
    dx = 50.0, dy = 50.0, dz = 50.0,
    .../

&runtime_parameters

    ...
    data_output_masks(1,:) = 'w',
    data_output_masks(2,:) = 'theta',
    data_output_masks(3,:) = 'u_av','v_av','w_av',
  
    mask_y(1,:) = 50., 200., 500., 1000.,
    mask_z_loop(1,:) = 0.,2000.,250.,

    mask_x(2,:) = 0.,
    mask_y(2,:) = 500.,
    mask_z_loop(2,:) = 200.,400.,100.,

    mask_y(3,:) = 400., 450., 600.,

    dt_domask = 1800.,
    skip_time_domask = 0.,3600.,3600.,/
```
Output starts from the beginning for **mask 1** and after 1h simulated time (see [skip_time_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--skip_time_domask)) for **mask 2** and **3**. It is done every 30 minutes for **mask 1** and at time intervals of ​dt_data_output for **mask 2** and **3** (see [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask)).


### 4. Output of three different masks with scaling length

**Mask 1**, **2** and **3** as in [example 3](#3.-output-of-three-different-masks).

```Fortran
&initialization_parameters

    nx = 39, ny = 39, nz = 40,
    dx = 50.0, dy = 50.0, dz = 50.0,
    .../

&runtime_parameters

    ...
    data_output_masks(1,:) = 'w',
    data_output_masks(2,:) = 'theta',
    data_output_masks(3,:) = 'u_av','v_av','w_av',

    mask_scale_x = 10.0,
    mask_scale_y = 10.0,
  
    mask_y(1,:) = 5., 20., 50., 100.,
    mask_z_loop(1,:) = 0.,2000.,250.,

    mask_x(2,:) = 0.,
    mask_y(2,:) = 50.,
    mask_z_loop(2,:) = 200.,400.,100.,

    mask_y(3,:) = 40., 45., 60.,

    dt_domask = 1800.,
    skip_time_domask = 0.,3600.,3600.,/
```
Output starts from the beginning for **mask 1** and after 1h simulated time (see [skip_time_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--skip_time_domask)) for **mask 2** and **3**. It is done every 30 minutes for **mask 1** and at time intervals of ​dt_data_output for **mask 2** and **3** (see [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask)).

Since [mask_scale_x](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_x) and [mask_scale_y](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_y) are assigned, the parameters for steering the mask locations along x- and y-direction (here [mask_x](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_x) and [mask_y](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_y)) has to be given in multiples of [mask_scale_x](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_x) and [mask_scale_y](../../../../Reference/LES_Model/Namelists/#runtime_parameters--mask_scale_y).


### 5. Output of four different masks with user-defined quantities

**Mask 1**: 1d data along x-direction at y = *50* m, *200* m, *500* m and *1000* m and every fifth grid point in z-dircetion. Output of instantaneous data of the w-velocity component.  
**Mask 2**: Output of the potential temperature at position x = *0* m, y = *500* m, z = *200* m, *300* m and *400* m.  
**Mask 3**: xz cross section at y = *400* m, *450* m, *600* m. Output of time averaged data of the three wind components and the user-defined quantity *u2*.  
**Mask 4**: xy cross section at z = *100* m, *250* m, *500* m. Output of the user-defined quantity `u2`.

```Fortran
&initialization_parameters

    nx = 39, ny = 39, nz = 40,
    dx = 50.0, dy = 50.0, dz = 50.0,
    .../

&runtime_parameters

    ...
    data_output_masks(1,:) = 'w',
    data_output_masks(2,:) = 'theta',
    data_output_masks(3,:) = 'u_av','v_av','w_av',
  
    mask_y(1,:) = 50., 200., 500., 1000.,
    mask_z_loop(1,:) = 0.,2000.,250.,

    mask_x(2,:) = 0.,
    mask_y(2,:) = 500.,
    mask_z_loop(2,:) = 200.,400.,100.,

    mask_y(3,:) = 400., 450., 600.,

    mask_z(4,:) = 100., 250., 500.,

    dt_domask = 1800.,
    skip_time_domask = 0.,3600.,3600.,/

&user_parameters

    data_output_masks_user(3,:) = 'u2',
    data_output_masks_user(4,:) = 'u2',
```
Output starts from the beginning for **mask 1** and after 1h simulated time (see [skip_time_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--skip_time_domask)) for **mask 2**, **3** and **4**. It is done every 30 minutes for **mask 1** and at time intervals of [dt_data_output](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_data_output) for **mask 2**, **3** and **4** (see [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask)). 

### 6. Output of terrain-following mask

**Mask 1**: xy cross section from x = *500* m to x = *1000* m and from y = *1000* m to *2000* m at the second, fith and tenth grid point above the surface. Output of instantaneous and time averaged data of the three wind components.

```Fortran
&initialization_parameters

    nx = 39, ny = 39, nz = 40,
    dx = 50.0, dy = 50.0, dz = 50.0,
    .../

&runtime_parameters

    ...
    data_output_masks(1,:) = 'u','v','w','u_av','v_av','w_av',
   
    mask_x_loop(1,:) = 500.,1000.,50. ,
    mask_y_loop(1,:) = 1000.,2000.,50. ,
    mask_k_over_surface(1,:) = 2, 5, 10,

    dt_domask = 600.,
    skip_time_domask = 3600.,/
```
Output starts after 1h simulated time (see [skip_time_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--skip_time_domask)) for both masks and is done every 10 minutes for **mask 1** and every 30 minutes for **mask 2** (see [dt_domask](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_domask)). 
