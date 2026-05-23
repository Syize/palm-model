---
title: Overview
---
# Mesoscale Nesting

---

## Purpose

The mesoscale nesting module enables the simulation of spatially multi-scale atmospheric scenarios. It provides time-dependent boundary conditions at the lateral and upper boundaries of the LES model domain using data from an external source, typically a larger-scale numerical weather prediction model. This allows e.g. for the simulation of multi-day setups under evolving synoptic conditions. The approach effectively performs a downscaling from the mesoscale to the microscale, enabling the study of microscale phenomena while retaining the influence of broader synoptic patterns.


## General Information

The mesoscale nesting is activated by adding the namelist [`&nesting_offl_parameters`](../../../../Reference/LES_Model/Namelists/#offline-nesting-parameters) to the `_p3d` namelist file. The namelist is generally empty, except that the parameter [switch_off_module](../../../../Reference/LES_Model/Namelists/#nesting_offl_parameters--switch_off_module) has been set.

Since mesoscale model information typically lacks turbulence information, turbulence must first develop downstream of the boundary to inside the LES model. This development requires a certain fetch length, which depends on factors such as wind speed, boundary-layer depth, and atmospheric stability (Kadasch et al., 2021).
To accelerate turbulence development and reduce the required fetch length, it is highly recommended to use mesoscale nesting in conjunction with the [synthetic turbulence generator](../Synthetic_turbulence/index.md), which itself is enabled by adding the namelist [`&stg_par`](../../../../Reference/LES_Model/Namelists/#synthetic-turbulence-generator-parameters) to the `_p3d` namelist file. This configuration imposes spatially and temporally varying random perturbations onto the boundary values.
However, even with synthetic turbulence, the fetch length can still be considerable, typically around 10–15 km. To reduce computational costs, it is therefore recommended to combine mesoscale nesting with [grid nesting](../Nesting/index.md). In such a setup, a parent domain with coarse spatial resolution covers a larger area, while higher resolved nested child domains focus on specific regions of interest.

**Always** ensure that the domains of interest are located sufficiently far from the lateral boundaries, and avoid performing any analysis within the flow adjustment zone. Under convective conditions, according to Kadasch et al. (2021), the required fetch length $d$ can be estimated by

$$
\begin{equation}
d = 2\,\frac{u_\textrm{h} * z_\textrm{i}}{w_{*}} \quad ,
\end{equation}
$$

with $u_\textrm{h}$ being the mean horizontal wind speed within the boundary layer, $z_\textrm{i}$ being the boundary-layer depth and $w_{*}$ the convective velocity.

Switching on the mesoscale nesting has several effects:

- Time-dependent boundary values of the $u-$, $v-$ and $w-$component of the wind velocity, potential temperature $\theta$, and the water vapor mixing ratio $q$ are read from the [dynamic driver file](../../../../Reference/LES_Model/IO-Files/Drivers/dynamic). These values are provided for all variables labeled as `ls_forcing_#bound>_#var`, where `#bound` is one of `left`, `right`, `north`, `south`, and `top`, and `#var` is the name of the respective variable (`u`, `v`, `w`, `pt`, etc.). The boundary values are linearly interpolated to the model time level at each LES timestep.  For a full list of input variables, see the [dynamic driver file](../../../../Reference/LES_Model/IO-Files/Drivers/dynamic) documentation.
- If the [chemistry](../Chemistry/index.md) and/or [SALSA](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/salsatoc) modules are enabled, corresponding boundary values for prognostic variables can also be supplied via the dynamic driver. If no chemistry values are present in the dynamic driver, boundary conditions fall back to the settings of [bc_cs_l](../../../../Reference/LES_Model/Namelists/#chemistry_parameters--bc_cs_l), [bc_cs_r](../../../../Reference/LES_Model/Namelists/#chemistry_parameters--bc_cs_r), [bc_cs_n](../../../../Reference/LES_Model/Namelists/#chemistry_parameters--bc_cs_n), [bc_cs_s](../../../../Reference/LES_Model/Namelists/#chemistry_parameters--bc_cs_s), and [bc_cs_t](../../../../Reference/LES_Model/Namelists/#chemistry_parameters--bc_cs_t). Similarly, in the absence of aerosol concentration values, the fallback parameters are [bc_aer_l](../../../../Reference/LES_Model/Namelists/#salsa_parameters--bc_aer_l), [bc_aer_r](../../../../Reference/LES_Model/Namelists/#salsa_parameters--bc_aer_r), [bc_aer_n](../../../../Reference/LES_Model/Namelists/#salsa_parameters--bc_aer_n), [bc_aer_s](../../../../Reference/LES_Model/Namelists/#salsa_parameters--bc_aer_n), and [bc_aer_t](../../../../Reference/LES_Model/Namelists/#salsa_parameters--bc_aer_t).
- If the synthetic turbulence generator is active, random perturbations are superimposed on the boundary values of the $u$-, $v$-, and $w$-components.
- Internally, a Dirichlet boundary condition is assumed for all of the above mentioned quantities at the lateral boundaries (see [bc_lr](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_lr), [bc_ns](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_ns)) and the top boundary (see also [bc_pt_t](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_pt_t) and [bc_q_t](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_q_t)), because fixed values are prescribed. As for the grid nesting, a zero-gradient Neumann boundary condition is assumed for the perturbation pressure at the lateral and top boundaries (please see also [bc_p_t](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_p_t)).
- As PALM is an incompressible model but mesoscale models often do not necessarily satisfy the divergence-free condition, a mass-flux correction is employed each timestep. A correction velocity is calculated from the integrated mass flux through the lateral model boundaries and imposed onto the $w-$component at the model-top boundary. For further details, please see Kadasch et al. (2021).

Boundary values at the lateral and top boundaries can be provided in two formats:

- **lod = 1** (level-of-detail): as a vertical profile (at lateral boundaries) and a single value (at the top boundary), or
- **lod = 2**: as spatially resolved cross-sections: *y-z* for the left and right, *x-z* for the north and south, and *x-y* for the top boundary.
In general, nesting in different models is possible. Multiple tools are available to convert mesoscale model data into a dynamic driver file: INIFOR (deprecated; targeted towards the COSMO model), PALM-METEO (targeted towards ICON, WRF and CAMX), and PROMET (targeted towards ICON, ICON-ART and WRF, can be easily extended towards other models).
For more information on the dynamic driver and how to create it, please see the [dynamic driver file](../../../../Reference/LES_Model/IO-Files/Drivers/dynamic) documentation.

## Basic Usage / Settings

An example of a simple setup using mesoscale nesting is presented below, followed by detailed explanations.

```Fortran
 &nesting_parameters
    domain_layouts = 'parent',  1,  -1,     5612,     0.0,     0.0,
                     'child1',  2,   1,     1350,  15856.0, 19504.0,
                     'child2',  3,   1,     1350,  15344.0, 13264.0,
    nesting_mode = 'one-way',
    nesting_datatransfer_mode  = 'mixed',

    homogeneous_initialization_child = .T.,
 /

 &initialization_parameters
    dx = 40.0,
    dy = 40.0,
    dz = 40.0,

    initializing_actions = 'read_from_file',
    ...
 /

 &radiation_parameters
    radiation_scheme = 'external',
    ...
 /

  &stg_par
    dt_stg_adjust    = 1800.0,
    dt_stg_call      =    0.0,
    compute_velocity_seeds_local = .F.,
 /


 &nesting_offl_parameters
 /
```
- The mesoscale nesting is activated by simply placing the namelist [&nesting_offl_parameters](../../../../Reference/LES_Model/Namelists/#offline-nesting-parameters) to the `_p3d` namelist file.
- The mesoscale nesting module must be activated **only** in the outermost root domain and **not** within any embedded child domains.
- Since the outermost parent domain should cover a larger area, relatively coarse grid spacings [dx](../../../../Reference/LES_Model/Namelists/#initialization_parameters--dx), [dy](../../../../Reference/LES_Model/Namelists/#initialization_parameters--dy), [dz](../../../../Reference/LES_Model/Namelists/#initialization_parameters--dz) are typically used. Please note that the recommended values apply to convective boundary layers and are not sufficient to resolve turbulence under neutral or stable conditions. In such cases, finer grid resolutions are required.
- Using the mesoscale nesting requires that the corresponding prognostic variables are also initialized with data from the [dynamic driver file](../../../../Reference/LES_Model/IO-Files/Drivers/dynamic). This is achieved by setting [initializing_actions](../../../../Reference/LES_Model/Namelists/#initialization_parameters--initializing_actions) = *'read_from_file'*. This case, either lod = 1 or lod = 2 data is read and applied to initialize the respective prognostic variables.
- The synthetic turbulence generator must be used in setups with mesoscale nesting. The generation of random synthetic turbulence is periodically adjusted, e.g., every 30 minutes, to the current atmospheric conditions. This adjustment is controlled by the parameter [dt_stg_adjust](../../../../Reference/LES_Model/Namelists/#stg_par--dt_stg_adjust). Dynamic adjustment is necessary because turbulent length and time scales vary throughout the diurnal cycle. Note, it is recommended to set [dt_stg_adjust](../../../../Reference/LES_Model/Namelists/#stg_par--dt_stg_adjust) to a value at least twice the typical eddy turnover time in the atmospheric boundary layer.
- The synthetic turbulence generator is invoked at each LES timestep, as specified by setting parameter [dt_stg_call](../../../../Reference/LES_Model/Namelists/#stg_par--dt_stg_call) = *0.0*. Using larger values for this parameter means the generator is no longer called at every timestep, which typically increases the turbulence adjustment fetch length.
- The calculation of random numbers used in the synthetic turbulence generator is distributed among multiple cores, i.e. they are not calculated locally on the core that is assigned to the inflow subdomain. Especially for large integral length scales (see documentation of the [synthetic turbulence generator](../Synthetic_turbulence/index.md)), this approach is significantly faster compared to local computations. However, with small integral length scales (e.g. in neutral or stable flows), local computation of the random numbers via setting [compute_velocity_seeds_local](../../../../Reference/LES_Model/Namelists/#stg_par--compute_velocity_seeds_local) = *.T.* is faster.
- In some setups, it may be beneficial to provide external radiation fluxes, e.g., such as during cloud evolution that reduces incoming solar radiation. In this case, setting the radiation namelist parameter [radiation_scheme](../../../../Reference/LES_Model/Namelists/#radiation_parameters--radiation_scheme) = *'external'* instructs the model to read incoming shortwave and longwave radiation from the dynamic driver file rather than computing the solar radiation internally.
This setting is optional and should be chosen based on the specific physical characteristics of the simulation.


## Limitations

The following limitations apply to the mesoscale nesting module:

- Mesoscale nesting cannot be used in combination with the [ocean mode](../Ocean/index.md).
- Mesoscale nesting is currently not realized for prognostic quantities used in the [bulk cloud model](../Cloud_Microphysics/bulk_microphysics.md).

**Note:** A mesoscale nesting using `lod = 2` boundary data in the dynamic driver generally does not make sense when small model domains with horizontal sizes of less than a few kilometers (or vertical sizes smaller than the boundary layer height of the respective szenario) are used, because the horizontal variation of data from the large-scale model is usually very small on such scales. Also, `lod = 2` data may not really be representative. For such cases it is recommended to provide just `lod = 1` vertical profile data in the dynamic driver.

## Further Information

For more detailed scientific and technical information about the mesoscale nesting interface, please refer to [Kadasch et al. (2021)](https://doi.org/10.5194/gmd-14-5435-2021).

## References

- **Kadasch, E., Sühring, M., Gronemeier, T., Raasch, S.** 2021: Mesoscale nesting interface of the PALM model system 6.0. Geoscientific Model Development, 14: 5435–5465,. [10.5194/gmd-14-5435-2021](https://doi.org/10.5194/gmd-14-5435-2021)
