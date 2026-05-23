---
title: Overview
---
# Sythentic Turbulence Generator
---

## Purpose

The synthetic turbulence generator (STG) offers the possibility to generate a turbulent inflow condition in non-cyclic setups. The STG can be useful in idealized setups with fixed inflow profiles at the left model boundary ([bc_lr](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_lr) =*'dirichlet/radiation'*) in case where no other [turbulent inflow method](../../../../Reference/LES_Model/Namelists/#turbulent_inflow_parameters--turbulent_inflow_method) is applicable. Moreover, the STG is required in [mesoscale-nested](../Mesoscale_nesting/index.md) simulations, where time-dependent boundary conditions from mesoscale model output are imposed onto the lateral and top domain boundaries.

## General Information


The mesoscale nesting is activated by adding the namelist [`&stg_par`](../../../../Reference/LES_Model/Namelists/#synthetic-turbulence-generator-parameters) to the `_p3d` namelist file.
Due to the nature of mesoscale RANS models, turbulence is parameterized and thus the boundary values are free of any turbulent fluctuations. Mirocha et al. (2014) showed that without adding perturbations the turbulent flow needs several tens of kilometres to sufficiently develop. In order to accelerate the spatial development of turbulence in [mesoscale-nested](../Mesoscale_nesting/index.md) setups, the STG approach by Xie and Castro (2008) is employed, where perturbations are added onto the three velocity components imposed at the lateral boundaries.

It is **noted** that the STG imposed spatially and temporarily correlated random numbers which do not represent a full turbulence spectrum. Even if the turbulence profiles match with the actual conditions, the LES needs to generate a full turbulence spectrum. The STG can be used to accelerate this process but it still features a certain adjustment fetch, which **always** needs to be checked.

## Basic Usage / Settings

In the following, we refer to example setups using the STG. For an example of a mesoscale nested setup where the turbulence profiles for the synthetic turbulence generation are parametrized, we refer to the [mesoscale nesting](../Mesoscale_nesting/index.md) documentation, where also parameters for the STG are explained.
An example of a simple setup using the STG in case of idealized turbulent inflows is given in the following:

```Fortran
 &initialization_parameters
    dx = 4.0,
    dy = 4.0,
    dz = 4.0,

    bc_lr = 'dirichlet/radiation',
    bc_nr = 'cyclic',
    ...
 /


  &stg_par
    dt_stg_call = 0.0,
 /

 &turbulent_inflow_parameters
   switch_off_module = .T.,
/

```

- In the case of an inflow boundary condition on the left model domain boundary (given by [bc_lr](../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_lr) =*'dirichlet/radiation'*), the STG is simply activated by adding the namelist [`&stg_par`](../../../../Reference/LES_Model/Namelists/#synthetic-turbulence-generator-parameters) to the `_p3d` namelist file. In this example, the STG is invoked every timestep, as specified by the parameter [dt_stg_call](../../../../Reference/LES_Model/Namelists/#stg_par--dt_stg_call), which is the default. Using larger values for this parameter means the generator is no longer called at every timestep, which typically increases the turbulence adjustment fetch length.
- The STG is only applicable in case no other turbulent inflow boundary conditions is chosen, i.e. the [turbulent inflow module](../../../../Reference/LES_Model/Namelists/#turbulent-inflow-parameters) needs to be switched-off.
- In case of a turbulent inflow boundary condition, the Reynolds stress as well as the length and time scales are parametrized by default (see [reference section](../../../../Reference/LES_Model/Modules/Synthetic_turbulence/index.md)). In case of available turbulence information, as e.g., inferred from precursor simulations, these can be written to an ASCII file with the suffix `_stg`. If this file is present in the `INPUT` folder, turbulence data will be taken from this file rather than being parametrized. An example of the `_stg` file can be download [here](STG_PROFILES).
- Optionally, temperature perturbations can be imposed by setting [disturb_theta](../../../../Reference/LES_Model/Namelists/#stg_par--disturb_theta) = *.T.* and prescribing a time-constant maximum temperature amplitude by [theta_amplitude](../../../../Reference/LES_Model/Namelists/#stg_par--theta_amplitude). Temperature perturbations correlate with the perturbations added onto the w-component, i.e. they feature the same length and time-scales. Please note, this is still an experimental feature which will be improved in the future.

## Limitations

The STG cannot be used in combination with another turbulent inflow boundary condition as given by [turbulent inflow module](../../../../Reference/LES_Model/Namelists/#turbulent-inflow-parameters).

At this point we emphasize that using the implemented STG from Xie and Castro (2008) only generates turbulence which is correlated in space and time but not necessarily generate realistic turbulent structures. Large coherent structures like e.g. hexagonal pattern as typically observed in a convective boundary layer, however, cannot be generated by this method. Further, we want to add that turbulence is only added to the three wind components, and optionally to the potential temperature as an experimental model feature. In the latter case, however, the imposed temperature amplitude is fixed to a time-constant value and is not yet scaled according to the actual atmospheric conditions. Furthermore, no perturbations are added to the subgrid-scale turbulent-kinetic energy.

## Reference

For more detailed scientific and technical information about the STG see the [reference section](../../../../Reference/LES_Model/Modules/Synthetic_turbulence/index.md) as well as [Kadasch et al. (2021)](https://doi.org/10.5194/gmd-14-5435-2021).

## References

- **Kadasch, E., Sühring, M., Gronemeier, T., Raasch, S.** 2021: Mesoscale nesting interface of the PALM model system 6.0. Geoscientific Model Development, 14: 5435–5465,. [10.5194/gmd-14-5435-2021](https://doi.org/10.5194/gmd-14-5435-2021)

- **Mirocha, J., Kosović, B., Kirkil, G.** 2014: Resolved Turbulence Characteristics in Large-Eddy Simulations Nested within Mesoscale Simulations Using the Weather Research and Forecasting Model, Mon. Weather Rev., 142, 806–831, [doi.org/10.1175/MWR-D-13-00064.1](https://doi.org/10.1175/MWR-D-13-00064.1)

- **Xie, Z. and Castro, I.** 2008: Efficient Generation of Inflow Conditions for Large Eddy Simulation of Street-Scale Flows, Flow Turbul. Combust., 81, 449–470, [doi.org/10.1007/s10494-008-9151-5](https://doi.org/10.1007/s10494-008-9151-5)
