---
title: Pressure solver usage
---

# Pressure solver usage
<br>

The pressure solver removes any divergence from the velocity field. This is required in models that use the Boussinesq or anelastic approximation. Moreover, numerical schemes for advection are often based on flux-form formulations. Numerical errors in these schemes increase if divergence is only poorly removed from the velocity field; scalar conservation errors are a typical result.

In general, PALM uses a predictor–corrector method to remove the divergence. In the predictor step, the momentum equation is advanced without the dynamic pressure gradient term. The divergence is computed from the provisional velocity field, and solving a Poisson equation for pressure with the divergence on the right-hand side yields the so-called perturbation pressure. This is done by the pressure solver. In the corrector step, the final velocity field is obtained by performing another time step based on the provisional velocity field and the pressure gradient term.

The computational time for solving the Poisson equation and removing the divergence can exceed 70% for larger grids.

In the following, the term “total divergence” is defined as the sum over all grid points of the magnitude of the local divergence. The quality of the solver is defined by the ratio of the total divergences of the provisional and final velocity fields. The total divergence for the provisional and final velocity field (i.e. before and after the pressure solver has been applied) is displayed in the time step output of the run-control file, columns `DIVOLD` and `DIVNEW`.

```
Run-control output:
------------------

RUN  ITER.  HHH:MM:SS.SS    DT(E)     UMAX     VMAX     WMAX   ...   DIVOLD     DIVNEW   ...  MGCYC
-------------------------------------------------------------- ... --------------------- ... ------
  0      0  000:00:00.00  16.2000A   0.2544D -0.2772D  0.1603  ...  0.708E-02  0.101E-09 ...      2
  0      1  000:00:16.19  10.5000D  -0.2471  -0.2634   0.1500  ...  0.286E-03  0.122E-11 ...      2
```



## Available pressure solvers in PALM

PALM offers two solvers:

- A direct solver based on 2D fast Fourier transform (FFT) of the Poisson equation (hereafter called the FFT solver). It is switched on via [psolver](../../../Reference/LES_Model/Namelists/#initialization_parameters--topography) = *'poisfft'* (default setting).

- An iterative solver based on a multigrid method, using Gauss–Seidel iteration with red/black ordering on each grid level (MG solver), and switched on via [psolver](../../../Reference/LES_Model/Namelists/#initialization_parameters--topography) = *'multigrid'*.

The solver features vary based on the setting of parameter [topography](../../../Reference/LES_Model/Namelists/#initialization_parameters--topography).


## Solver features for setups using flat terrain 

Solver features for [topography](../../../Reference/LES_Model/Namelists/#initialization_parameters--topography) = *'flat'* are:

**FFT-solver:**

- The FFT solver typically reduces the total divergence (final velocity divergence relative to the provisional velocity) by more than 10 orders of magnitude.

- It is very fast in CPU terms but does not scale well beyond about 10,000 cores.

- It requires that the number of grid points along `x` and `y` can be factored into primes up to 13.

**MG-solver:**

- For the MG solver, the reduction of divergence depends on several factors. Better reduction often requires more CPU time. Several parameters can be used to control solver efficiency.

- With reasonable effort, total divergence can be reduced by about 8 orders of magnitude.

**For flat terrain without buildings, using the FFT solver is recommended. For runs with any topography or obstacles present, the MG-solver is recommended.**
    

## Solver features for setups using orography/buildings

**FFT-Solver:**

- The FFT solver always computes the pressure solution for all grid points, including those below the orography surface and inside buildings. At these points, divergence is masked to zero before calling the solver (the masking method). Boundary conditions for pressure at surfaces are not considered (except at the bottom and top of the domain). Applying the solver this way generates nonzero velocities at surfaces and within buildings. Wall-normal velocities at surface grid points are set to zero after the pressure correction to enforce impermeable walls.

**MG-Solver:**

- By default, the MG solver calculates the pressure for fluid grid points only and explicitly applies boundary conditions for the perturbation pressure at walls. Typically, Neumann (zero-gradient) conditions are used to ensure that velocities at walls do not change due to pressure gradients. Wall-normal velocities are zero and should remain zero after the solver.

- By setting [masking_method](../../../Reference/LES_Model/Namelists/#initialization_parameters--masking_method) = *.T.*, the MG solver can be forced to use the masking method.


### Disadvantage of the masking method

Pressure correction with masking method can change velocities at the walls, and divergence is properly corrected only when using those velocities. Since velocities are reset to zero at walls after the pressure solver has been applied, the divergence reduction near walls is less effective; often only one or two orders of magnitude. This can severely affect scalar conservation near surfaces.

**For this reason, for setups with orography/buildings it is generally recommended to use the MG-solver without masking.**


## Required / recommended setups and settings for the solvers

**FFT-solver:**

- The number of grid points along x and y must be composed of prime factors. Allowed primes depend on the FFT method used. Details are given in the [fft_method](../../../Reference/LES_Model/Namelists/#initialization_parameters--fft_method) parameter description. The most flexible option is [fft_method](../../../Reference/LES_Model/Namelists/#initialization_parameters--fft_method) = *'fftw'*. Note that the prime-factor criterion depends on the horizontal boundary conditions: for cyclic boundaries ([bc_lr](../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_lr) = [bc_ns](../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_ns) = *'cyclic'*), [nx](../../../Reference/LES_Model/Namelists/#initialization_parameters--nx)+1 and [ny](../../../Reference/LES_Model/Namelists/#initialization_parameters--ny)+1 are used; for non-cyclic boundaries, [nx](../../../Reference/LES_Model/Namelists/#initialization_parameters--nx) and [ny](../../../Reference/LES_Model/Namelists/#initialization_parameters--ny) are used. The FFT solver does not allow mixed boundaries (cyclic in one direction and non-cyclic in the other).

**MG-solver:**

- Subdomains must be equally sized (i.e., non-uniform subdomains are not allowed), and the number of subdomain grid points along `x`, `y`, and `z` should be a multiple of 2, e.g., `nx`(subdomain) = *40* = *2* × *2* × *2* × *5*, which allows the use of four multigrid levels (provided the same holds for `ny` and `nz`). Fewer than three levels may significantly compromise divergence reduction on larger spatial scales. The subdomain size can be calculated manually based on the total domain grid points and the virtual processor grid (npex and npey), but it is also printed in the header of the run-control (`_rc`) file:
```
 Domain size:       x =   1280.000 m     y =   1280.000 m  z(u) =   2562.500 m

 Number of gridpoints (x,y,z):  (0: 255, 0: 255, 0: 513)
 Subdomain size (x,y,z):        (    32,     32,    514)
```
The above example would allow 5 levels.

- Various parameters are available to steer the multigrid solver (see [ngsrb](../../../Reference/LES_Model/Namelists/#initialization_parameters--ngsrb), [ngsrb_initial](../../../Reference/LES_Model/Namelists/#initialization_parameters--ngsrb_initial), [ngsrb_initial_timesteps](../../../Reference/LES_Model/Namelists/#initialization_parameters--ngsrb_initial_timesteps), [mg_cycles](../../../Reference/LES_Model/Namelists/#initialization_parameters--mg_cycles), [cycle_mg](../../../Reference/LES_Model/Namelists/#initialization_parameters--cycle_mg), [residual_limit](../../../Reference/LES_Model/Namelists/#initialization_parameters--residual_limit), [mg_switch_to_pe0_level](../../../Reference/LES_Model/Namelists/#initialization_parameters--mg_switch_to_pe0_level), [max_mg_grid_levels](../../../Reference/LES_Model/Namelists/#initialization_parameters--max_mg_grid_levels)). The default settings of these parameters should serve as a good first choice for many setups. Depending on the user’s requirements (better divergence reduction vs. less CPU time), the settings should be optimized by conducting test runs. Settings that are used in a run are listed in the header section of the run-control (`_rc`) file:
```
 --> Solve perturbation pressure via multigrid_noopt method (w-cycle)
     number of grid levels:                     5
     Gauss-Seidel red/black iterations:         2
                      at initial start:       100   until timestep   0
     fixed number of multigrid cycles:          2
```

## Special settings for the MG-solver

- By default, the default number of Gauss-Seidel iterations as given by [ngsrb](../../../Reference/LES_Model/Namelists/#initialization_parameters--ngsrb) is significantly increased for the case that the solver is called for the first time in an initial run. The number of iterations to be used then can be set via [ngsrb_initial](../../../Reference/LES_Model/Namelists/#initialization_parameters--ngsrb_initial). This setting often makes sense, e.g. in case of setups with topography, because the initial wind field is not adjusted to the topography, and therefore can generate huge divergence. More iterations help to further reduce the divergence. It may even be required to use more iterations for the first time steps, too, which can be steered via [ngsrb_initial_timesteps](../../../Reference/LES_Model/Namelists/#initialization_parameters--ngsrb_initial_timesteps). Good evidence that more iterations are required is, if a jump in the velocity magnitude is observed in the run-control file after the first time step.

- The number of so-called multigrid-cylces (see [cycle_mg](../../../Reference/LES_Model/Namelists/#initialization_parameters--cycle_mg) and [mg_cycles](../../../Reference/LES_Model/Namelists/#initialization_parameters--mg_cycles)) to reduce the divergence by a given order of magnitude often varies throughout the simulation. With a fixed number of cycles (which is the default), the reduction of divergence can vary from time step to time step, and may sometimes be so poor, that it badly affects the scalar conservation. In such cases, parameter [residual_limit](../../../Reference/LES_Model/Namelists/#initialization_parameters--residual_limit) can be set to guarantee, that the total divergence is always reduced by a certain order of magnitude as given by the parameter value. It further requires to set [mg_cycles](../../../Reference/LES_Model/Namelists/#initialization_parameters--mg_cycles) = *-1*. Note that with this setting the CPU requirement of a run is difficult to estimate, because the number of required cycles are not known in advance. The actual number of MG cycles for each time step is listed in the run-control file in the last column (`MGCYC`)of the time step output (see listing further above).


## Special optimized FFT-solver

An optimized FFT-solver can be selected via [psolver](../../../Reference/LES_Model/Namelists/#initialization_parameters--psolver) = *'poisfft_sm'*. It uses a 1d-domain-decomposition for pressure only, reducing the number of transpositions required between the forward/backward FFTs and the tridiagonal solver from 6 to 2, which for larger setups may improve the overall model performance by more than 20-30%.


## Pressure output

The pressure that is output via [data_output](../../../Reference/LES_Model/Namelists/#runtime_parameters--data_output) = *'p'* is the dynamic (perturbation) pressure calculated by the pressure solver. As in every incompressible or anelastic CFD model, PALM does not use (and provide) static pressure.

The main task of the pressure solver in incompressible or anelastic CFD codes is to make the velocity field free of divergence. This does not work in case of the specific solvers used in PALM, if Neumann boundary conditions for pressure are set both at the bottom and the top of the domain (see [bc_p_b](../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_p_b) and [bc_p_t](../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_p_t)). Under such conditions, any mean (horizontally averaged) vertical velocity has to be removed from the velocity field before the solver is called. This changes the pressure calculated by the solver, because it removes the effect that buoyancy should have on the pressure. If buoyancy effects should be accounted for in the pressure output, parameters [bc_p_t](../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_p_t) = *'dirichlet'* and [reference_state](../../../Reference/LES_Model/Namelists/#initialization_parameters--reference_state) = *'horizontal_average'* should be used.

Furthermore, PALM's pressure output does not contain the effect of the isotropic term of the stress tensor (internally, the so-called deviatoric kinemtaic stress tensor is used in the momentum equation). To get the correct dynamic pressure, 2/3 of the SGS turbulent kinetic energy `e` has to be added to the pressure in a post-processing step. The SGS-TKE can be obtained via [data_output](../../../Reference/LES_Model/Namelists/#runtime_parameters--data_output) = *'e'*.
