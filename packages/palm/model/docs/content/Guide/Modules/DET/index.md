---
title: Overview
---
# Dust Emission and Transport

---

!!! warning
    This site is  Work in Progress.

## Theoretical Concept

### Introduction
With the Dust Emission and Transport (DET) module, PALM is able to simulate saltation-induced dust release and transport for different dust-sized particles. The implementation and physical background is based on the publications of [Klamt et al. (2024)](https://doi.org/10.1029/2023JD040058) and the AFWA scheme of [LeGrand et al. (2019)](https://doi.org/10.5194/gmd-12-131-2019).

**Attention:** In contrast to [Klamt et al. (2024)](https://doi.org/10.1029/2023JD040058), where the molecular viscosity of air has been calculated at each grid point and for each time step, a constant viscosity for an air temperature defined via namelist parameter [pt_surface](../../../../Reference/LES_Model/Namelists/#initialization_parameters--pt_surface) is used. This way, the performance of the DET is improved by one order of magnitude, while the effect on the results is only marginal. Furthermore, the calculated settling velocity is also based on [pt_surface](../../../../Reference/LES_Model/Namelists/#initialization_parameters--pt_surface), e.g. the height dependence of temperature is ignored.

All particles with a diameter less than $63\,\mu$m are interpreted as dust, i.e., clay ($\leq 4\,\mu$m) and silt ($4\,\mu$m $< D \leq 63\,\mu$m) particles belong to this group ([Klose, 2014](https://kups.ub.uni-koeln.de/5826/)). In general, the implementation follows an Eulerian approach, i.e., the simulation of dust transport is realized using a Cartesian 3D grid, and all calculations are performed at the grid points of this grid. The following dust-related processes are considered in the implementation: dust release from the (flat/smooth) ground due to saltation of particles, gravitational settling, dry deposition, large-scale subsidence, passive advection with the resolved turbulent flow, and subgrid-scale turbulent transport. Prognostic variables in DET are the mass concentrations (in kg m^-3^) of different dust size bins that can be defined by the user. 

#### Dust Emission
The parameterization of dust emission is based on the Air Force Weather Agency (AFWA) dust emission scheme of the GOCART aerosol model [LeGrand et al. (2019)](https://doi.org/10.5194/gmd-12-131-2019), which is part of WRF-Chem (the Weather Research and Forecasting (WRF) model coupled with Chemistry). In the AFWA scheme, dust emission is handled as a two-part process, wherein large particle saltation from coarser dust- and sand-sized particles is triggered by wind shear and leads to a fine-particle (dust-sized) bulk emission flux by saltation bombardment and aggregate disintegration. The bulk dust emission flux is then further distributed among different size bins. 

For estimating the vertical bulk dust emission flux, particle-size dependent horizontal saltation fluxes are calculated for each given saltation size bin. These saltation fluxes are then weighted depending on the soil configuration to form the so-called total vertically-integrated streamwise (horizontal) saltation flux. The default configuration of the saltation size bins is given in Table 1. The soil category 1 ("sand") of the STATSGO-FAO Database ([Pérez et al., 2011](https://doi.org/10.5194/acpd-11-17551-2011)) is used here. This category assumes 92% sand particles, 5% silt particles, and 3% clay particles. The mass weighting factors are the result of the product of these percentage values and the mass fractions given in Table 1.

| Saltation size bin          | 1    | 2      | 3      | 4      | 5      | 6     | 7     | 8     | 9     | 10     |
|-----------------------------|------|--------|--------|--------|--------|-------|-------|-------|-------|--------|
| Effective diameter ($\mu$m) | 1.42 | 8      | 20     | 32     | 44     | 70    | 130   | 200   | 620   | 1500   |
| Soil separate class         | Clay | Silt   | Silt   | Silt   | Silt   | Sand  | Sand  | Sand  | Sand  | Sand   |
| Particle density (kg m^-3^) | 2500 | 2650 | 2650 | 2650 | 2650 | 2650 | 2650 | 2650 | 2650 | 2650   |
| Mass fraction               | 1    | 0.25   | 0.25   | 0.25   | 0.25   | 0.0205| 0.0410| 0.0359| 0.3897| 0.5128 |
| Mass weighting factor       | 0.03 | 0.0125 | 0.0125 | 0.0125 | 0.0125 | 0.0189| 0.0377| 0.0330| 0.3585| 0.4718 |

_Table 1: Configuration of Saltation Size Bins and Associated Attributes for the AFWA Scheme Assuming the Soil Category 1 ("sand") of the STATSGO-FAO Database ([Pérez et al., 2011](https://doi.org/10.5194/acpd-11-17551-2011); [LeGrand et al. (2019)](https://doi.org/10.5194/gmd-12-131-2019))._

Once the vertical bulk dust emission flux is calculated, it is distributed among different dust size bins by using the so-called brittle fragmentation theory ([Kok, 2011](https://doi.org/10.1073/pnas.1014798108)). This theory enables the calculation of dust distribution weighting factors that are multiplied with the vertical bulk dust emission flux to form the bin-specific emission flux. This procedure requires information about the configuration of dust size bins, similar to the configuration of saltation size bins above. Table 2 shows the default configuration derived from [LeGrand et al. (2019)](https://doi.org/10.5194/gmd-12-131-2019).

| Dust size bin                 | 1     | 2     | 3     | 4     | 5     |
|-----------------------------------------------|-------|-------|-------|-------|-------|
| Lower bound diameter ($\mu$m) | 0.2   | 2     | 3.6   | 6     | 12    |
| Upper bound diameter ($\mu$m) | 2     | 3.6   | 6     | 12    | 20    |
| Effective diameter ($\mu$m)   | 1.46  | 2.8   | 4.8   | 9     | 16    |
| Particle density (kg m^-3^)   | 2500  | 2650  | 2650  | 2650  | 2650  |

_Table 2: Configuration of Dust Size Bins and Associated Attributes for the AFWA Scheme ([LeGrand et al. (2019)](https://doi.org/10.5194/gmd-12-131-2019))._

#### Gravitational Settling
The vertical divergence of the gravitational settling flux causes a change in the dust mass concentration which is why this process is considered in DET. The gravitational settling flux at a given position \((x,y,z)\) and time \(t\) is given in the Stoke's regime (laminar regime) by

\[
F_{\text{g,}j}(t,x,y,z) = - v_{\text{g,}j} \cdot c_{j}(t,x,y,z) \quad ,
\]

with the dust mass concentration \(c\) and the gravitational settling velocity \(v_\text{g}\), which is also known as the terminal velocity. This velocity is calculated based on [Jacobson (2005)](https://doi.org/10.1017/CBO9781139165389) (Eq. (20.4)) via:

\[
v_\text{g} = \frac{\bigl(\rho_{\text{d},j}-\rho_a\bigr) D_{\text{d},j}^2 g C_j}{18 \eta_\text{a}} \approx \frac{\rho_{\text{d},j} D_{\text{d},j}^2 g C_j}{18 \eta_\text{a}} \quad ,
\]

where \(\rho_{\text{a}}\) is the air density, \(D_{\text{d},j}\) is the dust particle effective diameter, \(\rho_{\text{d},j}\) is the dust particle density, \(\eta_\text{a}\) is the molecular (dynamic) viscosity of air, \(g\) is the gravitational acceleration, \(C_j\) is the Cunningham slip-flow correction for small particles, and the index \(j\) marks a specific dust bin. According to [Farrell and Sherman (2015)](https://doi.org/10.1177/0309133314562442), the assumption of the Stoke's regime (now commonly referred to as *Stoke's law*) is valid for particles with diameters less than $100\,\mu$m, i.e., the whole range of dust-sized particles. Thus, it is not recommended to use DET for particles other than dust.

#### Deposition
Instead of only considering the bin-specific dust emission flux, a net surface flux is eventually calculated which is the sum of \(F_{\text{e},j}\), the emission flux, and the deposition flux \(F_{\text{d},j}\):

\[
F_{\text{net},j} = F_{\text{e},j} + F_{\text{d},j} \quad .
\]

By definition, \(F_{\text{d},j}\) is negative at the surface (directed downwards). So far, wet deposition is not considered. The deposition flux is directly proportional to the local concentration of the dust size bin \(j\) (\(c_j\)) at a reference height (\(z_\text{r}\)) at which the bin-specific dry deposition velocity \(v_{\text{d},j}\) is evaluated:

\[
F_{\text{d},j} = - v_{\text{d},j} \cdot c_j \quad .
\]

The dry deposition velocity tries to capture all physical processes governing dry deposition. The reference height at which \(c_j\) and \(v_{\text{d},j}\) are evaluated is \(zu(1)=dz/2\) for default surfaces at the ground with \(dz\) being the vertical grid spacing. It is equal to the height of the constant flux layer where Monin–Obukhov similarity theory is assumed.

To calculate \(v_{\text{d},j}\) (being downward positive), we follow the dry deposition scheme Z01 of WRF-Chem ([Zeng et al., 2020](https://doi.org/10.5194/gmd-13-2125-2020)), which is based on electrical resistance models of [Slinn (1982)](https://doi.org/10.1016/0004-6981(82)90271-2) and [Zhang et al. (2001)](https://doi.org/10.1016/S1352-2310(00)00326-5). The dry deposition velocity is calculated for each dust size bin as

\[
v_\text{d} = \frac{1}{R_\text{a}+R_\text{s}+R_\text{a}R_\text{s}v_\text{g}}+v_\text{g} \quad .
\]

The overall transfer resistance is divided into an aerodynamic resistance above the canopy \(R_\text{a}\), accounting for the particle turbulent diffusion in the constant flux layer, and a surface resistance \(R_\text{s}\). The latter summarizes the quasi-laminar sublayer (or boundary) resistance and the canopy resistance (or bulk surface resistance). Brownian diffusion, impaction, and the rebound effect are considered here. Interception is neglected.

## Usage
PALM can be started using the following example setup:

```fortran
!-------------------------------------------------------------------------------
!-- INITIALIZATION PARAMETER NAMELIST
!-------------------------------------------------------------------------------
&initialization_parameters
!
!-- grid parameters
!-------------------------------------------------------------------------------
    nx                = 399, ! Number of gridboxes in x-direction (nx+1)
    ny                = 399, ! Number of gridboxes in y-direction (ny+1)
    nz                = 144, ! Number of gridboxes in z-direction (nz)

    dx                = 10.0, ! Size of single gridbox in x-direction
    dy                = 10.0, ! Size of single gridbox in y-direction
    dz                = 10.0, ! Size of single gridbox in z-direction

    dz_stretch_level  = 1200.0, ! Height (in m) where stretching starts
    dz_stretch_factor = 1.08,   ! dz(k+1) = dz(k) * dz_stretch_factor
!
!-- mode
!-------------------------------------------------------------------------------
    large_scale_subsidence       = .TRUE.,     ! Parameter to enable large-scale subsidence
    subs_vertical_gradient       = -0.0023, 0.0, ! Gradient(s) of the profile for the
                                                ! large-scale subsidence/ascent velocity
                                                ! (in (m/s) / 100 m)
    subs_vertical_gradient_level = 0.0, 1000.0, ! Height level from which on the gradient
                                                ! for the subsidence velocity is effective
                                                ! (in m, uv-grid)
!
!-- initialization
!-------------------------------------------------------------------------------
    initializing_actions       = 'set_constant_profiles', ! initial conditions

    ug_surface                 = 0.0, ! u-comp of geostrophic wind at surface
    vg_surface                 = 0.0, ! v-comp of geostrophic wind at surface

    pt_surface                 = 300.0, ! initial surface potential temp

    pt_vertical_gradient       =    0.0,
                                    2.0, ! piecewise temp gradients
    pt_vertical_gradient_level =    0.0,
                                  1000.0, ! height level of temp gradients

    reference_state            = 'horizontal_average', ! The instantaneous horizontally
                                                       ! averaged potential temperature
                                                       ! profile will be used as reference
                                                       ! state in the buoyancy term. Doing
                                                       ! so the pressure pertubation
                                                       ! will fluctuate around zero in a
                                                       ! xy-plane instead of an increasing
                                                       ! p* value of -25Pa/-50Pa/...
!
!-- boundary conditions
!-------------------------------------------------------------------------------
    bc_lr                      = 'cyclic',
    bc_ns                      = 'cyclic',

    bc_uv_b                    = 'dirichlet', ! (no-slip condition)
    bc_uv_t                    = 'neumann',   ! (free-slip condition)

    surface_heatflux           = 0.24,      ! sensible heat flux at the bottom surface
    bc_pt_b                    = 'neumann', ! required with surface_heatflux

    roughness_length           = 0.1,    ! Roughness length (in m)
    constant_flux_layer        = .TRUE., ! Parameter to switch on a constant flux layer
                                         ! at the bottom boundary.
!
!-- numerics
!-------------------------------------------------------------------------------
    rayleigh_damping_height    = 1300.0,  ! Well above inversion height
    rayleigh_damping_factor    = 0.1,     ! 0.0 <= rayleigh_damping_factor <= 1.0

    fft_method                 = 'fftw', ! http://www.fftw.org/
    momentum_advec             = 'ws-scheme', ! 5th order upwind scheme
    scalar_advec               = 'ws-scheme', ! 5th order upwind scheme

!
!-- physics
!-------------------------------------------------------------------------------
    omega     = 7.29212E-5, ! Angular velocity of the rotating system (Earth, in rad/s).
    latitude  = 52.37,      ! Geographical latitude (in degrees north) of Hannover

/ ! end of initialization parameter namelist

!-------------------------------------------------------------------------------
!-- RUNTIME PARAMETER NAMELIST
!-------------------------------------------------------------------------------
&runtime_parameters
!
!-- run steering
!-------------------------------------------------------------------------------
    end_time                   = 3600.0, ! simulation time of the 3D model

    npex = 16, npey = 16, ! Virtual processor topology.

    create_disturbances        = .TRUE., ! randomly perturbate horiz. velocity
    dt_disturb                 = 150.0,  ! interval for random perturbations
    disturbance_energy_limit   = 0.01,   ! upper limit for perturbation energy

    data_output_2d_on_each_pe  = .FALSE., ! don't do 2D output on each MPI rank
!
!-- data output
!-------------------------------------------------------------------------------
    netcdf_data_format         = 2, ! netCDF 64-bit-offset format

    dt_run_control             = 0.0,   ! output interval for run control
    dt_data_output             = 900.0, ! output interval for general data
    dt_data_output_av          = 900.0, ! output interval for averaged data
    dt_dopr                    = 900.0, ! output interval for profile data

    nz_do3d                    = 60, ! 3d output up to 600m

    data_output                = 'u_xy', 'u_xz', 'u_xy_av', 'u_xz_av',
                                 'v_xy', 'v_xz', 'v_xy_av', 'v_xz_av',
                                 'w_xy', 'w_xz', 'w_xy_av', 'w_xz_av',
                                 'p_xy', 'p_xz', 'p_xy_av', 'p_xz_av',
                                 'theta_xy', 'theta_xz',
                                 'wspeed', 'wdir',
                                 'dust_mc_bin1', 'dust_mc_bin2', 'dust_mc_bin3',
                                 'dust_mc_bin4', 'dust_mc_bin5', 
                                 'dust_mc_bin2_xy', 'dust_mc_bin2_xz', 'dust_mc_bin2_yz', 
                                 'dust_mc_bin2_av',
                                 'dust_mc_bin3_xy', 'dust_mc_bin3_xz', 'dust_mc_bin3_yz',
                                 'dust_mc_bin3_av',
                                 'dust_emis_flux*_bin1_xy', 'dust_emis_flux*_bin2_xy',
                                 'dust_emis_flux*_bin3_xy', 'dust_emis_flux*_bin4_xy', 
                                 'dust_emis_flux*_bin5_xy',
                                 'dust_emis_flux*_bin1_xy_av',
                                 'dust_depo_flux*_bin1_xy', 'dust_depo_flux*_bin2_xy',
                                 'dust_depo_flux*_bin3_xy', 'dust_depo_flux*_bin4_xy', 
                                 'dust_depo_flux*_bin5_xy',
                                 'dust_depo_flux*_bin1_xy_av',
                                 'clay', 'clay_av', 'clay_xy', 'clay_xy_av', 
                                 'dust', 'dust_av', 'dust_xz', 'dust_xz_av', 
                                 'silt', 'silt_av', 'silt_yz', 'silt_yz_av', 

    data_output_pr             = '#u', '#v', 'w', 'w_subs',
                                 'wtheta', 'w"theta"', 'w*theta*',
                                 'e', 'e*',
                                 'u*2', 'v*2', 'w*2', 'theta*2',
                                 '#theta', '#km', '#kh',
                                 'hyp', 'p', 'rho',
                                 'clay', 'silt', 'dust',

    section_xy                 = 1, 5, 10, ! grid index for 2D XY cross sections
    section_xz                 = 200, -1,  ! grid index for 2D XZ cross sections
    section_yz                 = 200, -1,  ! grid index for 2D YZ cross sections

    averaging_interval         = 900.0, ! averaging interval general data
    dt_averaging_input         = 10.0,  ! averaging general data sampling rate

    averaging_interval_pr      = 900.0, ! averaging interval profile data
    dt_averaging_input_pr      = 10.0,  ! averaging profile data sampling rate

/ ! end of runtime parameter namelist

!-------------------------------------------------------------------------------
!-- DUST EMISSION AND TRANSPORT PARAMETER NAMELIST
!-------------------------------------------------------------------------------
&det_parameters
!
!-- general parameters
!-------------------------------------------------------------------------------
    deposition_scheme = 'Z01',   ! scheme for calculating deposition
    det_start_time    = 0.0,     ! simulation time after det is active
    switch_off_module = .FALSE., ! switch of (de)activating det
!
!-- parameters for dry deposition
!-------------------------------------------------------------------------------
    alpha_imp                      = 50.0, ! land use dependent parameter for impaction efficiency calculation
    brownian_diffusion_coefficient = 0.54, ! parameter for Brownian diffusion
!
!-- boundary conditions
!-------------------------------------------------------------------------------
    bc_dm_b   = 'neumann', ! bottom boundary condition for prognostic det variables
    bc_dm_l   = 'cyclic',  ! west/left boundary condition for prognostic det variables
    bc_dm_n   = 'cyclic',  ! north boundary condition for prognostic det variables
    bc_dm_r   = 'cyclic',  ! east/right boundary condition for prognostic det variables
    bc_dm_s   = 'cyclic',  ! south boundary condition for prognostic det variables
    bc_dm_t   = 'neumann', ! top boundary condition for prognostic det variables
!
!-- saltation bin parameters
!-------------------------------------------------------------------------------
    bin_mass_fraction_ssc      = 1.0, 0.25, 0.25, 0.25, 0.25, 0.0205, 0.0410,
                                 0.0359, 0.3897, 0.5128,         ! bin-specific mass fraction of corresponding soil separate class
    diameter_saltation         = 1.42E-6, 8.0E-6, 20.0E-6, 32.0E-6, 44.0E-6,
                                 70.0E-6, 130.0E-6, 200.0E-6, 620.0E-6,
                                 1500.0E-6,                      ! effective diameter of a saltation size bin
    mass_fraction_ssc          = 0.03, 0.05, 0.05, 0.05, 0.05, 0.92, 0.92, 0.92,
                                 0.92, 0.92,                     ! mass fraction of soil separate class
    n_saltation_bins           = 10,                             ! number of saltation size bins
    particle_density_saltation = 2500.0, 2650.0, 2650.0, 2650.0, 2650.0, 2650.0,
                                 2650.0, 2650.0, 2650.0, 2650.0, ! particle density of a saltation size bin
!
!-- dust bin parameters
!-------------------------------------------------------------------------------
    diameter_dust         = 1.46E-6, 2.8E-6, 4.8E-6, 9.0E-6, 16.0E-6, ! effective diameter of a dust size bin
    lower_bound_diameter  = 0.2E-6, 2.0E-6, 3.6E-6, 6.0E-6, 12.0E-6,  ! minimum effective diameters represented by the dust size bin
    n_dust_bins           = 5,                                        ! number of dust size bins
    particle_density_dust = 2500.0, 2650.0, 2650.0, 2650.0, 2650.0,   ! particle density of a dust size bin
    upper_bound_diameter  = 2.0E-6, 3.6E-6, 6.0E-6, 12.0E-6, 20.0E-6, ! maximum effective diameters represented by the dust size bin

/ ! end of det parameter namelist
```

**Note:** DET can also be used with an empty namelist because the default configuration sets every namelist parameter appropriately. Output of time series of horizontally averaged bulk deposition flux and bulk emission flux (where bulk means the integral over all size bins), as well as maximum values of the bulk emission flux and minimum values of the bulk deposition flux, are added by default to the standard time series output.

## Notes, shortcommings and open issues

1. In contrast to [Klamt et al. (2024)](https://doi.org/10.1029/2023JD040058), where the molecular viscosity of air has been calculated at each grid point and for each time step, a constant viscosity for an air temperature defined via namelist parameter [pt_surface](../../../../Reference/LES_Model/Namelists/#initialization_parameters--pt_surface) is used. Furthermore, the calculated settling velocity is also based on [pt_surface](../../../../Reference/LES_Model/Namelists/#initialization_parameters--pt_surface), e.g. the height dependence of temperature is ignored.

2. It is not recommended to simulate sand particles ($63\,\mu$m $< D \leq 2\,$mm) with the provided implementation, especially not the larger fraction ($\gg 100\,\mu$m), because, on the one hand, the implemented emission algorithm is specially designed and tested for dust-sized particles and on the other hand, the gravitational settling is calculated from *Stoke's law*, which can only be regarded as valid for particle diameters $< 100\,\mu$m.

3. DET is realized for horizontal upward-facing default surfaces only. Land and urban surfaces are not considered. All soil parameters given in DET are assumed to be horizontally homogeneous.

4. You may use the LSM together with the DET, but other LSM surface types than flat sandy surfaces should be simulated with caution, especially rougher surfaces with plants. Also, heterogeneous LSM surfaces may not work appropriately. If you plan such simulations, check every setting of the parameters carefully. The equations of the deposition scheme might also required to be adjusted in such cases.

5. DET can not be used with topography.

6. The effects of soil moisture and wet deposition are not considered.

7. The configuration of the saltation size bins in Table 1 must start with at least one saltation size bin that is assigned to the soil separate class clay. Otherwise, the sandblasting efficiency, which connects the total horizontal saltation flux with the vertical bulk dust flux, can not be calculated ([Marticorena and Bergametti, 1995](https://doi.org/10.1029/95JD00690)).

8. Other, probably more accurate deposition schemes, should be implemented. [Zhang and Shao (2014)](https://doi.org/10.5194/acp-14-12429-2014) and [Zhang and He (2014)](https://doi.org/10.5194/acp-14-3729-2014) are two promising options. They were evaluated by [Bergametti et al., 2018](https://doi.org/10.1029/2018JD028964) and [Khan and Perlinger, 2017](https://doi.org/10.5194/gmd-10-3861-2017).  

9. So far not implemented:
    - horizontal boundary conditions for the prognostic quantities of DET that are not cyclic,
    - masked output.
