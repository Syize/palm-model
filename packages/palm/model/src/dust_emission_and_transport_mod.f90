!> @file dust_emission_and_transport.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 2024 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Sebastian Giersch (Leibniz Universitaet Hannover)
!
!
! Description:
! ------------
!> Calculation of dust-sized particle emission and transport over desert-like (sandy) surfaces. The
!> emission is caused by saltating particles. Transporting mechanisms are gravitational settling,
!> dry deposition, passive advection with the resolved-scale flow, and subgrid-scale turbulent
!> transport. Only default and natural type (upward-facing) surfaces are able to release dust.
!>
!> This routine is based on code developed by Jannis Klamt as part of his master thesis, and which
!> is published as Klamt, J., S. Giersch, and S. Raasch (2024): Geophys. Res. Atmos., 129(7),
!> DOI: 10.1029/2023JD040058.
!>
!> ATTENTION:
!> Deviating from that paper, this routine is using a constant (dynamic) air viscosity
!> at a temperature given by pt_surface and pressure of 1000 hPa), and instead of the air
!> temperature at a specific height, the potential temperature given via pt_surface is always used.
!> Results depend only very marginal on these changes (0.3%), but the code performance increases
!> drastically, since e.g. the Exner function isn't required.
!>
!> @note Functionality must be tested in detail
!> @todo Consider humidity
!> @todo Enable horizontal (lateral) boundary conditions different than cyclic and that are 
!> independent of boundary conditions for other scalars, similar to SALSA or the chemistry model
!> @todo Remove basic constants the mass of one air molecule to basic_constants_and_equations_mod
!> @todo Unify Cunningham correction factor and dynamic viscosity of air with salsa and include
!> equations into basic_constants_and_equations_mod
!> @todo Unify calculation of the aerodynamic resistance in the whole code
!> @todo Enable to set initial values for dust mass concentration
!> @todo Compute masked output
!> @todo Add alternative deposition scheme of Zhang and Shao (2014)
!> @todo Introduce debug output
!--------------------------------------------------------------------------------------------------!
 MODULE dust_emission_and_transport_mod

    USE arrays_3d,                                                                                 &
        ONLY:  pt,                                                                                 &
               rho_air,                                                                            &
               rho_air_zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  k_boltzmann,                                                                        &
               g,                                                                                  &
               pi

    USE control_parameters,                                                                        &
        ONLY:  cyclic_fill_initialization,                                                         &
               debug_output,                                                                       &
               det_enabled,                                                                        &
               max_pr_cs,                                                                          &
               max_pr_det,                                                                         &
               message_string,                                                                     &
               pt_surface,                                                                         &
               restart_data_format_output,                                                         &
               time_since_reference_point

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rrd_mpi_io,                                                                         &
               rd_mpi_io_check_array,                                                              &
               wrd_mpi_io

    USE surface_mod,                                                                               &
        ONLY:  surf_type,                                                                          &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_top,                                                                           &
               surf_usm


    IMPLICIT NONE

    CHARACTER(LEN=4) ::  deposition_scheme = 'Z01' !< namelist parameter

    CHARACTER(LEN=20) ::  bc_dm_b = 'neumann'    !< namelist parameter
    CHARACTER(LEN=20) ::  bc_dm_l = 'undefined'  !< namelist parameter
    CHARACTER(LEN=20) ::  bc_dm_n = 'undefined'  !< namelist parameter
    CHARACTER(LEN=20) ::  bc_dm_r = 'undefined'  !< namelist parameter
    CHARACTER(LEN=20) ::  bc_dm_s = 'undefined'  !< namelist parameter
    CHARACTER(LEN=20) ::  bc_dm_t = 'neumann'    !< namelist parameter

    INTEGER(iwp) ::  det_pr_count = 0          !< counter for det profiles
    INTEGER(iwp) ::  dots_num_det = 0          !< number of det time series
    INTEGER(iwp) ::  dots_start_index_det = 0  !< start index for time series of this module
    INTEGER(iwp) ::  ibc_dm_b                  !< index for the bottom boundary condition
    INTEGER(iwp) ::  ibc_dm_t                  !< index for the top boundary condition
    INTEGER(iwp) ::  n_dust_bins = 5           !< namelist parameter
    INTEGER(iwp) ::  n_saltation_bins = 10     !< namelist parameter

    INTEGER(iwp), PARAMETER ::  n_dust_bins_max = 10       !< maximum number of dust size bins
    INTEGER(iwp), PARAMETER ::  n_saltation_bins_max = 20  !< maximum number of saltation size bins

    INTEGER(iwp), DIMENSION(99) ::  det_pr_index = 0   !< index for det profiles

    LOGICAL ::  clay_calculated = .FALSE.  !< flag to indicated if diagnostic quantity has been calculated
    LOGICAL ::  silt_calculated = .FALSE.  !< flag to indicated if diagnostic quantity has been calculated
    LOGICAL ::  dust_calculated = .FALSE.  !< flag to indicated if diagnostic quantity has been calculated

    REAL(wp) ::  air_viscosity                             !< dynamic air viscosity in kg/(ms)
    REAL(wp) ::  alpha_imp = 50.0_wp                       !< namelist parameter
    REAL(wp) ::  alpha_s                                   !< sandblasting efficiency
    REAL(wp) ::  brownian_diffusion_coefficient = 0.54_wp  !< namelist parameter
    REAL(wp) ::  det_start_time = 0.0_wp                   !< namelist parameter

    REAL(wp), PARAMETER ::  am_airmol = 4.8096E-26_wp   !< average mass of an air molecule

    REAL(wp), PARAMETER ::  not_set = HUGE( 1.0_wp )  !< fill value to identify arrays that were not set by the user

    REAL(wp), DIMENSION(1:n_saltation_bins_max) ::  bin_mass_fraction_ssc = not_set       !< bin-specific mass fraction of corresponding soil separate class, namelist parameter
    REAL(wp), DIMENSION(1:n_saltation_bins_max) ::  diameter_saltation = not_set          !< effective diameter of a saltation size bin, namelist parameter
    REAL(wp), DIMENSION(1:n_saltation_bins_max) ::  mass_fraction_ssc = not_set           !< mass fraction of soil separate class, namelist parameter
    REAL(wp), DIMENSION(1:n_saltation_bins_max) ::  particle_density_saltation = not_set  !< particle density of a saltation size bin, namelist parameter

    REAL(wp), DIMENSION(1:n_dust_bins_max) ::  lower_bound_diameter = not_set   !< minimum effective diameters represented by the dust size bin, namelist parameter
    REAL(wp), DIMENSION(1:n_dust_bins_max) ::  upper_bound_diameter = not_set   !< maximum effective diameters represented by the dust size bin, namelist parameter
    REAL(wp), DIMENSION(1:n_dust_bins_max) ::  diameter_dust = not_set          !< effective diameter of a dust size bin, namelist parameter
    REAL(wp), DIMENSION(1:n_dust_bins_max) ::  particle_density_dust = not_set  !< particle density of a dust size bin, namelist parameter

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dm_rel  !< relative saltation bin-specific mass distribution of particles in the surface soil
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ds_rel  !< saltation bin-specific wheighting factors for the saltation flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ds_sfc  !< saltation bin-specific basal surface coverage fractions
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dv      !< normalized volume distributions of each dust size bin
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dw      !< suspended dust distribution weighting factors

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  v_grav  !< gravitational settling velocity

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  clay  !< clay mass concentration
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dust  !< total dust mass concentration
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  silt  !< silt mass concentration

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  dm_1  !< target for swapping of timelevels
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  dm_2  !< target for swapping of timelevels
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  dm_3  !< target for swapping of timelevels

    TYPE dust_bin_properties_type
      REAL(wp) :: density               !< particle density of a dust size bin, namelist parameter
      REAL(wp) :: diameter              !< effective diameter of a dust size bin, namelist parameter
      REAL(wp) :: lower_bound_diameter  !< minimum effective diameters represented by the dust size bin, namelist parameter
      REAL(wp) :: upper_bound_diameter  !< maximum effective diameters represented by the dust size bin,, namelist parameter
    END TYPE dust_bin_properties_type

    TYPE dust_flux_type
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  depo_flux_av  
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  emis_flux_av  
    END TYPE dust_flux_type

    TYPE progn_dust_type
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  init  !< initial concentration at given height

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  diss_s  !< for WS advection scheme, discretized artificial dissipation at southward-side
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  flux_s  !< for WS advection scheme, discretized 6th-order flux at northward-side

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  conc_av  !< time-averaged concentrations
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  diss_l   !< for WS advection scheme, discretized artificial dissipation at leftward-side
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  flux_l   !< for WS advection scheme, discretized 6th-order flux at leftward-side

       REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  conc     !< dust mass concentration at current Runge-Kutta step
       REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  conc_p   !< dust mass concentration at prognostic time level
       REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tconc_m  !< weighted concentration tendencies of previous Runge-Kutta steps
    END TYPE progn_dust_type

    TYPE salatation_bin_properties_type
       REAL(wp) ::  bin_mass_fraction_ssc  !< bin-specific mass fraction of corresponding soil separate class, namelist parameter
       REAL(wp) ::  density                !< particle density of a saltation size bin, namelist parameter
       REAL(wp) ::  diameter               !< effective diameter of a saltation size bin, namelist parameter
       REAL(wp) ::  mass_fraction_ssc      !< mass fraction of soil separate class, namelist parameter
    END TYPE salatation_bin_properties_type

    TYPE(dust_bin_properties_type), DIMENSION(:), ALLOCATABLE ::  db_prop  !< contains all properties of a dust bin

    TYPE(dust_flux_type), DIMENSION(:), ALLOCATABLE ::  dust_fluxes  !< array that stores all kind of fluxes

    TYPE(progn_dust_type), DIMENSION(:), ALLOCATABLE, TARGET :: dm  !< prognostic variable mass concentration (kg/m³)

    TYPE(salatation_bin_properties_type), DIMENSION(:), ALLOCATABLE ::  sb_prop  !< contains all properties of a saltation bin

    TYPE(surf_type), POINTER ::  surf  !< surf-type array for generalization purpose

    SAVE

    PRIVATE

    PUBLIC det_3d_data_averaging,                                                                  &
           det_actions,                                                                            &
           det_boundary_conditions,                                                                &
           det_check_data_output,                                                                  &
           det_check_data_output_pr,                                                               &
           det_check_data_output_ts,                                                               &
           det_check_parameters,                                                                   &
           det_data_output_2d,                                                                     &
           det_data_output_3d,                                                                     &
           det_define_netcdf_grid,                                                                 &
           det_exchange_horiz,                                                                     &
           det_header,                                                                             &
           det_init,                                                                               &
           det_init_arrays,                                                                        &
           det_non_advective_processes,                                                            &
           det_parin,                                                                              &
           det_prognostic_equations,                                                               &
           det_rrd_global,                                                                         &
           det_rrd_local,                                                                          &
           det_start_time,                                                                         &
           det_statistics,                                                                         &
           det_swap_timelevel,                                                                     &
           det_wrd_global,                                                                         &
           det_wrd_local,                                                                          &
           dm,                                                                                     &
           dm_2,                                                                                   &
           n_dust_bins


    INTERFACE cunningham_slip_flow_correction
       MODULE PROCEDURE cunningham_slip_flow_correction
    END INTERFACE cunningham_slip_flow_correction

    INTERFACE det_3d_data_averaging
       MODULE PROCEDURE det_3d_data_averaging
    END INTERFACE det_3d_data_averaging

    INTERFACE det_actions
       MODULE PROCEDURE det_actions
       MODULE PROCEDURE det_actions_ij
    END INTERFACE det_actions

    INTERFACE det_boundary_conditions
       MODULE PROCEDURE det_boundary_conditions
    END INTERFACE det_boundary_conditions

    INTERFACE det_calculate_dry_deposition_z01
       MODULE PROCEDURE det_calculate_dry_deposition_z01
    END INTERFACE det_calculate_dry_deposition_z01

    INTERFACE det_calculate_emission
       MODULE PROCEDURE det_calculate_emission
    END INTERFACE det_calculate_emission

    INTERFACE det_check_data_output
       MODULE PROCEDURE det_check_data_output
    END INTERFACE det_check_data_output

    INTERFACE det_check_data_output_pr
       MODULE PROCEDURE det_check_data_output_pr
    END INTERFACE det_check_data_output_pr

    INTERFACE det_check_data_output_ts
       MODULE PROCEDURE det_check_data_output_ts
    END INTERFACE det_check_data_output_ts

    INTERFACE det_check_parameters
       MODULE PROCEDURE det_check_parameters
    END INTERFACE det_check_parameters

    INTERFACE det_data_output_2d
       MODULE PROCEDURE det_data_output_2d
    END INTERFACE det_data_output_2d

    INTERFACE det_data_output_3d
       MODULE PROCEDURE det_data_output_3d
    END INTERFACE det_data_output_3d

    INTERFACE det_define_netcdf_grid
       MODULE PROCEDURE det_define_netcdf_grid
    END INTERFACE det_define_netcdf_grid

    INTERFACE det_exchange_horiz
       MODULE PROCEDURE det_exchange_horiz
    END INTERFACE det_exchange_horiz

    INTERFACE det_gravitational_settling
       MODULE PROCEDURE det_gravitational_settling
       MODULE PROCEDURE det_gravitational_settling_ij
    END INTERFACE det_gravitational_settling

    INTERFACE det_header
       MODULE PROCEDURE det_header
    END INTERFACE det_header

    INTERFACE det_init
       MODULE PROCEDURE det_init
    END INTERFACE det_init

    INTERFACE det_init_arrays
       MODULE PROCEDURE det_init_arrays
    END INTERFACE det_init_arrays

    INTERFACE det_non_advective_processes
       MODULE PROCEDURE det_non_advective_processes
       MODULE PROCEDURE det_non_advective_processes_ij
    END INTERFACE det_non_advective_processes

    INTERFACE det_parin
       MODULE PROCEDURE det_parin
    END INTERFACE det_parin

    INTERFACE det_prognostic_equations
       MODULE PROCEDURE det_prognostic_equations
       MODULE PROCEDURE det_prognostic_equations_ij
    END INTERFACE det_prognostic_equations

    INTERFACE det_rrd_global
       MODULE PROCEDURE det_rrd_global_ftn
       MODULE PROCEDURE det_rrd_global_mpi
    END INTERFACE det_rrd_global

    INTERFACE det_rrd_local
       MODULE PROCEDURE det_rrd_local_ftn
       MODULE PROCEDURE det_rrd_local_mpi
    END INTERFACE det_rrd_local

    INTERFACE det_statistics
      MODULE PROCEDURE det_statistics
    END INTERFACE det_statistics

    INTERFACE det_swap_timelevel
      MODULE PROCEDURE det_swap_timelevel
    END INTERFACE det_swap_timelevel

    INTERFACE det_wrd_global
       MODULE PROCEDURE det_wrd_global
    END INTERFACE det_wrd_global

    INTERFACE det_wrd_local
       MODULE PROCEDURE det_wrd_local
    END INTERFACE det_wrd_local


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute Cunningham slip-flow correction according to Jacobson (2005): Fundamentals of Atmospheric
!> Modeling, 2nd Edition, Eq. 15.30
!> @todo Move function to basic_constants_and_equations_mod since it can be used at different places
!> in PALM, e.g., in SALSA.
!--------------------------------------------------------------------------------------------------!
 FUNCTION cunningham_slip_flow_correction( density, particle_diameter )

    REAL(wp) ::  cunningham_slip_flow_correction  !< dynamic viscosity of air
    REAL(wp) ::  kn                               !< Knudsen number
    REAL(wp) ::  lambda_f                         !< molecular mean free path (m)
    REAL(wp) ::  v_th                             !< thermal speed of an air molecule

    REAL(wp), INTENT(IN) ::  density            !< air density
    REAL(wp), INTENT(IN) ::  particle_diameter  !< particle diameter


!
!-- Thermal velocity of an air molecule, Eq. 15.32.
    v_th = SQRT( 8.0_wp * k_boltzmann * pt_surface / ( pi * am_airmol ) )
!
!-- Mean free path, Eq. 15.24
    lambda_f = 2.0_wp * air_viscosity / ( density * v_th )
!
!-- Knudsen number, Eq. 15.23
    kn = lambda_f / ( particle_diameter * 0.5_wp )
!
!-- Cunningham slip-flow correction, Eq. 15.30
    cunningham_slip_flow_correction = 1.0_wp + kn * ( 1.249_wp + 0.42_wp * EXP( -0.87_wp / kn ) )

 END FUNCTION cunningham_slip_flow_correction


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average output quantities as well as allocate the array necessary for storing the
!> average. Note, if you just specify an averaged output quantity in the _p3dr file during restarts
!> the first output includes the time between the beginning of the restart run and the first output
!> time (not necessarily the whole averaging_interval you have specified in your_p3d/_p3dr file ).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_3d_data_averaging( mode, variable )

    USE averaging,                                                                                 &
        ONLY:  clay_av,                                                                            &
               dust_av,                                                                            &
               silt_av

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d


    CHARACTER(LEN=*), INTENT(IN) ::  mode      !< averaging interface mode
    CHARACTER(LEN=*), INTENT(IN) ::  variable  !< variable name

    INTEGER(iwp) ::  char_to_int  !< for converting character to integer and index for dust size bin
    INTEGER(iwp) ::  i            !< loop index x direction
    INTEGER(iwp) ::  j            !< loop index y direction
    INTEGER(iwp) ::  k            !< loop index z direction
    INTEGER(iwp) ::  m            !< running index for surface elements


    IF ( time_since_reference_point < det_start_time )  RETURN

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'clay' )
             IF ( .NOT. ALLOCATED( clay_av ) )  THEN
                ALLOCATE( clay_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             clay_av = 0.0_wp

          CASE ( 'dust' )
             IF ( .NOT. ALLOCATED( dust_av ) )  THEN
                ALLOCATE( dust_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             dust_av = 0.0_wp

          CASE ( 'silt' )
             IF ( .NOT. ALLOCATED( silt_av ) )  THEN
                ALLOCATE( silt_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             silt_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

       IF ( variable(6:11) ==  'mc_bin' )  THEN
           READ( variable(12:),* ) char_to_int
           IF ( .NOT. ALLOCATED( dm(char_to_int)%conc_av ) )  THEN
              ALLOCATE( dm(char_to_int)%conc_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
           ENDIF
           dm(char_to_int)%conc_av(nzb:nzt+1,nys:nyn,nxl:nxr) = 0.0_wp
       ENDIF

       IF ( variable(6:19) ==  'emis_flux*_bin' )  THEN
          READ( variable(20:),* ) char_to_int
          IF ( .NOT. ALLOCATED( dust_fluxes(char_to_int)%emis_flux_av ) )  THEN
             ALLOCATE( dust_fluxes(char_to_int)%emis_flux_av(nys:nyn,nxl:nxr) )
          ENDIF
          dust_fluxes(char_to_int)%emis_flux_av(nys:nyn,nxl:nxr) = 0.0_wp
       ENDIF

       IF ( variable(6:19) ==  'depo_flux*_bin' )  THEN
          READ( variable(20:),* ) char_to_int
          IF ( .NOT. ALLOCATED( dust_fluxes(char_to_int)%depo_flux_av ) )  THEN
             ALLOCATE( dust_fluxes(char_to_int)%depo_flux_av(nys:nyn,nxl:nxr) )
          ENDIF
          dust_fluxes(char_to_int)%depo_flux_av(nys:nyn,nxl:nxr) = 0.0_wp
       ENDIF

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'clay' )
             IF ( ALLOCATED( clay_av ) )  THEN
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. clay_calculated )  CALL det_actions( 'calculate clay' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         clay_av(k,j,i) = clay_av(k,j,i) + clay(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'dust' )
             IF ( ALLOCATED( dust_av ) )  THEN
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. dust_calculated )  CALL det_actions( 'calculate dust' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         dust_av(k,j,i) = dust_av(k,j,i) + dust(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'silt' )
             IF ( ALLOCATED( silt_av ) )  THEN
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. silt_calculated )  CALL det_actions( 'calculate silt' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         silt_av(k,j,i) = silt_av(k,j,i) + silt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

       IF ( variable(6:11) ==  'mc_bin' )  THEN
          READ( variable(12:),* ) char_to_int
          IF ( ALLOCATED( dm(char_to_int)%conc_av ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      dm(char_to_int)%conc_av(k,j,i) = dm(char_to_int)%conc_av(k,j,i) +            &
                                                       dm(char_to_int)%conc(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       IF ( variable(6:19) ==  'emis_flux*_bin' )  THEN
          READ( variable(20:),* ) char_to_int
          IF ( ALLOCATED( dust_fluxes(char_to_int)%emis_flux_av ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                      IF ( surf_def%upward_top(m) )  THEN
                         dust_fluxes(char_to_int)%emis_flux_av(j,i) =                              &
                                                      dust_fluxes(char_to_int)%emis_flux_av(j,i) + &
                                                      surf_def%dm_emis_flux(char_to_int,m)
                      ENDIF
                   ENDDO

                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      IF ( surf_lsm%upward_top(m) )  THEN
                         dust_fluxes(char_to_int)%emis_flux_av(j,i) =                              &
                                                      dust_fluxes(char_to_int)%emis_flux_av(j,i) + &
                                                      surf_lsm%dm_emis_flux(char_to_int,m)
                      ENDIF
                   ENDDO

                ENDDO
             ENDDO
          ENDIF
       ENDIF

       IF ( variable(6:19) ==  'depo_flux*_bin' )  THEN
          READ( variable(20:),* ) char_to_int
          IF ( ALLOCATED( dust_fluxes(char_to_int)%depo_flux_av ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                      IF ( surf_def%upward_top(m) )  THEN
                         dust_fluxes(char_to_int)%depo_flux_av(j,i) =                              &
                                                      dust_fluxes(char_to_int)%depo_flux_av(j,i) + &
                                                      surf_def%dm_depo_flux(char_to_int,m)
                      ENDIF
                   ENDDO

                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      IF ( surf_lsm%upward_top(m) )  THEN
                         dust_fluxes(char_to_int)%depo_flux_av(j,i) =                              &
                                                      dust_fluxes(char_to_int)%depo_flux_av(j,i) + &
                                                      surf_lsm%dm_depo_flux(char_to_int,m)
                      ENDIF
                   ENDDO

                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'clay' )
             IF ( ALLOCATED( clay_av ) )  THEN
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. clay_calculated )  CALL det_actions( 'calculate clay' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         clay_av(k,j,i) = clay_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'dust' )
             IF ( ALLOCATED( dust_av ) )  THEN
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. dust_calculated )  CALL det_actions( 'calculate dust' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         dust_av(k,j,i) = dust_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'silt' )
             IF ( ALLOCATED( silt_av ) )  THEN
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. silt_calculated )  CALL det_actions( 'calculate silt' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         silt_av(k,j,i) = silt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

       IF ( variable(6:11) ==  'mc_bin' )  THEN
          READ( variable(12:),* ) char_to_int
          IF ( ALLOCATED( dm(char_to_int)%conc_av ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      dm(char_to_int)%conc_av(k,j,i) = dm(char_to_int)%conc_av(k,j,i) /            &
                                                       REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       IF ( variable(6:19) ==  'emis_flux*_bin' )  THEN
          READ( variable(20:),* ) char_to_int
          IF ( ALLOCATED( dust_fluxes(char_to_int)%emis_flux_av ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   dust_fluxes(char_to_int)%emis_flux_av(j,i) =                                    &
                                                      dust_fluxes(char_to_int)%emis_flux_av(j,i) / &
                                                      REAL( average_count_3d, KIND = wp )
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       IF ( variable(6:19) ==  'depo_flux*_bin' )  THEN
          READ( variable(20:),* ) char_to_int
          IF ( ALLOCATED( dust_fluxes(char_to_int)%depo_flux_av ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   dust_fluxes(char_to_int)%depo_flux_av(j,i) =                                    &
                                                      dust_fluxes(char_to_int)%depo_flux_av(j,i) / &
                                                      REAL( average_count_3d, KIND = wp )
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ENDIF

 END SUBROUTINE det_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_actions( location )

    CHARACTER(LEN=*), INTENT(IN) ::  location  !< location string, describes position of the call

    INTEGER(iwp) ::  i   !< grid index in x-direction
    INTEGER(iwp) ::  id  !< loop variable for dust size bins
    INTEGER(iwp) ::  j   !< grid index in y-direction
    INTEGER(iwp) ::  k   !< grid index in z-direction


    IF ( time_since_reference_point < det_start_time )  RETURN

!
!-- No calls for single grid points are allowed at locations before and after the timestep, since
!-- these calls are not within an i,j-loop.
    SELECT CASE ( location )

       CASE ( 'calculate clay' )

          IF ( clay_calculated )  RETURN

          IF ( .NOT. ALLOCATED( clay ) )  ALLOCATE( clay(nzb:nzt+1,nys:nyn,nxl:nxr) )

          clay = 0.0_wp
          DO  id = 1, n_dust_bins
             IF ( diameter_dust(id) <= 4.0E-6_wp )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         clay(k,j,i) = clay(k,j,i) + dm(id)%conc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

          clay_calculated = .TRUE.

       CASE ( 'calculate dust' )

          IF ( dust_calculated )  RETURN

          IF ( .NOT. ALLOCATED( dust ) )  ALLOCATE( dust(nzb:nzt+1,nys:nyn,nxl:nxr) )

          dust = 0.0_wp
          DO  id = 1, n_dust_bins
             IF ( diameter_dust(id) <= 63.0E-6_wp )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         dust(k,j,i) = dust(k,j,i) + dm(id)%conc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

          dust_calculated = .TRUE.

       CASE ( 'calculate silt' )

          IF ( silt_calculated )  RETURN

          IF ( .NOT. ALLOCATED( silt ) )  ALLOCATE( silt(nzb:nzt+1,nys:nyn,nxl:nxr) )

          silt = 0.0_wp
          DO  id = 1, n_dust_bins
             IF ( diameter_dust(id) > 4.0E-6_wp  .AND.  diameter_dust(id) <= 63.0E-6_wp )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         silt(k,j,i) = silt(k,j,i) + dm(id)%conc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

          silt_calculated = .TRUE.

       CASE DEFAULT

          CONTINUE

    END SELECT

 END SUBROUTINE det_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j. So far, this is just a framework for performing cache-optimized
!> module-specific actions while in time-integration.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_actions_ij( i, j, location )

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp), INTENT(IN) ::  i  !< loop variable in x-direction
    INTEGER(iwp), INTENT(IN) ::  j  !< loop variable in y-direction


    IF ( time_since_reference_point < det_start_time )  RETURN

!
!-- Next line is to avoid warning about unused variables- Please remove.
    IF ( i == 0 .OR. j == 0 )  CONTINUE

    SELECT CASE ( location )

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE det_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set boundary conditions for prognostic variables in the dust emission and transport model (DET).
!> @todo Implement initial gradient (ibc_dm_t == 2) as top boundary condition
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_boundary_conditions

    USE surface_mod,                                                                               &
        ONLY:  bc_hv


    INTEGER(iwp) ::  i   !< grid index x direction.
    INTEGER(iwp) ::  id  !< index for dust size bins
    INTEGER(iwp) ::  j   !< grid index y direction.
    INTEGER(iwp) ::  k   !< grid index z direction.
    INTEGER(iwp) ::  m   !< running index surface elements.


    IF ( time_since_reference_point < det_start_time )  RETURN

!
!-- Boundary conditions for dm%conc.
!-- Bottom boundary: Neumann condition because dust mass flux is always given. Run loop over all
!-- non-natural and natural walls. Note, in wall-datatype the k,j,i coordinate belong to the
!-- atmospheric grid point, therefore, set s_p at k+koff, j+joff, i+ioff, respectively.
    IF ( ibc_dm_b == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE(i, j, k, m, id)
       DO  m = 1, bc_hv%ns
          i = bc_hv%i(m)
          j = bc_hv%j(m)
          k = bc_hv%k(m)
          DO  id = 1, n_dust_bins
             dm(id)%conc_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = dm(id)%conc_p(k,j,i)
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

!
!-- Top boundary conditions: Dirichlet (ibc_dm_t == 0) or Neumann (ibc_dm_t == 1).
    IF ( ibc_dm_t == 0 )  THEN
       DO  id = 1, n_dust_bins
          dm(id)%conc_p(nzt+1,:,:) = dm(id)%conc(nzt+1,:,:)
       ENDDO
    ELSEIF ( ibc_dm_t == 1 )  THEN
       DO  id = 1, n_dust_bins
          dm(id)%conc_p(nzt+1,:,:) = dm(id)%conc_p(nzt,:,:)
       ENDDO
    ENDIF

 END SUBROUTINE det_boundary_conditions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The parameterization of the surface deposition emission flux is based on the dry deposition
!> scheme of Zhang et al. (2001): A size-segregated particle dry deposition scheme for an
!> atmospheric aerosol module (Z01). In this scheme, a dry deposition velocity is estimated based on
!> electrical resistance analogy. Aerodynamic and surface resistances are calculated. The surface
!> resistance considers Brownian diffusion, impaction and rebound of particles. Interception is
!> neglected. Relevant equations can also be found in Seinfeld & Pandis (2016): Atmospheric
!> Chemistry and Physics - From Air Pollution to Climate Change, 3rd Edition, Jacobson (2005):
!> Fundamentals of Atmospheric Modeling, 2nd Edition, and Slinn (1982): Predictions for particle
!> deposition to vegetative canopies. The routine is called for a certain surface element and dust
!> bin.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_calculate_dry_deposition_z01( i, j, k, m, depo_flux )

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  psi_h


    INTEGER(iwp), INTENT(IN) ::  i   !< loop index x direction
    INTEGER(iwp), INTENT(IN) ::  j   !< loop index y direction
    INTEGER(iwp), INTENT(IN) ::  k   !< height index
    INTEGER(iwp), INTENT(IN) ::  m   !< index for surface element

    INTEGER(iwp) ::  id  !< loop index for dust size bin

    REAL(wp) ::  diff_brown       !< Brownian diffusion for particles
    REAL(wp) ::  e_brownian_diff  !< collection efficiency from Brownian diffusion
    REAL(wp) ::  e_impaction      !< collection efficiency from impaction
    REAL(wp) ::  e_rebound        !< inverse of rebound efficiency
    REAL(wp) ::  cc               !< Cunningham slip-flow correction factor
    REAL(wp) ::  r_aero           !< aerodynamic resistance at current surface element
    REAL(wp) ::  r_surf           !< surface resistance at current surface element
    REAL(wp) ::  schmidt_number   !< Schmidt number
    REAL(wp) ::  stokes_num       !< Stokes number
    REAL(wp) ::  v_depo           !< deposition velocity

    REAL(wp), DIMENSION(n_dust_bins), INTENT(OUT) ::  depo_flux  !<deposition flux in kg/(m²s) for a dust bin


    DO  id = 1, n_dust_bins
!
!--    Calculate aerodynamic resistance, Zhang et al. (2001) Eq. 4. Reference height z is z_mo.
!--    ATTENTION: The below formula has been used in Klamt et al. (2024), but it occasionally
!--               gives negative resistances. Therefore, the resistance as calculated for the LSM
!--               is used here (see routine calc_aerodynamic_resistance in surface_layer_fluxes_mod.
!       r_aero = ( surf%ln_z_z0(m) - psi_h( surf%z_mo(m)/ surf%ol(m) ) ) / ( kappa * surf_def%us(m) )
       r_aero = surf_def%r_a(m)
!
!--    Cunningham slip-flow correction.
       cc = cunningham_slip_flow_correction( v_grav(k,id), db_prop(id)%diameter )
!
!--    Calculate collection efficiency from Brownian diffusion, Zhang et al. (2001) Eq. 6 and
!--    Seinfeld & Pandis (2006) Eq. 19.20, 19.21.
       diff_brown      = k_boltzmann * pt_surface * cc /                                           &
                         ( 3.0_wp * pi * air_viscosity * db_prop(id)%diameter )
       schmidt_number  = air_viscosity / rho_air(k) / diff_brown
       e_brownian_diff = schmidt_number**( -brownian_diffusion_coefficient )
!
!--    Calculate collection efficiency from impaction, Zhang et al. (2001) Eq. 7c and Seinfeld &
!--    Pandis (2006) Eq. 19.22, 19.23.
       stokes_num = v_grav(k,id) * surf%us(m)**2 / ( g * air_viscosity / rho_air(k) )
       e_impaction = ( stokes_num / ( alpha_imp + stokes_num ) )**2
!
!--    Calculate reduction in collection caused by rebound (fraction of particles that sticks to
!--    ground), Zhang et al. (2001) Eq. 9 and Slinn (1982) Eq. 29
       e_rebound = EXP( -SQRT( stokes_num ) )
!
!--    Calculate surface resistance:
       r_surf = 1.0_wp / ( 3.0_wp * MAX( 1.0E-8_wp, surf%us(m) ) *                                 &
                           ( e_brownian_diff + e_impaction ) * e_rebound )
!
!--    Calculate deposition velocity according to Seinfeld and Pandis (2006) Eq 19.7
       v_depo = 1.0_wp / ( r_aero + r_surf + r_aero * r_surf * v_grav(k,id) ) + v_grav(k,id)
!
!--    Calculate deposition flux ("kinematic" units, i.e., kg/(m²s)) for a dust bin, Seinfeld and
!--    Pandis (2006) Eq 19.1
       depo_flux(id) = -v_depo * dm(id)%conc(1,j,i)

    ENDDO

 END SUBROUTINE det_calculate_dry_deposition_z01


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The parameterization of the surface dust emission flux is based on the AFWA dust emission scheme.
!> In the AFWA scheme, dust emission is handled as a two-part process, wherein large particle
!> saltation from coarser dust- and sand-sized particles is triggered by wind shear and leads to a
!> fine-particle (dust-sized) bulk emission flux by saltation bombardment and aggregate
!> disintegration. The bulk dust emission flux is then further distributed among different size
!> bins (LeGrand et al., 2019: The afwa dust emission scheme for the gocart aerosol model in
!> wrf-chem v3.8.1). The routine is called for a certain surface element and dust bin.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_calculate_emission( density, m, em_flux )

    INTEGER(iwp) ::  is  !< loop index for saltation bins
    INTEGER(iwp) ::  id  !< loop index for dust size bin

    INTEGER(iwp), INTENT(IN) ::  m   !< index for surface element

    REAL(wp) ::  g_salt  !< total horizontal streamwise saltation flux in kg/(m*s)

    REAL(wp), INTENT(IN) ::  density  !< air density

    REAL(wp), DIMENSION(n_dust_bins), INTENT(OUT) ::  em_flux  !< vertical surface dust flux in kg/(m²s)

    REAL(wp), DIMENSION(n_saltation_bins) ::  h_salt  !< bin-specific horizontal streamwise saltation flux in kg/(m*s)


!
!-- Calculate horizontal saltation flux, LeGrand et al. (2019) Eq. 10 and 13.
    g_salt = 0.0_wp
    DO  is = 1, n_saltation_bins
       h_salt(is) = MAX( 0.0_wp, density * surf%us(m)**3  / g *                                    &
                                 ( 1.0_wp + surf%us_t(m,is) / surf%us(m) ) *                       &
                                 ( 1.0_wp - ( surf%us_t(m,is) / surf%us(m) )**2 ) )
       g_salt = g_salt +  h_salt(is) * ds_rel(is)
    ENDDO
!
!-- Calculation of vertical dust emission flux ("kinematic" units, i.e., kg/(m²s)) for a dust bin.
    DO  id = 1, n_dust_bins
       em_flux(id) = alpha_s * dw(id) * g_salt
    ENDDO

 END SUBROUTINE det_calculate_emission


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for DET.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_check_data_output( var, unit )

    CHARACTER(LEN=*), INTENT(IN)    ::  var  !< variable name

    CHARACTER(LEN=*), INTENT(INOUT) ::  unit  !< physical unit

    INTEGER(iwp) ::  char_to_int  !< for converting character to integer


    SELECT CASE ( TRIM( var ) )

       CASE ( 'clay', 'dust', 'silt' )
          unit = 'kg/m3'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT
!
!-- Treat bin-specific output, e.g., dust_mc_bin1, dust_mc_bin2, etc.
    IF ( var(6:11) ==  'mc_bin' )  THEN
       READ( var(12:),* )  char_to_int
       IF ( char_to_int >= 1  .AND.  char_to_int <= n_dust_bins )  THEN
          unit = 'kg/m3'
       ELSE
          unit = 'illegal'
          RETURN
       ENDIF
    ELSEIF ( var(6:19) ==  'emis_flux*_bin' )  THEN
       READ( var(20:),* )  char_to_int
       IF ( char_to_int >= 1  .AND.  char_to_int <= n_dust_bins )  THEN
          unit = 'kg m-2 s-1'
       ELSE
          unit = 'illegal'
          RETURN
       ENDIF
    ELSEIF ( var(6:19) ==  'depo_flux*_bin' )  THEN
       READ( var(20:),* )  char_to_int
       IF ( char_to_int >= 1  .AND.  char_to_int <= n_dust_bins )  THEN
          unit = 'kg m-2 s-1'
       ELSE
          unit = 'illegal'
          RETURN
       ENDIF
    ENDIF

 END SUBROUTINE det_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check profile data output for DET. So far, only a framework is given and no profile output is
!> implemented.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_check_data_output_pr( var, var_count, unit, dopr_unit )

    USE arrays_3d,                                                                                 &
        ONLY:  zu

    USE profil_parameter,                                                                          &
        ONLY:  dopr_index

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               pr_palm,                                                                            &
               statistic_regions


    CHARACTER (LEN=*), INTENT(IN)    ::  var  !< variable name

    CHARACTER (LEN=*), INTENT(INOUT) ::  dopr_unit  !< physical unit
    CHARACTER (LEN=*), INTENT(INOUT) ::  unit       !< physical unit

    INTEGER(iwp), INTENT(IN) ::  var_count  !< number of data-output quantity


    SELECT CASE ( TRIM( var ) )

       CASE ( 'clay' )
          det_pr_count = det_pr_count + 1
          det_pr_index(det_pr_count) = 1
          dopr_index(var_count) = pr_palm + max_pr_cs + det_pr_count
          dopr_unit = 'kg/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )
       CASE ( 'dust' )
          det_pr_count = det_pr_count + 1
          det_pr_index(det_pr_count) = 2
          dopr_index(var_count) = pr_palm + max_pr_cs + det_pr_count
          dopr_unit = 'kg/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )
       CASE ( 'silt' )
          det_pr_count = det_pr_count + 1
          det_pr_index(det_pr_count) = 3
          dopr_index(var_count) = pr_palm + max_pr_cs + det_pr_count
          dopr_unit = 'kg/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE det_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    INTEGER(iwp), INTENT(IN)     ::  dots_max  !< maximum number of timeseries data output

    INTEGER(iwp), INTENT(INOUT)  ::  dots_num  !< index number of timeseries data output

    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_label  !< label of timeseries
    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_unit   !< unit of timeseries

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( dots_label(1)(1:1) == ' '  .OR.  dots_unit(1)(1:1) == ' ' )  CONTINUE

!
!-- For each time series quantity a label and a unit is given, which will be used for the
!-- NetCDF file. The value of dots_num has to be increased by the number of new time series
!-- quantities. Its old value has to be stored in dots_num_palm. See routine user_statistics on how
!-- to calculate and output these quantities.
    dots_start_index_det = dots_num + 1

!
!-- Mean bulk deposition flux.
!-- Here and in the following "bulk" means the sum over all dust bins.
    dots_num = dots_num + 1
    dots_num_det = dots_num_det + 1
    dots_label(dots_num) = 'dm_df'
    dots_unit(dots_num)  = 'kg m-2 s-1'

!
!-- Mean bulk emission flux.
    dots_num = dots_num + 1
    dots_num_det = dots_num_det + 1
    dots_label(dots_num) = 'dm_ef'
    dots_unit(dots_num)  = 'kg m-2 s-1'

!
!-- Maximum of bulk deposition flux.
    dots_num = dots_num + 1
    dots_num_det = dots_num_det + 1
    dots_label(dots_num) = 'dm_df_max'
    dots_unit(dots_num)  = 'kg m-2 s-1'

!
!-- Maximum of bulk emission flux.
    dots_num = dots_num + 1
    dots_num_det = dots_num_det + 1
    dots_label(dots_num) = 'dm_ef_max'
    dots_unit(dots_num)  = 'kg m-2 s-1'


 END SUBROUTINE det_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for the dust emission and transport model (DET).
!> @todo Check for other top boundary conditions like nested or initial_gradient and other lateral 
!> boundary conditions than cyclic as soon as they are implemented.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  bc_lr,                                                                              &
               bc_ns,                                                                              &
               data_output_pr,                                                                     &
               nesting_offline,                                                                    &
               ocean_mode,                                                                         &
               topography,                                                                         &
               urban_surface


    INTEGER(iwp) ::  count  !< counter variable
    INTEGER(iwp) ::  i      !< running index to determine number of det output profiles

!
!-- Set bottom boundary condition flag. So far, only Neumann condition allowed because dust mass
!-- flux is always active.
    IF ( bc_dm_b == 'neumann' )  THEN
       ibc_dm_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_dm_b = "' // TRIM( bc_dm_b ) // '"'
       CALL message( 'det_check_parameters', 'DET0001', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set top boundary conditions flag.
    IF ( bc_dm_t == 'dirichlet' )  THEN
       ibc_dm_t = 0
    ELSEIF ( bc_dm_t == 'neumann' )  THEN
       ibc_dm_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_dm_t = "' // TRIM( bc_dm_t ) // '"'
       CALL message( 'det_check_parameters', 'DET0001', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check left and right boundary conditions. First set default value if not set by user.
    IF ( bc_dm_l == 'undefined' )  bc_dm_l = bc_lr
    IF ( bc_dm_r == 'undefined' )  bc_dm_l = bc_lr
!
!-- Check boundary conditions that are set by the user.
    IF ( bc_dm_l /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_dm_l = "' // TRIM( bc_dm_l ) // '"'
       CALL message( 'det_check_parameters', 'DET0001', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( bc_dm_r /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_dm_r = "' // TRIM( bc_dm_r ) // '"'
       CALL message( 'det_check_parameters', 'DET0001', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check north and south boundary conditions. First set default value if not set by user.
    IF ( bc_dm_n == 'undefined' )  bc_dm_n = bc_ns
    IF ( bc_dm_s == 'undefined' )  bc_dm_s = bc_ns
!
!-- Check boundary conditions that are set by the user.
    IF ( bc_dm_n /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_dm_n = "' // TRIM( bc_dm_n ) // '"'
       CALL message( 'det_check_parameters', 'DET0001', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( bc_dm_s /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_dm_s = "' // TRIM( bc_dm_s ) // '"'
       CALL message( 'det_check_parameters', 'DET0001', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Cyclic conditions must be set identically at opposing boundaries.
    IF ( ( bc_dm_l == 'cyclic' .AND. bc_dm_r /= 'cyclic' )  .OR.                                   &
         ( bc_dm_r == 'cyclic' .AND. bc_dm_l /= 'cyclic' ) )  THEN
       message_string = 'inconsistent left and right boundary conditions'
       CALL message( 'det_check_parameters', 'DET0002', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ( bc_dm_n == 'cyclic' .AND. bc_dm_s /= 'cyclic' )  .OR.                                   &
         ( bc_dm_s == 'cyclic' .AND. bc_dm_n /= 'cyclic' ) )  THEN
       message_string = 'inconsistent north and south boundary conditions'
       CALL message( 'det_check_parameters', 'DET0003', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check that number of bins is not larger than the limit.
    IF ( n_dust_bins > n_dust_bins_max )  THEN
       WRITE( message_string, '(A,I2,A,I2)' )  'number of specified dust bins ', n_dust_bins,      &
                                               ' exceeds the limit of ', n_dust_bins_max
       CALL message( 'det_check_parameters', 'DET0004', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( n_saltation_bins > n_saltation_bins_max )  THEN
       WRITE( message_string, '(A,I2,A,I2)' )  'number of specified saltation bins ',              &
                                               n_saltation_bins, ' exceeds the limit of ',         &
                                               n_saltation_bins_max
       CALL message( 'det_check_parameters', 'DET0005', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check that only dust-sized particles are specified by the user. Otherwise print a warning.
    IF ( ANY( diameter_dust(1:n_dust_bins) > 63.0E-6_wp )  .AND.  ANY( diameter_dust /= not_set ) )&
    THEN
       message_string = 'effective diameters of some/all specified dust bins exceed 63.0E-6'
       CALL message( 'det_check_parameters', 'DET0006', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Check if the first given diameter is less than or equal to 4microns (clay) because the 
!-- calculation of the sandblasting efficiency requires at least one size bin of clay particles.
    IF ( diameter_saltation(1) > 4.0E-6_wp  .AND.  ANY( diameter_dust /= not_set ) )  THEN
       WRITE( message_string, '(A,E11.4,A)' )  'diameter_saltation(1) = ',diameter_saltation(1),   &
                                               ' is larger than 4.0E-6'
       CALL message( 'det_check_parameters', 'DET0007', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check that the number of saltation size bins the user has specified matches the number of the
!-- given values for each saltation bin property
    count = 0
    DO  i = 1, n_saltation_bins_max
       IF ( bin_mass_fraction_ssc(i) /= not_set  .AND.  diameter_saltation(i) /= not_set  .AND.    &
            mass_fraction_ssc(i) /= not_set .AND. particle_density_saltation(i) /= not_set )  THEN
          count = count + 1
       ENDIF
    ENDDO
    IF ( count /= n_saltation_bins .AND. count /= 0 )  THEN
       message_string = 'inconsistent number of specified saltation bins'
       CALL message( 'det_check_parameters', 'DET0008', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check that the number of dust size bins the user has specified matches the number of the
!-- given values for each dust bin property
    count = 0
    DO  i = 1, n_dust_bins_max
       IF ( diameter_dust(i) /= not_set  .AND.  lower_bound_diameter(i) /= not_set  .AND.          &
            upper_bound_diameter(i) /= not_set .AND. particle_density_dust(i) /= not_set )  THEN
          count = count + 1
       ENDIF
    ENDDO
    IF ( count /= n_dust_bins .AND. count /= 0 )  THEN
       message_string = 'inconsistent number of specified dust bins'
       CALL message( 'det_check_parameters', 'DET0009', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for illegal combinations.
    message_string = ''
    IF ( topography /= 'flat' )  message_string = 'topography = "' // TRIM( topography ) // '", '
    IF ( nesting_offline )  THEN
       message_string =  TRIM( message_string ) // ' nesting (self/offline),'
    ENDIF
    IF ( ocean_mode )  message_string =  TRIM( message_string ) // 'ocean-mode, '
    IF ( urban_surface )  message_string =  TRIM( message_string ) // 'urban surfaces, '
    IF ( TRIM( message_string ) /= '' )  THEN
       message_string = 'dust module does not allow to use ' // TRIM( message_string )
       CALL message( 'det_check_parameters', 'DET0010', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Determine the number of det profiles and append them to the standard data output.
    i = 1
    DO  WHILE ( data_output_pr(i)  /= ' '  .AND.  i <= SIZE( data_output_pr ) )
       IF ( TRIM( data_output_pr(i) ) == 'clay' .OR.                                               &
            TRIM( data_output_pr(i) ) == 'dust' .OR.                                               &
            TRIM( data_output_pr(i) ) == 'silt'  )  THEN
          max_pr_det = max_pr_det + 1
       ENDIF
       i = i + 1
    ENDDO

 END SUBROUTINE det_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k) and sets the grid on which it is defined. Allowed values for grid are "zu" and "zw".
!> "mode" from argument list is not used.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )

    USE averaging,                                                                                 &
        ONLY:  clay_av,                                                                            &
               dust_av,                                                                            &
               silt_av


    CHARACTER(LEN=30) ::  dust_bin_name  !< name of dust size bin, e.g. dust_mc_bin10

    CHARACTER(LEN=*), INTENT(IN) ::  mode      !< either 'xy', 'xz' or 'yz'
    CHARACTER(LEN=*), INTENT(IN) ::  variable  !< variable name

    CHARACTER(LEN=*), INTENT(INOUT) ::  grid  !< name of vertical grid, either zu or zw

    INTEGER(iwp) ::  char_len     !< length of a character string
    INTEGER(iwp) ::  char_to_int  !< for converting character to integer and index for dust size bin
    INTEGER(iwp) ::  i            !< grid index along x-direction
    INTEGER(iwp) ::  j            !< grid index along y-direction
    INTEGER(iwp) ::  k            !< grid index along z-direction
    INTEGER(iwp) ::  m            !< running index for surface elements

    INTEGER(iwp), INTENT(IN) ::  av      !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp), INTENT(IN) ::  nzb_do  !< lower limit of the domain (usually nzb)
    INTEGER(iwp), INTENT(IN) ::  nzt_do  !< upper limit of the domain (usually nzt+1)

    LOGICAL, INTENT(OUT) ::  two_d  !< flag indicating 2D variables (horizontal cross sections)

    LOGICAL, INTENT(INOUT) ::  found  !< flag if output variable is found

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf  !< local array to which output data is resorted to


    IF ( time_since_reference_point < det_start_time )  RETURN

    found         = .FALSE.
    two_d         = .FALSE.
    char_len      = LEN_TRIM( variable )
    dust_bin_name = TRIM( variable(1:char_len-3) )

    IF ( variable(6:11)  == 'mc_bin' .AND. ( ( variable(char_len-2:) == '_xy' )  .OR.              &
                                             ( variable(char_len-2:) == '_xz')   .OR.              &
                                             ( variable(char_len-2:) == '_yz') ) )                 &
    THEN
       READ( dust_bin_name(12:), * ) char_to_int
       IF (av == 0)  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_do, nzt_do
                   IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                    local_pf(i,j,k) = dm(char_to_int)%conc(k,j,i)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ELSE
          IF ( .NOT. ALLOCATED( dm(char_to_int)%conc_av ) )  THEN
             ALLOCATE( dm(char_to_int)%conc_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             dm(char_to_int)%conc_av = 0.0_wp
          ENDIF
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_do, nzt_do
                   IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                      local_pf(i,j,k) = dm(char_to_int)%conc_av(k,j,i)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       IF ( mode == 'xy' ) grid = 'zu'
       found = .TRUE.
    ENDIF

    IF ( variable(6:19)  == 'emis_flux*_bin' .AND. variable(char_len-2:) == '_xy' )                &
    THEN
       READ( dust_bin_name(20:), * ) char_to_int
       IF (av == 0)  THEN
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             local_pf(i,j,nzb+1) = MERGE( surf_def%dm_emis_flux(char_to_int,m),                    &
                                          local_pf(i,j,nzb+1), surf_def%upward(m) )
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             local_pf(i,j,nzb+1) = MERGE( surf_lsm%dm_emis_flux(char_to_int,m),                    &
                                         local_pf(i,j,nzb+1), surf_lsm%upward(m) )
          ENDDO
       ELSE
          IF ( .NOT. ALLOCATED( dust_fluxes(char_to_int)%emis_flux_av ) )  THEN
             ALLOCATE( dust_fluxes(char_to_int)%emis_flux_av(nys:nyn,nxl:nxr) )
             dust_fluxes(char_to_int)%emis_flux_av = 0.0_wp
          ENDIF
          DO  i = nxl, nxr
             DO  j = nys, nyn
                local_pf(i,j,nzb+1) = dust_fluxes(char_to_int)%emis_flux_av(j,i)
             ENDDO
          ENDDO
       ENDIF
       grid  = 'zu1'
       found = .TRUE.
       two_d = .TRUE.
    ENDIF

    IF ( variable(6:19)  == 'depo_flux*_bin' .AND. variable(char_len-2:) == '_xy' )                &
    THEN
       READ( dust_bin_name(20:), * ) char_to_int
       IF (av == 0)  THEN
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             local_pf(i,j,nzb+1) = MERGE( surf_def%dm_depo_flux(char_to_int,m),                    &
                                          local_pf(i,j,nzb+1), surf_def%upward(m) )
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             local_pf(i,j,nzb+1) = MERGE( surf_lsm%dm_depo_flux(char_to_int,m),                    &
                                         local_pf(i,j,nzb+1), surf_lsm%upward(m) )
          ENDDO
       ELSE
          IF ( .NOT. ALLOCATED( dust_fluxes(char_to_int)%depo_flux_av ) )  THEN
             ALLOCATE( dust_fluxes(char_to_int)%depo_flux_av(nys:nyn,nxl:nxr) )
             dust_fluxes(char_to_int)%depo_flux_av = 0.0_wp
          ENDIF
          DO  i = nxl, nxr
             DO  j = nys, nyn
                local_pf(i,j,nzb+1) = dust_fluxes(char_to_int)%depo_flux_av(j,i)
             ENDDO
          ENDDO
       ENDIF
       grid  = 'zu1'
       found = .TRUE.
       two_d = .TRUE.
    ENDIF

    IF ( found ) RETURN

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'clay_xy', 'clay_xz', 'clay_yz' )

!
!--       Calculate diagnostic quantities, if not done so far.
          IF ( .NOT. clay_calculated )  CALL det_actions( 'calculate clay' )

          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = clay(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( clay_av ) )  THEN
                ALLOCATE( clay_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                clay_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = clay_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          found = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'dust_xy', 'dust_xz', 'dust_yz' )
!
!--       Calculate diagnostic quantities, if not done so far.
          IF ( .NOT. dust_calculated )  CALL det_actions( 'calculate dust' )


          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dust(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( dust_av ) )  THEN
                ALLOCATE( dust_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                dust_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dust_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          found = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'silt_xy', 'silt_xz', 'silt_yz' )
!
!--       Calculate diagnostic quantities, if not done so far.
          IF ( .NOT. silt_calculated )  CALL det_actions( 'calculate silt' )

          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = silt(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( silt_av ) )  THEN
                ALLOCATE( silt_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                silt_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = silt_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          found = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

 END SUBROUTINE det_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_data_output_3d( av, variable, found, local_pf, resorted, nzb_do, nzt_do )

    USE averaging,                                                                                 &
        ONLY:  clay_av,                                                                            &
               dust_av,                                                                            &
               silt_av


    CHARACTER (LEN=*), INTENT(IN) ::  variable  !<  variable name

    INTEGER(iwp) ::  char_to_int  !< for converting character to integer and index for dust size bin
    INTEGER(iwp) ::  i            !< grid index along x-direction
    INTEGER(iwp) ::  j            !< grid index along y-direction
    INTEGER(iwp) ::  k            !< grid index along z1-direction

    INTEGER(iwp), INTENT(IN) ::  av      !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do  !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN) ::  nzt_do  !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL, INTENT(INOUT) ::  found     !< flag if output variable is found
    LOGICAL, INTENT(INOUT) ::  resorted  !< flag if output is resorted

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf


    IF ( time_since_reference_point < det_start_time )  RETURN

    found    = .TRUE.
    resorted = .TRUE.

    IF ( variable(6:11) ==  'mc_bin' )  THEN
       READ( variable(12:),* ) char_to_int
       IF ( av == 0 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = dm(char_to_int)%conc(k,j,i)
                 ENDDO
              ENDDO
           ENDDO
       ELSE
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = dm(char_to_int)%conc_av(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       found    = .TRUE.
       resorted = .TRUE.
    ENDIF

    IF ( found ) RETURN

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'clay' )
!
!--       Calculate diagnostic quantities, if not done so far.
          IF ( .NOT. clay_calculated )  CALL det_actions( 'calculate clay' )

          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = clay(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( clay_av ) )  THEN
                ALLOCATE( clay_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                clay_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = clay_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'dust' )
!
!--       Calculate diagnostic quantities, if not done so far.
          IF ( .NOT. dust_calculated )  CALL det_actions( 'calculate dust' )

          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dust(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( dust_av ) )  THEN
                ALLOCATE( dust_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                dust_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = dust_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'silt' )
!
!--       Calculate diagnostic quantities, if not done so far.
          IF ( .NOT. silt_calculated )  CALL det_actions( 'calculate silt' )

          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = silt(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( silt_av ) )  THEN
                ALLOCATE( silt_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                silt_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = silt_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          found    = .FALSE.
          resorted = .FALSE.

    END SELECT

 END SUBROUTINE det_data_output_3d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining the grid on which netcdf output variables are defined. Same grid as for other
!> scalars (see netcdf_interface_mod.f90)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    CHARACTER(LEN=*), INTENT(IN)  ::  var     !<

    CHARACTER(LEN=*), INTENT(OUT) ::  grid_x  !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_y  !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_z  !<

    LOGICAL, INTENT(OUT) ::  found   !<


    found  = .TRUE.

    SELECT CASE ( TRIM( var ) )

       CASE ( 'clay', 'clay_xy', 'clay_xz', 'clay_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
          RETURN

       CASE ( 'dust', 'dust_xy', 'dust_xz', 'dust_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
          RETURN

       CASE ( 'silt', 'silt_xy', 'silt_xz', 'silt_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
          RETURN

       CASE DEFAULT
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT
!
!-- Check for dust variables.
    IF ( var(1:7) == 'dust_mc' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu'
    ELSEIF ( var(1:15) == 'dust_emis_flux*' ) THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu1'
    ELSEIF ( var(1:15) == 'dust_depo_flux*' ) THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu1'    
    ELSE
       found  = .FALSE.
       grid_x = 'none'
       grid_y = 'none'
       grid_z = 'none'
    ENDIF

 END SUBROUTINE det_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting of cyclic lateral
!> boundary conditions for the total domain.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_exchange_horiz( location )

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< location string, describes position of the call

    INTEGER(iwp) ::  id  !< loop variable for dust size bins


    IF ( time_since_reference_point < det_start_time )  RETURN

    SELECT CASE ( location )

       CASE ( 'before_prognostic_equation' )

       CASE ( 'after_prognostic_equation' )

          DO  id = 1, n_dust_bins
             CALL exchange_horiz( dm(id)%conc, nbgp )
          ENDDO

    END SELECT

END SUBROUTINE det_exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate divergence of the gravitational settling flux causing a change in concentration
!> (kg/m³). The gravitational settling flux is defined on the zw-grid. Also v_grav was calculated
!> for the zw-grid. For k=1, the flux at the bottom is always zero because v_grav is zero there.
!> Thus, the concentration always increases at zu(1) due to the gravitational settling flux. Only
!> the deposition can decrease concentration at zu(1).
!> Vector version.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_gravitational_settling( conc, id )

    USE arrays_3d,                                                                                 &
        ONLY:  ddzw,                                                                               &
               tend

    INTEGER(iwp) ::  i  !< loop index x direction
    INTEGER(iwp) ::  j  !< loop index y direction
    INTEGER(iwp) ::  k  !< loop index z direction

    INTEGER(iwp), INTENT(IN) ::  id  !< dust bin index

    REAL(wp) ::  fluxdiv_gs  !< divergence of gravitational settling flux (sedimentation flux)

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) :: conc  !< concentration field


    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             fluxdiv_gs = ( 0.5_wp * ( conc(k+1,j,i) + conc(k,j,i)   ) * v_grav(k,id) -            &
                            0.5_wp * ( conc(k,j,i)   + conc(k-1,j,i) ) * v_grav(k-1,id) ) *        &
                          ddzw(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--          Update of tendency-term.
             tend(k,j,i) = tend(k,j,i) + fluxdiv_gs
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE det_gravitational_settling


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate divergence of the gravitational settling flux causing a change in concentration
!> (kg/m³). The gravitational settling flux is defined on the zw-grid. Also v_grav was calculated
!> for the zw-grid. For k=1, the flux at the bottom is always zero because v_grav is zero there.
!> Thus, the concentration always increases at zu(1) due to the gravitational settling flux. Only
!> the deposition can decrease concentration at zu(1).
!> Cache version.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_gravitational_settling_ij( i, j, conc, id )

    USE arrays_3d,                                                                                 &
        ONLY:  ddzw,                                                                               &
               tend


    INTEGER(iwp) ::  k   !< loop index z direction

    INTEGER(iwp), INTENT(IN) ::  i   !< loop index x direction
    INTEGER(iwp), INTENT(IN) ::  id  !< dust bin index
    INTEGER(iwp), INTENT(IN) ::  j   !< loop index y direction

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr) ::  fluxdiv_gs !< divergence of gravitational settling flux (sedimentation flux)

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) :: conc  !< concentration field


    DO  k = nzb+1, nzt
       fluxdiv_gs(k,j,i) = ( 0.5_wp * ( conc(k+1,j,i) + conc(k,j,i)   ) * v_grav(k,id) -           &
                             0.5_wp * ( conc(k,j,i)   + conc(k-1,j,i) ) * v_grav(k-1,id) ) *       &
                           ddzw(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--    Update of tendency-term.
       tend(k,j,i) = tend(k,j,i) + fluxdiv_gs(k,j,i)
    ENDDO

 END SUBROUTINE det_gravitational_settling_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for det parameters
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_header( io )

    INTEGER(iwp) ::  id  !< loop variable for dust size bins
    INTEGER(iwp) ::  is  !< loop variable for saltation size bins

    INTEGER(iwp), INTENT(IN) ::  io   !< Unit of the output file


!
!-- Write det header
    WRITE( io, 1 )
    WRITE( io, 2 ) det_start_time
    WRITE( io, 3 ) deposition_scheme
    WRITE( io, 4, ADVANCE='NO' ) n_saltation_bins
    DO is = 1, n_saltation_bins
       WRITE( io, 5, ADVANCE='NO' ) bin_mass_fraction_ssc(is)
    END DO
    WRITE( io, 6, ADVANCE='NO' )
    DO is = 1, n_saltation_bins
       WRITE( io, 7, ADVANCE='NO' ) 1.0E6*diameter_saltation(is)
    END DO
    WRITE( io, 8, ADVANCE='NO' )
    DO is = 1, n_saltation_bins
       WRITE( io, 9, ADVANCE='NO' ) mass_fraction_ssc(is)
    END DO
    WRITE( io, 10, ADVANCE='NO' )
    DO is = 1, n_saltation_bins
       WRITE( io, 11, ADVANCE='NO' ) particle_density_saltation(is)
    END DO
    WRITE( io, 12, ADVANCE='NO' ) n_dust_bins
    DO id = 1, n_dust_bins
       WRITE( io, 13, ADVANCE='NO' ) 1.0E6*lower_bound_diameter(id)
    END DO
    WRITE( io, 14, ADVANCE='NO' )
    DO id = 1, n_dust_bins
       WRITE( io, 15, ADVANCE='NO' ) 1.0E6*diameter_dust(id)
    END DO
    WRITE( io, 16, ADVANCE='NO' )
    DO id = 1, n_dust_bins
       WRITE( io, 17, ADVANCE='NO' ) 1.0E6*upper_bound_diameter(id)
    END DO
    WRITE( io, 18, ADVANCE='NO' )
    DO id = 1, n_dust_bins
       WRITE( io, 19, ADVANCE='NO' ) particle_density_dust(id)
    END DO
    WRITE( io, 20 ) TRIM( bc_dm_b), TRIM( bc_dm_t), TRIM( bc_dm_s),            &
                    TRIM( bc_dm_n), TRIM( bc_dm_l), TRIM( bc_dm_r)

!
!-- Format specifications
1   FORMAT (//' DET settings:'/' ------------'/)
2   FORMAT ('    Module starts at:  det_start_time = ', F10.2, ' s')
3   FORMAT ('    Deposition_scheme: deposition_scheme = ', A4 /)
4   FORMAT ('    Saltation bin settings:' //                                   &
            '       Number of bins: n_saltation_bins = ', I2, //               &
            '       Mass fractions for each soil separate class:   ')
5   FORMAT (F6.4, ',', 1X)
6   FORMAT (/'       Effective diameters (µm):                   ')
7   FORMAT (F7.2, ',', 1X)
8   FORMAT (/'       Mass fractions of soil separate classes:       ')
9   FORMAT (F4.2, ',', 1X)
10  FORMAT (/'       Particle density (kg m-3):                     ')
11  FORMAT (F6.1, ',', 1X)
12  FORMAT (//'    Dust bin settings:' //                                      &
            '       Number of bins: n_dust_bins = ', I2, //                    &
            '       Lower diameter for each bin: ')
13  FORMAT (F5.2, ',', 1X)
14  FORMAT (/'       Effective diameters (µm):    ')
15  FORMAT (F5.2, ',', 1X)
16  FORMAT (/'       Upper diameter for each bin: ')
17  FORMAT (F5.2, ',', 1X)
18  FORMAT (/'       Particle density (kg m-3):    ')
19  FORMAT (F6.1, ',', 1X)
20  FORMAT (//'    Boundary conditions for dust mass:' //                      &
              '       bottom/top:  ', 2(A20) /                                 &
              '       north/south: ', 2(A20) /                                 &
              '       left/right:  ', 2(A20))

 END SUBROUTINE det_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the dust emission and transport model (DET).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_init

    INTEGER(iwp) ::  i   !< grid index in x-direction
    INTEGER(iwp) ::  id  !< loop variable for dust size bins
    INTEGER(iwp) ::  is  !< loop variable for saltation size bins
    INTEGER(iwp) ::  j   !< grid index in y-direction
    INTEGER(iwp) ::  k   !< grid index in z-direction
    INTEGER(iwp) ::  m   !< running index for surface elements

    REAL(wp) ::  cc              !< Cunningham slip-flow correction factor
    REAL(wp) ::  kn              !< Knudsen number
    REAL(wp) ::  lambda_f        !< molecular mean free path (m)
    REAL(wp) ::  ln_arg          !< natural logarithm argument in normalized volume distribution
    REAL(wp) ::  n_sfc = 0.0_wp  !< total basal surface area of soil bed
    REAL(wp) ::  n_v = 0.0_wp    !< total normalized volume distribution of emitted dust
    REAL(wp) ::  v_th            !< thermal speed of an air molecule


    IF ( debug_output )  CALL debug_message( 'det_init', 'start' )

!
!-- If parameters of salatation bins are not defined by the user they are set according to LeGrand
!-- et al. (2019), Table A1. The mass fractions of the soil separate classes clay, silt, and sand
!-- (mass_fraction_ssc) are based on the soil type sand (Pérez et al., 2011).
    IF ( ALL( bin_mass_fraction_ssc == not_set ) )  THEN
       bin_mass_fraction_ssc(1:n_saltation_bins) = (/ 1.0_wp, 0.25_wp, 0.25_wp, 0.25_wp, 0.25_wp,  &
                                                      0.0205_wp, 0.0410_wp, 0.0359_wp, 0.3897_wp,  &
                                                      0.5128_wp /)
    ENDIF

    IF ( ALL( diameter_saltation == not_set ) )  THEN
       diameter_saltation(1:n_saltation_bins) = (/   1.42E-6_wp,   8.0E-6_wp,  20.0E-6_wp,         &
                                                     32.0E-6_wp,  44.0E-6_wp,  70.0E-6_wp,         &
                                                    130.0E-6_wp, 200.0E-6_wp, 620.0E-6_wp,         &
                                                   1500.0E-6_wp /)
    ENDIF

    IF ( ALL( mass_fraction_ssc == not_set ) )  THEN
       mass_fraction_ssc(1:n_saltation_bins) = (/ 0.03_wp, 0.05_wp, 0.05_wp, 0.05_wp, 0.05_wp,     &
                                                  0.92_wp, 0.92_wp, 0.92_wp, 0.92_wp, 0.92_wp /)
    ENDIF

    IF ( ALL( particle_density_saltation == not_set ) )  THEN
       particle_density_saltation(1:n_saltation_bins) = (/ 2500.0_wp, 2650.0_wp, 2650.0_wp,        &
                                                           2650.0_wp, 2650.0_wp, 2650.0_wp,        &
                                                           2650.0_wp, 2650.0_wp, 2650.0_wp,        &
                                                           2650.0_wp /)
    ENDIF

    sb_prop(:)%diameter              = diameter_saltation(1:n_saltation_bins)
    sb_prop(:)%bin_mass_fraction_ssc = bin_mass_fraction_ssc(1:n_saltation_bins)
    sb_prop(:)%mass_fraction_ssc     = mass_fraction_ssc(1:n_saltation_bins)
    sb_prop(:)%density               = particle_density_saltation(1:n_saltation_bins)
!
!-- If parameters of dust bins are not defined by the user they are set according to LeGrand
!-- et al. (2019), Table 2.
    IF ( ALL( lower_bound_diameter == not_set ) )  THEN
       lower_bound_diameter(1:n_dust_bins) = (/ 0.2E-6_wp, 2.0E-6_wp, 3.6E-6_wp, 6.0E-6_wp,         &
                                               12.0E-6_wp /)
     ENDIF

     IF ( ALL( diameter_dust == not_set ) )  THEN
        diameter_dust(1:n_dust_bins) = (/ 1.46E-6_wp, 2.8E-6_wp, 4.8E-6_wp, 9.0E-6_wp, 16.0E-6_wp /)
     ENDIF

     IF ( ALL( upper_bound_diameter == not_set ) )  THEN
        upper_bound_diameter(1:n_dust_bins) = (/ 2.0E-6_wp, 3.6E-6_wp, 6.0E-6_wp, 12.0E-6_wp,      &
                                                20.0E-6_wp /)
     ENDIF

     IF ( ALL( particle_density_dust == not_set ) )  THEN
        particle_density_dust(1:n_dust_bins) = (/ 2500.0_wp, 2650.0_wp, 2650.0_wp, 2650.0_wp,      &
                                                 2650.0_wp /)
     ENDIF

     db_prop(:)%diameter             = diameter_dust(1:n_dust_bins)
     db_prop(:)%lower_bound_diameter = lower_bound_diameter(1:n_dust_bins)
     db_prop(:)%upper_bound_diameter = upper_bound_diameter(1:n_dust_bins)
     db_prop(:)%density              = particle_density_dust(1:n_dust_bins)

!
!-- Calculate size-resolved basal surface coverage fractions in m²/kg, LeGrand et al. (2019) Eq. 11.
    DO  is = 1, n_saltation_bins
       dm_rel(is) = sb_prop(is)%bin_mass_fraction_ssc * sb_prop(is)%mass_fraction_ssc
       ds_sfc(is) = dm_rel(is) / ( 2.0_wp / 3.0_wp * sb_prop(is)%density * sb_prop(is)%diameter )
       n_sfc      = n_sfc + ds_sfc(is)
    ENDDO
!
!-- Calculate salatation bin weighting factors, LeGrand et al. (2019) Eq. 12.
    DO  is = 1, n_saltation_bins
       ds_rel(is) = ds_sfc(is) / n_sfc
    ENDDO
!
!-- Calculate sandblasting efficiency, Marticorena and Bergametti (1995) Eq. 47.
!-- Factor 100 due to conversion from 1/cm to 1/m.
    alpha_s = 100.0_wp * ( 10.0_wp ** (0.134_wp * sb_prop(1)%mass_fraction_ssc - 6.0_wp) )
!
!-- Calculate dynamic viscosity of air in kg/(ms). Deviating from Klamt et al. (2024), a constant
!-- value of the viscosity is used, based on pt_surface.
    air_viscosity = 1.8325E-5_wp * ( 416.16_wp / ( pt_surface + 120.0_wp ) ) *                     &
                                   ( pt_surface / 296.16_wp )**1.5
!
!-- Calculate gravitational settling velocities for each dust bin:
!-- Thermal velocity of an air molecule, Eq. 15.32 (Cunningham slip-flow correction).
    v_th = SQRT( 8.0_wp * k_boltzmann * pt_surface / ( pi * am_airmol ) )
    DO  id = 1, n_dust_bins

       v_grav(nzb,id) = 0.0_wp

       DO  k = nzb+1, nzt
!
!--       Mean free path in m, Eq. 15.24
          lambda_f = 2.0_wp * air_viscosity / ( rho_air_zw(k) * v_th )
!
!--       Knudsen number, Eq. 15.23
          kn = lambda_f / ( db_prop(id)%diameter * 0.5_wp )
!
!--       Cunningham slip-flow correction, Eq. 15.30
          cc = 1.0_wp + kn * ( 1.249_wp + 0.42_wp * EXP( -0.87_wp / kn ) )
!
!--       Critical fall speed, i.e., gravitational settling velocity, Eq. 20.4
          v_grav(k,id) = db_prop(id)%diameter**2 * ( db_prop(id)%density - rho_air_zw(k) ) * g *   &
                         cc / ( 18.0_wp * air_viscosity )
       ENDDO

    ENDDO

!
!-- Calculate bin-specific threshold friction velocity. Because friction velocity is defined in the 
!-- grid box center and not at the surface location, the index k of the surface element is used, 
!-- which directly refers to the first atmospheric grid point above the upward-facing surface.
    DO  is = 1, n_saltation_bins
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Default type surfaces.
             DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                IF ( surf_def%upward(m) )  THEN
                   k = surf_def%k(m)
                   surf_def%us_t(m,is) = 0.129_wp * ( SQRT( sb_prop(is)%density * g *              &
                                                            sb_prop(is)%diameter / rho_air(k) )    &
                                                    * SQRT( 1.0_wp + 6.0E-7_wp /                   &
                                                            ( sb_prop(is)%density * g *            &
                                                              sb_prop(is)%diameter**2.5_wp ) )     &
                                                    ) /                                            &
                                         SQRT( 1.928_wp *                                          &
                                               ( 1.75E6_wp * sb_prop(is)%diameter**1.56_wp         &
                                               + 0.38_wp )**0.092 - 1.0_wp                         &
                                             )
                ENDIF
             ENDDO
!
!--          Natural type surfaces. surf_usm and surf_top are not realized so far.
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m) )  THEN
                   k = surf_lsm%k(m)
                   surf_lsm%us_t(m,is) = 0.129_wp * ( SQRT( sb_prop(is)%density * g *              &
                                                            sb_prop(is)%diameter / rho_air(k) )    &
                                                    * SQRT( 1.0_wp + 6.0E-7_wp /                   &
                                                            ( sb_prop(is)%density * g *            &
                                                              sb_prop(is)%diameter**2.5_wp ) )     &
                                                    ) /                                            &
                                         SQRT( 1.928_wp *                                          &
                                               ( 1.75E6_wp * sb_prop(is)%diameter**1.56_wp         &
                                               + 0.38_wp )**0.092 - 1.0_wp                         &
                                             )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate suspended dust distribution weighting factors, LeGrand et al. (2019) Eq. 15.
!-- First, calculate normalized volume size distributions for each dust size bin.
!-- For empirical constants see LeGrand et al. (2019).
    DO  id = 1, n_dust_bins
       ln_arg = LOG( db_prop(id)%diameter / 3.4E-6_wp ) / ( SQRT( 2.0_wp ) * LOG( 3.0_wp ) )
       dv(id) = db_prop(id)%diameter / 12.62E-6_wp * ( 1.0_wp + ERF( ln_arg ) ) *                  &
                EXP( -( db_prop(id)%diameter / 12.0E-6_wp )**3 ) *                                 &
                LOG( db_prop(id)%upper_bound_diameter / db_prop(id)%lower_bound_diameter )
       n_v = n_v + dv(id)
    ENDDO
!
!-- Second, calculate dust distribution weighting factors,
    DO  id = 1, n_dust_bins
       dw(id) = dv(id) / n_v
    ENDDO

    IF ( debug_output )  CALL debug_message( 'det_init', 'end' )

 END SUBROUTINE det_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate and initialize arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_init_arrays

    USE pegrid,                                                                                    &
        ONLY:  threads_per_task


    INTEGER(iwp) ::  id  !< loop variable for dust size bins


    IF ( debug_output )  CALL debug_message( 'det_init_arrays', 'start' )
!
!-- General array allocation and initialization.
    ALLOCATE( dm_rel(n_saltation_bins) )
    ALLOCATE( ds_sfc(n_saltation_bins) )
    ALLOCATE( ds_rel(n_saltation_bins) )
    ALLOCATE( db_prop(n_dust_bins) )
    ALLOCATE( dv(n_dust_bins) )
    ALLOCATE( dw(n_dust_bins) )
    ALLOCATE( sb_prop(n_saltation_bins) )

    dm_rel(:) = 0.0_wp
    ds_sfc(:) = 0.0_wp
    ds_rel(:) = 0.0_wp
    dv(:)     = 0.0_wp
    dw(:)     = 0.0_wp

!
!-- Array allocation and initialization for prognostic variables and fluxes.
    ALLOCATE( dust_fluxes(n_dust_bins) )
    ALLOCATE( dm(n_dust_bins) )
    ALLOCATE( v_grav(nzb:nzt,n_dust_bins) )
    ALLOCATE( dm_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,n_dust_bins) )
    ALLOCATE( dm_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,n_dust_bins) )
    ALLOCATE( dm_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,n_dust_bins) )

    dm_1 = 0.0_wp
    dm_2 = 0.0_wp
    dm_3 = 0.0_wp

    DO  id = 1, n_dust_bins
       dm(id)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => dm_1(:,:,:,id)
       dm(id)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => dm_2(:,:,:,id)
       dm(id)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => dm_3(:,:,:,id)

       ALLOCATE( dm(id)%flux_s(nzb+1:nzt,0:threads_per_task-1) )
       ALLOCATE( dm(id)%diss_s(nzb+1:nzt,0:threads_per_task-1) )
       ALLOCATE( dm(id)%flux_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
       ALLOCATE( dm(id)%diss_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
       ALLOCATE( dm(id)%init(nzb:nzt+1) )

       dm(id)%flux_s  = 0.0_wp
       dm(id)%diss_s  = 0.0_wp
       dm(id)%flux_l  = 0.0_wp
       dm(id)%diss_l  = 0.0_wp
       dm(id)%init    = 0.0_wp
    ENDDO
!
!-- Surface-related data: dm = dust mass.
!-- Default type surfaces.
    ALLOCATE( surf_def%dmsws(surf_def%ns,n_dust_bins) )
    ALLOCATE( surf_def%dm_depo_flux(n_dust_bins,surf_def%ns) )
    ALLOCATE( surf_def%dm_emis_flux(n_dust_bins,surf_def%ns) )
    ALLOCATE( surf_def%us_t(surf_def%ns,n_saltation_bins) )

    surf_def%dm_depo_flux = 0.0_wp
    surf_def%dm_emis_flux = 0.0_wp
    surf_def%dmsws = 0.0_wp
    surf_def%us_t  = 0.0_wp
!
!-- Natural and urban type surfaces. Note, urban type surface and model top 
!-- surfaces are not able to release dust so far. However, they need to be 
!-- allocated to enable the CALL of diffusion_s
    ALLOCATE( surf_lsm%dmsws(surf_lsm%ns,n_dust_bins) )
    ALLOCATE( surf_lsm%dm_depo_flux(n_dust_bins,surf_lsm%ns) )
    ALLOCATE( surf_lsm%dm_emis_flux(n_dust_bins,surf_lsm%ns) )
    ALLOCATE( surf_lsm%us_t(surf_lsm%ns,n_saltation_bins) )

    surf_def%dm_depo_flux = 0.0_wp
    surf_def%dm_emis_flux = 0.0_wp
    surf_lsm%dmsws = 0.0_wp
    surf_lsm%us_t  = 0.0_wp

    ALLOCATE( surf_usm%dmsws(surf_usm%ns,n_dust_bins) )

    surf_usm%dmsws = 0.0_wp
!
!-- Model top surfaces
    ALLOCATE( surf_top%dmsws(surf_top%ns,n_dust_bins) )

    surf_top%dmsws = 0.0_wp

    IF ( debug_output )  CALL debug_message( 'det_init_arrays', 'end' )

 END SUBROUTINE det_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Contains processes like deposition that are calculated before the actual prognostic equations
!> (see det_prognostic_equations_ij) similar to salsa_non_advective_processes_ij and
!> chem_non_advective_processes_ij. The aim is to calculate a surface dust mass net flux for
!> the prognostic equations: surf%dmsws(m,id).
!> Vector version.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_non_advective_processes

    INTEGER(iwp) ::  i     !< grid index in x-direction
    INTEGER(iwp) ::  j     !< grid index in y-direction
    INTEGER(iwp) ::  k     !< index for z-direction
    INTEGER(iwp) ::  koff  !< offset index for z-direction
    INTEGER(iwp) ::  m     !< running index for surface elements


    IF ( time_since_reference_point < det_start_time )  RETURN

!
!-- Start the calculations after a certain simulation time
    CALL cpu_log( log_point_s(106), 'dust emission', 'start' )

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Start with default type surfaces.
          IF ( surf_def%ns >= 1 )  THEN
             surf => surf_def
             DO  m = surf%start_index(j,i), surf%end_index(j,i)
!
!--             Only upward-facing surfaces are considered for dust release and deposition.
                IF ( surf%upward(m) )  THEN
!
!--                Initialize local loop variables. k(m) refers to first atmospheric grid point
!--                above the upward-facing surface at zu(k). Here, koff(m) is -1 , i.e, k + koff
!--                refers to the surface position index. The surface height is at zw(k + koff).
                   k    = surf%k(m)
                   koff = surf%koff(m)
!
!--                Calculation of bin-specific dust emission.
                   CALL det_calculate_emission( rho_air_zw(k+koff), m, surf%dm_emis_flux(:,m) )
!
!--                Calculation of bin-specific deposition flux.
                   IF ( deposition_scheme == 'Z01' )  THEN
                      CALL det_calculate_dry_deposition_z01( i, j, k, m, surf%dm_depo_flux(:,m) )
                   ENDIF
!
!--                Calculate surface flux of dust mass. Because later on the dynamic flux is
!--                required, e.g., in the diffusion scheme, multiply with the density
!--                ( kg/(m2*s) * kg/m3 ), see also chem_emissions_mod.f90.
                   surf%dmsws(m,:) = ( surf%dm_emis_flux(:,m) + surf%dm_depo_flux(:,m) ) *         &
                                     rho_air_zw(k+koff)
                ENDIF
             ENDDO
          ENDIF
!
!--       Natural type surfaces.
          IF ( surf_lsm%ns >= 1 )  THEN
             surf => surf_lsm
             DO  m = surf%start_index(j,i), surf%end_index(j,i)
!
!--             Only upward-facing surfaces are considered for dust release and deposition.
                IF ( surf%upward(m) )  THEN
!
!--                Initialize local loop variables. k(m) refers to first atmospheric grid point
!--                above the upward-facing surface at zu(k). Here, koff(m) is -1 , i.e, k + koff
!--                refers to the surface position index. The surface height is at zw(k + koff).
                   k    = surf%k(m)
                   koff = surf%koff(m)
!
!--                Calculation of bin-specific dust emission.
                   CALL det_calculate_emission( rho_air_zw(k+koff), m, surf%dm_emis_flux(:,m) )
!
!--                Calculation of bin-specific deposition flux.
                   IF ( deposition_scheme == 'Z01' )  THEN
                      CALL det_calculate_dry_deposition_z01( i, j, k, m, surf%dm_depo_flux(:,m) )
                   ENDIF
!
!--                Calculate surface flux of dust mass. Because later on the dynamic flux is
!--                require, e.g., in the diffusion scheme, multiply with the density
!--                ( kg/(m2*s) * kg/m3 ), see also chem_emissions_mod.f90.
                   surf%dmsws(m,:) = ( surf%dm_emis_flux(:,m) + surf%dm_depo_flux(:,m) ) *         &
                                     rho_air_zw(k+koff)
                ENDIF
             ENDDO
          ENDIF

       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(106), 'dust emission', 'stop' )

 END SUBROUTINE det_non_advective_processes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Contains processes like deposition that are calculated before the actual prognostic equations
!> (see det_prognostic_equations_ij) similar to salsa_non_advective_processes_ij and
!> chem_non_advective_processes_ij. The aim is to calculate a surface dust mass net flux for
!> the prognostic equations: surf%dmsws(m,id).
!> Cache version.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_non_advective_processes_ij( i, j )

    INTEGER(iwp) ::  k     !< index for z-direction
    INTEGER(iwp) ::  koff  !< offset index for z-direction
    INTEGER(iwp) ::  m     !< running index for surface elements

    INTEGER(iwp), INTENT(IN) ::  i  !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j  !< grid index in y-direction


    IF ( time_since_reference_point < det_start_time )  RETURN

!
!-- Start the calculations after a certain simulation time
!
!-- Start with default type surfaces.
    IF ( surf_def%ns >= 1 )  THEN
       surf => surf_def
       DO  m = surf%start_index(j,i), surf%end_index(j,i)
!
!--       Only upward-facing surfaces are considered for dust release and deposition.
          IF ( surf%upward(m) )  THEN
!
!--          Initialize local loop variables. k(m) refers to first atmospheric grid point
!--          above the upward-facing surface at zu(k). Here, koff(m) is -1 , i.e, k + koff
!--          refers to the surface position index. The surface height is at zw(k + koff).
             k    = surf%k(m)
             koff = surf%koff(m)
!
!--          Calculation of bin-specific dust emission.
             CALL det_calculate_emission( rho_air_zw(k+koff), m, surf%dm_emis_flux(:,m) )
!
!--          Calculation of bin-specific deposition flux.
             IF ( deposition_scheme == 'Z01' )  THEN
                CALL det_calculate_dry_deposition_z01( i, j, k, m, surf%dm_depo_flux(:,m) )
             ENDIF
!
!--          Calculate surface flux of dust mass. Because later on the dynamic flux is
!--          required, e.g., in the diffusion scheme, multiply with the density ( kg/(m2*s) *
!--          kg/m3 ), see also chem_emissions_mod.f90.
             surf%dmsws(m,:) = ( surf%dm_emis_flux(:,m) + surf%dm_depo_flux(:,m) ) *               &
                               rho_air_zw(k+koff)
          ENDIF
       ENDDO
    ENDIF
!
!-- Natural type surfaces.
    IF ( surf_lsm%ns >= 1 )  THEN
       surf => surf_lsm
       DO  m = surf%start_index(j,i), surf%end_index(j,i)
!
!--       Only upward-facing surfaces are considered for dust release and deposition.
          IF ( surf%upward(m) )  THEN
!
!--          Initialize local loop variables. k(m) refers to first atmospheric grid point
!--          above the upward-facing surface at zu(k). Here, koff(m) is -1 , i.e, k + koff
!--          refers to the surface position index. The surface height is at zw(k + koff).
             k    = surf%k(m)
             koff = surf%koff(m)
!
!--          Calculation of bin-specific dust emission.
             CALL det_calculate_emission( rho_air_zw(k+koff), m, surf%dm_emis_flux(:,m) )
!
!--          Calculation of bin-specific deposition flux.
             IF ( deposition_scheme == 'Z01' )  THEN
                CALL det_calculate_dry_deposition_z01( i, j, k, m, surf%dm_depo_flux(:,m) )
             ENDIF
!
!--          Calculate surface flux of dust mass. Because later on the dynamic flux is
!--          require, e.g., in the diffusion scheme, multiply with the density ( kg/(m2*s) *
!--          kg/m3 ), see also chem_emissions_mod.f90.
             surf%dmsws(m,:) = ( surf%dm_emis_flux(:,m) + surf%dm_depo_flux(:,m) ) *               &
                               rho_air_zw(k+koff)
          ENDIF
       ENDDO
    ENDIF

 END SUBROUTINE det_non_advective_processes_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read namelist &det_parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file


    NAMELIST /det_parameters/  alpha_imp,                                                          &
                               bc_dm_b,                                                            &
                               bc_dm_l,                                                            &
                               bc_dm_n,                                                            &
                               bc_dm_r,                                                            &
                               bc_dm_s,                                                            &
                               bc_dm_t,                                                            &
                               bin_mass_fraction_ssc,                                              &
                               brownian_diffusion_coefficient,                                     &
                               deposition_scheme,                                                  &
                               det_start_time,                                                     &
                               diameter_dust,                                                      &
                               diameter_saltation,                                                 &
                               lower_bound_diameter,                                               &
                               mass_fraction_ssc,                                                  &
                               n_dust_bins,                                                        &
                               n_saltation_bins,                                                   &
                               particle_density_dust,                                              &
                               particle_density_saltation,                                         &
                               switch_off_module,                                                  &
                               upper_bound_diameter

!
!-- Move to the beginning of the namelist file and try to find and read the user-defined namelist
!-- det_parameters.
    REWIND( 11 )
    READ( 11, det_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    det_parameters namelist was found and read correctly. Set flag that indicates that
!--    the dust emission and transport module (DET) is switched on.
       IF ( .NOT. switch_off_module )  det_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    det_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'det_parameters', line )

    ENDIF

 END SUBROUTINE det_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate prognostic equation for dust mass concentration. Vector optimized version, i.e., call
!> for all grid points. Not implemented yet.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_prognostic_equations

    USE advec_s_pw_mod,                                                                            &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                                            &
        ONLY:  advec_s_up

    USE advec_ws,                                                                                  &
        ONLY:  advec_s_ws

    USE arrays_3d,                                                                                 &
        ONLY:  rdf_sc,                                                                             &
               tend

    USE control_parameters,                                                                        &
        ONLY:  advanced_div_correction,                                                            &
               allow_negative_scalar_values,                                                       &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               dt_3d,                                                                              &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               use_subsidence_tendencies,                                                          &
               timestep_scheme,                                                                    &
               tsc,                                                                                &
               ws_scheme_sca,                                                                      &
               large_scale_forcing,                                                                &
               large_scale_subsidence

    USE diffusion_s_mod,                                                                           &
        ONLY:  diffusion_s

    USE indices,                                                                                   &
        ONLY:  advc_flags_s

    USE subsidence_mod,                                                                            &
        ONLY:  subsidence


    INTEGER(iwp) ::  i   !< loop index x direction
    INTEGER(iwp) ::  id  !< loop variable for dust size bins
    INTEGER(iwp) ::  j   !< loop index y direction
    INTEGER(iwp) ::  k   !< loop index z direction


    IF ( time_since_reference_point < det_start_time )  RETURN

    CALL cpu_log( log_point(68), 'dm-equation', 'start' )

!
!-- Reset flag that diagnostic det quantities have not been calculated so far.
    clay_calculated = .FALSE.
    dust_calculated = .FALSE.
    silt_calculated = .FALSE.
!
!-- Advective processes for all dust size bins id.
    DO  id = 1, n_dust_bins
!
!--    Tendency-terms for dust mass inside size bin.
       tend = 0.0_wp
!
!--    Advection terms.
!--    TODO: Standard flags indicating dirichlet/radiation boundary conditions are used here.
!--          However, as soons as (horizontal) boundary conditions for prognostic dust variables
!--          can be different to PALM's standard quantities, new flags need to be introduced
!--          (see SALSA or CHEM).
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( ws_scheme_sca )  THEN
             IF ( .NOT. advanced_div_correction )  THEN
                CALL advec_s_ws( advc_flags_s, dm(id)%conc, 'dm',                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_ws( advc_flags_s, dm(id)%conc, 'dm',                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s,                             &
                                 advanced_div_correction )
             ENDIF
          ELSE
             CALL advec_s_pw( dm(id)%conc )
          ENDIF
       ELSE
          CALL advec_s_up( dm(id)%conc )
       ENDIF
!
!--    Diffusion terms. Note, because urban and top surfaces can not emit dust so far, these
!--    surface fluxes do not have an effect.
       CALL diffusion_s( dm(id)%conc, surf_top%dmsws(:,id),                                        &
                         surf_def%dmsws(:,id), surf_lsm%dmsws(:,id),                               &
                         surf_usm%dmsws(:,id) )
!
!--    If required compute influence of large-scale subsidence/ascent. Note, the last argument
!--    (lsf_index) is of no meaning in this case, as it is only used in conjunction with
!--    large_scale_forcing, which is to date not implemented for scalars like the dust mass,
!--    similar to the nudging. Because lsf_index will be an optional parameter in the future,
!--    subsidence can then be called without the last integer argument and the check for .NOT.
!--    large_scale_forcing can then be removed.
       IF ( large_scale_subsidence  .AND.  .NOT. use_subsidence_tendencies  .AND.                  &
            .NOT. large_scale_forcing )                                                            &
       THEN
          CALL subsidence( tend, dm(id)%conc, dm(id)%init )
       ENDIF
!
!--    Change in concentration due to the gravitational settling flux.
       CALL det_gravitational_settling( dm(id)%conc, id )
!
!--    Prognostic equation for mass concentration.
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                dm(id)%conc_p(k,j,i) = dm(id)%conc(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +     &
                                                                  tsc(3) * dm(id)%tconc_m(k,j,i)   &
                                                                      ) -                          &
                                                              tsc(5) * rdf_sc(k) *                 &
                                                           ( dm(id)%conc(k,j,i) - dm(id)%init(k) ) &
                                                            ) *                                    &
                                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 )  )
                IF ( dm(id)%conc_p(k,j,i) < 0.0_wp  .AND.  .NOT. allow_negative_scalar_values )    &
                THEN
                   dm(id)%conc_p(k,j,i) = 0.1_wp * dm(id)%conc(k,j,i)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step.
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      dm(id)%tconc_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      dm(id)%tconc_m(k,j,i) = -9.5625_wp * tend(k,j,i) +                           &
                                               5.3125_wp * dm(id)%tconc_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ENDDO

    CALL cpu_log( log_point(68), 'dm-equation', 'stop' )

END SUBROUTINE det_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate prognostic equation for dust mass concentration. Cache optimized version, i.e., call
!> for grid points i,j.
!> @todo string identifier dm (duss mass), used to assign fluxes to the correct dimension in the 
!> analysis array must be implemented in advec_ws.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_prognostic_equations_ij( i, j, i_omp_start, tn )

    USE advec_s_pw_mod,                                                                            &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                                            &
        ONLY:  advec_s_up

    USE advec_ws,                                                                                  &
        ONLY:  advec_s_ws

    USE arrays_3d,                                                                                 &
        ONLY:  rdf_sc,                                                                             &
               tend

    USE control_parameters,                                                                        &
        ONLY:  allow_negative_scalar_values,                                                       &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               dt_3d,                                                                              &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               monotonic_limiter_z,                                                                &
               use_subsidence_tendencies,                                                          &
               timestep_scheme,                                                                    &
               tsc,                                                                                &
               ws_scheme_sca,                                                                      &
               large_scale_forcing,                                                                &
               large_scale_subsidence

    USE diffusion_s_mod,                                                                           &
        ONLY:  diffusion_s

    USE indices,                                                                                   &
        ONLY:  advc_flags_s

    USE subsidence_mod,                                                                            &
        ONLY:  subsidence


    INTEGER(iwp), INTENT(IN) ::  i            !< loop index x direction
    INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in calling routine prognostic_equations
    INTEGER(iwp), INTENT(IN) ::  j            !< loop index y direction
    INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

    INTEGER(iwp) ::  id  !< loop variable for dust size bins
    INTEGER(iwp) ::  k   !< loop index z direction


    IF ( time_since_reference_point < det_start_time )  RETURN

!
!-- Reset flag that diagnostic det quantities have not been calculated so far.
    clay_calculated = .FALSE.
    dust_calculated = .FALSE.
    silt_calculated = .FALSE.
!
!-- Advective processes for all dust size bins id.
    DO  id = 1, n_dust_bins
!
!--    Tendency-terms for dust mass inside size bin.
       tend(:,j,i) = 0.0_wp
!
!--    Advection terms.
!--    TODO: Standard flags indicating dirichlet/radiation boundary conditions are used here.
!--          However, as soons as (horizontal) boundary conditions for prognostic dust variables
!--          can be different to PALM's standard quantities, new flags need to be introduced
!--          (see SALSA or CHEM).
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( ws_scheme_sca )  THEN
             CALL advec_s_ws( advc_flags_s, i, j, dm(id)%conc, 'dm',                               &
                              dm(id)%flux_s, dm(id)%diss_s,                                        &
                              dm(id)%flux_l, dm(id)%diss_l,                                        &
                              i_omp_start, tn,                                                     &
                              bc_dirichlet_l  .OR.  bc_radiation_l,                                &
                              bc_dirichlet_n  .OR.  bc_radiation_n,                                &
                              bc_dirichlet_r  .OR.  bc_radiation_r,                                &
                              bc_dirichlet_s  .OR.  bc_radiation_s,                                &
                              monotonic_limiter_z )
          ELSE
             CALL advec_s_pw( i, j, dm(id)%conc )
          ENDIF
       ELSE
          CALL advec_s_up( i, j, dm(id)%conc )
       ENDIF
!
!--    Diffusion terms. Note, because urban and top surfaces can not emit dust so far, these
!--    surface fluxes do not have an effect.
       CALL diffusion_s( i, j, dm(id)%conc, surf_top%dmsws(:,id),                                  &
                         surf_def%dmsws(:,id), surf_lsm%dmsws(:,id),                               &
                         surf_usm%dmsws(:,id) )
!
!--    If required compute influence of large-scale subsidence/ascent. Note, the last argument
!--    (lsf_index) is of no meaning in this case, as it is only used in conjunction with
!--    large_scale_forcing, which is to date not implemented for scalars like the dust mass,
!--    similar to the nudging. Because lsf_index will be an optional parameter in the future,
!--    subsidence can then be called without the last integer argument and the check for .NOT.
!--    large_scale_forcing can then be removed.
       IF ( large_scale_subsidence  .AND.  .NOT. use_subsidence_tendencies  .AND.                  &
            .NOT. large_scale_forcing )                                                            &
       THEN
          CALL subsidence( i, j, tend, dm(id)%conc, dm(id)%init )
       ENDIF
!
!--    Change in concentration due to the gravitational settling flux.
       CALL det_gravitational_settling( i, j, dm(id)%conc, id )
!
!--    Prognostic equation for mass concentration.
       DO  k = nzb+1, nzt
          dm(id)%conc_p(k,j,i) = dm(id)%conc(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +           &
                                                                  tsc(3) * dm(id)%tconc_m(k,j,i)   &
                                                                ) -                                &
                                                        tsc(5) * rdf_sc(k) *                       &
                                                        ( dm(id)%conc(k,j,i) - dm(id)%init(k) )    &
                                                      ) *                                          &
                                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 )  )
          IF ( dm(id)%conc_p(k,j,i) < 0.0_wp  .AND.  .NOT. allow_negative_scalar_values )   &
          THEN
             dm(id)%conc_p(k,j,i) = 0.1_wp * dm(id)%conc(k,j,i)
          ENDIF
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step.
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                dm(id)%tconc_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                dm(id)%tconc_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * dm(id)%tconc_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF
    ENDDO

 END SUBROUTINE det_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_rrd_global_ftn( found )

    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string

    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'n_dust_bins' )
          READ ( 13 )  n_dust_bins

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE det_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_rrd_global_mpi

    CALL rrd_mpi_io( 'n_dust_bins', n_dust_bins )

 END SUBROUTINE det_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,    &
                               nyn_on_file, nysf, nysc, nys_on_file, tmp_3d, found )

    USE control_parameters

    CHARACTER(LEN=5) ::  dust_bin_name  !< name of dust mass bin on restart file

    INTEGER(iwp) ::  id              !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT) :: found

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) &
                 :: tmp_3d   !< 3D array to temp store data


    found = .FALSE.

    IF ( ALLOCATED( dm ) )  THEN

       DO  id = 1, n_dust_bins

          WRITE( dust_bin_name, '(A3,I2.2)' )  'dm_', id

          IF ( restart_string(1:length) == dust_bin_name )  THEN

             IF ( k == 1 )  READ ( 13 )  tmp_3d
             dm(id)%conc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          ELSEIF ( restart_string(1:length) == dust_bin_name // '_av' )  THEN

             IF ( .NOT. ALLOCATED( dm(id)%conc_av ) )  THEN
                ALLOCATE( dm(id)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             dm(id)%conc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          ENDIF

       ENDDO

    ENDIF

 END SUBROUTINE det_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_rrd_local_mpi

    IMPLICIT NONE

    CHARACTER(LEN=5) ::  dust_bin_name  !< name of dust mass bin on restart file

    INTEGER(iwp) ::  id  !<

    LOGICAL      ::  array_found  !<


    DO  id = 1, n_dust_bins

       WRITE( dust_bin_name, '(A3,I2.2)' )  'dm_', id

       CALL rrd_mpi_io( dust_bin_name, dm(id)%conc )

!
!--    Restart input of time-averaged quantities is skipped in case of cyclic-fill initialization.
!--    This case, input of time-averaged data is useless and can lead to faulty averaging.
       IF ( .NOT. cyclic_fill_initialization )  THEN

          CALL rd_mpi_io_check_array( dust_bin_name // '_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( dm(id)%conc_av ) )  THEN
                ALLOCATE( dm(id)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( dust_bin_name // '_av', dm(id)%conc_av )
          ENDIF

       ENDIF

    ENDDO

 END SUBROUTINE det_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of det statistics, i.e. horizontally averaged profiles and time series.
!> This routine is called for every statistic region sr defined by the user, but at least for the
!> region "total domain" (sr=0).
!> @todo Implement time series output.
!> @todo Implement profile output for each dust-sized bin
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_statistics( mode, sr, tn )

    USE indices,                                                                                   &
        ONLY:  ngp_2dh

#if defined( __parallel )
    USE pegrid,                                                                                    &
        ONLY:  collective_wait,                                                                    &
               comm2d,                                                                             &
               ierr
#endif

    USE statistics,                                                                                &
        ONLY:  pr_palm,                                                                            &
               rmask,                                                                              &
               sums_l,                                                                             &
               ts_value

#if defined( __parallel )
    USE MPI
#endif


    CHARACTER(LEN=*), INTENT(IN) ::  mode  !<

    INTEGER(iwp) ::  i    !< grid index in x-direction
    INTEGER(iwp) ::  ii   !< loop index over det profiles
    INTEGER(iwp) ::  ind  !< index for statistical output
    INTEGER(iwp) ::  j    !< grid index in y-direction
    INTEGER(iwp) ::  k    !< grid index in y-direction
    INTEGER(iwp) ::  m    !< running index for surface elements

    INTEGER(iwp), INTENT(IN) ::  sr  !<  statistical region
    INTEGER(iwp), INTENT(IN) ::  tn  !<  thread number

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts_value_l  !< to store local maxima

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bulk_depo  !< bulk deposition flux (over all dust bins)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  bulk_emis  !< bulk emission flux (over all dust bins)


    IF ( time_since_reference_point < det_start_time )  RETURN

    IF ( mode == 'profiles' )  THEN
!
!--    Calculate horizontally averaged profiles. Each quantity is identified by the index
!--    "det_pr_index". These profile numbers must also be assigned to the respective strings
!--    given by data_output_pr in routine det_check_data_output_pr.
       DO  ii = 1, det_pr_count

          ind = pr_palm + max_pr_cs + ii

          SELECT CASE( det_pr_index(ii) )

             CASE( 1 )
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. clay_calculated )  CALL det_actions( 'calculate clay' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + clay(k,j,i)  * rmask(j,i,sr) *      &
                                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 2 )
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. dust_calculated )  CALL det_actions( 'calculate dust' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + dust(k,j,i)  * rmask(j,i,sr) *      &
                                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 3 )
!
!--             Calculate diagnostic quantities, if not done so far.
                IF ( .NOT. silt_calculated )  CALL det_actions( 'calculate silt' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + silt(k,j,i)  * rmask(j,i,sr) *      &
                                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

          END SELECT
       ENDDO

    ELSEIF ( mode == 'time_series' )  THEN

       IF ( .NOT. ALLOCATED( ts_value_l ) )  THEN
          ALLOCATE ( ts_value_l(dots_num_det) )
          ts_value_l = 0.0_wp
       ENDIF

!
!--    Calculate mean and maximum dust emission flux and dust deposition flux (only from bulk
!--    values of all dust bins).
!--    First calculate the bulk fluxes for each horizontal grid point.
       ALLOCATE( bulk_depo(nys:nyn,nxl:nxr) )
       ALLOCATE( bulk_emis(nys:nyn,nxl:nxr) )
       bulk_depo = 0.0_wp
       bulk_emis = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                IF ( surf_def%upward(m) )  THEN
                   bulk_depo(j,i) = bulk_depo(j,i) +                                               &
                                    SUM( surf_def%dm_depo_flux(:,m) ) * rmask(j,i,sr)
                   bulk_emis(j,i) = bulk_emis(j,i) +                                               &
                                    SUM( surf_def%dm_emis_flux(:,m) ) * rmask(j,i,sr)
                ENDIF
             ENDDO
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m) )  THEN
                   bulk_depo(j,i) = bulk_depo(j,i) +                                               &
                                    SUM( surf_lsm%dm_depo_flux(:,m) ) * rmask(j,i,sr)
                   bulk_emis(j,i) = bulk_emis(j,i) +                                               &
                                    SUM( surf_lsm%dm_emis_flux(:,m) ) * rmask(j,i,sr)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Store local sums and maxima. Deposition flux is negative, therefore MINVAL is used.
       ts_value_l(1) = SUM( bulk_depo )
       ts_value_l(2) = SUM( bulk_emis )
       ts_value_l(3) = MINVAL( bulk_depo )
       ts_value_l(4) = MAXVAL( bulk_emis )
!
!--    Collect / send values to PE0, because only PE0 outputs the time series.
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( ts_value_l(1), ts_value(dots_start_index_det,sr), 2, MPI_REAL,          &
                           MPI_SUM, comm2d, ierr )

       CALL MPI_ALLREDUCE( ts_value_l(3), ts_value(dots_start_index_det+2,sr), 2, MPI_REAL,        &
                           MPI_MAX, comm2d, ierr )
#else
       ts_value(dots_start_index_det:dots_start_index_det+dots_num_det,sr) = ts_value_l
#endif
!
!--    Normalize the sums.
       ts_value(dots_start_index_det:dots_start_index_det+1,sr) =                                  &
                            ts_value(dots_start_index_det:dots_start_index_det+1,sr) / ngp_2dh(sr)

       DEALLOCATE( bulk_depo, bulk_emis )

    ENDIF

 END SUBROUTINE det_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels for a user-defined prognostic quantity.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_swap_timelevel( mod_count )

    INTEGER(iwp) ::  id  !< index for dust size bin

    INTEGER, INTENT(IN) ::  mod_count  !< flag defining where pointers point to


    IF ( time_since_reference_point < det_start_time )  RETURN

    SELECT CASE ( mod_count )

       CASE ( 0 )

          DO  id = 1, n_dust_bins
             dm(id)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => dm_1(:,:,:,id)
             dm(id)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => dm_2(:,:,:,id)
          ENDDO

       CASE ( 1 )

          DO  id = 1, n_dust_bins
             dm(id)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => dm_2(:,:,:,id)
             dm(id)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => dm_1(:,:,:,id)
          ENDDO

    END SELECT

 END SUBROUTINE det_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the module-specific global restart data.
!> The number of dust bins is the only parameter which is not allowed to be changed during a
!> restart.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_wrd_global

    IMPLICIT NONE


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'n_dust_bins' )
       WRITE ( 14 )  n_dust_bins

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       CALL wrd_mpi_io( 'n_dust_bins', n_dust_bins )

    ENDIF

 END SUBROUTINE det_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the module-specific local (domain dependent) restart data arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE det_wrd_local

    CHARACTER(LEN=5) ::  dust_bin_name  !< name of dust mass bin on restart file

    INTEGER(iwp) ::  id  !< loop variable for dust size bins


    DO  id = 1, n_dust_bins

       WRITE( dust_bin_name, '(A3,I2.2)' )  'dm_', id

       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          CALL wrd_write_string( dust_bin_name )
          WRITE ( 14 )  dm(id)%conc
          IF ( ALLOCATED( dm(id)%conc_av ) )  THEN
             CALL wrd_write_string( dust_bin_name // '_av' )
             WRITE ( 14 )  dm(id)%conc_av
          ENDIF

       ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

          CALL wrd_mpi_io( dust_bin_name, dm(id)%conc )
          IF ( ALLOCATED( dm(id)%conc_av ) )  THEN
             CALL wrd_mpi_io( dust_bin_name // '_av', dm(id)%conc_av )
          ENDIF

       ENDIF
    ENDDO

 END SUBROUTINE det_wrd_local

 END MODULE dust_emission_and_transport_mod
