!> @file time_integration_spinup.f90
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
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Integration in time of the non-atmospheric model components such as land surface model and urban
!> surface model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE time_integration_spinup

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  pt,                                                                                 &
               pt_p,                                                                               &
               u,                                                                                  &
               u_init,                                                                             &
               v,                                                                                  &
               v_init

    USE control_parameters,                                                                        &
        ONLY:  averaging_interval_pr,                                                              &
               constant_diffusion,                                                                 &
               constant_flux_layer,                                                                &
               coupling_start_time,                                                                &
               data_output_during_spinup,                                                          &
               dcep,                                                                               &
               debug_output_timestep,                                                              &
               debug_string,                                                                       &
               dopr_n,                                                                             &
               do_sum,                                                                             &
               dt_averaging_input_pr,                                                              &
               dt_dopr,                                                                            &
               dt_dots,                                                                            &
               dt_do2d_xy,                                                                         &
               dt_do3d,                                                                            &
               dt_run_control_spinup,                                                              &
               dt_3d,                                                                              &
               humidity,                                                                           &
               indoor_model,                                                                       &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               land_surface,                                                                       &
               simulated_time,                                                                     &
               simulated_time_chr,                                                                 &
               skip_time_dopr,                                                                     &
               skip_time_do2d_xy,                                                                  &
               skip_time_do3d,                                                                     &
               slurb,                                                                              &
               spinup_pt_amplitude,                                                                &
               spinup_phase,                                                                       &
               spinup_pt_mean,                                                                     &
               spinup_time,                                                                        &
               surface_output,                                                                     &
               timestep_count,                                                                     &
               time_dopr,                                                                          &
               time_dopr_av,                                                                       &
               time_dots,                                                                          &
               time_do2d_xy,                                                                       &
               time_do3d,                                                                          &
               time_since_reference_point,                                                         &
               urban_surface
#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  nested_run
#endif


    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE dcep_mod,                                                                                  &
        ONLY:  dcep_main

    USE diagnostic_output_quantities_mod,                                                          &
        ONLY:  doq_calculate

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               nysg,                                                                               &
               nyng,                                                                               &
               nxlg,                                                                               &
               nxrg

    USE indoor_model_mod,                                                                          &
        ONLY:  dt_indoor,                                                                          &
               im_main_heatcool,                                                                   &
               indoor_during_spinup,                                                               &
               time_indoor

    USE kinds

    USE land_surface_model_mod,                                                                    &
        ONLY:  dt_lsm,                                                                             &
               lsm_energy_balance,                                                                 &
               lsm_swap_timelevel

    USE module_interface,                                                                          &
        ONLY:  module_interface_actions,                                                           &
               module_interface_timestep

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time,                                                                      &
               seconds_per_hour

    USE pegrid

    USE radiation_model_mod,                                                                       &
        ONLY:  dt_radiation,                                                                       &
               force_radiation_call,                                                               &
               radiation,                                                                          &
               radiation_called,                                                                   &
               radiation_control,                                                                  &
               radiation_interaction,                                                              &
               radiation_interactions,                                                             &
               time_radiation

    USE slurb_mod,                                                                                 &
        ONLY :  slurb_model,                                                                       &
                slurb_swap_timelevel

    USE statistics,                                                                                &
        ONLY:  flow_statistics_called

    USE surface_data_output_mod,                                                                   &
        ONLY:  dt_dosurf,                                                                          &
               skip_time_dosurf,                                                                   &
               surface_data_output,                                                                &
               time_dosurf

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  surface_layer_fluxes

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    USE urban_surface_mod,                                                                         &
        ONLY:  dt_usm,                                                                             &
               usm_energy_balance,                                                                 &
               usm_swap_timelevel

    IMPLICIT NONE

    CHARACTER(LEN=10) ::  time_since_reference_point_chr  !< time since reference point, i.e., negative during spinup
    CHARACTER(LEN=10) ::  time_to_string                  !< function to transform time from real to character string

    INTEGER(iwp) ::  current_timestep_number_spinup = 0  !< number if timestep during spinup
    INTEGER(iwp) ::  day_of_year                         !< day of the year

    INTEGER(iwp) ::  i  !< running index
    INTEGER(iwp) ::  j  !< running index
    INTEGER(iwp) ::  k  !< running index
    INTEGER(iwp) ::  m  !< running index

    LOGICAL ::  run_control_header_spinup = .FALSE.  !< flag parameter for steering whether the header information must be output

    REAL(wp) ::  dt_save                           !< temporary storage for time step
    REAL(wp) ::  pt_spinup                         !< temporary storage of temperature
    REAL(wp) ::  second_of_day                     !< second of the day
    REAL(wp) ::  time_run_control_spinup = 0.0_wp  !< time since last spinup run control output

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_save  !< temporary storage of temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_save   !< temporary storage of u wind component
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_save   !< temporary storage of v wind component


!
!-- Indicate other modules that we are in spinup phase.
    spinup_phase = .TRUE.

!
!-- Save 3D arrays because they are to be changed for spinup purpose
    ALLOCATE( pt_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( u_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( v_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    CALL exchange_horiz( pt, nbgp )
    CALL exchange_horiz( u, nbgp )
    CALL exchange_horiz( v, nbgp )

    pt_save = pt
    u_save  = u
    v_save  = v

!
!-- Set the same wall-adjacent velocity to all grid points. The sign of the original velocity field
!-- must be preserved because the surface schemes crash otherwise. The precise reason is still
!-- unknown. A minimum velocity of 0.1 m/s is used to maintain turbulent transfer at the surface.
    IF ( land_surface )  THEN
       DO  m = 1, surf_lsm%ns
          i = surf_lsm%i(m)
          j = surf_lsm%j(m)
          k = surf_lsm%k(m)
          u(k,j,i) = SIGN( 1.0_wp, u_init(k) ) * MAX( ABS( u_init(k) ), 0.1_wp )
          v(k,j,i) = SIGN( 1.0_wp, v_init(k) ) * MAX( ABS( v_init(k) ), 0.1_wp )
       ENDDO
    ENDIF

    IF ( urban_surface )  THEN
       DO  m = 1, surf_usm%ns
          i = surf_usm%i(m)
          j = surf_usm%j(m)
          k = surf_usm%k(m)
          u(k,j,i) = SIGN( 1.0_wp, u_init(k) ) * MAX( ABS( u_init(k) ), 0.1_wp )
          v(k,j,i) = SIGN( 1.0_wp, v_init(k) ) * MAX( ABS( v_init(k) ), 0.1_wp )
       ENDDO
    ENDIF

    CALL exchange_horiz( u, nbgp )
    CALL exchange_horiz( v, nbgp )

    dt_save = dt_3d

    CALL location_message( 'surface model spinup time-stepping', 'start' )

!
!-- Start of the time loop
    DO  WHILE ( simulated_time < spinup_time )

       CALL cpu_log( log_point_s(15), 'timesteps spinup', 'start' )

       IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'time_integration_spinup', simulated_time
           CALL debug_message( debug_string, 'start' )
       ENDIF

!
!--    Calculate timestep based on diffsuion criterion for lsm and usm.
       CALL module_interface_timestep
       dt_3d = MIN( dt_lsm, dt_usm )

!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend on the timestep scheme
          CALL timestep_scheme_steering


!
!--       Estimate a near-surface air temperature based on the position of the sun and user input
!--       about mean temperature and amplitude. The time is shifted by one hour to simulate a lag
!--       between air temperature and incoming radiation.
          CALL get_date_time( simulated_time - spinup_time - seconds_per_hour,                     &
                              day_of_year = day_of_year, second_of_day = second_of_day )

          pt_spinup = spinup_pt_mean + spinup_pt_amplitude *                                       &
                      solar_angle( day_of_year, second_of_day )

!
!--       Map air temperature to all grid points in the vicinity of a surface element
          IF ( land_surface )  THEN
             DO  m = 1, surf_lsm%ns
                i = surf_lsm%i(m)
                j = surf_lsm%j(m)
                k = surf_lsm%k(m)
                pt(k,j,i) = pt_spinup
             ENDDO
          ENDIF

          IF ( urban_surface )  THEN
             DO  m = 1, surf_usm%ns
                i = surf_usm%i(m)
                j = surf_usm%j(m)
                k = surf_usm%k(m)
                pt(k,j,i) = pt_spinup
                !!!!!!!!!!!!!!!!HACK!!!!!!!!!!!!!
                surf_usm%pt1(m) = pt_spinup
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ENDDO
          ENDIF

          IF ( .NOT. land_surface  .AND.  slurb )  THEN
             DO  m = 1, surf_def%ns
                i = surf_def%i(m)
                j = surf_def%j(m)
                k = surf_def%k(m)
                pt(k,j,i) = pt_spinup
             ENDDO
          ENDIF

          CALL exchange_horiz( pt, nbgp )

!
!--       Swap the time levels in preparation for the next time step.
          timestep_count = timestep_count + 1

          IF ( land_surface )  THEN
              CALL lsm_swap_timelevel( 0 )
             CALL lsm_swap_timelevel( MOD( timestep_count, 2 ) )
          ENDIF

          IF ( slurb )  THEN
             CALL slurb_swap_timelevel( 0 )
             CALL slurb_swap_timelevel( MOD( timestep_count, 2 ) )
          ENDIF

          IF ( urban_surface )  THEN
             CALL usm_swap_timelevel( 0 )
             CALL usm_swap_timelevel( MOD( timestep_count, 2 ) )
          ENDIF

!
!--       If required, compute virtual potential temperature
          IF ( humidity )  THEN
             CALL compute_vpt
          ENDIF

!
!--       Compute the diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN
!
!--          First the vertical (and horizontal) fluxes in the surface (constant flux) layer are
!--          computed
             IF ( constant_flux_layer )  THEN
                CALL surface_layer_fluxes
             ENDIF
!
!--          If required, solve the energy balance for the surface and run soil model. Call for
!--          horizontal as well as vertical surfaces. The prognostic equation for soil moisure is
!--          switched off
             IF ( land_surface )  THEN
                CALL lsm_energy_balance
             ENDIF
!
!--          If required, solve the energy balance for SLUrb surfaces.
             IF ( slurb )  THEN
                CALL slurb_model
             ENDIF
!
!--          If required, solve the energy balance for urban surfaces and run the material heat
!--          model.
             IF ( urban_surface )  THEN
                CALL usm_energy_balance
             ENDIF

          ENDIF

       ENDDO   ! Intermediate step loop

!
!--    If required, calculate radiative fluxes and heating rates
       IF ( radiation )  THEN

          time_radiation = time_radiation + dt_3d

          IF ( time_radiation >= dt_radiation .OR. force_radiation_call )  THEN

             IF ( .NOT. force_radiation_call )  THEN
                time_radiation = time_radiation - dt_radiation
             ENDIF

             CALL radiation_control

             IF ( radiation_interactions )  THEN
                CALL radiation_interaction
             ENDIF
!
!--          Reset forcing of radiation call
             force_radiation_call = .FALSE.
!
!--          Set a flag that radiation was called this timestep.
             radiation_called = .TRUE.

          ELSE
!
!--          Set a flag that radiation wasn't called this timestep.
             radiation_called = .FALSE.
          ENDIF
       ENDIF
!
!--    If DCEP set, call dcep main routine.
       IF ( dcep )  CALL dcep_main

!
!--    If required, calculate indoor temperature, waste heat, heat flux
!--    through wall, etc.
!--    dt_indoor steers the frequency of the indoor model calculations.
!--    Note, at first timestep indoor model is called, in order to provide
!--    a waste heat flux.
       IF ( indoor_model  .AND.  indoor_during_spinup )  THEN

          time_indoor = time_indoor + dt_3d

          IF ( time_indoor >= dt_indoor  .OR.  current_timestep_number_spinup == 0 )  THEN
             IF ( time_indoor >= dt_indoor )  time_indoor = time_indoor - dt_indoor
             CALL im_main_heatcool
          ENDIF

       ENDIF
!
!--    Increase simulation time and output times.
       current_timestep_number_spinup = current_timestep_number_spinup + 1
       simulated_time                 = simulated_time   + dt_3d
       simulated_time_chr             = time_to_string( simulated_time )
       time_run_control_spinup        = time_run_control_spinup + dt_3d
       time_since_reference_point     = simulated_time - coupling_start_time
       time_since_reference_point_chr = time_to_string( time_since_reference_point )
!
!--    After last timestep, set the simulated time so that it exactly matches the given spinup time.
!--    Otherwise, the simulation after spinup will not start at exactly t = 0. Small differences
!--    between the simulated time during spinup and the given spinup time itself are ignored here,
!--    because they are less equal the timestep.
!--    The time after spinup is based on time_since_reference_point.
       IF ( simulated_time >= spinup_time )  THEN
         simulated_time                 = spinup_time
         time_since_reference_point     = simulated_time - coupling_start_time
         time_since_reference_point_chr = time_to_string( time_since_reference_point )
       ENDIF

       IF ( data_output_during_spinup )  THEN

          IF ( simulated_time >= skip_time_do2d_xy )  THEN
             time_do2d_xy      = time_do2d_xy     + dt_3d
          ENDIF
          IF ( simulated_time >= skip_time_do3d    )  THEN
             time_do3d         = time_do3d        + dt_3d
          ENDIF
          time_dots            = time_dots        + dt_3d
          IF ( simulated_time >= skip_time_dopr )  THEN
             time_dopr         = time_dopr        + dt_3d
          ENDIF
          IF ( surface_output )  THEN
             IF ( simulated_time >= skip_time_dosurf )  THEN
                time_dosurf    = time_dosurf + dt_3d
             ENDIF
          ENDIF
!
!--       Carry out statistical analysis and output at the requested output times.
!--       The MOD function is used for calculating the output time counters (like time_dopr) in
!--       order to regard a possible decrease of the output time interval in case of restart runs.
!--       Set a flag indicating that so far no statistics have been created for this time step
          flow_statistics_called = .FALSE.
!
!--       Calculate module specific output.
          CALL module_interface_actions( 'do_integration_spinup' )
!
!--       If required, call flow_statistics for averaging in time
          IF ( averaging_interval_pr /= 0.0_wp  .AND.                                              &
               ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.                             &
               simulated_time >= skip_time_dopr )                                                  &
          THEN
             time_dopr_av = time_dopr_av + dt_3d
             IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
                do_sum = .TRUE.
                time_dopr_av = MOD( time_dopr_av, MAX( dt_averaging_input_pr, dt_3d ) )
             ENDIF
          ENDIF
          IF ( do_sum )  CALL flow_statistics
!
!--       Output of profiles
          IF ( time_dopr >= dt_dopr )  THEN
             IF ( dopr_n /= 0 )  CALL data_output_profiles
             time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
             time_dopr_av = 0.0_wp    ! Due to averaging (see above)
          ENDIF

!--       Output of time series
          IF ( time_dots >= dt_dots )  THEN
             CALL data_output_tseries
             time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
          ENDIF
!
!--       2d-data output (cross-sections)
          IF ( time_do2d_xy >= dt_do2d_xy )  THEN
             CALL doq_calculate
             CALL data_output_2d( 'xy', 0 )
             time_do2d_xy = MOD( time_do2d_xy, MAX( dt_do2d_xy, dt_3d ) )
          ENDIF
!
!--       3d-data output (volume data)
          IF ( time_do3d >= dt_do3d )  THEN
             CALL doq_calculate
             CALL data_output_3d( 0 )
             time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
          ENDIF
!
!--       Output of surface data
          IF ( surface_output )  THEN
             IF ( time_dosurf >= dt_dosurf )  THEN
                CALL surface_data_output( 0 )
                time_dosurf = MOD( time_dosurf, MAX( dt_dosurf, dt_3d ) )
             ENDIF
          ENDIF

       ENDIF

       CALL cpu_log( log_point_s(15), 'timesteps spinup', 'stop' )
!
!--    Run control output
       IF ( myid == 0  .AND. ( time_run_control_spinup >= dt_run_control_spinup  .OR.              &
                               current_timestep_number_spinup == 1  .OR.                           &
                               simulated_time >= spinup_time ) )                                   &
       THEN
!
!--       If necessary, write header
          IF ( .NOT. run_control_header_spinup )  THEN
             CALL check_open( 15 )
             WRITE ( 15, 100 )
             run_control_header_spinup = .TRUE.
          ENDIF
!
!--       Write some general information about the spinup in run control file
          WRITE ( 15, 101 )  current_timestep_number_spinup, time_since_reference_point_chr,       &
                             dt_3d, pt_spinup
!
!--       Write buffer contents to disc immediately
          FLUSH( 15 )

          IF ( time_run_control_spinup >= dt_run_control_spinup )  THEN
             time_run_control_spinup = MOD( time_run_control_spinup,                               &
                                            MAX( dt_run_control_spinup, dt_3d ) )
          ENDIF
       ENDIF

       IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'time_integration_spinup', simulated_time, pt_spinup
           CALL debug_message( debug_string, 'end' )
       ENDIF


    ENDDO   ! Time loop

!
!-- Write back saved arrays to the 3D arrays
    pt   = pt_save
    pt_p = pt_save
    u    = u_save
    v    = v_save

!
!-- Reset time step
    dt_3d = dt_save

!-- Force radiation step in time zero
!-- It is performed at the end of init_radiation in case of run without spinup
    time_radiation = dt_radiation

    DEALLOCATE( pt_save )
    DEALLOCATE( u_save )
    DEALLOCATE( v_save )

#if defined( __parallel )
!
!>  TODO: This barrier shouldn't be required. Test removing it.
    IF ( nested_run )  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif

    spinup_phase = .FALSE.

    CALL location_message( 'surface model spinup time-stepping', 'finished' )

!
!-- Formats
100 FORMAT (///'Spinup control output:'/                                                           &
               '----------------------'//                                                          &
            'ITER.   HHH:MM:SS    DT  PT(z_MO)'/                                                   &
            '---------------------------------')
101 FORMAT (I5,2X,A10,1X,F6.2,3X,F6.2,2X,F6.2)

 CONTAINS

!
!-- Returns the cosine of the solar zenith angle at a given time. This routine is similar to that
!-- for calculation zenith (see radiation_model_mod.f90)
    !> @todo Load function calc_zenith of radiation model instead of rewrite the function here.
    FUNCTION solar_angle( day_of_year, second_of_day )

       USE basic_constants_and_equations_mod,                                                      &
           ONLY:  pi

       USE kinds

       USE radiation_model_mod,                                                                    &
           ONLY:  decl_1,                                                                          &
                  decl_2,                                                                          &
                  decl_3,                                                                          &
                  lat,                                                                             &
                  lon

       IMPLICIT NONE


       INTEGER(iwp), INTENT(IN) ::  day_of_year  !< day of the year

       REAL(wp)             ::  declination    !< solar declination angle
       REAL(wp)             ::  hour_angle     !< solar hour angle
       REAL(wp), INTENT(IN) ::  second_of_day  !< current time of the day in UTC
       REAL(wp)             ::  solar_angle    !< cosine of the solar zenith angle
!
!--    Calculate solar declination and hour angle
       declination = ASIN( decl_1 * SIN( decl_2 * REAL( day_of_year, KIND = wp) - decl_3 ) )
       hour_angle  = 2.0_wp * pi * ( second_of_day / 86400.0_wp ) + lon - pi

!
!--    Calculate cosine of solar zenith angle
       solar_angle = SIN( lat ) * SIN( declination ) + COS( lat ) * COS( declination ) *           &
                     COS( hour_angle )

    END FUNCTION solar_angle


 END SUBROUTINE time_integration_spinup
