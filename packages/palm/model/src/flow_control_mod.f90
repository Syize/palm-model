!> @file flow_control_mod.f90
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
! Copyright 2025 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Oliver Maas (ENERCON)
!
!
! Description:
! ------------
!> Flow control module. It currently only contains the geostrophic wind controller, but is intended 
!> to contain other flow control features like large-scale forcing, nudging and volume flow 
!> conservation.
!--------------------------------------------------------------------------------------------------!
 MODULE flow_control_mod

    USE arrays_3d,                                                                                 &
        ONLY:  u,                                                                                  &
               v,                                                                                  &
               ug,                                                                                 &
               vg,                                                                                 &
               zu

    USE calc_mean_profile_mod

    USE control_parameters,                                                                        &
        ONLY:  averaging_interval_pr,                                                              &
               bc_lr,                                                                              &
               bc_ns,                                                                              &
               current_timestep_number,                                                            &
               dt_dopr,                                                                            &
               end_time,                                                                           &
               f,                                                                                  &
               fct_enabled,                                                                        &
               debug_output,                                                                       &
               dt_3d,                                                                              &
               initializing_actions,                                                               &
               latitude,                                                                           &
               message_string,                                                                     &
               rayleigh_damping_factor,                                                            &
               restart_data_format_output,                                                         &
               simulated_time,                                                                     &
               time_dopr,                                                                          &
               time_since_reference_point,                                                         &
               topography

    USE damping_mod,                                                                               &
        ONLY:  inertial_oscillation_damping_factor

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rrd_mpi_io,                                                                         &
               wrd_mpi_io,                                                                         &
               rrd_mpi_io_global_array,                                                            &
               wrd_mpi_io_global_array

    USE statistics,                                                                                &
        ONLY: hom


    IMPLICIT NONE

    LOGICAL ::  geostrophic_wind_controller = .FALSE.  !< namelist parameter
    LOGICAL ::  do_sum_ug                   = .FALSE.  !< switch to determine whether summation of ug and vg shall be done

    INTEGER(iwp) ::  average_count_ug = 0  !< counter for summation of geostrophic wind
    INTEGER(iwp) ::  k1 = 0                !< lower interpolation index
    INTEGER(iwp) ::  k2 = 0                !< upper interpolation index

    REAL(wp) ::  hom_u_interp = 0.0_wp      !< horizontal mean of u interpolated to target height
    REAL(wp) ::  hom_v_interp = 0.0_wp      !< horizontal mean of v interpolated to target height
    REAL(wp) ::  kp = 10.0_wp               !< proportional gain factor
    REAL(wp) ::  ki = 20.0_wp               !< integral gain factor
    REAL(wp) ::  ki_f = 0.0_wp              !< integral gain factor multiplied by Coriolis parameter
    REAL(wp) ::  target_height = 0.0_wp     !< height at which target_u and target_v shall be reached
    REAL(wp) ::  target_u = 0.0_wp          !< target wind speed in x-direction
    REAL(wp) ::  target_v = 0.0_wp          !< target wind speed in y-direction
    REAL(wp) ::  error_u_p = 0.0_wp         !< proportional controller error for u
    REAL(wp) ::  error_u_i = 0.0_wp         !< integral controller error for u
    REAL(wp) ::  error_v_p = 0.0_wp         !< proportional controller error for v
    REAL(wp) ::  error_v_i = 0.0_wp         !< integral controller error for v
    REAL(wp) ::  interp_weighting = 0.0_wp  !< weighting factor for interpolation
    REAL(wp) ::  skip_time_av_ug = 0.0_wp   !< time after which averaging of ug and vg starts

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ug_sum  !< sum array for temporally averaged ug
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vg_sum  !< sum array for temporally averaged vg


    SAVE

    PRIVATE

!
!-- Public subroutines.
    PUBLIC fct_actions,                                                                            &
           fct_check_parameters,                                                                   &
           fct_header,                                                                             &
           fct_init,                                                                               &
           fct_init_arrays,                                                                        &
           fct_parin,                                                                              &
           fct_rrd_global,                                                                         &
           fct_wrd_global

    INTERFACE fct_actions
       MODULE PROCEDURE fct_actions
    END INTERFACE fct_actions

    INTERFACE fct_check_parameters
       MODULE PROCEDURE fct_check_parameters
    END INTERFACE fct_check_parameters

    INTERFACE fct_header
       MODULE PROCEDURE fct_header
    END INTERFACE fct_header

    INTERFACE fct_init
       MODULE PROCEDURE fct_init
    END INTERFACE fct_init

    INTERFACE fct_init_arrays
       MODULE PROCEDURE fct_init_arrays
    END INTERFACE fct_init_arrays

    INTERFACE fct_parin
       MODULE PROCEDURE fct_parin
    END INTERFACE fct_parin

    INTERFACE fct_rrd_global
       MODULE PROCEDURE fct_rrd_global_ftn
       MODULE PROCEDURE fct_rrd_global_mpi
    END INTERFACE fct_rrd_global

    INTERFACE fct_wrd_global
       MODULE PROCEDURE fct_wrd_global
    END INTERFACE fct_wrd_global

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute flow control actions (vector-optimized).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_actions( location )

    CHARACTER(LEN=*), INTENT(IN) ::  location  !< call location string


    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( geostrophic_wind_controller )  THEN

!
!--          Calculate horizontally-averaged profiles of u and v and store the
             CALL calc_mean_profile( u, 1 )
             CALL calc_mean_profile( v, 2 )
!
!--          Interpolate to target height
             hom_u_interp = hom(k1,1,1,0) + interp_weighting * ( hom(k2,1,1,0) - hom(k1,1,1,0) )
             hom_v_interp = hom(k1,1,2,0) + interp_weighting * ( hom(k2,1,2,0) - hom(k1,1,2,0) )
!
!--          Calculate geostrophic wind controller error values for u and v.
!--          p: proportional error
!--          i: integral error
             error_u_p = target_u - hom_u_interp
             error_u_i = error_u_i + error_u_p*dt_3d
             error_v_p = target_v - hom_v_interp
             error_v_i = error_v_i + error_v_p*dt_3d
!
!--          Calculate geostrophic wind components with PI-controller.
             ug = kp * error_u_p + ki_f * error_u_i
             vg = kp * error_v_p + ki_f * error_v_i

          ENDIF

       CASE ( 'after_integration' )
!
!--       Determine whether summation of ug and vg shall be done.
          IF ( averaging_interval_pr /= 0.0_wp  .AND.                                              &
               ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.                             &
               time_since_reference_point >= skip_time_av_ug )                                     &
          THEN
             do_sum_ug = .TRUE.
          ENDIF
!
!--       Sum up ug and vg for calculation of temporal average.
          IF ( do_sum_ug  .AND.  simulated_time /= 0.0_wp )  THEN
             IF ( average_count_ug == 0 )  THEN
                ug_sum = 0.0_wp
                vg_sum = 0.0_wp
             ENDIF
             ug_sum(:) = ug_sum(:) + ug(:)
             vg_sum(:) = vg_sum(:) + vg(:)
             average_count_ug = average_count_ug + 1
             do_sum_ug = .FALSE.
          ENDIF

       CASE ( 'after_time_integration' )
!
!--       Calculate temporal average of ug and vg and write it to ug and vg so that they are 
!--       written to restart file for main run. Only do this if end_time is reached and not for
!--       other cases for which 'after_time_integration' is reached (e.g. during restarts).
          IF ( simulated_time >= end_time )  THEN
             IF ( averaging_interval_pr == 0.0_wp )  THEN
                ug_sum(:) = ug(:)
                vg_sum(:) = vg(:)
             ENDIF
             ug(:) = ug_sum(:) / REAL( average_count_ug, KIND=wp )
             vg(:) = vg_sum(:) / REAL( average_count_ug, KIND=wp )
          ENDIF

       CASE DEFAULT

          CONTINUE

    END SELECT

 END SUBROUTINE fct_actions

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for the flow_control module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_check_parameters

!
!-- In case of geostrophic wind controller, check whether bc_lr and bc_ns are cyclic,
!-- whether topography is flat, whether inertial oscillation damping is activated and whether
!-- target_height is positive and below the uppermost gridpoint.
    IF ( geostrophic_wind_controller )  THEN
       IF ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic' )  THEN
          message_string = 'illegal setting of bc_lr / bc_ns for ' //                              &
                           'geostrophic wind controller'
          CALL message( 'check_parameters', 'FCT0001', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( topography /= 'flat' )  THEN
          message_string = 'illegal setting of topography for geostrophic wind controller'
          CALL message( 'check_parameters', 'FCT0002', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inertial_oscillation_damping_factor == 0.0_wp )  THEN
          message_string = 'illegal value for inertial_oscillation_damping_factor ' //             &
                           'for geostrophic wind controller'
          CALL message( 'check_parameters', 'FCT0003', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( target_height > zu(nzt)  .OR.  target_height <=0.0_wp )  THEN
          message_string = 'illegal value for target_height ' //             &
                           'for geostrophic wind controller'
          CALL message( 'check_parameters', 'FCT0004', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( rayleigh_damping_factor /= 0.0_wp )  THEN
          message_string = 'illegal value for rayleigh_damping_factor ' //             &
                           'for geostrophic wind controller'
          CALL message( 'check_parameters', 'FCT0005', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE fct_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for fct parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_header( io )

    INTEGER(iwp), INTENT(IN) ::  io  !< unit of the output file


!
!-- Write fct header.
    WRITE ( io, 1 )
    IF ( geostrophic_wind_controller )  THEN
       WRITE ( io, 2 )
       WRITE ( io, 3 )  target_height
       WRITE ( io, 4 )  target_u
       WRITE ( io, 5 )  target_v
       WRITE ( io, 6 )  kp
       WRITE ( io, 7 )  ki
    ENDIF

!
!-- Format specifications.
1   FORMAT (//' Flow control settings:'/' ------------'/)
2   FORMAT ('    Geostrophic wind controller is activated')
3   FORMAT ('    target_height = ', F10.2)
4   FORMAT ('    target_u      = ', F10.2)
5   FORMAT ('    target_v      = ', F10.2)
6   FORMAT ('    kp            = ', F10.2)
7   FORMAT ('    ki            = ', F10.2)

 END SUBROUTINE fct_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the flow control module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_init

    INTEGER(iwp) ::  k  !< loop index in z-direction


    IF ( debug_output )  CALL debug_message( 'fct_init', 'start' )

!
!-- Initialize geostrophic wind controller.
    IF ( geostrophic_wind_controller )  THEN

!
!--    Calculate weighting factor for interpolation between two grid levels zu(k1) and zu(k2).
       DO  k = nzb, nzt
          IF ( zu(k+1) >= target_height )  THEN
             k1 = k
             k2 = k + 1
             EXIT
          ENDIF
       ENDDO

       interp_weighting = ( target_height - zu(k1) ) / ( zu(k2) - zu(k1) )

!
!--    Multiply ki with abs(f) to ensure same controller behavior for all latitudes.
       ki_f = ki * ABS( f )

       skip_time_av_ug = end_time - averaging_interval_pr

    ENDIF

    IF ( debug_output )  CALL debug_message( 'fct_init', 'end' )

 END SUBROUTINE fct_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate and initialize arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_init_arrays


    IF ( debug_output )  CALL debug_message( 'fct_init_arrays', 'start' )

!
!-- Allocate arrays for geostrophic wind controller
    IF ( .NOT. ALLOCATED( ug_sum ) )  ALLOCATE( ug_sum(0:nzt+1) )
    IF ( .NOT. ALLOCATED( vg_sum ) )  ALLOCATE( vg_sum(0:nzt+1) )

    IF ( debug_output )  CALL debug_message( 'fct_init_arrays', 'end' )

 END SUBROUTINE fct_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read namelist &fct_parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file


    NAMELIST /flow_control_parameters/                                                             &
       geostrophic_wind_controller,                                                                &
       ki,                                                                                         &
       kp,                                                                                         &
       target_height,                                                                              &
       target_u,                                                                                   &
       target_v

!
!-- Move to the beginning of the namelist file and try to find and read the user-defined namelist
!-- flow_control_parameters.
    REWIND( 11 )
    READ( 11, flow_control_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    flow_control_parameters namelist was found and read correctly. Set flag that indicates that
!--    the flow control module (FCT) is switched on.
       IF ( .NOT. switch_off_module )  fct_enabled = .TRUE.


    ELSEIF ( io_status > 0 )  THEN
!
!--    flow_control_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'flow_control_parameters', line )

    ENDIF

 END SUBROUTINE fct_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_rrd_global_ftn( found )

    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string

    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found  !< switch to indicate if variable has been found on restart file


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'average_count_ug' )
          READ ( 13 )  average_count_ug

       CASE ( 'error_u_i' )
          READ ( 13 )  error_u_i

       CASE ( 'error_v_i' )
          READ ( 13 )  error_v_i

       CASE ( 'ug_sum' )
          IF ( .NOT. ALLOCATED( ug_sum ) )  ALLOCATE( ug_sum(0:nz+1) )
          READ ( 13 )  ug_sum

       CASE ( 'vg_sum' )
          IF ( .NOT. ALLOCATED( vg_sum ) )  ALLOCATE( vg_sum(0:nz+1) )
          READ ( 13 )  vg_sum

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE fct_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_rrd_global_mpi

    IF ( .NOT. ALLOCATED( ug_sum ) )  ALLOCATE( ug_sum(0:nz+1) )
    IF ( .NOT. ALLOCATED( vg_sum ) )  ALLOCATE( vg_sum(0:nz+1) )
    CALL rrd_mpi_io( 'average_count_ug', average_count_ug )
    CALL rrd_mpi_io( 'error_u_i', error_u_i )
    CALL rrd_mpi_io( 'error_v_i', error_v_i )
    CALL rrd_mpi_io_global_array( 'ug_sum', ug_sum )
    CALL rrd_mpi_io_global_array( 'vg_sum', vg_sum )

 END SUBROUTINE fct_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the module-specific global restart data.
!> The number of dust bins is the only parameter which is not allowed to be changed during a
!> restart.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fct_wrd_global

    IMPLICIT NONE

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'average_count_ug' )
       WRITE ( 14 )  average_count_ug

       CALL wrd_write_string( 'error_u_i' )
       WRITE ( 14 )  error_u_i

       CALL wrd_write_string( 'error_v_i' )
       WRITE ( 14 )  error_v_i

       CALL wrd_write_string( 'ug_sum' )
       WRITE ( 14 )  ug_sum

       CALL wrd_write_string( 'vg_sum' )
       WRITE ( 14 )  vg_sum

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       CALL wrd_mpi_io( 'average_count_ug', average_count_ug )
       CALL wrd_mpi_io( 'error_u_i', error_u_i )
       CALL wrd_mpi_io( 'error_v_i', error_v_i )
       CALL wrd_mpi_io_global_array( 'ug_sum', ug_sum )
       CALL wrd_mpi_io_global_array( 'vg_sum', vg_sum )

    ENDIF

 END SUBROUTINE fct_wrd_global


 END MODULE flow_control_mod
