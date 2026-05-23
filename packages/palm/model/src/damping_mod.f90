!> @file damping_mod.f90
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
! @author Oliver Maas (ENERCON)
!
!
! Description:
! ------------
!> Damping module. It currently only contains inertial oscillation damping, but is intended to 
!> contain other damping features like Rayleigh-damping in the future.
!--------------------------------------------------------------------------------------------------!
 MODULE damping_mod

    USE arrays_3d,                                                                                 &
        ONLY:  u,                                                                                  &
               v

    USE calc_mean_profile_mod

    USE control_parameters,                                                                        &
        ONLY:  bc_lr,                                                                              &
               bc_ns,                                                                              &
               current_timestep_number,                                                            &
               dmp_enabled,                                                                        &
               debug_output,                                                                       &
               initializing_actions,                                                               &
               latitude,                                                                           &
               message_string,                                                                     &
               restart_data_format_output,                                                         &
               topography

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


    REAL(wp) ::  inertial_oscillation_damping_factor = 0.0_wp  !< namelist parameter

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hom_u       !< horizontal mean of u of current timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hom_u_m     !< horizontal mean of u of previous timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hom_v       !< horizontal mean of v of current timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hom_v_m     !< horizontal mean of v of previous timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tend_u_iod  !< u-tendency of inertial oscillation damping
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tend_v_iod  !< v-tendency of inertial oscillation damping
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tend_dmp    !< tendency field for damping (prognostic equations)


    SAVE

    PRIVATE

!
!-- Public variables
    PUBLIC                                                                                         &
           dmp_enabled,                                                                            &
           inertial_oscillation_damping_factor,                                                    &
           tend_dmp
!
!-- Public subroutines.
    PUBLIC dmp_actions,                                                                            &
           dmp_check_parameters,                                                                   &
           dmp_header,                                                                             &
           dmp_init_arrays,                                                                        &
           dmp_parin,                                                                              &
           dmp_rrd_global,                                                                         &
           dmp_wrd_global

    INTERFACE dmp_actions
       MODULE PROCEDURE dmp_actions
       MODULE PROCEDURE dmp_actions_ij
    END INTERFACE dmp_actions

    INTERFACE dmp_check_parameters
       MODULE PROCEDURE dmp_check_parameters
    END INTERFACE dmp_check_parameters

    INTERFACE dmp_header
       MODULE PROCEDURE dmp_header
    END INTERFACE dmp_header

    INTERFACE dmp_init_arrays
       MODULE PROCEDURE dmp_init_arrays
    END INTERFACE dmp_init_arrays

    INTERFACE dmp_parin
       MODULE PROCEDURE dmp_parin
    END INTERFACE dmp_parin

    INTERFACE dmp_rrd_global
       MODULE PROCEDURE dmp_rrd_global_ftn
       MODULE PROCEDURE dmp_rrd_global_mpi
    END INTERFACE dmp_rrd_global

    INTERFACE dmp_wrd_global
       MODULE PROCEDURE dmp_wrd_global
    END INTERFACE dmp_wrd_global

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate tendencies for inertial oscillation damping (vector-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_actions( location )

    CHARACTER(LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp) ::  k       !< grid index in z-direction


    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( inertial_oscillation_damping_factor /= 0.0_wp ) THEN

             CALL calc_mean_profile( u, 1 )
             CALL calc_mean_profile( v, 2 )

             hom_u(nzb+1:nzt) = hom(nzb+1:nzt,1,1,0)
             hom_v(nzb+1:nzt) = hom(nzb+1:nzt,1,2,0)

!        
!--          Calculate inertial oscillation damping tendencies
             IF ( current_timestep_number > 2 )  THEN
                DO  k = nzb+1, nzt
                   tend_u_iod(k) =   SIGN( inertial_oscillation_damping_factor , latitude ) *      &
                                     ( hom_v(k) - hom_v_m(k) )
                   tend_v_iod(k) = - SIGN( inertial_oscillation_damping_factor , latitude ) *      &
                                     ( hom_u(k) - hom_u_m(k) )
                ENDDO
             ENDIF
          ENDIF

!
!--       Save hom profiles of u and v of this timestep to be able to compute temporal derivatives
!--       at next time step.
          hom_u_m(nzb+1:nzt) = hom_u(nzb+1:nzt)
          hom_v_m(nzb+1:nzt) = hom_v(nzb+1:nzt)

       CASE ( 'u-tendency' )
          tend_dmp = tend_u_iod

       CASE ( 'v-tendency' )
          tend_dmp = tend_v_iod

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE dmp_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate tendencies for inertial oscillation damping (cache-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_actions_ij( i, j, location )


    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp) ::  i       !< grid index in x-direction
    INTEGER(iwp) ::  j       !< grid index in y-direction


!
!-- Next line is to avoid compiler warning about unused variables.
    IF ( i == 0  .OR.  j == 0 )  CONTINUE

    SELECT CASE ( location )

       CASE ( 'u-tendency' )
             tend_dmp = tend_u_iod

       CASE ( 'v-tendency' )
             tend_dmp = tend_v_iod

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE dmp_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for the damping module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_check_parameters

!
!-- In case of inertial oscillation damping, check whether bc_lr and bc_ns are cyclic and 
!-- whether topography is flat.
    IF ( inertial_oscillation_damping_factor > 0.0_wp )  THEN
       IF ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic' )  THEN
          message_string = 'illegal setting of bc_lr / bc_ns for ' //                              &
                           'inertial_oscillation_damping_factor'
          CALL message( 'check_parameters', 'DMP0001', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( topography /= 'flat' )  THEN
          message_string = 'illegal setting of topography for inertial_oscillation_damping_factor'
          CALL message( 'check_parameters', 'DMP0002', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE dmp_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for dmp parameters
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_header( io )

    INTEGER(iwp), INTENT(IN) ::  io   !< Unit of the output file

!
!-- Write dmp header
    WRITE( io, 1 )
    WRITE( io, 2 ) inertial_oscillation_damping_factor

!
!-- Format specifications
1   FORMAT (//' Damping settings:'/' ------------'/)
2   FORMAT ('    inertial_oscillation_damping_factor = ', F4.2)

 END SUBROUTINE dmp_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate and initialize arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_init_arrays

    IF ( debug_output )  CALL debug_message( 'dmp_init_arrays', 'start' )

!
!-- Allocate arrays for inertial oscillation damping.
    IF ( .NOT. ALLOCATED( hom_u      ) )  ALLOCATE( hom_u(nzb+1:nzt) )
    IF ( .NOT. ALLOCATED( hom_u_m    ) )  ALLOCATE( hom_u_m(nzb+1:nzt) )
    IF ( .NOT. ALLOCATED( hom_v      ) )  ALLOCATE( hom_v(nzb+1:nzt) )
    IF ( .NOT. ALLOCATED( hom_v_m    ) )  ALLOCATE( hom_v_m(nzb+1:nzt) )
    IF ( .NOT. ALLOCATED( tend_u_iod ) )  ALLOCATE( tend_u_iod(nzb+1:nzt) )
    IF ( .NOT. ALLOCATED( tend_v_iod ) )  ALLOCATE( tend_v_iod(nzb+1:nzt) )

    ALLOCATE( tend_dmp(nzb+1:nzt) )

    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
       tend_u_iod(:) = 0.0_wp
       tend_v_iod(:) = 0.0_wp
       hom_u(:)      = 0.0_wp
       hom_u_m(:)    = 0.0_wp
       hom_v(:)      = 0.0_wp
       hom_v_m(:)    = 0.0_wp
    ENDIF

    IF ( debug_output )  CALL debug_message( 'dmp_init_arrays', 'end' )

 END SUBROUTINE dmp_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read namelist &dmp_parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file


    NAMELIST /damping_parameters/  inertial_oscillation_damping_factor

!
!-- Move to the beginning of the namelist file and try to find and read the user-defined namelist
!-- damping_parameters.
    REWIND( 11 )
    READ( 11, damping_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    damping_parameters namelist was found and read correctly. Set flag that indicates that
!--    the damping module (DMP) is switched on.
       IF ( .NOT. switch_off_module )  dmp_enabled = .TRUE.


    ELSEIF ( io_status > 0 )  THEN
!
!--    damping_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'damping_parameters', line )

    ENDIF

 END SUBROUTINE dmp_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_rrd_global_ftn( found )

    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string

    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found  !< switch to indicate if variable has been found on restart file


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'hom_u_m' )
          IF ( .NOT. ALLOCATED( hom_u_m ) )  ALLOCATE( hom_u_m(1:nz) )
          READ ( 13 )  hom_u_m

       CASE ( 'hom_v_m' )
          IF ( .NOT. ALLOCATED( hom_v_m ) )  ALLOCATE( hom_v_m(1:nz) )
          READ ( 13 )  hom_v_m

       CASE ( 'inertial_oscillation_damping_factor' )
          READ ( 13 )  inertial_oscillation_damping_factor

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE dmp_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_rrd_global_mpi

    IF ( .NOT. ALLOCATED( hom_u_m ) )  ALLOCATE( hom_u_m(1:nz) )
    IF ( .NOT. ALLOCATED( hom_v_m ) )  ALLOCATE( hom_v_m(1:nz) )
    CALL rrd_mpi_io_global_array( 'hom_u_m', hom_u_m )
    CALL rrd_mpi_io_global_array( 'hom_v_m', hom_v_m )
    CALL rrd_mpi_io( 'inertial_oscillation_df', inertial_oscillation_damping_factor )

 END SUBROUTINE dmp_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the module-specific global restart data.
!> The number of dust bins is the only parameter which is not allowed to be changed during a
!> restart.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dmp_wrd_global

    IMPLICIT NONE

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'hom_u_m' )
       WRITE ( 14 )  hom_u_m

       CALL wrd_write_string( 'hom_v_m' )
       WRITE ( 14 )  hom_v_m

       CALL wrd_write_string( 'inertial_oscillation_damping_factor' )
       WRITE ( 14 )  inertial_oscillation_damping_factor

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       CALL wrd_mpi_io_global_array( 'hom_u_m', hom_u_m )
       CALL wrd_mpi_io_global_array( 'hom_v_m', hom_v_m )
       CALL wrd_mpi_io( 'inertial_oscillation_df', inertial_oscillation_damping_factor )

    ENDIF

 END SUBROUTINE dmp_wrd_global


 END MODULE damping_mod
