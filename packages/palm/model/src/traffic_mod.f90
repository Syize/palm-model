!> @traffic_mod.f90
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
! Copyright 2024-2025 Institute of Computer Science, Academy of Sciences, Prague, Czech Republic
!--------------------------------------------------------------------------------------------------!
!
!
! Authors:
! --------
! @author Jaroslav Resler
! @author Pavel Krc
!
! Description:
! ------------
!> The PALM module Traffic.
!> This module parametrizes the additional wind and potential temperature tendencies caused
!> by the vehicles running in the area (e.g. street canyon). It calculates the wind tendencies
!> from the relative car and air speed and from the car drag coefficient and potential tendencies
!> from the consumed petrol/gas oil and from the combustion heat. These tendencies help to avoid
!> high concentrations of the pollutants caused by the insufficient ventilation of the street
!> canyon during very stable episodes.
!> (see e.g. Resler et al.: Challenges of high-fidelity air quality modeling in urban environments
!> - PALM sensitivity study during stable conditions, GMD 2024).
!> @todo Consider adding of the imposed SGS-TKE tendency.
!--------------------------------------------------------------------------------------------------!
 MODULE traffic_mod

    USE arrays_3d

    USE basic_constants_and_equations_mod,                                                         &
       ONLY:                                                                                       &
          c_p

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  fetch_species_source_by_indices

    USE control_parameters

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices

    USE kinds

#if defined ( __netcdf )
    USE NETCDF
#endif

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               open_read_file

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time,                                                                      &
               date_time_str_len

    USE pegrid

    IMPLICIT NONE


    SAVE

!
!-- Configuration traffic variables (see module documentation)
    LOGICAL                                    ::  continuous_tendency = .FALSE.
    LOGICAL                                    ::  theta_tendency = .TRUE.
    LOGICAL                                    ::  wind_tendency = .TRUE.

    PRIVATE

!
!-- Type trm for one grid
    TYPE trm_type
!
!--    Grid coordinates
       INTEGER(iwp)                            ::  i                !< i-coordinate
       INTEGER(iwp)                            ::  j                !< j-coordinate
       INTEGER(iwp)                            ::  k                !< k-coordinate

!
!--    parameters of the street lines
       REAL(wp)                                ::  dirx             !< x-coord. of street direction vector
       REAL(wp)                                ::  diry             !< y-coord. of street direction vector
       REAL(wp)                                ::  frac             !< fraction of the traffic line
       REAL(wp)                                ::  slope            !< street slope (%)
       REAL(wp)                                ::  width            !< width of the street (m)

!
!--    mapping of the local index in mpi process to the global index
       INTEGER(iwp)                            ::  l2g              !< local index local to global array

!
!--    parameters of the traffic
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  heat             !< heat produced by car (Wats per car)
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  intensity        !< intensity of the traffic in the street
                                                                    !< (cars per hour)
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  speed            !< average speed of the traffic flow (m/s)

    END TYPE trm_type

!
!-- Type car_type - parameters of the car types
    TYPE car_type
       REAL(wp)                                ::  cd               !< average car aerodynamic resistance
       REAL(wp)                                ::  length           !< average car length (10 m)
       REAL(wp)                                ::  width            !< average car width (PC~2m,TC~3m,Bus~3m)
       REAL(wp)                                ::  height           !< average car height (PC~2m,TC~4m,Bus~3m)
    END TYPE car_type

!
!-- Ncfile
    INTEGER(iwp)                               ::  ncid             !< input nc file handle
!
!-- Dimensions
    INTEGER(iwp)                               ::  nct              !< number of car_types
    INTEGER(iwp)                               ::  nts              !< number of timestamps
!
!-- Constants for the module
    INTEGER(iwp), PARAMETER                    ::  tend_wind_u = 0  !< index of u wind tendency
    INTEGER(iwp), PARAMETER                    ::  tend_wind_v = 1  !< index of v wind tendency
    INTEGER(iwp), PARAMETER                    ::  tend_wind_w = 2  !< index of w wind tendency
    INTEGER(iwp), PARAMETER                    ::  tend_pt = 3      !< index of pt wind tendency
!
!-- Local indices and variables
    INTEGER(iwp)                               ::  trm_ns           !< number of s
    INTEGER(iwp)                               ::  trm_nsl = 0      !< number of s in local mpi subdomain
    INTEGER(iwp)                               ::  trm_ts_index     !< time step index

    CHARACTER(LEN=date_time_str_len)           ::  trm_current_ts   !< simulation time(time)
    CHARACTER(LEN=date_time_str_len)           ::  trm_next_ts      !< simulation time(time)
!
!-- Local variable arrays
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:)  ::  trm_indl         !< local index array
    CHARACTER(LEN=date_time_str_len),                                                              &
              ALLOCATABLE, DIMENSION(:)        ::  trm_ts           !< individual timestamps

    TYPE(trm_type), ALLOCATABLE, DIMENSION(:)  ::  trm              !< traffic data for one gridbox and one traffic stream

    TYPE(car_type), ALLOCATABLE, DIMENSION(:)  ::  car              !< properties of the car types

!
!-- Output variable arrays
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_theta          !< front theta tendency (instant)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_theta_av       !< front theta tendency (average)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_wind_front     !< front wind tendency (instant)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_wind_front_av  !< front wind tendency (average)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_wind_top       !< top wind tendency (instant)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_wind_top_av    !< top wind tendency (average)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_wind_side      !< side wind tendency (instant)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)    ::  tend_wind_side_av   !< side wind tendency (average)

!
!-- Public functions
    PUBLIC                                                                                         &
       trm_3d_data_averaging,                                                                      &
       trm_actions,                                                                                &
       trm_check_parameters,                                                                       &
       trm_check_data_output,                                                                      &
       trm_data_output_3d,                                                                         &
       trm_data_output_mask,                                                                       &
       trm_define_netcdf_grid,                                                                     &
       trm_init_arrays,                                                                            &
       trm_last_actions,                                                                           &
       trm_parin,                                                                                  &
       trm_rrd_local,                                                                              &
       trm_wrd_local


    INTERFACE trm_3d_data_averaging
       MODULE PROCEDURE trm_3d_data_averaging
    END INTERFACE trm_3d_data_averaging

    INTERFACE trm_actions
       MODULE PROCEDURE trm_actions
       MODULE PROCEDURE trm_actions_ij
    END INTERFACE trm_actions

    INTERFACE trm_check_parameters
       MODULE PROCEDURE trm_check_parameters
    END INTERFACE trm_check_parameters

    INTERFACE trm_check_data_output
       MODULE PROCEDURE trm_check_data_output
    END INTERFACE trm_check_data_output

     INTERFACE trm_data_output_3d
        MODULE PROCEDURE trm_data_output_3d
     END INTERFACE trm_data_output_3d

     INTERFACE trm_define_netcdf_grid
        MODULE PROCEDURE trm_define_netcdf_grid
     END INTERFACE trm_define_netcdf_grid

     INTERFACE trm_last_actions
        MODULE PROCEDURE trm_last_actions
     END INTERFACE trm_last_actions

     INTERFACE trm_parin
        MODULE PROCEDURE trm_parin
     END INTERFACE trm_parin

     INTERFACE trm_rrd_local
        MODULE PROCEDURE trm_rrd_local_ftn
        MODULE PROCEDURE trm_rrd_local_mpi
     END INTERFACE trm_rrd_local

     INTERFACE trm_wrd_local
        MODULE PROCEDURE trm_wrd_local
     END INTERFACE trm_wrd_local

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average traffic-defined output quantities as well as allocate the array
!> necessary for storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_3d_data_averaging( mode, variable )

    CHARACTER (LEN=*) ::  mode     !< Averaging mode: allocate, sum, or average
    CHARACTER (LEN=*) ::  variable !< The variable in question

    INTEGER(iwp) ::  i             !< Running index, x-direction
    INTEGER(iwp) ::  j             !< Running index, y-direction
    INTEGER(iwp) ::  k             !< Running index, z-direction

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'trm_tend_wind_front' )
             IF ( .NOT. ALLOCATED( tend_wind_front_av ) )  THEN
                ALLOCATE( tend_wind_front_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             tend_wind_front_av = 0.0_wp

          CASE ( 'trm_tend_wind_side' )
             IF ( .NOT. ALLOCATED( tend_wind_side_av ) )  THEN
                ALLOCATE( tend_wind_side_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             tend_wind_side_av = 0.0_wp

          CASE ( 'trm_tend_wind_top' )
             IF ( .NOT. ALLOCATED( tend_wind_top_av ) )  THEN
                ALLOCATE( tend_wind_top_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             tend_wind_top_av = 0.0_wp

          CASE ( 'trm_tend_theta' )
             IF ( .NOT. ALLOCATED( tend_theta_av ) )  THEN
                ALLOCATE( tend_theta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             tend_theta_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'trm_tend_wind_front' )
             IF ( ALLOCATED( tend_wind_front_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_wind_front_av(k,j,i) =                                               &
                            tend_wind_front_av(k,j,i) + tend_wind_front(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'trm_tend_wind_side' )
             IF ( ALLOCATED( tend_wind_side_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_wind_side_av(k,j,i) =                                                &
                            tend_wind_side_av(k,j,i) + tend_wind_side(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'trm_tend_wind_top' )
             IF ( ALLOCATED( tend_wind_top_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                        tend_wind_top_av(k,j,i) = tend_wind_top_av(k,j,i) + tend_wind_top(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'trm_tend_theta' )
             IF ( ALLOCATED( tend_theta_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_theta_av(k,j,i) = tend_theta_av(k,j,i) + tend_theta(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'trm_tend_wind_front' )
             IF ( ALLOCATED( tend_wind_front_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_wind_front_av(k,j,i) =                                               &
                            tend_wind_front_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'trm_tend_wind_side' )
             IF ( ALLOCATED( tend_wind_side_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_wind_side_av(k,j,i) =                                                &
                            tend_wind_side_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'trm_tend_wind_top' )
             IF ( ALLOCATED( tend_wind_top_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_wind_top_av(k,j,i) =                                                 &
                            tend_wind_top_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'trm_tend_theta' )
             IF ( ALLOCATED( tend_theta_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         tend_theta_av(k,j,i) =                                                    &
                            tend_theta_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE trm_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Traffic actions - call for all grid points.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_actions( location )

!
!-- input parameters
    CHARACTER(LEN=*) ::  location   !< action type
!
!-- local variables
    INTEGER(iwp)     ::  i          !< running index
    INTEGER(iwp)     ::  j          !< running index

!
!-- Check the traffic forcing is enabled.
    IF ( .NOT. traffic )  RETURN

!
!-- Check existence of the traffic in the local mpi subdomain.
    IF ( trm_nsl == 0 )  RETURN

!
!-- Log cpu time.
    CALL cpu_log( log_point_s(37), 'traffic', 'start' )

!
!-- Traffic module actions.
    SELECT CASE ( location )

       CASE ( 'before_timestep' )
!
!--       Actions before every timestep - zero traffic tendencies.
          tend_wind_front = 0.0_wp
          tend_wind_side = 0.0_wp
          tend_wind_top = 0.0_wp
          tend_theta = 0.0_wp
!
!--       Check and read car intensities for given timestamp.
          CALL read_trm_intensity_ts()

       CASE ( 'u-tendency' )
!
!--       Calculate and apply traffic u-tendency.
          DO i = nxl, nxr
             DO  j = nys, nyn
                CALL trm_apply_tend_car(i, j, 0)
             ENDDO
          ENDDO

       CASE ( 'v-tendency' )
!
!--       Calculate and apply traffic v-tendency.
          DO i = nxl, nxr
             DO  j = nys, nyn
                CALL trm_apply_tend_car(i, j, 1)
             ENDDO
          ENDDO

       CASE ( 'w-tendency' )
!
!--       Calculate and apply traffic w-tendency.
          DO i = nxl, nxr
             DO  j = nys, nyn
                CALL trm_apply_tend_car(i, j, 2)
             ENDDO
          ENDDO

       CASE ( 'pt-tendency' )
!
!--       Calculate and apply traffic pt-tendency.
          DO i = nxl, nxr
             DO  j = nys, nyn
                CALL trm_apply_tend_car(i, j, 3)
             ENDDO
          ENDDO

       CASE DEFAULT
          CONTINUE

    END SELECT
!
!-- Finish logging of cpu time.
    CALL cpu_log( log_point_s(37), 'traffic', 'stop' )

 END SUBROUTINE trm_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Traffic actions - call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_actions_ij( i, j, location )

!
!-- input parameters
    CHARACTER(LEN=*) ::  location   !< action type

    INTEGER(iwp)     ::  i          !< running index
    INTEGER(iwp)     ::  j          !< running index

!
!-- Check the traffic forcing is enabled.
    IF ( .NOT. traffic )  RETURN

!
!-- Check existence of the traffic in the local mpi subdomain.
    IF ( trm_nsl == 0 )  RETURN

!
!-- Log cpu time.
    CALL cpu_log( log_point_s(37), 'traffic', 'start' )

!
!-- Traffic module actions.
    SELECT CASE ( location )

       CASE ( 'before_timestep' )
!
!--       Actions before every timestep - zero traffic tendencies.
          tend_wind_front(:,j,i) = 0.0_wp
          tend_wind_side(:,j,i) = 0.0_wp
          tend_wind_top(:,j,i) = 0.0_wp
          tend_theta(:,j,i) = 0.0_wp
!
!--       Check and read car intensities for given timestamp.
          CALL read_trm_intensity_ts()

       CASE ( 'u-tendency' )
!
!--       Calculate and apply traffic u-tendency.
          CALL trm_apply_tend_car(i, j, 0)

       CASE ( 'v-tendency' )
!
!--       Calculate and apply traffic v-tendency.
          CALL trm_apply_tend_car(i, j, 1)


       CASE ( 'w-tendency' )
!
!--       Calculate and apply traffic w-tendency.
          CALL trm_apply_tend_car(i, j, 2)

       CASE ( 'pt-tendency' )
!
!--       Calculate and apply traffic pt-tendency.
          CALL trm_apply_tend_car(i, j, 3)

       CASE DEFAULT
          CONTINUE

    END SELECT

!
!-- Finish logging of cpu time.
    CALL cpu_log( log_point_s(37), 'traffic', 'stop' )

 END SUBROUTINE trm_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check Traffic parameters control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_check_parameters

!
!-- Theta tendency cannot be applied in neutral stratification simulation.
    IF ( neutral  .AND.  theta_tendency )  THEN
       message_string = 'Theta tendency cannot be applied for neutral stratification.' //          &
                        '&Check the configuration.'
       CALL message( 'trm_check_parameters', 'TRM0004', 1, 2, 0, 6, 0 )
       theta_tendency = .FALSE.
    ENDIF

 END SUBROUTINE trm_check_parameters



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of Traffic defined output quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_check_data_output( variable, unit )

    CHARACTER (LEN=*), INTENT(IN)  ::  variable  !< variable name
    CHARACTER (LEN=*), INTENT(OUT) ::  unit      !< unit of the quantity

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'trm_tend_wind_front' )
          unit = 'm/s2'

       CASE ( 'trm_tend_wind_side' )
          unit = 'm/s2'

       CASE ( 'trm_tend_wind_top' )
          unit = 'm/s2'

       CASE ( 'trm_tend_theta' )
          unit = 'K/s'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE trm_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Traffic output quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_data_output_3d( av, variable, found, local_pf, resorted, nzb_do, nzt_do )


    CHARACTER(LEN=*), INTENT(IN) ::  variable  !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av            !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do        !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN) ::  nzt_do        !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL, INTENT(INOUT)   ::  found         !< flag if output variable is found
    LOGICAL, INTENT(INOUT)   ::  resorted      !< flag if output is resorted

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf  !< local array
                                                            !< to which output data is resorted to
    INTEGER(iwp) ::  i      !< running index
    INTEGER(iwp) ::  j      !< running index
    INTEGER(iwp) ::  k      !< running index

    found = .TRUE.
    resorted = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'trm_tend_wind_front' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_wind_front(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_wind_front_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'trm_tend_wind_side' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_wind_side(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_wind_side_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'trm_tend_wind_top' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_wind_top(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_wind_top_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'trm_tend_theta' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_theta(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = tend_theta_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          found    = .FALSE.
          resorted = .FALSE.

    END SELECT

 END SUBROUTINE trm_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Output of the Traffic values into masked output file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_data_output_mask( av, variable, found, local_pf, mid )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp)         ::  av              !< averaging output
    INTEGER(iwp)         ::  mid             !< masked output running index
    INTEGER(iwp)         ::  i               !< running index
    INTEGER(iwp)         ::  im              !< running index
    INTEGER(iwp)         ::  j               !< running index
    INTEGER(iwp)         ::  jm              !< running index
    INTEGER(iwp)         ::  k               !< running index
    INTEGER(iwp)         ::  kk              !< running index

    LOGICAL              ::  found           !< variable found

    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  local_pf  !<

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )


       CASE ( 'trm_tend_wind_front' )
          IF ( av == 0 )  THEN
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_wind_front(mask_k(mid,k),mask_j(mid,j),            &
                                           mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_wind_front(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ELSE
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_wind_front_av(mask_k(mid,k),mask_j(mid,j),         &
                                           mask_i(mid,i) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_wind_front_av(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       CASE ( 'trm_tend_wind_side' )
          IF ( av == 0 )  THEN
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_wind_side(mask_k(mid,k),mask_j(mid,j),             &
                                           mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_wind_side(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ELSE
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_wind_side_av(mask_k(mid,k),mask_j(mid,j),          &
                                           mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_wind_side_av(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       CASE ( 'trm_tend_wind_top' )
          IF ( av == 0 )  THEN
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_wind_top(mask_k(mid,k),mask_j(mid,j),              &
                                            mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_wind_top(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ELSE
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_wind_top_av(mask_k(mid,k),mask_j(mid,j),           &
                                           mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_wind_top_av(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       CASE ( 'trm_tend_theta' )
          IF ( av == 0 )  THEN
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_theta(mask_k(mid,k),mask_j(mid,j),                 &
                                           mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_theta(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ELSE
             IF ( .NOT. mask_surface(mid) )  THEN
!
!--             Default masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend_theta_av(mask_k(mid,k),mask_j(mid,j),              &
                                           mask_i(mid,i))
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
!
!--             Terrain-following masked output.
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)

                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)

                      DO  k = 1, mask_size_l(mid,3)
!
!--                      Calculate the vertical index.
                         kk = MIN( topo_top_ind(jm,im,6) + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
!
!--                         Save output array.
                            local_pf(i,j,k) =  tend_theta_av(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE trm_data_output_mask


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which Traffic output quantities are defined.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

!
!-- Input variables
    CHARACTER (LEN=*), INTENT(IN)  ::  variable  !< Name of output variable.
!
!-- Output variables
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x    !< x grid of output variable.
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y    !< y grid of output variable.
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z    !< z grid of output variable.

    LOGICAL, INTENT(OUT)           ::  found     !< Flag if output var is found.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'trm_tend_wind_front', 'trm_tend_wind_side',                                         &
              'trm_tend_wind_top', 'trm_tend_theta' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE trm_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize parameters, variables, and arrays for Traffic module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_init_arrays

!
!-- General characteristics.
    CHARACTER(LEN=*), PARAMETER ::  input_nc_file         = 'PIDS_TRAFFIC'    !< traffic input file
!
!-- Dimensions in the input file.
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_car_type       = 'car_type'        !< car type dimension
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_s              = 's'               !< s dimension
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_time           = 'time'            !< time dimension

!
!-- Variables in the input file.
    CHARACTER(LEN=*), PARAMETER ::  nc_var_car_cd         = 'car_cd'          !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_car_height     = 'car_height'      !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_car_length     = 'car_length'      !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_car_width      = 'car_width'       !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_dirx           = 'dirx'            !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_diry           = 'diry'            !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_frac           = 'frac'            !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_i              = 'i'               !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_j              = 'j'               !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_k              = 'k'               !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_slope          = 'slope'           !<  variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_timestamp      = 'timestamp'       !< timestamp variable
    CHARACTER(LEN=*), PARAMETER ::  nc_var_width          = 'width'           !<  variable
!
!-- Local indices and variables.
    INTEGER(iwp)                             ::  l             !< loop index (local)
    INTEGER(iwp)                             ::  m             !< loop index (global)
    INTEGER(iwp)                             ::  nc_stat       !< result of nc operation
    INTEGER(iwp)                             ::  varid         !< variable ID
!
!-- Status of traffic input reading.
    LOGICAL                                  ::  input_pids_traffic  = .FALSE.   !< input reading flag
!
!-- Temporary buffers for timestamps, coordinates i,j,k and variables.
    CHARACTER, ALLOCATABLE, DIMENSION(:,:)   ::  buf_timestamp !< buffer variable
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)  ::  ibuf          !< buffer for reading i variable
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)  ::  jbuf          !< buffer for reading j variable
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)  ::  kbuf          !< buffer for reading k variable
    REAL(wp),     ALLOCATABLE, DIMENSION(:)  ::  realbuf       !< buffer for reading real variables

    IF ( .NOT. traffic )  RETURN

#if defined ( __netcdf )

    message_string = 'Allocation and init of traffic intput arrays.'
    CALL location_message(message_string, 'start' )
!
!-- Check if traffic PIDS file exists.
    INQUIRE( FILE = TRIM( input_nc_file )  // TRIM( coupling_char ), EXIST = input_pids_traffic  )
    IF ( .NOT. input_pids_traffic )  THEN
       message_string = 'Traffic file does not exists for this domain.' //                         &
                  '&Continue without vehicle induced momentum and heat.'
       CALL message( 'trm_check_parameters', 'TRM0001', 1, 2, 0, 6, 0 )
       traffic = .FALSE.
       RETURN
    ENDIF

!
!-- Open file and read timestamps.
    CALL open_read_file( TRIM( input_nc_file ) // TRIM( coupling_char ), ncid )

!
!-- Grab timestamps from netCDF file into a buffer.
!-- Note: buffer will be allocated in get_variable().
    CALL get_dimension_length( ncid, nts, nc_dim_time )
    nc_stat = NF90_INQ_VARID( ncid, TRIM(nc_var_timestamp), varid )
    ALLOCATE(buf_timestamp(date_time_str_len,nts))
    nc_stat = NF90_GET_VAR( ncid, varid, buf_timestamp,                                            &
                           start=(/1,1/), count=(/date_time_str_len,nts/) )

    ALLOCATE( trm_ts( nts ) )
    DO  l = 1, nts
       DO m = 1, date_time_str_len
          trm_ts(l)(m:m) = buf_timestamp(m,l)
       ENDDO
    ENDDO

!
!-- Deallocate array allocated in get_variable.
    DEALLOCATE( buf_timestamp )

!
!-- Car_type dimension and variable.
    CALL get_dimension_length( ncid, nct, nc_dim_car_type )

!
!-- Other car type parameters.
    ALLOCATE(car(nct))
    ALLOCATE(realbuf(nct))
    CALL get_variable( ncid, nc_var_car_cd, realbuf )
    DO l = 1, nct
       car(l)%cd = realbuf(l)
    ENDDO
    CALL get_variable( ncid, nc_var_car_length, realbuf )
    DO l = 1, nct
       car(l)%length = realbuf(l)
    ENDDO
    CALL get_variable( ncid, nc_var_car_width, realbuf )
    DO l = 1, nct
       car(l)%width = realbuf(l)
    ENDDO
    CALL get_variable( ncid, nc_var_car_height, realbuf )
    DO l = 1, nct
       car(l)%height = realbuf(l)
    ENDDO
    DEALLOCATE(realbuf)

!
!-- Get grid traffic index dimension and its size.
    CALL get_dimension_length( ncid, trm_ns, nc_dim_s )

!
!-- Allocate temporary i,j,k buffers.
    ALLOCATE(ibuf(trm_ns))
    ALLOCATE(jbuf(trm_ns))
    ALLOCATE(kbuf(trm_ns))

!
!-- Read i,j,k.
    CALL get_variable( ncid, nc_var_i, ibuf )
    CALL get_variable( ncid, nc_var_j, jbuf )
    CALL get_variable( ncid, nc_var_k, kbuf )

!
!-- Calculate number of elements inside the mpi subdomain.
    l = 0
    DO m = 1, trm_ns
       IF ( nxl <= ibuf(m)  .AND.  nxr >= ibuf(m)  .AND.  nys <= jbuf(m)  .AND.  nyn >= jbuf(m) )  THEN
          l = l + 1
       ENDIF
    ENDDO
    trm_nsl = l

!
!-- Allocate and assign trm array and array indexing from local to global indexex.
    ALLOCATE( trm(trm_nsl) )
    l = 0
    DO m = 1, trm_ns
       IF ( nxl <= ibuf(m)  .AND.  nxr >= ibuf(m)  .AND.  nys <= jbuf(m)  .AND.  nyn >= jbuf(m) )  THEN
!
!--       Increase index and add to index array.
          l = l + 1
          trm(l)%l2g = m
          trm(l)%i = ibuf(m)
          trm(l)%j = jbuf(m)
          trm(l)%k = kbuf(m)
!
!--       Allocate arrays for car intensiy, speed and heat
          ALLOCATE( trm(l)%intensity(nct) )
          ALLOCATE( trm(l)%speed(nct) )
          ALLOCATE( trm(l)%heat(nct) )
       ENDIF
    ENDDO

!
!-- Sort array trm by i,j
    CALL trm_quicksort( trm, 1, trm_nsl )

!
!-- Allocate index array and store position of the first element of the i,j grid block of the traffic.
    ALLOCATE( trm_indl(nys:nyn,nxl:nxr) )
    trm_indl = 0
    DO l = 1, trm_nsl
          IF ( trm_indl(trm(l)%j,trm(l)%i) == 0 )  THEN
             trm_indl(trm(l)%j,trm(l)%i) = l
          ENDIF
    ENDDO

!
!-- Read and assign all variables.
    ALLOCATE(realbuf(trm_ns))
    CALL get_variable( ncid, nc_var_width, realbuf )
    DO l = 1, trm_nsl
       trm(l)%width = realbuf(trm(l)%l2g)
    ENDDO
    CALL get_variable( ncid, nc_var_slope, realbuf )
    DO l = 1, trm_nsl
       trm(l)%slope = realbuf(trm(l)%l2g)
    ENDDO
    CALL get_variable( ncid, nc_var_dirx, realbuf )
    DO l = 1, trm_nsl
       trm(l)%dirx = realbuf(trm(l)%l2g)
    ENDDO
    CALL get_variable( ncid, nc_var_diry, realbuf )
    DO l = 1, trm_nsl
       trm(l)%diry = realbuf(trm(l)%l2g)
    ENDDO
    CALL get_variable( ncid, nc_var_frac, realbuf )
    DO l = 1, trm_nsl
       trm(l)%frac = realbuf(trm(l)%l2g)
    ENDDO
    DEALLOCATE(realbuf)

    CALL location_message(message_string, 'finished' )
!
!-- Initialize timestamp indices.
    trm_current_ts = ''
    trm_ts_index = 0
    trm_next_ts = ''

!
!-- Allocate output 3d arrays.
    message_string = 'Allocation and init of traffic output arrays.'
    CALL location_message(message_string, 'start' )
    IF ( .NOT. ALLOCATED(tend_wind_front) )  THEN
       ALLOCATE( tend_wind_front(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_wind_front = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED(tend_wind_side) )  THEN
       ALLOCATE( tend_wind_side(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_wind_side = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED(tend_wind_top) )  THEN
       ALLOCATE( tend_wind_top(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_wind_top = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED(tend_theta) )  THEN
       ALLOCATE( tend_theta(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_theta = 0.0_wp
    ENDIF
    CALL location_message(message_string, 'finished' )

#else
    message_string = 'Vehicle induced momentum and heat requires compilation with netcdf.' //      &
                     '&Continue without vehicle induced momentum and heat.'
    CALL message( 'trm_init_arrays', 'TRM0002', 1, 2, 0, 6, 0 )
    traffic = .FALSE.
    RETURN
#endif

 END SUBROUTINE trm_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Traffic actions at the end of a job.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_last_actions

    INTEGER          ::  nc_stat  !<

#if defined ( __netcdf )
!
!-- Close the traffic input file.
    nc_stat = NF90_CLOSE( ncid )
#endif

 END SUBROUTINE trm_last_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for parameters of Traffic module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_parin

    CHARACTER (LEN=80) ::  line       !< String containing the last line read from namelist file.
    INTEGER(iwp)       ::  io_status  !< Status after reading the namelist file.

    LOGICAL ::  switch_off_module = .FALSE.  !< Local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file.

    NAMELIST /traffic_parameters/     continuous_tendency,                                         &
                                      switch_off_module,                                           &
                                      theta_tendency,                                              &
                                      wind_tendency

!
!-- Position the namelist-file at the beginning (it has already been opened in parin), and try to
!-- read (find) a namelist named "traffic_parameters".
    REWIND ( 11 )
    READ( 11, traffic_parameters, IOSTAT=io_status )

!
!-- Actions depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    Traffic namelist found and correctly read. Set default module switch to true.
!--    This activates calls of the Traffic interface subroutines.
       IF ( switch_off_module  .OR.  ( .NOT. wind_tendency  .AND.  .NOT. theta_tendency ) )  THEN
         traffic = .FALSE.
       ELSE
          traffic = .TRUE.
       ENDIF

    ELSEIF ( io_status > 0 )  THEN
!
!--    Traffic namelist was found, but contained errors. Print an error message containing the line
!--    that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'traffic_parameters', line )

    ENDIF

 END SUBROUTINE trm_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf,          &
                               nync, nyn_on_file, nysf, nysc, nys_on_file, tmp_3d, found )

    INTEGER(iwp) ::  k               !< running index over previous input files covering current local domain
    INTEGER(iwp) ::  nxlc            !< index of left boundary on current subdomain
    INTEGER(iwp) ::  nxlf            !< index of left boundary on former subdomain
    INTEGER(iwp) ::  nxl_on_file     !< index of left boundary on former local domain
    INTEGER(iwp) ::  nxrc            !< index of right boundary on current subdomain
    INTEGER(iwp) ::  nxrf            !< index of right boundary on former subdomain
    INTEGER(iwp) ::  nxr_on_file     !< index of right boundary on former local domain
    INTEGER(iwp) ::  nync            !< index of north boundary on current subdomain
    INTEGER(iwp) ::  nynf            !< index of north boundary on former subdomain
    INTEGER(iwp) ::  nyn_on_file     !< index of north boundary on former local domain
    INTEGER(iwp) ::  nysc            !< index of south boundary on current subdomain
    INTEGER(iwp) ::  nysf            !< index of south boundary on former subdomain
    INTEGER(iwp) ::  nys_on_file     !< index of south boundary on former local domain

    LOGICAL, INTENT(OUT)  ::  found  !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d  !<

!
!-- Reading of Traffic restart data:

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'trm_tend_wind_front_av' )
          IF ( .NOT. ALLOCATED( tend_wind_front_av ) )                                             &
             ALLOCATE( tend_wind_front_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          tend_wind_front_av = 0.0_wp
          IF ( k == 1 )  READ ( 13 )  tmp_3d(:,:,:)
          tend_wind_front_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'trm_tend_wind_side_av' )
          IF ( .NOT. ALLOCATED( tend_wind_side_av) )                                               &
             ALLOCATE( tend_wind_side_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          tend_wind_side_av = 0.0_wp
          IF ( k == 1 )  READ ( 13 )  tmp_3d(:,:,:)
          tend_wind_side_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'trm_tend_wind_top_av' )
          IF ( .NOT. ALLOCATED( tend_wind_top_av ) )                                               &
             ALLOCATE( tend_wind_top_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          tend_wind_top_av = 0.0_wp
          IF ( k == 1 )  READ ( 13 )  tmp_3d(:,:,:)
          tend_wind_top_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'trm_tend_theta_av' )
          IF ( .NOT. ALLOCATED( tend_theta_av ) )                                                  &
             ALLOCATE( tend_theta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          tend_theta_av = 0.0_wp
          IF ( k == 1 )  READ ( 13 )  tmp_3d(:,:,:)
          tend_theta_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE trm_rrd_local_ftn

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_rrd_local_mpi

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array, rrd_mpi_io

    LOGICAL  ::  array_found     !< flag

!
!-- Allocate average output 3d arrays.
    CALL rd_mpi_io_check_array( 'trm_tend_wind_front_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( tend_wind_front_av ) )                                                &
          ALLOCATE( tend_wind_front_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_wind_front_av = 0.0_wp
       CALL rrd_mpi_io( 'trm_tend_wind_front_av', tend_wind_front_av )
    ENDIF
    CALL rd_mpi_io_check_array( 'trm_tend_wind_side_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( tend_wind_side_av ) )                                                 &
          ALLOCATE( tend_wind_side_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_wind_side_av = 0.0_wp
       CALL rrd_mpi_io( 'trm_tend_wind_side_av', tend_wind_side_av )
    ENDIF
    CALL rd_mpi_io_check_array( 'trm_tend_wind_top_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( tend_wind_top_av ) )                                                  &
          ALLOCATE( tend_wind_top_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_wind_top_av = 0.0_wp
       CALL rrd_mpi_io( 'trm_tend_wind_top_av', tend_wind_top_av )
    ENDIF
    CALL rd_mpi_io_check_array( 'trm_tend_theta_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( tend_theta_av ) )                                                     &
          ALLOCATE( tend_theta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       tend_theta_av = 0.0_wp
       CALL rrd_mpi_io( 'trm_tend_theta_av', tend_theta_av )
    ENDIF

 END SUBROUTINE trm_rrd_local_mpi

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes Traffic defined restart data into binary/mpi file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_wrd_local

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  wrd_mpi_io

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( tend_wind_front_av ) )  THEN
          CALL wrd_write_string( 'trm_tend_wind_front_av' )
          WRITE ( 14 )  tend_wind_front_av
       ENDIF

       IF ( ALLOCATED( tend_wind_side_av ) )  THEN
          CALL wrd_write_string( 'trm_tend_wind_side_av' )
          WRITE ( 14 )  tend_wind_side_av
       ENDIF

       IF ( ALLOCATED( tend_wind_top_av ) )  THEN
          CALL wrd_write_string( 'trm_tend_wind_top_av' )
          WRITE ( 14 )  tend_wind_top_av
       ENDIF

       IF ( ALLOCATED( tend_theta_av ) )  THEN
          CALL wrd_write_string( 'trm_tend_theta_av' )
          WRITE ( 14 )  tend_theta_av
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       IF ( ALLOCATED( tend_wind_front_av ) )                                                      &
          CALL wrd_mpi_io( 'trm_tend_wind_front_av', tend_wind_front_av )
       IF ( ALLOCATED( tend_wind_side_av ) )                                                       &
          CALL wrd_mpi_io( 'trm_tend_wind_side_av', tend_wind_side_av )
       IF ( ALLOCATED( tend_wind_top_av ) )                                                        &
          CALL wrd_mpi_io( 'trm_tend_wind_top_av', tend_wind_top_av )
       IF ( ALLOCATED( tend_wind_top_av ) )                                                        &
          CALL wrd_mpi_io( 'trm_tend_wind_top_av', tend_wind_top_av )
       IF ( ALLOCATED( tend_theta_av ) )                                                           &
          CALL wrd_mpi_io( 'trm_tend_theta_av', tend_theta_av )

    ENDIF

 END SUBROUTINE trm_wrd_local


!--------------------------------------------------------------------------------------------------!
!> Section of the private module subroutines.
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads traffic intensities from file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE read_trm_intensity_ts()

!
!-- nc variable names
    CHARACTER(LEN=*), PARAMETER   ::  nc_var_intensity   = 'intensity'
    CHARACTER(LEN=*), PARAMETER   ::  nc_var_heat        = 'heat'
    CHARACTER(LEN=*), PARAMETER   ::  nc_var_speed       = 'speed'
!
!-- Local variables
    CHARACTER(LEN=date_time_str_len)  ::  timestamp  !< current time stamp

    INTEGER(iwp)                      ::  l          !< running index
    INTEGER(iwp)                      ::  ict        !< running index
    INTEGER(iwp)                      ::  nc_stat    !< results of nc operation
    INTEGER(iwp)                      ::  varid      !< variable ID

    LOGICAL                           ::  found

!
!-- Read buffer for time slice of the  2-d variable with dimension 'ntime', 's', 'car_type' (ntime fixed).
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)  ::  real2dbuf       !< read buffer

#if defined( __netcdf )
!
!-- Current timestamp.
    CALL get_date_time( time_since_reference_point, date_time_str=timestamp )

!
!-- Check timestamp
    IF ( timestamp < trm_next_ts  .OR.  trm_ts_index >= nts )  THEN
!
!--    Valid timestamp already read.
       RETURN
    ENDIF
!
!-- Find the new timestamp position and index (require ordered timestamp in nc file).
    WRITE( message_string, *) 'Reading new timestep from traffic intput data:', trm_ts_index, trm_next_ts
    CALL location_message(message_string, 'start' )
    found = .FALSE.
    DO l = trm_ts_index + 1, nts
       IF ( trm_ts(l) > timestamp )  THEN
          found = .TRUE.
          trm_ts_index = l - 1
          EXIT
       ENDIF
    ENDDO
    IF ( .NOT. found )  THEN
       trm_ts_index = nts
    ENDIF
    trm_current_ts = trm_ts(trm_ts_index)
    IF ( trm_ts_index < nts )  THEN
       trm_next_ts = trm_ts(trm_ts_index + 1)
    ELSE
       trm_next_ts = trm_ts(trm_ts_index)
    ENDIF
!
!-- Real 2d buffer for nc variable.
    ALLOCATE(real2dbuf(nct, trm_ns))
!
!-- Read intensity and transform to the local variable.
!-- PALM netcdf interface does not contain needed function - use netcdf interface directly.
    nc_stat = NF90_INQ_VARID( ncid, TRIM( nc_var_intensity ), varid )
    nc_stat = NF90_GET_VAR( ncid, varid, real2dbuf,                                                &
                            start=(/1,1,trm_ts_index/), count=(/nct,trm_ns,1/) )
    DO l = 1, trm_nsl
       DO ict = 1, nct
          trm(l)%intensity(ict) = real2dbuf(ict,trm(l)%l2g)
       ENDDO
    ENDDO

!
!-- Read heat and transform to the local variable.
    nc_stat = NF90_INQ_VARID( ncid, trim(nc_var_heat), varid )
    nc_stat = NF90_GET_VAR( ncid, varid, real2dbuf,                                                &
                           start=(/1,1,trm_ts_index/), count=(/nct,trm_ns,1/) )
    DO l = 1, trm_nsl
       DO ict = 1, nct
          trm(l)%heat(ict) = real2dbuf(ict,trm(l)%l2g)
       ENDDO
    ENDDO

!
!-- Read speed and transform to the local variable.
    nc_stat = NF90_INQ_VARID( ncid, trim(nc_var_speed), varid )
    nc_stat = NF90_GET_VAR( ncid, varid, real2dbuf,                                                &
                           start=(/1,1,trm_ts_index/), count=(/nct,trm_ns,1/) )
    DO l = 1, trm_nsl
       DO ict = 1, nct
          trm(l)%speed(ict) = real2dbuf(ict,trm(l)%l2g)
       ENDDO
    ENDDO

    DEALLOCATE(real2dbuf)

    WRITE( message_string, *) 'Reading new timestep from traffic intput data:',                    &
                              trm_ts_index, trm_next_ts
    CALL location_message(message_string, 'finished' )

#endif

 END SUBROUTINE read_trm_intensity_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates and applies wind velocity tendences from cars.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trm_apply_tend_car(i, j, tend_type)
!
    INTEGER(iwp), INTENT(in) :: i
    INTEGER(iwp), INTENT(in) :: j
    INTEGER(iwp), INTENT(in) :: tend_type

    INTEGER(iwp)   ::  car_layers       !< number of layers crossed by car
    INTEGER(iwp)   ::  ct               !< running index
    INTEGER(iwp)   ::  k                !< running index
    INTEGER(iwp)   ::  ktt              !< terrain height index
    INTEGER(iwp)   ::  l                !< running index

    REAL(wp)      :: car_interval       !< time interval between cars (s)
    REAL(wp)      :: nf_layer           !< relative percentage of the layer in the car
    REAL(wp)      :: nff                !< normalization factor (front)
    REAL(wp)      :: nfs                !< normalization factor (side)
    REAL(wp)      :: nfv                !< normalization factor (vertical)
    REAL(wp)      :: rcsf               !< relative car speed to wind (front component)
    REAL(wp)      :: rcss               !< relative car speed to wind (side component)
    REAL(wp)      :: rcsv               !< vertical wind tendency
    REAL(wp)      :: tend_car_front     !< velocity tendency induced by car in the front
    REAL(wp)      :: tend_car_heat      !< pt tendency induced by car
    REAL(wp)      :: tend_car_side      !< velocity tendency induced by car on the side
    REAL(wp)      :: tend_car_top       !< velocity tendency induced by car on the top
    REAL(wp)      :: time_coef          !< normalization time coefficient
    REAL(wp)      :: wsf                !< wind in the car coordinates (front component)
    REAL(wp)      :: wss                !< wind in the car coordinates (side component)

!
!-- Cycle over possible layers.
    l  = trm_indl(j,i)
    IF ( l == 0 )  RETURN
!
!-- Cycle over car streams going over traffic grids with [i,j] coordinates.
    DO WHILE  ( trm(l)%i == i  .AND.  trm(l)%j == j )
!
!--    Cycle over car types.
       DO ct = 1, nct
          IF ( trm(l)%intensity(ct) <= 1.0E-10_wp )  CYCLE
!
!--       Interval between cars of given type in i,j.
!--       Convert car_intensity (car/hour) -> car interval (seconds/car).
          car_interval = 3600.0_wp / trm(l)%intensity(ct)
          IF ( continuous_tendency )  THEN
!
!--          In case of the application of the tendency continuously, we need to apply folowing time_coef.
             time_coef = (car(ct)%length / trm(l)%speed(ct)) / car_interval
          ELSE
!
!--          Application of the tendency only during car presence.
!--          Check if car of given type is currently runing over i,j.
             IF ( .NOT. ( MODULO(time_since_reference_point, car_interval) <=                      &
                          car(ct)%length / trm(l)%speed(ct) ) )  CYCLE
!
!--          In case of the application of the tendency only during car presence, no time_coef is applied.
             time_coef = 1.0_wp
          ENDIF
!
!--       Terrain height
          ktt = topo_top_ind(j,i,0)
          car_layers = 1
          DO WHILE  ( ktt+car_layers <= nzt  .AND.                                                 &
                      zw(ktt+car_layers) - zw(ktt) < car(ct)%height )
             car_layers = car_layers + 1
          END DO
!
!--       Cycle over layers occupated by car (-1 in upper range as the indexing of trm(l)%k starts 1).
          DO k = ktt+trm(l)%k, ktt+trm(l)%k+car_layers-1
!
!--          Relative percentage of the layer in the car.
             nf_layer = MERGE( 1.0_wp, (car(ct)%height-(zw(ktt+car_layers-1)-zw(ktt))) /           &
                               dzw(ktt+car_layers) , k-ktt-trm(l)%k < car_layers-1 )
             IF ( tend_type == tend_pt )  THEN
                IF ( theta_tendency )  THEN
!
!--                Heat per car (W per car) * number of cars per second / rho / c_p /exner = pt/s (K/s)
!--                plus time and fraction normalizations.
                   tend_car_heat =                                                                 &
                      trm(l)%intensity(ct) * trm(l)%heat(ct) / 3600.0_wp /                         &
                      (rho_air(k) * trm(l)%width * car(ct)%height * car(ct)%length )               &
                      / c_p / exner(k) *                                                           &
                      time_coef * trm(l)%frac * nf_layer
!
!--                Apply tendency.
                   tend(k,j,i) = tend(k,j,i) + tend_car_heat
!
!--                Add to the output array.
                   tend_theta(k,j,i) = tend_theta(k,j,i) + tend_car_heat
                ENDIF
             ELSE
                IF ( wind_tendency )  THEN
!
!--                Fad = 0.5 · Cd · Scar · rhoair · (vcar - vwind)2
!--                tend_car = Fad / mair = Fad / rhoair / Scar / Lcar
!--                tend_car = 0.5 · Cd / Lcar * (vcar - vwind)2
                   IF ( tend_type == tend_wind_w )  THEN
!
!--                   Vertical wind tendency (calculation neglegts car inclining).
                      rcsv = trm(l)%speed(ct) * trm(l)%slope * 0.01_wp - w(k,j,i)
                      nfv = time_coef * trm(l)%frac * nf_layer *                                   &
                            car(ct)%width / trm(l)%width / car(ct)%height
                      tend_car_top = 0.5_wp * car(ct)%cd * rcsv * rcsv * nfv
!
!--                   Apply w-tend.
                      tend(k,j,i) = tend(k,j,i) + tend_car_top
!
!--                   Add to the output array.
                      tend_wind_top(k,j,i) = tend_wind_top(k,j,i) + tend_car_top
                   ELSE
!
!--                   Horizontal wind tendences.
!--                   Transformation of the wind into car coordinates.
                      wsf = u(k,j,i) * trm(l)%dirx + v(k,j,i) * trm(l)%diry
                      wss = u(k,j,i) * trm(l)%diry - v(k,j,i) * trm(l)%dirx
!
!--                   Relative car speed in front and side direction.
                      rcsf = trm(l)%speed(ct) - wsf
                      rcss = -wss
!
!--                   Normalization factors
!--                   (front/side area, width of cars/width of street, fraction and time normalization).
                      nff = time_coef * trm(l)%frac * nf_layer *                                   &
                            car(ct)%width / trm(l)%width / car(ct)%length
                      nfs = time_coef * trm(l)%frac * nf_layer / car(ct)%width
!
!--                   Front tendency in the car direction (= 1/2 c_d * cs^2 * Av/Al*lv).
                      tend_car_front = 0.5_wp * car(ct)%cd * rcsf * rcsf * nff
!
!--                   Side tendency in the direction orthogonal to car direction (90 deg. rotation).
                      tend_car_side = 0.5_wp * car(ct)%cd * rcss * rcss * nfs
!
!--                   Apply tend = tend + tend_car_front + tend_car_side.
!--                   It needs to transform back to the model coordinate system.
                      IF ( tend_type == tend_wind_u )  THEN
!
!--                      Apply u-tend front.
                         tend(k,j,i) = tend(k,j,i) + tend_car_front * trm(l)%dirx
!
!--                      Apply u-tend side.
                         tend(k,j,i) = tend(k,j,i) + tend_car_side * trm(l)%diry

                      ELSE IF ( tend_type == tend_wind_v )  THEN
!
!--                      Apply v-tend front.
                         tend(k,j,i) = tend(k,j,i) + tend_car_front * trm(l)%diry
!
!--                      Apply v-tend side.
                         tend(k,j,i) = tend(k,j,i) - tend_car_side * trm(l)%dirx

                      ENDIF
!
!--                   Add to the output array.
                      tend_wind_front(k,j,i) = tend_wind_front(k,j,i) + tend_car_front
                      tend_wind_side(k,j,i) = tend_wind_side(k,j,i) + tend_car_side
                   ENDIF
!
!--                Control debugging prints in case of suspicious tendencies.
                   IF ( ABS(tend_car_front) >= 2.0_wp  .OR.  ABS(tend_car_side) >= 2.0_wp )  THEN
                      WRITE( message_string, * )                                                   &
                         'Wind tendencies from traffic exceed limit.',                             &
                         '&i, j, k, tend_type, ct, tend_car_front, tend_car_side, wind_tend:',     &
                         i, j, k, tend_type, ct, tend_car_front, tend_car_side, tend(k,j,i),       &
                         tend(k,j,i) + tend_car_front + tend_car_side
                      CALL message( 'trm_apply_tend_car', 'TRM0003', 0, 1, 0, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       l = l + 1
       IF ( l > trm_nsl )  EXIT
    ENDDO


 END SUBROUTINE trm_apply_tend_car


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Quicksort algorithm for sorting trm array according to i,j coordinates.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE trm_quicksort( trm, first, last )

    IMPLICIT NONE
!
!-- parameters of the subroutine
    INTEGER(iwp), INTENT(IN)                     :: first, last   !<
    TYPE(trm_type),DIMENSION(:), INTENT(INOUT)   :: trm           !<
!
!-- local variables
    INTEGER(iwp)      :: i      !< left index
    INTEGER(iwp)      :: j      !< right index
    TYPE(trm_type)    :: t      !< pivot
    TYPE(trm_type)    :: x      !< temporary storage exchange of elements

    IF ( first >= last )  RETURN
    x = trm(( first + last ) / 2)
    i = first
    j = last
    DO
       DO WHILE  ( trm_lt( trm(i), x ) )
          i = i + 1
       ENDDO
       DO WHILE  ( trm_lt( x, trm(j) ) )
          j = j - 1
       ENDDO
       IF ( i >= j )  EXIT
       t = trm(i);  trm(i) = trm(j);  trm(j) = t
       i = i+1
       j = j-1
    ENDDO
    IF ( first < i-1  )  CALL trm_quicksort( trm,  first, i - 1 )
    IF ( j+1 < last )  CALL trm_quicksort( trm, j + 1, last  )

 END SUBROUTINE trm_quicksort

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Comparison function for quicksort algorithm, which returns whether this is less than target.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION trm_lt( this, target ) RESULT( res )
!
!-- parameters of the function
    TYPE(trm_type), INTENT(IN) :: this  !<
    TYPE(trm_type), INTENT(IN) :: target  !<
!
!-- local variables
    LOGICAL :: res  !<

    IF ( this%i < target%i  .OR. ( this%i == target%i  .AND.  this%j < target%j) )  THEN
       res = .TRUE.
    ELSE
       res = .FALSE.
    ENDIF

 END FUNCTION trm_lt


 END MODULE traffic_mod
