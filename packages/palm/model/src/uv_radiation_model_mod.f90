!> @file uv_radiation_model_mod.f90
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
! Copyright 2022-2023 Pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Computation of erythemally-weighted UV-irradiance from LIB-RADTRAN-provided inputs.
!> In flat environments without any obstacles and plants, the modelled erythemally-weighted
!> UV-irradiance equals the externally-provided irradiance, while in urban environments the
!> effect of directional shading obstacles and plants, as well as multiple reflections are
!> considered.
!--------------------------------------------------------------------------------------------------!
 MODULE uv_radiation_model_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  degrees_to_radiants,                                                                &
               pi,                                                                                 &
               radiants_to_degrees

    USE control_parameters,                                                                        &
        ONLY: average_count_3d,                                                                    &
              coupling_char,                                                                       &
              debug_output,                                                                        &
              debug_output_timestep,                                                               &
              debug_string,                                                                        &
              dt_3d,                                                                               &
              end_time,                                                                            &
              dt_do2d_xy,                                                                          &
              length,                                                                              &
              message_string,                                                                      &
              plant_canopy,                                                                        &
              restart_data_format_output,                                                          &
              restart_string,                                                                      &
              rotation_angle,                                                                      &
              skip_time_do2d_xy,                                                                   &
              time_do2d_xy,                                                                        &
              time_since_reference_point,                                                          &
              uv_radiation

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE general_utilities,                                                                         &
        ONLY:  interpolate_linear

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
               nzb

    USE, INTRINSIC ::  IEEE_ARITHMETIC

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  char_fill,                                                                          &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               vars_pids

    USE palm_date_time_mod,                                                                        &
        ONLY:  date_time_str_len,                                                                  &
               get_date_time,                                                                      &
               hours_per_day,                                                                      &
               seconds_per_hour

    USE pegrid

    USE radiation_model_mod,                                                                       &
        ONLY:  albedo_pars,                                                                        &
               albedo_surf,                                                                        &
               calc_zenith,                                                                        &
               cos_zenith,                                                                         &
               discr_azim_cent,                                                                    &
               discr_azim_bdry,                                                                    &
               discr_elev_cent,                                                                    &
               discr_elev_bdry,                                                                    &
               dt_radiation,                                                                       &
               dsidir_rev,                                                                         &
               dsitrans,                                                                           &
               idir,                                                                               &
               id,                                                                                 &
               ix,                                                                                 &
               iy,                                                                                 &
               iz,                                                                                 &
               jdir,                                                                               &
               kdir,                                                                               &
               nsurf_type,                                                                         &
               nsurfl,                                                                             &
               nsvfl,                                                                              &
               radiation,                                                                          &
               radiation_interactions,                                                             &
               radiation_interactions_on,                                                          &
               raytrace_discrete_azims,                                                            &
               raytrace_discrete_elevs,                                                            &
               surfoutsl,                                                                          &
               svf,                                                                                &
               svfsurf,                                                                            &
               skyvf,                                                                              &
               skyvft,                                                                             &
               spherical_view,                                                                     &
               sun_dir_lat,                                                                        &
               sun_dir_lon,                                                                        &
               sun_direction,                                                                      &
               surfins,                                                                            &
               surfinsw,                                                                           &
               surfinswdir,                                                                        &
               surfinswdif,                                                                        &
               surfl,                                                                              &
               surfstart,                                                                          &
               va_az,                                                                              &
               va_z

#if defined( __parallel )
    USE radiation_model_mod,                                                                       &
        ONLY:  isurf_send_radx,                                                                    &
               disp_send_radx,                                                                     &
               disp_recv_radx,                                                                     &
               nsend_radx,                                                                         &
               surfouts_recv,                                                                      &
               rtm_alltoallv,                                                                      &
               radx_send
#endif

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io,                                                                         &
               wrd_mpi_io

    USE surface_mod,                                                                               &
        ONLY:  ind_pav_green,                                                                      &
               ind_veg_wall,                                                                       &
               ind_wat_win,                                                                        &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_out,                                                                           &
               surf_usm,                                                                           &
               vertical_surfaces_exist


    IMPLICIT NONE

    CHARACTER(LEN=50) ::  uv_integration_method = 'from_irradiance'  !< namelist parameter

    INTEGER(iwp) ::  day_of_year              !< day of year
    INTEGER(iwp) ::  num_reflections = 0      !< namelist parameter

    LOGICAL ::  calc_uv_ewir1 = .FALSE.       !< flag to trigger calculation of erythemally-weighted UV irradiation from irradiation input
    LOGICAL ::  calc_uv_ewir2 = .FALSE.       !< flag to trigger calculation of erythemally-weighted UV irradiation from radiance input
    LOGICAL ::  calc_uv_ir1 = .FALSE.         !< flag to trigger calculation of UV irradiation from irradiation input
    LOGICAL ::  calc_uv_ir2 = .FALSE.         !< flag to trigger calculation of UV irradiation from radiation input
    LOGICAL ::  input_pids_uv                 !< flag indicating whether an UV-input file exists
    LOGICAL ::  uv_from_irradiance = .FALSE.  !< control flag to integrate UV-exposure from provided irradiance
    LOGICAL ::  uv_from_radiance   = .FALSE.  !< control flag to integrate UV-exposure from provided angle dependent radiance

    REAL(wp)            ::  second_of_day                        !< actual second of the day

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_ewir1  !< instantaneous erythemally-weighted UV-irradiance, integrated from irradiance
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_ewir2  !< instantaneous erythemally-weighted UV-irradiance, integrated from radiance
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_ir1    !< instantaneous UV-irradiance, integrated from irradiance
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_ir2    !< instantaneous UV-irradiance, integrated from radiance
!
!-- Define data structure that compasses the input from the external UV-radiation model.
    TYPE uv_ext

       CHARACTER(LEN=7)  ::  input_file_uv = 'PIDS_UV'            !< name of UV-input file
       CHARACTER(LEN=10) ::  char_sza = 'sun_zenith'              !< input char for sun zenith angle
       CHARACTER(LEN=10) ::  char_wl = 'wavelength'               !< input char for wavelenghts
       CHARACTER(LEN=12) ::  char_vaa = 'view_azimuth'            !< input char for view azimuth angle
       CHARACTER(LEN=11) ::  char_vza = 'view_zenith'             !< input char for view zenith angle
       CHARACTER(LEN=18) ::  char_ir_diff = 'uv_diff_irradiance'  !< input char for diffusive part of spectral irradiance
       CHARACTER(LEN=17) ::  char_ir_dir = 'uv_dir_irradiance'    !< input char for direct part of spectral irradiance
       CHARACTER(LEN=11) ::  char_rad = 'uv_radiance'             !< input char for spectral radiance

       INTEGER(iwp) ::  ind_sz_l  !< lower index of the two closest matching SZA's in the input file
       INTEGER(iwp) ::  ind_sz_u  !< lower index of the two closest matching SZA's in the input file
       INTEGER(iwp) ::  n_aa      !< number of provided spherical azimuth angles
       INTEGER(iwp) ::  n_sza     !< number of provided sun-zenith angles
       INTEGER(iwp) ::  n_wl      !< number of provided wavelengths
       INTEGER(iwp) ::  n_za      !< number of provided spherical zenith angles

       LOGICAL ::  file_read_initialized = .FALSE.  !< flag to check whether the file input has been already initialized

       REAL(wp) ::  altitude          !< assumed altitude for the considered scenario
       REAL(wp) ::  d_lambda          !< wavelength interval
       REAL(wp) ::  day_of_year       !< day of year
       REAL(wp) ::  fill_diff         !< fill value for diffuse irradiance
       REAL(wp) ::  fill_dir          !< fill value for direct irradiance
       REAL(wp) ::  fill_rad          !< fill value for radiance
       REAL(wp) ::  mean_albedo       !< assumed albedo for the considered scenario
       REAL(wp) ::  ozone_depth       !< assumed ozone depth for the considered scenario
       REAL(wp) ::  phi_sun = 0.0_wp  !< considered azimuth position of the sun (by default south)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  d_view_azimuth        !< increment between two given azimuth angles
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  d_view_zenith         !< increment between two given zenith angles
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  d_wavelength          !< increment between two given wavelengths
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  sun_zenith            !< sun zenith angles provided by external radiation model
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wavelength            !< wavelengths provided by external radiation model
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  view_azimuth          !< spherical azimuth angles provided by external radiation model
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  view_azimuth_rot      !< spherical azimuth angles provided by external radiation model
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  view_zenith           !< cos of spherical zenith angles provided by external radiation model
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  weighting_coeff_ery   !< erythema weighting coefficients
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  weighting_coeff_vitd  !< Vitamin-D weighting coefficients

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_diff_irradiance  !< input array of spectral diffuse irradiance
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_dir_irradiance   !< input array of spectral direct irradiance

       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  uv_radiance       !< input array of spectral radiance
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  uv_radiance_diff  !< diffuse part of spectral radiance

    END TYPE uv_ext
!
!-- Define data structure that compasses location-dependent information.
    TYPE uv_ij

       REAL(wp) ::  uv_ery  !< erythemally weighted UV irradiance
       REAL(wp) ::  uvi     !< UV-index
       REAL(wp) ::  svf     !< integral sky-view factor without plants
       REAL(wp) ::  svf_p   !< integral sky-view factor with plants considered

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ir_spectral1  !< spectral irradiance including shading, computed from irradiance values
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ir_spectral2  !< spectral irradiance including shading, computed from radiance values

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_reflect  !< reflected portion, angle dependent

    END TYPE uv_ij
!
!-- Define data structure for UV-input data.
    TYPE( uv_ext ) ::  uv_in  !< array of derived data structure to store required information
!
!-- Define 2D array for local UV data.
    TYPE( uv_ij ), DIMENSION(:,:), ALLOCATABLE ::  uv_rad  !< array of derived data structure to store required information

    SAVE

    PRIVATE

!
!-- Public functions.
    PUBLIC uv_radiation_actions,                                                                   &
           uv_radiation_calculate,                                                                 &
           uv_radiation_check_data_output,                                                         &
           uv_radiation_check_data_output_surf,                                                    &
           uv_radiation_check_parameters,                                                          &
           uv_radiation_data_output_2d,                                                            &
           uv_radiation_data_output_surf,                                                          &
           uv_radiation_define_netcdf_grid,                                                        &
           uv_radiation_header,                                                                    &
           uv_radiation_init,                                                                      &
           uv_radiation_parin

    INTERFACE uv_radiation_actions
       MODULE PROCEDURE uv_radiation_actions
    END INTERFACE uv_radiation_actions

    INTERFACE uv_radiation_calculate
       MODULE PROCEDURE uv_radiation_calculate
    END INTERFACE uv_radiation_calculate

    INTERFACE uv_radiation_check_data_output
       MODULE PROCEDURE uv_radiation_check_data_output
    END INTERFACE uv_radiation_check_data_output

    INTERFACE uv_radiation_check_data_output_surf
       MODULE PROCEDURE uv_radiation_check_data_output_surf
    END INTERFACE uv_radiation_check_data_output_surf

    INTERFACE uv_radiation_check_parameters
       MODULE PROCEDURE uv_radiation_check_parameters
    END INTERFACE uv_radiation_check_parameters

    INTERFACE uv_radiation_data_output_2d
       MODULE PROCEDURE uv_radiation_data_output_2d
    END INTERFACE uv_radiation_data_output_2d

    INTERFACE uv_radiation_data_output_surf
       MODULE PROCEDURE uv_radiation_data_output_surf
    END INTERFACE uv_radiation_data_output_surf

    INTERFACE uv_radiation_define_netcdf_grid
       MODULE PROCEDURE uv_radiation_define_netcdf_grid
    END INTERFACE uv_radiation_define_netcdf_grid

     INTERFACE uv_radiation_header
       MODULE PROCEDURE uv_radiation_header
    END INTERFACE uv_radiation_header

    INTERFACE uv_radiation_init
       MODULE PROCEDURE uv_radiation_init
    END INTERFACE uv_radiation_init

    INTERFACE uv_radiation_parin
       MODULE PROCEDURE uv_radiation_parin
    END INTERFACE uv_radiation_parin

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_actions( location )

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string


    SELECT CASE ( location )

        CASE ( 'after_integration' )
!
!--        Computation of UV-irradiances is only required at output timesteps.
           IF ( time_do2d_xy >= dt_do2d_xy  .AND.                                                  &
                time_since_reference_point >= skip_time_do2d_xy )  THEN
              CALL uv_radiation_calculate( time_since_reference_point )
           ENDIF

        CASE ( 'do_integration_spinup' )
!
!--        Computation of UV-irradiances is only required at output timesteps.
           IF ( time_do2d_xy >= dt_do2d_xy )  THEN
              CALL uv_radiation_calculate( time_since_reference_point )
           ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE uv_radiation_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of UV-radiation from externally LIB-RADTRAN-provided input. The incoming UV-radiance
!> at the top of the urban layer is taken as input for local integration over the non-obscured
!> parts of the sky (directional shading of UV-radiation).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_calculate( reference_time )

    INTEGER(iwp) ::  dum_sz_l  !< dummy index for lower sun-zenith angle in input dataset
    INTEGER(iwp) ::  dum_sz_u  !< dummy index for upper sun-zenith angle in input dataset
    INTEGER(iwp) ::  i         !< running index in x-direction
    INTEGER(iwp) ::  iaz       !< running index over obscured azimuth angles
    INTEGER(iwp) ::  ind_surf  !< reference index in the surfl array
    INTEGER(iwp) ::  ind_az_l  !< lower sun-azimuth range to extract direct radiation portion from radiance
    INTEGER(iwp) ::  ind_az_u  !< upper sun-azimuth range to extract direct radiation portion from radiance
    INTEGER(iwp) ::  ind_ze_l  !< lower sun-zenith range to extract direct radiation portion from radiance
    INTEGER(iwp) ::  ind_ze_u  !< upper sun-zenith range to extract direct radiation portion from radiance
    INTEGER(iwp) ::  isd       !< reference number of precalculated sun positions
    INTEGER(iwp) ::  ize       !< running index over obscured zenith angles
    INTEGER(iwp) ::  j         !< running index in y-direction
    INTEGER(iwp) ::  naa       !< running index for view azimuth angle
    INTEGER(iwp) ::  nwl       !< running index for wavelength
    INTEGER(iwp) ::  nza       !< running index for view zenith angle
    INTEGER(iwp) ::  sza       !< running index over loaded sun-zenith angles

    INTEGER(iwp), DIMENSION(0:1) ::  ind_sun_zenith  !< array to temporarily store corresponding loaded sun-zenith angles

    REAL(wp) ::  dir_shade_fac   !< shade factor for direct beam
    REAL(wp) ::  fac_sza         !< interpolation factor
    REAL(wp) ::  reference_time  !< actual time
    REAL(wp) ::  sa              !< dummy value for solar azimuth expressed in the coordinate system of the radiation model
    REAL(wp) ::  solar_azimuth   !< current solar azimuth angle
    REAL(wp) ::  solar_zenith    !< current solar zenith angle
    REAL(wp) ::  trans           !< plant transmission factor
    REAL(wp) ::  val1            !< dummy value used for linear interpolation between two sun-zenith angles
    REAL(wp) ::  val2            !< dummy value used for linear interpolation between two sun-zenith angles
    REAL(wp) ::  val_diff        !< interpolated diffuse spectral irradiation
    REAL(wp) ::  val_dir         !< interpolated direct spectral irradiation

    REAL(wp), DIMENSION(0:uv_in%n_wl-1,0:1) ::  dir_part  !< extracted direct portion from radiance


    CALL cpu_log( log_point_s(28), 'uv model', 'start' )

!
!-- Calculate current solar position.
    CALL get_date_time( reference_time, day_of_year = day_of_year, second_of_day = second_of_day )

    sun_direction = .TRUE.
    CALL calc_zenith( day_of_year, second_of_day )
!
!-- Skip all actions at night. However, before, set all integrated irradiances
!-- back to zero, in order to avoid non-zero output values at night.
    IF ( cos_zenith <= 0.0 )  THEN
       IF ( ALLOCATED( uv_ir1 )   )  uv_ir1   = 0.0_wp
       IF ( ALLOCATED( uv_ir2 )   )  uv_ir2   = 0.0_wp
       IF ( ALLOCATED( uv_ewir1 ) )  uv_ewir1 = 0.0_wp
       IF ( ALLOCATED( uv_ewir2 ) )  uv_ewir2 = 0.0_wp
       RETURN
    ENDIF
!
!-- Determine actual solar zenith angle that need to be read from the external dataset.
    solar_zenith = ACOS( cos_zenith ) * radiants_to_degrees
!
!-- Also, determine actual solar azimuth angle. This is required to properly rotate the dataset.
!-- Solar azimuth is defined as following:
!-- sun in the south, west, north, east gives an angle of 0, 90, 180, 270 degrees, respectively.
    solar_azimuth = pi * radiants_to_degrees +                                                     &
                    ATAN2( sun_dir_lon, sun_dir_lat ) * radiants_to_degrees - rotation_angle
!
!-- Determine reference index for current sun-zenith angle.
    dum_sz_l = MINLOC( ABS( uv_in%sun_zenith - solar_zenith ), DIM = 1 ) - 1
    IF ( solar_zenith < uv_in%sun_zenith(dum_sz_l)  .AND.  dum_sz_l > 0 )  THEN
       dum_sz_l = dum_sz_l - 1
    ENDIF
    dum_sz_u = MIN( dum_sz_l + 1, uv_in%n_sza - 1 )
!
!-- Check if the dataset for the required sun-zenith angles has been already loaded. If not, input
!-- the relevant dataset.
    IF ( uv_in%ind_sz_l /= dum_sz_l  .OR.  uv_in%ind_sz_u /= dum_sz_u )  THEN
       uv_in%ind_sz_l = dum_sz_l
       uv_in%ind_sz_u = dum_sz_u

       CALL uv_radiation_rff_input
    ENDIF
!
!-- Determine the interpolation factor between two sun-zenith angles.
    fac_sza = 0.0_wp
    IF ( uv_in%ind_sz_l < uv_in%ind_sz_u )  THEN
       fac_sza = ( solar_zenith - uv_in%sun_zenith(uv_in%ind_sz_l) ) /                             &
                 ( uv_in%sun_zenith(uv_in%ind_sz_u) - uv_in%sun_zenith(uv_in%ind_sz_l) )
    ENDIF

!
!-- Compute spectral UV-irradiation from direct and diffuse part of uvspec data.
!-- Interpolate between the two nearest sun-zenith angles.
    IF ( uv_from_irradiance )  THEN
!
!--    First, identify solar direction vector (discretized number). Therefore, use the solar
!--    azimuth in the reference system as it is used in the radiation model.
       sa  = ATAN2( sun_dir_lon, sun_dir_lat ) * radiants_to_degrees - rotation_angle
       naa = FLOOR( ACOS( cos_zenith ) / pi * REAL( raytrace_discrete_elevs, KIND = wp ) )
       nza = MODULO( NINT( sa / 360.0_wp * REAL( raytrace_discrete_azims, KIND = wp )              &
                           - 0.5_wp, iwp ), raytrace_discrete_azims )
       isd = dsidir_rev(naa,nza)
!
!--    First, determine the shading factor for the direct radiation. Therefore, take the
!--    already pre-calcualated shading directions. To determine the respective indicies, re-use
!--    the same computation from the RTM.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             dir_shade_fac = dsitrans(spherical_view(j,i)%ind_surfl,isd)
             DO  nwl = 0, uv_in%n_wl-1

                val1 = uv_in%uv_diff_irradiance(0,nwl) * uv_rad(j,i)%svf_p +                       &
                       uv_in%uv_dir_irradiance(0,nwl)  * dir_shade_fac
                val2 = uv_in%uv_diff_irradiance(1,nwl) * uv_rad(j,i)%svf_p +                       &
                       uv_in%uv_dir_irradiance(1,nwl)  * dir_shade_fac

                uv_rad(j,i)%ir_spectral1(nwl) = interpolate_linear( val1, val2, fac_sza )

             ENDDO
          ENDDO
       ENDDO

    ENDIF
!
!-- Compute spectral UV-irradiation from radiance input. The radiance coming from angles
!-- near sun-position is assumed to be the direct part.
    IF ( uv_from_radiance )  THEN
!
!--    Set possibly faulty uvspec output (NaNs) to zero.
       WHERE( IEEE_IS_NAN( uv_in%uv_radiance )  .OR.  uv_in%uv_radiance == uv_in%fill_rad  )       &
          uv_in%uv_radiance = 0.0_wp
!
!--    Separate direct and diffuse radiance. The direct part in the radiance array comes from
!--    zero azimuth and solar_zenith. Only required when data was loaded.
       ind_sun_zenith(0)  = uv_in%ind_sz_l
       ind_sun_zenith(1)  = uv_in%ind_sz_u

       uv_in%uv_radiance_diff = uv_in%uv_radiance
       DO  sza = 0, 1
!
!--       To compute the direct radiation portion from the radiance, mask an area of 15 x 15 degrees
!--       azimuth and zenith representative for the sun position. This area coverage is empirically
!--       derived from uvspec-delivered radiances.
          ind_ze_l = MINLOC( ABS( uv_in%view_zenith                                                &
                                - ( uv_in%sun_zenith(ind_sun_zenith(sza)) - 15.0_wp ) ),           &
                             DIM = 1 ) - 1
          ind_ze_u = MINLOC( ABS( uv_in%view_zenith                                                &
                                - ( uv_in%sun_zenith(ind_sun_zenith(sza)) + 15.0_wp ) ),           &
                             DIM = 1 ) - 1
!
!--       For azimuth, consider turn-over to correctly consider the range between 0-359 degrees.
          ind_az_l = MINLOC( ABS( MODULO( uv_in%view_azimuth - 15.0_wp, 360.0_wp ) ),              &
                             DIM = 1 ) - 1
          ind_az_u = MINLOC( ABS( MODULO( uv_in%view_azimuth + 15.0_wp, 360.0_wp ) ),              &
                             DIM = 1 ) - 1
!
!--       Extract the direct portion and correspondingly set the radiance at these spherical angles
!--       to zero.
          DO  nwl = 0, uv_in%n_wl-1
             dir_part(nwl,sza) = 0.0_wp
             DO  nza = ind_ze_l, ind_ze_u
                DO  naa = ind_az_u, uv_in%n_aa-1
                   dir_part(nwl,sza) = dir_part(nwl,sza) + uv_in%uv_radiance(sza,nwl,nza,naa)      &
                                                         * uv_in%d_view_zenith(nza)                &
                                                         * uv_in%d_view_azimuth(naa)
                   uv_in%uv_radiance_diff(sza,nwl,nza,naa) = 0.0_wp
                ENDDO
                DO  naa = 0, ind_az_l
                   dir_part(nwl,sza) = dir_part(nwl,sza) + uv_in%uv_radiance(sza,nwl,nza,naa)      &
                                                         * uv_in%d_view_zenith(nza)                &
                                                         * uv_in%d_view_azimuth(naa)
                   uv_in%uv_radiance_diff(sza,nwl,nza,naa) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!
!--    Rotate view-azimuth according to the actual sun azimuth position. This is required because
!--    the external uvspec dataset assumes the sun position always at the south.
       uv_in%view_azimuth_rot = uv_in%view_azimuth - solar_azimuth
       uv_in%view_azimuth_rot = MODULO( uv_in%view_azimuth_rot, 360.0_wp )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             dir_shade_fac = dsitrans(spherical_view(j,i)%ind_surfl,isd)
             DO  nwl = 0, uv_in%n_wl-1
                val1 = dir_part(nwl,0) * dir_shade_fac
                val2 = dir_part(nwl,1) * dir_shade_fac

                uv_rad(j,i)%ir_spectral2(nwl) = interpolate_linear( val1, val2, fac_sza )
             ENDDO
!
!--          In a first step, integrate over all spherical angles. Resolution of spherical angles
!--          is determined by the external radiation model.
             DO  naa = 0, uv_in%n_aa-1
                DO  nza = 0, uv_in%n_za-1
                   DO  nwl = 0, uv_in%n_wl-1

                      val1 = uv_in%uv_radiance_diff(0,nwl,nza,naa)
                      val2 = uv_in%uv_radiance_diff(1,nwl,nza,naa)

                      uv_rad(j,i)%ir_spectral2(nwl) = uv_rad(j,i)%ir_spectral2(nwl)                &
                                                    + interpolate_linear( val1, val2, fac_sza )    &
                                                    * uv_in%d_view_zenith(nza)                     &
                                                    * uv_in%d_view_azimuth(naa)
                   ENDDO
                ENDDO
             ENDDO
!
!--          In a second step, remove radiance at obscured angles from the already integrated
!--          irradiance. Resolution of spherical angles is determined by RTM settings.
!--          First, check if there are any obscured azimuth angles.
             IF ( spherical_view(j,i)%treat_az )  THEN
!
!--             Loop over all azimuth angles that compass obscured zenith angles.
                DO  iaz = 1, spherical_view(j,i)%n_az
!
!--                Determine the reference view azimuth angle in the rotated data set.
                   naa = MINLOC( ABS( uv_in%view_azimuth_rot -                                     &
                                      spherical_view(j,i)%az(iaz)%val_az ), DIM = 1 ) - 1

                   IF ( spherical_view(j,i)%az(iaz)%treat_ze )  THEN
!
!--                   Loop over all non-transparently obscured zenith angles.
                      DO  ize = 1, spherical_view(j,i)%az(iaz)%n_blocked
!
!--                      Determine reference zenith angle.
                         nza = MINLOC( ABS( uv_in%view_zenith -                                    &
                                            spherical_view(j,i)%az(iaz)%zenith_blocked(ize) ),     &
                                       DIM = 1 ) - 1
!
!--                      Remove irradiance coming from obscured angles.
                         DO  nwl = 0, uv_in%n_wl-1

                            val1 = uv_in%uv_radiance_diff(0,nwl,nza,naa)
                            val2 = uv_in%uv_radiance_diff(1,nwl,nza,naa)

                            uv_rad(j,i)%ir_spectral2(nwl) =                                        &
                                                   uv_rad(j,i)%ir_spectral2(nwl)                   &
                                                 - interpolate_linear( val1, val2, fac_sza )       &
                                                 *  uv_in%d_view_zenith(nza)                       &
                                                 *  uv_in%d_view_azimuth(naa)
                         ENDDO
                      ENDDO
!
!--                   Loop over all transparently obscured zenith angles.
                      DO  ize = 1, spherical_view(j,i)%az(iaz)%n_plant_affected
!
!--                      Determine reference zenith angle.
                         nza = MINLOC( ABS( uv_in%view_zenith -                                    &
                                            spherical_view(j,i)%az(iaz)%plant_affected(ize) ),     &
                                       DIM = 1) - 1
!
!--                      Remove plant-absorbed irradiation at obscured angles.
!--                      Absorbed portion is (1 - trans).
                         trans = spherical_view(j,i)%az(iaz)%transmitted_portion(ize)
                         DO  nwl = 0, uv_in%n_wl-1

                            val1 = uv_in%uv_radiance_diff(0,nwl,nza,naa)
                            val2 = uv_in%uv_radiance_diff(1,nwl,nza,naa)

                            uv_rad(j,i)%ir_spectral2(nwl) =                                        &
                                                   uv_rad(j,i)%ir_spectral2(nwl)                   &
                                                 - ( 1.0_wp - trans )                              &
                                                 * interpolate_linear( val1, val2, fac_sza )       &
                                                 * uv_in%d_view_zenith(nza)                        &
                                                 * uv_in%d_view_azimuth(naa)

                         ENDDO
                      ENDDO
                   ENDIF

                ENDDO
             ENDIF
          ENDDO
       ENDDO

    ENDIF

    IF ( num_reflections > 0 )  THEN
!
!--    Determine the reflected portion. Note, here it is strictly assumed that the diffuse
!--    radiation is isotropic, in constrast to the obstruction. Since the reflections depend on
!--    the ratio of diffuse and direct radiation, calculate the reflected portion for each
!--    wavelength separately.
       DO  nwl = 0, uv_in%n_wl-1
          val1 = uv_in%uv_diff_irradiance(0,nwl)
          val2 = uv_in%uv_diff_irradiance(1,nwl)
          val_diff = interpolate_linear( val1, val2, fac_sza )

          val1 = uv_in%uv_dir_irradiance(0,nwl)
          val2 = uv_in%uv_dir_irradiance(1,nwl)
          val_dir = interpolate_linear( val1, val2, fac_sza )

          CALL uv_radiation_interaction_sw_only( surf_lsm%albedo_uv, surf_usm%albedo_uv,           &
                                                 val_dir, val_diff, num_reflections )
!
!--       Add respective reflected portion onto the irradiation. This is the final incoming
!--       radiation at the surface minus the direct and diffuse radiation coming from the sky.
          IF ( uv_from_irradiance )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   ind_surf = spherical_view(j,i)%ind_surfl
                   uv_rad(j,i)%ir_spectral1(nwl) = uv_rad(j,i)%ir_spectral1(nwl)                   &
                                                 + surfinsw(ind_surf)                              &
                                                 - surfinswdir(ind_surf)                           &
                                                 - surfinswdif(ind_surf)
                ENDDO
             ENDDO
          ENDIF

          IF ( uv_from_radiance )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   ind_surf = spherical_view(j,i)%ind_surfl
                   uv_rad(j,i)%ir_spectral2(nwl) = uv_rad(j,i)%ir_spectral2(nwl)                   &
                                                 + surfinsw(ind_surf)                              &
                                                 - surfinswdir(ind_surf)                           &
                                                 - surfinswdif(ind_surf)
                ENDDO
             ENDDO
          ENDIF

       ENDDO

    ENDIF

!
!-- Finally, compute the UV irradiation from integration over the wavelength,
!-- with and/or without weighting.
    IF ( uv_from_irradiance )  THEN
!
!--    Irradiance using the the from_irradiance method.
       IF ( calc_uv_ir1 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                uv_ir1(j,i) = 0.0_wp
                DO  nwl = 0, uv_in%n_wl-1
                   uv_ir1(j,i) = uv_ir1(j,i) + uv_rad(j,i)%ir_spectral1(nwl)                       &
                                             * uv_in%d_wavelength(nwl)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
!
!--    Erythemally-weighted irradiance using the the from_irradiance method.
       IF ( calc_uv_ewir1 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                uv_ewir1(j,i) = 0.0_wp
                DO  nwl = 0, uv_in%n_wl-1
                   uv_ewir1(j,i) = uv_ewir1(j,i) + uv_rad(j,i)%ir_spectral1(nwl)                   &
                                                 * uv_in%d_wavelength(nwl)                         &
                                                 * uv_in%weighting_coeff_ery(nwl)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    IF ( uv_from_radiance )  THEN
!
!--    Irradiance using the the from_radiance method.
       IF ( calc_uv_ir2 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                uv_ir2(j,i) = 0.0_wp
                DO  nwl = 0, uv_in%n_wl-1
                   uv_ir2(j,i) = uv_ir2(j,i) + uv_rad(j,i)%ir_spectral2(nwl)                       &
                                             * uv_in%d_wavelength(nwl)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
!
!--    Erythemally-weighted irradiance using the the from_radiance method.
       IF ( calc_uv_ewir2 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                uv_ewir2(j,i) = 0.0_wp
                DO  nwl = 0, uv_in%n_wl-1
                   uv_ewir2(j,i) = uv_ewir2(j,i) + uv_rad(j,i)%ir_spectral2(nwl)                   &
                                                 * uv_in%d_wavelength(nwl)                         &
                                                 * uv_in%weighting_coeff_ery(nwl)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

    ENDIF

    CALL cpu_log( log_point_s(28), 'uv model', 'stop' )

 END SUBROUTINE uv_radiation_calculate


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for the UV model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_check_data_output( var, unit, av )

    CHARACTER (LEN=*) ::  unit  !<
    CHARACTER (LEN=*) ::  var   !<

    INTEGER(iwp) ::  av  !< indicates a time-averaged output

    LOGICAL ::  trigger_error_message  !< error flag


    trigger_error_message = .FALSE.

    SELECT CASE ( TRIM( var ) )
!
!--    Erythemally-weighted irradiance in the UV spectral range.
       CASE ( 'uv_ewir1*' )
          IF ( av == 1 )  trigger_error_message = .TRUE.
          unit = 'mW m-2'
          calc_uv_ewir1 = .TRUE.

       CASE ( 'uv_ewir2*' )
          IF ( av == 1 )  trigger_error_message = .TRUE.
          unit = 'mW m-2'
          calc_uv_ewir2 = .TRUE.
!
!--    Irradiance in the UV spectral range.
       CASE ( 'uv_ir1*' )
          IF ( av == 1 )  trigger_error_message = .TRUE.
          unit = 'mW m-2'
          calc_uv_ir1 = .TRUE.

       CASE ( 'uv_ir2*' )
          IF ( av == 1 )  trigger_error_message = .TRUE.
          unit = 'mW m-2'
          calc_uv_ir2 = .TRUE.

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

    IF ( trigger_error_message )  THEN
       message_string = 'averaging is not possible for variable ' // TRIM( var )
       CALL message( 'uv_radiation_check_data_output', 'UVM0001', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE uv_radiation_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check surface data output variables from the UV radiaton model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_check_data_output_surf( trimvar, unit, av )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    ::  trimvar  !< dummy for single output variable
    CHARACTER(LEN=*), INTENT(INOUT) ::  unit     !< dummy for unit of output variable

    INTEGER(iwp), INTENT(IN) ::  av  !< id indicating average or non-average data output


    SELECT CASE ( TRIM( trimvar ) )

       CASE ( 'uv_albedo' )   ! need to be specified
          IF ( av == 1 )  THEN
             message_string = 'time averaging of quantity "' // TRIM( trimvar ) //                 &
                              '" is not provided'
             CALL message( 'uv_check_data_output_surf', 'UVM0001', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( num_reflections < 1 )  THEN
             WRITE( message_string, '(A,I3,A)' )  'num_reflections = ', num_reflections,           &
                                                  ' out of range for quantity "' //                &
                                                  TRIM( trimvar ) // '"'
             CALL message( 'uv_check_data_output_surf', 'UVM0002', 1, 2, 0, 6, 0 )
          ENDIF

          unit = '1'

       CASE DEFAULT
           unit = 'illegal'

    END SELECT

 END SUBROUTINE uv_radiation_check_data_output_surf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for the UV-radiation model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_check_parameters

!
!-- Check for correct setting of the integration method.
    IF ( INDEX( uv_integration_method, 'from_irradiance' ) == 0  .AND.                             &
         INDEX( uv_integration_method, 'from_radiance' ) == 0 )  THEN
       message_string = 'unknown uv_integration_method = "' // TRIM( uv_integration_method ) // '"'
       CALL message( 'uv_radiation_check_parameters', 'UVM0003', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Initialize control flags.
    IF ( INDEX( uv_integration_method, 'from_irradiance' ) /= 0 )  uv_from_irradiance = .TRUE.
    IF ( INDEX( uv_integration_method, 'from_radiance'   ) /= 0 )  uv_from_radiance = .TRUE.

!
!-- Check if the plant-canopy model is switched-on. If not, shading by plants is not considered.
!-- Give a warning message in this case.
    IF ( .NOT. plant_canopy )  THEN
       message_string = 'plant canopy is switched-off and thus not considered in computation ' //  &
                        'of UV-radiation'
       CALL message( 'uv_radiation_check_parameters', 'UVM0004', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Check if an UV-input file exists. Each domain can have a separate input file, depending on
!-- the chosen spherical resolution in the RTM.
    INQUIRE( FILE = TRIM( uv_in%input_file_uv ) // TRIM( coupling_char ), EXIST = input_pids_uv )

    IF ( .NOT. input_pids_uv )  THEN
       message_string = 'no UV-input file found'
       CALL message( 'uv_radiation_check_parameters', 'UVM0005', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check if the number of spherical angles in the input file matches with the one considered for
!-- the sky-factor calculation. If not, give a warning message that this can lead to uncertainties.
!-- Therefore, call the NetCDF input routine for the UV model. Therein, first initializations of
!-- the dimension are carried out and checks if all variables and dimensions are present are
!-- performed. Variable input is skipped.
    CALL uv_radiation_rff_input( skip_variable_input = .TRUE. )

    IF ( uv_from_radiance )  THEN
       IF ( raytrace_discrete_azims /= uv_in%n_aa )  THEN
          WRITE( message_string, '(A,I3,A,I3)' )  'mismatch between raytrace_discrete_azims = ',   &
                                                  raytrace_discrete_azims,                         &
                                                  '&and value found in in UV input file = ',       &
                                                  uv_in%n_aa
          CALL message( 'uv_radiation_check_parameters', 'UVM0006', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( ( raytrace_discrete_elevs / 2 ) /= uv_in%n_za )  THEN
          WRITE( message_string, '(A,I3,A,I3)' )  'mismatch between raytrace_discrete_elevs/2 = ', &
                                                  raytrace_discrete_elevs / 2,                     &
                                                  '&and value found in in UV input file = ',       &
                                                  uv_in%n_za
          CALL message( 'uv_radiation_check_parameters', 'UVM0006', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- If reflections are switched-on, direct and diffuse irradiances are required, even if only
!-- the radiance method is used.
    IF ( num_reflections > 0  .AND.  .NOT. ( ALLOCATED( uv_in%uv_diff_irradiance )  .AND.          &
                                             ALLOCATED( uv_in%uv_dir_irradiance ) ) )              &
    THEN
       message_string = 'num_reflections > 0 but no variables "' // TRIM( uv_in%char_ir_diff ) //  &
                        '" and "' // TRIM( uv_in%char_ir_dir ) // '"&found in UV input file'
       CALL message( 'uv_radiation_check_parameters', 'UVM0007', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE uv_radiation_check_parameters


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do,       &
                                         nzt_do )

    CHARACTER (LEN=*), INTENT(OUT) ::  grid      !< grid type (always "zu1" for uv-radiation)
    CHARACTER (LEN=*)              ::  variable  !< treated variable

    INTEGER(iwp) ::  av      !< flag indicating instantaneous or averaged data output
    INTEGER(iwp) ::  i       !< grid index x-direction
    INTEGER(iwp) ::  j       !< grid index y-direction
    INTEGER(iwp) ::  nzb_do  !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do  !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL ::  found  !< flag indicating if variable is found
    LOGICAL ::  two_d  !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< data output array


    found = .TRUE.
    two_d = .FALSE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'uv_ewir1*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_ewir1(j,i)
                ENDDO
             ENDDO
          ENDIF
          grid = 'zu1'
          two_d = .TRUE.

       CASE ( 'uv_ewir2*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_ewir2(j,i)
                ENDDO
             ENDDO
          ENDIF
          grid = 'zu1'
          two_d = .TRUE.

       CASE ( 'uv_ir1*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_ir1(j,i)
                ENDDO
             ENDDO
          ENDIF
          grid = 'zu1'
          two_d = .TRUE.

       CASE ( 'uv_ir2*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_ir2(j,i)
                ENDDO
             ENDDO
          ENDIF
          grid = 'zu1'
          two_d = .TRUE.

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE uv_radiation_data_output_2d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> UV-radiation surface output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_data_output_surf( av, trimvar, found )

    CHARACTER (LEN=*), INTENT(IN) ::  trimvar  !< variable name

    INTEGER(iwp), INTENT(IN) ::  av  !< flag for (non-)average output

    INTEGER(iwp) ::  i       !< running index x-direction
    INTEGER(iwp) ::  j       !< running index y-direction
    INTEGER(iwp) ::  m       !< running index for surface elements
    INTEGER(iwp) ::  n_surf  !< running index for output surface elements


    LOGICAL, INTENT(INOUT) ::  found  !< flag if output variable is found


    found = .TRUE.

    n_surf = 0

    SELECT CASE ( TRIM( trimvar ) )

       CASE ( 'uv_albedo' )
          DO  i = nxl, nxr
             DO  j = nys, nyn

                DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                   n_surf = n_surf + 1
                   surf_out%var_out(n_surf) = SUM( surf_def%frac(m,:) * surf_def%albedo_uv(m,:) )
                ENDDO
                DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                   n_surf = n_surf + 1
                   surf_out%var_out(n_surf) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo_uv(m,:) )
                ENDDO
                DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                   n_surf = n_surf + 1
                   surf_out%var_out(n_surf) = SUM( surf_usm%frac(m,:) * surf_usm%albedo_uv(m,:) )
                ENDDO

             ENDDO
          ENDDO

       CASE DEFAULT
          found = .FALSE.

    END SELECT
!
!-- Silence the compiler warning about unused parameters.
    IF ( av == 0  .OR.  i == 0  .OR.  j == 0 )  CONTINUE

 END SUBROUTINE uv_radiation_data_output_surf


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    CHARACTER (LEN=*), INTENT(IN)  ::  var     !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x  !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y  !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z  !<

    LOGICAL, INTENT(OUT) ::  found  !<


    found  = .TRUE.

    SELECT CASE ( TRIM( var ) )

       CASE ( 'uv_ewir1*_xy', 'uv_ewir2*_xy', 'uv_ir1*_xy', 'uv_ir2*_xy' )

          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu1'

       CASE DEFAULT

          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE uv_radiation_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for the UV-radiation model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_header ( io )

    INTEGER(iwp),  INTENT(IN) ::  io  !< unit of the output file

    WRITE( io, 1 )
    WRITE( io, 2 )
    WRITE( io, 3 ) num_reflections

1   FORMAT (/ /' UV-radiation model switched-on.', A )
2   FORMAT (   ' -------------------------------', A )
3   FORMAT (   ' number of reflection steps: ', I1 )

 END SUBROUTINE uv_radiation_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the UV-radiation model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_init

    INTEGER(iwp) ::  i      !< running index in x-direction
    INTEGER(iwp) ::  iaz    !<
    INTEGER(iwp) ::  isurf  !< running index over sky-view factors
    INTEGER(iwp) ::  j      !< running index in y-direction
    INTEGER(iwp) ::  m      !< running index over type surface elements
    INTEGER(iwp) ::  n      !< running index over wavelength
    INTEGER(iwp) ::  n1     !< running index over view azimuth angles
    INTEGER(iwp) ::  n2     !< running index over view zenith angles

    REAL(wp) ::  check_sum  !< dummy used to sum up discrete spherical angles
    REAL(wp) ::  d_angle    !< angle increment
    REAL(wp) ::  theta_l    !< lower integration bound for view zenith angle
    REAL(wp) ::  theta_u    !< upper integration bound for view zenith angle


    IF ( debug_output )  CALL debug_message( 'uv_radiation_init', 'start' )
!
!-- The UV-radiation model requires the radiation model and the RTM. The control flag
!-- radiation_interactions, which controls the RTM calls, hasn't been set at this point in time.
!-- Hence, mimic the conditions that are equivalent to radiation_interactions. Please be aware
!-- that the UV-model can also run without any information from the radiation model, though this
!-- means that no shading is considered at all and the driver input is simply integrated and output.
    IF ( .NOT.  radiation  .OR.                                                                    &
         .NOT. ( radiation_interactions_on  .AND.  vertical_surfaces_exist ) )   THEN
       message_string = 'UV-radiation model runs without information from radiative transfer model'
       CALL message( 'uv_radiation_check_parameters', 'UVM0008', 0, 1, 0, 6, 0 )
    ENDIF

!
!-- Allocate output arrays.
    IF ( calc_uv_ewir1 )  THEN
       ALLOCATE( uv_ewir1(nys:nyn,nxl:nxr) )
       uv_ewir1 = 0.0_wp
    ENDIF
    IF ( calc_uv_ewir2 )  THEN
       ALLOCATE( uv_ewir2(nys:nyn,nxl:nxr) )
       uv_ewir2 = 0.0_wp
    ENDIF
    IF ( calc_uv_ir1 )  THEN
       ALLOCATE( uv_ir1(nys:nyn,nxl:nxr) )
       uv_ir1 = 0.0_wp
    ENDIF
    IF ( calc_uv_ir2 )  THEN
       ALLOCATE( uv_ir2(nys:nyn,nxl:nxr) )
       uv_ir2 = 0.0_wp
    ENDIF
!
!-- Initialize multiple reflections.
    IF ( num_reflections > 0 )  THEN
!
!--    Allocate surface arrays for UV-albedo (only required if multiple reflections are
!--    switched-on).
       IF ( surf_lsm%ns > 0 )  ALLOCATE( surf_lsm%albedo_uv(1:surf_lsm%ns,0:2) )
       IF ( surf_usm%ns > 0 )  ALLOCATE( surf_usm%albedo_uv(1:surf_usm%ns,0:2) )
!
!--    Initialize broadband UV albedo by bulk parameters (list need to be checked).
       DO  m = 1, surf_lsm%ns

          surf_lsm%albedo_uv(m,:) = 0.1_wp

          IF ( surf_lsm%albedo_type(m,ind_veg_wall) /= 0 )                                         &
             surf_lsm%albedo_uv(m,ind_veg_wall) =                                                  &
                                                albedo_pars(3,surf_lsm%albedo_type(m,ind_veg_wall))

          IF ( surf_lsm%albedo_type(m,ind_pav_green) /= 0 )                                        &
             surf_lsm%albedo_uv(m,ind_pav_green) =                                                 &
                                                albedo_pars(3,surf_lsm%albedo_type(m,ind_pav_green))

          IF ( surf_lsm%albedo_type(m,ind_wat_win) /= 0 )                                          &
             surf_lsm%albedo_uv(m,ind_wat_win) =                                                   &
                                                albedo_pars(3,surf_lsm%albedo_type(m,ind_wat_win))

       ENDDO

       DO  m = 1, surf_usm%ns

          surf_usm%albedo_uv(m,:) = 0.1_wp

          IF ( surf_usm%albedo_type(m,ind_veg_wall) /= 0 )                                         &
             surf_usm%albedo_uv(m,ind_veg_wall) =                                                  &
                                                albedo_pars(3,surf_usm%albedo_type(m,ind_veg_wall))

          IF ( surf_usm%albedo_type(m,ind_pav_green) /= 0 )                                        &
             surf_usm%albedo_uv(m,ind_pav_green) =                                                 &
                                                albedo_pars(3,surf_usm%albedo_type(m,ind_pav_green))

          IF ( surf_usm%albedo_type(m,ind_wat_win) /= 0 )                                          &
             surf_usm%albedo_uv(m,ind_wat_win) =                                                   &
                                                albedo_pars(3,surf_usm%albedo_type(m,ind_wat_win))

       ENDDO

    ENDIF
!
!-- Pre-compute the weighting coefficients for erythema weighting.
!-- Weighting function is computed according to ISO/CIE:17166 (2019).
    ALLOCATE( uv_in%weighting_coeff_ery(0:uv_in%n_wl-1) )

    DO  n = 0, uv_in%n_wl-1
       IF ( uv_in%wavelength(n) < 298.0_wp )  THEN
          uv_in%weighting_coeff_ery(n) = 1.0_wp
       ELSEIF ( uv_in%wavelength(n) >= 298.0_wp  .AND.  uv_in%wavelength(n) < 328.0_wp )  THEN
          uv_in%weighting_coeff_ery(n) = 10**( 0.094_wp * ( 298.0_wp - uv_in%wavelength(n) ) )
       ELSEIF ( uv_in%wavelength(n) >= 328.0_wp  .AND.  uv_in%wavelength(n) < 400.0_wp )  THEN
          uv_in%weighting_coeff_ery(n) = 10**( 0.015_wp * ( 139.0_wp - uv_in%wavelength(n) ) )
       ELSE
          uv_in%weighting_coeff_ery(n) = 0.0_wp
       ENDIF
    ENDDO
!
!-- Initialize the data structure for the actual computed radiation.
    ALLOCATE( uv_rad(nys:nyn,nxl:nxr) )
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Set initial values to zero.
          uv_rad(j,i)%uv_ery = 0.0_wp
          uv_rad(j,i)%uvi    = 0.0_wp
          uv_rad(j,i)%svf    = -HUGE( 1.0_wp )
          uv_rad(j,i)%svf_p  = -HUGE( 1.0_wp )
!
!--       Allocate array for spectral irradiation. Consider two arrays, one for
!--       calculation from irradiances, from for calculation from radiances.
          IF ( uv_from_irradiance )  ALLOCATE( uv_rad(j,i)%ir_spectral1(0:uv_in%n_wl-1) )
          IF ( uv_from_radiance   )  ALLOCATE( uv_rad(j,i)%ir_spectral2(0:uv_in%n_wl-1) )
       ENDDO
    ENDDO

!
!-- Store sky-view factors for horizontally upward facing surfaces. This is only required
!-- if UV radiation shall be computed from irradiance too.
    IF ( uv_from_irradiance )  THEN
       DO  isurf = 1, nsurfl
          i = surfl(ix,isurf)
          j = surfl(iy,isurf)
!
!--       If subsricpt id ( = 1 ) of dimension 0 equals 0, this is a horizontally upward facing
!--       surface.
          IF ( surfl(id,isurf) == 0  .AND.  spherical_view(j,i)%k_surf == surfl(iz,isurf) )  THEN
             uv_rad(j,i)%svf   = skyvf(isurf)
             uv_rad(j,i)%svf_p = skyvft(isurf)
          ENDIF
       ENDDO
    ENDIF
!
!-- In case of radiance integration, convert view factors (inferred from RTM raytracing) from
!-- radiants to degrees and bring them to the coordinate system used in the external file.
    IF( uv_from_radiance )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn

             IF ( ALLOCATED( spherical_view(j,i)%az ) )  THEN
                spherical_view(j,i)%treat_az = .TRUE.
                DO  iaz = 1, spherical_view(j,i)%n_az
                   spherical_view(j,i)%az(iaz)%val_az = spherical_view(j,i)%az(iaz)%val_az *       &
                                                        radiants_to_degrees

                   spherical_view(j,i)%az(iaz)%val_az =                                            &
                                 MODULO( spherical_view(j,i)%az(iaz)%val_az - 180.0_wp, 360.0_wp )

                   IF ( ALLOCATED( spherical_view(j,i)%az(iaz)%zenith_blocked )  .OR.              &
                        ALLOCATED( spherical_view(j,i)%az(iaz)%plant_affected ) )  THEN

                      IF ( spherical_view(j,i)%az(iaz)%n_blocked > 0 )  THEN
                         spherical_view(j,i)%az(iaz)%zenith_blocked =                              &
                                                      spherical_view(j,i)%az(iaz)%zenith_blocked * &
                                                      radiants_to_degrees
                      ENDIF

                      IF ( spherical_view(j,i)%az(iaz)%n_plant_affected > 0 )  THEN
                         spherical_view(j,i)%az(iaz)%plant_affected =                              &
                                                      spherical_view(j,i)%az(iaz)%plant_affected * &
                                                      radiants_to_degrees
                      ENDIF

                      spherical_view(j,i)%az(iaz)%treat_ze = .TRUE.

                   ELSE

                      spherical_view(j,i)%az(iaz)%treat_ze = .FALSE.

                   ENDIF
                ENDDO
             ELSE
                spherical_view(j,i)%treat_az = .FALSE.
             ENDIF

          ENDDO
       ENDDO
    ENDIF
!
!-- Determine integration increments for the wavelengths.
    ALLOCATE( uv_in%d_wavelength(0:uv_in%n_wl-1) )
    DO  n = 1, uv_in%n_wl-1
       uv_in%d_wavelength(n) = uv_in%wavelength(n) - uv_in%wavelength(n-1)
    ENDDO
    uv_in%d_wavelength(0) = uv_in%d_wavelength(1)
!
!-- Determine integration increments for the view azimuth and view zenith angles. Only required,
!-- if UV exposure shall be computed from radiances. At this step, already convert to
!-- radiants. d_view_azimuth times d_view_zenith gives the area increment for the given
!-- spherical angle. Furthermore, shift the azimuth and zenith angles by half a discrete
!-- distances. This is necessary as uvspec output is left-bounded.
    IF( uv_from_radiance )  THEN

       ALLOCATE( uv_in%view_azimuth_rot(0:uv_in%n_aa-1) )

       ALLOCATE( uv_in%d_view_azimuth(0:uv_in%n_aa-1) )
       DO  n1 = 1, uv_in%n_aa-1
          uv_in%d_view_azimuth(n1) = uv_in%view_azimuth(n1) - uv_in%view_azimuth(n1-1)
       ENDDO
       uv_in%d_view_azimuth(0) = uv_in%d_view_azimuth(1)
!
!--    Convert into radiants.
       uv_in%d_view_azimuth = uv_in%d_view_azimuth * degrees_to_radiants
!
!--    Zenith angles are converted into -cos(theta), with theta being the view zenith
!--    angle. Index 0 is for horizon, index n_za is for horizon.
       ALLOCATE( uv_in%d_view_zenith(0:uv_in%n_za-1) )
       DO  n2 = 0, uv_in%n_za-2
          theta_u = -COS( uv_in%view_zenith(n2+1) * degrees_to_radiants )
          theta_l = -COS( uv_in%view_zenith(n2)   * degrees_to_radiants )
          uv_in%d_view_zenith(n2) = theta_u - theta_l
       ENDDO
!
!--    Special treatment for horizon. This is necessary as the horizon angle is not
!--    part of the uvspec output.
       n2 = uv_in%n_za-1
       theta_u = -COS( 90.0_wp               * degrees_to_radiants )
       theta_l = -COS( uv_in%view_zenith(n2) * degrees_to_radiants )
       uv_in%d_view_zenith(n2) = theta_u - theta_l
!
!--    Finally, check if the given view angles have plausible values.
!--    Integration of d_view_zenith * d_view_azimuth over all elements should yield 2 x PI,
!--    i.e. the surface area of an upper unit half-sphere. Here, compute the check sum. If this
!--    deviates significantly from 2 x PI, give an error message.
       check_sum = 0.0_wp
       DO  n1 = 0, uv_in%n_aa-1
          DO  n2 = 0, uv_in%n_za-1
             check_sum = check_sum + uv_in%d_view_azimuth(n1) * uv_in%d_view_zenith(n2)
          ENDDO
       ENDDO

       IF ( ABS( check_sum - 2.0_wp * pi ) > 10E-1_wp )  THEN
          message_string = 'view_azimuth and view_zenith angles provided via UV input file ' //    &
                           '&do not cover an upper half-sphere'
          CALL message( 'uv_radiation_init', 'UVM0009', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Since view angles from uvspec output are left-sided but RTM view angles are centered for
!--    their representative spherical angle interval, shift view angles by half a spherical
!--    resolution. At this point, the view angles are only used to identify indices.
       d_angle = uv_in%view_zenith(1) - uv_in%view_zenith(0)
       DO  n2 = 0, uv_in%n_za-1
          uv_in%view_zenith(n2) = uv_in%view_zenith(n2) + 0.5_wp * d_angle
       ENDDO

       d_angle = uv_in%view_azimuth(1) - uv_in%view_azimuth(0)
       DO  n1 = 0, uv_in%n_aa-1
          uv_in%view_azimuth(n1) = uv_in%view_azimuth(n1) + 0.5_wp * d_angle
       ENDDO
!
!--    Allocate array for diffuse part of radiance.
       ALLOCATE( uv_in%uv_radiance_diff(0:1,0:uv_in%n_wl-1,0:uv_in%n_za-1,0:uv_in%n_aa-1) )

    ENDIF
!
!-- Initialize read index.
    uv_in%ind_sz_l = 0
    uv_in%ind_sz_u = 0

    IF ( debug_output )  CALL debug_message( 'uv_radiation_init', 'end' )

 END SUBROUTINE uv_radiation_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read namelist &uv_radiation_parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /uv_radiation_parameters/  num_reflections,                                           &
                                        switch_off_module,                                         &
                                        uv_integration_method


!
!-- Move to the beginning of the namelist file and try to find and read the user-defined namelist
!-- uv_radiation_parameters.
    REWIND( 11 )
    READ( 11, uv_radiation_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    uv_radiation_parameters namelist was found and read correctly. Set flag that indicates that
!--    the uv-radiation model is switched on.
       IF ( .NOT. switch_off_module )  uv_radiation = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    uv_radiation_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'uv_radiation_parameters', line )

    ENDIF

 END SUBROUTINE uv_radiation_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read NetCDF input data and initially allocate memory.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_rff_input( skip_variable_input )

#if defined ( __netcdf )
    CHARACTER(LEN=500) ::  variable_not_found = ''  !< string to gather information of missing dimensions and variables

    INTEGER(iwp) ::  pids_id_uv  !< NetCDF ID
#endif
    LOGICAL ::            skip                   !< flag to skip reading of variables
    LOGICAL, OPTIONAL ::  skip_variable_input    !< passed flag to skip reading of variables
#if defined ( __netcdf )
    LOGICAL ::            trigger_error_message  !< control flag to check whether all required input data is available
#endif

!
!-- Determine the skip flag. By default, reading of variables is not skipped.
    skip = .FALSE.
    IF ( PRESENT( skip_variable_input ) )  skip = skip_variable_input

    IF ( input_pids_uv )  THEN
#if defined ( __netcdf )
!
!--    Open file in read-only mode.
       CALL open_read_file( TRIM( uv_in%input_file_uv ) // TRIM( coupling_char ), pids_id_uv )
!
!--    If not yet initialized, inquire variable names, read dimensions, etc..
       IF ( .NOT. uv_in%file_read_initialized )  THEN
          trigger_error_message = .FALSE.
          variable_not_found = ''
!
!--       At first, inquire all variable names.
!--       This will be used to check whether an input variable exists or not.
          CALL inquire_num_variables( pids_id_uv, num_var_pids )
!
!--       Allocate memory to store variable names and read them.
          ALLOCATE( vars_pids(1:num_var_pids) )
          CALL inquire_variable_names( pids_id_uv, vars_pids )
!
!--       Read azimuth dimension.
          IF ( check_existence( vars_pids, uv_in%char_vaa ) )  THEN
             CALL get_dimension_length( pids_id_uv, uv_in%n_aa, uv_in%char_vaa )
             ALLOCATE( uv_in%view_azimuth(0:uv_in%n_aa-1) )
             CALL get_variable( pids_id_uv, uv_in%char_vaa, uv_in%view_azimuth )
          ELSE
             trigger_error_message = .TRUE.
             variable_not_found = TRIM( uv_in%char_vaa ) // ', '
          ENDIF
!
!--       Read zenith dimension.
          IF ( check_existence( vars_pids, uv_in%char_vza ) )  THEN
             CALL get_dimension_length( pids_id_uv, uv_in%n_za, uv_in%char_vza )
             ALLOCATE( uv_in%view_zenith(0:uv_in%n_za-1) )
             CALL get_variable( pids_id_uv, uv_in%char_vza, uv_in%view_zenith )
          ELSE
             trigger_error_message = .TRUE.
             variable_not_found = TRIM( uv_in%char_vza ) // ', '
          ENDIF
!
!--       Read sun zenith dimension.
          IF ( check_existence( vars_pids, uv_in%char_sza ) )  THEN
             CALL get_dimension_length( pids_id_uv, uv_in%n_sza, uv_in%char_sza )
             ALLOCATE( uv_in%sun_zenith(0:uv_in%n_sza-1) )
             CALL get_variable( pids_id_uv, uv_in%char_sza, uv_in%sun_zenith )
          ELSE
             trigger_error_message = .TRUE.
             variable_not_found = TRIM( uv_in%char_sza ) // ', '
          ENDIF
!
!--       Read wavelength dimension.
          IF ( check_existence( vars_pids, uv_in%char_wl ) )  THEN
             CALL get_dimension_length( pids_id_uv, uv_in%n_wl, uv_in%char_wl )
             ALLOCATE( uv_in%wavelength(0:uv_in%n_wl-1) )
             CALL get_variable( pids_id_uv, uv_in%char_wl, uv_in%wavelength )
          ELSE
             trigger_error_message = .TRUE.
             variable_not_found = TRIM( uv_in%char_wl ) // ', '
          ENDIF
!
!--       Check if any dimension is missing.
          IF ( trigger_error_message )  THEN
             message_string = 'dimension(s): &' // TRIM( variable_not_found ) //                   &
                              '& not found in UV input file'
             CALL message( 'uv_radiation_rff_input', 'UVM0010', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Now check for missing variables. Variables will be read later.
          variable_not_found = ''
!
!--       Direct irradiance.
          IF ( uv_from_irradiance  .OR.  num_reflections > 0 )  THEN
             IF ( .NOT. check_existence( vars_pids, uv_in%char_ir_dir ) )  THEN
                trigger_error_message = .TRUE.
                variable_not_found = TRIM( uv_in%char_ir_dir ) // ', '
             ENDIF
          ENDIF
!
!--       Diffuse irradiance.
          IF ( uv_from_irradiance  .OR.  num_reflections > 0 )  THEN
             IF ( .NOT. check_existence( vars_pids, uv_in%char_ir_diff ) )  THEN
                trigger_error_message = .TRUE.
                variable_not_found = TRIM( uv_in%char_ir_diff ) // ', '
             ENDIF
          ENDIF
!
!--       Radiance.
          IF ( uv_from_radiance )  THEN
             IF ( .NOT. check_existence( vars_pids, uv_in%char_rad ) )  THEN
                trigger_error_message = .TRUE.
                variable_not_found = TRIM( uv_in%char_rad ) // ', '
             ENDIF
          ENDIF

          IF ( trigger_error_message )  THEN
             message_string = 'variable(s): &' // TRIM( variable_not_found ) //                    &
                              '& not found in UV input file'
             CALL message( 'uv_radiation_rff_input', 'UVM0010', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Allocate memory for the input arrays. Only memory for the two closest sun-zeniths
!--       is allocated.
          IF ( uv_from_irradiance )  ALLOCATE( uv_in%uv_diff_irradiance(0:1,0:uv_in%n_wl-1) )
          IF ( uv_from_irradiance )  ALLOCATE( uv_in%uv_dir_irradiance(0:1,0:uv_in%n_wl-1) )

          IF ( uv_from_radiance )                                                                  &
             ALLOCATE( uv_in%uv_radiance(0:1,0:uv_in%n_wl-1,0:uv_in%n_za-1,0:uv_in%n_aa-1) )
!
!--       Read _FillValue attribute.
          IF ( uv_from_irradiance  .OR.  num_reflections > 0 )  THEN
             CALL get_attribute( pids_id_uv, char_fill, uv_in%fill_diff, .FALSE.,                  &
                                 uv_in%char_ir_diff, .FALSE. )
             CALL get_attribute( pids_id_uv, char_fill, uv_in%fill_dir, .FALSE.,                   &
                                 uv_in%char_ir_dir, .FALSE. )
          ENDIF
          IF ( uv_from_radiance )  THEN
             CALL get_attribute( pids_id_uv, char_fill, uv_in%fill_rad, .FALSE., uv_in%char_rad,   &
                                 .FALSE. )
          ENDIF

          uv_in%file_read_initialized = .TRUE.

          DEALLOCATE( vars_pids )
       ENDIF

!
!--    Read variables if not skipped. Reading is only skipped when the subroutine is called
!--    from check_parameters.
       IF ( .NOT. skip )  THEN

          IF ( uv_from_irradiance  .OR.  num_reflections > 0 )  THEN
!
!--          Read spectral diffuse irradiance.
             CALL get_variable( pids_id_uv, uv_in%char_ir_diff, uv_in%uv_diff_irradiance,          &
                                0, uv_in%n_wl-1, uv_in%ind_sz_l, uv_in%ind_sz_u, nbgp=0 )
!
!--          Read spectral diffuse irradiance.
             CALL get_variable( pids_id_uv, uv_in%char_ir_dir, uv_in%uv_dir_irradiance,            &
                                0, uv_in%n_wl-1, uv_in%ind_sz_l, uv_in%ind_sz_u, nbgp=0 )
          ENDIF

          IF ( uv_from_radiance )  THEN
!
!--          Read spectral radiance.
             CALL get_variable( pids_id_uv, uv_in%char_rad, uv_in%uv_radiance,                     &
                                0, uv_in%n_aa-1, 0, uv_in%n_za-1, 0, uv_in%n_wl-1,                 &
                                uv_in%ind_sz_l, uv_in%ind_sz_u, nbgp=0 )
          ENDIF

       ENDIF
!
!--    Close input file.
       CALL close_input_file( pids_id_uv )
#endif
    ENDIF

 END SUBROUTINE uv_radiation_rff_input


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine to simulate multiple reflections of shortwave radiation. See description of
!> radiation_interaction for more details. This routine is actually a slimmed-down version of
!> radiation_interaction in radiation_model_mod.f90
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uv_radiation_interaction_sw_only( alb_l, alb_u, dir_in, diff_in, num_reflections )

     IMPLICIT NONE

     INTEGER(iwp) ::  d                !< running index over surface orientations
     INTEGER(iwp) ::  i                !< grid index in x-direction
#if defined( __parallel )
     INTEGER(iwp) ::  surf_start_id    !< id of first surface in current processor
#endif
     INTEGER(iwp) ::  isd              !< corresponding index in pre-calculated direct solar shading array
     INTEGER(iwp) ::  isurf            !< index in RTM-related surface array for the target surface
     INTEGER(iwp) ::  isurfsrc         !< index in RTM-related surface array for the source surface
     INTEGER(iwp) ::  isvf             !< running index in RTM-related surface array
     INTEGER(iwp) ::  j                !< grid index in y-direction
     INTEGER(iwp) ::  m                !< running index surface elements
     INTEGER(iwp) ::  mm               !< index for aggregated surface properties
     INTEGER(iwp) ::  num_reflections  !< number of considered reflections
     INTEGER(iwp) ::  refstep          !< running index for reflections

     REAL(wp) ::  cos_rot                        !< cosine of rotation_angle
     REAL(wp) ::  diff_in                        !< diffuse sky radiation
     REAL(wp) ::  dir_in                         !< direct sky radiation
     REAL(wp) ::  min_stable_coszen = 0.0262_wp  !< 1.5 deg above horizon, eliminates most of circumsolar
     REAL(wp) ::  sin_rot                        !< sine of rotation_angle
     REAL(wp) ::  solar_azim                     !< solar azimuth in rotated model coordinates
     REAL(wp) ::  sun_direct_factor              !< factor for direct normal radiation from direct horizontal

     REAL(wp), DIMENSION(3) ::  sunorig  !< grid rotated solar direction unit vector (zyx)

     REAL(wp), DIMENSION(0:nsurf_type) ::  costheta  !< direct irradiance factor of solar angle

     REAL(wp), DIMENSION(3,3) ::  mrot  !< grid rotation matrix (zyx)

     REAL(wp), DIMENSION(3,0:nsurf_type) ::  vnorm  !< face direction normal vectors (zyx)

     REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  alb_l  !< land-surface albedo used from shortwave treatment
     REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  alb_u  !< urban-surface albedo used from shortwave treatment


     IF ( debug_output_timestep )  THEN
        WRITE( debug_string, * ) 'uv_radiation_interaction_sw_only', time_since_reference_point
        CALL debug_message( debug_string, 'start' )
     ENDIF

     sun_direction = .TRUE.
     CALL get_date_time( time_since_reference_point, day_of_year=day_of_year,                      &
                         second_of_day = second_of_day )
!
!--  Following data is also required for diffuse radiation.
     CALL calc_zenith( day_of_year, second_of_day )

!
!--  Prepare rotated normal vectors and irradiance factor.
     sin_rot = SIN( rotation_angle * pi / 180.0_wp )
     cos_rot = COS( rotation_angle * pi / 180.0_wp )
     vnorm(1,:) = kdir(:)
     vnorm(2,:) = jdir(:)
     vnorm(3,:) = idir(:)

     mrot(1,:) = (/ 1.0_wp,  0.0_wp,   0.0_wp /)
     mrot(2,:) = (/ 0.0_wp,  cos_rot, sin_rot /)
     mrot(3,:) = (/ 0.0_wp, -sin_rot, cos_rot /)
     sunorig = (/ cos_zenith, sun_dir_lat, sun_dir_lon /)
     sunorig = MATMUL( mrot, sunorig )

!
!--  Direct irradiance factor of solar angle, avoid negative value to prevent negative direct SW
!--  values.
     DO  d = 0, nsurf_type
        costheta(d) = MAX( DOT_PRODUCT( sunorig, vnorm(:,d) ), 0.0_wp )
     ENDDO

!
!--  Initialize relavant surface flux arrays and radiation energy sum.
     surfinswdir  = 0.0_wp
     surfins      = 0.0_wp
     surfoutsl(:) = 0.0_wp
!
!--  Set up thermal radiation from surfaces
     mm = 1
!--  Following code depends on the order of the execution. Do not parallelize by OpenMP!
     DO  i = nxl, nxr
        DO  j = nys, nyn
!
!--        Urban-type surfaces
           DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
              albedo_surf(mm) = SUM( surf_usm%frac(m,:) * alb_u(m,:) )
              mm = mm + 1
           ENDDO
!
!--        Land surfaces
           DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
              albedo_surf(mm) = SUM( surf_lsm%frac(m,:) * alb_l(m,:) )
              mm = mm + 1
           ENDDO
        ENDDO
     ENDDO
!
!--  Direct radiation
     IF ( cos_zenith > 0 )  THEN
!
!--     To avoid numerical instability near horizon depending on what direct radiation is used
!--     (slightly different zenith angle, considering circumsolar etc.), we use a minimum value for
!--     cos_zenith.
        sun_direct_factor = 1.0_wp / MAX( min_stable_coszen, cos_zenith )
!
!--     Identify solar direction vector (discretized number) (1).
        solar_azim = ATAN2( sun_dir_lon, sun_dir_lat ) * ( 180.0_wp / pi ) - rotation_angle
        j = FLOOR( ACOS( cos_zenith ) / pi * REAL( raytrace_discrete_elevs, KIND = wp ) )
        i = MODULO( NINT( solar_azim / 360.0_wp * REAL( raytrace_discrete_azims, KIND = wp )       &
                          - 0.5_wp, iwp ), raytrace_discrete_azims )
        isd = dsidir_rev(j,i)
!
!-- TODO: check if isd = -1 to report that this solar position is not precalculated
        DO  isurf = 1, nsurfl
           j = surfl(iy,isurf)
           i = surfl(ix,isurf)
           d = surfl(id,isurf)
           surfinswdir(isurf) = dir_in * costheta(surfl(id,isurf)) *                               &
                                dsitrans(isurf,isd) * sun_direct_factor
           surfinswdif(isurf) = diff_in * skyvft(isurf)
        ENDDO

     ENDIF

     surfins  = surfinswdir + surfinswdif
     surfinsw = surfins
!
!--  Next passes of radiation interactions: Radiation reflections.
     DO  refstep = 1, num_reflections

        surfoutsl = albedo_surf * surfins

#if defined( __parallel )
!
!--     Sending out flux, surfoutsl.
        surf_start_id = surfstart(myid)
        DO  i = 1, nsend_radx
           radx_send(i) = surfoutsl(isurf_send_radx(i) - surf_start_id)
        ENDDO

        CALL rtm_alltoallv( radx_send, disp_send_radx, surfouts_recv, disp_recv_radx )

#endif
!
!--     Reset for the input from next reflective pass.
        surfins = 0.0_wp
!
!--     Reflected radiation.
        DO  isvf = 1, nsvfl
           isurf = svfsurf(1,isvf)
           isurfsrc = svfsurf(2,isvf)
#if defined( __parallel )
           surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfouts_recv(isurfsrc)
#else
           surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutsl(isurfsrc)
#endif
        ENDDO

        surfinsw  = surfinsw  + surfins

     ENDDO

     IF ( debug_output_timestep )  CALL debug_message( 'uv_radiation_interaction_sw_only', 'end' )

 END SUBROUTINE uv_radiation_interaction_sw_only

 END MODULE uv_radiation_model_mod
