!> @file slurb_mod.f90
!--------------------------------------------------------------------------------------------------!
! This program is free software; you can redistribute it and/or modify it under the terms of the
! GNU General Public License as published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this program.
! If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2022-2024 University of Helsinki
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Sasu Karttunen <sasu.karttunen@helsinki.fi>
!
! Description:
! ------------
!> Single-layer urban surface model (SLUrb) to compute online fluxes for non-resolved urban
!> surfaces. The model is a resistance-based urban canopy model based on the physical formulation of
!> the Town Energy Balance (TEB) model of SURFEX (see V. Masson, 2000 and subsequent TEB papers),
!> with numerics similar to those of PALM's USM and LSM models.
!--------------------------------------------------------------------------------------------------!
 MODULE slurb_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  exner,                                                                              &
               dzw,                                                                                &
               dzw,                                                                                &
               ddzw,                                                                               &
               drho_air,                                                                           &
               drho_air_zw,                                                                        &
               d_exner,                                                                            &
               heatflux_input_conversion,                                                          &
               heatflux_output_conversion,                                                         &
               hyrho,                                                                              &
               prr,                                                                                &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               rho_air,                                                                            &
               rho_air,                                                                            &
               rho_air_zw,                                                                         &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               waterflux_input_conversion,                                                         &
               waterflux_output_conversion,                                                        &
               zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  c_p,                                                                                &
               g,                                                                                  &
               l_v,                                                                                &
               lv_d_cp,                                                                            &
               kappa,                                                                              &
               magnus,                                                                             &
               pi,                                                                                 &
               rd_d_rv,                                                                            &
               rho_l,                                                                              &
               sigma_sb

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d,                                                                   &
               cloud_droplets,                                                                     &
               cyclic_fill_initialization,                                                         &
               data_output_raw,                                                                    &
               debug_output,                                                                       &
               debug_output_timestep,                                                              &
               debug_string,                                                                       &
               dt_3d,                                                                              &
               flux_output_mode,                                                                   &
               humidity,                                                                           &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               latitude,                                                                           &
               longitude,                                                                          &
               loop_optimization,                                                                  &
               message_string,                                                                     &
               output_fill_value,                                                                  &
               pt_surface,                                                                         &
               read_spinup_data,                                                                   &
               restart_data_format_output,                                                         &
               rho_cp,                                                                             &
               rho_surface,                                                                        &
               rotation_angle,                                                                     &
               time_since_reference_point,                                                         &
               slurb,                                                                              &
               spinup,                                                                             &
               spinup_pt_mean,                                                                     &
               surface_pressure,                                                                   &
               timestep_scheme,                                                                    &
               tsc

    USE cpulog,                                                                                    &
       ONLY:  cpu_log,                                                                             &
              log_point,                                                                           &
              log_point_s

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_top_ind

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               precipitation

    USE kinds

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time,                                                                      &
               seconds_per_day

    USE pegrid

    USE radiation_model_mod,                                                                       &
        ONLY:  albedo,                                                                             &
               average_radiation,                                                                  &
               calc_zenith,                                                                        &
               cos_zenith,                                                                         &
               emissivity,                                                                         &
               radiation,                                                                          &
               radiation_calc_diffusion_radiation,                                                 &
               radiation_called,                                                                   &
               radiation_interactions,                                                             &
               radiation_scheme,                                                                   &
               rad_sw_in_diff,                                                                     &
               rad_sw_in_dir,                                                                      &
               rad_lw_in,                                                                          &
               sun_dir_lat,                                                                        &
               sun_dir_lon,                                                                        &
               unscheduled_radiation_calls

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io,                                                                         &
               wrd_mpi_io

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  calc_ol,                                                                            &
               calc_rib,                                                                           &
               psi_h,                                                                              &
               psi_m

    USE surface_mod,                                                                               &
        ONLY:  fr_urb,                                                                             &
               surf_lsm,                                                                           &
               surf => surf_slurb   ! renamed internally to shorten line lengths in the module


!
!-- Target arrays for timelevel switching.
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  m_liq_road_1  !< target array for liquid water reservoir on roads
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  m_liq_road_2  !< target array for liquid water reservoir on roads
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  m_liq_roof_1  !< target array for liquid water reservoir on roofs
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  m_liq_roof_2  !< target array for liquid water reservoir on roofs
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  q_can_1       !< target array for canyon water mixing ratio
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  q_can_2       !< target array for canyon water mixing ratio
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  t_can_1  !< target array for canyon temperature used to change timelevels
    REAL(wp), DIMENSION(:), TARGET, ALLOCATABLE ::  t_can_2  !< target array for canyon temperature

    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_road_1    !< target array for road temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_road_2    !< target array for road temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_roof_1    !< target array for roof temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_roof_2    !< target array for roof temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_wall_a_1  !< target array for wall A temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_wall_a_2  !< target array for wall A temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_wall_b_1  !< target array for wall B temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_wall_b_2  !< target array for wall B temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_win_a_1   !< target array for window A temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_win_a_2   !< target array for window A temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_win_b_1   !< target array for window B temperature
    REAL(wp), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_win_b_2   !< target array for window B temperature
!
!-- Arrays for output temporal averaging.
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  albedo_urb_av         !< road liquid water coverage
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq_road_av         !< road liquid water coverage
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq_roof_av         !< roof liquid water coverage
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  emiss_urb_av          !< road liquid water coverage
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_road_av           !< heat flux between the road bottom layer and soil
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_roof_av           !< heat flux between the roof bottom layer and indoor air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_wall_a_av         !< heat flux between the wall a bottom layer and indoor air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_wall_b_av         !< heat flux between the wall b bottom layer and indoor air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_win_a_av          !< heat flux between the window a bottom layer and indoor air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_win_b_av          !< heat flux between the window b bottom layer and indoor air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  m_liq_road_av         !< liquid water reservoir on roads
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  m_liq_roof_av         !< liquid water reservoir on roofs
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_can_av             !< street canyon top obukhov length
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_road_av            !< road obukhov length
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_roof_av            !< roof obukhov length
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_urb_av             !< urban obukhov length for momentum flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_can_av             !< street canyon air potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_road_av            !< road surface potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_roof_av            !< roof surface potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_wall_a_av          !< wall a surface potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_wall_b_av          !< wall b surface potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_win_a_av           !< window a surface potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_win_b_av           !< window b surface potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_can_av              !< street canyon water vapour mixing ratio
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_road_av             !< road surface mixing ratio
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_roof_av             !< roof surface mixing ratio
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs_road_av            !< road surface saturation mixing ratio
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs_roof_av            !< roof surface saturation mixing ratio
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_can_av           !< latent heat flux between the street canyon and the atmosphere
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_external_av      !< latent heat flux external to the model (e.g. from industry)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_road_av          !< latent heat flux between the road and the street canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_roof_av          !< latent heat flux between the roof and the atmosphere
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_urb_av           !< urban aggregated latent heat flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_road_av    !< road surface net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_roof_av    !< roof surface net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_urb_av     !< urban aggergated net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_wall_a_av  !< wall a surface net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_wall_b_av  !< wall b surface net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_win_a_av   !< window a surface net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_win_b_av   !< window b surface net longwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_road_av    !< road surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_roof_av    !< roof surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_urb_av     !< aggegated urban surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_wall_a_av  !< wall a surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_wall_b_av  !< wall b surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_win_a_av   !< window a surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_win_b_av   !< window b surface net shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_tr_win_a_av    !< window a surface transmitted shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_tr_win_b_av    !< window b surface transmitted shortwave radiative flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_can_av            !< street canyon aerodynamic resistance for heat
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_road_av           !< road aerodynamic resistance for heat
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_roof_av           !< roof aerodynamic resistance for heat
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_wall_a_av         !< wall A aerodynamic resistance for heat (DOE-2)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_wall_b_av         !< wall B aerodynamic resistance for heat (DOE-2)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_win_a_av          !< window A aerodynamic resistance for heat (DOE-2)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_win_b_av          !< window B aerodynamic resistance for heat (DOE-2)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_facade_av         !< wall and window aerodynamic resistance for heat (combined)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ram_urb_av            !< urban aerodynamic resistance for momentum
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_can_av            !< street canyon top bulk richardson number
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_road_av           !< road bulk richardson number
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_roof_av           !< roof bulk richardson number
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_can_av            !< sensible heat flux between the street canyon and the atmosphere
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_external_av       !< sensible heat flux external to the model (e.g. from industry)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_road_av           !< sensible heat flux between the road and the street canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_roof_av           !< sensible heat flux between the roof and the atmosphere
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_traffic_av        !< sensible heat flux from traffic to the canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_urb_av            !< urban aggregated sensible heat flux
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_wall_a_av         !< sensible heat flux between the wall a and the canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_wall_b_av         !< sensible heat flux between the wall b and the canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_win_a_av          !< sensible heat flux between the window a and the canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_win_b_av          !< sensible heat flux between the window b and the canyon air
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_2m_urb_av           !< extrapolated 2-metre urban surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_c_urb_av            !< complete urban surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_can_av              !< street canyon air temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_h_urb_av            !< effective urban surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_rad_urb_av          !< effective urban surface radiative temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_road_av        !< road surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_roof_av        !< roof surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_wall_a_av      !< wall a surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_wall_b_av      !< wall b surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_win_a_av       !< window a surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_win_b_av       !< window b surface temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_can_av             !< friction velocity for street canyons
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_road_av            !< friction velocity for roads
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_roof_av            !< friction velocity for roofs
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_urb_av             !< urban friction velocity
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  usws_urb_av           !< urban momentum flux (u-component)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  uv_abs_can_av         !< street canyon wind speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  uv_eff_can_av         !< street canyon effective wind speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt_can_av            !< street canyon air virtual potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt_road_av           !< road surface virtual potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt_roof_av           !< roof surface virtual potential temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsws_urb_av           !< urban momentum flux (v-component
!
!-- Arrays for output of unmodified LSM fluxes (2D).
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  shf_lsm_av   !< sensible heat flux from lsm surfaces
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  qsws_lsm_av  !< latent heat flux from lsm surfaces

!
!-- Arrays for output temporal averaging (2D).
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_road_av    !< road temperature (all layers)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_roof_av    !< roof temperature (all layers)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_wall_a_av  !< wall a temperature (all layers)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_wall_b_av  !< wall b temperature (all layers)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_win_a_av   !< window a temperature (all layers)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_win_b_av   !< window b temperature (all layers)

!
!-- Derived type for temporally dynamic SLUrb input variables.
    TYPE slurb_dynamic_var_type
       INTEGER(iwp) ::  lod = 0  !< 1 = no spatial dependency, 2 = spatial dependency

       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  var1d  !< 1D input (time)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var2d  !< 2D input (time,m)
    END TYPE slurb_dynamic_var_type

!
!-- Derived type for temporally dynamic SLUrb inputs.
    TYPE slurb_dynamic_type
       INTEGER(iwp) ::  ntime = 0  !< number of time steps in the input

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time  !< time dimension (seconds from init time)

       TYPE(slurb_dynamic_var_type) ::  qsws_external  !< latent heat flux from model external sources (e.g. industry)
       TYPE(slurb_dynamic_var_type) ::  shf_traffic    !< sensible heat flux from traffic
       TYPE(slurb_dynamic_var_type) ::  shf_external   !< sensible heat flux from model external sources (e.g. industry)
    END TYPE slurb_dynamic_type

    TYPE(slurb_dynamic_type) ::  slurb_dynamic

!
!-- Default subsurface layer configuration.
    INTEGER(iwp) ::  nzt_wall  !< top of the wall model (outer surface)
    INTEGER(iwp) ::  nzb_wall  !< bottom of the wall model (inside surface)
    INTEGER(iwp) ::  nzt_win   !< top of the window model (outer surface)
    INTEGER(iwp) ::  nzb_win   !< bottom of the window model (inside surface)
    INTEGER(iwp) ::  nzt_roof  !< top of the roof model (outer surface)
    INTEGER(iwp) ::  nzb_roof  !< bottom the roof model (inside surface)
    INTEGER(iwp) ::  nzt_road  !< top of the road model
    INTEGER(iwp) ::  nzb_road  !< bottom of the road model

    REAL(wp) ::  dt_slurb = HUGE( 1.0_wp )  !< maximum allowed timestep of SLUrb

!
!-- Model constants.
    REAL(wp) ::  drho_l_lv  !< (rho_l * l_v)**-1
    REAL(wp) ::  rho_lv     !< rho_surface * l_v

!
!-- Parameter defaults.
    REAL(wp), PARAMETER ::  m_liq_max_road = 1.0E-3_wp  !< maximum capacity of the liquid water reservoir on roads (m)
    REAL(wp), PARAMETER ::  m_liq_max_roof = 1.0E-3_wp  !< maximum capacity of the liquid water reservoir on roofs (m)
    REAL(wp), PARAMETER ::  rah_max   = 1.0E6_wp        !< maximum aerodynamic resistance for scalars
    REAL(wp), PARAMETER ::  rah_min   = 1.0_wp          !< minimum aerodynamic resistance for scalars
    REAL(wp), PARAMETER ::  ram_min   = 1.0_wp          !< minimum aerodynamic resistance for momentum
    REAL(wp), PARAMETER ::  urb_thres = 1.0E-2_wp       !< minimum urban fraction to consider (1%)
    REAL(wp), PARAMETER ::  us_min    = 1.0E-8_wp       !< minimum friction velocity
    REAL(wp), PARAMETER ::  zeta_min  = 1.0E-3_wp       !< minimum stability parameter absolute value (neutral limit)
!
!-- slurb_parameters namelist defaults.
    CHARACTER(LEN=20) ::  aero_roughness_heat = 'kanda'                !< SLURrb namelist parameter
    CHARACTER(LEN=20) ::  facade_resistance_parametrization = 'doe-2'  !< SLURrb namelist parameter
    CHARACTER(LEN=20) ::  street_canyon_wspeed_factor = 'surfex'       !< SLURrb namelist parameter

    INTEGER(iwp) ::  building_type = 2     !< SLURrb namelist parameter
    INTEGER(iwp) ::  n_layers_roads = 4    !< SLURrb namelist parameter
    INTEGER(iwp) ::  n_layers_roofs = 4    !< SLURrb namelist parameter
    INTEGER(iwp) ::  n_layers_walls = 4    !< SLURrb namelist parameter
    INTEGER(iwp) ::  n_layers_windows = 4  !< SLURrb namelist parameter
    INTEGER(iwp) ::  pavement_type = 2     !< SLURrb namelist parameter

    LOGICAL ::  anisotropic_street_canyons = .FALSE.  !< SLURrb namelist parameter
    LOGICAL ::  moist_physics = .TRUE.                !< SLURrb namelist parameter

    REAL(wp) ::  building_frontal_area_fraction = -9999.0_wp  !< SLURrb namelist parameter
    REAL(wp) ::  building_height = -9999.0_wp                 !< SLURrb namelist parameter
    REAL(wp) ::  building_indoor_temperature = 295.15_wp      !< SLURrb namelist parameter
    REAL(wp) ::  building_plan_area_fraction = -9999.0_wp     !< SLURrb namelist parameter
    REAl(wp) ::  deep_soil_temperature = -9999.0_wp           !< SLURrb namelist parameter
    REAL(wp) ::  qsws_external = 0.0_wp                       !< SLURrb namelist parameter
    REAL(wp) ::  shf_external = 0.0_wp                        !< SLURrb namelist parameter
    REAL(wp) ::  shf_traffic = 0.0_wp                         !< SLURrb namelist parameter
    REAL(wp) ::  street_canyon_aspect_ratio = -9999.0_wp      !< SLURrb namelist parameter
    REAL(wp) ::  street_canyon_orientation = -9999.0_wp       !< SLURrb namelist parameter
    REAL(wp) ::  urban_fraction = -9999.0_wp                  !< SLURrb namelist parameter
    REAL(wp) ::  urban_roughness_length = -9999.0_wp          !< SLURrb namelist parameter
    REAL(wp) ::  window_fraction = -9999.0_wp                 !< SLURrb namelist parameter
!
!-- Internal logical switches for character-based namelist settings.
    LOGICAL ::  facade_rah_doe       = .FALSE.  !< facade resistance parameterization using DOE-2
    LOGICAL ::  facade_rah_kray      = .FALSE.  !< facade resistance parameterization using Krayenhoff&Voogt (2007)
    LOGICAL ::  facade_rah_rowley    = .FALSE.  !< facade resistance parameterization using Rowley (1932)
    LOGICAL ::  roughness_kanda      = .FALSE.  !< roughness parameterization of horizontal surfaces using Kanda et al. (2007)
    LOGICAL ::  uv_can_factor_kray   = .FALSE.  !< street canyon wind speed factor following Krayenhoff&Voogt (2007)
    LOGICAL ::  uv_can_factor_masson = .FALSE.  !< street canyon wind speed factor following Masson (2000)
    LOGICAL ::  uv_can_factor_surfex = .FALSE.  !< street canyon wind speed factor following the SURFEX model
!
!-- Internal switches to enable computation of on-demand statistics if needed for output.
    LOGICAL ::  calc_t_2m = .FALSE.  !< enables computation of extrapolated 2-metre air temperature
    LOGICAL ::  calc_t_c  = .FALSE.  !< enables computation of complete surface temperature
    LOGICAL ::  calc_t_h  = .FALSE.  !< enables computation of effective surface temperature

!
!-- Internal switch to designate if that this is the first call of the model, which forces the call
!-- of the internal shortwave model and prevents usage of some variables which do not yet have a
!-- physical initialization.
    LOGICAL ::  first_call = .TRUE.

!
!-- Default surface description.
    REAL(wp), DIMENSION(0:45,1:6) ::  building_pars_slurb  !< building default parameters derived from USM
    REAL(wp), DIMENSION(0:14,1:5) ::  pavement_pars_slurb  !< pavement default parameters derived from LSM

    SAVE

    PRIVATE

!
!-  Public subroutines and funtions.
    PUBLIC  dt_slurb,                                                                              &
            n_layers_roofs,                                                                        &
            n_layers_walls,                                                                        &
            n_layers_windows,                                                                      &
            n_layers_roads,                                                                        &
            slurb_3d_data_averaging,                                                               &
            slurb_check_data_output,                                                               &
            slurb_check_parameters,                                                                &
            slurb_data_output_2d,                                                                  &
            slurb_data_output_3d,                                                                  &
            slurb_define_netcdf_grid,                                                              &
            slurb_header,                                                                          &
            slurb_init,                                                                            &
            slurb_init_arrays,                                                                     &
            slurb_model,                                                                           &
            slurb_parin,                                                                           &
            slurb_swap_timelevel,                                                                  &
            slurb_rrd_local,                                                                       &
            slurb_wrd_local,                                                                       &
            slurb_timestep

    INTERFACE slurb_3d_data_averaging
       MODULE PROCEDURE slurb_3d_data_averaging
    END INTERFACE slurb_3d_data_averaging

    INTERFACE slurb_atmospheric_model_coupler
       MODULE PROCEDURE slurb_atmospheric_model_coupler
    END INTERFACE slurb_atmospheric_model_coupler

    INTERFACE slurb_canyon_model
       MODULE PROCEDURE slurb_canyon_model
    END INTERFACE slurb_canyon_model

    INTERFACE slurb_check_data_output
       MODULE PROCEDURE slurb_check_data_output
    END INTERFACE slurb_check_data_output

    INTERFACE slurb_check_parameters
       MODULE PROCEDURE slurb_check_parameters
    END INTERFACE slurb_check_parameters

    INTERFACE slurb_data_output_2d
       MODULE PROCEDURE slurb_data_output_2d
    END INTERFACE slurb_data_output_2d

    INTERFACE slurb_data_output_3d
       MODULE PROCEDURE slurb_data_output_3d
    END INTERFACE slurb_data_output_3d

    INTERFACE slurb_define_netcdf_grid
       MODULE PROCEDURE slurb_define_netcdf_grid
    END INTERFACE slurb_define_netcdf_grid

    INTERFACE slurb_energy_balance_model
       MODULE PROCEDURE slurb_energy_balance_model
    END INTERFACE slurb_energy_balance_model

    INTERFACE slurb_header
       MODULE PROCEDURE slurb_header
    END INTERFACE slurb_header

    INTERFACE slurb_init
       MODULE PROCEDURE slurb_init
    END INTERFACE slurb_init

    INTERFACE slurb_init_arrays
       MODULE PROCEDURE slurb_init_arrays
    END INTERFACE slurb_init_arrays

    INTERFACE slurb_model
       MODULE PROCEDURE slurb_model
    END INTERFACE slurb_model

    INTERFACE slurb_parin
       MODULE PROCEDURE slurb_parin
    END INTERFACE slurb_parin

    INTERFACE slurb_resistance_model
       MODULE PROCEDURE slurb_resistance_model
    END INTERFACE slurb_resistance_model

     INTERFACE slurb_rrd_local
        MODULE PROCEDURE slurb_rrd_local_ftn
        MODULE PROCEDURE slurb_rrd_local_mpi
     END INTERFACE slurb_rrd_local

    INTERFACE slurb_swap_timelevel
       MODULE PROCEDURE slurb_swap_timelevel
    END INTERFACE slurb_swap_timelevel

    INTERFACE slurb_update_external_vars
       MODULE PROCEDURE slurb_update_external_vars
    END INTERFACE slurb_update_external_vars

    INTERFACE slurb_urban_aggregation_model
       MODULE PROCEDURE slurb_urban_aggregation_model
    END INTERFACE slurb_urban_aggregation_model

    INTERFACE slurb_wrd_local
       MODULE PROCEDURE slurb_wrd_local
    END INTERFACE slurb_wrd_local

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for averaging 3D data (and also 2D in reality).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_3d_data_averaging( mode, variable )

    CHARACTER (LEN=*) ::  mode      !< mode of function
    CHARACTER (LEN=*) ::  variable  !< variable name

    INTEGER(iwp) ::  i   !< loop index (x-direction)
    INTEGER(iwp) ::  j   !< loop index (y-direction)
    INTEGER(iwp) ::  k   !< loop index (z-direction)
    INTEGER(iwp) ::  m   !< running SLUrb tile index

    IF ( variable(1:6) /= 'slurb_' ) RETURN

    IF ( mode == 'allocate' )  THEN
       SELECT CASE ( TRIM( variable ) )

          CASE ( 'slurb_albedo_urb*' )
             IF ( .NOT. ALLOCATED( albedo_urb_av ) )  THEN
                ALLOCATE( albedo_urb_av(1:surf%ns) )
             ENDIF
             albedo_urb_av(:) = 0.0_wp

          CASE ( 'slurb_c_liq_road*' )
             IF ( .NOT. ALLOCATED( c_liq_road_av ) )  THEN
                ALLOCATE( c_liq_road_av(1:surf%ns) )
             ENDIF
              c_liq_road_av(:) = 0.0_wp

          CASE ( 'slurb_c_liq_roof*' )
             IF ( .NOT. ALLOCATED( c_liq_roof_av ) )  THEN
                ALLOCATE( c_liq_roof_av(1:surf%ns) )
             ENDIF
             c_liq_roof_av(:) = 0.0_wp

          CASE ( 'slurb_emiss_urb*' )
             IF ( .NOT. ALLOCATED( emiss_urb_av ) )  THEN
                ALLOCATE( emiss_urb_av(1:surf%ns) )
             ENDIF
             emiss_urb_av(:) = 0.0_wp

          CASE ( 'slurb_ghf_road*' )
             IF ( .NOT. ALLOCATED( ghf_road_av ) )  THEN
                ALLOCATE( ghf_road_av(1:surf%ns) )
             ENDIF
             ghf_road_av(:) = 0.0_wp

          CASE ( 'slurb_ghf_roof*' )
             IF ( .NOT. ALLOCATED( ghf_roof_av ) )  THEN
                ALLOCATE( ghf_roof_av(1:surf%ns) )
             ENDIF
             ghf_roof_av(:) = 0.0_wp

          CASE ( 'slurb_ghf_wall_a*' )
             IF ( .NOT. ALLOCATED( ghf_wall_a_av ) )  THEN
                ALLOCATE( ghf_wall_a_av(1:surf%ns) )
             ENDIF
             ghf_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_ghf_wall_b*' )
             IF ( .NOT. ALLOCATED( ghf_wall_b_av ) )  THEN
                ALLOCATE( ghf_wall_b_av(1:surf%ns) )
             ENDIF
             ghf_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_ghf_win_a*' )
             IF ( .NOT. ALLOCATED( ghf_win_a_av ) )  THEN
                ALLOCATE( ghf_win_a_av(1:surf%ns) )
             ENDIF
             ghf_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_ghf_win_b*' )
             IF ( .NOT. ALLOCATED( ghf_win_b_av ) )  THEN
                ALLOCATE( ghf_win_b_av(1:surf%ns) )
             ENDIF
             ghf_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_m_liq_road*' )
             IF ( .NOT. ALLOCATED( m_liq_road_av ) )  THEN
                ALLOCATE( m_liq_road_av(1:surf%ns) )
             ENDIF
             m_liq_road_av(:) = 0.0_wp

          CASE ( 'slurb_m_liq_roof*' )
             IF ( .NOT. ALLOCATED( m_liq_roof_av ) )  THEN
                ALLOCATE( m_liq_roof_av(1:surf%ns) )
             ENDIF
             m_liq_roof_av(:) = 0.0_wp

          CASE ( 'slurb_ol_canyon*' )
             IF ( .NOT. ALLOCATED( ol_can_av ) )  THEN
                ALLOCATE( ol_can_av(1:surf%ns) )
             ENDIF
             ol_can_av(:) = 0.0_wp

          CASE ( 'slurb_ol_road*' )
             IF ( .NOT. ALLOCATED( ol_road_av ) )  THEN
                ALLOCATE( ol_road_av(1:surf%ns) )
             ENDIF
             ol_road_av(:) = 0.0_wp

          CASE ( 'slurb_ol_roof*' )
             IF ( .NOT. ALLOCATED( ol_roof_av ) )  THEN
                ALLOCATE( ol_roof_av(1:surf%ns) )
             ENDIF
             ol_roof_av(:) = 0.0_wp

          CASE ( 'slurb_ol_urb*' )
             IF ( .NOT. ALLOCATED( ol_urb_av ) )  THEN
                ALLOCATE( ol_urb_av(1:surf%ns) )
             ENDIF
             ol_urb_av(:) = 0.0_wp

          CASE ( 'slurb_q_canyon*' )
             IF ( .NOT. ALLOCATED( q_can_av ) )  THEN
                ALLOCATE( q_can_av(1:surf%ns) )
             ENDIF
             q_can_av(:) = 0.0_wp

          CASE ( 'slurb_q_road*' )
             IF ( .NOT. ALLOCATED( q_road_av ) )  THEN
                ALLOCATE( q_road_av(1:surf%ns) )
             ENDIF
             q_road_av(:) = 0.0_wp

          CASE ( 'slurb_q_roof*' )
             IF ( .NOT. ALLOCATED( q_roof_av ) )  THEN
                ALLOCATE( q_roof_av(1:surf%ns) )
             ENDIF
             q_roof_av(:) = 0.0_wp

          CASE ( 'slurb_qs_road*' )
             IF ( .NOT. ALLOCATED( qs_road_av ) )  THEN
                ALLOCATE( qs_road_av(1:surf%ns) )
             ENDIF
             qs_road_av(:) = 0.0_wp

          CASE ( 'slurb_qs_roof*' )
             IF ( .NOT. ALLOCATED( qs_roof_av ) )  THEN
                ALLOCATE( qs_roof_av(1:surf%ns) )
             ENDIF
             qs_roof_av(:) = 0.0_wp

          CASE ( 'slurb_qsws_canyon*' )
             IF ( .NOT. ALLOCATED( qsws_can_av ) )  THEN
                ALLOCATE( qsws_can_av(1:surf%ns) )
             ENDIF
             qsws_can_av(:) = 0.0_wp

          CASE ( 'slurb_qsws_road*' )
             IF ( .NOT. ALLOCATED( qsws_road_av ) )  THEN
                ALLOCATE( qsws_road_av(1:surf%ns) )
             ENDIF
             qsws_road_av(:) = 0.0_wp

          CASE ( 'slurb_qsws_roof*' )
             IF ( .NOT. ALLOCATED( qsws_roof_av ) )  THEN
                ALLOCATE( qsws_roof_av(1:surf%ns) )
             ENDIF
             qsws_roof_av(:) = 0.0_wp

          CASE ( 'slurb_qsws_lsm*' )
             IF ( .NOT. ALLOCATED( qsws_lsm_av ) )  THEN
                ALLOCATE( qsws_lsm_av(nys:nyn,nxl:nxr) )
             ENDIF
             qsws_lsm_av(:,:) = 0.0_wp

          CASE ( 'slurb_qsws_urb*' )
             IF ( .NOT. ALLOCATED( qsws_urb_av ) )  THEN
                ALLOCATE( qsws_urb_av(1:surf%ns) )
             ENDIF
             qsws_urb_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_road*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_road_av ) )  THEN
                ALLOCATE( rad_lw_net_road_av(1:surf%ns) )
             ENDIF
             rad_lw_net_road_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_roof*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_roof_av ) )  THEN
                ALLOCATE( rad_lw_net_roof_av(1:surf%ns) )
             ENDIF
             rad_lw_net_roof_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_urb*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_urb_av ) )  THEN
                ALLOCATE( rad_lw_net_urb_av(1:surf%ns) )
             ENDIF
             rad_lw_net_urb_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_wall_a*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_wall_a_av ) )  THEN
                ALLOCATE( rad_lw_net_wall_a_av(1:surf%ns) )
             ENDIF
             rad_lw_net_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_wall_b*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_wall_b_av ) )  THEN
                ALLOCATE( rad_lw_net_wall_b_av(1:surf%ns) )
             ENDIF
             rad_lw_net_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_win_a*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_win_a_av ) )  THEN
                ALLOCATE( rad_lw_net_win_a_av(1:surf%ns) )
             ENDIF
             rad_lw_net_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_rad_lw_net_win_b*' )
             IF ( .NOT. ALLOCATED( rad_lw_net_win_b_av ) )  THEN
                ALLOCATE( rad_lw_net_win_b_av(1:surf%ns) )
             ENDIF
             rad_lw_net_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_road*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_road_av ) )  THEN
                ALLOCATE( rad_sw_net_road_av(1:surf%ns) )
             ENDIF
             rad_sw_net_road_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_roof*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_roof_av ) )  THEN
                ALLOCATE( rad_sw_net_roof_av(1:surf%ns) )
             ENDIF
             rad_sw_net_roof_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_urb*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_urb_av ) )  THEN
                ALLOCATE( rad_sw_net_urb_av(1:surf%ns) )
             ENDIF
             rad_sw_net_urb_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_wall_a*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_wall_a_av ) )  THEN
                ALLOCATE( rad_sw_net_wall_a_av(1:surf%ns) )
             ENDIF
             rad_sw_net_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_wall_b*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_wall_b_av ) )  THEN
                ALLOCATE( rad_sw_net_wall_b_av(1:surf%ns) )
             ENDIF
             rad_sw_net_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_win_a*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_win_a_av ) )  THEN
                ALLOCATE( rad_sw_net_win_a_av(1:surf%ns) )
             ENDIF
             rad_sw_net_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_net_win_b*' )
             IF ( .NOT. ALLOCATED( rad_sw_net_win_b_av ) )  THEN
                ALLOCATE( rad_sw_net_win_b_av(1:surf%ns) )
             ENDIF
             rad_sw_net_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_tr_win_a*' )
             IF ( .NOT. ALLOCATED( rad_sw_tr_win_a_av ) )  THEN
                ALLOCATE( rad_sw_tr_win_a_av(1:surf%ns) )
             ENDIF
             rad_sw_tr_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_rad_sw_tr_win_b*' )
             IF ( .NOT. ALLOCATED( rad_sw_tr_win_b_av ) )  THEN
                ALLOCATE( rad_sw_tr_win_b_av(1:surf%ns) )
             ENDIF
             rad_sw_tr_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_rah_canyon*' )
             IF ( .NOT. ALLOCATED( rah_can_av ) )  THEN
                ALLOCATE( rah_can_av(1:surf%ns) )
             ENDIF
             rah_can_av(:) = 0.0_wp

          CASE ( 'slurb_rah_road*' )
             IF ( .NOT. ALLOCATED( rah_road_av ) )  THEN
                ALLOCATE( rah_road_av(1:surf%ns) )
             ENDIF
             rah_road_av(:) = 0.0_wp

          CASE ( 'slurb_rah_roof*' )
             IF ( .NOT. ALLOCATED( rah_roof_av ) )  THEN
                ALLOCATE( rah_roof_av(1:surf%ns) )
             ENDIF
             rah_roof_av(:) = 0.0_wp

          CASE ( 'slurb_rah_wall_a*' )
             IF ( .NOT. ALLOCATED( rah_wall_a_av ) )  THEN
                ALLOCATE( rah_wall_a_av(1:surf%ns) )
             ENDIF
             rah_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_rah_wall_b*' )
             IF ( .NOT. ALLOCATED( rah_wall_b_av ) )  THEN
                ALLOCATE( rah_wall_b_av(1:surf%ns) )
             ENDIF
             rah_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_rah_win_a*' )
             IF ( .NOT. ALLOCATED( rah_win_a_av ) )  THEN
                ALLOCATE( rah_win_a_av(1:surf%ns) )
             ENDIF
             rah_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_rah_win_b*' )
             IF ( .NOT. ALLOCATED( rah_win_b_av ) )  THEN
                ALLOCATE( rah_win_b_av(1:surf%ns) )
             ENDIF
             rah_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_rah_facade*' )
             IF ( .NOT. ALLOCATED( rah_facade_av ) )  THEN
                ALLOCATE( rah_facade_av(1:surf%ns) )
             ENDIF
             rah_facade_av(:) = 0.0_wp

          CASE ( 'slurb_ram_urb*' )
             IF ( .NOT. ALLOCATED( ram_urb_av ) )  THEN
                ALLOCATE( ram_urb_av(1:surf%ns) )
             ENDIF
             ram_urb_av(:) = 0.0_wp

          CASE ( 'slurb_rib_canyon*' )
             IF ( .NOT. ALLOCATED( rib_can_av ) )  THEN
                ALLOCATE( rib_can_av(1:surf%ns) )
             ENDIF
             rib_can_av(:) = 0.0_wp

          CASE ( 'slurb_rib_road*' )
             IF ( .NOT. ALLOCATED( rib_road_av ) )  THEN
                ALLOCATE( rib_road_av(1:surf%ns) )
             ENDIF
             rib_road_av(:) = 0.0_wp

          CASE ( 'slurb_rib_roof*' )
             IF ( .NOT. ALLOCATED( rib_roof_av ) )  THEN
                ALLOCATE( rib_roof_av(1:surf%ns) )
             ENDIF
             rib_roof_av(:) = 0.0_wp

          CASE ( 'slurb_shf_canyon*' )
             IF ( .NOT. ALLOCATED( shf_can_av ) )  THEN
                ALLOCATE( shf_can_av(1:surf%ns) )
             ENDIF
             shf_can_av(:) = 0.0_wp

          CASE ( 'slurb_shf_external*' )
             IF ( .NOT. ALLOCATED( shf_external_av ) )  THEN
                ALLOCATE( shf_external_av(1:surf%ns) )
             ENDIF
             shf_external_av(:) = 0.0_wp

          CASE ( 'slurb_shf_road*' )
             IF ( .NOT. ALLOCATED( shf_road_av ) )  THEN
                ALLOCATE( shf_road_av(1:surf%ns) )
             ENDIF
             shf_road_av(:) = 0.0_wp

          CASE ( 'slurb_shf_roof*' )
             IF ( .NOT. ALLOCATED( shf_roof_av ) )  THEN
                ALLOCATE( shf_roof_av(1:surf%ns) )
             ENDIF
             shf_roof_av(:) = 0.0_wp

          CASE ( 'slurb_shf_traffic*' )
             IF ( .NOT. ALLOCATED( shf_traffic_av ) )  THEN
                ALLOCATE( shf_traffic_av(1:surf%ns) )
             ENDIF
             shf_traffic_av(:) = 0.0_wp

          CASE ( 'slurb_shf_lsm*' )
             IF ( .NOT. ALLOCATED( shf_lsm_av ) )  THEN
                ALLOCATE( shf_lsm_av(nys:nyn,nxl:nxr) )
             ENDIF
             shf_lsm_av(:,:) = 0.0_wp

          CASE ( 'slurb_shf_urb*' )
             IF ( .NOT. ALLOCATED( shf_urb_av ) )  THEN
                ALLOCATE( shf_urb_av(1:surf%ns) )
             ENDIF
             shf_urb_av(:) = 0.0_wp

          CASE ( 'slurb_shf_wall_a*' )
             IF ( .NOT. ALLOCATED( shf_wall_a_av ) )  THEN
                ALLOCATE( shf_wall_a_av(1:surf%ns) )
             ENDIF
             shf_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_shf_wall_b*' )
             IF ( .NOT. ALLOCATED( shf_wall_b_av ) )  THEN
                ALLOCATE( shf_wall_b_av(1:surf%ns) )
             ENDIF
             shf_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_shf_win_a*' )
             IF ( .NOT. ALLOCATED( shf_win_a_av ) )  THEN
                ALLOCATE( shf_win_a_av(1:surf%ns) )
             ENDIF
             shf_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_shf_win_b*' )
             IF ( .NOT. ALLOCATED( shf_win_b_av ) )  THEN
                ALLOCATE( shf_win_b_av(1:surf%ns) )
             ENDIF
             shf_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_t_canyon*' )
             IF ( .NOT. ALLOCATED( t_can_av ) )  THEN
                ALLOCATE( t_can_av(1:surf%ns) )
             ENDIF
             t_can_av(:) = 0.0_wp

          CASE ( 'slurb_t_rad_urb*' )
             IF ( .NOT. ALLOCATED( t_rad_urb_av ) )  THEN
                ALLOCATE( t_rad_urb_av(1:surf%ns) )
             ENDIF
             t_rad_urb_av(:) = 0.0_wp

          CASE ( 'slurb_t_surf_road*' )
             IF ( .NOT. ALLOCATED( t_surf_road_av ) )  THEN
                ALLOCATE( t_surf_road_av(1:surf%ns) )
             ENDIF
             t_surf_road_av(:) = 0.0_wp

          CASE ( 'slurb_t_surf_roof*' )
             IF ( .NOT. ALLOCATED( t_surf_roof_av ) )  THEN
                ALLOCATE( t_surf_roof_av(1:surf%ns) )
             ENDIF
             t_surf_roof_av(:) = 0.0_wp

          CASE ( 'slurb_t_surf_wall_a*' )
             IF ( .NOT. ALLOCATED( t_surf_wall_a_av ) )  THEN
                ALLOCATE( t_surf_wall_a_av(1:surf%ns) )
             ENDIF
             t_surf_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_t_surf_wall_b*' )
             IF ( .NOT. ALLOCATED( t_surf_wall_b_av ) )  THEN
                ALLOCATE( t_surf_wall_b_av(1:surf%ns) )
             ENDIF
             t_surf_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_t_surf_win_a*' )
             IF ( .NOT. ALLOCATED( t_surf_win_a_av ) )  THEN
                ALLOCATE( t_surf_win_a_av(1:surf%ns) )
             ENDIF
             t_surf_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_t_surf_win_b*' )
             IF ( .NOT. ALLOCATED( t_surf_win_b_av ) )  THEN
                ALLOCATE( t_surf_win_b_av(1:surf%ns) )
             ENDIF
             t_surf_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_t_c_urb*' )
             IF ( .NOT. ALLOCATED( t_c_urb_av ) )  THEN
                ALLOCATE( t_c_urb_av(1:surf%ns) )
             ENDIF
             t_c_urb_av(:) = 0.0_wp

          CASE ( 'slurb_t_h_urb*' )
             IF ( .NOT. ALLOCATED( t_h_urb_av ) )  THEN
                ALLOCATE( t_h_urb_av(1:surf%ns) )
             ENDIF
             t_h_urb_av(:) = 0.0_wp

          CASE ( 'slurb_t_2m_urb*' )
             IF ( .NOT. ALLOCATED( t_2m_urb_av ) )  THEN
                ALLOCATE( t_2m_urb_av(1:surf%ns) )
             ENDIF
             t_2m_urb_av(:) = 0.0_wp

          CASE ( 'slurb_theta_canyon*' )
             IF ( .NOT. ALLOCATED( pt_can_av ) )  THEN
                ALLOCATE( pt_can_av(1:surf%ns) )
             ENDIF
             pt_can_av(:) = 0.0_wp

          CASE ( 'slurb_theta_road*' )
             IF ( .NOT. ALLOCATED( pt_road_av ) )  THEN
                ALLOCATE( pt_road_av(1:surf%ns) )
             ENDIF
             pt_road_av(:) = 0.0_wp

          CASE ( 'slurb_theta_roof*' )
             IF ( .NOT. ALLOCATED( pt_roof_av ) )  THEN
                ALLOCATE( pt_roof_av(1:surf%ns) )
             ENDIF
             pt_roof_av(:) = 0.0_wp

          CASE ( 'slurb_theta_wall_a*' )
             IF ( .NOT. ALLOCATED( pt_wall_a_av ) )  THEN
                ALLOCATE( pt_wall_a_av(1:surf%ns) )
             ENDIF
             pt_wall_a_av(:) = 0.0_wp

          CASE ( 'slurb_theta_wall_b*' )
             IF ( .NOT. ALLOCATED( pt_wall_b_av ) )  THEN
                ALLOCATE( pt_wall_b_av(1:surf%ns) )
             ENDIF
             pt_wall_b_av(:) = 0.0_wp

          CASE ( 'slurb_theta_win_a*' )
             IF ( .NOT. ALLOCATED( pt_win_a_av ) )  THEN
                ALLOCATE( pt_win_a_av(1:surf%ns) )
             ENDIF
             pt_win_a_av(:) = 0.0_wp

          CASE ( 'slurb_theta_win_b*' )
             IF ( .NOT. ALLOCATED( pt_win_b_av ) )  THEN
                ALLOCATE( pt_win_b_av(1:surf%ns) )
             ENDIF
             pt_win_b_av(:) = 0.0_wp

          CASE ( 'slurb_thetav_canyon*' )
             IF ( .NOT. ALLOCATED( vpt_can_av ) )  THEN
                ALLOCATE( vpt_can_av(1:surf%ns) )
             ENDIF
             vpt_can_av(:) = 0.0_wp

          CASE ( 'slurb_thetav_road*' )
             IF ( .NOT. ALLOCATED( vpt_road_av ) )  THEN
                ALLOCATE( vpt_road_av(1:surf%ns) )
             ENDIF
             vpt_road_av(:) = 0.0_wp

          CASE ( 'slurb_thetav_roof*' )
             IF ( .NOT. ALLOCATED( vpt_roof_av ) )  THEN
                ALLOCATE( vpt_roof_av(1:surf%ns) )
             ENDIF
             vpt_roof_av(:) = 0.0_wp

          CASE ( 'slurb_wspeed_canyon*' )
             IF ( .NOT. ALLOCATED( uv_abs_can_av ) )  THEN
                ALLOCATE( uv_abs_can_av(1:surf%ns) )
             ENDIF
             uv_abs_can_av(:) = 0.0_wp

          CASE ( 'slurb_wspeed_eff_canyon*' )
             IF ( .NOT. ALLOCATED( uv_eff_can_av ) )  THEN
                ALLOCATE( uv_eff_can_av(1:surf%ns) )
             ENDIF
             uv_eff_can_av(:) = 0.0_wp

          CASE ( 'slurb_us_canyon*' )
             IF ( .NOT. ALLOCATED( us_can_av ) )  THEN
                ALLOCATE( us_can_av(1:surf%ns) )
             ENDIF
             us_can_av(:) = 0.0_wp

          CASE ( 'slurb_us_road*' )
             IF ( .NOT. ALLOCATED( us_road_av ) )  THEN
                ALLOCATE( us_road_av(1:surf%ns) )
             ENDIF
             us_road_av(:) = 0.0_wp

          CASE ( 'slurb_us_roof*' )
             IF ( .NOT. ALLOCATED( us_roof_av ) )  THEN
                ALLOCATE( us_roof_av(1:surf%ns) )
             ENDIF
             us_roof_av(:) = 0.0_wp

          CASE ( 'slurb_us_urb*' )
             IF ( .NOT. ALLOCATED( us_urb_av ) )  THEN
                ALLOCATE( us_urb_av(1:surf%ns) )
             ENDIF
             us_urb_av(:) = 0.0_wp

          CASE ( 'slurb_usws_urb*' )
             IF ( .NOT. ALLOCATED( usws_urb_av ) )  THEN
                ALLOCATE( usws_urb_av(1:surf%ns) )
             ENDIF
             usws_urb_av(:) = 0.0_wp

          CASE ( 'slurb_vsws_urb*' )
             IF ( .NOT. ALLOCATED( vsws_urb_av ) )  THEN
                ALLOCATE( vsws_urb_av(1:surf%ns) )
             ENDIF
             vsws_urb_av(:) = 0.0_wp

          CASE ( 'slurb_t_road' )
             IF ( .NOT. ALLOCATED( t_road_av ) )  THEN
                ALLOCATE( t_road_av(nzt_road:nzb_road,1:surf%ns) )
             ENDIF
             t_road_av(:,:) = 0.0_wp

          CASE ( 'slurb_t_roof' )
             IF ( .NOT. ALLOCATED( t_roof_av ) )  THEN
                ALLOCATE( t_roof_av(nzt_roof:nzb_roof,1:surf%ns) )
             ENDIF
             t_roof_av(:,:) = 0.0_wp

          CASE ( 'slurb_t_wall_a' )
             IF ( .NOT. ALLOCATED( t_wall_a_av ) )  THEN
                ALLOCATE( t_wall_a_av(nzt_wall:nzb_wall,1:surf%ns) )
             ENDIF
             t_wall_a_av(:,:) = 0.0_wp

          CASE ( 'slurb_t_wall_b' )
             IF ( .NOT. ALLOCATED( t_wall_b_av ) )  THEN
                ALLOCATE( t_wall_b_av(nzt_wall:nzb_wall,1:surf%ns) )
             ENDIF
             t_wall_b_av(:,:) = 0.0_wp

          CASE ( 'slurb_t_win_a' )
             IF ( .NOT. ALLOCATED( t_win_a_av ) )  THEN
                ALLOCATE( t_win_a_av(nzt_win:nzb_win,1:surf%ns) )
             ENDIF
             t_win_a_av(:,:) = 0.0_wp

          CASE ( 'slurb_t_win_b' )
             IF ( .NOT. ALLOCATED( t_win_b_av ) )  THEN
                ALLOCATE( t_win_b_av(nzt_win:nzb_win,1:surf%ns) )
             ENDIF
             t_win_b_av(:,:) = 0.0_wp

       CASE DEFAULT
!
!--       In case of missing or incorrect SLUrb variable, give a meaningful error.
          IF ( variable(1:6) == 'slurb_' )  THEN
             message_string = 'Unknown temporally averaged SLUrb output '                          &
                              // TRIM ( variable ) // ' requested.'
             CALL message( 'slurb_3d_data_averaging', 'SLU0035', 1, 2, 0, 6, 0 )
          ENDIF

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'slurb_albedo_urb*' )
             IF ( ALLOCATED( albedo_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   albedo_urb_av(m) = albedo_urb_av(m) + surf%albedo_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_c_liq_road*' )
             IF ( ALLOCATED( c_liq_road_av ) )  THEN
                DO  m = 1, surf%ns
                   c_liq_road_av(m) = c_liq_road_av(m) + surf%c_liq_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_c_liq_roof*' )
             IF ( ALLOCATED( c_liq_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   c_liq_roof_av(m) = c_liq_roof_av(m) + surf%c_liq_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_emiss_urb*' )
             IF ( ALLOCATED( emiss_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   emiss_urb_av(m) = emiss_urb_av(m) + surf%emiss_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_road*' )
             IF ( ALLOCATED( ghf_road_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_road_av(m) = ghf_road_av(m) + surf%ghf_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_roof*' )
             IF ( ALLOCATED( ghf_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_roof_av(m) = ghf_roof_av(m) + surf%ghf_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_wall_a*' )
             IF ( ALLOCATED( ghf_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_wall_a_av(m) = ghf_wall_a_av(m) + surf%ghf_wall_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_wall_b*' )
             IF ( ALLOCATED( ghf_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_wall_b_av(m) = ghf_wall_b_av(m) + surf%ghf_wall_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_win_a*' )
             IF ( ALLOCATED( ghf_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_win_a_av(m) = ghf_win_a_av(m) + surf%ghf_win_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_win_b*' )
             IF ( ALLOCATED( ghf_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_win_b_av(m) = ghf_win_b_av(m) + surf%ghf_win_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_m_liq_road*' )
             IF ( ALLOCATED( m_liq_road_av ) )  THEN
                DO  m = 1, surf%ns
                   m_liq_road_av(m) = m_liq_road_av(m) + surf%m_liq_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_m_liq_roof*' )
             IF ( ALLOCATED( m_liq_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   m_liq_roof_av(m) = m_liq_roof_av(m) + surf%m_liq_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_canyon*' )
             IF ( ALLOCATED( ol_can_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_can_av(m) = ol_can_av(m) + surf%ol_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_road*' )
             IF ( ALLOCATED( ol_road_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_road_av(m) = ol_road_av(m) + surf%ol_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_roof*' )
             IF ( ALLOCATED( ol_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_roof_av(m) = ol_roof_av(m) + surf%ol_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_urb*' )
             IF ( ALLOCATED( ol_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_urb_av(m) = ol_urb_av(m) + surf%ol_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_q_canyon*' )
             IF ( ALLOCATED( q_can_av ) )  THEN
                DO  m = 1, surf%ns
                   q_can_av(m) = q_can_av(m) + surf%q_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_q_road*' )
             IF ( ALLOCATED( q_road_av ) )  THEN
                DO  m = 1, surf%ns
                   q_road_av(m) = q_road_av(m) + surf%q_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_q_roof*' )
             IF ( ALLOCATED( q_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   q_roof_av(m) = q_roof_av(m) + surf%q_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_qs_road*' )
             IF ( ALLOCATED( qs_road_av ) )  THEN
                DO  m = 1, surf%ns
                   qs_road_av(m) = qs_road_av(m) + surf%qs_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_qs_roof*' )
             IF ( ALLOCATED( qs_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   qs_roof_av(m) = qs_roof_av(m) + surf%qs_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_canyon*' )
             IF ( ALLOCATED( qsws_can_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_can_av(m) = qsws_can_av(m) + surf%qsws_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_road*' )
             IF ( ALLOCATED( qsws_road_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_road_av(m) = qsws_road_av(m) + surf%qsws_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_roof*' )
             IF ( ALLOCATED( qsws_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_roof_av(m) = qsws_roof_av(m) + surf%qsws_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_lsm*' )
             IF ( ALLOCATED( qsws_lsm_av ) )  THEN
                DO  m = 1, surf_lsm%ns
                   i = surf_lsm%i(m)
                   j = surf_lsm%j(m)
                   qsws_lsm_av(j,i) = qsws_lsm_av(j,i) +                                           &
                                      MERGE( surf_lsm%qsws(m), 0.0_wp, surf_lsm%upward(m) )
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_urb*' )
             IF ( ALLOCATED( qsws_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_urb_av(m) = qsws_urb_av(m) + surf%qsws_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_road*' )
             IF ( ALLOCATED( rad_lw_net_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_road_av(m) = rad_lw_net_road_av(m) + surf%rad_lw_net_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_roof*' )
             IF ( ALLOCATED( rad_lw_net_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_roof_av(m) = rad_lw_net_roof_av(m) + surf%rad_lw_net_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_urb*' )
             IF ( ALLOCATED( rad_lw_net_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_urb_av(m) = rad_lw_net_urb_av(m) + surf%rad_lw_net_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_wall_a*' )
             IF ( ALLOCATED( rad_lw_net_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_wall_a_av(m) = rad_lw_net_wall_a_av(m) + surf%rad_lw_net_wall_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_wall_b*' )
             IF ( ALLOCATED( rad_lw_net_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_wall_b_av(m) = rad_lw_net_wall_b_av(m) + surf%rad_lw_net_wall_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_win_a*' )
             IF ( ALLOCATED( rad_lw_net_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_win_a_av(m) = rad_lw_net_win_a_av(m) + surf%rad_lw_net_win_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_win_b*' )
             IF ( ALLOCATED( rad_lw_net_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_win_b_av(m) = rad_lw_net_win_b_av(m) + surf%rad_lw_net_win_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_road*' )
             IF ( ALLOCATED( rad_sw_net_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_road_av(m) = rad_sw_net_road_av(m) + surf%rad_sw_net_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_roof*' )
             IF ( ALLOCATED( rad_sw_net_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_roof_av(m) = rad_sw_net_roof_av(m) + surf%rad_sw_net_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_urb*' )
             IF ( ALLOCATED( rad_sw_net_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_urb_av(m) = rad_sw_net_urb_av(m) + surf%rad_sw_net_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_wall_a*' )
             IF ( ALLOCATED( rad_sw_net_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_wall_a_av(m) = rad_sw_net_wall_a_av(m) + surf%rad_sw_net_wall_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_wall_b*' )
             IF ( ALLOCATED( rad_sw_net_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_wall_b_av(m) = rad_sw_net_wall_b_av(m) + surf%rad_sw_net_wall_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_win_a*' )
             IF ( ALLOCATED( rad_sw_net_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_win_a_av(m) = rad_sw_net_win_a_av(m) + surf%rad_sw_net_win_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_win_b*' )
             IF ( ALLOCATED( rad_sw_net_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_win_b_av(m) = rad_sw_net_win_b_av(m) + surf%rad_sw_net_win_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_tr_win_a*' )
             IF ( ALLOCATED( rad_sw_tr_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_tr_win_a_av(m) = rad_sw_tr_win_a_av(m) + surf%rad_sw_net_win_a(m) *      &
                                           ( 1.0_wp - SUM( surf%absorption_win(:,m) ) )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_tr_win_b*' )
             IF ( ALLOCATED( rad_sw_tr_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_tr_win_b_av(m) = rad_sw_tr_win_b_av(m) + surf%rad_sw_net_win_b(m) *      &
                                           ( 1.0_wp - SUM( surf%absorption_win(:,m) ) )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_canyon*' )
             IF ( ALLOCATED( rah_can_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_can_av(m) = rah_can_av(m) + surf%rah_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_road*' )
             IF ( ALLOCATED( rah_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_road_av(m) = rah_road_av(m) + surf%rah_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_roof*' )
             IF ( ALLOCATED( rah_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_roof_av(m) = rah_roof_av(m) + surf%rah_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_wall_a*' )
             IF ( ALLOCATED( rah_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_wall_a_av(m) = rah_wall_a_av(m) + surf%rah_wall_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_wall_b*' )
             IF ( ALLOCATED( rah_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_wall_b_av(m) = rah_wall_b_av(m) + surf%rah_wall_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_win_a*' )
             IF ( ALLOCATED( rah_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_win_a_av(m) = rah_win_a_av(m) + surf%rah_win_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_win_b*' )
             IF ( ALLOCATED( rah_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_win_b_av(m) = rah_win_b_av(m) + surf%rah_win_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_facade*' )
             IF ( ALLOCATED( rah_facade_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_facade_av(m) = rah_facade_av(m) + surf%rah_facade(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_ram_urb*' )
             IF ( ALLOCATED( ram_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   ram_urb_av(m) = ram_urb_av(m) + surf%ram_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rib_canyon*' )
             IF ( ALLOCATED( rib_can_av ) )  THEN
                DO  m = 1, surf%ns
                   rib_can_av(m) = rib_can_av(m) + surf%rib_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rib_road*' )
             IF ( ALLOCATED( rib_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rib_road_av(m) = rib_road_av(m) + surf%rib_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_rib_roof*' )
             IF ( ALLOCATED( rib_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rib_roof_av(m) = rib_roof_av(m) + surf%rib_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_canyon*' )
             IF ( ALLOCATED( shf_can_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_can_av(m) = shf_can_av(m) + surf%shf_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_external*' )
             IF ( ALLOCATED( shf_external_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_external_av(m) = shf_external_av(m) + surf%shf_external(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_road*' )
             IF ( ALLOCATED( shf_road_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_road_av(m) = shf_road_av(m) + surf%shf_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_roof*' )
             IF ( ALLOCATED( shf_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_roof_av(m) = shf_roof_av(m) + surf%shf_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_traffic*' )
             IF ( ALLOCATED( shf_traffic_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_traffic_av(m) = shf_traffic_av(m) + surf%shf_traffic(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_urb*' )
             IF ( ALLOCATED( shf_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_urb_av(m) = shf_urb_av(m) + surf%shf_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_lsm*' )
             IF ( ALLOCATED( shf_lsm_av ) )  THEN
                DO  m = 1, surf_lsm%ns
                   i = surf_lsm%i(m)
                   j = surf_lsm%j(m)
                   shf_lsm_av(j,i) = shf_lsm_av(j,i) +                                             &
                                     MERGE( surf_lsm%shf(m), 0.0_wp, surf_lsm%upward(m) )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_wall_a*' )
             IF ( ALLOCATED( shf_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_wall_a_av(m) = shf_wall_a_av(m) + surf%shf_wall_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_wall_b*' )
             IF ( ALLOCATED( shf_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_wall_b_av(m) = shf_wall_b_av(m) + surf%shf_wall_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_win_a*' )
             IF ( ALLOCATED( shf_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_win_a_av(m) = shf_win_a_av(m) + surf%shf_win_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_win_b*' )
             IF ( ALLOCATED( shf_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_win_b_av(m) = shf_win_b_av(m) + surf%shf_win_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_canyon*' )
             IF ( ALLOCATED( t_can_av ) )  THEN
                DO  m = 1, surf%ns
                   t_can_av(m) = t_can_av(m) + surf%t_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_rad_urb*' )
             IF ( ALLOCATED( t_rad_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_rad_urb_av(m) = t_rad_urb_av(m) + surf%t_rad_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_road*' )
             IF ( ALLOCATED( t_surf_road_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_road_av(m) = t_surf_road_av(m) + surf%t_road(nzt_road,m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_roof*' )
             IF ( ALLOCATED( t_surf_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_roof_av(m) = t_surf_roof_av(m) + surf%t_roof(nzt_roof,m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_wall_a*' )
             IF ( ALLOCATED( t_surf_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_wall_a_av(m) = t_surf_wall_a_av(m) + surf%t_wall_a(nzt_wall,m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_wall_b*' )
             IF ( ALLOCATED( t_surf_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_wall_b_av(m) = t_surf_wall_b_av(m) + surf%t_wall_b(nzt_wall,m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_win_a*' )
             IF ( ALLOCATED( t_surf_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_win_a_av(m) = t_surf_win_a_av(m) + surf%t_win_a(nzt_win,m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_win_b*' )
             IF ( ALLOCATED( t_surf_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_win_b_av(m) = t_surf_win_b_av(m) + surf%t_win_b(nzt_win,m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_c_urb*' )
             IF ( ALLOCATED( t_c_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_c_urb_av(m) = t_c_urb_av(m) + surf%t_c_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_h_urb*' )
             IF ( ALLOCATED( t_h_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_h_urb_av(m) = t_h_urb_av(m) + surf%t_h_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_2m_urb*' )
             IF ( ALLOCATED( t_2m_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_2m_urb_av(m) = t_2m_urb_av(m) + surf%t_2m_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_canyon*' )
             IF ( ALLOCATED( pt_can_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_can_av(m) = pt_can_av(m) + surf%pt_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_road*' )
             IF ( ALLOCATED( pt_road_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_road_av(m) = pt_road_av(m) + surf%pt_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_roof*' )
             IF ( ALLOCATED( pt_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_roof_av(m) = pt_roof_av(m) + surf%pt_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_wall_a*' )
             IF ( ALLOCATED( pt_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_wall_a_av(m) = pt_wall_a_av(m) + surf%pt_wall_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_wall_b*' )
             IF ( ALLOCATED( pt_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_wall_b_av(m) = pt_wall_b_av(m) + surf%pt_wall_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_win_a*' )
             IF ( ALLOCATED( pt_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_win_a_av(m) = pt_win_a_av(m) + surf%pt_win_a(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_win_b*' )
             IF ( ALLOCATED( pt_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_win_b_av(m) = pt_win_b_av(m) + surf%pt_win_b(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_thetav_canyon*' )
             IF ( ALLOCATED( vpt_can_av ) )  THEN
                DO  m = 1, surf%ns
                   vpt_can_av(m) = vpt_can_av(m) + surf%vpt_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_thetav_road*' )
             IF ( ALLOCATED( vpt_road_av ) )  THEN
                DO  m = 1, surf%ns
                   vpt_road_av(m) = vpt_road_av(m) + surf%vpt_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_thetav_roof*' )
             IF ( ALLOCATED( vpt_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   vpt_roof_av(m) = vpt_roof_av(m) + surf%vpt_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_wspeed_canyon*' )
             IF ( ALLOCATED( uv_abs_can_av ) )  THEN
                DO  m = 1, surf%ns
                   uv_abs_can_av(m) = uv_abs_can_av(m) + surf%uv_abs_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_wspeed_eff_canyon*' )
             IF ( ALLOCATED( uv_eff_can_av ) )  THEN
                DO  m = 1, surf%ns
                   uv_eff_can_av(m) = uv_eff_can_av(m) + surf%uv_eff_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_us_canyon*' )
             IF ( ALLOCATED( us_can_av ) )  THEN
                DO  m = 1, surf%ns
                   us_can_av(m) = us_can_av(m) + surf%us_can(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_us_road*' )
             IF ( ALLOCATED( us_road_av ) )  THEN
                DO  m = 1, surf%ns
                   us_road_av(m) = us_road_av(m) + surf%us_road(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_us_roof*' )
             IF ( ALLOCATED( us_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   us_roof_av(m) = us_roof_av(m) + surf%us_roof(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_us_urb*' )
             IF ( ALLOCATED( us_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   us_urb_av(m) = us_urb_av(m) + surf%us_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_usws_urb*' )
             IF ( ALLOCATED( usws_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   usws_urb_av(m) = usws_urb_av(m) + surf%usws_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_vsws_urb*' )
             IF ( ALLOCATED( vsws_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   vsws_urb_av(m) = vsws_urb_av(m) + surf%vsws_urb(m)
                ENDDO
             ENDIF

          CASE ( 'slurb_t_road' )
             IF ( ALLOCATED( t_road_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_road, nzb_road
                      t_road_av(k,m) = t_road_av(k,m) + surf%t_road(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_roof' )
             IF ( ALLOCATED( t_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_roof, nzb_roof
                      t_roof_av(k,m) = t_roof_av(k,m) + surf%t_roof(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_wall_a' )
             IF ( ALLOCATED( t_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_wall, nzb_wall
                      t_wall_a_av(k,m) = t_wall_a_av(k,m) + surf%t_wall_a(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_wall_b' )
             IF ( ALLOCATED( t_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_wall, nzb_wall
                      t_wall_b_av(k,m) = t_wall_b_av(k,m) + surf%t_wall_b(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_win_a' )
             IF ( ALLOCATED( t_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_win, nzb_win
                      t_win_a_av(k,m) = t_win_a_av(k,m) + surf%t_win_a(k,m)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_win_b' )
             IF ( ALLOCATED( t_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_win, nzb_win
                      t_win_b_av(k,m) = t_win_b_av(k,m) + surf%t_win_b(k,m)
                   ENDDO
                ENDDO
             ENDIF

       CASE DEFAULT
!
!--       In case of missing or incorrect SLUrb variable, give a meaningful error.
          IF ( variable(1:6) == 'slurb_' )  THEN
             message_string = 'Unknown temporally averaged SLUrb output ' //                       &
                              TRIM ( variable ) // ' requested.'
             CALL message( 'slurb_3d_data_averaging', 'SLU0036', 1, 2, 0, 6, 0 )
          ENDIF

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )
                    CASE ( 'slurb_albedo_urb*' )
             IF ( ALLOCATED( albedo_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   albedo_urb_av(m) = albedo_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_c_liq_road*' )
             IF ( ALLOCATED( c_liq_road_av ) )  THEN
                DO  m = 1, surf%ns
                   c_liq_road_av(m) = c_liq_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_c_liq_roof*' )
             IF ( ALLOCATED( c_liq_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   c_liq_roof_av(m) = c_liq_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_emiss_urb*' )
             IF ( ALLOCATED( emiss_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   emiss_urb_av(m) = emiss_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_road*' )
             IF ( ALLOCATED( ghf_road_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_road_av(m) = ghf_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_roof*' )
             IF ( ALLOCATED( ghf_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_roof_av(m) = ghf_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_wall_a*' )
             IF ( ALLOCATED( ghf_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_wall_a_av(m) = ghf_wall_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_wall_b*' )
             IF ( ALLOCATED( ghf_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_wall_b_av(m) = ghf_wall_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_win_a*' )
             IF ( ALLOCATED( ghf_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_win_a_av(m) = ghf_win_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ghf_win_b*' )
             IF ( ALLOCATED( ghf_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   ghf_win_b_av(m) = ghf_win_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_m_liq_road*' )
             IF ( ALLOCATED( m_liq_road_av ) )  THEN
                DO  m = 1, surf%ns
                   m_liq_road_av(m) = m_liq_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_m_liq_roof*' )
             IF ( ALLOCATED( m_liq_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   m_liq_roof_av(m) = m_liq_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_canyon*' )
             IF ( ALLOCATED( ol_can_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_can_av(m) = ol_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_road*' )
             IF ( ALLOCATED( ol_road_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_road_av(m) = ol_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_roof*' )
             IF ( ALLOCATED( ol_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_roof_av(m) = ol_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ol_urb*' )
             IF ( ALLOCATED( ol_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   ol_urb_av(m) = ol_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_q_canyon*' )
             IF ( ALLOCATED( q_can_av ) )  THEN
                DO  m = 1, surf%ns
                   q_can_av(m) = q_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_q_road*' )
             IF ( ALLOCATED( q_road_av ) )  THEN
                DO  m = 1, surf%ns
                   q_road_av(m) = q_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_q_roof*' )
             IF ( ALLOCATED( q_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   q_roof_av(m) = q_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_qs_road*' )
             IF ( ALLOCATED( qs_road_av ) )  THEN
                DO  m = 1, surf%ns
                   qs_road_av(m) = qs_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_qs_roof*' )
             IF ( ALLOCATED( qs_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   qs_roof_av(m) = qs_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_canyon*' )
             IF ( ALLOCATED( qsws_can_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_can_av(m) = qsws_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_road*' )
             IF ( ALLOCATED( qsws_road_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_road_av(m) = qsws_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_roof*' )
             IF ( ALLOCATED( qsws_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_roof_av(m) = qsws_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_lsm*' )
             IF ( ALLOCATED( qsws_lsm_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      qsws_lsm_av(j,i) = qsws_lsm_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_qsws_urb*' )
             IF ( ALLOCATED( qsws_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   qsws_urb_av(m) = qsws_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_road*' )
             IF ( ALLOCATED( rad_lw_net_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_road_av(m) = rad_lw_net_road_av(m)                                   &
                                           / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_roof*' )
             IF ( ALLOCATED( rad_lw_net_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_roof_av(m) = rad_lw_net_roof_av(m)                                   &
                                           / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_urb*' )
             IF ( ALLOCATED( rad_lw_net_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_urb_av(m) = rad_lw_net_urb_av(m)                                     &
                                          / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_wall_a*' )
             IF ( ALLOCATED( rad_lw_net_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_wall_a_av(m) = rad_lw_net_wall_a_av(m)                               &
                                             / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_wall_b*' )
             IF ( ALLOCATED( rad_lw_net_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_wall_b_av(m) = rad_lw_net_wall_b_av(m)                               &
                                             / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_win_a*' )
             IF ( ALLOCATED( rad_lw_net_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_win_a_av(m) = rad_lw_net_win_a_av(m)                                 &
                                            / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_lw_net_win_b*' )
             IF ( ALLOCATED( rad_lw_net_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_lw_net_win_b_av(m) = rad_lw_net_win_b_av(m)                                 &
                                            / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_road*' )
             IF ( ALLOCATED( rad_sw_net_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_road_av(m) = rad_sw_net_road_av(m)                                   &
                                           / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_roof*' )
             IF ( ALLOCATED( rad_sw_net_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_roof_av(m) = rad_sw_net_roof_av(m)                                   &
                                           / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_urb*' )
             IF ( ALLOCATED( rad_sw_net_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_urb_av(m) = rad_sw_net_urb_av(m)                                     &
                                          / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_wall_a*' )
             IF ( ALLOCATED( rad_sw_net_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_wall_a_av(m) = rad_sw_net_wall_a_av(m)                               &
                                             / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_wall_b*' )
             IF ( ALLOCATED( rad_sw_net_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_wall_b_av(m) = rad_sw_net_wall_b_av(m)                               &
                                             / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_win_a*' )
             IF ( ALLOCATED( rad_sw_net_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_win_a_av(m) = rad_sw_net_win_a_av(m)                                 &
                                            / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_net_win_b*' )
             IF ( ALLOCATED( rad_sw_net_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_net_win_b_av(m) = rad_sw_net_win_b_av(m)                                 &
                                            / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_tr_win_a*' )
             IF ( ALLOCATED( rad_sw_tr_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_tr_win_a_av(m) = rad_sw_tr_win_a_av(m)                                   &
                                           / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rad_sw_tr_win_b*' )
             IF ( ALLOCATED( rad_sw_tr_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rad_sw_tr_win_b_av(m) = rad_sw_tr_win_b_av(m)                                   &
                                           / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_canyon*' )
             IF ( ALLOCATED( rah_can_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_can_av(m) = rah_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_road*' )
             IF ( ALLOCATED( rah_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_road_av(m) = rah_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_roof*' )
             IF ( ALLOCATED( rah_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_roof_av(m) = rah_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_wall_a*' )
             IF ( ALLOCATED( rah_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_wall_a_av(m) = rah_wall_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_wall_b*' )
             IF ( ALLOCATED( rah_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_wall_b_av(m) = rah_wall_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_win_a*' )
             IF ( ALLOCATED( rah_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_win_a_av(m) = rah_win_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_win_b*' )
             IF ( ALLOCATED( rah_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_win_b_av(m) = rah_win_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rah_facade*' )
             IF ( ALLOCATED( rah_facade_av ) )  THEN
                DO  m = 1, surf%ns
                   rah_facade_av(m) = rah_facade_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_ram_urb*' )
             IF ( ALLOCATED( ram_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   ram_urb_av(m) = ram_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rib_canyon*' )
             IF ( ALLOCATED( rib_can_av ) )  THEN
                DO  m = 1, surf%ns
                   rib_can_av(m) = rib_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rib_road*' )
             IF ( ALLOCATED( rib_road_av ) )  THEN
                DO  m = 1, surf%ns
                   rib_road_av(m) = rib_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_rib_roof*' )
             IF ( ALLOCATED( rib_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   rib_roof_av(m) = rib_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_canyon*' )
             IF ( ALLOCATED( shf_can_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_can_av(m) = shf_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_external*' )
             IF ( ALLOCATED( shf_can_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_external_av(m) = shf_external_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_road*' )
             IF ( ALLOCATED( shf_road_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_road_av(m) = shf_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_roof*' )
             IF ( ALLOCATED( shf_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_roof_av(m) = shf_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_traffic*' )
             IF ( ALLOCATED( shf_traffic_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_traffic_av(m) = shf_traffic_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_lsm*' )
             IF ( ALLOCATED( shf_lsm_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      shf_lsm_av(j,i) = shf_lsm_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_urb*' )
             IF ( ALLOCATED( shf_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_urb_av(m) = shf_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_wall_a*' )
             IF ( ALLOCATED( shf_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_wall_a_av(m) = shf_wall_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_wall_b*' )
             IF ( ALLOCATED( shf_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_wall_b_av(m) = shf_wall_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_win_a*' )
             IF ( ALLOCATED( shf_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_win_a_av(m) = shf_win_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_shf_win_b*' )
             IF ( ALLOCATED( shf_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   shf_win_b_av(m) = shf_win_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_canyon*' )
             IF ( ALLOCATED( t_can_av ) )  THEN
                DO  m = 1, surf%ns
                   t_can_av(m) = t_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_rad_urb*' )
             IF ( ALLOCATED( t_rad_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_rad_urb_av(m) = t_rad_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_road*' )
             IF ( ALLOCATED( t_surf_road_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_road_av(m) = t_surf_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_roof*' )
             IF ( ALLOCATED( t_surf_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_roof_av(m) = t_surf_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_wall_a*' )
             IF ( ALLOCATED( t_surf_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_wall_a_av(m) = t_surf_wall_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_wall_b*' )
             IF ( ALLOCATED( t_surf_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_wall_b_av(m) = t_surf_wall_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_win_a*' )
             IF ( ALLOCATED( t_surf_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_win_a_av(m) = t_surf_win_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_surf_win_b*' )
             IF ( ALLOCATED( t_surf_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   t_surf_win_b_av(m) = t_surf_win_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_c_urb*' )
             IF ( ALLOCATED( t_c_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_c_urb_av(m) = t_c_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_h_urb*' )
             IF ( ALLOCATED( t_h_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_h_urb_av(m) = t_h_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_2m_urb*' )
             IF ( ALLOCATED( t_2m_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   t_2m_urb_av(m) = t_2m_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_canyon*' )
             IF ( ALLOCATED( pt_can_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_can_av(m) = pt_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_road*' )
             IF ( ALLOCATED( pt_road_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_road_av(m) = pt_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_roof*' )
             IF ( ALLOCATED( pt_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_roof_av(m) = pt_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_wall_a*' )
             IF ( ALLOCATED( pt_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_wall_a_av(m) = pt_wall_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_wall_b*' )
             IF ( ALLOCATED( pt_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_wall_b_av(m) = pt_wall_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_win_a*' )
             IF ( ALLOCATED( pt_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_win_a_av(m) = pt_win_a_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_theta_win_b*' )
             IF ( ALLOCATED( pt_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   pt_win_b_av(m) = pt_win_b_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_thetav_canyon*' )
             IF ( ALLOCATED( vpt_can_av ) )  THEN
                DO  m = 1, surf%ns
                   vpt_can_av(m) = vpt_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_thetav_road*' )
             IF ( ALLOCATED( vpt_road_av ) )  THEN
                DO  m = 1, surf%ns
                   vpt_road_av(m) = vpt_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_thetav_roof*' )
             IF ( ALLOCATED( vpt_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   vpt_roof_av(m) = vpt_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_wspeed_canyon*' )
             IF ( ALLOCATED( uv_abs_can_av ) )  THEN
                DO  m = 1, surf%ns
                   uv_abs_can_av(m) = uv_abs_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_wspeed_eff_canyon*' )
             IF ( ALLOCATED( uv_eff_can_av ) )  THEN
                DO  m = 1, surf%ns
                   uv_eff_can_av(m) = uv_eff_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_us_canyon*' )
             IF ( ALLOCATED( us_can_av ) )  THEN
                DO  m = 1, surf%ns
                   us_can_av(m) = us_can_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_us_road*' )
             IF ( ALLOCATED( us_road_av ) )  THEN
                DO  m = 1, surf%ns
                   us_road_av(m) = us_road_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_us_roof*' )
             IF ( ALLOCATED( us_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   us_roof_av(m) = us_roof_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_us_urb*' )
             IF ( ALLOCATED( us_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   us_urb_av(m) = us_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_usws_urb*' )
             IF ( ALLOCATED( usws_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   usws_urb_av(m) = usws_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_vsws_urb*' )
             IF ( ALLOCATED( vsws_urb_av ) )  THEN
                DO  m = 1, surf%ns
                   vsws_urb_av(m) = vsws_urb_av(m) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDIF

          CASE ( 'slurb_t_road' )
             IF ( ALLOCATED( t_road_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_road, nzb_road
                      t_road_av(k,m) = t_road_av(k,m) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_roof' )
             IF ( ALLOCATED( t_roof_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_roof, nzb_roof
                      t_roof_av(k,m) = t_roof_av(k,m) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_wall_a' )
             IF ( ALLOCATED( t_wall_a_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_wall, nzb_wall
                      t_wall_a_av(k,m) = t_wall_a_av(k,m) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_wall_b' )
             IF ( ALLOCATED( t_wall_b_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_wall, nzb_wall
                      t_wall_b_av(k,m) = t_wall_b_av(k,m) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_win_a' )
             IF ( ALLOCATED( t_win_a_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_win, nzb_win
                      t_win_a_av(k,m) = t_win_a_av(k,m) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'slurb_t_win_b' )
             IF ( ALLOCATED( t_win_b_av ) )  THEN
                DO  m = 1, surf%ns
                   DO  k = nzt_win, nzb_win
                      t_win_b_av(k,m) = t_win_b_av(k,m) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

       CASE DEFAULT
!
!--       In case of missing or incorrect SLUrb variable, give a meaningful error.
          IF ( variable(1:6) == 'slurb_' )  THEN
             message_string = 'Unknown temporally averaged SLUrb output ' //                       &
                              TRIM ( variable ) // ' requested.'
             CALL message( 'slurb_3d_data_averaging', 'SLU2001', 1, 2, 0, 6, 0 )
          ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE slurb_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Aggregate the urban+lsm fluxes which are used in surface-atmosphere coupling.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE slurb_atmospheric_model_coupler

    INTEGER(iwp) ::  i   !< loop index (x-direction)
    INTEGER(iwp) ::  j   !< loop index (y-direction)
    INTEGER(iwp) ::  m   !< running LSM tile index
    INTEGER(iwp) ::  mm  !< running SLUrb tile index


!
!-- The aggregated values for shf and qsws are stored in special arrays to prevent interference
!-- to MOST calculations in surface_layer_fluxes_mod. In the time integration scheme, the
!-- surface_layer_fluxes, which computes the scaling parameters for the LSM surfaces, is called
!-- before LSM itself. Computation of the scaling parameters depends on shf and qsws, thus modifying
!-- these directly here would lead to SLUrb modifying the LSM scheme. Therefore, additional arrays
!-- to store the aggregated values (shf_agg, qsws_agg) are introduced. These are used for surface-
!-- atmosphere coupling by the diffusion_s routine. LSM and surface_layer_fluxes continue to use
!-- the unaggregated values. For momentum fluxes these are not needed as surface_layer_fluxes
!-- has no dependency on them. The aggregation is done for upward surfaces only, as SLUrb doesn't
!-- have implementation for vertical ones. This is to be changed with cut-cell topography at some
!-- point in the future.
    DO  m = 1, surf_lsm%ns
       i = surf_lsm%i(m)
       j = surf_lsm%j(m)
       mm = surf%m(j,i)
!
!--    Jump to next surface element, if not a SLUrb cell.
       IF ( .NOT. ( mm > 0 ) )  CYCLE

       surf_lsm%usws(m) = MERGE( ( 1.0 - fr_urb(j,i) ) * surf_lsm%usws(m)                          &
                                 + fr_urb(j,i) * surf%usws_urb(mm),                                &
                                   surf_lsm%usws(m), surf_lsm%upward(m) )
       surf_lsm%vsws(m) = MERGE( ( 1.0 - fr_urb(j,i) ) * surf_lsm%vsws(m)                          &
                                 + fr_urb(j,i) * surf%vsws_urb(mm),                                &
                                   surf_lsm%vsws(m), surf_lsm%upward(m) )
       surf_lsm%shf_agg(m) = MERGE( ( 1.0 - fr_urb(j,i) ) * surf_lsm%shf_agg(m)                    &
                                    + fr_urb(j,i) * surf%shf_urb(mm),                              &
                                    surf_lsm%shf_agg(m), surf_lsm%upward(m) )
       IF ( moist_physics )  THEN
          surf_lsm%qsws_agg(m) = MERGE( ( 1.0 - fr_urb(j,i) ) * surf_lsm%qsws_agg(m)               &
                                        + fr_urb(j,i) * surf%qsws_urb(mm),                         &
                                        surf_lsm%qsws_agg(m), surf_lsm%upward(m) )
!
!--    Case for moist physics disabled only for SLUrb, not for the atmospheric simulation
!--    (zero flux).
       ELSEIF ( .NOT. moist_physics  .AND.  humidity )  THEN
          surf_lsm%qsws_agg(m) = MERGE( ( 1.0 - fr_urb(j,i) ) * surf_lsm%qsws_agg(m),              &
                                        surf_lsm%qsws_agg(m), surf_lsm%upward(m) )
       ENDIF
    ENDDO

END SUBROUTINE slurb_atmospheric_model_coupler


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the dynamical conditions (wind speed, pt, q, vpt) in the street canyon.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE slurb_canyon_model

    INTEGER(iwp) ::  i       !< loop index (x-direction)
    INTEGER(iwp) ::  j       !< loop index (y-direction)
    INTEGER(iwp) ::  k_topo  !< k index of topography
    INTEGER(iwp) ::  k_atm   !< k index of the first atmospheric level
    INTEGER(iwp) ::  m       !< running SLUrb tile index

    LOGICAL  ::  runge_l  !< timestep scheme switch for vectorization

    REAL(wp) ::  c          !< total heat capacity of canyon air column per square metre
    REAL(wp) ::  coef_1     !< coefficient A for the prognostic equation
    REAL(wp) ::  coef_2     !< coefficient B for the prognostic equation
    REAL(wp) ::  f_shf      !< factor for the sensible heat flux
    REAL(wp) ::  f_qsws     !< factor for the latent heat flux
    REAL(wp) ::  qsws_surf  !< aggregated latent heat flux from canyon surfaces per unit area
    REAL(wp) ::  shf_surf   !< aggregated sensible heat flux from canyon surfaces per unit area
    REAL(wp) ::  tq_new     !< mixing ratio tendency for the new RK3 time step
    REAL(wp) ::  tt_new     !< temperature tendency for the new RK3 time step
    REAL(wp) ::  vtws       !< buoyancy flux
    REAL(wp) ::  ws         !< free-convection scale


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_canyon_model'
       CALL debug_message( debug_string, 'start' )
    ENDIF

    runge_l = ( timestep_scheme(1:5) == 'runge' )

    DO  m = 1, surf%ns
!
!--    Index offset of surface element point with respect to adjoining atmospheric grid point.
       i = surf%i(m)
       j = surf%j(m)
       k_topo = topo_top_ind(j,i,0)
       k_atm  = topo_top_ind(j,i,0) + 1
       f_shf  = rho_cp / surf%rah_can(m)
!
!--    Consider total air mass column within the street canyon.
       c = rho_cp * surf%h_bld(m)
!
!--    In canyon temperature prognostic equation, we use already computed fluxes from surfaces
!--    in order to ensure consistency and conservation of energy. Thus, only the fluxes between
!--    canyon air and the atmosphere are linearized.
!
!--    Aggregated sensible heat flux from canyon surfaces (per unit area).
       shf_surf = surf%hw_can(m) * ( ( 1.0_wp - surf%f_win(m) ) *                                  &
                                     ( surf%shf_wall_a(m) + surf%shf_wall_b(m) ) +                 &
                                     surf%f_win(m) * ( surf%shf_win_a(m) + surf%shf_win_b(m) )     &
                                   ) + surf%shf_road(m)
!
!--    Aggregated flux doesn't contain c_p yet.
       shf_surf = shf_surf * c_p

!
!--    Coefficients for the prognostic equation of street canyon temperature.
       coef_1 = f_shf * surf%pt1(m) + shf_surf
       coef_2 = f_shf * d_exner(k_topo)

       surf%t_can_p(m) = ( coef_1 * dt_3d * tsc(2) + c * surf%t_can(m) ) /                         &
                         ( c + coef_2 * dt_3d * tsc(2) )

       surf%t_can_p(m) = surf%t_can_p(m) + dt_3d * tsc(3) * surf%tt_can(m)

       tt_new = ( surf%t_can_p(m) - surf%t_can(m) - dt_3d * tsc(3) * surf%tt_can(m) ) /            &
                ( dt_3d  * tsc(2) )
!
!--    Compute the weighted RK3 tendency to be used in next time step.
       IF ( runge_l )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             surf%tt_can(m) = tt_new
          ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
             surf%tt_can(m) = -9.5625_wp * tt_new + 5.3125_wp * surf%tt_can(m)
          ENDIF
       ENDIF

!
!--    Calculate new pt and shf from canyon to atmosphere.
       surf%pt_can(m) = surf%t_can_p(m) * d_exner(k_topo)
       surf%shf_can(m) = -f_shf * ( surf%pt1(m) - surf%pt_can(m) ) / c_p

!
!--    Compute prognostic street canyon mixing ratio.
       IF ( moist_physics )  THEN

          f_qsws = rho_lv / surf%rah_can(m)
!
!--       Here our "latent heat capacity" is the canyon air column total mass.
          c = rho_lv * surf%h_bld(m)
!
!--       Same for the latent heat flux. Currently only the roads, walls are always dry.
!--       This is a placeholder aggregation for street canyon vegetation,
!--       e.g. green walls, low vegetation etc.
          qsws_surf = surf%qsws_road(m)
!
!--       Aggregated flux doesn't contain l_v yet.
          qsws_surf = surf%qsws_road(m) * l_v
!
!--       Compute new prognostic canyon mixing ratio.
          coef_1 = f_qsws * surf%q1(m) + qsws_surf
          coef_2 = f_qsws

          surf%q_can_p(m) = ( coef_1 * dt_3d * tsc(2) + c * surf%q_can(m) ) /                      &
                            ( c + coef_2 * dt_3d * tsc(2) )

          surf%q_can_p(m) = surf%q_can_p(m) + dt_3d * tsc(3) * surf%tq_can(m)
!
!--       Prevent negative mixing ratios due to temporal discretization. This is done before
!--       the computation of tq_new in order to conserve energy.
          IF ( surf%q_can_p(m) < 0.0_wp )  surf%q_can_p(m) = 0.0_wp

          tq_new = ( surf%q_can_p(m) - surf%q_can(m) - dt_3d * tsc(3) * surf%tq_can(m) ) /         &
                   ( dt_3d  * tsc(2) )
!
!--       Compute the weighted RK3 tendency to be used in next time step.
          IF ( runge_l )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                surf%tq_can(m) = tq_new
             ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
                surf%tq_can(m) = -9.5625_wp * tq_new + 5.3125_wp * surf%tq_can(m)
             ENDIF
          ENDIF

          surf%q_can(m) = surf%q_can_p(m)
          surf%vpt_can(m) = surf%pt_can(m) * ( 1.0_wp + 0.61_wp * surf%q_can_p(m) )
          surf%qsws_can(m) = - f_qsws * ( surf%q1(m) - surf%q_can_p(m) ) / l_v

       ENDIF
!
!--    Compute the canyon horizontal wind speed, Eq. (9), Krayenhoff & Voogt (2007).
       surf%uv_abs_can(m) = surf%uv_abs_can_coef(m) * surf%uv_abs1(m)
!
!--    Calculate the canyon effective wind speed taking into account turbulent processes
!--    Lemonsu et al. (2004) Eqs. (2-3).
!
!--    Free convection scale (wstar) at building roof height (for unstable cases)
!--    In case of moist physics, use virtual temperature (buoyancy) flux in free-convection scale.
       IF ( moist_physics )  THEN
          vtws = surf%shf_can(m) + lv_d_cp * surf%qsws_can(m)
       ELSE
          vtws = surf%shf_can(m)
       ENDIF
!
!--    No scaling for stable cases:
       vtws = MERGE( vtws, 0.0_wp, vtws > 0.0_wp )
       ws = ( g / surf%pt_can(m) * surf%z_mo_can(m) * vtws )**( 1.0_wp / 3.0_wp )
!
!--    Canyon effective wind speed taking account both mean and turbulent wind.
       surf%uv_eff_can(m) = SQRT( surf%uv_abs_can(m)**2 + ( surf%us_can(m) + ws )**2 )

    ENDDO

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_canyon_model'
       CALL debug_message( debug_string, 'end' )
    ENDIF

 END SUBROUTINE slurb_canyon_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to check data output for SLurb model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_check_data_output( variable, unit )

    CHARACTER(LEN=*) ::  unit      !< unit of the variable
    CHARACTER(LEN=*) ::  variable  !< variable name

    CHARACTER(LEN=30) ::  var  !< trimmed variable name


    var = TRIM( variable )

    IF ( var(1:6) == 'slurb_'  )  THEN
       IF ( .NOT. slurb )  THEN
          message_string = 'Output of "' // TRIM( var ) // '" requires slurb = .T..'
          CALL message( 'slurb_check_data_output', 'SLU0037', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
       unit = 'illegal'
       RETURN
    ENDIF

!
!-- Prevent output of moist physical quantities if moist phyiscs is not enabled.
    IF ( var(1:9) == 'slurb_qs_'  .OR.  var(1:8) == 'slurb_q_'  .OR.  var(1:13) == 'slurb_thetav_' &
         .OR.  var(1:11) == 'slurb_qsws_'  .OR.  var(1:12) == 'slurb_c_liq_'  .OR.                 &
         var(1:12) == 'slurb_m_liq_' )                                                             &
    THEN
       IF ( .NOT. moist_physics )  THEN
          message_string = 'Output of "' // TRIM( var ) // '" requires moist_physics = .TRUE. in ' &
                           // 'slurb_parameters and humidity = .TRUE. in initialization_parameters.'
          CALL message( 'slurb_check_data_output', 'SLU0038', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Availability of wall and window resistance output depends on the parametrization.
    IF ( ( var(1:16) == 'slurb_rah_facade' )  .AND.  facade_rah_doe )  THEN
       message_string = 'Output of "' // TRIM( var )  // '" requires '                             &
                        // 'facade_resistance_parametrization = "krayenhoff&voogt" or "rowley".'
       CALL message( 'slurb_check_data_output', 'SLU0041', 1, 2, 0, 6, 0 )
    ELSEIF ( ( var(1:14) == 'slurb_rah_wall'  .OR.  var(1:14) == 'slurb_rah_wall')                 &
             .AND.  ( .NOT.  facade_rah_doe ) )                                                    &
    THEN
       message_string = 'Output of "' // TRIM( var )  // '" requires '                             &
                        // 'facade_resistance_parametrization = "doe-2".'
       CALL message( 'slurb_check_data_output', 'SLU0042', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Search for the variable.
    IF ( var(1:13) == 'slurb_albedo_' )  unit = ''
    IF ( var(1:12) == 'slurb_c_liq_'  )  unit = '%'
    IF ( var(1:12) == 'slurb_emiss_'  )  unit = ''
    IF ( var(1:10) == 'slurb_ghf_'    )  unit = 'W/m2'
    IF ( var(1:12) == 'slurb_m_liq_'  )  unit = 'm'
    IF ( var(1:9)  == 'slurb_ol_'     )  unit = 'm'
    IF ( var(1:8)  == 'slurb_q_'      )  unit = 'kg/kg'
    IF ( var(1:9)  == 'slurb_qs_'     )  unit = 'kg/kg'
    IF ( var(1:11) == 'slurb_qsws_' )  THEN
       IF ( TRIM( flux_output_mode ) == 'kinematic' )  unit = 'K m/s'
       IF ( TRIM( flux_output_mode ) == 'dynamic'   )  unit = 'W/m2'
    ENDIF
    IF ( var(1:10) == 'slurb_rad_'    )  unit = 'W/m2'
    IF ( var(1:10) == 'slurb_rah_'    )  unit = 's/m'
    IF ( var(1:10) == 'slurb_ram_'    )  unit = 's/m'
    IF ( var(1:10) == 'slurb_rib_'    )  unit = ''

    IF ( var(1:9)  == 'slurb_shf' )  THEN
       IF ( TRIM( flux_output_mode ) == 'kinematic' )  unit = 'K m/s'
       IF ( TRIM( flux_output_mode ) == 'dynamic'   )  unit = 'W/m2'
    ENDIF
    IF ( var(1:8)  == 'slurb_t_'      )  unit = 'K'
    IF ( var(1:12) == 'slurb_theta_'  )  unit = 'K'
    IF ( var(1:13) == 'slurb_thetav_' )  unit = 'K'
    IF ( var(1:9)  == 'slurb_us_'     )  unit = 'm/s'
    IF ( var(1:11) == 'slurb_usws_'   )  unit = 'm2/s2'
    IF ( var(1:11) == 'slurb_vsws_'   )  unit = 'm2/s2'
    IF ( var(1:13) == 'slurb_wspeed_' )  unit = 'm/s'
!
!-- Set flags to enable on-demand computations for output statistics.
    IF ( var(1:10) == 'slurb_t_h_' ) calc_t_h = .TRUE.
    IF ( var(1:11) == 'slurb_t_2m_' ) calc_t_2m = .TRUE.
    IF ( var(1:10) == 'slurb_t_c_' ) calc_t_c = .TRUE.

 END SUBROUTINE slurb_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the user-set parameters and general configuration is valid and SLUrb-compatible.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  constant_diffusion,                                                                 &
               neutral,                                                                            &
               plant_canopy,                                                                       &
               urban_surface,                                                                      &
               use_fixed_time

    USE radiation_model_mod,                                                                       &
        ONLY:  sun_direction


!
!-- Initial checks of the control parameters.
    IF ( neutral )  THEN
       message_string = 'SLUrb requires does not support strictly neutral flow.'
       CALL message( 'slurb_check_parameters', 'SLU0000', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( .NOT. radiation )  THEN
       message_string = 'SLUrb requires the radiation model to be switched on.'
       CALL message( 'slurb_check_parameters', 'SLU0001', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Currently this error is not even possible as the combination is possible only with DCEP.
!-- However, in the future a switch to enable average radiation separately might get implemented.
    IF ( average_radiation  .AND.  .NOT. radiation_interactions )  THEN
       message_string = 'SLUrb does not support average radiation if radiation interactions are '  &
                        // 'disabled ( radiation_interactions_on = .F. ).'
       CALL message( 'slurb_check_parameters', 'SLU0002', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( constant_diffusion )  THEN
       message_string = 'SLUrb requires constant_diffusion = .F.'
       CALL message( 'slurb_check_parameters', 'SLU0003', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( plant_canopy )  THEN
       message_string = 'Using 3D plant canopy module together with SLUrb is not supported.'
       CALL message( 'slurb_check_parameters', 'SLU0004', 1, 2, 0, 6, 0 )
    ENDIF


    IF ( radiation_scheme == 'constant'  .AND.  .NOT. use_fixed_time )  THEN
       message_string = 'Using radiation_scheme = "constant" with the SLUrb model requires '       &
                        // 'the setting of use_fixed_time = .T..'
       CALL message( 'slurb_check_parameters', 'SLU0005', 1, 2, 0, 6, 0 )
    ELSE
       sun_direction = .TRUE.
    ENDIF

!
!-- Ensure that USM is not enabled.
    IF ( urban_surface )  THEN
       message_string = 'Enabling the urban surface model for resolved buildings (USM) is not '    &
                     // 'allowed together with the single-layer urban model (SLUrb) enabled in '   &
                     // 'the same simulation domain or nest.'
       CALL message( 'slurb_check_parameters', 'SLU0006', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Turn off moist physical processes if humidity is off (default is on if humidity is enabled).
    IF ( moist_physics  .AND.  .NOT. humidity )  THEN
       moist_physics = .FALSE.
    ENDIF

    IF ( .NOT. moist_physics  .AND.  precipitation)  THEN
       message_string = 'Using the SLUrb model with precipitation enabled requires '               &
                        // 'the setting of moist_physics = .T..'
       CALL message( 'slurb_check_parameters', 'SLU0007', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set parametrization for roughness length for heat for horizontal surfaces (only for MOST).
    IF ( TRIM( aero_roughness_heat ) == 'kanda' )  THEN
       roughness_kanda = .TRUE.
    ELSEIF ( TRIM( aero_roughness_heat ) == 'fixed' )  THEN
       roughness_kanda = .FALSE.
    ELSE
       message_string = 'Invalid setting for the parametrization of the aerodynamic roughness '    &
                        // 'for heat aero_roughness_heat = ' // TRIM( aero_roughness_heat ) // '.'
       CALL message( 'slurb_check_parameters', 'SLU0008', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set parametrization for canyon wind speed.
    IF ( TRIM( street_canyon_wspeed_factor ) == 'krayenhoff&voogt' )  THEN
       uv_can_factor_kray = .TRUE.
    ELSEIF ( TRIM( street_canyon_wspeed_factor ) == 'masson' )  THEN
       uv_can_factor_masson = .TRUE.
    ELSEIF ( TRIM( street_canyon_wspeed_factor ) == 'surfex' )  THEN
       uv_can_factor_surfex = .TRUE.
    ELSE
       message_string = 'Invalid setting for the parametrization of the street canyon wind speed ' &
                        // 'street_canyon_wspeed_factor = ' // TRIM( street_canyon_wspeed_factor ) &
                        // '.'
       CALL message( 'slurb_check_parameters', 'SLU0009', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( facade_resistance_parametrization ) == 'doe-2' )  THEN
       facade_rah_doe = .TRUE.
    ELSEIF ( TRIM( facade_resistance_parametrization ) == 'krayenhoff&voogt' )  THEN
       facade_rah_kray = .TRUE.
    ELSEIF ( TRIM( facade_resistance_parametrization ) == 'rowley' )  THEN
       facade_rah_rowley = .TRUE.
    ELSE
       message_string = 'Invalid setting for the parametrization of aerodynamic resistance for '   &
                        // 'facades facade_resistance_parametrization = '                          &
                        // TRIM( facade_resistance_parametrization ) // '.'
       CALL message( 'slurb_check_parameters', 'SLU0040', 1, 2, 0, 6, 0 )
    ENDIF



!
!-- Check if minimum number of material layers is set for each facet.
    IF ( n_layers_roads < 3 )  THEN
       WRITE( message_string, * ) 'At least three layers are required to model road surfaces '     &
                       // 'in SLUrb. Current setting is n_layers_roads = ', n_layers_roads, '.'
       CALL message( 'slurb_check_parameters', 'SLU0010', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( n_layers_roofs < 3 )  THEN
       WRITE( message_string, * ) 'At least three layers are required to model roof surfaces '     &
                       // 'in SLUrb. Current setting is n_layers_roofs = ', n_layers_roofs, '.'
       CALL message( 'slurb_check_parameters', 'SLU0011', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( n_layers_walls < 3 )  THEN
       WRITE( message_string, * ) 'At least three layers are required to model wall surfaces '     &
                       // 'in SLUrb. Current setting is n_layers_walls = ', n_layers_walls, '.'
       CALL message( 'slurb_check_parameters', 'SLU0012', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( n_layers_windows < 3 )  THEN
       WRITE( message_string, * ) 'At least three layers are required to model window surfaces '   &
                       // 'in SLUrb. Current setting is n_layers_windows = ', n_layers_windows, '.'
       CALL message( 'slurb_check_parameters', 'SLU0013', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE slurb_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 2D output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_data_output_2d( av, variable, found, grid_type, mode, local_pf, two_d, nzb_do,   &
                                  nzt_do )

    CHARACTER(LEN=*), INTENT(IN)  ::  variable  !< variable name
    CHARACTER(LEN=*), INTENT(IN)  ::  mode      !< cross-section direction, always 'xy' for SLUrb
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_type   !< Grid type (always "zu1" for biom)

    INTEGER(iwp), INTENT(IN) ::  av      !< data averaging flag: 0 = no, 1 = yes
    INTEGER(iwp), INTENT(IN) ::  nzb_do  !< vertical output index (bottom)
    INTEGER(iwp), INTENT(IN) ::  nzt_do  !< = nzb_do always for surface output

    LOGICAL, INTENT(OUT) ::  found  !< flag to indicate the variable is a SLUrb output
    LOGICAL, INTENT(OUT) ::  two_d  !< flag to indicate the variable is a 2D variable (always true here)

    INTEGER(iwp) ::  i   !< loop index (x-direction)
    INTEGER(iwp) ::  j   !< loop index (y-direction)
    INTEGER(iwp) ::  k   !< layer index for heat/waterflux conversion
    INTEGER(iwp) ::  m   !< running SLUrb tile index

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< result grid to return


    found = .TRUE.

    IF ( variable(1:6) /= 'slurb_' )  THEN
       found = .FALSE.
       RETURN
    ENDIF

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'slurb_albedo_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%albedo_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( albedo_urb_av ) )  THEN
                ALLOCATE( albedo_urb_av(1:surf%ns) )
                albedo_urb_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = albedo_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_c_liq_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%c_liq_road(m) * 100.0_wp
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( c_liq_road_av ) )  THEN
                ALLOCATE( c_liq_road_av(1:surf%ns) )
                c_liq_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = c_liq_road_av(m) * 100.0_wp
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_c_liq_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%c_liq_roof(m) * 100.0_wp
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( c_liq_roof_av ) )  THEN
                ALLOCATE( c_liq_roof_av(1:surf%ns) )
                c_liq_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = c_liq_roof_av(m) * 100.0_wp
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_emiss_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%emiss_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( emiss_urb_av ) )  THEN
                ALLOCATE( emiss_urb_av(1:surf%ns) )
                emiss_urb_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = emiss_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ghf_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ghf_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ghf_road_av ) )  THEN
                ALLOCATE( ghf_road_av(1:surf%ns) )
                ghf_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ghf_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ghf_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ghf_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ghf_roof_av ) )  THEN
                ALLOCATE( ghf_roof_av(1:surf%ns) )
                ghf_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ghf_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ghf_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ghf_wall_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ghf_wall_a_av ) )  THEN
                ALLOCATE( ghf_wall_a_av(1:surf%ns) )
                ghf_wall_a_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ghf_wall_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ghf_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ghf_wall_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ghf_wall_b_av ) )  THEN
                ALLOCATE( ghf_wall_b_av(1:surf%ns) )
                ghf_wall_b_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ghf_wall_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ghf_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ghf_win_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ghf_win_a_av ) )  THEN
                ALLOCATE( ghf_win_a_av(1:surf%ns) )
                ghf_win_a_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ghf_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ghf_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ghf_win_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ghf_win_b_av ) )  THEN
                ALLOCATE( ghf_win_b_av(1:surf%ns) )
                ghf_win_b_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ghf_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_m_liq_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%m_liq_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( m_liq_road_av ) )  THEN
                ALLOCATE( m_liq_road_av(1:surf%ns) )
                m_liq_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = m_liq_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_m_liq_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%m_liq_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( m_liq_roof_av ) )  THEN
                ALLOCATE( m_liq_roof_av(1:surf%ns) )
                m_liq_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = m_liq_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ol_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ol_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ol_can_av ) )  THEN
                ALLOCATE( ol_can_av(1:surf%ns) )
                ol_can_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ol_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ol_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ol_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ol_road_av ) )  THEN
                ALLOCATE( ol_road_av(1:surf%ns) )
                ol_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ol_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ol_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ol_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ol_roof_av ) )  THEN
                ALLOCATE( ol_roof_av(1:surf%ns) )
                ol_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ol_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ol_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ol_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ol_urb_av ) )  THEN
                ALLOCATE( ol_urb_av(1:surf%ns) )
                ol_urb_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ol_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_q_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%q_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( q_can_av ) )  THEN
                ALLOCATE( q_can_av(1:surf%ns) )
                q_can_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = q_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_q_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%q_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( q_road_av ) )  THEN
                ALLOCATE( q_road_av(1:surf%ns) )
                q_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = q_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_q_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%q_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( q_roof_av ) )  THEN
                ALLOCATE( q_roof_av(1:surf%ns) )
                q_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = q_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qs_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%qs_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qs_road_av ) )  THEN
                ALLOCATE( qs_road_av(1:surf%ns) )
                qs_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = qs_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qs_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%qs_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qs_roof_av ) )  THEN
                ALLOCATE( qs_roof_av(1:surf%ns) )
                qs_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = qs_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qsws_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%qsws_can(m) * waterflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qsws_can_av ) )  THEN
                ALLOCATE( qsws_can_av(1:surf%ns) )
                qsws_can_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = qsws_can_av(m) * waterflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qsws_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%qsws_road(m) * waterflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qsws_road_av ) )  THEN
                ALLOCATE( qsws_road_av(1:surf%ns) )
                qsws_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = qsws_road_av(m) * waterflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qsws_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%qsws_roof(m) * waterflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qsws_roof_av ) )  THEN
                ALLOCATE( qsws_roof_av(1:surf%ns) )
                qsws_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = qsws_roof_av(m) * waterflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qsws_lsm*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm%ns
                i = surf_lsm%i(m)
                j = surf_lsm%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = local_pf(i,j,nzb+1)                                          &
                                      + MERGE( surf_lsm%qsws(m), 0.0_wp, surf_lsm%upward(m) )      &
                                      * waterflux_output_conversion(k)
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qsws_lsm_av ) )  THEN
                ALLOCATE( qsws_lsm_av(nys:nyn,nxl:nxr) )
                qsws_lsm_av(:,:) = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   k = topo_top_ind(j,i,0)
                   local_pf(i,j,nzb+1) = qsws_lsm_av(j,i) * waterflux_output_conversion(k)
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_qsws_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%qsws_urb(m) * waterflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qsws_urb_av ) )  THEN
                ALLOCATE( qsws_urb_av(1:surf%ns) )
                qsws_urb_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = qsws_urb_av(m) * waterflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_road_av ) )  THEN
                ALLOCATE( rad_lw_net_road_av(1:surf%ns) )
                rad_lw_net_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_roof_av ) )  THEN
                ALLOCATE( rad_lw_net_roof_av(1:surf%ns) )
                rad_lw_net_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_urb_av ) )  THEN
                ALLOCATE( rad_lw_net_urb_av(1:surf%ns) )
                rad_lw_net_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_wall_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_wall_a_av ) )  THEN
                ALLOCATE( rad_lw_net_wall_a_av(1:surf%ns) )
                rad_lw_net_wall_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_wall_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_wall_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_wall_b_av ) )  THEN
                ALLOCATE( rad_lw_net_wall_b_av(1:surf%ns) )
                rad_lw_net_wall_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_wall_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_win_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_win_a_av ) )  THEN
                ALLOCATE( rad_lw_net_win_a_av(1:surf%ns) )
                rad_lw_net_win_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_lw_net_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_lw_net_win_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_net_win_b_av ) )  THEN
                ALLOCATE( rad_lw_net_win_b_av(1:surf%ns) )
                rad_lw_net_win_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_lw_net_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_road_av ) )  THEN
                ALLOCATE( rad_sw_net_road_av(1:surf%ns) )
                rad_sw_net_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_roof_av ) )  THEN
                ALLOCATE( rad_sw_net_roof_av(1:surf%ns) )
                rad_sw_net_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_urb_av ) )  THEN
                ALLOCATE( rad_sw_net_urb_av(1:surf%ns) )
                rad_sw_net_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_wall_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_wall_a_av ) )  THEN
                ALLOCATE( rad_sw_net_wall_a_av(1:surf%ns) )
                rad_sw_net_wall_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_wall_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_wall_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_wall_b_av ) )  THEN
                ALLOCATE( rad_sw_net_wall_b_av(1:surf%ns) )
                rad_sw_net_wall_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_wall_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_win_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_win_a_av ) )  THEN
                ALLOCATE( rad_sw_net_win_a_av(1:surf%ns) )
                rad_sw_net_win_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_net_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_win_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_net_win_b_av ) )  THEN
                ALLOCATE( rad_sw_net_win_b_av(1:surf%ns) )
                rad_sw_net_win_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_net_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_tr_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_win_a(m)                                     &
                                      * ( 1.0_wp - SUM( surf%absorption_win(:,m) ) )
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_tr_win_a_av ) )  THEN
                ALLOCATE( rad_sw_tr_win_a_av(1:surf%ns) )
                rad_sw_tr_win_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_tr_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rad_sw_tr_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rad_sw_net_win_b(m)                                     &
                                      * ( 1.0_wp - SUM( surf%absorption_win(:,m) ) )
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_tr_win_b_av ) )  THEN
                ALLOCATE( rad_sw_tr_win_b_av(1:surf%ns) )
                rad_sw_tr_win_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rad_sw_tr_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_can_av ) )  THEN
                ALLOCATE( rah_can_av(1:surf%ns) )
                rah_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_road_av ) )  THEN
                ALLOCATE( rah_road_av(1:surf%ns) )
                rah_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_roof_av ) )  THEN
                ALLOCATE( rah_roof_av(1:surf%ns) )
                rah_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_wall_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_wall_a_av ) )  THEN
                ALLOCATE( rah_wall_a_av(1:surf%ns) )
                rah_wall_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_wall_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_wall_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_wall_b_av ) )  THEN
                ALLOCATE( rah_wall_b_av(1:surf%ns) )
                rah_wall_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_wall_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_win_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_win_a_av ) )  THEN
                ALLOCATE( rah_win_a_av(1:surf%ns) )
                rah_win_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_win_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_win_b_av ) )  THEN
                ALLOCATE( rah_win_b_av(1:surf%ns) )
                rah_win_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rah_facade*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rah_facade(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rah_facade_av ) )  THEN
                ALLOCATE( rah_facade_av(1:surf%ns) )
                rah_facade_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rah_facade_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_ram_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%ram_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ram_urb_av ) )  THEN
                ALLOCATE( ram_urb_av(1:surf%ns) )
                ram_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = ram_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rib_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rib_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rib_can_av ) )  THEN
                ALLOCATE( rib_can_av(1:surf%ns) )
                rib_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rib_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rib_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rib_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rib_road_av ) )  THEN
                ALLOCATE( rib_road_av(1:surf%ns) )
                rib_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rib_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_rib_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%rib_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rib_roof_av ) )  THEN
                ALLOCATE( rib_roof_av(1:surf%ns) )
                rib_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = rib_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_can(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_can_av ) )  THEN
                ALLOCATE( shf_can_av(1:surf%ns) )
                shf_can_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_can_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_external*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_external(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_external_av ) )  THEN
                ALLOCATE( shf_external_av(1:surf%ns) )
                shf_external_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_external_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_road(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_road_av ) )  THEN
                ALLOCATE( shf_road_av(1:surf%ns) )
                shf_road_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_road_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_roof(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_roof_av ) )  THEN
                ALLOCATE( shf_roof_av(1:surf%ns) )
                shf_roof_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_roof_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_traffic*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_traffic(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_traffic_av ) )  THEN
                ALLOCATE( shf_traffic_av(1:surf%ns) )
                shf_traffic_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_traffic_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_lsm*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf_lsm%ns
                i = surf_lsm%i(m)
                j = surf_lsm%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = local_pf(i,j,nzb+1)                                          &
                                      + MERGE( surf_lsm%shf(m), 0.0_wp, surf_lsm%upward(m) )       &
                                      * heatflux_output_conversion(k)
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_lsm_av ) )  THEN
                ALLOCATE( shf_lsm_av(nys:nyn,nxl:nxr) )
                shf_lsm_av(:,:) = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   k = topo_top_ind(j,i,0)
                   local_pf(i,j,nzb+1) = shf_lsm_av(j,i) * heatflux_output_conversion(k)
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_urb(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_urb_av ) )  THEN
                ALLOCATE( shf_urb_av(1:surf%ns) )
                shf_urb_av(:) = 0.0_wp
            ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_urb_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_wall_a(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_wall_a_av ) )  THEN
                ALLOCATE( shf_wall_a_av(1:surf%ns) )
                shf_wall_a_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_wall_a_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_wall_b(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_wall_b_av ) )  THEN
                ALLOCATE( shf_wall_b_av(1:surf%ns) )
                shf_wall_b_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_wall_b_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_win_a(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_win_a_av ) )  THEN
                ALLOCATE( shf_win_a_av(1:surf%ns) )
                shf_win_a_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_win_a_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_shf_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = surf%shf_win_b(m) * heatflux_output_conversion(k)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( shf_win_b_av ) )  THEN
                ALLOCATE( shf_win_b_av(1:surf%ns) )
                shf_win_b_av(:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = topo_top_ind(j,i,0)
                local_pf(i,j,nzb+1) = shf_win_b_av(m) * heatflux_output_conversion(k)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_can_av ) )  THEN
                ALLOCATE( t_can_av(1:surf%ns) )
                t_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_rad_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_rad_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_rad_urb_av ) )  THEN
                ALLOCATE( t_rad_urb_av(1:surf%ns) )
                t_rad_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_rad_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_surf_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_road(nzt_road,m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_surf_road_av ) )  THEN
                ALLOCATE( t_surf_road_av(1:surf%ns) )
                t_surf_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_surf_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_surf_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_roof(nzt_roof,m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_surf_roof_av ) )  THEN
                ALLOCATE( t_surf_roof_av(1:surf%ns) )
                t_surf_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_surf_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_surf_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_wall_a(nzt_wall,m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_surf_wall_a_av ) )  THEN
                ALLOCATE( t_surf_wall_a_av(1:surf%ns) )
                t_surf_wall_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_surf_wall_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_surf_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_wall_b(nzt_wall,m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_surf_wall_b_av ) )  THEN
                ALLOCATE( t_surf_wall_b_av(1:surf%ns) )
                t_surf_wall_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_surf_wall_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_surf_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_win_a(nzt_win,m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_surf_win_a_av ) )  THEN
                ALLOCATE( t_surf_win_a_av(1:surf%ns) )
                t_surf_win_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_surf_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_surf_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_win_b(nzt_win,m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_surf_win_b_av ) )  THEN
                ALLOCATE( t_surf_win_b_av(1:surf%ns) )
                t_surf_win_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_surf_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_c_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_c_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_c_urb_av ) )  THEN
                ALLOCATE( t_c_urb_av(1:surf%ns) )
                t_c_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_c_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_t_h_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_h_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_h_urb_av ) )  THEN
                ALLOCATE( t_h_urb_av(1:surf%ns) )
                t_h_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_h_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

          CASE ( 'slurb_t_2m_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%t_2m_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_2m_urb_av ) )  THEN
                ALLOCATE( t_2m_urb_av(1:surf%ns) )
                t_2m_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = t_2m_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_can_av ) )  THEN
                ALLOCATE( pt_can_av(1:surf%ns) )
                pt_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_road_av ) )  THEN
                ALLOCATE( pt_road_av(1:surf%ns) )
                pt_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_roof_av ) )  THEN
                ALLOCATE( pt_roof_av(1:surf%ns) )
                pt_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_wall_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_wall_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_wall_a_av ) )  THEN
                ALLOCATE( pt_wall_a_av(1:surf%ns) )
                pt_wall_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_wall_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_wall_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_wall_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_wall_b_av ) )  THEN
                ALLOCATE( pt_wall_b_av(1:surf%ns) )
                pt_wall_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_wall_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_win_a*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_win_a(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_win_a_av ) )  THEN
                ALLOCATE( pt_win_a_av(1:surf%ns) )
                pt_win_a_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_win_a_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_theta_win_b*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%pt_win_b(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_win_b_av ) )  THEN
                ALLOCATE( pt_win_b_av(1:surf%ns) )
                pt_win_b_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = pt_win_b_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_thetav_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%vpt_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vpt_can_av ) )  THEN
                ALLOCATE( vpt_can_av(1:surf%ns) )
                vpt_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = vpt_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_thetav_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%vpt_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vpt_road_av ) )  THEN
                ALLOCATE( vpt_road_av(1:surf%ns) )
                vpt_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = vpt_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_thetav_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%vpt_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vpt_roof_av ) )  THEN
                ALLOCATE( vpt_roof_av(1:surf%ns) )
                vpt_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = vpt_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_wspeed_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%uv_abs_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uv_abs_can_av ) )  THEN
                ALLOCATE( uv_abs_can_av(1:surf%ns) )
                uv_abs_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = uv_abs_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_wspeed_eff_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%uv_eff_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uv_eff_can_av ) )  THEN
                ALLOCATE( uv_eff_can_av(1:surf%ns) )
                uv_eff_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = uv_eff_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_us_canyon*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%us_can(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( us_can_av ) )  THEN
                ALLOCATE( us_can_av(1:surf%ns) )
                us_can_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = us_can_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_us_road*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%us_road(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( us_road_av ) )  THEN
                ALLOCATE( us_road_av(1:surf%ns) )
                us_road_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = us_road_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_us_roof*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%us_roof(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( us_roof_av ) )  THEN
                ALLOCATE( us_roof_av(1:surf%ns) )
                us_roof_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = us_roof_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_us_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%us_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( us_urb_av ) )  THEN
                ALLOCATE( us_urb_av(1:surf%ns) )
                us_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = us_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_usws_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%usws_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( usws_urb_av ) )  THEN
                ALLOCATE( usws_urb_av(1:surf%ns) )
                usws_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = usws_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE ( 'slurb_vsws_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = surf%vsws_urb(m)
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vsws_urb_av ) )  THEN
                ALLOCATE( vsws_urb_av(1:surf%ns) )
                vsws_urb_av(:) = 0.0_wp
             ENDIF
            DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                local_pf(i,j,nzb+1) = vsws_urb_av(m)
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid_type = 'zu1'
          two_d = .TRUE.

       CASE DEFAULT
!
!--       In case of missing or incorrect SLUrb variable, give a meaningful error.
          IF ( variable(1:6) == 'slurb_' )  THEN
             message_string = 'Unknown SLUrb output ' // TRIM ( variable ) // ' requested.'
             CALL message( 'slurb_data_output_2d', 'SLU2000', 1, 2, 0, 6, 0 )
          ENDIF

          found = .FALSE.
          grid_type = 'none'

    END SELECT

 END SUBROUTINE slurb_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables.
!--------------------------------------------------------------------------------------------------!
  SUBROUTINE slurb_data_output_3d( av, variable, found, local_pf, mask_topography, nzb_do, nzt_do, &
                                   fill_value )

    CHARACTER(LEN=*), INTENT(IN) ::  variable !< variable name

    INTEGER(iwp), INTENT(IN)  ::  av      !< data averaging flag: 0 = no, 1 = yes
    INTEGER(iwp), INTENT(OUT) ::  nzb_do  !< start index of the layers
    INTEGER(iwp), INTENT(OUT) ::  nzt_do  !< end index of the layers

    REAL(wp), INTENT(IN) :: fill_value  !< output fill value

    INTEGER(iwp) ::  i  !< loop index (x-direction)
    INTEGER(iwp) ::  j  !< loop index (y-direction)
    INTEGER(iwp) ::  k  !< layer index for heat/waterflux conversion
    INTEGER(iwp) ::  m  !< running SLUrb tile index

    LOGICAL, INTENT(OUT)   ::  found            !< flag to indicate the variable is a SLUrb output
    LOGICAL, INTENT(INOUT) ::  mask_topography  !< flag if topography points shall be masked

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf  !< result grid to return


    found = .TRUE.

    IF ( variable(1:6) /= 'slurb_' )  THEN
       found = .FALSE.
       RETURN
    ENDIF

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'slurb_t_road' )
          mask_topography = .FALSE.
          nzb_do = LBOUND(surf%t_road,1)
          nzt_do = UBOUND(surf%t_road,1)
          ALLOCATE ( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
          IF ( .NOT. data_output_raw )  local_pf(:,:,:) = fill_value
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = surf%t_road(k,m)
                ENDDO
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_road_av ) )  THEN
                ALLOCATE( t_road_av(nzb_do:nzt_do,1:surf%ns) )
                t_road_av(:,:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = t_road_av(k,m)
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'slurb_t_roof' )
          mask_topography = .FALSE.
          nzb_do = LBOUND(surf%t_roof,1)
          nzt_do = UBOUND(surf%t_roof,1)
          ALLOCATE ( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
          IF ( .NOT. data_output_raw )  local_pf(:,:,:) = fill_value
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = surf%t_roof(k,m)
                ENDDO
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_roof_av ) )  THEN
                ALLOCATE( t_roof_av(nzb_do:nzt_do,1:surf%ns) )
                t_roof_av(:,:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = t_roof_av(k,m)
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'slurb_t_wall_a' )
          mask_topography = .FALSE.
          nzb_do = LBOUND(surf%t_wall_a,1)
          nzt_do = UBOUND(surf%t_wall_a,1)
          ALLOCATE ( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
          IF ( .NOT. data_output_raw )  local_pf(:,:,:) = fill_value
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = surf%t_wall_a(k,m)
                ENDDO
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_wall_a_av ) )  THEN
                ALLOCATE( t_wall_a_av(nzb_do:nzt_do,1:surf%ns) )
                t_wall_a_av(:,:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = t_wall_a_av(k,m)
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'slurb_t_wall_b' )
          mask_topography = .FALSE.
          nzb_do = LBOUND(surf%t_wall_b,1)
          nzt_do = UBOUND(surf%t_wall_b,1)
          ALLOCATE ( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
          IF ( .NOT. data_output_raw )  local_pf(:,:,:) = fill_value
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = surf%t_wall_b(k,m)
                ENDDO
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_wall_b_av ) )  THEN
                ALLOCATE( t_wall_b_av(nzb_do:nzt_do,1:surf%ns) )
                t_wall_b_av(:,:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = t_wall_b_av(k,m)
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'slurb_t_win_a' )
          mask_topography = .FALSE.
          nzb_do = LBOUND(surf%t_win_a,1)
          nzt_do = UBOUND(surf%t_win_a,1)
          ALLOCATE ( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
          IF ( .NOT. data_output_raw )  local_pf(:,:,:) = fill_value
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = surf%t_win_a(k,m)
                ENDDO
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_win_a_av ) )  THEN
                ALLOCATE( t_win_a_av(nzb_do:nzt_do,1:surf%ns) )
                t_win_a_av(:,:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = t_win_a_av(k,m)
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'slurb_t_win_b' )
          mask_topography = .FALSE.
          nzb_do = LBOUND(surf%t_win_b,1)
          nzt_do = UBOUND(surf%t_win_b,1)
          ALLOCATE ( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
          IF ( .NOT. data_output_raw )  local_pf(:,:,:) = fill_value
          IF ( av == 0 )  THEN
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = surf%t_win_b(k,m)
                ENDDO
              ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( t_win_b_av ) )  THEN
                ALLOCATE( t_win_b_av(nzb_do:nzt_do,1:surf%ns) )
                t_win_b_av(:,:) = 0.0_wp
             ENDIF
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = t_win_b_av(k,m)
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
!
!--       In case of missing or incorrect SLUrb variable, give a meaningful error.
          IF ( variable(1:6) == 'slurb_' )  THEN
             message_string = 'Unknown SLUrb output ' // TRIM ( variable ) // ' requested.'
             CALL message( 'slurb_data_output_3d', 'SLU2000', 1, 2, 0, 6, 0 )
          ENDIF
          found = .FALSE.

    END SELECT

 END SUBROUTINE slurb_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    CHARACTER(LEN=*), INTENT(IN) ::  variable  !< variable name

    LOGICAL, INTENT(OUT) ::  found  !< flag to indicate the variable is a SLUrb output

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x  !< coordinate variable name for x-grid
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y  !< coordinate variable name for y-grid
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z  !< coordinate varaible name for z-grid

    INTEGER(iwp) ::  ilen  !< length of the variable name string after trimming


    found = .TRUE.

    ilen = LEN_TRIM( variable )

    IF ( variable(1:6) == 'slurb_'  .AND.  variable(ilen-2:ilen) == '_xy' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu1'
    ELSEIF ( variable(1:12) == 'slurb_t_road' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'nroad_3d'
    ELSEIF ( variable(1:12) == 'slurb_t_roof' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'nroof_3d'
    ELSEIF ( variable(1:12) == 'slurb_t_wall' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'nwall_3d'
    ELSEIF ( variable(1:11) == 'slurb_t_win' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'nwin_3d'
    ELSE
          found = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'
    ENDIF

 END SUBROUTINE slurb_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Surface and subsurface energy balance computations of roofs, walls, windows and roads.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_energy_balance_model

    INTEGER(iwp) ::  i       !< loop index (x-direction)
    INTEGER(iwp) ::  j       !< loop index (y-direction)
    INTEGER(iwp) ::  k_topo  !< k index of topography
    INTEGER(iwp) ::  k_atm   !< k index of the first atmospheric level
    INTEGER(iwp) ::  m       !< running SLUrb tile index

    LOGICAL ::  runge_l  !< flag to vectorize timestep scheme switch


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_energy_balance_model'
       CALL debug_message( debug_string, 'start' )
    ENDIF

    runge_l = ( timestep_scheme(1:5) == 'runge' )

    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k_topo = topo_top_ind(j,i,0)
       k_atm = topo_top_ind(j,i,0) + 1
!
!--    Call specific models for all the facets.
       CALL roof_model
       CALL wall_model
       IF ( surf%f_win(m) /= 0.0_wp )  CALL window_model
       CALL road_model
    ENDDO

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_energy_balance_model'
       CALL debug_message( debug_string, 'end' )
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the new surface prognostic temperature for current time step using RK3.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_surf_t_p ( t, t_p, tt_current, coef_1, coef_2, c )

    REAL(wp), INTENT(IN) ::  c       !< total layer heat capacity
    REAL(wp), INTENT(IN) ::  coef_1  !< coefficient A in the prognostic equation
    REAL(wp), INTENT(IN) ::  coef_2  !< coefficient B in the prognostic equation
    REAL(wp), INTENT(IN) ::  t       !< current layer temperature

    REAL(wp), INTENT(OUT) ::  t_p  !< new layer temperature

    REAL(wp), INTENT(INOUT) ::  tt_current  !< current temperature tendency

    REAL(wp) ::  tt_new  !< new temperature tendency


!
!-- Compute the prognostic temperature without RK weighting.
    t_p = ( coef_1 * dt_3d * tsc(2) + c * t )  / ( c + coef_2 * dt_3d * tsc(2) )

!
!-- Compute the RK3 tendency for next time step.
    IF ( c /= 0.0_wp )  THEN

       t_p    = t_p + dt_3d * tsc(3) * tt_current
       tt_new = ( t_p - t - dt_3d * tsc(3) * tt_current ) / ( dt_3d  * tsc(2) )

       CALL calc_rk3_tend( tt_current, tt_new )

    ENDIF

 END SUBROUTINE calc_surf_t_p


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the new layer prognostic temperature by solving the Fourier diffusion equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_heat_diffusion ( t, t_p, tt_current, c, lambda, t_bc, sw_in, phi )

    REAL(wp), INTENT(IN) ::  t_bc  !< temperature boundary condition

    REAL(wp), INTENT(IN), OPTIONAL ::  sw_in  !< incoming shortwave radiation for windows

    REAL(wp), DIMENSION(:), INTENT(IN) ::  c       !< total heat capacity of the layer
    REAL(wp), DIMENSION(:), INTENT(IN) ::  lambda  !< total heat conductivity between layers
    REAL(wp), DIMENSION(:), INTENT(IN) ::  t       !< current time level temperature

    REAL(wp), DIMENSION(:), INTENT(IN), OPTIONAL ::  phi  !< fraction of incoming shortwave radiation absorbed at window layer

    REAL(wp), DIMENSION(:), INTENT(OUT) ::  t_p  !< new layer temperature

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  tt_current  !< current temperature tendency

    INTEGER(iwp) ::  k  !< material layer loop index

    REAL(wp) ::  tt_new  !<  new temperature tendency


!
!-- Loop through non-boundary layers of the material.
!-- @todo Split loop into three to move IFs out for better vecotrization.
    DO  k = LBOUND( t, 1 ) + 1, UBOUND( t, 1 )
!
!--    New prognostic layer temperature.
!--    Compute the diffusion between neighbouring layers.
       IF ( k /= UBOUND( t , 1 ) )  THEN
          tt_new = ( 1.0_wp / c(k) ) * ( lambda(k) * ( t(k+1) - t(k) ) +                           &
                   lambda(k-1) * ( t(k-1) - t(k) ) )
       ELSE
!
!--    Use a constant value boundary condition (skin temperature) for the innermost layer.
          tt_new = ( 1.0_wp / c(k) ) * ( lambda(k) * ( t_bc - t(k) ) +                             &
                   lambda(k-1) * ( t(k-1) - t(k) ) )
       ENDIF
!
!--    Add tendency from absorbed shortwave radiation.
       IF ( PRESENT( sw_in ) )  THEN
          tt_new = tt_new + ( 1.0_wp / c(k) ) * sw_in * phi(k)
       ENDIF
!
!--    Compute the prognostic temperature and RK3 tendency for next time step.
       t_p(k) = t(k) + dt_3d * ( tsc(2) * tt_new + tsc(3) * tt_current(k) )

       CALL calc_rk3_tend( tt_current(k), tt_new )
    ENDDO

 END SUBROUTINE calc_heat_diffusion


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the weighted RK3 tendency for the next timestep
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rk3_tend ( tend_current, tend_new )

    REAL(wp), INTENT(IN) ::  tend_new  !< new RK3 tendency

    REAL(wp), INTENT(INOUT) ::  tend_current  !< current time level tendency


    IF ( runge_l )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          tend_current = tend_new
       ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
          tend_current = -9.5625_wp * tend_new + 5.3125_wp * tend_current
       ENDIF
    ENDIF

 END SUBROUTINE calc_rk3_tend


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance and subsurface heat diffusion for roofs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE roof_model

    REAL(wp) ::  coef_1      !< coefficient A of the prognostic equation
    REAL(wp) ::  coef_2      !< coefficient B of the prognostic equation
    REAL(wp) ::  dq_s_dt     !< water vapour mixing ratio tendency
    REAL(wp) ::  e_s         !< saturation water vapour pressure
    REAL(wp) ::  e_s_dt      !< saturation water vapour pressure tendency
    REAL(wp) ::  f_shf       !< factor for the roof sensible heat flux
    REAL(wp) ::  f_qsws_liq  !< factor for the latent heat flux from/to liquid water reservoir
    REAL(wp) ::  tm_new      !< new liquid water reservoir tendency


!
!-- Surface sensible heat flux factor.
    f_shf = rho_cp / surf%rah_roof(m)

!
!-- Compute the nominator and denominator coefficients in
!-- the prognostic equation for the moist case.
    IF ( moist_physics )  THEN
!
!--    Computation of factor for the latent heat flux due to
!--    liquid water reservoir evaporation/condensation.
       e_s = surface_pressure * surf%qs_roof(m) / ( surf%qs_roof(m) + rd_d_rv )

!
!--    In case of evaporation, evaporate only for the liquid water coverage area,
!--    in case of condensation, use the total surface.
       IF ( surf%qs_roof(m) > surf%q1(m) )  THEN
          f_qsws_liq = rho_lv * surf%c_liq_roof(m) / surf%rah_roof(m)
       ELSE
          f_qsws_liq = rho_lv / surf%rah_roof(m)
       ENDIF

       e_s_dt = e_s * ( 17.62_wp / ( surf%t_roof(nzt_roof,m) -  29.65_wp ) -                       &
                        17.62_wp * ( surf%t_roof(nzt_roof,m) - 273.15_wp ) /                       &
                        ( surf%t_roof(nzt_roof,m) - 29.65_wp )**2                                  &
                      )

       dq_s_dt = rd_d_rv * e_s_dt / ( surface_pressure - e_s_dt )

!
!--    The coefficients for the moist prognostic equation for temperature.
       coef_1 = surf%rad_sw_net_roof(m) + surf%rad_lw_net_roof(m)                                  &
                - 3.0_wp * surf%lw_roof_coef(1,m) * surf%t_roof(nzt_roof,m)**4                     &
                + f_shf * surf%pt1(m)                                                              &
                + f_qsws_liq * ( surf%q1(m) - surf%qs_roof(m)                                      &
                                 + dq_s_dt * surf%t_roof(nzt_roof,m) )                             &
                + surf%conductivity_roof(nzt_roof,m) * surf%t_roof(nzt_roof+1,m)

       coef_2 = -4.0_wp * surf%lw_roof_coef(1,m) * surf%t_roof(nzt_roof,m)**3                      &
                + f_shf * d_exner(k_topo)                                                          &
                + f_qsws_liq * dq_s_dt                                                             &
                + surf%conductivity_roof(nzt_roof,m)

    ELSE
!
!-- The coefficients for the dry prognostic equation for temperature.
       coef_1 = surf%rad_sw_net_roof(m) + surf%rad_lw_net_roof(m)                                  &
                -3.0_wp * surf%lw_roof_coef(1,m) * surf%t_roof(nzt_roof,m)**4                      &
                + f_shf * surf%pt1(m)                                                              &
                + surf%conductivity_roof(nzt_roof,m) * surf%t_roof(nzt_roof+1,m)

       coef_2 = -4.0_wp * surf%lw_roof_coef(1,m) * surf%t_roof(nzt_roof,m)**3                      &
                + f_shf * d_exner(k_topo)                                                          &
                + surf%conductivity_roof(nzt_roof,m)
    ENDIF

    CALL calc_surf_t_p( surf%t_roof(nzt_roof,m), surf%t_roof_p(nzt_roof,m),                        &
                        surf%tt_roof(nzt_roof,m), coef_1, coef_2, surf%c_roof(nzt_roof,m) )

!
!-- Explicit solution of the Fourier heat equation for the subsurface layers.
    CALL calc_heat_diffusion( surf%t_roof(:,m), surf%t_roof_p(:,m), surf%tt_roof(:,m),             &
                              surf%c_roof(:,m), surf%conductivity_roof(:,m), surf%t_indoor(m) )

!
!-- Compute the diagnostic fluxes for the roof surface.
    surf%ghf_roof(m) = surf%conductivity_roof(nzb_roof,m) *                                        &
                       ( surf%t_roof_p(nzb_roof,m) - surf%t_indoor(m) )

    surf%pt_roof(m) = surf%t_roof_p(nzt_roof,m) * d_exner(k_topo)

    surf%shf_roof(m) = -f_shf * ( surf%pt1(m) - surf%pt_roof(m) ) / c_p

!
!-- Update longwave radiative flux following linearization.
    surf%rad_lw_net_roof(m) = surf%rad_lw_net_roof(m)                                              &
                              + surf%lw_roof_coef(1,m) * surf%t_roof(nzt_roof,m)**4                &
                              - 4.0_wp * surf%lw_roof_coef(1,m) * surf%t_roof(nzt_roof,m)**3       &
                              * ( surf%t_roof(nzt_roof,m) - surf%t_roof_p(nzt_roof,m) )

!
!-- Compute the water vapor flux from/to liquid water reservoir and the prognostic reservoir level.
    IF ( moist_physics )  THEN
       surf%qsws_liq_roof(m) = -f_qsws_liq * ( surf%q1(m) - surf%qs_roof(m) +                      &
                                               dq_s_dt * surf%t_roof(nzt_roof,m) -                 &
                                               dq_s_dt * surf%t_roof_p(nzt_roof,m)                 &
                                             )

       surf%qsws_roof(m) = surf%qsws_liq_roof(m)
!
!
!--    Modification due to precipitiation. If the liquid reservoir is full, the liquid water
!--    is assumed to be drained into the drainage system (liquid water is not conserved).
!--    The precipitation flux is not included in the surface-atmosphere latent heat flux (qsws).
       IF ( precipitation )  THEN
          IF ( surf%m_liq_roof(m) < m_liq_max_roof )  THEN
             surf%qsws_liq_roof(m) = surf%qsws_roof(m) -                                           &
                                     prr(k_atm,j,i) * hyrho(k_atm) * 0.001_wp * rho_l * l_v
          ENDIF
       ENDIF

!
!--    Compute the total latent heat flux.
       surf%qsws_roof(m) = surf%qsws_roof(m) / l_v
!
!--    Compute the prognostic liquid water reservoir.
       tm_new = - surf%qsws_liq_roof(m) * drho_l_lv
       surf%m_liq_roof_p(m) = surf%m_liq_roof(m) +                                                 &
                              dt_3d * ( tsc(2) * tm_new + tsc(3) * surf%tm_liq_roof(m) )
!
!--    Check if the liquid water reservoir is overfull. If so, drain excess to the
!--    assumed drainage system (water is not conserved here).
       surf%m_liq_roof_p(m) = MIN( surf%m_liq_roof_p(m), m_liq_max_roof )
!
!--    Check for negative water reservoir. @todo store the removed water as runoff for output.
       surf%m_liq_roof_p(m) = MAX( surf%m_liq_roof_p(m), 0.0_wp )
!
!--    Compute RK3 tendency.
       CALL calc_rk3_tend( surf%tm_liq_roof(m), tm_new )
!
!--    Compute the new liquid water coverage.
       surf%c_liq_roof(m) = MIN( 1.0_wp, ( surf%m_liq_roof_p(m) / m_liq_max_roof )**0.67 )
!
!--    Compute new saturation mixing ratio.
       e_s = 0.01_wp * magnus( MIN( surf%t_roof_p(nzt_roof,m), 333.15_wp ) )
       surf%qs_roof(m) = rd_d_rv * e_s / ( surface_pressure - e_s )
!
!--    Calculate new mixing ratio and vpt at roof surface.
       surf%q_roof(m) = q_surf( surf%qs_roof(m), surf%rah_roof(m), surf%q1(m), f_qsws_liq )
       surf%vpt_roof(m) = surf%pt_roof(m) * ( 1.0_wp + 0.61_wp * surf%q_roof(m) )

    ENDIF

 END SUBROUTINE roof_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance and subsurface heat diffusion for roads.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE road_model

    REAL(wp) ::  coef_1      !< coefficient A of the prognostic equation
    REAL(wp) ::  coef_2      !< coefficient B of the prognostic equation
    REAL(wp) ::  dq_s_dt     !< water vapour mixing ratio tendency
    REAL(wp) ::  e_s         !< saturation water vapour pressure
    REAL(wp) ::  e_s_dt      !< saturation water vapour pressure tendency
    REAL(wp) ::  f_shf       !< factor for the sensible heat flux from roads
    REAL(wp) ::  f_qsws_liq  !< factor for the latent heat flux from/to liquid water reservoir
    REAL(wp) ::  tm_new      !< new liquid water reservoir tendency


!
!-- Surface sensible heat flux factor.
    f_shf = rho_cp / surf%rah_road(m)
!
!-- Compute the nominator and denominator coefficients in
!-- the prognostic equation for the moist case.
    IF ( moist_physics )  THEN
!
!--    Computation of factor for the latent heat flux due to
!--    liquid water reservoir evaporation/condensation.
       e_s = surface_pressure * surf%qs_road(m) / ( surf%qs_road(m) + rd_d_rv )

!
!--    In case of evaporation, evaporate only for the liquid water coverage area,
!--    in case of condensation, use the total surface.
       IF ( surf%qs_road(m) > surf%q_can(m) )  THEN
          f_qsws_liq = rho_lv * surf%c_liq_road(m) / surf%rah_road(m)
       ELSE
          f_qsws_liq = rho_lv / surf%rah_road(m)
       ENDIF

       e_s_dt = e_s * ( 17.62_wp / ( surf%t_road(nzt_road,m) - 29.65_wp ) -                        &
                        17.62_wp * ( surf%t_road(nzt_road,m) - 273.15_wp ) /                       &
                        ( surf%t_road(nzt_road,m) - 29.65_wp)**2                                   &
                      )

       dq_s_dt = rd_d_rv * e_s_dt / ( surface_pressure - e_s_dt )

!
!--    The coefficients for the moist prognostic equation for temperature. For the longwave balance,
!--    both direct emission and the effect of backreflection are linearized.
       coef_1 = surf%rad_sw_net_road(m) + surf%rad_lw_net_road(m)                                  &
                -3.0_wp * surf%lw_road_coef(1,m) * surf%t_road(nzt_road,m)**4                      &
                + f_shf * surf%t_can(m)                                                            &
                + f_qsws_liq * ( surf%q_can(m) - surf%qs_road(m)                                   &
                                 + dq_s_dt * surf%t_road(nzt_road,m) )                             &
                + surf%conductivity_road(nzt_road,m) * surf%t_road(nzt_road+1,m)

       coef_2 = -4.0_wp * surf%lw_road_coef(1,m) * surf%t_road(nzt_road,m)**3                      &
                + f_shf                                                                            &
                + f_qsws_liq * dq_s_dt                                                             &
                + surf%conductivity_road(nzt_road,m)

    ELSE
!
!--    The coefficients for the dry prognostic equation for temperature.
       coef_1 = surf%rad_sw_net_road(m) + surf%rad_lw_net_road(m)                                  &
                -3.0_wp * surf%lw_road_coef(1,m) * surf%t_road(nzt_road,m)**4                      &
                + f_shf * surf%t_can(m)                                                            &
                + surf%conductivity_road(nzt_road,m) * surf%t_road(nzt_road+1,m)

       coef_2 = -4.0_wp * surf%lw_road_coef(1,m) * surf%t_road(nzt_road,m)**3                      &
                + f_shf                                                                            &
                + surf%conductivity_road(nzt_road,m)
    ENDIF

    CALL calc_surf_t_p( surf%t_road(nzt_road,m), surf%t_road_p(nzt_road,m),                        &
                        surf%tt_road(nzt_road,m), coef_1, coef_2, surf%c_road(nzt_road,m) )

!
!-- Heat diffusion through subsurface layers.
    CALL calc_heat_diffusion( surf%t_road(:,m), surf%t_road_p(:,m), surf%tt_road(:,m),             &
                              surf%c_road(:,m), surf%conductivity_road(:,m), surf%t_soil(m) )

    surf%shf_road(m) = -f_shf * ( surf%t_can(m) - surf%t_road_p(nzt_road,m) ) / c_p

    surf%pt_road(m)  = surf%t_road_p(nzt_road,m) * d_exner(k_topo)

    surf%ghf_road(m) = surf%conductivity_road(nzb_road,m) *                                        &
                       ( surf%t_road_p(nzb_road,m) - surf%t_soil(m) )

!
!-- Update longwave radiative flux following linearization.
    surf%rad_lw_net_road(m) = surf%rad_lw_net_road(m)                                              &
                              + surf%lw_road_coef(1,m) * surf%t_road(nzt_road,m)**4                &
                              - 4.0_wp * surf%lw_road_coef(1,m) * surf%t_road(nzt_road,m)**3       &
                              * ( surf%t_road(nzt_road,m) - surf%t_road_p(nzt_road,m) )

!
!-- Compute the water vapor flux from/to liquid water reservoir and the prognostic reservoir level.
    IF ( moist_physics )  THEN
       surf%qsws_liq_road(m) = -f_qsws_liq * ( surf%q_can(m) - surf%qs_road(m) +                   &
                                               dq_s_dt * surf%t_road(nzt_road,m) -                 &
                                               dq_s_dt * surf%t_road_p(nzt_road,m)                 &
                                             )

       surf%qsws_road(m) = surf%qsws_liq_road(m)
!
!--    Modification due to precipitiation. If the liquid reservoir is full, the liquid water
!--    is assumed to be drained into the drainage system (liquid water is not conserved).
!--    The precipitation flux is not included in the surface-atmosphere latent heat flux (qsws).
       IF ( precipitation )  THEN
          IF ( surf%m_liq_road(m) < m_liq_max_road )  THEN
             surf%qsws_road(m) = surf%qsws_road(m) -                                               &
                                 prr(k_atm,j,i+i) * hyrho(k_atm) * 0.001_wp * rho_l * l_v
          ENDIF
       ENDIF
!
!--    Compute the total latent heat flux.
       surf%qsws_road(m) = surf%qsws_road(m) / l_v
!
!--    Compute the prognostic liquid water reservoir.
       tm_new = - surf%qsws_liq_road(m) * drho_l_lv
       surf%m_liq_road_p(m) = surf%m_liq_road(m) +                                                 &
                              dt_3d * ( tsc(2) * tm_new + tsc(3) * surf%tm_liq_road(m) )
!
!--    Check if the liquid water reservoir is overfull. If so, drain excess to the
!--    assumed drainage system (water is not conserved here).
       surf%m_liq_road_p(m) = MIN( surf%m_liq_road_p(m), m_liq_max_road )
!
!--    Check for negative water reservoir. Should we adjust qsws_road accordingly?
       surf%m_liq_road_p(m) = MAX( surf%m_liq_road_p(m), 0.0_wp )
!
!--    Compute RK3 tendency
       CALL calc_rk3_tend( surf%tm_liq_road(m), tm_new )
!
!--    Compute the new liquid water coverage.
       surf%c_liq_road(m) = MIN( 1.0_wp, ( surf%m_liq_road_p(m) / m_liq_max_road )**0.67 )
!
!--    Compute new saturation mixing ratio.
       e_s = 0.01_wp * magnus( MIN( surf%t_road_p(nzt_road,m), 333.15_wp ) )
       surf%qs_road(m) = rd_d_rv * e_s / ( surface_pressure - e_s )
!
!--    Calculate new mixing ratio and vpt at road surface.
       surf%q_road(m) = q_surf( surf%qs_road(m), surf%rah_road(m), surf%q1(m), f_qsws_liq )
       surf%vpt_road(m) = surf%pt_road(m) * ( 1.0_wp + 0.61_wp * surf%q_road(m) )

    ENDIF

 END SUBROUTINE road_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance and subsurface heat diffusion for both walls.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wall_model

   REAL(wp) ::  coef_1   !< coefficient A of the prognostic equation
   REAL(wp) ::  coef_2   !< coefficient B of the prognostic equation
   REAL(wp) ::  f_shf_a  !< factor for the wall surface heat flux
   REAL(wp) ::  f_shf_b  !< factor for the wall surface heat flux


    IF ( facade_rah_doe )  THEN
       f_shf_a = rho_cp / surf%rah_wall_a(m)
       IF ( surf%anisotropic_canyon(m) )  f_shf_b = rho_cp / surf%rah_wall_b(m)
    ELSE
       f_shf_a = rho_cp / surf%rah_facade(m)
       IF ( surf%anisotropic_canyon(m) )  f_shf_b = f_shf_a
    ENDIF

!
!-- The coefficients for the moist prognostic equation for temperature. For the longwave balance,
!-- both direct emission and the effect of backreflection are linearized. The linearization depends
!-- if the canyon is isotropic or not as an average backreflection is used for isotropic canyons.
!-- We consider the walls are dry in all cases, so moist physical processes are not considered.
    IF ( surf%anisotropic_canyon(m) )  THEN
       coef_1 = surf%rad_sw_net_wall_a(m) + surf%rad_lw_net_wall_a(m)                              &
                - 3.0_wp * surf%lw_wall_coef(1,m) * surf%t_wall_a(nzt_wall,m)**4                   &
                + f_shf_a * surf%t_can(m)                                                          &
                + surf%conductivity_wall(nzt_wall,m) * surf%t_wall_a(nzt_wall+1,m)

       coef_2 = -4.0_wp * surf%lw_wall_coef(1,m) * surf%t_wall_a(nzt_wall,m)**3                    &
                + f_shf_a                                                                          &
                + surf%conductivity_wall(nzt_wall,m)

       CALL calc_surf_t_p(surf%t_wall_a(nzt_wall,m), surf%t_wall_a_p(nzt_wall,m),                  &
                          surf%tt_wall_a(nzt_wall,m), coef_1, coef_2, surf%c_wall(nzt_wall,m) )

       coef_1 = surf%rad_sw_net_wall_b(m) + surf%rad_lw_net_wall_b(m)                              &
                - 3.0_wp * surf%lw_wall_coef(1,m) * surf%t_wall_b(nzt_wall,m)**4                   &
                + f_shf_b * surf%t_can(m)                                                          &
                + surf%conductivity_wall(nzt_wall,m) * surf%t_wall_b(nzt_wall+1,m)

       coef_2 = -4.0_wp * surf%lw_wall_coef(1,m) * surf%t_wall_b(nzt_wall,m)**3                    &
                + f_shf_b                                                                          &
                + surf%conductivity_wall(nzt_wall,m)

       CALL calc_surf_t_p( surf%t_wall_b(nzt_wall,m), surf%t_wall_b_p(nzt_wall,m),                 &
                           surf%tt_wall_b(nzt_wall,m), coef_1, coef_2, surf%c_wall(nzt_wall,m) )
    ELSE
!
!--    In case of isotropic canyon, wall A and B temperatures are averaged, and thus the prognostic
!--    equation for t_wall_a is representative of both of the walls. Thus both the terms for
!--    t_wall_a as well as for t_wall_b in the longwave radiation balance has a dependency on
!--    the surface temperature.
       coef_1 = surf%rad_sw_net_wall_a(m) + surf%rad_lw_net_wall_a(m)                              &
                - 3.0_wp * ( surf%lw_wall_coef(1,m) + surf%lw_wall_coef(3,m) )                     &
                   * surf%t_wall_a(nzt_wall,m)**4                                                  &
                + f_shf_a * surf%t_can(m)                                                          &
                + surf%conductivity_wall(nzt_wall,m) * surf%t_wall_a(nzt_wall+1,m)

       coef_2 = -4.0_wp * ( surf%lw_wall_coef(1,m) + surf%lw_wall_coef(3,m) )                      &
                   * surf%t_wall_a(nzt_wall,m)**3                                                  &
                + f_shf_a                                                                          &
                + surf%conductivity_wall(nzt_wall,m)

       CALL calc_surf_t_p( surf%t_wall_a(nzt_wall,m), surf%t_wall_a_p(nzt_wall,m),                 &
                           surf%tt_wall_a(nzt_wall,m), coef_1, coef_2, surf%c_wall(nzt_wall,m) )
    ENDIF

    surf%pt_wall_a(m)  = surf%t_wall_a_p(nzt_wall,m) * d_exner(k_topo)
    surf%shf_wall_a(m) = -f_shf_a * ( surf%t_can(m) - surf%t_wall_a_p(nzt_wall,m) ) / c_p

!
!-- Heat diffusion through subsurface layers.
    CALL calc_heat_diffusion( surf%t_wall_a(:,m), surf%t_wall_a_p(:,m), surf%tt_wall_a(:,m),       &
                              surf%c_wall(:,m), surf%conductivity_wall(:,m), surf%t_indoor(m) )

    surf%ghf_wall_a(m) = surf%conductivity_wall(nzb_wall,m) *                                      &
                         ( surf%t_wall_a_p(nzb_wall,m) - surf%t_indoor(m) )

!
!-- Same treatment for wall B if this is an anisotropic canyon, otherwise copy.
    IF ( surf%anisotropic_canyon(m) )  THEN
       surf%pt_wall_b(m)  = surf%t_wall_b_p(nzt_wall,m) * d_exner(k_topo)
       surf%shf_wall_b(m) = -f_shf_b * ( surf%t_can(m) - surf%t_wall_b_p(nzt_wall,m) ) / c_p

       CALL calc_heat_diffusion( surf%t_wall_b(:,m), surf%t_wall_b_p(:,m), surf%tt_wall_b(:,m),    &
                                 surf%c_wall(:,m), surf%conductivity_wall(:,m), surf%t_indoor(m) )

       surf%ghf_wall_b(m) = surf%conductivity_wall(nzb_wall,m) *                                   &
                            ( surf%t_wall_b_p(nzb_wall,m) - surf%t_indoor(m) )

!
!--    Update longwave radiative fluxes following linearization.
       surf%rad_lw_net_wall_a(m) = surf%rad_lw_net_wall_a(m)                                       &
                                   + surf%lw_wall_coef(1,m) * surf%t_wall_a(nzt_wall,m)**4         &
                                   - 4.0_wp * surf%lw_wall_coef(1,m)                               &
                                      * surf%t_wall_a(nzt_wall,m)**3                               &
                                   * ( surf%t_wall_a(nzt_wall,m) - surf%t_wall_a_p(nzt_wall,m) )

       surf%rad_lw_net_wall_b(m) = surf%rad_lw_net_wall_b(m)                                       &
                                   + surf%lw_wall_coef(1,m) * surf%t_wall_b(nzt_wall,m)**4         &
                                   - 4.0_wp * surf%lw_wall_coef(1,m)                               &
                                      * surf%t_wall_b(nzt_wall,m)**3                               &
                                   * ( surf%t_wall_b(nzt_wall,m) - surf%t_wall_b_p(nzt_wall,m) )
    ELSE
!
!--    Copy all layers including the surface for wall B.
       surf%t_wall_b_p(:,m) = surf%t_wall_a_p(:,m)
       surf%tt_wall_b(:,m)  = surf%tt_wall_a(:,m)
       surf%pt_wall_b(m)    = surf%pt_wall_a(m)
       surf%shf_wall_b(m)   = surf%shf_wall_a(m)
       surf%ghf_wall_b(m)   = surf%ghf_wall_a(m)

!
!--    For longwave radiative flux, we need to add terms for both the t_wall_a and the t_wall_b
!--    in the longwave balance, thus coefficients 1 and 3 are summed here.
       surf%rad_lw_net_wall_a(m) = surf%rad_lw_net_wall_a(m)                                       &
                                   + ( surf%lw_wall_coef(1,m) + surf%lw_wall_coef(3,m) )           &
                                      * surf%t_wall_a(nzt_wall,m)**4                               &
                                   - 4.0_wp * ( surf%lw_wall_coef(1,m)                             &
                                                + surf%lw_wall_coef(3,m) )                         &
                                   * surf%t_wall_a(nzt_wall,m)**3                                  &
                                   * ( surf%t_wall_a(nzt_wall,m) - surf%t_wall_a_p(nzt_wall,m) )
       surf%rad_lw_net_wall_b(m) = surf%rad_lw_net_wall_a(m)
    ENDIF

 END SUBROUTINE wall_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance, SW transmission and subsurface heat diffusion for windows.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE window_model

    REAL(wp) ::  coef_1   !< coefficient A of the prognostic equation
    REAL(wp) ::  coef_2   !< coefficient B of the prognostic equation
    REAL(wp) ::  f_shf_a  !< factor for the window surface heat flux
    REAL(wp) ::  f_shf_b  !< factor for the window surface heat flux


    IF ( facade_rah_doe )  THEN
       f_shf_a = rho_cp / surf%rah_win_a(m)
       IF ( surf%anisotropic_canyon(m) )  f_shf_b = rho_cp / surf%rah_win_b(m)
    ELSE
       f_shf_a = rho_cp / surf%rah_facade(m)
       IF ( surf%anisotropic_canyon(m) )  f_shf_b = f_shf_a
    ENDIF

!
!-- Computation of the prognostic equation similarly to the walls, with exception of added
!-- shortwave transmission component for surface and subsurface layers. Explanatory comments
!-- are not repeated from the wall model, comments reflect differences specific to windows.
    IF ( surf%anisotropic_canyon(m) )  THEN
!
!--    For windows, some of the incoming shortwave radiation is transmitted through the material.
       coef_1 = surf%rad_sw_net_win_a(m) * surf%absorption_win(nzt_win,m)                          &
                + surf%rad_lw_net_win_a(m)                                                         &
                - 3.0_wp * surf%lw_win_coef(1,m) * surf%t_win_a(nzt_win,m)**4                      &
                + f_shf_a * surf%t_can(m)                                                          &
                + surf%conductivity_win(nzt_win,m) * surf%t_win_a(nzt_win+1,m)

       coef_2 = -4.0_wp * surf%lw_win_coef(1,m) * surf%t_win_a(nzt_win,m)**3                       &
                + f_shf_a                                                                          &
                + surf%conductivity_win(nzt_win,m)

       CALL calc_surf_t_p( surf%t_win_a(nzt_win,m), surf%t_win_a_p(nzt_win,m),                     &
                           surf%tt_win_a(nzt_win,m), coef_1, coef_2, surf%c_win(nzt_win,m) )

       coef_1 = surf%rad_sw_net_win_b(m) * surf%absorption_win(nzt_win,m)                          &
                + surf%rad_lw_net_win_b(m)                                                         &
                - 3.0_wp * surf%lw_win_coef(1,m) * surf%t_win_b(nzt_win,m)**4                      &
                + f_shf_b * surf%t_can(m)                                                          &
                + surf%conductivity_win(nzt_win,m) * surf%t_win_b(nzt_win+1,m)

       coef_2 = -4.0_wp * surf%lw_win_coef(1,m) * surf%t_win_b(nzt_win,m)**3                       &
                + f_shf_b                                                                          &
                + surf%conductivity_win(nzt_win,m)

       CALL calc_surf_t_p( surf%t_win_b(nzt_win,m), surf%t_win_b_p(nzt_win,m),                     &
                           surf%tt_win_b(nzt_win,m), coef_1, coef_2, surf%c_win(nzt_win,m) )

    ELSE
       coef_1 = surf%rad_sw_net_win_a(m) * surf%absorption_win(nzt_win,m)                          &
                + surf%rad_lw_net_win_a(m)                                                         &
                - 3.0_wp * ( surf%lw_win_coef(1,m) + surf%lw_win_coef(3,m) )                       &
                   * surf%t_win_a(nzt_win,m)**4                                                    &
                + f_shf_a * surf%t_can(m)                                                          &
                + surf%conductivity_win(nzt_win,m) * surf%t_win_a(nzt_win+1,m)

       coef_2 = -4.0_wp * ( surf%lw_win_coef(1,m) + surf%lw_win_coef(3,m) )                        &
                   * surf%t_win_a(nzt_win,m)**3                                                    &
                + f_shf_a                                                                          &
                + surf%conductivity_win(nzt_win,m)

       CALL calc_surf_t_p( surf%t_win_a(nzt_win,m), surf%t_win_a_p(nzt_win,m),                     &
                           surf%tt_win_a(nzt_win,m), coef_1, coef_2, surf%c_win(nzt_win,m) )
    ENDIF

    surf%pt_win_a(m)  = surf%t_win_a_p(nzt_win,m) * d_exner(k_topo)
    surf%shf_win_a(m) = -f_shf_a * ( surf%t_can(m) - surf%t_win_a_p(nzt_win,m) ) / c_p

!
!-- The transmitted shortwave radiation is included also in the prognostic equations for material
!-- subsurface temperatures.
    CALL calc_heat_diffusion( surf%t_win_a(:,m), surf%t_win_a_p(:,m),                              &
                              surf%tt_win_a(:,m), surf%c_win(:,m),                                 &
                              surf%conductivity_win(:,m), surf%t_indoor(m),                        &
                              surf%rad_sw_net_win_a(m), surf%absorption_win(:,m) )

    surf%ghf_win_a(m) = surf%conductivity_win(nzb_win,m) *                                         &
                        ( surf%t_win_a_p(nzb_win,m) - surf%t_indoor(m) )

    IF ( surf%anisotropic_canyon(m) )  THEN
       surf%pt_win_b(m)  = surf%t_win_b_p(nzt_win,m) * d_exner(k_topo)
       surf%shf_win_b(m) = -f_shf_b * ( surf%t_can(m) - surf%t_win_b_p(nzt_win,m) ) / c_p

       CALL calc_heat_diffusion( surf%t_win_b(:,m), surf%t_win_b_p(:,m),                           &
                                 surf%tt_win_b(:,m), surf%c_win(:,m),                              &
                                 surf%conductivity_win(:,m), surf%t_indoor(m),                     &
                                 surf%rad_sw_net_win_b(m), surf%absorption_win(:,m) )

       surf%ghf_win_b(m) = surf%conductivity_win(nzb_win,m) *                                      &
                           ( surf%t_win_b_p(nzb_win,m) - surf%t_indoor(m) )

       surf%rad_lw_net_win_a(m) = surf%rad_lw_net_win_a(m)                                         &
                                  + surf%lw_win_coef(1,m) * surf%t_win_a(nzt_win,m)**4             &
                                  - 4.0_wp * surf%lw_win_coef(1,m) * surf%t_win_a(nzt_win,m)**3    &
                                     * ( surf%t_win_a(nzt_win,m) - surf%t_win_a_p(nzt_win,m) )

       surf%rad_lw_net_win_b(m) = surf%rad_lw_net_win_b(m)                                         &
                                  + surf%lw_win_coef(1,m) * surf%t_win_b(nzt_win,m)**4             &
                                  - 4.0_wp * surf%lw_win_coef(1,m) * surf%t_win_b(nzt_win,m)**3    &
                                     * ( surf%t_win_b(nzt_win,m) - surf%t_win_b_p(nzt_win,m) )
    ELSE
       surf%t_win_b_p(:,m) = surf%t_win_a_p(:,m)
       surf%tt_win_b(:,m)  = surf%tt_win_a(:,m)
       surf%pt_win_b(m)    = surf%pt_win_a(m)
       surf%shf_win_b(m)   = surf%shf_win_a(m)
       surf%ghf_win_b(m)   = surf%ghf_win_a(m)

       surf%rad_lw_net_win_a(m) = surf%rad_lw_net_win_a(m)                                         &
                                  + ( surf%lw_win_coef(1,m) + surf%lw_win_coef(3,m) )              &
                                     * surf%t_win_a(nzt_win,m)**4                                  &
                                  - 4.0_wp * ( surf%lw_win_coef(1,m) + surf%lw_win_coef(3,m) )     &
                                  * surf%t_win_a(nzt_win,m)**3                                     &
                                  * ( surf%t_win_a(nzt_win,m) - surf%t_win_a_p(nzt_win,m) )
       surf%rad_lw_net_win_b(m) = surf%rad_lw_net_win_a(m)
    ENDIF

 END SUBROUTINE window_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface mixing ratio using resistance weighting.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION q_surf( q_s, rah, q_a, f_qsws )

    REAL(wp), INTENT(IN) ::  f_qsws  !< factor for the latent heat flux
    REAL(wp), INTENT(IN) ::  q_a     !< mixing ratio of adjacent air
    REAL(wp), INTENT(IN) ::  q_s     !< saturation mixing ratio at the surface
    REAL(wp), INTENT(IN) ::  rah     !< aerodynamic resistance for heat (and for water vapor)

    REAL(wp) ::  q_surf  !< mixing ratio for the surface
    REAL(wp) ::  res     !< total surface resistance


!
!-- Total surface resistance.
    res = rah / ( rah + ABS( rho_lv / ( f_qsws + 1.0E-20_wp ) - rah ) )

!
!-- Use the newly calculated total surface resistance to compute weighted surface mixing ratio.
    IF ( bulk_cloud_model )  THEN
!
!--    Assume equal liquid water content in canyon as in air above.
       q_surf = res * q_s + ( 1.0_wp - res ) * ( q_a - ql(k_atm,j,i) )
    ELSE
       q_surf = res * q_s + ( 1.0_wp - res ) * q_a
    ENDIF

 END FUNCTION q_surf

 END SUBROUTINE slurb_energy_balance_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes information on SLUrb setup to the PALM header file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_header ( io )

    INTEGER(iwp), INTENT(IN) ::  io  !< unit id of the output file


    WRITE( io,  1 )
    WRITE( io,  2 )
    IF ( moist_physics )  WRITE( io, 3 )
    IF (  .NOT.  moist_physics )  WRITE( io, 4 )
    IF ( anisotropic_street_canyons )  WRITE( io, 5 )
    IF (  .NOT.  anisotropic_street_canyons )  WRITE( io, 6 )
    WRITE( io, 7 )  TRIM( aero_roughness_heat )
    WRITE( io, 8 )  TRIM( facade_resistance_parametrization )
    WRITE( io, 9 )  TRIM( street_canyon_wspeed_factor )

1   FORMAT (//' Single-layer urban surface model (PALM-SLUrb):'/                                   &
             ' ------------------------------'/)
2   FORMAT ('  --> Module enabled.')
3   FORMAT ('   Moist physics enabled.')
4   FORMAT ('   Moist physics disabled.')
5   FORMAT ('   Anisotropic street canyons allowed.')
6   FORMAT ('   Anisotropic street canyons not allowed.')
7   FORMAT ('   Parametrization for aerodynamic roughness for heat on horizontal surfaces: ',A)
8   FORMAT ('   Parametrization for aerodynamic resistance for heat vertical surfaces: ',A)
9   FORMAT ('   Parametrization for street canyon wind speed: ',A)

 END SUBROUTINE slurb_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes the SLUrb model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_init

    USE control_parameters,                                                                        &
       ONLY:  coupling_char

    USE netcdf_data_input_mod,                                                                     &
       ONLY:  check_existence,                                                                     &
              close_input_file,                                                                    &
              get_dimension_length,                                                                &
              get_attribute,                                                                       &
              get_variable,                                                                        &
              inquire_num_variables,                                                               &
              inquire_variable_names,                                                              &
              open_read_file

    CHARACTER(LEN=100) ::  input_file_slurb = 'PIDS_SLURB'  !< name of driver file which comprises SLUrb input data

    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< array of variable names in the input driver

    INTEGER(iwp) ::  dimlen    !< lenght of dimension
    INTEGER(iwp) ::  i         !< loop index along x
    INTEGER(iwp) ::  id_slurb  !< netCDF id of the input netCDF file
    INTEGER(iwp) ::  j         !< loop index along y
    INTEGER(iwp) ::  k         !< layer index for loops
    INTEGER(iwp) ::  k_atm     !< k index of the first atmospheric grid level
    INTEGER(iwp) ::  k_topo    !< k index of the topography top
    INTEGER(iwp) ::  m         !< SLUrb tile index for loops
    INTEGER(iwp) ::  num_vars  !< number of variables in the input netCDF file

    LOGICAL ::  input_file_present  !< flag to indicate that the driver file has been found

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  canyon_orientation_tmp  !< temporary array which doesn't have a corresponding model array, they are used during
                                                                    !< the initialization to fill model arrays for latent variables of variable with a different data type.


    IF ( debug_output )  CALL debug_message( 'slurb_init', 'start' )
!
!-- Begin reading input data.
    INQUIRE( FILE = TRIM( input_file_slurb ) //  TRIM( coupling_char ), EXIST = input_file_present )

    IF ( input_file_present )  THEN
       CALL open_read_file( TRIM( input_file_slurb ) // TRIM( coupling_char ), id_slurb )
       CALL inquire_num_variables( id_slurb, num_vars )

!
!--    Allocate memory to store variable names.
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_slurb, var_names )

!
!--    Verify correct dimension lengths in the user input. The layer ids are up to the user.
       IF ( check_existence( var_names, 'nroof_3d' ) ) THEN
          CALL get_dimension_length( id_slurb, dimlen, 'nroad_3d' )
          IF ( dimlen /= n_layers_roads ) THEN
             WRITE( message_string, * ) 'Mismatch in the number of road layers in the SLUrb input' &
                                     // ' driver and the SLUrb model configuration (', nroad_3d,   &
                                        ' and', n_layers_roads, ' respectively).'
             CALL message( 'slurb_init', 'SLU0014', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( check_existence( var_names, 'nroof_3d' ) ) THEN
          CALL get_dimension_length( id_slurb, dimlen, 'nroof_3d' )
          IF ( dimlen /= n_layers_roofs ) THEN
             WRITE( message_string, * ) 'Mismatch in the number of roof layers in the SLUrb input' &
                                     // ' driver and the SLUrb model configuration (', nroof_3d,   &
                                        ' and', n_layers_roofs, ' respectively).'
             CALL message( 'slurb_init', 'SLU0015', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( check_existence( var_names, 'nwall_3d' ) ) THEN
          CALL get_dimension_length( id_slurb, dimlen, 'nwall_3d' )
          IF ( dimlen /= n_layers_walls ) THEN
             WRITE( message_string, * ) 'Mismatch in the number of wall layers in the SLUrb input' &
                                     // ' driver and the SLUrb model configuration (', nwall_3d,   &
                                        ' and', n_layers_walls, ' respectively).'
             CALL message( 'slurb_init', 'SLU0016', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( check_existence( var_names, 'nwin_3d' ) ) THEN
          CALL get_dimension_length( id_slurb, dimlen, 'nwin_3d' )
          IF ( dimlen /= n_layers_windows ) THEN
             WRITE( message_string, * ) 'Mismatch in the number of win layers in the SLUrb input'  &
                                     // ' driver and the SLUrb model configuration (', nwin_3d,    &
                                        ' and', n_layers_windows, ' respectively).'
             CALL message( 'slurb_init', 'SLU0017', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Process variables related to urban form.
!
!-- Internally, f_bld refers to building plan area fraction of the urban surface. However, it is
!-- more common to report the building plan area fraction as a fraction of the total surface, e.g.
!-- in the case of LCZs. We want the user input correspond to the latter, thus the scaling.
    surf%f_bld(:) = -9999.0_wp
    CALL get_grid_variable_1d_real( 'building_plan_area_fraction', surf%f_bld,                     &
                                    building_plan_area_fraction )
    CALL check_grid_variable_1d_real( 'building_plan_area_fraction', surf%f_bld,                   &
                                      TINY( 1.0_wp ), 0.99_wp)
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       IF ( surf%f_bld(m) > fr_urb(j,i) ) THEN
          WRITE( message_string, * ) 'building_plan_area_fraction = ', surf%f_bld(m),              &
                                     ' is higher than urban_fraction = ', fr_urb(j,i),             &
                                     ' for grid cell (j,i) = ', surf%j(m), surf%i(m), '.'
          CALL message( 'slurb_init', 'SLU1020', 2, 2, 0, 6, 0 )
       ENDIF
       surf%f_bld(m) = surf%f_bld(m) / fr_urb(j,i)
    ENDDO

    surf%f_bld_frn(:) = -9999.0_wp
    CALL get_grid_variable_1d_real( 'building_frontal_area_fraction', surf%f_bld_frn,              &
                                    building_frontal_area_fraction )
    CALL check_grid_variable_1d_real( 'building_frontal_area_fraction', surf%f_bld_frn,            &
                                      0.0_wp, HUGE( 1.0_wp ) )

    surf%h_bld(:) = -9999.0_wp
    CALL get_grid_variable_1d_real( 'building_height', surf%h_bld, building_height )
    CALL check_grid_variable_1d_real( 'building_height', surf%h_bld, 0.0_wp, 1000.0_wp )

!
!-- Urban surface and street canyon MOST heights.
    DO  m = 1, surf%ns
       surf%z_mo(m) = 0.5_wp *  dzw(topo_top_ind(j,i,0)+1)
       surf%z_mo_can(m) = 0.5_wp * surf%h_bld(m)
    ENDDO

!
!-- Process canyon direction information if anisotropic street canyons are enabled.
    IF ( anisotropic_street_canyons )  THEN
       ALLOCATE( canyon_orientation_tmp(1:surf%ns) )
       CALL get_grid_variable_1d_real( 'street_canyon_orientation', canyon_orientation_tmp,        &
                                       street_canyon_orientation )
!
!--    In order to make it possible to have a mix of isotropic and anisotropic tiles, use
!--    anisotropic canyons only in the case a canyon orientation has been given either for all tiles
!--    in the namelist or per-patch basis in input file.
       DO  m = 1, surf%ns
          IF ( canyon_orientation_tmp(m) /= -9999.0_wp )  THEN
             surf%anisotropic_canyon(m) = .TRUE.
             surf%theta_can(m) = canyon_orientation_tmp(m)  * ( pi / 180.0_wp )
          ELSE
             surf%anisotropic_canyon(m) = .FALSE.
             surf%theta_can(m) = -9999.0_wp
          ENDIF
       ENDDO
       DEALLOCATE( canyon_orientation_tmp )
    ELSE
       DO  m = 1, surf%ns
          surf%anisotropic_canyon(m) = .FALSE.
          surf%theta_can(m) = -9999.0_wp
       ENDDO
    ENDIF

    CALL get_grid_variable_1d_real( 'street_canyon_aspect_ratio', surf%hw_can,                     &
                                    street_canyon_aspect_ratio )
    CALL check_grid_variable_1d_real( 'street_canyon_aspect_ratio', surf%hw_can,                   &
                                      TINY( 1.0_wp ), HUGE( 0.0_wp ) )

    CALL get_grid_variable_1d_real( 'z0_urb', surf%z0_urb, urban_roughness_length )
    CALL check_grid_variable_1d_real( 'z0_urb', surf%z0_urb, TINY( 1.0_wp ), MINVAL( surf%z_mo ) )

!
!-- Initialize material parameters such as heat capacities.
    CALL process_surface_parameters

!
!-- End reading input data and continue initialization by precomputing variables.
!
!-- Compute latent variables which can be inferred from the inputs and do not update runtime.
    CALL precompute_latent_variables

!
!-- Initialization of thermodynamical variables
    CALL get_grid_variable_1d_real( 'building_indoor_temperature', surf%t_indoor,                  &
                                    building_indoor_temperature )
    CALL check_grid_variable_1d_real( 'building_indoor_temperature', surf%t_indoor,                &
                                      TINY( 1.0_wp ), HUGE( 1.0_wp ) )

    CALL get_grid_variable_1d_real( 'deep_soil_temperature', surf%t_soil, deep_soil_temperature )
    CALL check_grid_variable_1d_real( 'deep_soil_temperature', surf%t_soil,                        &
                                      TINY( 1.0_wp ), HUGE( 1.0_wp ) )

    CALL process_dynamic_inputs

    CALL init_slurb_variables

    IF ( debug_output )  CALL debug_message( 'slurb_init', 'end' )

    CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE check_grid_variable_1d_int( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    INTEGER(iwp), INTENT(IN) ::  valid_max  !< upper bound of valid range for var
    INTEGER(iwp), INTENT(IN) ::  valid_min  !< lower bound of valid range for var

    INTEGER(iwp), DIMENSION(:), INTENT(IN) ::  var  !< target variable


!
!-- Check if the input is within allowed bounds.
    DO  m = 1, surf%ns
       IF ( ( var(m) < valid_min  .OR.  var(m) > valid_max ) )  THEN
          WRITE( message_string, * ) 'Input variable ' // TRIM( varname ) //                       &
                                     ' for grid cell (j,i) = ', surf%j(m), surf%i(m), ' set to',   &
                                     var(m), ' valid range is [', valid_min, valid_max, '].'
          CALL message( 'slurb_init', 'SLU0018', 2, 2, myid, 6, 0 )
       ENDIF
    ENDDO

 END SUBROUTINE check_grid_variable_1d_int


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE check_grid_variable_1d_real( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(wp), INTENT(IN) ::  valid_max  !< upper bound of valid range for var
    REAL(wp), INTENT(IN) ::  valid_min  !< lower bound of valid range for var

    REAL(wp), DIMENSION(:), INTENT(IN) ::  var  !< target variable


!
!-- Check if the input is within allowed bounds.
    DO  m = 1, surf%ns
       IF ( ( var(m) < valid_min  .OR.  var(m) > valid_max ) )  THEN
          WRITE( message_string, * ) 'input variable ' // TRIM( varname ) //                       &
                                     ' for grid cell (j,i) = ', surf%j(m), surf%i(m), ' set to',   &
                                     var(m), ' valid range is [', valid_min, valid_max, '].'
          CALL message( 'slurb_init', 'SLU0019', 2, 2, myid, 6, 0 )
       ENDIF
    ENDDO

 END SUBROUTINE check_grid_variable_1d_real


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE check_grid_variable_2d_real( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(wp), INTENT(IN) ::  valid_max   !< upper bound of valid range for var
    REAL(wp), INTENT(IN) ::  valid_min   !< lower bound of valid range for var

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::  var   !< target variable


!
!-- Check if the input is within allowed bounds
    DO  m = 1, surf%ns
       DO  k = LBOUND( var, 1 ), UBOUND( var, 1 )
          IF ( ( var(k,m) < valid_min  .OR.  var(k,m) > valid_max ) )  THEN
             WRITE( message_string, * ) 'Input variable ' // TRIM( varname ) //                    &
                                        ' for grid cell (k,j,i) = ', k, surf%j(m), surf%i(m),      &
                                        ' set to', var(k,m), ' valid range is [',                  &
                                        valid_min, valid_max, '].'
             CALL message( 'slurb_init', 'SLU0020', 2, 2, myid, 6, 0 )
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE check_grid_variable_2d_real


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes SLUrb model variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_slurb_variables

    REAL(wp) ::  bc_atm  !< initial atmospheric boundary condition for temperature
    REAL(wp) ::  e_s     !< initial water vapor saturation pressure


    DO  m = 1, surf%ns

       i = surf%i(m)
       j = surf%j(m)
       k_atm = topo_top_ind(j,i,0) + 1
       k_topo = topo_top_ind(j,i,0)
!
!--    Calculate the pt, vpt and q for atmosphere depending on what modules are enabled.
       IF ( bulk_cloud_model )  THEN
          surf%pt1(m)  = pt(k_atm,j,i) + lv_d_cp * d_exner(k_atm) * ql(k_atm,j,i)
          surf%q1(m)   = q(k_atm,j,i) - ql(k_atm,j,i)
          surf%vpt1(m) = surf%pt1(m) * ( 1.0_wp + 0.61_wp * surf%q1(m) )
       ELSEIF ( cloud_droplets )  THEN
          surf%pt1(m)  = pt(k_atm,j,i) + lv_d_cp * d_exner(k_atm) * ql(k_atm,j,i)
          surf%q1(m)   = q(k_atm,j,i)
          surf%vpt1(m) = surf%pt1(m) * ( 1.0_wp + 0.61_wp * surf%q1(m) )
       ELSE
          surf%pt1(m) = pt(k_atm,j,i)
          IF ( moist_physics )  THEN
             surf%q1(m)   = q(k_atm,j,i)
             surf%vpt1(m) = surf%pt1(m) * ( 1.0_wp + 0.61_wp * surf%q1(m) )
          ENDIF
       ENDIF

       surf%uv_abs1(m) = SQRT( ( 0.5 * ( u(k_atm,j,i) + u(k_atm,j,i+1) ) )**2 +                    &
                               ( 0.5 * ( v(k_atm,j,i) + v(k_atm,j+1,i) ) )**2 )
       surf%uv_eff1(m) = surf%uv_abs1(m)

!
!--    Initialize the prognostic model variables and all the other variables where the value from
!--    the previous time step is used. For restart runs or if spinup data is available, these are
!--    read from the restart files, and should not be overwritten here.
       IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.  .NOT. read_spinup_data )   &
       THEN
!
!--       If spinup is enabled for current run, use diurnal mean spinup pt as the initial
!--       atmospheric boundary condition. Otherwise, use the first atmospheric grid level.
          IF ( spinup )  THEN
             bc_atm = spinup_pt_mean * exner(k_topo)
          ELSE
             bc_atm = surf%pt1(m) * exner(k_topo)
          ENDIF

          surf%us_urb(m)  = 1.0_wp
          surf%us_can(m)  = 1.0_wp
          surf%us_roof(m) = 1.0_wp
          surf%us_road(m) = 1.0_wp

          surf%uv_abs_can(m) = surf%uv_abs_can_coef(m) * surf%uv_abs1(m)
          surf%uv_eff_can(m) = surf%uv_abs_can(m)

          surf%shf_urb(m)  = 0.0_wp
          surf%qsws_urb(m) = 0.0_wp

          surf%t_can(m)   = bc_atm
          surf%t_can_p(m) = surf%t_can(m)

!
!--       For subsurface temps, a steady-state 1D heat equation solution will be used as the
!--       initial temperature profile. This might or might not speed up the spinup process.
!--       In case of windowless facade, set window temps to fill value to prevent meaningless
!--       output values. Vice versa for the opposite case.
          IF ( surf%f_win(m) < 1.0_wp )  THEN
             surf%t_wall_a(:,m) = calc_1d_heat_equation( SIZE( surf%t_wall_a, 1 ), bc_atm,         &
                                                         surf%t_indoor(m),                         &
                                                         surf%conductivity_wall(:,m) )
             surf%t_wall_a_p(:,m) = surf%t_wall_a(:,m)
             surf%t_wall_b(:,m)   = surf%t_wall_a(:,m)
             surf%t_wall_b_p(:,m) = surf%t_wall_a(:,m)
          ELSE
             IF ( .NOT. data_output_raw )  THEN
                surf%t_wall_a(:,m)   = output_fill_value
                surf%t_wall_a_p(:,m) = output_fill_value
             ENDIF
          ENDIF

          IF ( surf%f_win(m) > 0.0_wp )  THEN
             surf%t_win_a(:,m) = calc_1d_heat_equation( SIZE( surf%t_win_a, 1 ), bc_atm,           &
                                                        surf%t_indoor(m),                          &
                                                        surf%conductivity_win(:,m) )
             surf%t_win_a_p(:,m) = surf%t_win_a(:,m)
             surf%t_win_b(:,m)   = surf%t_win_a(:,m)
             surf%t_win_b_p(:,m) = surf%t_win_a(:,m)
          ELSE
             IF ( .NOT. data_output_raw )  THEN
                surf%t_win_a(:,m)   = output_fill_value
                surf%t_win_a_p(:,m) = output_fill_value
             ENDIF
          ENDIF

          surf%t_roof(:,m) = calc_1d_heat_equation( SIZE( surf%t_roof, 1 ), bc_atm,                &
                                                    surf%t_indoor(m), surf%conductivity_roof(:,m) )
          surf%t_roof_p(:,m) = surf%t_roof(:,m)

          surf%t_road(:,m) = calc_1d_heat_equation( SIZE( surf%t_road, 1 ), bc_atm,                &
                                                    surf%t_soil(m), surf%conductivity_road(:,m) )
          surf%t_road_p(:,m) = surf%t_road(:,m)

          IF ( moist_physics )  THEN
             surf%vpt_can(m) = 0.0_wp

             surf%q_can(m)        = surf%q1(m)
             surf%q_can_p(m)      = surf%q_can(m)
             surf%m_liq_roof(m)   = 0.0_wp
             surf%m_liq_roof_p(m) = surf%m_liq_roof(m)
             surf%m_liq_road(m)   = 0.0_wp
             surf%m_liq_road_p(m) = surf%m_liq_road(m)

             surf%q_roof(m) = surf%q1(m)
             surf%q_road(m) = surf%q_can(m)
          ENDIF

          surf%ol_roof(m) = surf%z_mo(m)     / zeta_min
          surf%ol_road(m) = surf%z_mo_can(m) / zeta_min
          surf%ol_can(m)  = surf%z_mo(m)     / zeta_min
          surf%ol_urb(m)  = surf%z_mo(m)     / zeta_min
       ENDIF

!
!--    Init potential temperatures and virtual potential temperatures. These need to be computed
!--    also for the restart case, as d_exner is not yet available when rrd routines are called.
       surf%pt_can(m) = surf%t_can(m) * d_exner(k_topo)

       IF ( surf%f_win(m) < 1.0_wp )  THEN
          surf%pt_wall_a(m) = surf%t_wall_a(nzt_wall,m) * d_exner(k_topo)
          surf%pt_wall_b(m) = surf%t_wall_b(nzt_wall,m) * d_exner(k_topo)
       ELSE
          IF ( .NOT. data_output_raw )  THEN
             surf%pt_wall_a(m) = output_fill_value
             surf%pt_wall_b(m) = output_fill_value
          ENDIF
       ENDIF

       IF ( surf%f_win(m) > 0.0_wp )  THEN
          surf%pt_win_a(m) = surf%t_win_a(nzt_win,m) * d_exner(k_topo)
          surf%pt_win_b(m) = surf%t_win_b(nzt_win,m) * d_exner(k_topo)
       ELSE
          IF ( .NOT. data_output_raw )  THEN
             surf%pt_win_a(m) = output_fill_value
             surf%pt_win_b(m) = output_fill_value
          ENDIF
       ENDIF

       surf%pt_roof(m) = surf%t_roof(nzt_roof,m) * d_exner(k_topo)
       surf%pt_road(m) = surf%t_road(nzt_road,m) * d_exner(k_topo)

       IF ( moist_physics )  THEN
          surf%vpt_can(m)  = surf%pt_can(m)  * ( 1.0_wp + 0.61_wp * surf%q_can(m)  )
          surf%vpt_roof(m) = surf%pt_roof(m) * ( 1.0_wp + 0.61_wp * surf%q_roof(m) )
          surf%vpt_road(m) = surf%pt_road(m) * ( 1.0_wp + 0.61_wp * surf%q_road(m) )
       ENDIF

!
!--    Initialize tendencies to zero.
       surf%tt_can(m)      = 0.0_wp
       surf%tt_wall_a(:,m) = 0.0_wp
       surf%tt_wall_b(:,m) = 0.0_wp
       surf%tt_win_a(:,m)  = 0.0_wp
       surf%tt_win_b(:,m)  = 0.0_wp
       surf%tt_roof(:,m)   = 0.0_wp
       surf%tt_road(:,m)   = 0.0_wp
       IF ( moist_physics )  THEN
          surf%tq_can(m)      = 0.0_wp
          surf%tm_liq_roof(m) = 0.0_wp
          surf%tm_liq_road(m) = 0.0_wp
       ENDIF

!
!--    Initialize model variables which are not used prior to an assignment in the model itself.
!--    Thus these initializations should not end up being used in code, but as this is not
!--    guaranteed with e.g. future changes, initialize them nevertheless. For the same reason,
!--    these are not included in the restart data. But if in the future there is an usage prior to
!--    proper assignment by the model, the respective variable should be added to restart routines,
!--    and given a proper intialization.
       surf%albedo_urb(m)     = 0.0_wp
       surf%emiss_urb(m)      = 1.0_wp
       surf%rad_lw_in_urb(m)  = 0.0_wp
       surf%rad_lw_out_urb(m) = 0.0_wp
       surf%rad_sw_in_urb(m)  = 0.0_wp
       surf%rad_sw_out_urb(m) = 0.0_wp
       surf%ram_urb(m)        = 1E3_wp
       surf%rib_urb(m)        = 0.0_wp
       surf%t_2m_urb(m)       = 0.0_wp
       surf%t_c_urb(m)        = 0.0_wp
       surf%t_h_urb(m)        = 0.0_wp
       surf%t_rad_urb(m)      = 0.0_wp
       surf%usws_urb(m)       = 0.0_wp
       surf%vsws_urb(m)       = 0.0_wp

       surf%shf_can(m)    = 0.0_wp
       surf%shf_road(m)   = 0.0_wp
       surf%shf_roof(m)   = 0.0_wp
       surf%shf_wall_a(m) = 0.0_wp
       surf%shf_wall_b(m) = 0.0_wp
       surf%shf_win_a(m)  = 0.0_wp
       surf%shf_win_b(m)  = 0.0_wp

       surf%ghf_road(m)   = 0.0_wp
       surf%ghf_roof(m)   = 0.0_wp
       surf%ghf_wall_a(m) = 0.0_wp
       surf%ghf_wall_b(m) = 0.0_wp
       surf%ghf_win_a(m)  = 0.0_wp
       surf%ghf_win_b(m)  = 0.0_wp

       surf%rad_lw_net_can(m)    = 0.0_wp
       surf%rad_lw_net_road(m)   = 0.0_wp
       surf%rad_lw_net_roof(m)   = 0.0_wp
       surf%rad_lw_net_urb(m)    = 0.0_wp
       surf%rad_lw_net_wall_a(m) = 0.0_wp
       surf%rad_lw_net_wall_b(m) = 0.0_wp
       surf%rad_lw_net_win_a(m)  = 0.0_wp
       surf%rad_lw_net_win_b(m)  = 0.0_wp
       surf%rad_sw_in_road(m)    = 0.0_wp
       surf%rad_sw_in_win_a(m)   = 0.0_wp
       surf%rad_sw_in_win_b(m)   = 0.0_wp
       surf%rad_sw_net_road(m)   = 0.0_wp
       surf%rad_sw_net_roof(m)   = 0.0_wp
       surf%rad_sw_net_urb(m)    = 0.0_wp
       surf%rad_sw_net_wall_a(m) = 0.0_wp
       surf%rad_sw_net_wall_b(m) = 0.0_wp
       surf%rad_sw_net_win_a (m) = 0.0_wp
       surf%rad_sw_net_win_b (m) = 0.0_wp

       surf%rib_can(m) = 0.0_wp
       surf%rib_road(m) = 0.0_wp
       surf%rib_roof(m) = 0.0_wp

       surf%rah_can(m)  = 1E3_wp
       surf%rah_road(m) = 1E3_wp
       surf%rah_roof    = 1E3_wp

       IF ( facade_rah_doe )  THEN
          surf%rah_wall_a(m) = 1E3_wp
          surf%rah_wall_b(m) = 1E3_wp
          surf%rah_win_a(m)  = 1E3_wp
          surf%rah_win_b(m)  = 1E3_wp
       ELSE
          surf%rah_facade(m) = 1E3_wp
       ENDIF

       IF ( moist_physics )  THEN
          surf%qsws_can(m)      = 0.0_wp
          surf%qsws_liq_road(m) = 0.0_wp
          surf%qsws_liq_roof(m) = 0.0_wp
          surf%qsws_road(m)     = 0.0_wp
          surf%qsws_roof(m)     = 0.0_wp

          e_s = 0.01_wp * magnus( MIN( surf%t_road(nzt_road,m), 333.15_wp ) )
          surf%qs_road(m) = rd_d_rv * e_s / ( surface_pressure - e_s )
          e_s = 0.01_wp * magnus( MIN( surf%t_roof(nzt_roof,m), 333.15_wp ) )
          surf%qs_roof(m) = rd_d_rv * e_s / ( surface_pressure - e_s )

          surf%c_liq_road(m)    = MIN( 1.0_wp, ( surf%m_liq_road(m) / m_liq_max_road )**0.67 )
          surf%c_liq_roof(m)    = MIN( 1.0_wp, ( surf%m_liq_roof(m) / m_liq_max_roof )**0.67 )
       ENDIF


    ENDDO

 END SUBROUTINE init_slurb_variables


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a gridded input INTEGER variable either from NAMELIST or netCDF4 input.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_grid_variable_1d_int( varname, var, var_namelist )

    CHARACTER(len=*), INTENT(IN) ::  varname  !< variable name in file

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  var_namelist  !< variable in namelist

    INTEGER(iwp), DIMENSION(:), INTENT(INOUT) ::  var  !< target variable

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  tmp  !< temporary array for 2D -> 1D array mapping


!
!-- First, try a NAMELIST initialization or use default if available. If the variable is not set
!-- in the namelist and does not have a proper default, the NAMELIST variable will be
!-- initialized to -9999.
    IF ( PRESENT( var_namelist ) )  THEN
       IF ( var_namelist /= -9999 )  THEN
          DO  m = 1, surf%ns
             var(m) = var_namelist
          ENDDO
       ENDIF
    ENDIF

!
!-- Third, try to read from a netCDF file.
    IF ( input_file_present )  THEN
       IF ( check_existence( var_names, varname ) )  THEN
!
!--       Allocate a temporary array matching the 2D PALM grid.
          ALLOCATE( tmp(nys:nyn,nxl:nxr) )
!
!--       Read the 2D input array to the temporary array.
          CALL get_variable( id_slurb, varname, tmp, nxl, nxr, nys, nyn, nbgp=0 )
!
!--       Map the temporary array into the SLUrb surface type.
          DO  m = 1, surf%ns
             i = surf%i(m)
             j = surf%j(m)
!
!--          Set the input into target array.
             var(m) = tmp(j,i)
!
!--          Missing values not allowed for the array for urban tiles -> give an error.
             IF ( var(m) == -9999 )  THEN
                WRITE( message_string, * ) 'Missing value for ' // TRIM( varname ) //              &
                                           ' for grid cell (j,i) = ', j, i,                        &
                                           ' with non-zero urban fraction in file ' //             &
                                           TRIM( input_file_slurb ) // '.'
                CALL message( 'slurb_init', 'SLU0021', 2, 2, myid, 6, 0 )
             ENDIF
          ENDDO
          DEALLOCATE( tmp )
       ENDIF
    ENDIF

 END SUBROUTINE get_grid_variable_1d_int


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a gridded input INTEGER variable either from NAMELIST or netCDF4 input.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_grid_variable_1d_real( varname, var, var_namelist )

    CHARACTER(len=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(wp), INTENT(IN), OPTIONAL ::  var_namelist  !< variable in namelist

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var  !< target variable

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tmp  !< temporary array for 2D -> 1D array mapping


!
!-- First, try a NAMELIST initialization or use default if available. If the variable is not set
!-- in the namelist and does not have a proper default, the NAMELIST variable will be
!-- initialized to -9999.
    IF ( PRESENT( var_namelist ) )  THEN
       IF ( var_namelist /= -9999.0_wp )  THEN
          DO  m = 1, surf%ns
             var(m) = var_namelist
          ENDDO
       ENDIF
    ENDIF

!
!-- Third, try to read from a netCDF file.
    IF ( input_file_present )  THEN
       IF ( check_existence( var_names, varname ) )  THEN
!
!--       Allocate a temporary array matching the 2D PALM grid.
          ALLOCATE( tmp(nys:nyn,nxl:nxr) )
!
!--       Read the 2D input array to the temporary array
          CALL get_variable( id_slurb, varname, tmp, nxl, nxr, nys, nyn, nbgp=0 )
!
!--       Map the temporary array into the SLUrb surface type.
          DO  m = 1, surf%ns
             i = surf%i(m)
             j = surf%j(m)
!
!--          Set the input into target array.
             var(m) = tmp(j,i)
!
!--          Missing values not allowed for the array for urban tiles -> give an error.
             IF ( var(m) == -9999.0_wp )  THEN
                WRITE( message_string, * ) 'Missing value for ' // TRIM( varname ) //              &
                                           ' for grid cell (j,i)=', j, i,                          &
                                           ' with non-zero urban fraction in file '                &
                                           // TRIM( input_file_slurb ) // '.'
                CALL message( 'slurb_init', 'SLU0022', 2, 2, myid, 6, 0 )
             ENDIF
          ENDDO
          DEALLOCATE( tmp )
       ENDIF
    ENDIF

 END SUBROUTINE get_grid_variable_1d_real


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a gridded and layered input variable either from NAMELIST or netCDF4 input.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_grid_variable_2d_real( varname, var, var_namelist  )

    CHARACTER(len=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(wp), DIMENSION(:), INTENT(IN), OPTIONAL ::  var_namelist  !< variable in namelist

    REAL(wp), DIMENSION(:,:), INTENT(INOUT)::  var  !< target variable

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp  !< temporary array for 2D -> 1D array mapping

    INTEGER(iwp) ::  k_target  !< target k-index


!
!-- First, try a NAMELIST initialization or use default if available. If the variable is not set
!-- in the namelist and does not have a proper default, the NAMELIST variable will be
!-- initialized to -9999.0_wp.
    IF ( PRESENT( var_namelist ) )  THEN
       IF ( ALL( var_namelist /= -9999.0_wp ) )  THEN
          DO  m = 1, surf%ns
             k_target = LBOUND( var_namelist, 1 )
             DO  k = LBOUND( var, 1 ), UBOUND( var, 1 )
                var(k,m) = var_namelist(k_target)
                k_target = k_target + 1
             ENDDO
          ENDDO
       ENDIF
    ENDIF

!
!-- Third, try to read from a netCDF file.
    IF ( input_file_present )  THEN
       IF ( check_existence( var_names, varname ) )  THEN
!
!--       Allocate a temporary array matching the 2D PALM grid.
          ALLOCATE( tmp(LBOUND(var, 1):UBOUND(var,1),nys:nyn,nxl:nxr) )
!
!--       Read the 2D input array to the temporary array.
          CALL get_variable( id_slurb, varname, tmp, nxl, nxr, nys, nyn, 0, SIZE( var, 1 )-1,      &
                             nbgp=0 )
!
!--       Map the temporary array into the SLUrb surface type.
          DO  m = 1, surf%ns
             i = surf%i(m)
             j = surf%j(m)
!
!--          Set the input into target array.
             DO  k = LBOUND( tmp, 1 ), UBOUND( tmp, 1 )
                IF ( tmp(k,j,i) /= -9999.0_wp )  var(k,m) = tmp(k,j,i)
!
!--             Missing values not allowed for the array for urban tiles -> give an error.
                IF ( var(k,m) == -9999.0_wp )  THEN
                   WRITE( message_string, * ) 'missing value for ' // TRIM( varname ) //           &
                                              ' for grid cell (j,i)=', surf%j(m), surf%i(m),       &
                                              ' with non-zero urban fraction in file ' //          &
                                              TRIM( input_file_slurb ) // '.'
                   CALL message( 'slurb_init', 'SLU0023', 2, 2, myid, 6, 0 )
                ENDIF
             ENDDO
          ENDDO
          DEALLOCATE( tmp )
       ENDIF
    ENDIF

 END SUBROUTINE get_grid_variable_2d_real


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes latent variables which can be inferred from the inputs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE precompute_latent_variables

    REAL(wp) ::  emiss_facade       !< aggregated facade emissivity
    REAL(wp) ::  f_wall             !< wall fraction
    REAL(wp) ::  f_win              !< window fraction
    REAL(wp) ::  wake               !< wake parameter for U_can parametrization in SURFEX
    REAL(wp) ::  win_nonrefl_1side  !< 1-side nonreflected radiation (for windows)
    REAL(wp) ::  win_absorp         !< window absorption coefficient


!
!-- Precompute model constants.
    rho_lv = l_v * rho_surface
    drho_l_lv = 1.0_wp / (rho_l * l_v)

!
!-- Check street canyon height-to-width ratio.
    DO  m = 1, surf%ns
       IF ( surf%h_bld(m) == 0.0_wp  .OR.  surf%hw_can(m) == 0.0_wp )  THEN
          i = surf%i(m)
          j = surf%j(m)
          WRITE( message_string, * ) 'Building height or street canyon aspect ratio ' //           &
                                     'set to zero at (j,i)=', j, i, 'with non-zero ' //            &
                                     'urban_fraction = ', fr_urb(j,i), '.'
          CALL message( 'slurb_init', 'SLU0024', 2, 2, myid, 6, 0 )
       ENDIF
    ENDDO

!
!-- Precompute layer total conductivities from layer thicknesses and thermal conductivities.
    DO  m = 1, surf%ns
       DO  k = nzt_roof, nzb_roof-1
          surf%conductivity_roof(k,m) = 2.0_wp / ( surf%dz_roof(k,m)   / surf%lambda_roof(k,m) +   &
                                                   surf%dz_roof(k+1,m) / surf%lambda_roof(k+1,m) )
       ENDDO
       surf%conductivity_roof(nzb_roof,m) = 2.0_wp * surf%lambda_roof(nzb_roof,m) /                &
                                            surf%dz_roof(nzb_roof,m)

       DO  k = nzt_wall, nzb_wall-1
          surf%conductivity_wall(k,m) = 2.0_wp / ( surf%dz_wall(k,m)   / surf%lambda_wall(k,m) +   &
                                                   surf%dz_wall(k+1,m) / surf%lambda_wall(k+1,m) )
       ENDDO
       surf%conductivity_wall(nzb_wall,m) = 2.0_wp * surf%lambda_wall(nzb_wall,m) /                &
                                            surf%dz_wall(nzb_wall,m)

       DO  k = nzt_win, nzb_win-1
          surf%conductivity_win(k,m) = 2.0_wp / ( surf%dz_win(k,m)   / surf%lambda_win(k,m) +      &
                                                  surf%dz_win(k+1,m) / surf%lambda_win(k+1,m) )
       ENDDO
       surf%conductivity_win(nzb_wall,m) = 2.0_wp * surf%lambda_wall(nzb_win,m) /                  &
                                           surf%dz_wall(nzb_win,m)

!
!--    For the road, the last conductance depends on the soil conductance, so we need to
!--    compute it during time-stepping (as it depends on soil moisture).
       DO  k = nzt_road, nzb_road-1
          surf%conductivity_road(k,m) = 2.0_wp / ( surf%dz_road(k,m)   / surf%lambda_road(k,m) +   &
                                                   surf%dz_road(k+1,m) / surf%lambda_road(k+1,m) )
       ENDDO
       surf%conductivity_road(nzb_road,m) = 2.0_wp * surf%lambda_road(nzb_road,m) /                &
                                            surf%dz_road(nzb_road,m)
    ENDDO

!
!-- Precompute sky-view factors.
    surf%svf_road(:) = 0.0_wp
    surf%svf_wall(:) = 0.0_wp

    DO  m = 1, surf%ns
       surf%svf_road(m) = SQRT( (surf%hw_can(m))**2 + 1.0_wp ) - surf%hw_can(m)

       surf%svf_wall(m) = ( surf%hw_can(m) + 1.0_wp - SQRT( surf%hw_can(m)**2 + 1.0_wp ) ) /       &
                          ( 2.0_wp * surf%hw_can(m) )
    ENDDO

!
!-- Precompute urban emissivity based on SVFs.
    DO  m = 1, surf%ns
       surf%emiss_urb(m) = surf%f_bld(m) * surf%emiss_roof(m)                                      &
                           + ( 1.0_wp - surf%f_bld(m) )                                            &
                              * ( surf%svf_road(m) * surf%emiss_road(m)                            &
                                + surf%svf_wall(m) * surf%hw_can(m)                                &
                                  * ( ( 1.0_wp - surf%f_win(m) ) * 2.0_wp * surf%emiss_wall(m)     &
                                      + surf%f_win(m) * 2.0_wp * surf%emiss_win(m) ) )
    ENDDO

!
!-- Preompute the longwave interaction coefficients for surface elements as these are
!-- static in time. Based on Johnson et al. (1991) general formula. Absorption from reflected
!-- radiation is taken into account only after first reflection. The first reflections contribute
!-- around 5% of the total LW budget, while higher order reflections would contribute only <0.5%.
!-- Coefficients are grouped per variable, so they can be effectively used in time-stepping
!-- without wasting too much computational time or memory.
    DO  m = 1, surf%ns
!
!--    Compute aggregated facade emissivity for simplification of reflections.
       f_win  = surf%f_win(m)
       f_wall = ( 1.0_wp - f_win )
       emiss_facade = f_wall * surf%emiss_wall(m) + f_win * surf%emiss_win(m)
!
!--    Roof.
!--    To be multiplied by t_roof**4 in the LW budget:
       surf%lw_roof_coef(1,m) = -surf%emiss_roof(m) * sigma_sb
!
!--    To be multiplied by lw_rad_in_urb in the LW budget:
       surf%lw_roof_coef(2,m) = surf%emiss_roof(m)
!
!--    Roads.
!--    To be multiplied by t_road**4 in the LW budget:
       surf%lw_road_coef(1,m) = ( - surf%emiss_road(m)                                             &
                                  + surf%emiss_road(m)**2 * ( 1 - emiss_facade )                   &
                                    * ( 1.0_wp - surf%svf_road(m) ) * surf%svf_wall(m)             &
                                ) * sigma_sb
!
!--    To be multiplied by lw_rad_in in the LW budget:
       surf%lw_road_coef(2,m) = surf%emiss_road(m) * surf%svf_road(m)                              &
                               - surf%emiss_road(m) * ( 1.0_wp - emiss_facade ) * surf%svf_wall(m) &
                                 * ( 1.0_wp - surf%svf_road(m) )
!
!--    To be multiplied by (t_wall_a**4 + t_wall_b**4) in the LW budget:
       surf%lw_road_coef(3,m) = ( surf%emiss_road(m) * surf%emiss_wall(m)                          &
                                    * ( 1.0_wp - surf%svf_road(m) )                                &
                                  + surf%emiss_road(m) * surf%emiss_wall(m)                        &
                                    * ( 1.0_wp - emiss_facade )                                    &
                                    * ( 1.0_wp - surf%svf_road(m) )                                &
                                    * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                       &
                                ) * 0.5_wp * f_wall * sigma_sb
       surf%lw_road_coef(4,m) = ( surf%emiss_road(m) * surf%emiss_win(m)                           &
                                    * ( 1.0_wp - surf%svf_road(m) )                                &
                                  + surf%emiss_road(m) * surf%emiss_win(m)                         &
                                    * ( 1.0_wp - emiss_facade )                                    &
                                    * ( 1.0_wp - surf%svf_road(m) )                                &
                                    * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                       &
                                ) * 0.5_wp * f_win * sigma_sb
!
!--    Walls.
!--    To be multiplied by t_wall_a**4 in the LW budget:
       surf%lw_wall_coef(1,m) = ( - surf%emiss_wall(m)                                             &
                                 + 0.5_wp * f_wall * surf%emiss_wall(m)**2                         &
                                   * ( 1.0_wp - surf%emiss_road(m) ) * surf%svf_wall(m)            &
                                   * ( 1.0_wp - surf%svf_road(m) )                                 &
                                 + f_wall * surf%emiss_wall(m)**2                                  &
                                   * ( 1.0_wp - emiss_facade )                                     &
                                   * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )**2                     &
                                ) * sigma_sb
!
!--    To be multiplied by lw_rad_in in the LW budget:
       surf%lw_wall_coef(2,m) = surf%emiss_wall(m) * surf%svf_wall(m)                              &
                                + surf%emiss_wall(m) * ( 1.0_wp - surf%emiss_road(m) )             &
                                  * surf%svf_wall(m) * surf%svf_road(m)                            &
                                + surf%emiss_wall(m)                                               &
                                  * ( 1.0_wp - emiss_facade )                                      &
                                  * surf%svf_wall(m) * surf%svf_road(m)                            &
                                + surf%emiss_wall(m)                                               &
                                  * ( 1.0_wp - emiss_facade )                                      &
                                  * surf%svf_wall(m)                                               &
                                  * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )
!
!--    To be multiplied by t_wall_b**4 in the LW budget:
       surf%lw_wall_coef(3,m) = ( 0.5_wp * f_wall * surf%emiss_wall(m)**2                          &
                                  * ( 1.0_wp - surf%emiss_road(m) )                                &
                                  * surf%svf_wall(m)                                               &
                                  * ( 1.0_wp - surf%svf_road(m) )                                  &
                                    + f_wall * surf%emiss_wall(m)**2                               &
                                      * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                     &
                                ) * sigma_sb
!
!--    To be multiplied by t_win_a**4 in the LW budget:
       surf%lw_wall_coef(4,m) = ( 0.5_wp * f_win * surf%emiss_wall(m)                              &
                                    * surf%emiss_win(m)                                            &
                                    * ( 1.0_wp - surf%emiss_road(m) )                              &
                                    * surf%svf_wall(m)                                             &
                                    * ( 1.0_wp - surf%svf_road(m) )                                &
                                  + f_win * surf%emiss_wall(m)                                     &
                                     * surf%emiss_win(m)                                           &
                                     * ( 1.0_wp - emiss_facade )                                   &
                                     * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )**2                   &
                                ) * sigma_sb
!
!--    To be multiplied by t_win_b**4 in the LW budget:
       surf%lw_wall_coef(5,m) = ( 0.5_wp * f_win * surf%emiss_wall(m)                              &
                                     * surf%emiss_win(m)                                           &
                                     * ( 1.0_wp - surf%emiss_road(m) )                             &
                                     * surf%svf_wall(m)                                            &
                                     * ( 1.0_wp - surf%svf_road(m) )                               &
                                  + f_win * surf%emiss_wall(m)                                     &
                                     * surf%emiss_win(m)                                           &
                                     * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                      &
                                ) * sigma_sb
!
!--    To be multiplied by t_road**4 in the LW budget:
       surf%lw_wall_coef(6,m) = ( surf%emiss_wall(m) * surf%emiss_road(m)                          &
                                     * surf%svf_wall(m)                                            &
                                  + surf%emiss_wall(m) * surf%emiss_road(m)                        &
                                     * ( 1.0_wp - emiss_facade )                                   &
                                     * surf%svf_wall(m)                                            &
                                     * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                      &
                                ) * sigma_sb
!
!--    Windows.
!--    To be multiplied by t_wall_a**4 in the LW budget:
       surf%lw_win_coef(1,m) = ( - surf%emiss_win(m)                                               &
                                 + 0.5_wp * f_win * surf%emiss_win(m)**2                           &
                                   * ( 1.0_wp - surf%emiss_road(m) )                               &
                                   * surf%svf_wall(m)                                              &
                                   * ( 1.0_wp - surf%svf_road(m) )                                 &
                                 + f_win * surf%emiss_win(m)**2                                    &
                                   * ( 1.0_wp - emiss_facade )                                     &
                                   * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )**2                     &
                               ) * sigma_sb
!
!--    To be multiplied by lw_rad_in in the LW budget:
       surf%lw_win_coef(2,m) = surf%emiss_win(m) * surf%svf_wall(m)                                &
                               + surf%emiss_win(m)                                                 &
                                 * ( 1.0_wp - surf%emiss_road(m) )                                 &
                                 * surf%svf_wall(m) * surf%svf_road(m)                             &
                               + surf%emiss_win(m)                                                 &
                                 * ( 1.0_wp - emiss_facade )                                       &
                                 * surf%svf_wall(m) * surf%svf_road(m)                             &
                               + surf%emiss_win(m)                                                 &
                                 * ( 1.0_wp - emiss_facade )                                       &
                                 * surf%svf_wall(m)                                                &
                                 * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )
!
!--    To be multiplied by t_win_b**4 in the LW budget:
       surf%lw_win_coef(3,m) = ( 0.5_wp * f_win * surf%emiss_win(m)**2                             &
                                   * ( 1.0_wp - surf%emiss_road(m) )                               &
                                   * surf%svf_wall(m)                                              &
                                   * ( 1.0_wp - surf%svf_road(m) )                                 &
                                 + f_win * surf%emiss_win(m)**2                                    &
                                   * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                        &
                               ) * sigma_sb
!
!--    To be multiplied by t_wall_a**4 in the LW budget:
       surf%lw_win_coef(4,m) = ( 0.5_wp * f_wall * surf%emiss_win(m)                               &
                                   * surf%emiss_wall(m)                                            &
                                   * ( 1.0_wp - surf%emiss_road(m) )                               &
                                   * surf%svf_wall(m)                                              &
                                   * ( 1.0_wp - surf%svf_road(m) )                                 &
                                 + f_wall * surf%emiss_win(m)                                      &
                                   * surf%emiss_wall(m)                                            &
                                   * ( 1.0_wp - emiss_facade )                                     &
                                   * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )**2                     &
                               ) * sigma_sb
!
!--    To be multiplied by t_wall_b**4 in the LW budget:
       surf%lw_win_coef(5,m) = ( 0.5_wp * f_wall * surf%emiss_win(m)                               &
                                   * surf%emiss_wall(m)                                            &
                                   * ( 1.0_wp - surf%emiss_road(m) )                               &
                                   * surf%svf_wall(m)                                              &
                                   * ( 1.0_wp - surf%svf_road(m) )                                 &
                                 + f_wall * surf%emiss_win(m)                                      &
                                   * surf%emiss_wall(m)                                            &
                                   * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                        &
                               ) * sigma_sb
!
!--    To be multiplied by t_road**4 in the LW budget:
       surf%lw_win_coef(6,m) = ( surf%emiss_win(m) * surf%emiss_road(m)                            &
                                   * surf%svf_wall(m)                                              &
                                 + surf%emiss_win(m) * surf%emiss_road(m)                          &
                                   * ( 1.0_wp - emiss_facade )                                     &
                                   * surf%svf_wall(m)                                              &
                                   * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )                        &
                               ) * sigma_sb

    ENDDO

!
!-- Precompute shortwave radiation reflection denominator.
    DO  m = 1, surf%ns
       surf%sw_ref_denom(m) = 1.0_wp - surf%albedo_road(m)                                         &
                             * surf%albedo_wall_win(m)                                             &
                             * surf%svf_wall(m) * ( 1.0_wp - surf%svf_road(m) )                    &
                             - surf%albedo_wall_win(m) * ( 1.0_wp - 2.0_wp * surf%svf_wall(m) )
    ENDDO

!
!-- Compute window layer shortwave absorption based on USM documentation.
!-- @todo This computation needs checking. Now sw_transmitted is not simply equal to
!-- sw_net_win*transmissivity, as 1.0_wp - SUM(absorption(:,m)) != transmissivity(m). This is
!-- mitigated for now at the output side. No side effects for the model, as the transmitted
!-- radiation is purely an output.
    DO  m = 1, surf%ns
       win_nonrefl_1side = 1.0 - (surf%albedo_win(m) + surf%transmissivity_win(m)                  &
                           + 1.0_wp  - SQRT( ( surf%albedo_win(m)                                  &
                           + surf%transmissivity_win(m) + 1.0_wp )**2                              &
                           - 4.0_wp * surf%albedo_win(m) ) ) / 2.0_wp

       win_absorp = -LOG( ( surf%transmissivity_win(m) + surf%albedo_win(m)                        &
                            - 1.0_wp + win_nonrefl_1side ) / win_nonrefl_1side                     &
                        ) / surf%zw_win(nzb_win,m)

       DO  k = nzt_win, nzb_win
          IF ( k /= nzt_win)  THEN
!
!--          The absorbed fraction is difference between cumulative absorption over the layer.
             surf%absorption_win(k,m) = win_nonrefl_1side                                          &
                                        * ( EXP( -win_absorp * surf%zw_win(k-1,m) ) -              &
                                            EXP( -win_absorp * surf%zw_win(k,m)   ) )
          ELSE
!
!--          For the first layer, it is the cumulative absorption so far.
             surf%absorption_win(k,m) = win_nonrefl_1side *                                        &
                                        ( 1.0_wp - EXP( -win_absorp * surf%zw_win(k,m) ) )
          ENDIF
       ENDDO
    ENDDO

!
!-- Coefficient for the canyon wind speed Krayenhoff & Voogt (2007) Eq. (9).
    IF ( uv_can_factor_kray )  THEN
       DO  m = 1, surf%ns
          surf%uv_abs_can_coef(m) = LOG( surf%h_bld(m)  / ( 3.0_wp * surf%z0_urb(m) ) ) /          &
                               LOG( ( surf%z_mo(m) + surf%h_bld(m) / 3.0_wp ) / surf%z0_urb(m) ) * &
                               EXP( -surf%f_bld_frn(m) / ( 2.0_wp * ( 1.0_wp - surf%f_bld(m) ) ) )
       ENDDO
!
!-- Coefficient for the canyon windspeed as derived in Masson (2000) (original TEB)
    ELSEIF ( uv_can_factor_masson )  THEN
       DO  m = 1, surf%ns
          surf%uv_abs_can_coef(m) = LOG( surf%h_bld(m)  / ( 3.0_wp * surf%z0_urb(m) ) ) /          &
                               LOG( ( surf%z_mo(m) + surf%h_bld(m) / 3.0_wp ) / surf%z0_urb(m) ) * &
                               EXP( -surf%hw_can(m) / 4.0_wp )
       ENDDO
!
!-- Coefficient for the canyon windspeed as implemented in SURFEX v8.1
    ELSEIF ( uv_can_factor_surfex )  THEN
       DO  m = 1, surf%ns
          wake = 1.0_wp + ( 2.0_wp / pi - 1.0_wp ) * 2.0 * ( surf%hw_can(m) - 0.5_wp )
          wake = MAX( MIN( wake, 1.0_wp ), 2.0_wp / pi )
          surf%uv_abs_can_coef(m) = wake * EXP( - surf%hw_can(m) / 4.0_wp ) *                      &
                                    LOG( 2.0_wp * surf%h_bld(m) / ( 3.0_wp * surf%z0_urb(m) ) ) /  &
                                    LOG( ( surf%z_mo(m) + 2.0_wp * surf%h_bld(m) ) /               &
                                         ( 3.0_wp * surf%z0_urb(m) ) )
       ENDDO
    ENDIF

!
!-- Compute minimum timestep based on SLUrb internal diffusivities.
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)

!
!--    Criterion based on subsurface heat diffusion. Heat capacities are already multiplied with
!--    layer thickness and conductivity is already divided with it. Thus, no need to multiply with
!--    dz**2 like in usm and lsm. Dimension analysis:
!--    c [J m-2 K-1] and layer conductivity [W m-2 K-1] -> [J/W] -> [s]
       surf%dt_max(m) = MIN( surf%dt_max(m),                                                       &
                             MINVAL( surf%c_roof(:,m) / surf%conductivity_roof(:,m) ) )

       surf%dt_max(m) = MIN( surf%dt_max(m),                                                       &
                             MINVAL( surf%c_road(:,m) / surf%conductivity_road(:,m) ) )

       surf%dt_max(m) = MIN( surf%dt_max(m),                                                       &
                             MINVAL( surf%c_wall(:,m) / surf%conductivity_wall(:,m) ) )

       IF ( surf%f_win(m) > 0.0_wp )  THEN
          surf%dt_max(m) = MIN( surf%dt_max(m),                                                    &
                                MINVAL( surf%c_win(:,m) / surf%conductivity_win(:,m) ) )
       ENDIF

    ENDDO
!
!-- Consider a pre-factor (1/8) for the diffusion criterion.
    dt_slurb = MINVAL( surf%dt_max ) * 0.125_wp
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_slurb, 1, MPI_REAL, MPI_MIN, comm2d, ierr )
#endif

 END SUBROUTINE precompute_latent_variables


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Processes inputs that are dynamic in time. @todo There is room for improvement, e.g. adapt
!> the pre-existing get_grid_variable_* routines to handle time dimension as well. However, at
!> this point the benefit is not that big as there are only three dynamic inputs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE process_dynamic_inputs

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp  !<


!
!-- First, apply a default NAMELIST initialization (lod=0).
!
!-- Anthropogenic shf from industry/heating systems/etc. that enters the urban aggregated flux.
    surf%shf_external(:) = shf_external

!
!-- Anthropogenic qsws from industry/heating systems/etc. that enters the urban aggregated flux.
    IF ( moist_physics )  surf%qsws_external(:) = qsws_external

!
!-- Anthropogenic shf from traffic to the street canyon.
    surf%shf_traffic(:) = shf_traffic

    IF ( input_file_present )  THEN
!
!--    Check if there is time dimension in the input driver.
       IF ( check_existence( var_names, 'time' ) ) THEN
          CALL get_dimension_length( id_slurb, slurb_dynamic%ntime, 'time' )
          CALL get_variable( id_slurb, 'time', slurb_dynamic%time )

!
!--       If the input doesn't cover the whole simulation period, give an error.
          IF ( end_time - spinup_time > slurb_dynamic%time(slurb_dynamic%ntime) )  THEN
             WRITE( message_string, * ) 'Time dimension of the SLUrb input driver does not ' //    &
                                        'cover the entire simulation period.'
             CALL message( 'slurb_init', 'SLU0025', 2, 2, 0, 6, 0 )
!
!--       Warn if the input doesn't cover the spinup period.
          ELSEIF ( -spinup_time < slurb_dynamic%time(1) )  THEN
             WRITE( message_string, * ) 'Time dimension of the SLUrb input driver does not ' //    &
                            'cover the spinup period (negative reference time). The first input'// &
                            'values will be used for temporally dynamic inputs during the spinup.'
             CALL message( 'slurb_init', 'SLU0026', 0, 1, 0, 6, 0 )
!
!--       Raise an error if the temporal dimension doesn't cover the reference time.
          ELSEIF ( 0.0_wp < slurb_dynamic%time(1) )  THEN
             WRITE( message_string, * ) 'Time dimension of the SLUrb input driver does not ' //    &
                                        'cover the reference initialization time (0.0 seconds).'
             CALL message( 'slurb_init', 'SL0027', 2, 2, 0, 6, 0 )
          ENDIF

!
!--       Process shf_external if found from the driver.
          IF ( check_existence( var_names, 'shf_external' ) )  THEN
!
!--          Check lod of the input for dimensionality: 1 = (time), 2 = (time,j,i)
             CALL get_attribute( id_slurb, 'lod', slurb_dynamic%shf_external%lod,                  &
                                 .FALSE., 'shf_external' )
             IF ( slurb_dynamic%shf_external%lod == 1 )  THEN
!
!--             Allocate and read the 1D time series.
                ALLOCATE( slurb_dynamic%shf_external%var1d(1:slurb_dynamic%ntime) )
                CALL get_variable( id_slurb, 'shf_external', slurb_dynamic%shf_external%var1d )
!
!--             Check for missing values.
                IF ( ANY( slurb_dynamic%shf_external%var1d(:) == -9999.0_wp ) ) THEN
                   WRITE( message_string, * ) 'Missing value found in input time series' //        &
                                              'shf_external in file ' // TRIM( input_file_slurb )  &
                                              // '.'
                   CALL message( 'slurb_init', 'SLU0028', 2, 2, 0, 6, 0 )
                ENDIF

             ELSEIF ( slurb_dynamic%shf_external%lod == 2 )  THEN
!
!--             Allocate and read the 2D spatiotemporal data.
                ALLOCATE( slurb_dynamic%shf_external%var2d(1:slurb_dynamic%ntime,1:surf%ns) )
!
!--             Allocate a temporary array matching the 2D PALM grid.
                ALLOCATE( tmp(1:slurb_dynamic%ntime,nys:nyn,nxl:nxr) )
!
!--             Read the 2D input array to the temporary array.
                CALL get_variable( id_slurb, 'shf_external', tmp, nxl, nxr, nys, nyn )

                DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
!
!--             Set the input into target array.
                   DO  k = 1, slurb_dynamic%ntime
                      IF ( tmp(k,j,i) /= -9999.0_wp )  THEN
                         slurb_dynamic%shf_external%var2d(k,m) = tmp(k,j,i)
                      ELSE
                         WRITE( message_string, * ) 'Missing value for shf_external ' //           &
                                                    ' for grid cell (j,i) = ',                     &
                                                    surf%j(m), surf%i(m),                          &
                                                    ' with non-zero urban fraction in file ' //    &
                                                    TRIM( input_file_slurb ) // '.'
                         CALL message( 'slurb_init', 'SL0029', 2, 2, myid, 6, 0 )
                      ENDIF
                   ENDDO
                ENDDO
                DEALLOCATE( tmp )
             ELSE
                WRITE( message_string, * ) 'Invalid lod attribute for input variable ' //          &
                                           'shf_external in file '// TRIM( input_file_slurb ) //   &
                                           ': lod=', slurb_dynamic%shf_external%lod,               &
                                           ', expected lod=1 or lod=2.'
                CALL message( 'slurb_init', 'SLU0030', 2, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Process qsws_external if found from the driver.
          IF ( moist_physics )  THEN
             IF ( check_existence( var_names, 'qsws_external' ) )  THEN
!
!--             Check lod of the input for dimensionality: 1 = (time), 2 = (time,j,i)
                CALL get_attribute( id_slurb, 'lod', slurb_dynamic%qsws_external%lod,              &
                                    .FALSE., 'qsws_external' )
                IF ( slurb_dynamic%qsws_external%lod == 1 )  THEN
!
!--                Allocate and read the 1D time series.
                   ALLOCATE( slurb_dynamic%qsws_external%var1d(1:slurb_dynamic%ntime) )
                   CALL get_variable( id_slurb, 'qsws_external', slurb_dynamic%qsws_external%var1d )
!
!--                Check for missing values.
                   IF ( ANY( slurb_dynamic%qsws_external%var1d(:) == -9999.0_wp ) ) THEN
                      WRITE( message_string, * ) 'Missing value found in input time series' //     &
                                                 'qsws_external in file ' //                       &
                                                 TRIM( input_file_slurb ) // '.'
                      CALL message( 'slurb_init', 'SLU0028', 2, 2, 0, 6, 0 )
                   ENDIF
                ELSEIF ( slurb_dynamic%qsws_external%lod == 2 )  THEN
!
!--                Allocate and read the 2D spatiotemporal data.
                   ALLOCATE( slurb_dynamic%qsws_external%var2d(1:slurb_dynamic%ntime,1:surf%ns) )
!
!--                Allocate a temporary array matching the 2D PALM grid.
                   ALLOCATE( tmp(1:slurb_dynamic%ntime,nys:nyn,nxl:nxr) )
!
!--                Read the 2D input array to the temporary array.
                   CALL get_variable( id_slurb, 'qsws_external', tmp, nxl, nxr, nys, nyn )

                   DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
!
!--                Set the input into target array.
                      DO  k = 1, slurb_dynamic%ntime
                         IF ( tmp(k,j,i) /= -9999.0_wp )  THEN
                            slurb_dynamic%qsws_external%var2d(k,m) = tmp(k,j,i)
                         ELSE
                            WRITE( message_string, * ) 'missing value for qsws_external ' //       &
                                                       ' for grid cell (j,i)=',                    &
                                                       surf%j(m), surf%i(m),                       &
                                                       ' with non-zero urban fraction in file '    &
                                                     // TRIM( input_file_slurb )
                            CALL message( 'slurb_init', 'SLU0029', 2, 2, myid, 6, 0 )
                         ENDIF
                      ENDDO
                   ENDDO
                   DEALLOCATE( tmp )
                ELSE
                   WRITE( message_string, * ) 'Invalid lod attribute for input variable ' //       &
                                              'qsws_external in file '//                           &
                                              TRIM( input_file_slurb ) //                          &
                                              ': lod=', slurb_dynamic%qsws_external%lod,           &
                                              ', expected lod=1 or lod=2.'
                   CALL message( 'slurb_init', 'SLU0030', 2, 2, 0, 6, 0 )
                ENDIF
             ENDIF
          ENDIF
!
!--       Process shf_traffic if found from the driver.
          IF ( check_existence( var_names, 'shf_traffic' ) )  THEN
!
!--          Check lod of the input for dimensionality: 1 = (time), 2 = (time,j,i)
             CALL get_attribute( id_slurb, 'lod', slurb_dynamic%shf_traffic%lod,                   &
                                 .FALSE., 'shf_traffic' )
             IF ( slurb_dynamic%shf_traffic%lod == 1 )  THEN
!
!--             Allocate and read the 1D time series.
                ALLOCATE( slurb_dynamic%shf_traffic%var1d(1:slurb_dynamic%ntime) )
                CALL get_variable( id_slurb, 'shf_traffic', slurb_dynamic%shf_traffic%var1d )
!
!--             Check for missing values.
                IF ( ANY( slurb_dynamic%shf_traffic%var1d(:) == -9999.0_wp ) ) THEN
                   WRITE( message_string, * ) 'Missing value found in input time series' //        &
                                              'shf_traffic in file ' // TRIM( input_file_slurb ) //&
                                              '.'
                   CALL message( 'slurb_init', 'SLU0028', 2, 2, 0, 6, 0 )
                ENDIF
             ELSEIF ( slurb_dynamic%shf_traffic%lod == 2 )  THEN
!
!--             Allocate and read the 2D spatiotemporal data.
                ALLOCATE( slurb_dynamic%shf_traffic%var2d(1:slurb_dynamic%ntime,1:surf%ns) )
!
!--             Allocate a temporary array matching the 2D PALM grid.
                ALLOCATE( tmp(1:slurb_dynamic%ntime,nys:nyn,nxl:nxr) )
!
!--             Read the 2D input array to the temporary array.
                CALL get_variable( id_slurb, 'shf_traffic', tmp, nxl, nxr, nys, nyn )

                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
!
!--                Set the input into target array.
                   DO  k = 1, slurb_dynamic%ntime
                      IF ( tmp(k,j,i) /= -9999.0_wp )  THEN
                         slurb_dynamic%shf_traffic%var2d(k,m) = tmp(k,j,i)
                      ELSE
                         WRITE( message_string, * ) 'missing value for shf_traffic ' //            &
                                                    ' for grid cell (j,i)=', surf%j(m), surf%i(m), &
                                                    ' with non-zero urban fraction in file ' //    &
                                                    TRIM( input_file_slurb ) // '.'
                         CALL message( 'slurb_init', 'SLU0029', 2, 2, myid, 6, 0 )
                      ENDIF
                   ENDDO
                ENDDO
                DEALLOCATE( tmp )
             ELSE
                WRITE( message_string, * ) 'Invalid lod attribute for input variable ' //          &
                                           'shf_traffic in file '// TRIM( input_file_slurb ) //    &
                                           ': lod=', slurb_dynamic%shf_traffic%lod,                &
                                           ', expected lod=1 or lod=2.'
                CALL message( 'slurb_init', 'SLU0030', 2, 2, 0, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
    ENDIF
!
!-- Convert units for all inputs and scale the shf_traffic to street area.
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k_topo = topo_top_ind(j,i,0)
       surf%shf_external(m) = surf%shf_external(m) / fr_urb(j,i) * heatflux_input_conversion(k_topo)
       surf%shf_traffic(m)  = surf%shf_traffic(m) / ( ( 1.0_wp - surf%f_bld(m) ) * fr_urb(j,i) ) * &
                              heatflux_input_conversion(k_topo)
    ENDDO
    IF ( moist_physics )  THEN
       DO  m = 1, surf%ns
          surf%qsws_external(m) = surf%qsws_external(m) * waterflux_input_conversion(k_topo)
       ENDDO
    ENDIF

 END SUBROUTINE process_dynamic_inputs


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Process parameters dependent on the building/pavement type and properties.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE process_surface_parameters

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_tmp  !< array to contain building type temporarily


    ALLOCATE( type_tmp(1:surf%ns) )
    type_tmp(:) = -9999

    CALL slurb_default_pars
    CALL get_grid_variable_1d_int( 'building_type', type_tmp, building_type )
    CALL check_grid_variable_1d_int( 'building_type', type_tmp, 1, 6 )

    DO  m = 1, surf%ns
       surf%f_win(m) = building_pars_slurb(0,type_tmp(m))

       IF ( n_layers_roofs == 4 )  THEN
          surf%dz_roof(1,m) = building_pars_slurb(1,type_tmp(m))
          surf%dz_roof(2,m) = building_pars_slurb(2,type_tmp(m))
          surf%dz_roof(3,m) = building_pars_slurb(3,type_tmp(m))
          surf%dz_roof(4,m) = building_pars_slurb(4,type_tmp(m))

          surf%c_roof(1,m) = building_pars_slurb(5,type_tmp(m))
          surf%c_roof(2,m) = building_pars_slurb(6,type_tmp(m))
          surf%c_roof(3,m) = building_pars_slurb(7,type_tmp(m))
          surf%c_roof(4,m) = building_pars_slurb(8,type_tmp(m))

          surf%lambda_roof(1,m) = building_pars_slurb(9,type_tmp(m))
          surf%lambda_roof(2,m) = building_pars_slurb(10,type_tmp(m))
          surf%lambda_roof(3,m) = building_pars_slurb(11,type_tmp(m))
          surf%lambda_roof(4,m) = building_pars_slurb(12,type_tmp(m))
       ENDIF

       surf%z0_roof(m)     = building_pars_slurb(13,type_tmp(m))
       surf%z0h_roof(m)    = building_pars_slurb(13,type_tmp(m)) * 1.0E-2
       surf%albedo_roof(m) = building_pars_slurb(14,type_tmp(m))
       surf%emiss_roof(m)  = building_pars_slurb(15,type_tmp(m))

       IF ( n_layers_walls == 4 )  THEN
          surf%dz_wall(1,m) = building_pars_slurb(16,type_tmp(m))
          surf%dz_wall(2,m) = building_pars_slurb(17,type_tmp(m))
          surf%dz_wall(3,m) = building_pars_slurb(18,type_tmp(m))
          surf%dz_wall(4,m) = building_pars_slurb(19,type_tmp(m))

          surf%c_wall(1,m) = building_pars_slurb(20,type_tmp(m))
          surf%c_wall(2,m) = building_pars_slurb(21,type_tmp(m))
          surf%c_wall(3,m) = building_pars_slurb(22,type_tmp(m))
          surf%c_wall(4,m) = building_pars_slurb(23,type_tmp(m))

          surf%lambda_wall(1,m) = building_pars_slurb(24,type_tmp(m))
          surf%lambda_wall(2,m) = building_pars_slurb(25,type_tmp(m))
          surf%lambda_wall(3,m) = building_pars_slurb(26,type_tmp(m))
          surf%lambda_wall(4,m) = building_pars_slurb(27,type_tmp(m))
       ENDIF

       surf%z0_wall(m)     = building_pars_slurb(28,type_tmp(m))
       surf%albedo_wall(m) = building_pars_slurb(29,type_tmp(m))
       surf%emiss_wall(m)  = building_pars_slurb(30,type_tmp(m))

       IF ( n_layers_windows == 4 )  THEN
          surf%dz_win(1,m) = building_pars_slurb(31,type_tmp(m))
          surf%dz_win(2,m) = building_pars_slurb(32,type_tmp(m))
          surf%dz_win(3,m) = building_pars_slurb(33,type_tmp(m))
          surf%dz_win(4,m) = building_pars_slurb(34,type_tmp(m))

          surf%c_win(1,m) = building_pars_slurb(35,type_tmp(m))
          surf%c_win(2,m) = building_pars_slurb(36,type_tmp(m))
          surf%c_win(3,m) = building_pars_slurb(37,type_tmp(m))
          surf%c_win(4,m) = building_pars_slurb(38,type_tmp(m))

          surf%lambda_win(1,m) = building_pars_slurb(39,type_tmp(m))
          surf%lambda_win(2,m) = building_pars_slurb(40,type_tmp(m))
          surf%lambda_win(3,m) = building_pars_slurb(41,type_tmp(m))
          surf%lambda_win(4,m) = building_pars_slurb(42,type_tmp(m))
       ENDIF

       surf%transmissivity_win(m) = building_pars_slurb(43,type_tmp(m))

       surf%albedo_win(m) = building_pars_slurb(44,type_tmp(m))
       surf%emiss_win(m)  = building_pars_slurb(45,type_tmp(m))

    ENDDO

!
!-- Process pavement type.
    type_tmp(:) = -9999
    CALL get_grid_variable_1d_int( 'pavement_type', type_tmp, pavement_type )
    CALL check_grid_variable_1d_int( 'pavement_type', type_tmp, 1, 5 )
    DO  m = 1, surf%ns
       IF ( n_layers_roads == 4 )  THEN
          surf%dz_road(1,m) = pavement_pars_slurb(0,type_tmp(m))
          surf%dz_road(2,m) = pavement_pars_slurb(1,type_tmp(m))
          surf%dz_road(3,m) = pavement_pars_slurb(2,type_tmp(m))
          surf%dz_road(4,m) = pavement_pars_slurb(3,type_tmp(m))

          surf%c_road(1,m) = pavement_pars_slurb(4,type_tmp(m))
          surf%c_road(2,m) = pavement_pars_slurb(5,type_tmp(m))
          surf%c_road(3,m) = pavement_pars_slurb(6,type_tmp(m))
          surf%c_road(4,m) = pavement_pars_slurb(7,type_tmp(m))

          surf%lambda_road(1,m) = pavement_pars_slurb(8,type_tmp(m))
          surf%lambda_road(2,m) = pavement_pars_slurb(9,type_tmp(m))
          surf%lambda_road(3,m) = pavement_pars_slurb(10,type_tmp(m))
          surf%lambda_road(4,m) = pavement_pars_slurb(11,type_tmp(m))

          surf%z0_road(m)     = pavement_pars_slurb(12,type_tmp(m))
          surf%z0h_road(m)    = pavement_pars_slurb(12,type_tmp(m)) * 1.0E-2
          surf%albedo_road(m) = pavement_pars_slurb(13,type_tmp(m))
          surf%emiss_road(m)  = pavement_pars_slurb(14,type_tmp(m))
       ENDIF
    ENDDO
    DEALLOCATE( type_tmp )

!
!-- Process material layer information such as thickness, heat capacities, if given.
!-- By default, use information provided on building type.
    CALL get_grid_variable_1d_real( 'albedo_roof', surf%albedo_roof )
    CALL check_grid_variable_1d_real( 'albedo_roof', surf%albedo_roof, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'dz_roof', surf%dz_roof )
    CALL check_grid_variable_2d_real( 'dz_roof', surf%dz_roof, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
    CALL get_grid_variable_1d_real( 'emiss_roof', surf%emiss_roof )
    CALL check_grid_variable_1d_real( 'emiss_roof', surf%emiss_roof, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'c_roof', surf%c_roof )
    CALL check_grid_variable_2d_real( 'c_roof', surf%c_roof, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
!-- so multiply c_roof by dz_roof.
    surf%c_roof = surf%c_roof * surf%dz_roof

    CALL get_grid_variable_1d_real( 'z0_roof', surf%z0_roof )
    CALL check_grid_variable_1d_real( 'z0_roof', surf%z0_roof,                                     &
                                      TINY( 1.0_wp ), 0.5_wp * MINVAL( surf%z_mo ) )
    CALL get_grid_variable_1d_real( 'z0h_roof', surf%z0h_roof )
    CALL check_grid_variable_1d_real( 'z0h_roof', surf%z0h_roof,                                   &
                                      TINY( 1.0_wp ), 0.5_wp * MINVAL( surf%z_mo ) )
    CALL get_grid_variable_2d_real( 'lambda_roof', surf%lambda_roof )
    CALL check_grid_variable_2d_real( 'lambda_roof', surf%lambda_roof,                             &
                                      TINY( 1.0_wp ), HUGE( 1.0_wp )  )
    CALL get_grid_variable_1d_real( 'albedo_wall', surf%albedo_wall )
    CALL check_grid_variable_1d_real( 'albedo_wall', surf%albedo_wall, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'dz_wall', surf%dz_wall )
    CALL check_grid_variable_2d_real( 'dz_wall', surf%dz_wall, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
    CALL get_grid_variable_1d_real( 'emiss_wall', surf%emiss_wall )
    CALL check_grid_variable_1d_real( 'emiss_wall', surf%emiss_wall, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'c_wall', surf%c_wall )
    CALL check_grid_variable_2d_real( 'c_wall', surf%c_wall, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
!-- so multiply c_wall by dz_wall.
    surf%c_wall = surf%c_wall * surf%dz_wall
    CALL get_grid_variable_1d_real( 'z0_wall', surf%z0_wall )
    CALL check_grid_variable_1d_real( 'z0_wall', surf%z0_wall, TINY( 1.0_wp ), 1.0_wp )
    CALL get_grid_variable_2d_real( 'lambda_wall', surf%lambda_wall )
    CALL check_grid_variable_2d_real( 'lambda_wall', surf%lambda_wall,                             &
                                      TINY( 1.0_wp ), HUGE( 1.0_wp )  )
    CALL get_grid_variable_1d_real( 'albedo_window', surf%albedo_win )
    CALL check_grid_variable_1d_real( 'albedo_window', surf%albedo_win, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'dz_window', surf%dz_win )
    CALL check_grid_variable_2d_real( 'dz_window', surf%dz_win, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
    CALL get_grid_variable_1d_real( 'emiss_window', surf%emiss_win )
    CALL check_grid_variable_1d_real( 'emiss_window', surf%emiss_win, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'c_window', surf%c_win )
    CALL check_grid_variable_2d_real( 'c_window', surf%c_win, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity, so multiply c_win
!-- by dz_win.
    surf%c_win = surf%c_win * surf%dz_win
    CALL get_grid_variable_2d_real( 'lambda_window', surf%lambda_win )
    CALL check_grid_variable_2d_real( 'lambda_window', surf%lambda_win,                            &
                                      TINY( 1.0_wp ), HUGE( 1.0_wp )  )
    CALL get_grid_variable_1d_real( 'transmissivity_window', surf%transmissivity_win )
    CALL check_grid_variable_1d_real( 'transmissivity_window', surf%transmissivity_win,            &
                                      0.0_wp, 1.0_wp )

    CALL get_grid_variable_1d_real( 'window_fraction', surf%f_win, window_fraction )
    CALL check_grid_variable_1d_real( 'window_fraction', surf%f_win, 0.0_wp, 1.0_wp )


    CALL get_grid_variable_1d_real( 'albedo_road', surf%albedo_road )
    CALL check_grid_variable_1d_real( 'albedo_road', surf%albedo_road, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'dz_road', surf%dz_road )
    CALL check_grid_variable_2d_real( 'dz_road', surf%dz_road, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
    CALL get_grid_variable_1d_real( 'emiss_road', surf%emiss_road )
    CALL check_grid_variable_1d_real( 'emiss_road', surf%emiss_road, 0.0_wp, 1.0_wp )
    CALL get_grid_variable_2d_real( 'c_road', surf%c_road )
    CALL check_grid_variable_2d_real( 'c_road', surf%c_road, TINY( 1.0_wp ), HUGE( 1.0_wp ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
!-- so multiply c_road by dz_road.
    surf%c_road = surf%c_road * surf%dz_road
    CALL get_grid_variable_1d_real( 'z0_road', surf%z0_road )
    CALL check_grid_variable_1d_real( 'z0_road', surf%z0_road, TINY( 1.0_wp ), 1.0_wp )
    CALL check_grid_variable_1d_real( 'z0h_road', surf%z0h_road, TINY( 1.0_wp ), 1.0_wp )
    CALL get_grid_variable_2d_real( 'lambda_road', surf%lambda_road )
    CALL check_grid_variable_2d_real( 'lambda_road', surf%lambda_road,                             &
                                      TINY( 1.0_wp ), HUGE( 1.0_wp )  )
!
!-- Compute weighted wall-window albedo.
    DO  m = 1, surf%ns
       surf%albedo_wall_win(m) = ( 1.0_wp - surf%f_win(m) ) * surf%albedo_wall(m) +                &
                                 surf%f_win(m) * surf%albedo_win(m)
    ENDDO
!
!-- Compute the cumulative layer thickness zw for windows.
    DO  m = 1, surf%ns
       surf%zw_win(nzt_win, m) = surf%dz_win(nzt_win, m)
       DO  k = nzt_win+1, nzb_win
          surf%zw_win(k,m) = surf%zw_win(k-1,m) + surf%dz_win(k,m)
       ENDDO
    ENDDO

 END SUBROUTINE process_surface_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes a steady-state solution for 1D heat equation using Gauss-Seidel iteration.
!> This is used to initialize the material temperatures for roofs, walls, windows and roads,
!> shortening the time required for the spinup. For windows, SW absorption is not considered.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION calc_1d_heat_equation( result_size , t_bc_1, t_bc_2, lambda ) RESULT( t_result )

    INTEGER(iwp), INTENT(IN) ::  result_size  !< output target size

    REAL(wp), INTENT(IN) ::  t_bc_1  !< outer t boundary condition
    REAL(wp), INTENT(IN) ::  t_bc_2  !< inner t boundary condition

    REAL(wp), DIMENSION(:), INTENT(IN) ::  lambda  !< total layer heat conductivity

    INTEGER(iwp) ::  ix  !< iteration counter
    INTEGER(iwp) ::  kx  !< layer running index

    REAL(wp), PARAMETER ::  omega = 1.0_wp   !< relaxation to control convergence
    REAL(wp), PARAMETER ::  tol = 1.0E-6_wp  !< maximum residual for convergence

    REAL(wp) ::  res    !< iteration residual for convergence check
    REAL(wp) ::  t_old  !< previous t of layer for convergence check

    REAL(wp), DIMENSION(result_size) ::  t_result  !< result t profile
    REAL(wp), DIMENSION(1:result_size+1) ::  t     !< intermediate t array containing also the BCs


!
!-- Set boundary conditions for the iteration array. The layer against the atmosphere will have a
!-- constant boundary condition and the other boundary is treated similarly to the inner layer
!-- boundary condition in prognostic equations. Thus, an extra layer is neeeded for the inner
!-- temperature array for iteration.
    t(LBOUND( t, 1 )) = t_bc_1
    t(UBOUND( t, 1 )) = t_bc_2

!
!-- Set initial guess for temperature for subsurface layers.
    t(LBOUND( t, 1 )+1:UBOUND( t, 1 )-1) = (t_bc_1 + t_bc_2) / 2.0_wp

!
!-- Gauss-Seidel iteration.
    DO  ix = 1, 1000
       DO  kx = LBOUND( t, 1 )+1, UBOUND( t, 1 )-1
          t_old = t(kx)
          res = 0.0_wp
          t(kx) = ( lambda(kx) * ( t(kx+1) - t(kx) ) + lambda(kx-1) * ( t(kx-1) - t(kx) ) ) /      &
                  ( lambda(kx) + lambda(kx-1) ) * omega + t(kx)
          res = MAX( res, ABS( t(kx) - t_old ) )
       ENDDO
!
!--    Check for convergence using the stored maximum residual.
       IF ( res < tol )  EXIT
    ENDDO
!
!-- Return the solution.
    t_result(:) = t(LBOUND( t, 1 ):UBOUND( t, 1 )-1)

 END FUNCTION calc_1d_heat_equation

 END SUBROUTINE slurb_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocates SLUrb model arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_init_arrays

    USE indices,                                                                                   &
       ONLY:  nx,                                                                                  &
              ny

    USE netcdf_data_input_mod,                                                                     &
       ONLY:  check_existence,                                                                     &
              close_input_file,                                                                    &
              get_dimension_length,                                                                &
              get_attribute,                                                                       &
              get_variable,                                                                        &
              inquire_num_variables,                                                               &
              inquire_variable_names,                                                              &
              open_read_file

    USE control_parameters,                                                                        &
       ONLY:  coupling_char


    CHARACTER(LEN=100) ::  input_file_slurb = 'PIDS_SLURB'  !< name of driver file which comprises SLUrb input data

    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< array of variable names in the input driver

    INTEGER(iwp) ::  i         !< loop index for x in 2D loops
    INTEGER(iwp) ::  id_slurb  !< netCDF id of the input netCDF file
    INTEGER(iwp) ::  j         !< loop index for y in 2D loops
    INTEGER(iwp) ::  num_vars  !< number of variables in the input netCDF file
    INTEGER(iwp) ::  nx_f      !< number of grid points in x direction in the netCDF file
    INTEGER(iwp) ::  ny_f      !< number of grid points in y direction in the netCDF file

    LOGICAL ::  input_file_present  !< flag to indicate that a netCDF input has been found


    IF ( debug_output )  CALL debug_message( 'slurb_init_arrays', 'start' )

!
!-- Allocate urban fraction based on the PALM domain grid.
    IF ( .NOT. ALLOCATED( fr_urb ) )  ALLOCATE( fr_urb(nys:nyn,nxl:nxr) )

!
!-- We need to read the possible netCDF input in order to find the correct grid mapping. Otherwise
!-- the netCDF input is read in slurb_init, not here.
    INQUIRE( FILE = TRIM( input_file_slurb ) //  TRIM( coupling_char ), EXIST = input_file_present )

    IF ( input_file_present )  THEN
!
!--    Open the input file
       CALL open_read_file( TRIM( input_file_slurb ) // TRIM( coupling_char ), id_slurb )

       CALL get_dimension_length( id_slurb, nx_f, 'x'  )
       CALL get_dimension_length( id_slurb, ny_f, 'y'  )

       IF ( nx_f-1 /= nx  .OR.  ny_f-1 /= ny )  THEN
          message_string = 'Grid dimensions of the SLUrb input driver does not match ' //          &
                           'the PALM grid configuration in initialization_parameters.'
          CALL message( 'slurb_init_arrays', 'SLU00312', 2, 2, myid, 6, 0 )
       ENDIF

!
!--    Check the extistence of urban fraction.
       CALL inquire_num_variables( id_slurb, num_vars )
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_slurb, var_names )

!
!--    Check that urban_fraction exists somewhere.
       IF ( check_existence( var_names, 'urban_fraction' ) )  THEN
          CALL get_variable( id_slurb, 'urban_fraction', fr_urb, nxl, nxr, nys, nyn, nbgp=0 )
          DO  i = nxl, nxr
             DO  j = nys, nyn
!
!--             Check that urban fraction is on valid range.
                IF ( ( fr_urb(j,i) < 0.0_wp )  .OR.   ( fr_urb(j,i) > 1.0_wp ) )  THEN
                   WRITE( message_string, * ) 'urban_fraction = ', fr_urb(j,i),                    &
                                              ' is out of the valid range at (j,i) = ', j, i,      &
                                              ' in file ' // TRIM( input_file_slurb ) // '.'
                   CALL message( 'slurb_init_arrays', 'SLU0032', 2, 2, myid, 6, 0 )
                ENDIF
             ENDDO
          ENDDO

       ELSE
          IF ( urban_fraction == -9999.0_wp)  THEN
             WRITE( message_string, * ) 'urban_fraction not set in slurb_parameters or ' //        &
                                        TRIM( input_file_slurb ) // '.'
             CALL message( 'slurb_init_arrays', 'SLU0033', 2, 2, myid, 6, 0 )
          ENDIF
          DO  i = nxl, nxr
             DO  j = nys, nyn
                fr_urb(j,i) = urban_fraction
             ENDDO
          ENDDO
       ENDIF
!
!--    Close input file.
       CALL close_input_file( id_slurb )
    ELSE
!
!--    Homogeneous namelist-based initialization.
       IF ( urban_fraction == -9999.0_wp)  THEN
          WRITE( message_string, * ) 'Urban fraction not set in slurb_parameters or '              &
                                     // TRIM( input_file_slurb ) // '.'
          CALL message( 'slurb_init_arrays', 'SLU0034', 2, 2, myid, 6, 0 )
       ENDIF
       DO  i = nxl, nxr
          DO  j = nys, nyn
             fr_urb(j,i) = urban_fraction
          ENDDO
       ENDDO
    ENDIF

!
!-- Count the number of tiles for the 1D SLUrb grid mapping.
    surf%ns = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( fr_urb(j,i) > urb_thres ) surf%ns = surf%ns + 1
       ENDDO
    ENDDO

!
!-- Init SLUrb 1D grid mapping
    ALLOCATE( surf%i(1:surf%ns) )
    ALLOCATE( surf%j(1:surf%ns) )
    ALLOCATE( surf%m(nys:nyn,nxl:nxr) )
    ALLOCATE( surf%dt_max(1:surf%ns) )
    CALL init_grid

!
!-- Initialize bounds for subsurface layers.
    nzt_road = 1
    nzb_road = n_layers_roads
    nzt_roof = 1
    nzb_roof = n_layers_roofs
    nzt_wall = 1
    nzb_wall = n_layers_walls
    nzt_win  = 1
    nzb_win  = n_layers_windows

!
!-- Bulk allocation
    ALLOCATE( surf%dz_roof(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( surf%dz_wall(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( surf%dz_road(nzt_road:nzb_road,1:surf%ns) )
    ALLOCATE( surf%dz_win(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%zw_win(nzt_win:nzb_win,1:surf%ns) )

    ALLOCATE( surf%t_c_urb(1:surf%ns) )
    ALLOCATE( surf%t_rad_urb(1:surf%ns) )
    ALLOCATE( surf%t_h_urb(1:surf%ns) )
    ALLOCATE( surf%t_2m_urb(1:surf%ns) )
    ALLOCATE( surf%shf_urb(1:surf%ns) )
    ALLOCATE( surf%qsws_urb(1:surf%ns) )
    ALLOCATE( surf%ol_urb(1:surf%ns) )
    ALLOCATE( surf%rib_urb(1:surf%ns) )
    ALLOCATE( surf%ram_urb(1:surf%ns) )
    ALLOCATE( surf%usws_urb(1:surf%ns) )
    ALLOCATE( surf%vsws_urb(1:surf%ns) )

    ALLOCATE( surf%albedo_urb(1:surf%ns) )
    ALLOCATE( surf%emiss_urb(1:surf%ns) )

    ALLOCATE( surf%t_indoor(1:surf%ns) )
    ALLOCATE( surf%t_soil(1:surf%ns) )

    ALLOCATE( surf%tt_can(1:surf%ns) )
    ALLOCATE( surf%tt_wall_a(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( surf%tt_wall_b(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( surf%tt_win_a(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%tt_win_b(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%tt_roof(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( surf%tt_road(nzt_road:nzb_road,1:surf%ns) )

    ALLOCATE( surf%pt_wall_a(1:surf%ns) )
    ALLOCATE( surf%pt_wall_b(1:surf%ns) )
    ALLOCATE( surf%pt_win_a(1:surf%ns) )
    ALLOCATE( surf%pt_win_b(1:surf%ns) )
    ALLOCATE( surf%pt_roof(1:surf%ns) )
    ALLOCATE( surf%pt_road(1:surf%ns) )

    ALLOCATE( surf%shf_can(1:surf%ns) )
    ALLOCATE( surf%shf_roof(1:surf%ns) )
    ALLOCATE( surf%shf_road(1:surf%ns) )
    ALLOCATE( surf%shf_wall_a(1:surf%ns) )
    ALLOCATE( surf%shf_wall_b(1:surf%ns) )
    ALLOCATE( surf%shf_win_a(1:surf%ns) )
    ALLOCATE( surf%shf_win_b(1:surf%ns) )

    ALLOCATE( surf%shf_external(1:surf%ns) )
    ALLOCATE( surf%shf_traffic(1:surf%ns) )

    ALLOCATE( surf%ghf_road(1:surf%ns) )
    ALLOCATE( surf%ghf_roof(1:surf%ns) )
    ALLOCATE( surf%ghf_wall_a(1:surf%ns) )
    ALLOCATE( surf%ghf_wall_b(1:surf%ns) )
    ALLOCATE( surf%ghf_win_a(1:surf%ns) )
    ALLOCATE( surf%ghf_win_b(1:surf%ns) )

    ALLOCATE( surf%rad_lw_in_urb(1:surf%ns) )
    ALLOCATE( surf%rad_sw_in_urb(1:surf%ns) )
    ALLOCATE( surf%rad_lw_out_urb(1:surf%ns) )
    ALLOCATE( surf%rad_sw_out_urb(1:surf%ns) )

    ALLOCATE( surf%rad_lw_net_urb(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_urb(1:surf%ns) )

    ALLOCATE( surf%rad_lw_net_can(1:surf%ns) )

    ALLOCATE( surf%rad_lw_net_roof(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_roof(1:surf%ns) )
    ALLOCATE( surf%rad_lw_net_road(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_road(1:surf%ns) )
    ALLOCATE( surf%rad_sw_in_road(1:surf%ns) )
    ALLOCATE( surf%rad_lw_net_wall_a(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_wall_a(1:surf%ns) )
    ALLOCATE( surf%rad_lw_net_wall_b(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_wall_b(1:surf%ns) )
    ALLOCATE( surf%rad_lw_net_win_a(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_win_a(1:surf%ns) )
    ALLOCATE( surf%rad_sw_in_win_a(1:surf%ns) )
    ALLOCATE( surf%rad_lw_net_win_b(1:surf%ns) )
    ALLOCATE( surf%rad_sw_net_win_b(1:surf%ns) )
    ALLOCATE( surf%rad_sw_in_win_b(1:surf%ns) )

    ALLOCATE( surf%pt_can(1:surf%ns) )
    ALLOCATE( surf%uv_abs_can(1:surf%ns) )
    ALLOCATE( surf%uv_eff_can(1:surf%ns) )
    ALLOCATE( surf%us_can(1:surf%ns) )
    ALLOCATE( surf%rib_can(1:surf%ns) )
    ALLOCATE( surf%ol_can(1:surf%ns) )

    ALLOCATE( surf%rib_roof(1:surf%ns) )
    ALLOCATE( surf%ol_roof(1:surf%ns) )
    ALLOCATE( surf%rib_road(1:surf%ns) )
    ALLOCATE( surf%ol_road(1:surf%ns) )

    ALLOCATE( surf%us_roof(1:surf%ns) )
    ALLOCATE( surf%us_road(1:surf%ns) )

    ALLOCATE( surf%hw_can(1:surf%ns) )
    ALLOCATE( surf%anisotropic_canyon(1:surf%ns) )
    ALLOCATE( surf%theta_can(1:surf%ns) )
    ALLOCATE( surf%h_bld(1:surf%ns) )
    ALLOCATE( surf%f_bld(1:surf%ns) )
    ALLOCATE( surf%f_bld_frn(1:surf%ns) )
    ALLOCATE( surf%f_win(1:surf%ns) )
    ALLOCATE( surf%svf_road(1:surf%ns) )
    ALLOCATE( surf%svf_wall(1:surf%ns) )
    ALLOCATE( surf%z0_urb(1:surf%ns) )

    ALLOCATE( surf%rah_roof(1:surf%ns) )
    ALLOCATE( surf%rah_road(1:surf%ns) )
    ALLOCATE( surf%rah_can(1:surf%ns) )

    IF ( facade_rah_doe )  THEN
       ALLOCATE( surf%rah_wall_a(1:surf%ns) )
       ALLOCATE( surf%rah_wall_b(1:surf%ns) )
       ALLOCATE( surf%rah_win_a(1:surf%ns) )
       ALLOCATE( surf%rah_win_b(1:surf%ns) )
    ELSE
       ALLOCATE( surf%rah_facade(1:surf%ns) )
    ENDIF

    ALLOCATE( surf%lambda_roof(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( surf%c_roof(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( surf%albedo_roof(1:surf%ns) )
    ALLOCATE( surf%emiss_roof(1:surf%ns) )
    ALLOCATE( surf%z0_roof(1:surf%ns) )
    ALLOCATE( surf%z0h_roof(1:surf%ns) )
    ALLOCATE( surf%lambda_wall(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( surf%c_wall(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( surf%albedo_wall(1:surf%ns) )
    ALLOCATE( surf%emiss_wall(1:surf%ns) )
    ALLOCATE( surf%z0_wall(1:surf%ns) )
    ALLOCATE( surf%lambda_win(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%c_win(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%albedo_wall_win(1:surf%ns) )
    ALLOCATE( surf%albedo_win(1:surf%ns) )
    ALLOCATE( surf%emiss_win(1:surf%ns) )
    ALLOCATE( surf%transmissivity_win(1:surf%ns) )
    ALLOCATE( surf%absorption_win(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%lambda_road(nzt_road:nzb_road,1:surf%ns) )
    ALLOCATE( surf%c_road(nzt_road:nzb_road,1:surf%ns) )
    ALLOCATE( surf%albedo_road(1:surf%ns) )
    ALLOCATE( surf%emiss_road(1:surf%ns) )
    ALLOCATE( surf%z0_road(1:surf%ns) )
    ALLOCATE( surf%z0h_road(1:surf%ns) )

    ALLOCATE( surf%conductivity_roof(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( surf%conductivity_wall(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( surf%conductivity_win(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( surf%conductivity_road(nzt_road:nzb_road,1:surf%ns) )

    ALLOCATE( surf%z_mo(1:surf%ns) )
    ALLOCATE( surf%z_mo_can(1:surf%ns) )
    ALLOCATE( surf%uv_abs_can_coef(1:surf%ns) )
    ALLOCATE( surf%wall_hor_a_ratio(1:surf%ns) )


    ALLOCATE( surf%lw_roof_coef(1:2,1:surf%ns) )
    ALLOCATE( surf%lw_road_coef(1:4,1:surf%ns) )
    ALLOCATE( surf%lw_wall_coef(1:6,1:surf%ns) )
    ALLOCATE( surf%lw_win_coef(1:6,1:surf%ns) )
    ALLOCATE( surf%sw_ref_denom(1:surf%ns) )

    ALLOCATE( surf%us_urb(1:surf%ns) )
    ALLOCATE( surf%uv_eff1(1:surf%ns) )
    ALLOCATE( surf%uv_abs1(1:surf%ns) )
    ALLOCATE( surf%pt1(1:surf%ns) )

    ALLOCATE( t_can_1(1:surf%ns) )
    ALLOCATE( t_can_2(1:surf%ns) )
    ALLOCATE( t_wall_a_1(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( t_wall_a_2(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( t_wall_b_1(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( t_wall_b_2(nzt_wall:nzb_wall,1:surf%ns) )
    ALLOCATE( t_win_a_1(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( t_win_a_2(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( t_win_b_1(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( t_win_b_2(nzt_win:nzb_win,1:surf%ns) )
    ALLOCATE( t_roof_1(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( t_roof_2(nzt_roof:nzb_roof,1:surf%ns) )
    ALLOCATE( t_road_1(nzt_road:nzb_road,1:surf%ns) )
    ALLOCATE( t_road_2(nzt_road:nzb_road,1:surf%ns) )

    IF ( moist_physics )  THEN
       ALLOCATE( surf%tq_can(1:surf%ns))
       ALLOCATE( surf%tm_liq_roof(1:surf%ns) )
       ALLOCATE( surf%tm_liq_road(1:surf%ns) )

       ALLOCATE( surf%vpt_roof(1:surf%ns) )
       ALLOCATE( surf%vpt_road(1:surf%ns) )

       ALLOCATE( surf%q_roof(1:surf%ns) )
       ALLOCATE( surf%q_road(1:surf%ns) )
       ALLOCATE( surf%qs_roof(1:surf%ns) )
       ALLOCATE( surf%qs_road(1:surf%ns) )

       ALLOCATE( surf%qsws_can(1:surf%ns) )
       ALLOCATE( surf%qsws_roof(1:surf%ns) )
       ALLOCATE( surf%qsws_road(1:surf%ns) )
       ALLOCATE( surf%qsws_liq_roof(1:surf%ns) )
       ALLOCATE( surf%qsws_liq_road(1:surf%ns) )

       ALLOCATE( surf%c_liq_roof(1:surf%ns) )
       ALLOCATE( surf%c_liq_road(1:surf%ns) )

       ALLOCATE( surf%vpt_can(1:surf%ns) )

       ALLOCATE( surf%q1(1:surf%ns) )
       ALLOCATE( surf%vpt1(1:surf%ns) )

       ALLOCATE( surf%qsws_external(1:surf%ns) )

       ALLOCATE( q_can_1(1:surf%ns) )
       ALLOCATE( q_can_2(1:surf%ns) )
       ALLOCATE( m_liq_roof_1(1:surf%ns) )
       ALLOCATE( m_liq_roof_2(1:surf%ns) )
       ALLOCATE( m_liq_road_1(1:surf%ns) )
       ALLOCATE( m_liq_road_2(1:surf%ns) )

    ENDIF

!
!-- Set the initial timelevel.
    CALL slurb_swap_timelevel(0)

    IF ( debug_output )  CALL debug_message( 'slurb_init_arrays', 'end' )

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes the variables of the internal SLUrb grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_grid

    INTEGER(iwp) ::  m  !< running SLUrb tile index


!
!-- Fill the SLUrb grid meta arrays.
    m = 0
!
!-- Initialize m(:,:) such that zero indicates non-SLUrb tile. This is helpful when checking
!-- if the necessary inputs have been prorvided.
    surf%m(:,:) = -9999
    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( fr_urb(j,i) > 0.01_wp )  THEN
             m = m + 1
             surf%m(j,i) = m
             surf%i(m) = i
             surf%j(m) = j
             surf%dt_max(m) = HUGE( 1.0_wp )
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE init_grid

 END SUBROUTINE slurb_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Run the SLUrb model for the time step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_model

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_model'
       CALL debug_message( debug_string, 'start' )
    ENDIF

    CALL slurb_update_external_vars

    CALL slurb_radiation_model

    CALL slurb_resistance_model

    CALL slurb_energy_balance_model

    CALL slurb_canyon_model

    CALL slurb_urban_aggregation_model

    CALL slurb_atmospheric_model_coupler

    IF ( first_call )  first_call = .FALSE.

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_model'
       CALL debug_message( debug_string, 'end' )
    ENDIF

 END SUBROUTINE slurb_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &slurb_parameters for SLUrb model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter

    INTEGER(iwp) ::  io_status !< status after reading the nameslist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module


    NAMELIST /slurb_parameters/ aero_roughness_heat,                                               &
                                anisotropic_street_canyons,                                        &
                                building_frontal_area_fraction,                                    &
                                building_height,                                                   &
                                building_indoor_temperature,                                       &
                                building_plan_area_fraction,                                       &
                                building_type,                                                     &
                                deep_soil_temperature,                                             &
                                facade_resistance_parametrization,                                 &
                                moist_physics,                                                     &
                                n_layers_roads,                                                    &
                                n_layers_roofs,                                                    &
                                n_layers_walls,                                                    &
                                n_layers_windows,                                                  &
                                pavement_type,                                                     &
                                qsws_external,                                                     &
                                shf_external,                                                      &
                                shf_traffic,                                                       &
                                street_canyon_aspect_ratio,                                        &
                                street_canyon_orientation,                                         &
                                street_canyon_wspeed_factor,                                       &
                                switch_off_module,                                                 &
                                urban_fraction,                                                    &
                                urban_roughness_length,                                            &
                                window_fraction

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, slurb_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    slurb_parameters namelist was found and read correctly. Set flag that indicates that
!--    the land surface model is switched on.
       IF ( .NOT. switch_off_module )  slurb = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    slurb_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'slurb_parameters', line )

    ENDIF

  END SUBROUTINE slurb_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Swap timelevel of the SLUrb model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_swap_timelevel ( mod_count )

    INTEGER, INTENT(IN) ::  mod_count

    SELECT CASE ( mod_count )

       CASE ( 0 )
          surf%q_can => q_can_1; surf%q_can_p => q_can_2
          surf%t_can => t_can_1; surf%t_can_p => t_can_2
          surf%m_liq_road => m_liq_road_1; surf%m_liq_road_p => m_liq_road_2
          surf%m_liq_roof => m_liq_roof_1; surf%m_liq_roof_p => m_liq_roof_2
          surf%t_wall_a => t_wall_a_1; surf%t_wall_a_p => t_wall_a_2
          surf%t_wall_b => t_wall_b_1; surf%t_wall_b_p => t_wall_b_2
          surf%t_win_a => t_win_a_1; surf%t_win_a_p => t_win_a_2
          surf%t_win_b => t_win_b_1; surf%t_win_b_p => t_win_b_2
          surf%t_road => t_road_1; surf%t_road_p => t_road_2
          surf%t_roof => t_roof_1; surf%t_roof_p => t_roof_2

       CASE ( 1 )
          surf%q_can => q_can_2; surf%q_can_p => q_can_1
          surf%t_can => t_can_2; surf%t_can_p => t_can_1
          surf%m_liq_road => m_liq_road_2; surf%m_liq_road_p => m_liq_road_1
          surf%m_liq_roof => m_liq_roof_2; surf%m_liq_roof_p => m_liq_roof_1
          surf%t_wall_a => t_wall_a_2; surf%t_wall_a_p => t_wall_a_1
          surf%t_wall_b => t_wall_b_2; surf%t_wall_b_p => t_wall_b_1
          surf%t_win_a => t_win_a_2; surf%t_win_a_p => t_win_a_1
          surf%t_win_b => t_win_b_2; surf%t_win_b_p => t_win_b_1
          surf%t_roof => t_roof_2; surf%t_roof_p => t_roof_1
          surf%t_road => t_road_2; surf%t_road_p => t_road_1

    END SELECT

 END SUBROUTINE slurb_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Shortwave and longwave radiation parametrisations of the model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_radiation_model

    INTEGER(iwp) ::  day_of_year  !< day of year for the current day
    INTEGER(iwp) ::  i            !< loop index
    INTEGER(iwp) ::  j            !< loop index
    INTEGER(iwp) ::  k_topo       !< k index of topography top
    INTEGER(iwp) ::  k_atm        !< k index of the first atmospheric level
    INTEGER(iwp) ::  m            !< running index of surface tiles

    REAL(wp) ::  azimuth        !< solar azimuth angle
    REAL(wp) ::  second_of_day  !< second of the current day
    REAL(wp) ::  tan_zenith     !< tangent of the solar zenith angle
    REAL(wp) ::  zenith         !< solar zenith angle


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_radiation_model'
       CALL debug_message( debug_string, 'start' )
    ENDIF

!
!-- Calculate solar angles if not already done by RTM.
    IF ( .NOT. radiation_interactions )  THEN
       CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                  &
                           second_of_day = second_of_day )
       CALL calc_zenith( day_of_year, second_of_day )
    ENDIF

    azimuth = ATAN2( sun_dir_lon, sun_dir_lat )
    zenith = ACOS( cos_zenith )

!
!-- Split the incoming SW radiation into direct and diffuse parts.
!-- Direct-diffuse SW split is quite weirdly done in the radiation mod if radiation
!-- interactions are enabled. However, we do need it here even without interactions.
    IF ( cos_zenith > 0.0_wp )  CALL radiation_calc_diffusion_radiation

    DO  m = 1, surf%ns

       i = surf%i(m)
       j = surf%j(m)
       k_topo = topo_top_ind(j,i,0)
       k_atm = topo_top_ind(j,i,0) + 1

!
!--    Update SLUrb internal radiative fluxes based on the new surface temperatures
!--    Compute the internal longwave radiation interactions at every timestep.
       CALL calc_rad_lw

!
!--    Compute the SW radiation fluxesd.
!--    Do this only if the radiation model has updated SW fluxes at previous timestep,
!--    as otherwise the computation would just yield the same fluxes.
       IF ( radiation_called  .OR.  first_call )  CALL calc_rad_sw
    ENDDO

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_radiation_model'
       CALL debug_message( debug_string, 'end' )
    ENDIF

!
!-- Private functions and subroutines of slurb_radiation_model.
    CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the LW radiative fluxes and their differentials for the time step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rad_lw

    REAL(wp) ::  t_rad_sky  !< Radiative temperature of the sky


!
!-- Compute the effective radiative temperature of the incoming LW radiation.
    surf%rad_lw_in_urb(m) = rad_lw_in(0,j,i)
    t_rad_sky = SQRT( SQRT( surf%rad_lw_in_urb(m) / sigma_sb ) )

!
!-- Computation of net LW fluxes based on Lemonsu et al. 2012 Eqs. (1-3) (+ windows).
!-- Note that these are NOT YET the final net longwave fluxes for the surfaces, as the term
!-- dependent on the surface's own surface temperature (coef=1) is omitted at this stage.
!-- This term is added after computing the prognostic equation for the surface temperature,
!-- as it is included in the prognostic equations in an linearized form.
    surf%rad_lw_net_roof(m) = surf%lw_roof_coef(2,m) * surf%rad_lw_in_urb(m)


    surf%rad_lw_net_road(m) = surf%lw_road_coef(2,m) * surf%rad_lw_in_urb(m) +                     &
                              surf%lw_road_coef(3,m) * surf%t_wall_a(nzt_wall,m)**4 +              &
                              surf%lw_road_coef(3,m) * surf%t_wall_b(nzt_wall,m)**4 +              &
                              surf%lw_road_coef(4,m) * surf%t_win_a(nzt_win,m)**4 +                &
                              surf%lw_road_coef(4,m) * surf%t_win_b(nzt_win,m)**4

!
!-- The term dependent on t_wall_b is omitted at this stage, as for isotropic canyons the mean wall
!-- temperature is used, including both wall A and B interactions. Thus, the terms for both
!-- t_wall_a and t_wall_b have to be included in linearization. For anisotropic canyons there is
!-- no direct dependence, so it can be directly added (see below).
    surf%rad_lw_net_wall_a(m) = surf%lw_wall_coef(2,m) * surf%rad_lw_in_urb(m) +                   &
                                surf%lw_wall_coef(4,m) * surf%t_win_a(nzt_win,m)**4 +              &
                                surf%lw_wall_coef(5,m) * surf%t_win_b(nzt_win,m)**4 +              &
                                surf%lw_wall_coef(6,m) * surf%t_road(nzt_road,m)**4

    IF ( surf%f_win(m) > 0.0_wp )  THEN
       surf%rad_lw_net_win_a(m) = surf%lw_win_coef(2,m) * surf%rad_lw_in_urb(m) +                  &
                                  surf%lw_win_coef(4,m) * surf%t_wall_a(nzt_wall,m)**4 +           &
                                  surf%lw_win_coef(5,m) * surf%t_wall_b(nzt_wall,m)**4 +           &
                                  surf%lw_win_coef(6,m) * surf%t_road(nzt_road,m)**4
    ENDIF

!
!-- Inverse for facade B, if anisotropic canyons are used. If not, copy.
    IF ( surf%anisotropic_canyon(m) )  THEN
!
!--    In case of anisotropic canyons, t_wall_b doesn't have dependency on t_wall_a in the
!--    prognostic equation, and thus it's contribution to longwave balance can be directly added
!--    to the net longwave radiation before prognostic equations. Vice versa for t_wall_b.
       surf%rad_lw_net_wall_a(m) = surf%rad_lw_net_wall_a(m) +                                     &
                                   surf%lw_wall_coef(3,m) * surf%t_wall_b(nzt_wall,m)**4

!
!--    Note that for wall (and window) B the coefficients 4 and 5 are also swapped.
       surf%rad_lw_net_wall_b(m) = surf%lw_wall_coef(2,m) * surf%rad_lw_in_urb(m) +                &
                                   surf%lw_wall_coef(3,m) * surf%t_wall_a(nzt_wall,m)**4 +         &
                                   surf%lw_wall_coef(4,m) * surf%t_win_b(nzt_win,m)**4 +           &
                                   surf%lw_wall_coef(5,m) * surf%t_win_a(nzt_win,m)**4 +           &
                                   surf%lw_wall_coef(6,m) * surf%t_road(nzt_road,m)**4

       IF ( surf%f_win(m) > 0.0_wp )  THEN
          surf%rad_lw_net_win_a(m) = surf%rad_lw_net_win_a(m) *                                    &
                                     surf%lw_win_coef(3,m) * surf%t_win_b(nzt_win,m)**4

          surf%rad_lw_net_win_b(m) = surf%lw_win_coef(2,m) * surf%rad_lw_in_urb(m) +               &
                                     surf%lw_win_coef(3,m) * surf%t_win_a(nzt_win,m)**4 +          &
                                     surf%lw_win_coef(4,m) * surf%t_wall_b(nzt_wall,m)**4 +        &
                                     surf%lw_win_coef(5,m) * surf%t_wall_a(nzt_wall,m)**4 +        &
                                     surf%lw_win_coef(6,m) * surf%t_road(nzt_road,m)**4
       ENDIF
    ELSE
       surf%rad_lw_net_wall_b(m) = surf%rad_lw_net_wall_a(m)
       surf%rad_lw_net_win_b(m)  = surf%rad_lw_net_win_a(m)
    ENDIF

 END SUBROUTINE calc_rad_lw


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the SW radiative fluxes for the time step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rad_sw

    REAL(wp) ::  rad_sw_diff_road      !< incoming diffuse shortwave radiation on road
    REAL(wp) ::  rad_sw_diff_wall_a    !< incoming diffuse shortwave radiation on wall A
    REAL(wp) ::  rad_sw_diff_wall_b    !< incoming diffuse shortwave radiation on wall B
    REAL(wp) ::  rad_sw_dir_road       !< incoming direct shortwave radiation on road
    REAL(wp) ::  rad_sw_dir_wall_a     !< incoming direct shortwave radiation on wall A
    REAL(wp) ::  rad_sw_dir_wall_b     !< incoming direct shortwave radiation on wall B
    REAL(wp) ::  rad_sw_ref_nomin      !< nominator of the sum of reflections at infinity.
    REAL(wp) ::  rad_sw_wall_modifier  !< modifier term for anisotropic walls
    REAL(wp) ::  theta0                !< critical canyon orientation for road illumination
    REAL(wp) ::  w_inf                 !< mean wall reflection at infinity



!
!-- Check if there is any shortwave radiation to take care of in the first place.
    IF ( .NOT. ( cos_zenith > 0.0_wp ) )  THEN
       surf%rad_sw_in_urb(m)     = 0.0_wp
       surf%rad_sw_net_urb(m)    = 0.0_wp
       surf%rad_sw_net_roof(m)   = 0.0_wp
       surf%rad_sw_net_road(m)   = 0.0_wp
       surf%rad_sw_net_wall_a(m) = 0.0_wp
       surf%rad_sw_net_wall_b(m) = 0.0_wp
       surf%albedo_urb(m)        = 0.1_wp
       RETURN
    ENDIF

    surf%rad_sw_in_urb(m) = rad_sw_in_dir(j,i) + rad_sw_in_diff(j,i)

!
!-- Compute the net shortwave radiation for roofs, which is the simplest case.
    surf%rad_sw_net_roof(m) = ( 1.0_wp - surf%albedo_roof(m) ) * surf%rad_sw_in_urb(m)

!
!-- Next, compute then et shortwave radiation within the street canyon. This is quite complex,
!-- including the effect of shading and within-canyon reflections. See Lemonsu et al. (2012)
!-- for reference.

!
!-- Calculate tangent of the zenith angle, with limiters and safety margins applied to prevent
!-- floating point overflows and division by zero. Shouldn't affect the physics too much.
    IF ( ABS( 0.5_wp * pi - zenith ) < 1.0E-6_wp )  THEN
       IF ( 0.5_wp * pi - zenith >  0.0_wp )  tan_zenith = TAN( 0.5_wp * pi - 1.0E-6_wp )
       IF ( 0.5_wp * pi - zenith <= 0.0_wp )  tan_zenith = TAN( 0.5_wp * pi + 1.0E-6_wp )
    ELSEIF ( ABS( zenith ) < 1.0E-6_wp )  THEN
       tan_zenith = SIGN(1.0_wp, zenith) * TAN( 1.0E-6_wp )
    ELSE
       tan_zenith = TAN( zenith )
    ENDIF

!
!-- Direct SW radiation received by the walls (and windows), the road and vegetation.
    IF ( surf%anisotropic_canyon(m) )  THEN
!
!--    Lemonsu et al. (2012) Eq. (A1)
!--    @note There is an error in this equation in the article. It should be that
!--    the direct radiation on road should decrease when difference between the sun azimuth
!--    angles increase, not vice versa.
       rad_sw_dir_road = rad_sw_in_dir(j,i) * MAX( 0.0_wp, 1.0_wp - surf%hw_can(m) *               &
                         tan_zenith *  SIN( ABS( azimuth - surf%theta_can(m) ) ) )

!
!--    Lemonsu et al. (2012) Eqs. (A2-A4)
       rad_sw_dir_wall_a = ( rad_sw_in_dir(j,i) - rad_sw_dir_road ) * 0.5_wp / surf%hw_can(m)

       IF ( SIN( azimuth - surf%theta_can(m) ) > 0.0_wp )  THEN
          rad_sw_dir_wall_a = 2.0_wp * rad_sw_dir_wall_a
          rad_sw_dir_wall_b = 0.0_wp
       ELSE
          rad_sw_dir_wall_b = 2.0_wp * rad_sw_dir_wall_a
          rad_sw_dir_wall_a = 0.0_wp
       ENDIF

    ELSE
!
!--    Revert to the anisotropic integrated solution by Masson (2000).
!
!--    Calculate the critical canyon orientation theta0 for anisotropic street canyons.
       theta0 = ASIN( MIN( 1.0_wp / ( tan_zenith * surf%hw_can(m) ), 1.0_wp ) )

!
!--    Masson (2000) Eqs. (13-15)
       rad_sw_dir_road = rad_sw_in_dir(j,i) * ( 2.0_wp * theta0 / pi -                             &
                         2.0_wp * tan_zenith / pi * surf%hw_can(m) * ( 1.0_wp - COS( theta0 ) ) )

       rad_sw_dir_wall_a = ( rad_sw_in_dir(j,i) - rad_sw_dir_road ) * 0.5_wp / surf%hw_can(m)

       rad_sw_dir_wall_b = rad_sw_dir_wall_a

   ENDIF

!
!-- Diffuse (from sky) solar radiation received by the surfaces.
    rad_sw_diff_road   = rad_sw_in_diff(j,i) * surf%svf_road(m)
    rad_sw_diff_wall_a = rad_sw_in_diff(j,i) * surf%svf_wall(m)
    rad_sw_diff_wall_b = rad_sw_diff_wall_a

!
!-- Canyon internal scattering based on both Masson (2000) Eqs. (16-20) and
!-- Lemonsu et al. (2012) Appendix A2. This has been modified to include windows: the weighted
!-- average reflection from walls and windows is taken into account by using weighted average
!-- albedo. The wall and window surfaces are assumed to be uniformly distributed.

!
!-- Nominator of the sum of reflections at infinity.
    rad_sw_ref_nomin = surf%albedo_wall_win(m) * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a +        &
                       rad_sw_dir_wall_b + rad_sw_diff_wall_b ) / 2.0_wp +                         &
                       surf%albedo_wall_win(m) * surf%svf_wall(m) * surf%albedo_road(m) *          &
                       rad_sw_dir_road

!
!-- Sum of refelctions at infinity.
    w_inf = rad_sw_ref_nomin / surf%sw_ref_denom(m)

!
!-- Total solar radiation absorbed after infinite reflections.
    surf%rad_sw_in_road(m) = rad_sw_dir_road + rad_sw_diff_road +                                  &
                             ( 1.0_wp - surf%svf_road(m) ) * w_inf
    surf%rad_sw_net_road(m) = ( 1.0_wp - surf%albedo_road(m) ) * surf%rad_sw_in_road(m)

    surf%rad_sw_net_wall_a(m) = ( 1.0_wp - surf%albedo_wall(m) ) *                                 &
                                ( 0.5_wp * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a +              &
                                             rad_sw_dir_wall_b + rad_sw_diff_wall_b )              &
                                + surf%albedo_road(m) * surf%svf_wall(m) *                         &
                                  ( rad_sw_dir_road + rad_sw_diff_road )                           &
                                + surf%albedo_road(m) * surf%svf_wall(m) *                         &
                                  ( 1.0_wp - surf%svf_road(m) ) * w_inf                            &
                                + ( 1.0_wp - 2.0_wp * surf%svf_wall(m) ) * w_inf                   &
                                )

    surf%rad_sw_net_wall_b(m) = surf%rad_sw_net_wall_a(m)

    IF ( surf%f_win(m) /= 0.0_wp  )  THEN
       surf%rad_sw_in_win_a(m) =   0.5_wp * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a               &
                                            + rad_sw_dir_wall_b + rad_sw_diff_wall_b )             &
                                 + surf%albedo_road(m) * surf%svf_wall(m) *                        &
                                   ( rad_sw_dir_road + rad_sw_diff_road )                          &
                                 + surf%albedo_road(m) * surf%svf_wall(m) *                        &
                                   ( 1.0_wp - surf%svf_road(m) ) * w_inf                           &
                                 + ( 1.0_wp - 2.0_wp * surf%svf_wall(m) ) * w_inf

       surf%rad_sw_net_win_a(m) = ( 1.0_wp - surf%albedo_win(m) ) * surf%rad_sw_in_win_a(m)

       surf%rad_sw_in_win_b(m)  = surf%rad_sw_in_win_a(m)
       surf%rad_sw_net_win_b(m) = surf%rad_sw_net_win_a(m)
    ENDIF

!
!-- Modification of reflected solar radiation for anisotropic street canyons.
    IF ( surf%anisotropic_canyon(m) )  THEN
       rad_sw_wall_modifier = ( 1.0_wp + surf%albedo_wall_win(m) *                                 &
                                ( 1.0_wp - 2.0_wp * surf%svf_wall(m) ) /                           &
                                ( 1.0_wp + surf%albedo_wall_win(m) *                               &
                                  ( 1.0_wp - 2.0_wp * surf%svf_wall(m) ) )                         &
                              ) *                                                                  &
                              0.5_wp * ( ( rad_sw_dir_wall_a + rad_sw_diff_wall_a )                &
                                       - ( rad_sw_dir_wall_b + rad_sw_diff_wall_b ) )

       surf%rad_sw_net_wall_a(m) = surf%rad_sw_net_wall_a(m) +                                     &
                                   ( 1.0_wp - surf%albedo_wall(m) ) * rad_sw_wall_modifier

       surf%rad_sw_net_wall_b(m) = surf%rad_sw_net_wall_b(m) -                                     &
                                   ( 1.0_wp - surf%albedo_wall(m) ) * rad_sw_wall_modifier

       IF ( surf%f_win(m) /= 0.0_wp )  THEN
          surf%rad_sw_in_win_a(m)  = surf%rad_sw_in_win_a(m) + rad_sw_wall_modifier
          surf%rad_sw_net_win_a(m) = surf%rad_sw_in_win_a(m) * ( 1.0_wp - surf%albedo_win(m) )
          surf%rad_sw_in_win_b(m)  = surf%rad_sw_in_win_b(m) - rad_sw_wall_modifier
          surf%rad_sw_net_win_b(m) = surf%rad_sw_in_win_b(m) * ( 1.0_wp - surf%albedo_win(m) )
       ENDIF
    ENDIF

!
!-- The upward shortwave radiation is computed as residual of absorbed radiation per uniturban
!-- area. Aggregated effective albedo of urban surface is computed so that the raditaiton models end
!-- up with the same figure for outgoing shortwave radiation.
    surf%rad_sw_out_urb(m) = surf%rad_sw_in_urb(m) -                                               &
                             ( ( 1.0_wp - surf%f_bld(m) ) *                                        &
                               ( surf%hw_can(m) * ( ( 1.0_wp - surf%f_win(m) ) *                   &
                                       ( surf%rad_sw_net_wall_a(m) + surf%rad_sw_net_wall_b(m) )   &
                                       + surf%f_win(m) *                                           &
                                       ( surf%rad_sw_net_win_a(m)  + surf%rad_sw_net_win_b(m)  )   &
                                                  )                                                &
                               + surf%rad_sw_net_road(m)                                           &
                               )                                                                   &
                             + surf%f_bld(m) * surf%rad_sw_net_roof(m)                             &
                             )

!
!-- Compute the net SW flux for diagnostics and output.
    surf%rad_sw_net_urb(m) = surf%rad_sw_in_urb(m) - surf%rad_sw_out_urb(m)

!
!-- Save effective albedo for the radiation model.
    surf%albedo_urb(m) = surf%rad_sw_out_urb(m) / surf%rad_sw_in_urb(m)

 END SUBROUTINE calc_rad_sw

 END SUBROUTINE slurb_radiation_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Model to compute aerodynamic resistances for individual surfaces and urban form.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_resistance_model


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_resistance_model'
       CALL debug_message( debug_string, 'start' )
    ENDIF
!
!-- Compute the resistances between the atmosphere and urban surface.
    CALL calc_urban_resistances
!
!-- Compute aerodynamical variables shared by within-canyon surfaces.
    CALL calc_canyon_resistances

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_resistance_model'
       CALL debug_message( debug_string, 'end' )
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the heat and momentum fluxes between the atmosphere and the urban surface.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_urban_resistances

    INTEGER(iwp) ::  m  !< running index of surface tiles

    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0_roof   !< temporary array to store logarithm
    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0h_roof  !< temporary array to store logarithm
    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0_urb    !< temporary array to store logarithm
    REAL(wp), DIMENSION(surf%ns) ::  pt_surface     !< temporary array to store weighted temperature


!
!-- Calculate logarithms of ratio z/z0.
!>  TODO: Since the ratios do not change during the simulation, they can be stored once at the
!>        and stored in surf_slurb, like it is done for the other surface types, too.
    DO  m = 1, surf%ns
       ln_z_z0_roof(m)  = LOG( surf%z_mo(m) / surf%z0_roof(m)  )
       ln_z_z0h_roof(m) = LOG( surf%z_mo(m) / surf%z0h_roof(m) )
       ln_z_z0_urb(m)   = LOG( surf%z_mo(m) / surf%z0_urb(m)   )
    ENDDO

!
!-- Compute friction velocity and aerodynamic resistance for momentum for the whole urban surface.
!-- As SLUrb doesn't explicitly compute the momentum flux for each individual surface, and as
!-- pressure drag needs to be included in the total urban drag, us_urb and rah_urb are computed
!-- using given roughness length for whole urban fabric (z0_urb, user input). To compute the MOST
!-- stability corrections, we follow the SURFEX implementation where weighted pt/vpt from canyons
!-- and roofs is used to represent the pt/vpt at roof level. For urban heat fluxes, aggregated
!-- values from roofs and canyons are directly used, so rah_urb is not needed.
    IF ( moist_physics )  THEN

       DO  m = 1, surf%ns
          pt_surface(m) = surf%f_bld(m)              * surf%vpt_roof(m) +                          &
                          ( 1.0_wp - surf%f_bld(m) ) * surf%vpt_can(m)
       ENDDO
       CALL calc_rib( surf%ns, surf%vpt1, pt_surface, surf%rib_urb, surf%uv_eff1, surf%z_mo, slurb )

    ELSE

       DO  m = 1, surf%ns
          pt_surface(m) = surf%f_bld(m)              * surf%pt_roof(m) +                           &
                          ( 1.0_wp - surf%f_bld(m) ) * surf%pt_can(m)
       ENDDO
       CALL calc_rib( surf%ns, surf%pt1, pt_surface, surf%rib_urb, surf%uv_eff1, surf%z_mo, slurb )

    ENDIF

    CALL calc_ol( surf%ns, ln_z_z0_urb, ln_z_z0_urb, surf%ol_urb, surf%rib_urb, surf%z0_urb,       &
                  surf%z0_urb, surf%z_mo )

    DO  m = 1, surf%ns
       surf%us_urb(m) = kappa * surf%uv_eff1(m) /                                                  &
                        ( LOG( surf%z_mo(m) / surf%z0_urb(m) ) -                                   &
                          psi_m( surf%z_mo(m) / surf%ol_urb(m) ) +                                 &
                          psi_m( surf%z0_urb(m) / surf%ol_urb(m) ) )
    ENDDO

!
!-- Ensure physical friction velocity (might be needed due to instabilities in e.g. initialization)
    DO  m = 1, surf%ns
       IF ( surf%us_urb(m) <= us_min ) surf%us_urb(m) = us_min

       surf%ram_urb(m) = 1.0_wp / ( kappa * surf%us_urb(m) ) *                                     &
                         ( LOG( surf%z_mo(m) / surf%z0_urb(m) ) -                                  &
                           psi_m( surf%z_mo(m) / surf%ol_urb(m) ) +                                &
                           psi_m( surf%z0_urb(m) / surf%ol_urb(m) ) )

       IF ( surf%ram_urb(m) < ram_min )  surf%ram_urb(m) = ram_min
    ENDDO
!
!-- For street canyons, effective mixing between canyon half-height and roof height is assumed, thus
!-- z_mo is used as as reference height when considering atmosphere-street canyon air mixing.
!-- This is equivalent to mixing of canyon air between the roof top level and the first atm grid
!-- level. For canyons, roughness length for the whole urban fabric (z0_urb) is used instead of
!-- local z0m/z0h. This is based on an assumption that turbulence can effectively mix the two air
!-- masses (canyon air and the atmospheric air). The same assumption is used in TEB/SURFEX.
!-- Using e.g. the Kanda et al. (2007) parametrization or any other surface parametrization for
!-- canyon z0h would yield unrealistically low mixing.
!
!-- Update z0h for roofs following Kanda et al. (2007) parametrization if enabled.
    IF ( roughness_kanda )  THEN
       DO  m = 1, surf%ns
          surf%z0h_roof(m) = surf%z0_roof(m) * 7.4_wp *                                            &
                             EXP( -1.29_wp * SQRT( SQRT( surf%z0_roof(m) * surf%us_roof(m) /       &
                                                         1.461E-5_wp) ) )
          ln_z_z0h_roof(m) = LOG( surf%z_mo(m) / surf%z0h_roof(m) )
       ENDDO
    ENDIF

    IF ( moist_physics )  THEN
       CALL calc_rib( surf%ns, surf%vpt1, surf%vpt_roof, surf%rib_roof, surf%uv_eff1, surf%z_mo,   &
                      slurb )
       CALL calc_rib( surf%ns, surf%vpt1, surf%vpt_can,  surf%rib_can,  surf%uv_eff1, surf%z_mo,   &
                      slurb )
    ELSE
       CALL calc_rib( surf%ns, surf%pt1, surf%pt_roof, surf%rib_roof, surf%uv_eff1, surf%z_mo,     &
                      slurb )
       CALL calc_rib( surf%ns, surf%pt1, surf%pt_can,  surf%rib_can,  surf%uv_eff1, surf%z_mo,     &
                      slurb )
    ENDIF

    CALL calc_ol( surf%ns, ln_z_z0_roof, ln_z_z0h_roof, surf%ol_roof, surf%rib_roof, surf%z0_roof, &
                  surf%z0h_roof, surf%z_mo )
    CALL calc_ol( surf%ns, ln_z_z0_urb, ln_z_z0_urb, surf%ol_can, surf%rib_can, surf%z0_urb,       &
                  surf%z0_urb, surf%z_mo )

!
!-- Compute the local friction velocity for roof and canyon.
    DO  m = 1, surf%ns
       surf%us_roof(m) = kappa * surf%uv_eff1(m) /                                                 &
                         ( LOG( surf%z_mo(m) / surf%z0_roof(m) ) -                                 &
                           psi_m( surf%z_mo(m) / surf%ol_roof(m) ) +                               &
                           psi_m( surf%z0_roof(m) / surf%ol_roof(m) ) )

!
!--    For canyons, use urban roughness length (assume the air mixes efficiently
!--    between the canyon air and atmosphere).
       surf%us_can(m) = kappa * surf%uv_eff1(m) /                                                  &
                        ( LOG( surf%z_mo(m) / surf%z0_urb(m) ) -                                   &
                          psi_m( surf%z_mo(m) / surf%ol_can(m) ) +                                 &
                          psi_m( surf%z0_urb(m) / surf%ol_can(m) ) )

!
!--    Ensure physical friction velocity.
       IF ( surf%us_roof(m) <= us_min )  surf%us_roof(m) = us_min
    ENDDO

!
!-- Compute the aerodynamic resistances for heat.
    DO  m = 1, surf%ns
       surf%rah_roof(m) = 1.0_wp / ( kappa * surf%us_roof(m) ) *                                   &
                          ( LOG( surf%z_mo(m) / surf%z0h_roof(m) ) -                               &
                            psi_h( surf%z_mo(m) / surf%ol_roof(m) ) +                              &
                            psi_h( surf%z0h_roof(m) / surf%ol_roof(m) ) )

       surf%rah_can(m) = 1.0_wp / ( kappa * surf%us_can(m) ) *                                     &
                         ( LOG( surf%z_mo(m) / surf%z0_urb(m) ) -                                  &
                           psi_h( surf%z_mo(m) / surf%ol_can(m) ) +                                &
                           psi_h( surf%z0_urb(m) / surf%ol_can(m) ) )

       IF ( surf%rah_roof(m) < rah_min )  surf%rah_roof(m) = rah_min
       IF ( surf%rah_roof(m) > rah_max )  surf%rah_roof(m) = rah_max
!
!--    Use ram_min for canyon air as turbulence is able to mix the air.
       IF ( surf%rah_can(m) < ram_min )  surf%rah_can(m) = ram_min
    ENDDO

 END SUBROUTINE calc_urban_resistances


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Model for the surface resistances within the street canyon.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_canyon_resistances

    INTEGER(iwp) ::  i       !< loop index x-direction
    INTEGER(iwp) ::  j       !< loop index y-direction
    INTEGER(iwp) ::  k_topo  !< k-index of topography
    INTEGER(iwp) ::  m       !< running index of surface tiles

    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0_road   !< temporary array to store logarithm
    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0h_road  !< temporary array to store logarithm
    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0_roof   !< temporary array to store logarithm
    REAL(wp), DIMENSION(surf%ns) ::  ln_z_z0h_roof  !< temporary array to store logarithm


!
!-- Calculate logarithms of ratio z/z0.
!>  TODO: Since the ratios do not change during the simulation, they can be stored once at the
!>        and stored in surf_slurb, like it is done for the other surface types, too.
    DO  m = 1, surf%ns
       ln_z_z0_road(m)  = LOG( surf%z_mo_can(m) / surf%z0_road(m)  )
       ln_z_z0h_road(m) = LOG( surf%z_mo_can(m) / surf%z0h_road(m) )
       ln_z_z0_roof(m)  = LOG( surf%z_mo(m)     / surf%z0_roof(m)  )
       ln_z_z0h_roof(m) = LOG( surf%z_mo(m)     / surf%z0h_roof(m) )
    ENDDO

!
!-- Update z0h for roads following Kanda et al. (2007) parametrization if necessary.
    IF ( roughness_kanda )  THEN
       DO  m = 1, surf%ns
          surf%z0h_road(m) = surf%z0_road(m) * 7.4_wp *                                            &
                             EXP( -1.29_wp * SQRT( SQRT( surf%z0_road(m) * surf%us_road(m) /       &
                                                         1.461E-5_wp) ) )
          ln_z_z0h_road(m) = LOG( surf%z_mo_can(m) / surf%z0h_road(m) )
       ENDDO
    ENDIF

!
!-- Compute the new Obukhov length for road.
    IF ( moist_physics )  THEN
       CALL calc_rib( surf%ns, surf%vpt_can, surf%vpt_road, surf%rib_road, surf%uv_eff_can,        &
                      surf%z_mo_can, slurb )
    ELSE
       CALL calc_rib( surf%ns, surf%pt_can, surf%pt_road, surf%rib_road, surf%uv_eff_can,          &
                      surf%z_mo_can, slurb )
    ENDIF

    CALL calc_ol( surf%ns, ln_z_z0_road, ln_z_z0h_road, surf%ol_road, surf%rib_road, surf%z0_road, &
                  surf%z0h_road, surf%z_mo_can )

!
!-- Compute the local friction velocity for roads.
    DO  m = 1, surf%ns
       surf%us_road(m) = kappa * surf%uv_eff_can(m) /                                              &
                         ( LOG( surf%z_mo_can(m) / surf%z0_road(m) ) -                             &
                           psi_m( surf%z_mo_can(m) / surf%ol_road(m) ) +                           &
                           psi_m( surf%z0_road(m) / surf%ol_road(m) ) )
    ENDDO

!
!-- The resistance between the street canyon air and facades (walls and windows).
    IF ( facade_rah_doe )  THEN

       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k_topo = topo_top_ind(j,i,0)
          surf%rah_wall_a(m) = rah_doe2( k_topo, surf%t_can(m), surf%t_wall_a(nzt_wall,m),         &
                                         surf%uv_eff_can(m), .TRUE. )
          IF ( surf%rah_wall_a(m) < rah_min )  surf%rah_wall_a(m) = rah_min
          IF ( surf%rah_wall_a(m) > rah_max )  surf%rah_wall_a(m) = rah_max
          IF ( surf%f_win(m) /= 0.0_wp )  THEN
             surf%rah_win_a(m) = rah_doe2( k_topo, surf%t_can(m), surf%t_win_a(nzt_win,m),         &
                                           surf%uv_eff_can(m), .FALSE. )
             IF ( surf%rah_win_a(m) < rah_min )  surf%rah_win_a(m) = rah_min
             IF ( surf%rah_win_a(m) > rah_max )  surf%rah_win_a(m) = rah_max
          ENDIF

          IF ( surf%anisotropic_canyon(m) )  THEN
             surf%rah_wall_b(m) = rah_doe2( k_topo, surf%t_can(m), surf%t_wall_b(nzt_wall,m),      &
                                            surf%uv_eff_can(m), .TRUE. )
             IF ( surf%rah_wall_b(m) < rah_min )  surf%rah_wall_b(m) = rah_min
             IF ( surf%rah_wall_b(m) > rah_max )  surf%rah_wall_b(m) = rah_max
             IF ( surf%f_win(m) /= 0.0_wp )  THEN
                surf%rah_win_b(m) = rah_doe2( k_topo, surf%t_can(m), surf%t_win_b(nzt_win,m),      &
                                              surf%uv_eff_can(m), .FALSE. )
                IF ( surf%rah_win_b(m) < rah_min )  surf%rah_win_b(m) = rah_min
                IF ( surf%rah_win_b(m) > rah_max )  surf%rah_win_b(m) = rah_max
             ENDIF
          ENDIF
       ENDDO

    ELSEIF ( facade_rah_kray )  THEN

       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k_topo = topo_top_ind(j,i,0)
          surf%rah_facade(m) = rah_kray( k_topo, surf%z0_wall(m), surf%uv_eff_can(m) )
          IF ( surf%rah_facade(m) < rah_min )  surf%rah_facade(m) = rah_min
          IF ( surf%rah_facade(m) > rah_max )  surf%rah_facade(m) = rah_max
       ENDDO

    ELSEIF ( facade_rah_rowley )  THEN
!
!--    Rowley et al. (1930) , Cole and Sturrock (1977)  Mills (1993).
       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k_topo = topo_top_ind(j,i,0)
          surf%rah_facade(m) = c_p * rho_air_zw(k_topo) / ( 11.8_wp + 4.2_wp * surf%uv_eff_can(m) )
          IF ( surf%rah_facade(m) < rah_min )  surf%rah_facade(m) = rah_min
          IF ( surf%rah_facade(m) > rah_max )  surf%rah_facade(m) = rah_max
       ENDDO
    ENDIF

    DO  m = 1, surf%ns
       surf%rah_road(m) = 1.0_wp / ( kappa * surf%us_can(m) ) *                                    &
                          ( ln_z_z0h_road(m) -                                                     &
                            psi_h( surf%z_mo_can(m) / surf%ol_road(m) ) +                          &
                            psi_h( surf%z0h_road(m) / surf%ol_road(m) ) )

       IF ( surf%rah_road(m) < rah_min )  surf%rah_road(m) = rah_min
       IF ( surf%rah_road(m) > rah_max )  surf%rah_road(m) = rah_max
    ENDDO

 END SUBROUTINE calc_canyon_resistances


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute aerodynamic resistance for heat for vertical surfaces following DOE-2 parametrization,
!> which takes natural convection into account. Average of leeward and windward sides.
!> Source: EnegyPlus 23.2.0 Engineering Reference p.68.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION rah_doe2( k_topo, t_air, t_surf, u_eff, rough )

    LOGICAL, INTENT(IN) ::  rough  !< flag for rough surface, true for walls, false for windows

    INTEGER(iwp), INTENT(IN) ::  k_topo  !< k-index of topography

    REAL(wp), INTENT(IN) ::  t_air   !< temperature of adjacent air
    REAL(wp), INTENT(IN) ::  t_surf  !< surface temperature
    REAL(wp), INTENT(IN) ::  u_eff   !< effective wind speed

    REAL(wp), PARAMETER ::  r_f = 1.52_wp  !< surface roughness multiplier

    REAL(wp) ::  chtcn       !< convective heat transfer coefficient for natural convection
    REAL(wp) ::  chtcs       !< convective heat transfer coefficient for smooth surface
    REAL(wp) ::  chtcs_lee   !< convective heat transfer coefficient for smooth surface (leeward)
    REAL(wp) ::  chtcs_wind  !< convective heat transfer coefficient for smooth surface (windward)
    REAL(wp) ::  rah_doe2    !< resulting resistance


    chtcn = 1.31_wp * ABS( t_air - t_surf )**0.33333_wp

    chtcs_lee  = SQRT( chtcn**2 + ( 2.86_wp * u_eff**0.617_wp )**2 )
    chtcs_wind = SQRT( chtcn**2 + ( 2.38_wp * u_eff**0.89_wp  )**2 )

    chtcs = 0.5 * ( chtcs_lee + chtcs_wind )

    IF ( rough )  THEN
       rah_doe2 = c_p * rho_air_zw(k_topo) / ( chtcn + r_f * ( chtcs - chtcn ) )
    ELSE
       rah_doe2 = c_p * rho_air_zw(k_topo) / chtcs
    ENDIF

 END FUNCTION rah_doe2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute aerodynamic resistance for heat for vertical surfaces following
!> Krayenhoff & Voogt (2007).
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION rah_kray( k_topo, z0, u_eff )

    INTEGER(iwp), INTENT(IN) ::  k_topo  !< k-index of topography

    REAL(wp), INTENT(IN) ::  u_eff  !< effective wind speed
    REAL(wp), INTENT(IN) ::  z0     !< roughness length for momentum

    REAL(wp) ::  kray_coeff  !< denominator for the parametrization
    REAL(wp) ::  rah_kray    !< resulting resistance


!
!-- Compute denominator first, ensuring it is a positive number.
    kray_coeff = MAX( z0 * 1000.0_wp * ( 11.8_wp + 4.2_wp * u_eff ) - 4.0_wp, 1.0E-3_wp )

    rah_kray = c_p * rho_air_zw(k_topo) / kray_coeff

 END FUNCTION

 END SUBROUTINE slurb_resistance_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read restart data using standard Fortran I/O.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,  &
                                 nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, found )

    USE control_parameters

    USE indices

    INTEGER(iwp) ::  k                 !<
    INTEGER(iwp) ::  nxlc              !<
    INTEGER(iwp) ::  nxlf              !<
    INTEGER(iwp) ::  nxl_on_file       !< index of left boundary on former local domain
    INTEGER(iwp) ::  nxrc              !<
    INTEGER(iwp) ::  nxrf              !<
    INTEGER(iwp) ::  nxr_on_file       !< index of right boundary on former local domain
    INTEGER(iwp) ::  nync              !<
    INTEGER(iwp) ::  nynf              !<
    INTEGER(iwp) ::  nyn_on_file       !< index of north boundary on former local domain
    INTEGER(iwp) ::  nysc              !<
    INTEGER(iwp) ::  nysf              !<
    INTEGER(iwp) ::  nys_on_file       !< index of south boundary on former local domain

    LOGICAL, INTENT(OUT) ::  found

    INTEGER(iwp) ::  i   !< loop index in the x-direction
    INTEGER(iwp) ::  j   !< loop index in the y-direction
    INTEGER(iwp) ::  m   !< running SLUrb tile index

    REAL(wp),                                                                                      &
       DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::  tmp_2d   !< temporary array for 2D vars

    REAL(wp), DIMENSION(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) ::  tmp_2d_tgt   !< temporary target array for 2D vars

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )
!
!--    Prognostic 3D arrays.
       CASE ( 'surf_slurb%t_wall_a' )
          CALL rrd_local_3d( surf%t_wall_a, nzt_wall, nzb_wall )
       CASE ( 'surf_slurb%t_wall_b' )
          CALL rrd_local_3d( surf%t_wall_b, nzt_wall, nzb_wall )
       CASE ( 'surf_slurb%t_win_a' )
          CALL rrd_local_3d( surf%t_win_a, nzt_win, nzb_win )
       CASE ( 'surf_slurb%t_win_b' )
          CALL rrd_local_3d( surf%t_win_b, nzt_win, nzb_win )
       CASE ( 'surf_slurb%t_roof' )
          CALL rrd_local_3d( surf%t_roof, nzt_roof, nzb_roof )
       CASE ( 'surf_slurb%t_road' )
          CALL rrd_local_3d( surf%t_road, nzt_road, nzb_road )
!
!--    Prognostic 2D arrays.
       CASE ( 'surf_slurb%us_urb' )
          CALL rrd_local_2d( surf%us_urb )
       CASE ( 'surf_slurb%shf_urb' )
          CALL rrd_local_2d( surf%shf_urb )
       CASE ( 'surf_slurb%t_can' )
          CALL rrd_local_2d( surf%t_can )
       CASE ( 'surf_slurb%m_liq_roof' )
          CALL rrd_local_2d( surf%m_liq_roof )
       CASE ( 'surf_slurb%m_liq_road' )
          CALL rrd_local_2d( surf%m_liq_road )
       CASE ( 'surf_slurb%q_can' )
          CALL rrd_local_2d( surf%q_can )
       CASE ( 'surf_slurb%q_road' )
          CALL rrd_local_2d( surf%q_road )
       CASE ( 'surf_slurb%q_roof' )
          CALL rrd_local_2d( surf%q_roof )
       CASE ( 'surf_slurb%ol_urb' )
          CALL rrd_local_2d( surf%ol_urb )
       CASE ( 'surf_slurb%ol_roof' )
          CALL rrd_local_2d( surf%ol_roof )
       CASE ( 'surf_slurb%ol_road' )
          CALL rrd_local_2d( surf%ol_road )
       CASE ( 'surf_slurb%ol_can' )
          CALL rrd_local_2d( surf%ol_can )
       CASE ( 'surf_slurb%us_roof' )
          CALL rrd_local_2d( surf%us_roof)
       CASE ( 'surf_slurb%us_road' )
          CALL rrd_local_2d( surf%us_road )
       CASE ( 'surf_slurb%uv_eff_can' )
          CALL rrd_local_2d( surf%uv_eff_can )
!
!-- Arrays for time averaging (2D). We need to check the allocation status before reading.
       CASE ( 'slurb_albedo_urb_av' )
          IF ( .NOT.  ALLOCATED( albedo_urb_av ) )  ALLOCATE( albedo_urb_av(1:surf%ns) )
          CALL rrd_local_2d( albedo_urb_av )
       CASE ( 'slurb_c_liq_road_av' )
          IF ( .NOT.  ALLOCATED( c_liq_road_av ) )  ALLOCATE( c_liq_road_av(1:surf%ns) )
          CALL rrd_local_2d( c_liq_road_av )
       CASE ( 'slurb_c_liq_roof_av' )
          IF ( .NOT.  ALLOCATED( c_liq_roof_av ) )  ALLOCATE( c_liq_roof_av(1:surf%ns) )
          CALL rrd_local_2d( c_liq_roof_av )
       CASE ( 'slurb_emiss_urb_av' )
          IF ( .NOT.  ALLOCATED( emiss_urb_av ) )  ALLOCATE( emiss_urb_av(1:surf%ns) )
          CALL rrd_local_2d( emiss_urb_av )
       CASE ( 'slurb_ghf_road_av' )
          IF ( .NOT.  ALLOCATED( ghf_road_av ) )  ALLOCATE( ghf_road_av(1:surf%ns) )
          CALL rrd_local_2d( ghf_road_av )
       CASE ( 'slurb_ghf_roof_av' )
          IF ( .NOT.  ALLOCATED( ghf_roof_av ) )  ALLOCATE( ghf_roof_av(1:surf%ns) )
          CALL rrd_local_2d( ghf_roof_av )
       CASE ( 'slurb_ghf_wall_a_av' )
          IF ( .NOT.  ALLOCATED( ghf_wall_a_av ) )  ALLOCATE( ghf_wall_a_av(1:surf%ns) )
          CALL rrd_local_2d( ghf_wall_a_av )
       CASE ( 'slurb_ghf_wall_b_av' )
          IF ( .NOT.  ALLOCATED( ghf_wall_b_av ) )  ALLOCATE( ghf_wall_b_av(1:surf%ns) )
          CALL rrd_local_2d( ghf_wall_b_av )
       CASE ( 'slurb_ghf_win_a_av' )
          IF ( .NOT.  ALLOCATED( ghf_win_a_av ) )  ALLOCATE( ghf_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( ghf_win_a_av )
       CASE ( 'slurb_ghf_win_b_av' )
          IF ( .NOT.  ALLOCATED( ghf_win_b_av ) )  ALLOCATE( ghf_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( ghf_win_b_av )
       CASE ( 'slurb_m_liq_road_av' )
          IF ( .NOT.  ALLOCATED( m_liq_road_av ) )  ALLOCATE( m_liq_road_av(1:surf%ns) )
          CALL rrd_local_2d( m_liq_road_av )
       CASE ( 'slurb_m_liq_roof_av' )
          IF ( .NOT.  ALLOCATED( m_liq_roof_av ) )  ALLOCATE( m_liq_roof_av(1:surf%ns) )
          CALL rrd_local_2d( m_liq_roof_av )
       CASE ( 'slurb_ol_can_av' )
          IF ( .NOT.  ALLOCATED( ol_can_av ) )  ALLOCATE( ol_can_av(1:surf%ns) )
          CALL rrd_local_2d( ol_can_av )
       CASE ( 'slurb_ol_road_av' )
          IF ( .NOT.  ALLOCATED( ol_road_av ) )  ALLOCATE( ol_road_av(1:surf%ns) )
          CALL rrd_local_2d( ol_road_av )
       CASE ( 'slurb_ol_roof_av' )
          IF ( .NOT.  ALLOCATED( ol_roof_av ) )  ALLOCATE( ol_roof_av(1:surf%ns) )
          CALL rrd_local_2d( ol_roof_av )
       CASE ( 'slurb_pt_can_av' )
          IF ( .NOT.  ALLOCATED( pt_can_av ) )  ALLOCATE( pt_can_av(1:surf%ns) )
          CALL rrd_local_2d( pt_can_av )
       CASE ( 'slurb_pt_road_av' )
          IF ( .NOT.  ALLOCATED( pt_road_av ) )  ALLOCATE( pt_road_av(1:surf%ns) )
          CALL rrd_local_2d( pt_road_av )
       CASE ( 'slurb_pt_roof_av' )
          IF ( .NOT.  ALLOCATED( pt_roof_av ) )  ALLOCATE( pt_roof_av(1:surf%ns) )
          CALL rrd_local_2d( pt_roof_av )
       CASE ( 'slurb_pt_wall_a_av' )
          IF ( .NOT.  ALLOCATED( pt_wall_a_av ) )  ALLOCATE( pt_wall_a_av(1:surf%ns) )
          CALL rrd_local_2d( pt_wall_a_av )
       CASE ( 'slurb_pt_wall_b_av' )
          IF ( .NOT.  ALLOCATED( pt_wall_b_av ) )  ALLOCATE( pt_wall_b_av(1:surf%ns) )
          CALL rrd_local_2d( pt_wall_b_av )
       CASE ( 'slurb_pt_win_a_av' )
          IF ( .NOT.  ALLOCATED( pt_win_a_av ) )  ALLOCATE( pt_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( pt_win_a_av )
       CASE ( 'slurb_pt_win_b_av' )
          IF ( .NOT.  ALLOCATED( pt_win_b_av ) )  ALLOCATE( pt_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( pt_win_b_av )
       CASE ( 'slurb_q_can_av' )
          IF ( .NOT.  ALLOCATED( q_can_av ) )  ALLOCATE( q_can_av(1:surf%ns) )
          CALL rrd_local_2d( q_can_av )
       CASE ( 'slurb_q_road_av' )
          IF ( .NOT.  ALLOCATED( q_road_av ) )  ALLOCATE( q_road_av(1:surf%ns) )
          CALL rrd_local_2d( q_road_av )
       CASE ( 'slurb_q_roof_av' )
          IF ( .NOT.  ALLOCATED( q_roof_av ) )  ALLOCATE( q_roof_av(1:surf%ns) )
          CALL rrd_local_2d( q_roof_av )
       CASE ( 'slurb_qs_road_av' )
          IF ( .NOT.  ALLOCATED( qs_road_av ) )  ALLOCATE( qs_road_av(1:surf%ns) )
          CALL rrd_local_2d( qs_road_av )
       CASE ( 'slurb_qs_roof_av' )
          IF ( .NOT.  ALLOCATED( qs_roof_av ) )  ALLOCATE( qs_roof_av(1:surf%ns) )
          CALL rrd_local_2d( qs_roof_av )
       CASE ( 'slurb_qsws_can_av' )
          IF ( .NOT.  ALLOCATED( qsws_can_av ) )  ALLOCATE( qsws_can_av(1:surf%ns) )
          CALL rrd_local_2d( qsws_can_av )
       CASE ( 'slurb_qsws_external_av' )
          IF ( .NOT.  ALLOCATED( qsws_external_av ) )  ALLOCATE( qsws_external_av(1:surf%ns) )
          CALL rrd_local_2d( qsws_external_av )
       CASE ( 'slurb_qsws_road_av' )
          IF ( .NOT.  ALLOCATED( qsws_road_av ) )  ALLOCATE( qsws_road_av(1:surf%ns) )
          CALL rrd_local_2d( qsws_road_av )
       CASE ( 'slurb_qsws_roof_av' )
          IF ( .NOT.  ALLOCATED( qsws_roof_av ) )  ALLOCATE( qsws_roof_av(1:surf%ns) )
          CALL rrd_local_2d( qsws_roof_av )
       CASE ( 'slurb_qsws_urb_av' )
          IF ( .NOT.  ALLOCATED( qsws_urb_av ) )  ALLOCATE( qsws_urb_av(1:surf%ns) )
          CALL rrd_local_2d( qsws_urb_av )
       CASE ( 'slurb_rad_lw_net_road_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_road_av ) )  ALLOCATE( rad_lw_net_road_av(1:surf%ns) )
          CALL rrd_local_2d( rad_lw_net_road_av )
       CASE ( 'slurb_rad_lw_net_roof_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_roof_av ) )  ALLOCATE( rad_lw_net_roof_av(1:surf%ns) )
          CALL rrd_local_2d( rad_lw_net_roof_av )
       CASE ( 'slurb_rad_lw_net_urb_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_urb_av ) )  ALLOCATE( rad_lw_net_urb_av(1:surf%ns) )
          CALL rrd_local_2d( rad_lw_net_urb_av )
       CASE ( 'slurb_rad_lw_net_wall_a_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_wall_a_av ) )  ALLOCATE(rad_lw_net_wall_a_av(1:surf%ns))
          CALL rrd_local_2d( rad_lw_net_wall_a_av )
       CASE ( 'slurb_rad_lw_net_wall_b_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_wall_b_av ) )  ALLOCATE(rad_lw_net_wall_b_av(1:surf%ns))
          CALL rrd_local_2d( rad_lw_net_wall_b_av )
       CASE ( 'slurb_rad_lw_net_win_a_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_win_a_av ) )  ALLOCATE( rad_lw_net_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( rad_lw_net_win_a_av )
       CASE ( 'slurb_rad_lw_net_win_b_av' )
          IF ( .NOT.  ALLOCATED( rad_lw_net_win_b_av ) )  ALLOCATE( rad_lw_net_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( rad_lw_net_win_b_av )
       CASE ( 'slurb_rad_sw_net_road_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_road_av ) )  ALLOCATE( rad_sw_net_road_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_net_road_av )
       CASE ( 'slurb_rad_sw_net_roof_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_roof_av ) )  ALLOCATE( rad_sw_net_roof_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_net_roof_av )
       CASE ( 'slurb_rad_sw_net_urb_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_urb_av ) )  ALLOCATE( rad_sw_net_urb_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_net_urb_av )
       CASE ( 'slurb_rad_sw_net_wall_a_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_wall_a_av ) )  ALLOCATE(rad_sw_net_wall_a_av(1:surf%ns))
          CALL rrd_local_2d( rad_sw_net_wall_a_av )
       CASE ( 'slurb_rad_sw_net_wall_b_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_wall_b_av ) )  ALLOCATE(rad_sw_net_wall_b_av(1:surf%ns))
          CALL rrd_local_2d( rad_sw_net_wall_b_av )
       CASE ( 'slurb_rad_sw_net_win_a_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_win_a_av ) )  ALLOCATE( rad_sw_net_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_net_win_a_av )
       CASE ( 'slurb_rad_sw_net_win_b_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_net_win_b_av ) )  ALLOCATE( rad_sw_net_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_net_win_b_av )
       CASE ( 'slurb_rad_sw_tr_win_a_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_tr_win_a_av ) )  ALLOCATE( rad_sw_tr_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_tr_win_a_av )
       CASE ( 'slurb_rad_sw_tr_win_b_av' )
          IF ( .NOT.  ALLOCATED( rad_sw_tr_win_b_av ) )  ALLOCATE( rad_sw_tr_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( rad_sw_tr_win_b_av )
       CASE ( 'slurb_rah_can_av' )
          IF ( .NOT.  ALLOCATED( rah_can_av ) )  ALLOCATE( rah_can_av(1:surf%ns) )
          CALL rrd_local_2d( rah_can_av )
       CASE ( 'slurb_rah_road_av' )
          IF ( .NOT.  ALLOCATED( rah_road_av ) )  ALLOCATE( rah_road_av(1:surf%ns) )
          CALL rrd_local_2d( rah_road_av )
       CASE ( 'slurb_rah_roof_av' )
          IF ( .NOT.  ALLOCATED( rah_roof_av ) )  ALLOCATE( rah_roof_av(1:surf%ns) )
          CALL rrd_local_2d( rah_roof_av )
       CASE ( 'slurb_rah_wall_a_av' )
          IF ( .NOT.  ALLOCATED( rah_wall_a_av ) )  ALLOCATE( rah_wall_a_av(1:surf%ns) )
          CALL rrd_local_2d( rah_wall_a_av )
       CASE ( 'slurb_rah_wall_b_av' )
          IF ( .NOT.  ALLOCATED( rah_wall_b_av ) )  ALLOCATE( rah_wall_b_av(1:surf%ns) )
          CALL rrd_local_2d( rah_wall_b_av )
       CASE ( 'slurb_rah_win_a_av' )
          IF ( .NOT.  ALLOCATED( rah_win_a_av ) )  ALLOCATE( rah_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( rah_win_a_av )
       CASE ( 'slurb_rah_win_b_av' )
          IF ( .NOT.  ALLOCATED( rah_win_b_av ) )  ALLOCATE( rah_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( rah_win_b_av )
       CASE ( 'slurb_rah_facade_av' )
          IF ( .NOT.  ALLOCATED( rah_facade_av ) )  ALLOCATE( rah_facade_av(1:surf%ns) )
          CALL rrd_local_2d( rah_facade_av )
       CASE ( 'slurb_ram_urb_av' )
          IF ( .NOT.  ALLOCATED( ram_urb_av ) )  ALLOCATE( ram_urb_av(1:surf%ns) )
          CALL rrd_local_2d( ram_urb_av )
       CASE ( 'slurb_rib_can_av' )
          IF ( .NOT.  ALLOCATED( rib_can_av ) )  ALLOCATE( rib_can_av(1:surf%ns) )
          CALL rrd_local_2d( rib_can_av )
       CASE ( 'slurb_rib_road_av' )
          IF ( .NOT.  ALLOCATED( rib_road_av ) )  ALLOCATE( rib_road_av(1:surf%ns) )
          CALL rrd_local_2d( rib_road_av )
       CASE ( 'slurb_rib_roof_av' )
          IF ( .NOT.  ALLOCATED( rib_roof_av ) )  ALLOCATE( rib_roof_av(1:surf%ns) )
          CALL rrd_local_2d( rib_roof_av )
       CASE ( 'slurb_shf_can_av' )
          IF ( .NOT.  ALLOCATED( shf_can_av ) )  ALLOCATE( shf_can_av(1:surf%ns) )
          CALL rrd_local_2d( shf_can_av )
       CASE ( 'slurb_shf_external_av' )
          IF ( .NOT.  ALLOCATED( shf_external_av ) )  ALLOCATE( shf_external_av(1:surf%ns) )
          CALL rrd_local_2d( shf_external_av )
       CASE ( 'slurb_shf_road_av' )
          IF ( .NOT.  ALLOCATED( shf_road_av ) )  ALLOCATE( shf_road_av(1:surf%ns) )
          CALL rrd_local_2d( shf_road_av )
       CASE ( 'slurb_shf_roof_av' )
          IF ( .NOT.  ALLOCATED( shf_roof_av ) )  ALLOCATE( shf_roof_av(1:surf%ns) )
          CALL rrd_local_2d( shf_roof_av )
       CASE ( 'slurb_shf_traffic_av' )
          IF ( .NOT.  ALLOCATED( shf_traffic_av ) )  ALLOCATE( shf_traffic_av(1:surf%ns) )
          CALL rrd_local_2d( shf_traffic_av )
       CASE ( 'slurb_shf_urb_av' )
          IF ( .NOT.  ALLOCATED( shf_urb_av ) )  ALLOCATE( shf_urb_av(1:surf%ns) )
          CALL rrd_local_2d( shf_urb_av )
       CASE ( 'slurb_shf_wall_a_av' )
          IF ( .NOT.  ALLOCATED( shf_wall_a_av ) )  ALLOCATE( shf_wall_a_av(1:surf%ns) )
          CALL rrd_local_2d( shf_wall_a_av )
       CASE ( 'slurb_shf_wall_b_av' )
          IF ( .NOT.  ALLOCATED( shf_wall_b_av ) )  ALLOCATE( shf_wall_b_av(1:surf%ns) )
          CALL rrd_local_2d( shf_wall_b_av )
       CASE ( 'slurb_shf_win_a_av' )
          IF ( .NOT.  ALLOCATED( shf_win_a_av ) )  ALLOCATE( shf_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( shf_win_a_av )
       CASE ( 'slurb_shf_win_b_av' )
          IF ( .NOT.  ALLOCATED( shf_win_b_av ) )  ALLOCATE( shf_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( shf_win_b_av )
       CASE ( 'slurb_t_2m_urb_av' )
          IF ( .NOT.  ALLOCATED( t_2m_urb_av ) )  ALLOCATE( t_2m_urb_av(1:surf%ns) )
          CALL rrd_local_2d( t_2m_urb_av )
       CASE ( 'slurb_t_c_urb_av' )
          IF ( .NOT.  ALLOCATED( t_c_urb_av ) )  ALLOCATE( t_c_urb_av(1:surf%ns) )
          CALL rrd_local_2d( t_c_urb_av )
       CASE ( 'slurb_t_can_av' )
          IF ( .NOT.  ALLOCATED( t_can_av ) )  ALLOCATE( t_can_av(1:surf%ns) )
          CALL rrd_local_2d( t_can_av )
       CASE ( 'slurb_t_h_urb_av' )
          IF ( .NOT.  ALLOCATED( t_h_urb_av ) )  ALLOCATE( t_h_urb_av(1:surf%ns) )
          CALL rrd_local_2d( t_h_urb_av )
       CASE ( 'slurb_t_rad_urb_av' )
          IF ( .NOT.  ALLOCATED( t_rad_urb_av ) )  ALLOCATE( t_rad_urb_av(1:surf%ns) )
          CALL rrd_local_2d( t_rad_urb_av )
       CASE ( 'slurb_t_surf_road_av' )
          IF ( .NOT.  ALLOCATED( t_surf_road_av ) )  ALLOCATE( t_surf_road_av(1:surf%ns) )
          CALL rrd_local_2d( t_surf_road_av )
       CASE ( 'slurb_t_surf_roof_av' )
          IF ( .NOT.  ALLOCATED( t_surf_roof_av ) )  ALLOCATE( t_surf_roof_av(1:surf%ns) )
          CALL rrd_local_2d( t_surf_roof_av )
       CASE ( 'slurb_t_surf_wall_a_av' )
          IF ( .NOT.  ALLOCATED( t_surf_wall_a_av ) )  ALLOCATE( t_surf_wall_a_av(1:surf%ns) )
          CALL rrd_local_2d( t_surf_wall_a_av )
       CASE ( 'slurb_t_surf_wall_b_av' )
          IF ( .NOT.  ALLOCATED( t_surf_wall_b_av ) )  ALLOCATE( t_surf_wall_b_av(1:surf%ns) )
          CALL rrd_local_2d( t_surf_wall_b_av )
       CASE ( 'slurb_t_surf_win_a_av' )
          IF ( .NOT.  ALLOCATED( t_surf_win_a_av ) )  ALLOCATE( t_surf_win_a_av(1:surf%ns) )
          CALL rrd_local_2d( t_surf_win_a_av )
       CASE ( 'slurb_t_surf_win_b_av' )
          IF ( .NOT.  ALLOCATED( t_surf_win_b_av ) )  ALLOCATE( t_surf_win_b_av(1:surf%ns) )
          CALL rrd_local_2d( t_surf_win_b_av )
       CASE ( 'slurb_us_can_av' )
          IF ( .NOT.  ALLOCATED( us_can_av ) )  ALLOCATE( us_can_av(1:surf%ns) )
          CALL rrd_local_2d( us_can_av )
       CASE ( 'slurb_us_road_av' )
          IF ( .NOT.  ALLOCATED( us_road_av ) )  ALLOCATE( us_road_av(1:surf%ns) )
          CALL rrd_local_2d( us_road_av )
       CASE ( 'slurb_us_roof_av' )
          IF ( .NOT.  ALLOCATED( us_roof_av ) )  ALLOCATE( us_roof_av(1:surf%ns) )
          CALL rrd_local_2d( us_roof_av )
       CASE ( 'slurb_us_urb_av' )
          IF ( .NOT.  ALLOCATED( us_urb_av ) )  ALLOCATE( us_urb_av(1:surf%ns) )
          CALL rrd_local_2d( us_urb_av )
       CASE ( 'slurb_usws_urb_av' )
          IF ( .NOT.  ALLOCATED( usws_urb_av ) )  ALLOCATE( usws_urb_av(1:surf%ns) )
          CALL rrd_local_2d( usws_urb_av )
       CASE ( 'slurb_uv_abs_can_av' )
          IF ( .NOT.  ALLOCATED( uv_abs_can_av ) )  ALLOCATE( uv_abs_can_av(1:surf%ns) )
          CALL rrd_local_2d( uv_abs_can_av )
       CASE ( 'slurb_uv_eff_can_av' )
          IF ( .NOT.  ALLOCATED( uv_eff_can_av ) )  ALLOCATE( uv_eff_can_av(1:surf%ns) )
          CALL rrd_local_2d( uv_eff_can_av )
       CASE ( 'slurb_vpt_can_av' )
          IF ( .NOT.  ALLOCATED( vpt_can_av ) )  ALLOCATE( vpt_can_av(1:surf%ns) )
          CALL rrd_local_2d( vpt_can_av )
       CASE ( 'slurb_vpt_road_av' )
          IF ( .NOT.  ALLOCATED( vpt_road_av ) )  ALLOCATE( vpt_road_av(1:surf%ns) )
          CALL rrd_local_2d( vpt_road_av )
       CASE ( 'slurb_vpt_roof_av' )
          IF ( .NOT.  ALLOCATED( vpt_roof_av ) )  ALLOCATE( vpt_roof_av(1:surf%ns) )
          CALL rrd_local_2d( vpt_roof_av )
       CASE ( 'slurb_vsws_urb_av' )
          IF ( .NOT.  ALLOCATED( vsws_urb_av ) )  ALLOCATE( vsws_urb_av(1:surf%ns) )
          CALL rrd_local_2d( vsws_urb_av )
!
!--    Arrays for time averaging (2D, ji-grid).
       CASE ( 'slurb_shf_lsm_av' )
          IF ( .NOT.  ALLOCATED( shf_lsm_av ) )  ALLOCATE( shf_lsm_av(nys:nyn,nxl:nxr) )
          CALL rrd_local_2d_ji( shf_lsm_av )
       CASE ( 'slurb_qsws_lsm_av' )
          IF ( .NOT.  ALLOCATED( qsws_lsm_av ) )  ALLOCATE( qsws_lsm_av(nys:nyn,nxl:nxr) )
          CALL rrd_local_2d_ji( qsws_lsm_av )
!
!--    Arrays for time averaging (3D).
       CASE ( 'slurb_t_road_av' )
          IF ( .NOT.  ALLOCATED( t_road_av ) )  ALLOCATE( t_road_av(nzt_road:nzb_road,1:surf%ns) )
          CALL rrd_local_3d( t_road_av, nzt_road, nzb_road )
       CASE ( 'slurb_t_roof_av' )
          IF ( .NOT.  ALLOCATED( t_roof_av ) )  ALLOCATE( t_roof_av(nzt_roof:nzb_roof,1:surf%ns) )
          CALL rrd_local_3d( t_roof_av, nzt_roof, nzb_roof )
       CASE ( 'slurb_t_wall_a_av' )
          IF ( .NOT.  ALLOCATED( t_wall_a_av ) )  ALLOCATE(t_wall_a_av(nzt_wall:nzb_wall,1:surf%ns))
          CALL rrd_local_3d( t_wall_a_av, nzt_wall, nzb_wall )
       CASE ( 'slurb_t_wall_b_av' )
          IF ( .NOT.  ALLOCATED( t_wall_b_av ) )  ALLOCATE(t_wall_b_av(nzt_wall:nzb_wall,1:surf%ns))
          CALL rrd_local_3d( t_wall_b_av, nzt_wall, nzb_wall )
       CASE ( 'slurb_t_win_a_av' )
          IF ( .NOT.  ALLOCATED( t_win_a_av ) )  ALLOCATE( t_win_a_av(nzt_win:nzb_win,1:surf%ns) )
          CALL rrd_local_3d( t_win_a_av, nzt_win, nzb_win )
       CASE ( 'slurb_t_win_b_av' )
          IF ( .NOT.  ALLOCATED( t_win_b_av ) )  ALLOCATE( t_win_b_av(nzt_win:nzb_win,1:surf%ns) )
          CALL rrd_local_3d( t_win_b_av, nzt_win, nzb_win )

       CASE DEFAULT
          found = .FALSE.

    END SELECT
!
!-- Update the prognostic arrays.
    DO  m = 1, surf%ns
       surf%t_wall_a_p(:,m) = surf%t_wall_a(:,m)
       surf%t_wall_b_p(:,m) = surf%t_wall_b(:,m)
       surf%t_win_a_p(:,m) = surf%t_win_a(:,m)
       surf%t_win_b_p(:,m) = surf%t_win_b(:,m)
       surf%t_roof_p(:,m) = surf%t_roof(:,m)
       surf%t_road_p(:,m) = surf%t_road(:,m)
       surf%t_can_p(m) = surf%t_can(m)
    ENDDO

    IF ( moist_physics )  THEN
       DO  m = 1, surf%ns
          surf%m_liq_roof_p(m) = surf%m_liq_roof(m)
          surf%m_liq_road_p(m) = surf%m_liq_road(m)
          surf%q_can_p(m) = surf%q_can(m)
       ENDDO
    ENDIF

 CONTAINS


 SUBROUTINE rrd_local_2d( var_tgt )

    REAL(wp), DIMENSION(1:surf%ns), INTENT(INOUT) ::  var_tgt  !< targed SLUrb model variable


!
!-- Read and map the former local domain to the new domain (note the inclusion of ghost points).
    IF ( k == 1 )  READ ( 13 ) tmp_2d
    tmp_2d_tgt(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)                                            &
                  = tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!-- Map the new array to SLUrb tiles.
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       IF ( tmp_2d_tgt(j,i) == -9999.0_wp )  CYCLE
       var_tgt(m) = tmp_2d_tgt(j,i)
    ENDDO
    found = .TRUE.

 END SUBROUTINE rrd_local_2d


 SUBROUTINE rrd_local_2d_ji( var_tgt )

    REAL(wp), DIMENSION(nys:nyn,nxl:nxr), INTENT(INOUT) ::  var_tgt  !< targed SLUrb model variable


!
!-- Read and map the former local domain to the new domain.
    IF ( k == 1 )  READ ( 13 ) tmp_2d
    tmp_2d_tgt(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)                                            &
                  = tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!-- Map (j,i) grid with ghost points to (j,i) grid without them.
    var_tgt(:,:) = tmp_2d_tgt(nys:nyn,nxl:nxr)

    found = .TRUE.

 END SUBROUTINE rrd_local_2d_ji


 SUBROUTINE rrd_local_3d( var_tgt, nztl, nzbl )

    INTEGER(iwp) ::  nzbl  !< layer bottom index
    INTEGER(iwp) ::  nztl  !< layer top index

    REAL(wp), DIMENSION(nztl:nzbl,1:surf%ns), INTENT(INOUT) ::  var_tgt  !< targed SLUrb model variable

    REAL(wp), DIMENSION(nztl:nzbl,nys_on_file-nbgp:nyn_on_file+nbgp,                               &
                        nxl_on_file-nbgp:nxr_on_file+nbgp) ::  tmp_3d  !< temporary array for 3D vars

    REAL(wp), DIMENSION(nztl:nzbl,nysg:nyng,nxlg:nxrg) ::  tmp_3d_tgt  !< temporary target array for 3D vars


!
!-- Read and map the former local domain to the new domain.
    IF ( k == 1 )  READ ( 13 ) tmp_3d
    tmp_3d_tgt(nztl:nzbl,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)                                  &
                     = tmp_3d(nztl:nzbl,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
!-- Map the new local array to SLUrb tiles.
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       IF ( ANY ( tmp_3d_tgt(:,j,i) == -9999.0_wp ) )  CYCLE
       var_tgt(:,m) = tmp_3d_tgt(:,j,i)
    ENDDO
    found = .TRUE.

 END SUBROUTINE rrd_local_3d

 END SUBROUTINE slurb_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read restart data using the MPI I/O routines.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_rrd_local_mpi
!
!-- The PALM restart routines assume 2D arrays including ghost points.
    USE indices,                                                                                   &
        ONLY:  nxlg,                                                                               &
               nxrg,                                                                               &
               nysg,                                                                               &
               nyng

    INTEGER(iwp)  ::  i  !< loop index x-direction
    INTEGER(iwp)  ::  j  !< loop index y-direction
    INTEGER(iwp)  ::  k  !< loop index z-direction
    INTEGER(iwp)  ::  m  !< loop index for surface elements on individual surface array

!
!-- No implementation for cyclic fill initialization.
    IF ( cyclic_fill_initialization )  THEN
       WRITE( message_string, * ) 'Cyclic fill initialization is not implemented in SLUrb. ' //    &
                                  'Restart data is not used for initialization of surfaces.'
       CALL message( 'slurb_rrd_local_mpi', 'SLU0039', 0, 1, 0, 6, 0 )
       RETURN
    ENDIF

!
!-- Prognostic 3D arrays
    CALL rrd_local_3d( 'surf_slurb%t_wall_a', surf%t_wall_a, nzt_wall, nzb_wall )
    CALL rrd_local_3d( 'surf_slurb%t_wall_b', surf%t_wall_b, nzt_wall, nzb_wall )
    CALL rrd_local_3d( 'surf_slurb%t_win_a', surf%t_win_a, nzt_win, nzb_win )
    CALL rrd_local_3d( 'surf_slurb%t_win_b', surf%t_win_b, nzt_win, nzb_win )
    CALL rrd_local_3d( 'surf_slurb%t_roof', surf%t_roof, nzt_roof, nzb_roof )
    CALL rrd_local_3d( 'surf_slurb%t_road', surf%t_road, nzt_road, nzb_road )

!
!-- Prognostic 2D arrays.
    CALL rrd_local_2d( 'surf_slurb%us_urb', surf%us_urb )
    CALL rrd_local_2d( 'surf_slurb%shf_urb', surf%shf_urb )
    CALL rrd_local_2d( 'surf_slurb%t_can', surf%t_can )

    IF ( moist_physics )  THEN
       CALL rrd_local_2d( 'surf_slurb%m_liq_roof', surf%m_liq_roof )
       CALL rrd_local_2d( 'surf_slurb%m_liq_road', surf%m_liq_road )
       CALL rrd_local_2d( 'surf_slurb%qsws_road', surf%qsws_road )
       CALL rrd_local_2d( 'surf_slurb%q_can', surf%q_can )
       CALL rrd_local_2d( 'surf_slurb%q_road', surf%q_road )
       CALL rrd_local_2d( 'surf_slurb%q_roof', surf%q_roof )
    ENDIF

    CALL rrd_local_2d( 'surf_slurb%ol_urb', surf%ol_urb )
    CALL rrd_local_2d( 'surf_slurb%ol_roof', surf%ol_roof )
    CALL rrd_local_2d( 'surf_slurb%ol_road', surf%ol_road )
    CALL rrd_local_2d( 'surf_slurb%ol_can', surf%ol_can )
    CALL rrd_local_2d( 'surf_slurb%us_roof', surf%us_roof )
    CALL rrd_local_2d( 'surf_slurb%us_road', surf%us_road )
    CALL rrd_local_2d( 'surf_slurb%uv_eff_can', surf%uv_eff_can )

!
!-- Arrays for time averaging (2D).
    CALL rrd_local_2d_av( 'slurb_albedo_urb_av', albedo_urb_av )
    CALL rrd_local_2d_av( 'slurb_c_liq_road_av', c_liq_road_av )
    CALL rrd_local_2d_av( 'slurb_c_liq_roof_av', c_liq_roof_av )
    CALL rrd_local_2d_av( 'slurb_emiss_urb_av', emiss_urb_av )
    CALL rrd_local_2d_av( 'slurb_ghf_road_av', ghf_road_av )
    CALL rrd_local_2d_av( 'slurb_ghf_roof_av', ghf_roof_av )
    CALL rrd_local_2d_av( 'slurb_ghf_wall_a_av', ghf_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_ghf_wall_b_av', ghf_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_ghf_win_a_av', ghf_win_a_av )
    CALL rrd_local_2d_av( 'slurb_ghf_win_b_av', ghf_win_b_av )
    CALL rrd_local_2d_av( 'slurb_m_liq_road_av', m_liq_road_av )
    CALL rrd_local_2d_av( 'slurb_m_liq_roof_av', m_liq_roof_av )
    CALL rrd_local_2d_av( 'slurb_ol_can_av', ol_can_av )
    CALL rrd_local_2d_av( 'slurb_ol_road_av', ol_road_av )
    CALL rrd_local_2d_av( 'slurb_ol_roof_av', ol_roof_av )
    CALL rrd_local_2d_av( 'slurb_pt_can_av', pt_can_av )
    CALL rrd_local_2d_av( 'slurb_pt_road_av', pt_road_av )
    CALL rrd_local_2d_av( 'slurb_pt_roof_av', pt_roof_av )
    CALL rrd_local_2d_av( 'slurb_pt_wall_a_av', pt_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_pt_wall_b_av', pt_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_pt_win_a_av', pt_win_a_av )
    CALL rrd_local_2d_av( 'slurb_pt_win_b_av', pt_win_b_av )
    CALL rrd_local_2d_av( 'slurb_q_can_av', q_can_av )
    CALL rrd_local_2d_av( 'slurb_q_road_av', q_road_av )
    CALL rrd_local_2d_av( 'slurb_q_roof_av', q_roof_av )
    CALL rrd_local_2d_av( 'slurb_qs_road_av', qs_road_av )
    CALL rrd_local_2d_av( 'slurb_qs_roof_av', qs_roof_av )
    CALL rrd_local_2d_av( 'slurb_qsws_can_av', qsws_can_av )
    CALL rrd_local_2d_av( 'slurb_qsws_external_av', qsws_external_av )
    CALL rrd_local_2d_av( 'slurb_qsws_road_av', qsws_road_av )
    CALL rrd_local_2d_av( 'slurb_qsws_roof_av', qsws_roof_av )
    CALL rrd_local_2d_av( 'slurb_qsws_urb_av', qsws_urb_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_road_av', rad_lw_net_road_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_roof_av', rad_lw_net_roof_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_urb_av', rad_lw_net_urb_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_wall_a_av', rad_lw_net_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_wall_b_av', rad_lw_net_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_win_a_av', rad_lw_net_win_a_av )
    CALL rrd_local_2d_av( 'slurb_rad_lw_net_win_b_av', rad_lw_net_win_b_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_road_av', rad_sw_net_road_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_roof_av', rad_sw_net_roof_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_urb_av', rad_sw_net_urb_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_wall_a_av', rad_sw_net_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_wall_b_av', rad_sw_net_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_win_a_av', rad_sw_net_win_a_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_net_win_b_av', rad_sw_net_win_b_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_tr_win_a_av', rad_sw_tr_win_a_av )
    CALL rrd_local_2d_av( 'slurb_rad_sw_tr_win_b_av', rad_sw_tr_win_b_av )
    CALL rrd_local_2d_av( 'slurb_rah_can_av', rah_can_av )
    CALL rrd_local_2d_av( 'slurb_rah_road_av', rah_road_av )
    CALL rrd_local_2d_av( 'slurb_rah_roof_av', rah_roof_av )
    CALL rrd_local_2d_av( 'slurb_rah_wall_a_av', rah_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_rah_wall_b_av', rah_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_rah_win_a_av', rah_win_a_av )
    CALL rrd_local_2d_av( 'slurb_rah_win_b_av', rah_win_b_av )
    CALL rrd_local_2d_av( 'slurb_rah_facade_av', rah_facade_av )
    CALL rrd_local_2d_av( 'slurb_ram_urb_av', ram_urb_av )
    CALL rrd_local_2d_av( 'slurb_rib_can_av', rib_can_av )
    CALL rrd_local_2d_av( 'slurb_rib_road_av', rib_road_av )
    CALL rrd_local_2d_av( 'slurb_rib_roof_av', rib_roof_av )
    CALL rrd_local_2d_av( 'slurb_shf_can_av', shf_can_av )
    CALL rrd_local_2d_av( 'slurb_shf_external_av', shf_external_av )
    CALL rrd_local_2d_av( 'slurb_shf_road_av', shf_road_av )
    CALL rrd_local_2d_av( 'slurb_shf_roof_av', shf_roof_av )
    CALL rrd_local_2d_av( 'slurb_shf_traffic_av', shf_traffic_av )
    CALL rrd_local_2d_av( 'slurb_shf_urb_av', shf_urb_av )
    CALL rrd_local_2d_av( 'slurb_shf_wall_a_av', shf_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_shf_wall_b_av', shf_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_shf_win_a_av', shf_win_a_av )
    CALL rrd_local_2d_av( 'slurb_shf_win_b_av', shf_win_b_av )
    CALL rrd_local_2d_av( 'slurb_t_2m_urb_av', t_2m_urb_av )
    CALL rrd_local_2d_av( 'slurb_t_c_urb_av', t_c_urb_av )
    CALL rrd_local_2d_av( 'slurb_t_can_av', t_can_av )
    CALL rrd_local_2d_av( 'slurb_t_h_urb_av', t_h_urb_av )
    CALL rrd_local_2d_av( 'slurb_t_rad_urb_av', t_rad_urb_av )
    CALL rrd_local_2d_av( 'slurb_t_surf_road_av', t_surf_road_av )
    CALL rrd_local_2d_av( 'slurb_t_surf_roof_av', t_surf_roof_av )
    CALL rrd_local_2d_av( 'slurb_t_surf_wall_a_av', t_surf_wall_a_av )
    CALL rrd_local_2d_av( 'slurb_t_surf_wall_b_av', t_surf_wall_b_av )
    CALL rrd_local_2d_av( 'slurb_t_surf_win_a_av', t_surf_win_a_av )
    CALL rrd_local_2d_av( 'slurb_t_surf_win_b_av', t_surf_win_b_av )
    CALL rrd_local_2d_av( 'slurb_us_can_av', us_can_av )
    CALL rrd_local_2d_av( 'slurb_us_road_av', us_road_av )
    CALL rrd_local_2d_av( 'slurb_us_roof_av', us_roof_av )
    CALL rrd_local_2d_av( 'slurb_us_urb_av', us_urb_av )
    CALL rrd_local_2d_av( 'slurb_usws_urb_av', usws_urb_av )
    CALL rrd_local_2d_av( 'slurb_uv_abs_can_av', uv_abs_can_av )
    CALL rrd_local_2d_av( 'slurb_uv_eff_can_av', uv_eff_can_av )
    CALL rrd_local_2d_av( 'slurb_vpt_can_av', vpt_can_av )
    CALL rrd_local_2d_av( 'slurb_vpt_road_av', vpt_road_av )
    CALL rrd_local_2d_av( 'slurb_vpt_roof_av', vpt_roof_av )
    CALL rrd_local_2d_av( 'slurb_vsws_urb_av', vsws_urb_av )

!
!-- Arrays for time averaging (2D, ji-grid).
    CALL rrd_local_2d_ji_av( 'slurb_shf_lsm_av', shf_lsm_av )
    CALL rrd_local_2d_ji_av( 'slurb_qsws_lsm_av', qsws_lsm_av )

!
!-- Arrays for time averaging (3D).
    CALL rrd_local_3d_av( 'slurb_t_road_av', t_road_av, nzt_road, nzb_road )
    CALL rrd_local_3d_av( 'slurb_t_roof_av', t_roof_av, nzt_roof, nzb_roof)
    CALL rrd_local_3d_av( 'slurb_t_wall_a_av', t_wall_a_av, nzt_wall, nzb_wall )
    CALL rrd_local_3d_av( 'slurb_t_wall_b_av', t_wall_b_av, nzt_wall, nzb_wall )
    CALL rrd_local_3d_av( 'slurb_t_win_a_av', t_win_a_av, nzt_win, nzb_win )
    CALL rrd_local_3d_av( 'slurb_t_win_b_av', t_win_b_av, nzt_win, nzb_win )

!
!-- Update the prognostic arrays.
    DO  m = 1, surf%ns
       surf%t_wall_a_p(:,m) = surf%t_wall_a(:,m)
       surf%t_wall_b_p(:,m) = surf%t_wall_b(:,m)
       surf%t_win_a_p(:,m) = surf%t_win_a(:,m)
       surf%t_win_b_p(:,m) = surf%t_win_b(:,m)
       surf%t_roof_p(:,m) = surf%t_roof(:,m)
       surf%t_road_p(:,m) = surf%t_road(:,m)
       surf%t_can_p(m) = surf%t_can(m)
    ENDDO

    IF ( moist_physics )  THEN
       DO  m = 1, surf%ns
          surf%m_liq_roof_p(m) = surf%m_liq_roof(m)
          surf%m_liq_road_p(m) = surf%m_liq_road(m)
          surf%q_can_p(m) = surf%q_can(m)
       ENDDO
    ENDIF

 CONTAINS


 SUBROUTINE rrd_local_2d( varname, tgt )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), DIMENSION(1:surf%ns), INTENT(INOUT) ::  tgt  !< target array of the variable

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_2d  !<


    CALL rrd_mpi_io( varname, tmp_2d )
!
!-- Map the variable back from the (j,i) grid to SLUrb grid.
    tgt(:) = -9999.0_wp
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       tgt(m) = tmp_2d(j,i)
    ENDDO

 END SUBROUTINE rrd_local_2d


 SUBROUTINE rrd_local_2d_av( varname, tgt )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) ::  tgt  !< target array of the variable

    LOGICAL ::  found

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_2d  !<


!
!-- Check if the averaged variable exists, if not, simply return.
    CALL rd_mpi_io_check_array( varname , found = found )
    IF ( .NOT.  found )  RETURN
!
!-  If found, allocate if needed.
    IF ( .NOT. ALLOCATED( tgt ) )  ALLOCATE( tgt(1:surf%ns) )

    CALL rrd_mpi_io( varname, tmp_2d )
!
!-- Map the variable back from the (j,i) grid to SLUrb grid.
    tgt(:) = -9999.0_wp
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       tgt(m) = tmp_2d(j,i)
    ENDDO

 END SUBROUTINE rrd_local_2d_av


 SUBROUTINE rrd_local_2d_ji_av( varname, tgt )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  tgt  !< target array of the variable

    LOGICAL ::  found

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_2d  !<


!
!-- Check if the averaged variable exists, if not, simply return.
    CALL rd_mpi_io_check_array( varname , found = found )
    IF ( .NOT.  found )  RETURN
!
!-  If found, allocate if needed.
    IF ( .NOT. ALLOCATED( tgt ) )  ALLOCATE( tgt(nys:nyn,nxl:nxr) )

    CALL rrd_mpi_io( varname, tmp_2d )
!
!-- Map the variable back from the (j,i) grid to SLUrb grid.
    tgt(:,:) = tmp_2d(nys:nyn,nxl:nxr)

 END SUBROUTINE rrd_local_2d_ji_av


 SUBROUTINE rrd_local_3d( varname, tgt, nztl, nzbl )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), DIMENSION(nztl:nzbl,1:surf%ns), INTENT(OUT) ::  tgt  !< target array of the variable

    CHARACTER(LEN=3) ::  id  !< layer identifier

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_3dto2d  !< array to temporarily map 3D array from 2D layers

    REAL(wp), DIMENSION(nztl:nzbl,nysg:nyng,nxlg:nxrg) :: tmp_3d


!
!-- Read the restart variable layer-by-layer.
    DO  k = nztl, nzbl
       WRITE( id, '(I3.3)') k
       CALL rrd_mpi_io( TRIM( varname ) // '_' // id, tmp_3dto2d )
       tmp_3d(k,:,:) = tmp_3dto2d(:,:)
    ENDDO
!
!-- Map the variable back from the (j,i) grid to SLUrb grid.
    tgt(:,:) = -9999.0_wp
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       tgt(:,m) = tmp_3d(:,j,i)
    ENDDO

 END SUBROUTINE rrd_local_3d


 SUBROUTINE rrd_local_3d_av( varname, tgt, nztl, nzbl )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) ::  tgt  !< target array of the variable

    CHARACTER(LEN=3) ::  id  !< layer identifier

    LOGICAL ::  found  !<

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) :: tmp_3dto2d  !< array to temporarily map 3D array from 2D layers

    REAL(wp), DIMENSION(nztl:nzbl,nysg:nyng,nxlg:nxrg) :: tmp_3d  !<


!
!-- Check if the averaged variable exists, if not, simply return.
    CALL rd_mpi_io_check_array( varname , found = found )
    IF ( .NOT.  found )  RETURN
!
!-  If found, allocate if needed.
    IF ( .NOT. ALLOCATED( tgt ) )  ALLOCATE( tgt(nztl:nzbl,1:surf%ns) )

!
!-- Read the restart variable layer-by-layer.
    DO  k = nztl, nzbl
       WRITE( id, '(I3.3)') k
       CALL rrd_mpi_io( TRIM( varname ) // '_' // id, tmp_3dto2d )
       tmp_3d(k,:,:) = tmp_3dto2d(:,:)
    ENDDO
!
!-- Map the variable back from the (j,i) grid to SLUrb grid.
    tgt(:,:) = -9999.0_wp
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       tgt(:,m) = tmp_3d(:,j,i)
    ENDDO

 END SUBROUTINE rrd_local_3d_av

 END SUBROUTINE slurb_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write restart data using either the legacy Fortran I/O or or MPI I/O. Only minimum necessary
!> set of variables are stored for the restart, namely those which are used at next time step
!> before re-computation and cannot be inferred from other variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_wrd_local

    USE control_parameters,                                                                        &
        ONLY:  restart_data_format_output
!
!-- The PALM restart routines assume 2D arrays including ghost points.
    USE indices,                                                                                   &
        ONLY:  nxlg,                                                                               &
               nxrg,                                                                               &
               nysg,                                                                               &
               nyng

    INTEGER(iwp)  ::  i  !< loop index x-direction
    INTEGER(iwp)  ::  j  !< loop index y-direction
    INTEGER(iwp)  ::  k  !< loop index z-direction
    INTEGER(iwp)  ::  m  !< loop index for surface elements on individual surface array

!
!-- Prognostic 3D arrays.
    CALL wrd_local_3d( 'surf_slurb%t_wall_a', surf%t_wall_a, nzt_wall, nzb_wall )
    CALL wrd_local_3d( 'surf_slurb%t_wall_b', surf%t_wall_b, nzt_wall, nzb_wall )
    CALL wrd_local_3d( 'surf_slurb%t_win_a', surf%t_wall_b, nzt_win, nzb_win )
    CALL wrd_local_3d( 'surf_slurb%t_win_b', surf%t_wall_b, nzt_win, nzb_win )
    CALL wrd_local_3d( 'surf_slurb%t_roof', surf%t_roof, nzt_roof, nzb_roof )
    CALL wrd_local_3d( 'surf_slurb%t_road', surf%t_road, nzt_road, nzb_road )

!
!-- Prognostic 2D arrays.
    CALL wrd_local_2d( 'surf_slurb%us_urb', surf%us_urb )
    CALL wrd_local_2d( 'surf_slurb%shf_urb', surf%shf_urb )
    CALL wrd_local_2d( 'surf_slurb%t_can', surf%t_can )

    IF ( moist_physics )  THEN
       CALL wrd_local_2d( 'surf_slurb%m_liq_roof', surf%m_liq_roof )
       CALL wrd_local_2d( 'surf_slurb%m_liq_road', surf%m_liq_road )
       CALL wrd_local_2d( 'surf_slurb%qsws_road', surf%qsws_road )
       CALL wrd_local_2d( 'surf_slurb%q_can', surf%q_can )
       CALL wrd_local_2d( 'surf_slurb%q_road', surf%q_road )
       CALL wrd_local_2d( 'surf_slurb%q_roof', surf%q_roof )
    ENDIF

!
!-- Previous Obukhov lengths are used as the initial guess in slurb_resistance_model when computing
!-- one for the next time step, so save them to ensure continuity from run to run.
    CALL wrd_local_2d( 'surf_slurb%ol_urb', surf%ol_urb )
    CALL wrd_local_2d( 'surf_slurb%ol_roof', surf%ol_roof )
    CALL wrd_local_2d( 'surf_slurb%ol_road', surf%ol_road )
    CALL wrd_local_2d( 'surf_slurb%ol_can', surf%ol_can )

!
!-- Previous local friction veloctiy is used for Kanda et al. (2007) z0h parametrizations
!-- in slurb_resistance_model before updated in the same routine.
    CALL wrd_local_2d( 'surf_slurb%us_roof', surf%us_roof )
    CALL wrd_local_2d( 'surf_slurb%us_road', surf%us_road )
!
!-- Previous effective canyon wind speed is used in slurb_resistance_model before being computed in
!-- slurb_canyon_model for the next time step.
    CALL wrd_local_2d( 'surf_slurb%uv_eff_can', surf%uv_eff_can )

!
!-- Arrays for time averaging (2D). Prefix slurb_ is added to prevent possible name conflicts.
    IF ( ALLOCATED( albedo_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_albedo_urb_av', albedo_urb_av )
    ENDIF
    IF ( ALLOCATED( c_liq_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_c_liq_road_av', c_liq_road_av )
    ENDIF
    IF ( ALLOCATED( c_liq_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_c_liq_roof_av', c_liq_roof_av )
    ENDIF
    IF ( ALLOCATED( emiss_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_emiss_urb_av', emiss_urb_av )
    ENDIF
    IF ( ALLOCATED( ghf_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ghf_road_av', ghf_road_av )
    ENDIF
    IF ( ALLOCATED( ghf_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ghf_roof_av', ghf_roof_av )
    ENDIF
    IF ( ALLOCATED( ghf_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ghf_wall_a_av', ghf_wall_a_av )
    ENDIF
    IF ( ALLOCATED( ghf_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ghf_wall_b_av', ghf_wall_b_av )
    ENDIF
    IF ( ALLOCATED( ghf_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ghf_win_a_av', ghf_win_a_av )
    ENDIF
    IF ( ALLOCATED( ghf_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ghf_win_b_av', ghf_win_b_av )
    ENDIF
    IF ( ALLOCATED( m_liq_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_m_liq_road_av', m_liq_road_av )
    ENDIF
    IF ( ALLOCATED( m_liq_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_m_liq_roof_av', m_liq_roof_av )
    ENDIF
    IF ( ALLOCATED( ol_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ol_can_av', ol_can_av )
    ENDIF
    IF ( ALLOCATED( ol_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ol_road_av', ol_road_av )
    ENDIF
    IF ( ALLOCATED( ol_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ol_roof_av', ol_roof_av )
    ENDIF
    IF ( ALLOCATED( pt_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_can_av', pt_can_av )
    ENDIF
    IF ( ALLOCATED( pt_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_road_av', pt_road_av )
    ENDIF
    IF ( ALLOCATED( pt_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_roof_av', pt_roof_av )
    ENDIF
    IF ( ALLOCATED( pt_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_wall_a_av', pt_wall_a_av )
    ENDIF
    IF ( ALLOCATED( pt_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_wall_b_av', pt_wall_b_av )
    ENDIF
    IF ( ALLOCATED( pt_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_win_a_av', pt_win_a_av )
    ENDIF
    IF ( ALLOCATED( pt_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_pt_win_b_av', pt_win_b_av )
    ENDIF
    IF ( ALLOCATED( q_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_q_can_av', q_can_av )
    ENDIF
    IF ( ALLOCATED( q_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_q_road_av', q_road_av )
    ENDIF
    IF ( ALLOCATED( q_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_q_roof_av', q_roof_av )
    ENDIF
    IF ( ALLOCATED( qs_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qs_road_av', qs_road_av )
    ENDIF
    IF ( ALLOCATED( qs_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qs_roof_av', qs_roof_av )
    ENDIF
    IF ( ALLOCATED( qsws_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qsws_can_av', qsws_can_av )
    ENDIF
    IF ( ALLOCATED( qsws_external_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qsws_external_av', qsws_external_av )
    ENDIF
    IF ( ALLOCATED( qsws_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qsws_road_av', qsws_road_av )
    ENDIF
    IF ( ALLOCATED( qsws_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qsws_roof_av', qsws_roof_av )
    ENDIF
    IF ( ALLOCATED( qsws_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_qsws_urb_av', qsws_urb_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_road_av', rad_lw_net_road_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_roof_av', rad_lw_net_roof_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_urb_av', rad_lw_net_urb_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_wall_a_av', rad_lw_net_wall_a_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_wall_b_av', rad_lw_net_wall_b_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_win_a_av', rad_lw_net_win_a_av )
    ENDIF
    IF ( ALLOCATED( rad_lw_net_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_lw_net_win_b_av', rad_lw_net_win_b_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_road_av', rad_sw_net_road_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_roof_av', rad_sw_net_roof_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_urb_av', rad_sw_net_urb_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_wall_a_av', rad_sw_net_wall_a_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_wall_b_av', rad_sw_net_wall_b_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_win_a_av', rad_sw_net_win_a_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_net_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_net_win_b_av', rad_sw_net_win_b_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_tr_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_tr_win_a_av', rad_sw_tr_win_a_av )
    ENDIF
    IF ( ALLOCATED( rad_sw_tr_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rad_sw_tr_win_b_av', rad_sw_tr_win_b_av )
    ENDIF
    IF ( ALLOCATED( rah_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_can_av', rah_can_av )
    ENDIF
    IF ( ALLOCATED( rah_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_road_av', rah_road_av )
    ENDIF
    IF ( ALLOCATED( rah_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_roof_av', rah_roof_av )
    ENDIF
    IF ( ALLOCATED( rah_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_wall_a_av', rah_wall_a_av )
    ENDIF
    IF ( ALLOCATED( rah_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_wall_b_av', rah_wall_b_av )
    ENDIF
    IF ( ALLOCATED( rah_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_win_a_av', rah_win_a_av )
    ENDIF
    IF ( ALLOCATED( rah_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_win_b_av', rah_win_b_av )
    ENDIF
    IF ( ALLOCATED( rah_facade_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rah_facade_av', rah_facade_av )
    ENDIF
    IF ( ALLOCATED( ram_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_ram_urb_av', ram_urb_av )
    ENDIF
    IF ( ALLOCATED( rib_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rib_can_av', rib_can_av )
    ENDIF
    IF ( ALLOCATED( rib_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rib_road_av', rib_road_av )
    ENDIF
    IF ( ALLOCATED( rib_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_rib_roof_av', rib_roof_av )
    ENDIF
    IF ( ALLOCATED( shf_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_can_av', shf_can_av )
    ENDIF
    IF ( ALLOCATED( shf_external_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_external_av', shf_external_av )
    ENDIF
    IF ( ALLOCATED( shf_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_road_av', shf_road_av )
    ENDIF
    IF ( ALLOCATED( shf_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_roof_av', shf_roof_av )
    ENDIF
    IF ( ALLOCATED( shf_traffic_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_traffic_av', shf_traffic_av )
    ENDIF
    IF ( ALLOCATED( shf_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_urb_av', shf_urb_av )
    ENDIF
    IF ( ALLOCATED( shf_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_wall_a_av', shf_wall_a_av )
    ENDIF
    IF ( ALLOCATED( shf_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_wall_b_av', shf_wall_b_av )
    ENDIF
    IF ( ALLOCATED( shf_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_win_a_av', shf_win_a_av )
    ENDIF
    IF ( ALLOCATED( shf_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_shf_win_b_av', shf_win_b_av )
    ENDIF
    IF ( ALLOCATED( t_2m_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_2m_urb_av', t_2m_urb_av )
    ENDIF
    IF ( ALLOCATED( t_c_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_c_urb_av', t_c_urb_av )
    ENDIF
    IF ( ALLOCATED( t_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_can_av', t_can_av )
    ENDIF
    IF ( ALLOCATED( t_h_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_h_urb_av', t_h_urb_av )
    ENDIF
    IF ( ALLOCATED( t_rad_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_rad_urb_av', t_rad_urb_av )
    ENDIF
    IF ( ALLOCATED( t_surf_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_surf_road_av', t_surf_road_av )
    ENDIF
    IF ( ALLOCATED( t_surf_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_surf_roof_av', t_surf_roof_av )
    ENDIF
    IF ( ALLOCATED( t_surf_wall_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_surf_wall_a_av', t_surf_wall_a_av )
    ENDIF
    IF ( ALLOCATED( t_surf_wall_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_surf_wall_b_av', t_surf_wall_b_av )
    ENDIF
    IF ( ALLOCATED( t_surf_win_a_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_surf_win_a_av', t_surf_win_a_av )
    ENDIF
    IF ( ALLOCATED( t_surf_win_b_av ) )  THEN
       CALL wrd_local_2d( 'slurb_t_surf_win_b_av', t_surf_win_b_av )
    ENDIF
    IF ( ALLOCATED( us_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_us_can_av', us_can_av )
    ENDIF
    IF ( ALLOCATED( us_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_us_road_av', us_road_av )
    ENDIF
    IF ( ALLOCATED( us_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_us_roof_av', us_roof_av )
    ENDIF
    IF ( ALLOCATED( us_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_us_urb_av', us_urb_av )
    ENDIF
    IF ( ALLOCATED( usws_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_usws_urb_av', usws_urb_av )
    ENDIF
    IF ( ALLOCATED( uv_abs_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_uv_abs_can_av', uv_abs_can_av )
    ENDIF
    IF ( ALLOCATED( uv_eff_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_uv_eff_can_av', uv_eff_can_av )
    ENDIF
    IF ( ALLOCATED( vpt_can_av ) )  THEN
       CALL wrd_local_2d( 'slurb_vpt_can_av', vpt_can_av )
    ENDIF
    IF ( ALLOCATED( vpt_road_av ) )  THEN
       CALL wrd_local_2d( 'slurb_vpt_road_av', vpt_road_av )
    ENDIF
    IF ( ALLOCATED( vpt_roof_av ) )  THEN
       CALL wrd_local_2d( 'slurb_vpt_roof_av', vpt_roof_av )
    ENDIF
    IF ( ALLOCATED( vsws_urb_av ) )  THEN
       CALL wrd_local_2d( 'slurb_vsws_urb_av', vsws_urb_av )
    ENDIF

!
!-- Arrays for time averaging (2D, ji-grid).
    IF ( ALLOCATED( shf_lsm_av ) )  THEN
       CALL wrd_local_2d_ji( 'slurb_shf_lsm_av', shf_lsm_av )
    ENDIF
    IF ( ALLOCATED( qsws_lsm_av ) )  THEN
       CALL wrd_local_2d_ji( 'slurb_qsws_lsm_av', qsws_lsm_av )
    ENDIF

!
!-- Arrays for time averaging (3D).
    IF ( ALLOCATED( t_road_av ) ) THEN
       CALL wrd_local_3d( 'slurb_t_road_av', t_road_av, nzt_road, nzb_road )
    ENDIF
    IF ( ALLOCATED( t_roof_av ) ) THEN
       CALL wrd_local_3d( 'slurb_t_roof_av', t_roof_av, nzt_roof, nzb_roof )
    ENDIF
    IF ( ALLOCATED( t_wall_a_av ) ) THEN
       CALL wrd_local_3d( 'slurb_t_wall_a_av', t_wall_a_av, nzt_wall, nzb_wall )
    ENDIF
    IF ( ALLOCATED( t_wall_b_av ) ) THEN
       CALL wrd_local_3d( 'slurb_t_wall_b_av', t_wall_b_av, nzt_wall, nzb_wall )
    ENDIF
    IF ( ALLOCATED( t_win_a_av ) ) THEN
       CALL wrd_local_3d( 'slurb_t_win_a_av', t_win_a_av, nzt_win, nzb_win )
    ENDIF
    IF ( ALLOCATED( t_win_b_av ) ) THEN
       CALL wrd_local_3d( 'slurb_t_win_b_av', t_win_b_av, nzt_win, nzb_win )
    ENDIF

 CONTAINS


 SUBROUTINE wrd_local_2d( varname, src )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), INTENT(IN), DIMENSION(:) ::  src  !<
!
!-- For restarts, the arrays are mapped to (nysg:nyng,nxlg:nxrg) grid (includes ghost points),
!-- as there is no implementation for wrd_mpi_io for (nys:nyn,nxl:nxr) grid. Technically, the
!-- grid without the ghost points could be used for Fortran I/O, but this would require splitting
!-- the mapping routines (i.e. this and following subroutines) between the restart data formats.
    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_2d  !<


!
!-- Map the SLUrb grid on a regular (j,i) grid.
    tmp_2d(:,:) = -9999.0_wp
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       tmp_2d(j,i) = src(m)
    ENDDO

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
       CALL wrd_write_string( TRIM( varname ) )
       WRITE ( 14 ) tmp_2d
    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
       CALL wrd_mpi_io( TRIM( varname ), tmp_2d )
    ENDIF

 END SUBROUTINE wrd_local_2d


 SUBROUTINE wrd_local_2d_ji( varname, src )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< name of the variable to be written

    REAL(wp), INTENT(IN), DIMENSION(:,:) ::  src  !<

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_2d  !<


!
!-- Map to (j,i) grid to (j,i) grid including ghost points.
    tmp_2d(:,:) = -9999.0_wp
    tmp_2d(nys:nyn,nxl:nxr) = src(:,:)

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
       CALL wrd_write_string( TRIM( varname ) )
       WRITE ( 14 ) tmp_2d
    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
       CALL wrd_mpi_io( TRIM( varname ), tmp_2d )
    ENDIF

 END SUBROUTINE wrd_local_2d_ji


 SUBROUTINE wrd_local_3d( varname, src, nztl, nzbl )

    CHARACTER(LEN=*), INTENT(IN) ::  varname   !< name of the variable to be written

    INTEGER(iwp) ::  nzbl  !< layer bottom index
    INTEGER(iwp) ::  nztl  !< layer top index

    REAL(wp), INTENT(IN), DIMENSION(:,:) ::  src  !< source array of the variable

    CHARACTER(LEN=3) ::  id  !< layer identifier

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp_3dto2d  !< array to temporarily map 3D array layers into 2D

    REAL(wp), DIMENSION(nztl:nzbl,nysg:nyng,nxlg:nxrg) ::  tmp_3d  !<


!
!-- Map the SLUrb grid on a regular (j,i) grid.
    tmp_3d(:,:,:) = -9999.0_wp
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       tmp_3d(:,j,i) = src(:,m)
    ENDDO

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
       CALL wrd_write_string( TRIM( varname ) )
       WRITE ( 14 ) tmp_3d
!
!-- With MPI, 3D variables are written as 2D arrays layer-by-layer with layer id as suffix.
!-- This is due to the fact that wrd_mpi_io is implemented only for 3D arrays covering the whole
!-- vertical extent of the domain, not arbitrary number of layers.
    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
       DO  k = nztl, nzbl
          WRITE( id, '(I3.3)') k
          tmp_3dto2d(:,:) = tmp_3d(k,:,:)
          CALL wrd_mpi_io( TRIM( varname ) // '_' // id, tmp_3dto2d )
       ENDDO
    ENDIF

 END SUBROUTINE wrd_local_3d

 END SUBROUTINE slurb_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates model external variables, e.g. the variables defined at the first atmospheric level based
!> on atmospheric simulation state as well as the temporally dynamic SLUrb input variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_update_external_vars

    INTEGER(iwp) ::  i      !< loop index
    INTEGER(iwp) ::  j      !< loop index
    INTEGER(iwp) ::  k_atm  !< k index of the first atmospheric level
    INTEGER(iwp) ::  m      !< loop index of surface tiles
    INTEGER(iwp) ::  t      !< current timestep index
    INTEGER(iwp) ::  tm     !< previous timestep index

    REAL(wp) ::  fac_dt  !< factor for linear interpolation between timesteps
    REAL(wp) ::  vtws    !< buoyancy flux
    REAL(wp) ::  ws      !< free convection velocity scale


    DO  m = 1, surf%ns

       i = surf%i(m)
       j = surf%j(m)
       k_atm = topo_top_ind(j,i,0) + 1
!
!--    Calculate the pt, vpt and q for atmosphere depending on what modules are enabled.
       IF ( bulk_cloud_model )  THEN
          surf%pt1(m) = pt(k_atm,j,i) + lv_d_cp * d_exner(k_atm) * ql(k_atm,j,i)
          surf%q1(m) = q(k_atm,j,i) - ql(k_atm,j,i)
          surf%vpt1(m) = surf%pt1(m) * ( 1.0_wp + 0.61_wp * surf%q1(m) )
       ELSEIF ( cloud_droplets )  THEN
          surf%pt1(m) = pt(k_atm,j,i) + lv_d_cp * d_exner(k_atm) * ql(k_atm,j,i)
          surf%q1(m) = q(k_atm,j,i)
          surf%vpt1(m) = surf%pt1(m) * ( 1.0_wp + 0.61_wp * surf%q1(m) )
       ELSE
          surf%pt1(m) = pt(k_atm,j,i)
          IF ( moist_physics )  THEN
             surf%q1(m) = q(k_atm,j,i)
             surf%vpt1(m) = surf%pt1(m) * ( 1.0_wp + 0.61_wp * surf%q1(m) )
          ENDIF
       ENDIF

       surf%uv_abs1(m) = SQRT( ( 0.5 * ( u(k_atm,j,i) + u(k_atm,j,i+1) ) )**2 +                    &
                               ( 0.5 * ( v(k_atm,j,i) + v(k_atm,j+1,i) ) )**2 )

!
!--    Calculate surface-parallel absolute velocity uv_eff1 at cell center using
!--    free convection scale (w_star, for unstable cases).
       vtws = surf%shf_urb(m) + MERGE( lv_d_cp * surf%qsws_urb(m), 0.0_wp, moist_physics )
!
!--    No scaling for stable cases:
       vtws = MERGE( vtws, 0.0_wp, vtws > 0.0_wp )
       ws = ( g / surf%pt1(m) * surf%z_mo(m) * vtws )**( 1.0_wp / 3.0_wp )

       surf%uv_eff1(m) = SQRT( ( 0.5 * ( u(k_atm,j,i) + u(k_atm,j,i+1) ) )**2 +                    &
                               ( 0.5 * ( v(k_atm,j,i) + v(k_atm,j+1,i) ) )**2 + ws**2 )
    ENDDO

    IF ( slurb_dynamic%ntime > 0 )  CALL update_dynamic_inputs

 CONTAINS


 SUBROUTINE update_dynamic_inputs


!
!-- Update temporally dynamic input variables using linear interpolation.
    IF ( time_since_reference_point  <= slurb_dynamic%time(1) )  THEN
!
!--    If the current time is before the first time step of the input, which could happen for
!--    example in the case the input doesn't cover spinup, use the first time step.
       t      = 1
       tm     = 1
       fac_dt = 0.0_wp
    ELSEIF ( time_since_reference_point >= slurb_dynamic%time(slurb_dynamic%ntime) )  THEN
!
!--    Similarly, use constant extrapolation in the case the input file doesn't cover the whole
!--    simulation period. The user has been warned about this in initialization.
       t      = slurb_dynamic%ntime
       tm     = slurb_dynamic%ntime
       fac_dt = 0.0_wp
    ELSE
!
!--    Compute factor for linear interpolation (weighted average).
       t = 0
       DO WHILE ( slurb_dynamic%time(t) <= time_since_reference_point )
          t = t + 1
       ENDDO

       tm = MAX( t-1, 0 )

       fac_dt = ( time_since_reference_point - slurb_dynamic%time(tm) + dt_3d ) /                  &
                  MAX( TINY( 1.0_wp ), ( slurb_dynamic%time(t)  - slurb_dynamic%time(tm) ) )
       fac_dt = MIN( 1.0_wp, fac_dt )
    ENDIF

!
!-- Update model fields to interpolated values.
    IF ( slurb_dynamic%shf_external%lod == 1 )  THEN
       surf%shf_external(:) = ( 1.0_wp - fac_dt ) * slurb_dynamic%shf_external%var1d(tm) +         &
                              fac_dt * slurb_dynamic%shf_external%var1d(t)
    ELSEIF ( slurb_dynamic%shf_external%lod == 2 ) THEN
       DO  m = 1, surf%ns
          surf%shf_external(m) = ( 1.0_wp - fac_dt ) * slurb_dynamic%shf_external%var2d(tm,m) +    &
                                 fac_dt * slurb_dynamic%shf_external%var2d(t,m)
       ENDDO
    ENDIF

    IF ( moist_physics )  THEN
       IF ( slurb_dynamic%qsws_external%lod == 1 )  THEN
          surf%qsws_external(:) = ( 1.0_wp - fac_dt ) * slurb_dynamic%qsws_external%var1d(tm) +    &
                                  fac_dt * slurb_dynamic%qsws_external%var1d(t)
       ELSEIF ( slurb_dynamic%qsws_external%lod == 2 ) THEN
          DO  m = 1, surf%ns
             surf%qsws_external(m) = ( 1.0_wp - fac_dt ) *                                         &
                                     slurb_dynamic%qsws_external%var2d(tm,m) +                     &
                                     fac_dt * slurb_dynamic%qsws_external%var2d(t,m)
          ENDDO
       ENDIF
    ENDIF

    IF ( slurb_dynamic%shf_traffic%lod == 1 )  THEN
       surf%shf_traffic(:) = ( 1.0_wp - fac_dt ) * slurb_dynamic%shf_traffic%var1d(tm) +           &
                             fac_dt * slurb_dynamic%shf_traffic%var1d(t)
    ELSEIF ( slurb_dynamic%shf_traffic%lod == 2 ) THEN
       DO  m = 1, surf%ns
          surf%shf_traffic(m) = ( 1.0_wp - fac_dt ) * slurb_dynamic%shf_traffic%var2d(tm,m) +      &
                                fac_dt * slurb_dynamic%shf_traffic%var2d(t,m)
       ENDDO
    ENDIF

 END SUBROUTINE update_dynamic_inputs

 END SUBROUTINE slurb_update_external_vars


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> SLUrb's internal model to model urban surface - atmosphere coupling.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_urban_aggregation_model

    INTEGER(iwp) ::  i       !< running index
    INTEGER(iwp) ::  j       !< running index
    INTEGER(iwp) ::  k_topo  !< k index of topography
    INTEGER(iwp) ::  k_atm   !< k index of the first atmospheric level
    INTEGER(iwp) ::  m       !< running index of surface tiles

    LOGICAL ::  runge_l  !< flag for timestep scheme to allow vectorization


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_urban_aggregation_model'
       CALL debug_message( debug_string, 'start' )
    ENDIF

    runge_l = ( timestep_scheme(1:5) == 'runge' )

    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k_topo = topo_top_ind(j,i,0)
       k_atm = topo_top_ind(j,i,0) + 1

!
!--    For shf and qsws, use direct aggregation.
       surf%shf_urb(m) = surf%f_bld(m) * surf%shf_roof(m) +                                        &
                         ( 1.0_wp - surf%f_bld(m) ) * surf%shf_can(m) + surf%shf_external(m)

       IF ( moist_physics )  THEN
          surf%qsws_urb(m) = surf%f_bld(m) * surf%qsws_roof(m) +                                   &
                             ( 1.0_wp - surf%f_bld(m) ) * surf%qsws_can(m) + surf%qsws_external(m)
       ENDIF

!
!--    Calculate momentum flux for horizontal wind components.
       surf%usws_urb(m) = -u(k_atm,j,i) / surf%ram_urb(m) * rho_air_zw(k_topo)
       surf%vsws_urb(m) = -v(k_atm,j,i) / surf%ram_urb(m) * rho_air_zw(k_topo)

!
!--    Aggregate radiative fluxes. Note that this aggregation is done here rather than in the
!--    slurb_radiation_model on purpose to include the longwave term dependent on
!--    the surface temperature of given surface.
!
!--    Compute the net LW radiation flux at canyon top and for urban surface.
       surf%rad_lw_net_can(m) = surf%rad_lw_net_road(m) + surf%hw_can(m) *                         &
                                (   ( 1.0_wp - surf%f_win(m) ) *                                   &
                                    ( surf%rad_lw_net_wall_a(m) + surf%rad_lw_net_wall_b(m) )      &
                                  + surf%f_win(m) *                                                &
                                    ( surf%rad_lw_net_win_a(m)  + surf%rad_lw_net_win_b(m) )       &
                                )

       surf%rad_lw_net_urb(m) = surf%f_bld(m) * surf%rad_lw_net_roof(m) +                          &
                                ( 1.0_wp - surf%f_bld(m) ) * surf%rad_lw_net_can(m)
!
!--    Outgoing LW flux.
       surf%rad_lw_out_urb(m) = surf%rad_lw_in_urb(m) - surf%rad_lw_net_urb(m)
!
!--    Calculate urban aggregated surface temperatures.
       CALL calc_urban_aggregated_temperatures

    ENDDO

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'slurb_urban_aggregation_model'
       CALL debug_message( debug_string, 'end' )
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Update the aggregated urban surface temperatures. These are diagnostic outputs and are not
!> prognostic model variables. Four aggregated urban surface temperatures are computed:
!> 1) Effective surface temperature T_H derived from conservation of heat flux contributions
!> 2) Radiative surface temperature T_rad derived from the outgoing LW radiation
!> 3) Complete surface temperature T_C which is an area-weighted temperature of all facets
!> 4) Theoretical temperature at 2 m height extrapolated using stability-corrected log profile
!> For 1-3 formulations of Kanda et al. 2005, adapted for SLUrb configuration, are used.
!< Note that prognostic temperatures (suffix _p) are used. These are the ones that are output
!> for current time step, as the timelevel is swapped right after the prognostic equation calls.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_urban_aggregated_temperatures

    REAL(wp) ::  c_h_roof    !< bulk heat transfer coefficient for roof
    REAL(wp) ::  c_h_wall_a  !< bulk heat transfer coefficient for wall a
    REAL(wp) ::  c_h_wall_b  !< bulk heat transfer coefficient for wall b
    REAL(wp) ::  c_h_win_a   !< bulk heat transfer coefficient for window a
    REAL(wp) ::  c_h_win_b   !< bulk heat transfer coefficient for window b
    REAL(wp) ::  c_h_road    !< bulk heat transfer coefficient for road
    REAL(wp) ::  ts          !< scaling temperature
    REAL(wp) ::  vtws        !< virtual potential temperature flux (buoyancy flux)

!
!-- 1) Effective surface temperature T_H.
!-- First, compute the bulk heat transfer coefficients.
    IF ( calc_t_h )  THEN
       c_h_roof = ABS( surf%shf_roof(m) / ( rho_cp * surf%uv_eff1(m) *                             &
                       ( surf%t_roof_p(nzt_roof,m) - surf%pt1(m) * exner(k_atm) ) ) )

       c_h_wall_a = ABS( surf%shf_wall_a(m) / ( rho_cp * surf%uv_eff1(m) *                         &
                         ( surf%t_wall_a_p(nzt_wall,m) - surf%pt1(m) * exner(k_atm) ) ) )

       c_h_wall_b = ABS( surf%shf_wall_b(m) / ( rho_cp * surf%uv_eff1(m) *                         &
                         ( surf%t_wall_b_p(nzt_wall,m) - surf%pt1(m) * exner(k_atm) ) ) )

       c_h_win_a = ABS( surf%shf_win_a(m) / ( rho_cp * surf%uv_eff1(m) *                           &
                        ( surf%t_win_a_p(nzt_win,m) - surf%pt1(m) * exner(k_atm) ) ) )

       c_h_win_b = ABS( surf%shf_win_b(m) / ( rho_cp * surf%uv_eff1(m) *                           &
                        ( surf%t_win_b_p(nzt_win,m) - surf%pt1(m) * exner(k_atm) ) ) )

       c_h_road = ABS( surf%shf_road(m) / ( rho_cp * surf%uv_eff1(m) *                             &
                       ( surf%t_road_p(nzt_road,m) - surf%pt1(m) * exner(k_atm) ) ) )

       surf%t_h_urb(m) = ( ( 1.0_wp - surf%f_bld(m) ) *                                            &
                           ( surf%hw_can(m) * (                                                    &
                                                ( 1.0_wp - surf%f_win(m) ) *                       &
                                                ( c_h_wall_a * surf%t_wall_a_p(nzt_wall,m)         &
                                                + c_h_wall_b * surf%t_wall_b_p(nzt_wall,m) )       &
                                              + surf%f_win(m) *                                    &
                                                ( c_h_win_a * surf%t_win_a_p(nzt_win,m)            &
                                                + c_h_win_b * surf%t_win_b_p(nzt_win,m) )          &
                                              )                                                    &
                           + c_h_road * surf%t_road_p(nzt_road,m)                                  &
                           )                                                                       &
                         + surf%f_bld(m) * c_h_roof * surf%t_roof_p(nzt_roof,m)                    &
                         ) /                                                                       &
                         ( ( 1.0_wp - surf%f_bld(m) ) *                                            &
                           ( surf%hw_can(m) * (                                                    &
                                                ( 1.0_wp - surf%f_win(m) ) *                       &
                                                ( c_h_wall_a + c_h_wall_b )                        &
                                              + surf%f_win(m) *                                    &
                                                ( c_h_win_a + c_h_win_b )                          &
                                              )                                                    &
                           + c_h_road                                                              &
                           )                                                                       &
                         + surf%f_bld(m) * c_h_roof + 1E-10_wp                                     &
                         )
    ENDIF

!
!-- 2) Radiative surface temperature T_rad.
    surf%t_rad_urb(m) = SQRT( SQRT( surf%rad_lw_out_urb(m) / ( surf%emiss_urb(m) * sigma_sb ) ) )

!
!-- 3) Complete surface temperature T_C, similarly to T_H but without the C_h weighting.
    IF ( calc_t_c )  THEN
       surf%t_c_urb(m) = ( ( 1.0_wp - surf%f_bld(m) ) *                                            &
                           ( surf%hw_can(m) * ( ( 1.0_wp - surf%f_win(m) ) *                       &
                                     ( surf%t_wall_a_p(nzt_wall,m) + surf%t_wall_b_p(nzt_wall,m) ) &
                                     + surf%f_win(m) *                                             &
                                     ( surf%t_win_a_p(nzt_win,m)   + surf%t_win_b_p(nzt_win,m)   ) &
                                              )                                                    &
                           + surf%t_road_p(nzt_road,m)                                             &
                           )                                                                       &
                           + surf%f_bld(m) * surf%t_roof_p(nzt_roof,m)                             &
                         ) /                                                                       &
                         ( ( 1.0_wp - surf%f_bld(m) ) * ( 2.0_wp * surf%hw_can(m) + 1.0_wp )       &
                           + surf%f_bld(m)                                                         &
                         )
    ENDIF

!
!-- 4) Theoretical 2 m temperature extrapolated using MOST.
    IF ( calc_t_2m )  THEN
       IF ( moist_physics )  THEN
          vtws = surf%shf_can(m) + lv_d_cp * surf%qsws_can(m)
       ELSE
          vtws = surf%shf_can(m)
       ENDIF
       ts = -vtws * drho_air(k_atm) / surf%us_urb(m)

       surf%t_2m_urb(m) = ts / kappa *                                                             &
                          ( LOG( 2.0_wp / ( surf%z_mo(m) + surf%h_bld(m) ) ) -                     &
                            psi_h( 2.0_wp / surf%ol_urb(m) ) +                                     &
                            psi_h( ( surf%z_mo(m) + surf%h_bld(m) ) / surf%ol_urb(m) )             &
                          ) + surf%pt1(m) * exner(k_atm)
    ENDIF

 END SUBROUTINE calc_urban_aggregated_temperatures

 END SUBROUTINE slurb_urban_aggregation_model


!--------------------------------------------------------------------------------------------------!
!   MODLULE PREDEFINED PARAMETERS
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Default parameters for the building types. These are based on the urban surface mod.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_default_pars

!
!-- Residential, < 1950.
    building_pars_slurb(:,1) = (/                                                                  &
       0.18_wp,        &   !< parameter 0   - [-] window fraction
       0.02_wp,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_wp,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.02_wp,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_wp,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.51200E6_wp,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_wp,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.70965E6_wp,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_wp,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_wp,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_wp,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.120_wp,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_wp,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_wp,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_wp,        &   !< parameter 14  - [-] albedo
       0.90_wp,        &   !< parameter 15  - [-] emissivity
       0.02_wp,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.18_wp,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.18_wp,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_wp,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_wp,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       1.5120E6_wp,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.5120E6_wp,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_wp,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_wp,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.810_wp,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.810_wp,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_wp,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_wp,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_wp,        &   !< parameter 29  - [-] albedo
       0.93_wp,        &   !< parameter 30  - [-] emissivity
       0.02_wp,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_wp,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_wp,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_wp,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_wp,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_wp,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_wp,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_wp,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.45_wp,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.45_wp,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.45_wp,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.45_wp,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.70_wp,        &   !< parameter 43  - [-] transmissivity
       0.12_wp,        &   !< parameter 44  - [-] albedo
       0.91_wp         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Residential, 1950 - 2000.
    building_pars_slurb(:,2) = (/                                                                  &
       0.25_wp,        &   !< parameter 0   - [-] window fraction
       0.02_wp,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.15_wp,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.20_wp,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_wp,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.70000E6_wp,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.07920E6_wp,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       2.11200E6_wp,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_wp,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.160_wp,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.046_wp,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       2.100_wp,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_wp,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_wp,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.10_wp,        &   !< parameter 14  - [-] albedo
       0.95_wp,        &   !< parameter 15  - [-] emissivity
       0.02_wp,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.06_wp,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.24_wp,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_wp,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_wp,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_wp,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       2.1120E6_wp,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_wp,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_wp,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.046_wp,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       2.100_wp,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_wp,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_wp,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_wp,        &   !< parameter 29  - [-] albedo
       0.93_wp,        &   !< parameter 30  - [-] emissivity
       0.02_wp,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_wp,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_wp,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_wp,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_wp,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_wp,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_wp,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_wp,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.18_wp,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.18_wp,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.18_wp,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.18_wp,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.65_wp,        &   !< parameter 43  - [-] transmissivity
       0.15_wp,        &   !< parameter 44  - [-] albedo
       0.87_wp         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Residential, > 2000.
    building_pars_slurb(:,3) = (/                                                                  &
       0.29_wp,        &   !< parameter 0   - [-] window fraction
       0.02_wp,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_wp,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.30_wp,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_wp,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       3.75360E6_wp,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_wp,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.07920E6_wp,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_wp,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_wp,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_wp,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.035_wp,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_wp,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_wp,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_wp,        &   !< parameter 14  - [-] albedo
       0.92_wp,        &   !< parameter 15  - [-] emissivity
       0.02_wp,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.20_wp,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.36_wp,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_wp,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_wp,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_wp,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.3400E6_wp,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_wp,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_wp,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.035_wp,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.680_wp,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_wp,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_wp,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.37_wp,        &   !< parameter 29  - [-] albedo
       0.93_wp,        &   !< parameter 30  - [-] emissivity
       0.02_wp,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_wp,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_wp,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_wp,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_wp,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_wp,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_wp,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_wp,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.11_wp,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.11_wp,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.11_wp,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.11_wp,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.57_wp,        &   !< parameter 43  - [-] transmissivity
       0.18_wp,        &   !< parameter 44  - [-] albedo
       0.80_wp         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, < 1950.
    building_pars_slurb(:,4) = (/                                                                  &
       0.18_wp,        &   !< parameter 0   - [-] window fraction
       0.02_wp,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_wp,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.02_wp,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_wp,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.51200E6_wp,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_wp,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.70965E6_wp,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_wp,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_wp,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_wp,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.120_wp,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_wp,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_wp,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_wp,        &   !< parameter 14  - [-] albedo
       0.90_wp,        &   !< parameter 15  - [-] emissivity
       0.02_wp,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.18_wp,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.18_wp,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_wp,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_wp,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       1.5120E6_wp,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.5120E6_wp,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_wp,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_wp,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.810_wp,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.810_wp,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_wp,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_wp,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_wp,        &   !< parameter 29  - [-] albedo
       0.93_wp,        &   !< parameter 30  - [-] emissivity
       0.02_wp,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_wp,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_wp,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_wp,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_wp,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_wp,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_wp,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_wp,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.45_wp,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.45_wp,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.45_wp,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.45_wp,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.70_wp,        &   !< parameter 43  - [-] transmissivity
       0.12_wp,        &   !< parameter 44  - [-] albedo
       0.91_wp         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, 1950 - 2000.
    building_pars_slurb(:,5) = (/                                                                  &
       0.25_wp,        &   !< parameter 0   - [-] window fraction
       0.02_wp,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.15_wp,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.20_wp,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_wp,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.70000E6_wp,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.07920E6_wp,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       2.11200E6_wp,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_wp,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.160_wp,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.046_wp,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       2.100_wp,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_wp,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_wp,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.10_wp,        &   !< parameter 14  - [-] albedo
       0.95_wp,        &   !< parameter 15  - [-] emissivity
       0.02_wp,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.06_wp,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.24_wp,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_wp,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_wp,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_wp,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       2.1120E6_wp,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_wp,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_wp,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.046_wp,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       2.100_wp,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_wp,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_wp,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_wp,        &   !< parameter 29  - [-] albedo
       0.93_wp,        &   !< parameter 30  - [-] emissivity
       0.02_wp,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_wp,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_wp,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_wp,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_wp,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_wp,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_wp,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_wp,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.18_wp,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.18_wp,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.18_wp,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.18_wp,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.65_wp,        &   !< parameter 43  - [-] transmissivity
       0.15_wp,        &   !< parameter 44  - [-] albedo
       0.87_wp         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, > 2000.
    building_pars_slurb(:,6) = (/                                                                  &
       0.29_wp,        &   !< parameter 0   - [-] window fraction
       0.02_wp,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_wp,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.30_wp,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_wp,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       3.75360E6_wp,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_wp,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.07920E6_wp,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_wp,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_wp,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_wp,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.035_wp,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_wp,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_wp,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_wp,        &   !< parameter 14  - [-] albedo
       0.92_wp,        &   !< parameter 15  - [-] emissivity
       0.02_wp,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.20_wp,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.36_wp,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_wp,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_wp,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_wp,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.3400E6_wp,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_wp,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_wp,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.035_wp,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.680_wp,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_wp,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_wp,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.37_wp,        &   !< parameter 29  - [-] albedo
       0.93_wp,        &   !< parameter 30  - [-] emissivity
       0.02_wp,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_wp,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_wp,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_wp,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_wp,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_wp,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_wp,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_wp,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.11_wp,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.11_wp,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.11_wp,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.11_wp,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.57_wp,        &   !< parameter 43  - [-] transmissivity
       0.18_wp,        &   !< parameter 44  - [-] albedo
       0.80_wp         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Asphalt concrete mix (I-II), stone aggregate(III), gravel and soil(IV), PALM-LSM default.
    pavement_pars_slurb(:,1) = (/                                                                  &
       0.01_wp,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_wp,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_wp,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_wp,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.00E6_wp,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.00E6_wp,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_wp,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_wp,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       1.00_wp,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       1.00_wp,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_wp,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_wp,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_wp,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_wp,      &   !< parameter 13  - [-] albedo
       0.93_wp       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Asphalt concrete (I-II), stone aggregate (III), gravel and soil (IV), Masson et al. (2002).
    pavement_pars_slurb(:,2) = (/                                                                  &
       0.01_wp,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_wp,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_wp,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_wp,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       1.74E6_wp,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       1.74E6_wp,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_wp,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_wp,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       0.82_wp,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       0.82_wp,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_wp,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_wp,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_wp,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.10_wp,      &   !< parameter 13  - [-] albedo
       0.95_wp       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Concrete (Portland concrete, I-II), stone aggregate (III), gravel and soil (IV),
!-- Masson et al. (2002) and Yaghoobian et al. (2009).
    pavement_pars_slurb(:,3) = (/                                                                  &
       0.01_wp,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_wp,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_wp,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_wp,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.11E6_wp,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.11E6_wp,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_wp,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_wp,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       1.51_wp,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       1.51_wp,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_wp,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_wp,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_wp,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.30_wp,      &   !< parameter 13  - [-] albedo
       0.90_wp       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Sett (I-II), stone aggregate (III), gravel and soil (IV),Masson et al. (2002), Oke (1987)
!-- and Mandanici et al. (2016).
    pavement_pars_slurb(:,4) = (/                                                                  &
       0.01_wp,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_wp,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_wp,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_wp,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.25E6_wp,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.25E6_wp,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_wp,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_wp,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       2.19_wp,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       2.19_wp,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_wp,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_wp,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_wp,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_wp,      &   !< parameter 13  - [-] albedo
       0.95_wp       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Pavement stones (I-II), stone aggregate (III), gravel and soil (IV),
!-- Masson et al. (2002), Oke (1987) and Göttsche & Hulley (2012).
    pavement_pars_slurb(:,5) = (/                                                                  &
       0.01_wp,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_wp,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_wp,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_wp,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.25E6_wp,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.25E6_wp,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_wp,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_wp,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       2.19_wp,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       2.19_wp,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_wp,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_wp,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_wp,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_wp,      &   !< parameter 13  - [-] albedo
       0.93_wp       &   !< parameter 14  - [-] emissivity
    /)

 END SUBROUTINE slurb_default_pars

END MODULE slurb_mod
