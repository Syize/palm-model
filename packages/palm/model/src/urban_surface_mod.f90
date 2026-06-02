!> @file urban_surface_mod.f90
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
! Copyright 2015-2021 Czech Technical University in Prague
! Copyright 2015-2021 Institute of Computer Science of the Czech Academy of Sciences, Prague
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
! 2016/6/9 - Initial version of the USM (Urban Surface Model)
!            authors: Jaroslav Resler, Pavel Krc (Czech Technical University in Prague and Institute
!            of Computer Science of the Czech Academy of Sciences, Prague)
!            with contributions: Michal Belda, Nina Benesova, Ondrej Vlcek
!            partly inspired by PALM LSM (B. Maronga)
!            parameterizations of Ra checked with TUF3D (E. S. Krayenhoff)
!> Module for Urban Surface Model (USM)
!> The module includes:
!>    1. Radiation model with direct/diffuse radiation, shading, reflections and integration with
!>       plant canopy
!>    2. Wall and wall surface model
!>    3. Surface layer energy balance
!>    4. Anthropogenic heat (only from transportation so far)
!>    5. Necessary auxiliary subroutines (reading inputs, writing outputs, restart simulations, ...)
!> It also makes use of standard radiation and integrates it into urban surface model.
!>
!> Further work:
!> -------------
!> @todo Revise initialization when building_pars / building_surface_pars are provided -
!>       intialization is not consistent to building_pars
!> @todo Revise flux conversion in energy-balance solver
!--------------------------------------------------------------------------------------------------!
module urban_surface_mod

#if defined( __parallel )
   use MPI
#endif

   use arrays_3d, &
      only: exner, &
            hyp, &
            hyrho, &
            p, &
            prr, &
            pt, &
            q, &
            ql, &
            tend, &
            u, &
            v, &
            vpt, &
            w, &
            zu

   use calc_mean_profile_mod, &
      only: calc_mean_profile

   use basic_constants_and_equations_mod, &
      only: c_p, &
            degc_to_k, &
            g, &
            kappa, &
            l_v, &
            magnus_tl, &
            pi, &
            r_d, &
            rho_l, &
            sigma_sb

   use control_parameters, &
      only: allow_roughness_limitation, &
            average_count_3d, &
            coupling_char, &
            coupling_start_time, &
            cut_cell_topography, &
            cyclic_fill_initialization, &
            data_output_raw, &
            debug_output, &
            debug_output_timestep, &
            debug_string, &
            dt_do3d, &
            dt_3d, &
            dz, &
            end_time, &
            humidity, &
            indoor_model, &
            initializing_actions, &
            intermediate_timestep_count, &
            intermediate_timestep_count_max, &
            io_blocks, &
            io_group, &
            large_scale_forcing, &
            lsf_surf, &
            message_string, &
            output_fill_value, &
            pt_surface, &
            read_spinup_data, &
            restart_data_format_output, &
            spinup_phase, &
            surface_pressure, &
            time_since_reference_point, &
            timestep_scheme, &
            topography, &
            tsc, &
            urban_surface, &
            varnamelength

   use bulk_cloud_model_mod, &
      only: bulk_cloud_model, &
            precipitation

   use cpulog, &
      only: cpu_log, &
            log_point, &
            log_point_s

   use grid_variables, &
      only: ddx, &
            ddx2, &
            ddy, &
            ddy2, &
            dx, &
            dy

   use indices, &
      only: nbgp, &
            nnx, &
            nny, &
            nnz, &
            nx, &
            nxl, &
            nxlg, &
            nxr, &
            nxrg, &
            ny, &
            nyn, &
            nyng, &
            nys, &
            nysg, &
            nzb, &
            nzt, &
            topo_top_ind

   use, intrinsic :: iso_c_binding

   use kinds

   use netcdf_data_input_mod, &
      only: albedo_type_f, &
            building_surface_pars_f, &
            building_type_f, &
            input_file_static, &
            netcdf_data_input_parameter_lists, &
            pars, &
            terrain_height_f

   use palm_date_time_mod, &
      only: get_date_time, &
            seconds_per_hour

   use pegrid

   use radiation_model_mod, &
      only: albedo_type, &
            dirname, &
            diridx, &
            dirint, &
            force_radiation_call, &
            id, &
            idown, &
            ieast, &
            inorth, &
            isouth, &
            iup, &
            iwest, &
            nd, &
            nz_urban_b, &
            nz_urban_t, &
            radiation_interaction, &
            radiation, &
            rad_lw_in, &
            rad_lw_out, &
            rad_sw_in, &
            rad_sw_out, &
            unscheduled_radiation_calls

   use restart_data_mpi_io_mod, &
      only: rd_mpi_io_check_array, &
            rd_mpi_io_surface_filetypes, &
            rrd_mpi_io, &
            rrd_mpi_io_surface, &
            wrd_mpi_io, &
            wrd_mpi_io_surface

   use statistics, &
      only: hom, &
            statistic_regions

   use surface_data_handling, &
      only: surface_restore_elements

   use surface_mod, &
      only: ind_pav_green, &
            ind_veg_wall, &
            ind_wat_win, &
            surf_type, &
            surf_usm

   implicit none

!
!-- Indices of input attributes in building_surface_pars (except for radiation-related, which are in
!-- radiation_model_mod)
   character(37), dimension(0:7), parameter ::  building_type_name = (/ &
                                               'user-defined                         ', &  !< type 0
                                               'residential - 1950                   ', &  !< type  1
                                               'residential 1951 - 2000              ', &  !< type  2
                                               'residential 2001 -                   ', &  !< type  3
                                               'office - 1950                        ', &  !< type  4
                                               'office 1951 - 2000                   ', &  !< type  5
                                               'office 2001 -                        ', &  !< type  6
                                               'bridges                              ' &  !< type  7
                                               /)
!
!-- USM model constants
   real(wp), parameter ::  b_ch = 6.04_wp    !< Clapp & Hornberger exponent
   real(wp), parameter ::  lambda_h_green_dry = 0.19_wp    !< heat conductivity for dry soil
   real(wp), parameter ::  lambda_h_green_sm = 3.44_wp    !< heat conductivity of the soil matrix
   real(wp), parameter ::  lambda_h_water = 0.57_wp    !< heat conductivity of water
   real(wp), parameter ::  psi_sat = -0.388_wp  !< soil matrix potential at saturation
   real(wp), parameter ::  rho_c_soil = 2.19e6_wp  !< volumetric heat capacity of soil
   real(wp), parameter ::  rho_c_water = 4.20e6_wp  !< volumetric heat capacity of water
!    REAL(wp), PARAMETER ::  m_max_depth        = 0.0002_wp  !< Maximum capacity of the water reservoir (m)

!
!-- Soil parameters I           alpha_vg,      l_vg_green,    n_vg, gamma_w_green_sat
   real(wp), dimension(0:3, 1:7), parameter ::  soil_pars = reshape((/ &
                                                                    3.83_wp, 1.250_wp, 1.38_wp, 6.94e-6_wp, &  !< soil 1
                                                                    3.14_wp, -2.342_wp, 1.28_wp, 1.16e-6_wp, &  !< soil 2
                                                                    0.83_wp, -0.588_wp, 1.25_wp, 0.26e-6_wp, &  !< soil 3
                                                                    3.67_wp, -1.977_wp, 1.10_wp, 2.87e-6_wp, &  !< soil 4
                                                                    2.65_wp, 2.500_wp, 1.10_wp, 1.74e-6_wp, &  !< soil 5
                                                                    1.30_wp, 0.400_wp, 1.20_wp, 0.93e-6_wp, &  !< soil 6
                                                                    0.00_wp, 0.00_wp, 0.00_wp, 0.57e-6_wp &  !< soil 7
                                                                    /), (/4, 7/))

!
!-- Soil parameters II              swc_sat,     fc,   wilt,    swc_res
   real(wp), dimension(0:3, 1:7), parameter ::  m_soil_pars = reshape((/ &
                                                                      0.403_wp, 0.244_wp, 0.059_wp, 0.025_wp, &  !< soil 1
                                                                      0.439_wp, 0.347_wp, 0.151_wp, 0.010_wp, &  !< soil 2
                                                                      0.430_wp, 0.383_wp, 0.133_wp, 0.010_wp, &  !< soil 3
                                                                      0.520_wp, 0.448_wp, 0.279_wp, 0.010_wp, &  !< soil 4
                                                                      0.614_wp, 0.541_wp, 0.335_wp, 0.010_wp, &  !< soil 5
                                                                      0.766_wp, 0.663_wp, 0.267_wp, 0.010_wp, &  !< soil 6
                                                                      0.472_wp, 0.323_wp, 0.171_wp, 0.000_wp &  !< soil 7
                                                                      /), (/4, 7/))
!
!-- Value 9999999.9_wp -> Generic available or user-defined value must be set otherwise
!-- -> No generic variable and user setting is optional
   real(wp) ::  alpha_vangenuchten = 9999999.9_wp      !< NAMELIST alpha_vg
   real(wp) ::  field_capacity = 9999999.9_wp          !< NAMELIST fc
   real(wp) ::  hydraulic_conductivity = 9999999.9_wp  !< NAMELIST gamma_w_green_sat
   real(wp) ::  l_vangenuchten = 9999999.9_wp          !< NAMELIST l_vg
   real(wp) ::  n_vangenuchten = 9999999.9_wp          !< NAMELIST n_vg
   real(wp) ::  residual_moisture = 9999999.9_wp       !< NAMELIST m_res
   real(wp) ::  saturation_moisture = 9999999.9_wp     !< NAMELIST m_sat
   real(wp) ::  wilting_point = 9999999.9_wp           !< NAMELIST m_wilt

!
!-- Configuration parameters (they can be setup in PALM config)
   logical ::  force_radiation_call_l = .false.   !< flag parameter for unscheduled radiation model calls

   integer(iwp) ::  building_type = 1               !< default building type (preleminary setting)
   integer(iwp) ::  roof_category = 2               !< default category for root surface
   integer(iwp) ::  wall_category = 2               !< default category for wall surface over pedestrian zone

   integer(iwp), parameter ::  nzb_wall = 0  !< inner side of the wall model (to be switched)
   integer(iwp), parameter ::  nzt_wall = 3  !< outer side of the wall model (to be switched)
   integer(iwp), parameter ::  nzw = 4  !< number of wall layers (fixed for now)

   integer(iwp)            ::  soil_type     !<

!
!-- Indicies of input attributes for urban surface types.
   integer(iwp) ::  ind_green_agfl = 7    !< index in input list for green on wall, above ground floor level
   integer(iwp) ::  ind_green_gfl = 6    !< index in input list for green on wall, ground floor level
   integer(iwp) ::  ind_green_roof = 8    !< index in input list for green on roof
   integer(iwp) ::  ind_wall_agfl = 1    !< index in input list for wall, above ground floor level
   integer(iwp) ::  ind_wall_gfl = 0    !< index in input list for wall, ground floor level
   integer(iwp) ::  ind_wall_roof = 2    !< index in input list for wall, roof
   integer(iwp) ::  ind_win_agfl = 4    !< index in input list for window, above ground floor level
   integer(iwp) ::  ind_win_gfl = 3    !< index in input list for window, ground floor level
   integer(iwp) ::  ind_win_roof = 5    !< index in input list for window, roof
!
!-- Indices of input attributes for urban surface level.
   integer(iwp) ::  ind_agfl = 1    !< index in input list above ground floor level
   integer(iwp) ::  ind_gfl = 0    !< index in input list ground floor level
   integer(iwp) ::  ind_roof = 2    !< index in input list roof
!
!-- Indicies of input attributes in building_gen.
   integer(iwp) ::  ind_gflh = 0    !< index in input list for ground floor level height
   integer(iwp) ::  ind_green_type_roof = 1    !< index in input list for type of green roof
!
!-- Indicies of input attributes in building_indoor.
   integer(iwp) ::  ind_theta_int_c_set = 0   !< index in input list for indoor target summer temperature
   integer(iwp) ::  ind_theta_int_h_set = 1   !< index in input list for indoor target winter temperature
   integer(iwp) ::  ind_f_c_win = 2   !< index in input list for shading factor
   integer(iwp) ::  ind_g_value_win = 3   !< index in input list for g-value windows
   integer(iwp) ::  ind_u_value_win = 4   !< index in input list for u-value windows
   integer(iwp) ::  ind_airflow_unocc = 5   !< index in input list for basic airflow without occupancy of the room
   integer(iwp) ::  ind_airflow_occ = 6   !< index in input list for additional airflow dependent on occupancy of the room
   integer(iwp) ::  ind_eta_ve = 7   !< index in input list for heat recovery efficiency
   integer(iwp) ::  ind_factor_a = 8   !< index in input list for dynamic parameter specific effective surface
   integer(iwp) ::  ind_factor_c = 9   !< index in input list for dynamic parameter innner heat storage
   integer(iwp) ::  ind_lambda_at = 10  !< index in input list for ratio internal surface/floor area
   integer(iwp) ::  ind_q_h_max = 11  !< index in input list for maximal heating capacity
   integer(iwp) ::  ind_q_c_max = 12  !< index in input list for maximal cooling capacity
   integer(iwp) ::  ind_qint_high = 13
   !< index in input list for additional internal heat gains dependent on occupancy of the room
   integer(iwp) ::  ind_qint_low = 14  !< index in input list for basic internal heat gains without occupancy of the room
   integer(iwp) ::  ind_height_storey = 15  !< index in input list for storey height
   integer(iwp) ::  ind_height_cei_con = 16  !< index in input list for ceiling construction height
   integer(iwp) ::  ind_params_waste_heat_h = 17  !< index in input list for anthropogenic heat output factor for heating
   integer(iwp) ::  ind_params_waste_heat_c = 18  !< index in input list for anthropogenic heat output factor for cooling

   integer(iwp) ::  ind_s_emis_green = 14  !< index for emissivity of green fraction (0-1)
   integer(iwp) ::  ind_s_emis_wall = 13  !< index for emissivity of wall fraction (0-1)
   integer(iwp) ::  ind_s_emis_win = 15  !< index for emissivity o f window fraction (0-1)
   integer(iwp) ::  ind_s_green_frac_r = 3   !< index for green fraction on roof (0-1)
   integer(iwp) ::  ind_s_green_frac_w = 2   !< index for green fraction on wall (0-1)
   integer(iwp) ::  ind_s_hc1 = 5   !< index for heat capacity of wall layer 1
   integer(iwp) ::  ind_s_hc2 = 6   !< index for heat capacity of wall layer 2
   integer(iwp) ::  ind_s_hc3 = 7   !< index for heat capacity of wall layer 3
   integer(iwp) ::  ind_s_indoor_target_temp_summer = 11  !< index for indoor target summer temperature
   integer(iwp) ::  ind_s_indoor_target_temp_winter = 12  !< index for indoor target winter temperature
   integer(iwp) ::  ind_s_lai_r = 4   !< index for leaf area index of green fraction
   integer(iwp) ::  ind_s_tc1 = 8   !< index for thermal conducivity of wall layer 1
   integer(iwp) ::  ind_s_tc2 = 9   !< index for thermal conducivity of wall layer 2
   integer(iwp) ::  ind_s_tc3 = 10  !< index for thermal conducivity of wall layer 3
   integer(iwp) ::  ind_s_trans = 16  !< index for transmissivity of window fraction (0-1)
   integer(iwp) ::  ind_s_wall_frac = 0   !< index for wall fraction (0-1)
   integer(iwp) ::  ind_s_win_frac = 1   !< index for window fraction (0-1)
   integer(iwp) ::  ind_s_z0 = 17  !< index for roughness length for momentum (m)
   integer(iwp) ::  ind_s_z0qh = 18  !< index for roughness length for heat (m)

   real(wp) ::  d_roughness_concrete                 !< inverse roughness length of average concrete surface
   real(wp) ::  dt_usm = huge(1.0_wp)              !< maximum allowed timestep of the urban-surface model
   real(wp) ::  ground_floor_level = 4.0_wp          !< default ground floor level
   real(wp) ::  m_total = 0.0_wp    !< weighted total water content of the soil (m3/m3)
   real(wp) ::  roof_inner_temperature = 295.0_wp  !< temperature of the inner roof surface (~22 degrees C) (K)
   real(wp) ::  roughness_concrete = 0.001_wp        !< roughness length of average concrete surface
   real(wp) ::  wall_inner_temperature = 295.0_wp  !< temperature of the inner wall surface (~22 degrees C) (K)
   real(wp) ::  window_inner_temperature = 295.0_wp  !< temperature of the inner window surface (~22 degrees C) (K)

!
!-- Building properties. Different kind of properties are listed in different tables.
!-- Initialization of building parameters is outsourced to usm_init_pars. This is needed because of the
!-- huge number of attributes given in building_pars (>700), while intel and gfortran compiler have
!-- hard limit of continuation lines of 511.
   integer(iwp), dimension(0:8, 1:7) ::  building_alb_type   !< albedo type of building surfaces

   real(wp), dimension(0:8, 1:7)   ::  building_emis    !< emissivity of building surfaces
   real(wp), dimension(0:8, 1:7)   ::  building_frac    !< wall, window and green fractions
   real(wp), dimension(0:1, 1:7)   ::  building_gen     !< general building parameters
   real(wp), dimension(0:18, 1:7)  ::  building_indoor  !< indoor building parameters
   real(wp), dimension(0:2, 1:7)   ::  building_lai     !< leaf area index of building surfaces
   real(wp), dimension(0:149, 1:7) ::  building_pars    !< list containing all building parameters (deprecated)
   real(wp), dimension(0:2, 1:7)   ::  building_trans   !< window transmissivities
   real(wp), dimension(0:2, 1:7)   ::  building_z0      !< roughness length of building surfaces
   real(wp), dimension(0:2, 1:7)   ::  building_z0qh    !< roughness length for moisture and heat of building surfaces

   real(wp), dimension(0:8, 0:nzw - 1, 1:7) ::  building_hcap   !< heat capacities at different layers
   real(wp), dimension(0:8, 0:nzw - 1, 1:7) ::  building_hcond  !< heat conductivities at different layers
   real(wp), dimension(0:8, 0:nzw - 1, 1:7) ::  building_depth  !< layer depths

!
!-- Define the corresponding NetCDF input variables.
   type(pars) ::  building_pars_f      !< input variable for building parameters (deprecated)
   type(pars) ::  building_alb_type_f  !< input variable for building surface albedo type
   type(pars) ::  building_emis_f      !< input variable for building surface emissivity
   type(pars) ::  building_frac_f      !< input variable for building surface fractions
   type(pars) ::  building_gen_f       !< input variable for general building parameters
   type(pars) ::  building_hcap_f      !< input variable for building surface layer heat capacities
   type(pars) ::  building_hcond_f     !< input variable for building surface layer conductivities
   type(pars) ::  building_indoor_f    !< input variable for building indoor parameters
   type(pars) ::  building_lai_f       !< input variable for building surface leaf area index
   type(pars) ::  building_thick_f     !< input variable for building surface layer thicknesses
   type(pars) ::  building_trans_f     !< input variable for building window transmissivities
   type(pars) ::  building_z0_f        !< input variable for building roughness length
   type(pars) ::  building_z0qh_f      !< input variable for building roughness length for moisture and heat
!
!-- Type for 1d surface variables as surface temperature and liquid water reservoir
   type surf_type_1d_usm
      real(wp), dimension(:), allocatable ::  val  !<
   end type surf_type_1d_usm
!
!-- Type for 2d surface variables as wall temperature
   type surf_type_2d_usm
      real(wp), dimension(:, :), allocatable ::  val  !<
   end type surf_type_2d_usm
!
!-- Surface and material model variables for walls, ground, roofs
   type(surf_type_1d_usm), pointer ::  t_surf_green      !< prognostic array for green surface temperature
   type(surf_type_1d_usm), pointer ::  t_surf_green_p    !< prognostic array for green surface temperature
   type(surf_type_1d_usm), pointer ::  t_surf_wall       !< prognostic array for wall surface temperature
   type(surf_type_1d_usm), pointer ::  t_surf_wall_p     !< prognostic array for wall surface temperature
   type(surf_type_1d_usm), pointer ::  t_surf_window     !< prognostic array for window surface temperature
   type(surf_type_1d_usm), pointer ::  t_surf_window_p   !< prognostic array for window surface temperature

   type(surf_type_1d_usm), target ::  t_surf_green_1    !<
   type(surf_type_1d_usm), target ::  t_surf_green_2    !<
   type(surf_type_1d_usm), target ::  t_surf_wall_1     !<
   type(surf_type_1d_usm), target ::  t_surf_wall_2     !<
   type(surf_type_1d_usm), target ::  t_surf_window_1   !<
   type(surf_type_1d_usm), target ::  t_surf_window_2   !<

!
!-- Energy balance variables
!-- Parameters of the land, roof and wall surfaces (only for horizontally upward surfaces)
   type(surf_type_1d_usm), pointer ::  m_liq_usm   !< liquid water reservoir (m), horizontal surface elements
   type(surf_type_1d_usm), pointer ::  m_liq_usm_p !< progn. liquid water reservoir (m), horizontal surface elements

   type(surf_type_1d_usm), target ::  m_liq_usm_1  !<
   type(surf_type_1d_usm), target ::  m_liq_usm_2  !<
   type(surf_type_1d_usm), target ::  tm_liq_usm_m !< liquid water reservoir tendency (m), horizontal surface elements

   type(surf_type_2d_usm), pointer ::  fc          !<
   type(surf_type_2d_usm), pointer ::  rootfr      !<
   type(surf_type_2d_usm), pointer ::  swc         !<
   type(surf_type_2d_usm), pointer ::  swc_p       !<
   type(surf_type_2d_usm), pointer ::  swc_res     !<
   type(surf_type_2d_usm), pointer ::  swc_sat     !<
   type(surf_type_2d_usm), pointer ::  t_green     !<
   type(surf_type_2d_usm), pointer ::  t_green_p   !<
   type(surf_type_2d_usm), pointer ::  t_wall      !<
   type(surf_type_2d_usm), pointer ::  t_wall_p    !<
   type(surf_type_2d_usm), pointer ::  wilt        !<
   type(surf_type_2d_usm), pointer ::  t_window    !<
   type(surf_type_2d_usm), pointer ::  t_window_p  !<

   type(surf_type_2d_usm), target ::  fc_1        !<
   type(surf_type_2d_usm), target ::  rootfr_1    !<
   type(surf_type_2d_usm), target ::  swc_1       !<
   type(surf_type_2d_usm), target ::  swc_2       !<
   type(surf_type_2d_usm), target ::  swc_res_1   !<
   type(surf_type_2d_usm), target ::  swc_sat_1   !<
   type(surf_type_2d_usm), target ::  t_green_1   !<
   type(surf_type_2d_usm), target ::  t_green_2   !<
   type(surf_type_2d_usm), target ::  t_wall_1    !<
   type(surf_type_2d_usm), target ::  t_wall_2    !<
   type(surf_type_2d_usm), target ::  wilt_1      !<
   type(surf_type_2d_usm), target ::  t_window_1  !<
   type(surf_type_2d_usm), target ::  t_window_2  !<

!
!-- Arrays for time averages
   type(surf_type_1d_usm) ::  wghf_eb_av          !< average of wghf_eb
   type(surf_type_1d_usm) ::  wghf_eb_window_av   !< average of wghf_eb window
   type(surf_type_1d_usm) ::  wghf_eb_green_av    !< average of wghf_eb window
   type(surf_type_1d_usm) ::  iwghf_eb_av         !< indoor average of wghf_eb
   type(surf_type_1d_usm) ::  iwghf_eb_window_av  !< indoor average of wghf_eb window
   type(surf_type_1d_usm) ::  wshf_eb_av          !< average of wshf_eb
   type(surf_type_1d_usm) ::  qsws_av             !< average of qsws
   type(surf_type_1d_usm) ::  qsws_veg_av         !< average of qsws_veg_eb
   type(surf_type_1d_usm) ::  qsws_liq_av         !< average of qsws_liq_eb
   type(surf_type_1d_usm) ::  t_surf_wall_av      !< average of wall surface temperature (K)
   type(surf_type_1d_usm) ::  t_surf_window_av    !< average of window surface temperature (K)
   type(surf_type_1d_usm) ::  t_surf_green_av     !< average of green wall surface temperature (K)

   type(surf_type_2d_usm) ::  t_wall_av           !< average of t_wall
   type(surf_type_2d_usm) ::  t_window_av         !< average of t_window
   type(surf_type_2d_usm) ::  t_green_av          !< average of t_green
   type(surf_type_2d_usm) ::  swc_av              !< average of swc

!
!-- Interfaces of subroutines accessed from outside of this module
   interface usm_3d_data_averaging
      module procedure usm_3d_data_averaging
   end interface usm_3d_data_averaging

   interface usm_boundary_condition
      module procedure usm_boundary_condition
   end interface usm_boundary_condition

   interface usm_check_data_output
      module procedure usm_check_data_output
   end interface usm_check_data_output

   interface usm_check_parameters
      module procedure usm_check_parameters
   end interface usm_check_parameters

   interface usm_data_output_3d
      module procedure usm_data_output_3d
   end interface usm_data_output_3d

   interface usm_define_netcdf_grid
      module procedure usm_define_netcdf_grid
   end interface usm_define_netcdf_grid

   interface usm_init
      module procedure usm_init
   end interface usm_init

   interface usm_init_arrays
      module procedure usm_init_arrays
   end interface usm_init_arrays

   interface usm_parin
      module procedure usm_parin
   end interface usm_parin

   interface usm_rrd_local
      module procedure usm_rrd_local_ftn
      module procedure usm_rrd_local_mpi
   end interface usm_rrd_local

   interface usm_energy_balance
      module procedure usm_energy_balance
   end interface usm_energy_balance

   interface usm_swap_timelevel
      module procedure usm_swap_timelevel
   end interface usm_swap_timelevel

   interface usm_timestep
      module procedure usm_timestep
   end interface usm_timestep

   interface usm_vm_sampling
      module procedure usm_vm_sampling
   end interface usm_vm_sampling

   interface usm_wrd_local
      module procedure usm_wrd_local
   end interface usm_wrd_local

   save

   private

!
!-- Public functions
   public usm_boundary_condition, &
      usm_check_data_output, &
      usm_check_parameters, &
      usm_data_output_3d, &
      usm_define_netcdf_grid, &
      usm_init, &
      usm_init_arrays, &
      usm_parin, &
      usm_rrd_local, &
      usm_energy_balance, &
      usm_swap_timelevel, &
      usm_timestep, &
      usm_vm_sampling, &
      usm_wrd_local, &
      usm_3d_data_averaging

!
!-- Public parameters, constants and initial values
   public building_type, &
      building_hcond, &
      building_gen, &
      building_indoor, &
      building_depth, &
      dt_usm, &
      ind_wall_gfl, &
      ind_theta_int_c_set, &
      ind_theta_int_h_set, &
      ind_f_c_win, &
      ind_g_value_win, &
      ind_u_value_win, &
      ind_airflow_unocc, &
      ind_airflow_occ, &
      ind_eta_ve, &
      ind_factor_a, &
      ind_factor_c, &
      ind_lambda_at, &
      ind_q_h_max, &
      ind_q_c_max, &
      ind_qint_high, &
      ind_qint_low, &
      ind_height_storey, &
      ind_height_cei_con, &
      ind_params_waste_heat_h, &
      ind_params_waste_heat_c, &
      nzb_wall, &
      nzt_wall, &
      t_green, &
      t_wall, &
      t_window

contains

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates the necessary indices of the urban surfaces and plant canopy and it
!> allocates the needed arrays for USM
!--------------------------------------------------------------------------------------------------!
   subroutine usm_init_arrays

      implicit none

      if (debug_output) call debug_message('usm_init_arrays', 'start')

!
!-- Allocate radiation arrays which are part of the new data type.
!-- For horizontal surfaces.
      allocate (surf_usm%surfhf(1:surf_usm%ns))
      allocate (surf_usm%rad_net_l(1:surf_usm%ns))

!
!-- Wall surface model
!-- Allocate arrays for wall surface model and define pointers
!-- Allocate array of wall types and wall parameters
      allocate (surf_usm%surface_types(1:surf_usm%ns))
      allocate (surf_usm%building_type(1:surf_usm%ns))
      allocate (surf_usm%building_type_name(1:surf_usm%ns))
      surf_usm%building_type = 0
      surf_usm%building_type_name = 'none'
!
!-- Allocate albedo_type and albedo. Each surface element has 3 values, 0: wall fraction,
!-- 1: green fraction, 2: window fraction.
      allocate (surf_usm%albedo_type(1:surf_usm%ns, 0:2))
      allocate (surf_usm%albedo(1:surf_usm%ns, 0:2))
      surf_usm%albedo_type = albedo_type
!
!-- Allocate indoor target temperature for summer and winter
      allocate (surf_usm%target_temp_summer(1:surf_usm%ns))
      allocate (surf_usm%target_temp_winter(1:surf_usm%ns))
!
!-- Allocate flag indicating ground floor level surface elements
      allocate (surf_usm%gfl(1:surf_usm%ns))
!
!-- Allocate arrays for relative surface fraction.
!-- 0 - wall fraction, 1 - green fraction, 2 - window fraction
      allocate (surf_usm%frac(1:surf_usm%ns, 0:2))
      surf_usm%frac = 0.0_wp
!
!-- Wall and roof surface parameters.
      allocate (surf_usm%isroof_surf(1:surf_usm%ns))
      allocate (surf_usm%lambda_surf(1:surf_usm%ns))
      allocate (surf_usm%lambda_surf_window(1:surf_usm%ns))
      allocate (surf_usm%lambda_surf_green(1:surf_usm%ns))
      allocate (surf_usm%c_surface(1:surf_usm%ns))
      allocate (surf_usm%c_surface_window(1:surf_usm%ns))
      allocate (surf_usm%c_surface_green(1:surf_usm%ns))
      allocate (surf_usm%transmissivity(1:surf_usm%ns))
      allocate (surf_usm%lai(1:surf_usm%ns))
      allocate (surf_usm%emissivity(1:surf_usm%ns, 0:2))
      allocate (surf_usm%r_a_green(1:surf_usm%ns))
      allocate (surf_usm%r_a_window(1:surf_usm%ns))
      allocate (surf_usm%green_type_roof(1:surf_usm%ns))
      allocate (surf_usm%r_s(1:surf_usm%ns))
!
!-- Allocate wall and roof material parameters.
      allocate (surf_usm%thickness_wall(1:surf_usm%ns))
      allocate (surf_usm%thickness_window(1:surf_usm%ns))
      allocate (surf_usm%thickness_green(1:surf_usm%ns))
      allocate (surf_usm%lambda_h(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%lambda_h_layer(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%rho_c_wall(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%lambda_h_window(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%lambda_h_window_layer(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%rho_c_window(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%lambda_h_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%rho_c_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%rho_c_total_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%n_vg_green(1:surf_usm%ns))
      allocate (surf_usm%alpha_vg_green(1:surf_usm%ns))
      allocate (surf_usm%l_vg_green(1:surf_usm%ns))
      allocate (surf_usm%gamma_w_green_sat(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      allocate (surf_usm%lambda_w_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%lambda_w_green_layer(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%gamma_w_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%gamma_w_green_layer(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%tswc_m(nzb_wall:nzt_wall, 1:surf_usm%ns))

!
!-- Allocate green wall and roof vegetation and soil parameters.
      allocate (surf_usm%g_d(1:surf_usm%ns))
      allocate (surf_usm%c_liq(1:surf_usm%ns))
      allocate (surf_usm%qsws_liq(1:surf_usm%ns))
      allocate (surf_usm%qsws_veg(1:surf_usm%ns))
      allocate (surf_usm%r_canopy(1:surf_usm%ns))
      allocate (surf_usm%r_canopy_min(1:surf_usm%ns))
!
!-- Allocate wall and roof layers sizes.
      allocate (surf_usm%dz_wall(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%dz_window(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%dz_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%ddz_wall(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%dz_wall_center(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%ddz_wall_center(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%zw(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%ddz_window(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%dz_window_center(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%ddz_window_center(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%zw_window(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%ddz_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%dz_green_center(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%ddz_green_center(nzb_wall:nzt_wall, 1:surf_usm%ns))
      allocate (surf_usm%zw_green(nzb_wall:nzt_wall, 1:surf_usm%ns))
!
!-- Allocate wall and roof temperature arrays.
!-- Allocate if required. Note, in case of restarts, some of these arrays might be already allocated.
      if (.not. allocated(t_surf_wall_1%val)) allocate (t_surf_wall_1%val(1:surf_usm%ns))
      if (.not. allocated(t_surf_wall_2%val)) allocate (t_surf_wall_2%val(1:surf_usm%ns))
      if (.not. allocated(t_wall_1%val)) allocate (t_wall_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(t_wall_2%val)) allocate (t_wall_2%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(t_surf_window_1%val)) allocate (t_surf_window_1%val(1:surf_usm%ns))
      if (.not. allocated(t_surf_window_2%val)) allocate (t_surf_window_2%val(1:surf_usm%ns))
      if (.not. allocated(t_window_1%val)) allocate (t_window_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(t_window_2%val)) allocate (t_window_2%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(t_surf_green_1%val)) allocate (t_surf_green_1%val(1:surf_usm%ns))
      if (.not. allocated(t_surf_green_2%val)) allocate (t_surf_green_2%val(1:surf_usm%ns))
      if (.not. allocated(t_green_1%val)) allocate (t_green_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(t_green_2%val)) allocate (t_green_2%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(swc_1%val)) allocate (swc_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(swc_sat_1%val)) allocate (swc_sat_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(swc_res_1%val)) allocate (swc_res_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(swc_2%val)) allocate (swc_2%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(rootfr_1%val)) allocate (rootfr_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(wilt_1%val)) allocate (wilt_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(fc_1%val)) allocate (fc_1%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      if (.not. allocated(m_liq_usm_1%val)) allocate (m_liq_usm_1%val(1:surf_usm%ns))
      if (.not. allocated(m_liq_usm_2%val)) allocate (m_liq_usm_2%val(1:surf_usm%ns))

!
!-- Initial assignment of the pointers
      t_wall => t_wall_1; t_wall_p => t_wall_2
      t_window => t_window_1; t_window_p => t_window_2
      t_green => t_green_1; t_green_p => t_green_2
      t_surf_wall => t_surf_wall_1; t_surf_wall_p => t_surf_wall_2
      t_surf_window => t_surf_window_1; t_surf_window_p => t_surf_window_2
      t_surf_green => t_surf_green_1; t_surf_green_p => t_surf_green_2
      m_liq_usm => m_liq_usm_1; m_liq_usm_p => m_liq_usm_2
      swc => swc_1; swc_p => swc_2
      swc_sat => swc_sat_1
      swc_res => swc_res_1
      rootfr => rootfr_1
      wilt => wilt_1
      fc => fc_1

!
!-- Allocate intermediate timestep arrays. For horizontal surfaces.
      allocate (surf_usm%tt_surface_wall_m(1:surf_usm%ns))
      allocate (surf_usm%tt_wall_m(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      allocate (surf_usm%tt_surface_window_m(1:surf_usm%ns))
      allocate (surf_usm%tt_window_m(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      allocate (surf_usm%tt_green_m(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
      allocate (surf_usm%tt_surface_green_m(1:surf_usm%ns))
!
!-- Allocate intermediate timestep arrays
      allocate (tm_liq_usm_m%val(1:surf_usm%ns))
      tm_liq_usm_m%val = 0.0_wp
!
!-- Set inital values for prognostic quantities
      if (allocated(surf_usm%tt_surface_wall_m)) surf_usm%tt_surface_wall_m = 0.0_wp
      if (allocated(surf_usm%tt_wall_m)) surf_usm%tt_wall_m = 0.0_wp
      if (allocated(surf_usm%tt_surface_window_m)) surf_usm%tt_surface_window_m = 0.0_wp
      if (allocated(surf_usm%tt_window_m)) surf_usm%tt_window_m = 0.0_wp
      if (allocated(surf_usm%tt_green_m)) surf_usm%tt_green_m = 0.0_wp
      if (allocated(surf_usm%tt_surface_green_m)) surf_usm%tt_surface_green_m = 0.0_wp
!
!-- Allocate wall heat flux output arrays and set initial values.
      allocate (surf_usm%ghf(1:surf_usm%ns))
      allocate (surf_usm%wshf_eb(1:surf_usm%ns))
      allocate (surf_usm%wghf_eb(1:surf_usm%ns))
      allocate (surf_usm%wghf_eb_window(1:surf_usm%ns))
      allocate (surf_usm%wghf_eb_green(1:surf_usm%ns))
      allocate (surf_usm%iwghf_eb(1:surf_usm%ns))
      allocate (surf_usm%iwghf_eb_window(1:surf_usm%ns))
      if (allocated(surf_usm%ghf)) surf_usm%ghf = 0.0_wp
      if (allocated(surf_usm%wshf_eb)) surf_usm%wshf_eb = 0.0_wp
      if (allocated(surf_usm%wghf_eb)) surf_usm%wghf_eb = 0.0_wp
      if (allocated(surf_usm%wghf_eb_window)) surf_usm%wghf_eb_window = 0.0_wp
      if (allocated(surf_usm%wghf_eb_green)) surf_usm%wghf_eb_green = 0.0_wp
      if (allocated(surf_usm%iwghf_eb)) surf_usm%iwghf_eb = 0.0_wp
      if (allocated(surf_usm%iwghf_eb_window)) surf_usm%iwghf_eb_window = 0.0_wp
!
!-- Initialize building-surface properties, which are also required by other modules, e.g. the
!-- indoor model.
      call usm_define_pars

      if (debug_output) call debug_message('usm_init_arrays', 'end')

   end subroutine usm_init_arrays

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average urban surface output quantities as well as allocate the array necessary
!> for storing the average.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_3d_data_averaging(mode, variable)

      implicit none

      character(LEN=*), intent(IN) ::  variable  !<
      character(LEN=*), intent(IN) ::  mode      !<
      character(LEN=varnamelength) ::  var       !< trimmed variable

      integer(iwp) ::  i, j, k, m, ids, idsint, iwl, istat  !< running indices

      logical ::  downward  !< control flag indicating the output of downward-facing surfaces
      logical ::  eastward  !< control flag indicating the output of east-facing surfaces
      logical ::  northward !< control flag indicating the output of northward-facing surfaces
      logical ::  southward !< control flag indicating the output of southward-facing surfaces
      logical ::  upward    !< control flag indicating the output of upward-facing surfaces
      logical ::  westward  !< control flag indicating the output of westward-facing surfaces

      if (.not. variable(1:4) == 'usm_') return  ! Is such a check really required?

!
!-- Find the real name of the variable
      ids = -1
      var = trim(variable)
      do i = 0, nd - 1
         k = len(trim(var))
         j = len(trim(dirname(i)))
         if (trim(var(k - j + 1:k)) == trim(dirname(i))) then
            ids = i
            idsint = dirint(ids)
            var = var(:k - j)
            exit
         end if
      end do
      if (ids == -1) then
         var = trim(variable)
      else
!
!--    Set direction control flags
         downward = .false.
         eastward = .false.
         northward = .false.
         southward = .false.
         upward = .false.
         westward = .false.
         if (idsint == iup) then
            upward = .true.
         elseif (idsint == idown) then
            downward = .true.
         elseif (idsint == ieast) then
            eastward = .true.
         elseif (idsint == iwest) then
            westward = .true.
         elseif (idsint == inorth) then
            northward = .true.
         elseif (idsint == isouth) then
            southward = .true.
         end if
      end if
      if (var(1:11) == 'usm_t_wall_' .and. len(trim(var)) >= 12) then
!
!--    Wall layers
         read (var(12:12), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:10)
         else
!
!--       Wrong wall layer index
            return
         end if
      end if
      if (var(1:13) == 'usm_t_window_' .and. len(trim(var)) >= 14) then
!
!--      Wall layers
         read (var(14:14), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:12)
         else
!
!--         Wrong window layer index
            return
         end if
      end if
      if (var(1:12) == 'usm_t_green_' .and. len(trim(var)) >= 13) then
!
!--      Wall layers
         read (var(13:13), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:11)
         else
!
!--         Wrong green layer index
            return
         end if
      end if
      if (var(1:8) == 'usm_swc_' .and. len(trim(var)) >= 9) then
!
!--      Swc layers
         read (var(9:9), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:7)
         else
!
!--         Wrong swc layer index
            return
         end if
      end if

      if (mode == 'allocate') then

         select case (trim(var))

         case ('usm_wshf')
!
!--          Sensible heat flux
            if (.not. allocated(wshf_eb_av%val)) then
               allocate (wshf_eb_av%val(1:surf_usm%ns))
               wshf_eb_av%val = 0.0_wp
            end if

         case ('usm_qsws')
!
!--          Latent heat flux
            if (.not. allocated(qsws_av%val)) then
               allocate (qsws_av%val(1:surf_usm%ns))
               qsws_av%val = 0.0_wp
            end if

         case ('usm_qsws_veg')
!
!--          Latent heat flux from vegetation surfaces
            if (.not. allocated(qsws_veg_av%val)) then
               allocate (qsws_veg_av%val(1:surf_usm%ns))
               qsws_veg_av%val = 0.0_wp
            end if

         case ('usm_qsws_liq')
!
!--          Latent heat flux from surfaces with liquid
            if (.not. allocated(qsws_liq_av%val)) then
               allocate (qsws_liq_av%val(1:surf_usm%ns))
               qsws_liq_av%val = 0.0_wp
            end if
!
!--       Please note, the following output quantities belongs to the individual tile fractions -
!--       ground heat flux at wall-, window-, and green fraction. Aggregated ground-heat flux is
!--       treated accordingly in average_3d_data, sum_up_3d_data, etc..
         case ('usm_wghf')
!
!--          Heat flux from ground (wall)
            if (.not. allocated(wghf_eb_av%val)) then
               allocate (wghf_eb_av%val(1:surf_usm%ns))
               wghf_eb_av%val = 0.0_wp
            end if

         case ('usm_wghf_window')
!
!--          Heat flux from window ground
            if (.not. allocated(wghf_eb_window_av%val)) then
               allocate (wghf_eb_window_av%val(1:surf_usm%ns))
               wghf_eb_window_av%val = 0.0_wp
            end if

         case ('usm_wghf_green')
!
!--          Heat flux from green ground
            if (.not. allocated(wghf_eb_green_av%val)) then
               allocate (wghf_eb_green_av%val(1:surf_usm%ns))
               wghf_eb_green_av%val = 0.0_wp
            end if

         case ('usm_iwghf')
!
!--          Heat flux from indoor ground
            if (.not. allocated(iwghf_eb_av%val)) then
               allocate (iwghf_eb_av%val(1:surf_usm%ns))
               iwghf_eb_av%val = 0.0_wp
            end if

         case ('usm_iwghf_window')
!
!--          Heat flux from indoor window ground
            if (.not. allocated(iwghf_eb_window_av%val)) then
               allocate (iwghf_eb_window_av%val(1:surf_usm%ns))
               iwghf_eb_window_av%val = 0.0_wp
            end if

         case ('usm_t_surf_wall')
!
!--          Surface temperature for wall surfaces
            if (.not. allocated(t_surf_wall_av%val)) then
               allocate (t_surf_wall_av%val(1:surf_usm%ns))
               t_surf_wall_av%val = 0.0_wp
            end if

         case ('usm_t_surf_window')
!
!--          Surface temperature for window surfaces
            if (.not. allocated(t_surf_window_av%val)) then
               allocate (t_surf_window_av%val(1:surf_usm%ns))
               t_surf_window_av%val = 0.0_wp
            end if

         case ('usm_t_surf_green')
!
!--          Surface temperature for green surfaces
            if (.not. allocated(t_surf_green_av%val)) then
               allocate (t_surf_green_av%val(1:surf_usm%ns))
               t_surf_green_av%val = 0.0_wp
            end if

         case ('usm_t_wall')
!
!--          Wall temperature for iwl layer
            if (.not. allocated(t_wall_av%val)) then
               allocate (t_wall_av%val(nzb_wall:nzt_wall, 1:surf_usm%ns))
               t_wall_av%val = 0.0_wp
            end if

         case ('usm_t_window')
!
!--          Window temperature for iwl layer
            if (.not. allocated(t_window_av%val)) then
               allocate (t_window_av%val(nzb_wall:nzt_wall, 1:surf_usm%ns))
               t_window_av%val = 0.0_wp
            end if

         case ('usm_t_green')
!
!--          Green temperature for iwl layer
            if (.not. allocated(t_green_av%val)) then
               allocate (t_green_av%val(nzb_wall:nzt_wall, 1:surf_usm%ns))
               t_green_av%val = 0.0_wp
            end if

         case ('usm_swc')
!
!--          Soil water content for iwl layer
            if (.not. allocated(swc_av%val)) then
               allocate (swc_av%val(nzb_wall:nzt_wall, 1:surf_usm%ns))
               swc_av%val = 0.0_wp
            end if

         case DEFAULT
            continue

         end select

      elseif (mode == 'sum') then

         select case (trim(var))

         case ('usm_wshf')
!
!--          Sensible heat flux
            call average_surfaces(wshf_eb_av%val, surf_usm%wshf_eb)

         case ('usm_qsws')
!
!--          Latent heat flux
            call average_surfaces(qsws_av%val, surf_usm%qsws)

         case ('usm_qsws_veg')
!
!--          Latent heat flux from vegetation surfaces
            call average_surfaces(qsws_veg_av%val, surf_usm%qsws_veg)

         case ('usm_qsws_liq')
!
!--          Latent heat flux from surfaces with liquid
            call average_surfaces(qsws_liq_av%val, surf_usm%qsws_liq)

         case ('usm_wghf')
!
!--           Heat flux from ground (wall)
            call average_surfaces(wghf_eb_av%val, surf_usm%wghf_eb)

         case ('usm_wghf_window')
!
!--          Heat flux from window ground
            call average_surfaces(wghf_eb_window_av%val, surf_usm%wghf_eb_window)

         case ('usm_wghf_green')
!
!--           Heat flux from green ground
            call average_surfaces(wghf_eb_green_av%val, surf_usm%wghf_eb_green)

         case ('usm_iwghf')
!
!--          Heat flux from indoor ground
            call average_surfaces(iwghf_eb_av%val, surf_usm%iwghf_eb)

         case ('usm_iwghf_window')
!
!--          Heat flux from indoor window ground
            call average_surfaces(iwghf_eb_window_av%val, surf_usm%iwghf_eb_window)

         case ('usm_t_surf_wall')
!
!--          Surface temperature of wall surfaces
            call average_surfaces(t_surf_wall_av%val, t_surf_wall%val)

         case ('usm_t_surf_window')
!
!--          Surface temperature for window surfaces
            call average_surfaces(t_surf_window_av%val, t_surf_window%val)

         case ('usm_t_surf_green')
!
!--           Surface temperature for green surfaces
            call average_surfaces(t_surf_green_av%val, t_surf_green%val)

         case ('usm_t_wall')
!
!--          Wall temperature for iwl layer
            call average_surfaces(t_wall_av%val(iwl, :), t_wall%val(iwl, :))

         case ('usm_t_window')
!
!--          Window temperature for iwl layer
            call average_surfaces(t_window_av%val(iwl, :), t_window%val(iwl, :))

         case ('usm_t_green')
!
!--          Green temperature for iwl layer
            call average_surfaces(t_green_av%val(iwl, :), t_green%val(iwl, :))

         case ('usm_swc')
!
!--          Soil water content for iwl layer
            call average_surfaces(swc_av%val(iwl, :), swc%val(iwl, :))

         case DEFAULT
            continue

         end select

      elseif (mode == 'average') then

         select case (trim(var))

         case ('usm_wshf')
!
!--          Sensible heat flux
            call average_surfaces(wshf_eb_av%val)

         case ('usm_qsws')
!
!--          Latent heat flux
            call average_surfaces(qsws_av%val)

         case ('usm_qsws_veg')
!
!--          Latent heat flux from vegetation surfaces
            call average_surfaces(qsws_veg_av%val)

         case ('usm_qsws_liq')
!
!--          Latent heat flux from surfaces with liquid
            call average_surfaces(qsws_liq_av%val)

         case ('usm_wghf')
!
!--          Heat flux from ground
            call average_surfaces(wghf_eb_av%val)

         case ('usm_wghf_window')
!
!--          Heat flux from window ground
            call average_surfaces(wghf_eb_window_av%val)

         case ('usm_wghf_green')
!
!--          Heat flux from green ground
            call average_surfaces(wghf_eb_green_av%val)

         case ('usm_iwghf')
!
!--          Heat flux from indoor ground
            call average_surfaces(iwghf_eb_av%val)

         case ('usm_iwghf_window')
!
!--          Heat flux from indoor window ground
            call average_surfaces(iwghf_eb_window_av%val)

         case ('usm_t_surf_wall')
!
!--          Surface temperature for wall surfaces
            call average_surfaces(t_surf_wall_av%val)

         case ('usm_t_surf_window')
!
!--          Surface temperature for window surfaces
            call average_surfaces(t_surf_window_av%val)

         case ('usm_t_surf_green')
!
!--          Surface temperature for green surfaces
            call average_surfaces(t_surf_green_av%val)

         case ('usm_t_wall')
!
!--          Wall temperature for iwl layer
            call average_surfaces(t_wall_av%val(iwl, :))

         case ('usm_t_window')
!
!--          Window temperature for iwl layer
            call average_surfaces(t_window_av%val(iwl, :))

         case ('usm_t_green')
!
!--          Green temperature for iwl layer
            call average_surfaces(t_green_av%val(iwl, :))

         case ('usm_swc')
!
!--          Soil water content for iwl layer
            call average_surfaces(swc_av%val(iwl, :))

         end select

      end if

   contains

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Average surface data accorrding to its facing.
!--------------------------------------------------------------------------------------------------!
      subroutine average_surfaces(mean_array, input_array)

         real(wp), dimension(1:surf_usm%ns), optional ::  input_array !< array to be averaged and summed-up
         real(wp), dimension(1:surf_usm%ns)           ::  mean_array  !< averaged and summed-up array

         if (mode == 'sum') then
!
!--       Sum-up surface array. Thereby, distinguish between different facings. This is
!--       necessary since the routine for a given quantity can be called several times
!--       (for each facing separately). If a surface element does not belong to the currently
!--       treated facing, just add a zero.
            if (upward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m) + merge(input_array(m), 0.0_wp, surf_usm%upward(m))
               end do
            elseif (downward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m) + merge(input_array(m), 0.0_wp, surf_usm%downward(m))
               end do
            elseif (eastward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m) + merge(input_array(m), 0.0_wp, surf_usm%eastward(m))
               end do
            elseif (westward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m) + merge(input_array(m), 0.0_wp, surf_usm%westward(m))
               end do
            elseif (northward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m) + merge(input_array(m), 0.0_wp, surf_usm%northward(m))
               end do
            elseif (southward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m) + merge(input_array(m), 0.0_wp, surf_usm%southward(m))
               end do
            end if

         elseif (mode == 'average') then
!
!--       Average the surface array. Thereby, distinguish between different facings. This is
!--       necessary since the routine for a given quantity can be called several times
!--       (for each facing separately). If a surface element does not belong to the currently
!--       treated facing, just divide it by one.
            if (upward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m)/ &
                                  merge(real(average_count_3d, KIND=wp), 1.0_wp, surf_usm%upward(m))
               end do
            elseif (downward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m)/ &
                                  merge(real(average_count_3d, KIND=wp), 1.0_wp, surf_usm%downward(m))
               end do
            elseif (eastward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m)/ &
                                  merge(real(average_count_3d, KIND=wp), 1.0_wp, surf_usm%eastward(m))
               end do
            elseif (westward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m)/ &
                                  merge(real(average_count_3d, KIND=wp), 1.0_wp, surf_usm%westward(m))
               end do
            elseif (northward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m)/ &
                                  merge(real(average_count_3d, KIND=wp), 1.0_wp, surf_usm%northward(m))
               end do
            elseif (southward) then
               do m = 1, surf_usm%ns
                  mean_array(m) = mean_array(m)/ &
                                  merge(real(average_count_3d, KIND=wp), 1.0_wp, surf_usm%southward(m))
               end do
            end if
         end if

      end subroutine average_surfaces

   end subroutine usm_3d_data_averaging

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set internal Neumann boundary condition at outer soil grid points for temperature and humidity.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_boundary_condition

      implicit none

      integer(iwp) ::  i      !< grid index x-direction
      integer(iwp) ::  ioff   !< offset index x-direction indicating location of soil grid point
      integer(iwp) ::  j      !< grid index y-direction
      integer(iwp) ::  joff   !< offset index x-direction indicating location of soil grid point
      integer(iwp) ::  k      !< grid index z-direction
      integer(iwp) ::  koff   !< offset index x-direction indicating location of soil grid point
      integer(iwp) ::  m      !< running index surface elements

      do m = 1, surf_usm%ns
         ioff = surf_usm%ioff(m)
         joff = surf_usm%joff(m)
         koff = surf_usm%koff(m)
         i = surf_usm%i(m)
         j = surf_usm%j(m)
         k = surf_usm%k(m)
         pt(k + koff, j + joff, i + ioff) = pt(k, j, i)
      end do

   end subroutine usm_boundary_condition

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine checks variables and assigns units.
!> It is called out from subroutine check_parameters.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_check_data_output(variable, unit)

      implicit none

      character(LEN=*), intent(IN)    ::  variable   !<
      character(LEN=*), intent(OUT)   ::  unit       !<

      character(LEN=2)                              ::  ls            !<

      character(LEN=varnamelength)                  ::  var           !< TRIM(variable)

      integer(iwp)                                  ::  i, j, l         !< index

      integer(iwp), parameter                       ::  nl1 = 14      !< number of directional usm variables
      character(LEN=varnamelength), dimension(nl1)  ::  varlist1 = &  !< list of directional usm variables
                                                       (/'usm_wshf                      ', &
                                                         'usm_wghf                      ', &
                                                         'usm_wghf_window               ', &
                                                         'usm_wghf_green                ', &
                                                         'usm_iwghf                     ', &
                                                         'usm_iwghf_window              ', &
                                                         'usm_surfz                     ', &
                                                         'usm_surfwintrans              ', &
                                                         'usm_surfcat                   ', &
                                                         'usm_t_surf_wall               ', &
                                                         'usm_t_surf_window             ', &
                                                         'usm_t_surf_green              ', &
                                                         'usm_t_green                   ', &
                                                         'usm_qsws                      '/)

      integer(iwp), parameter                       ::  nl2 = 3       !< number of directional layer usm variables
      character(LEN=varnamelength), dimension(nl2)  ::  varlist2 = &  !< list of directional layer usm variables
                                                       (/'usm_t_wall                    ', &
                                                         'usm_t_window                  ', &
                                                         'usm_t_green                   '/)

      logical                                       ::  lfound     !< flag if the variable is found

      lfound = .false.

      var = trim(variable)

!
!-- Check if variable exists
!-- Directional variables
      do i = 1, nl1
         do j = 0, nd - 1
            if (trim(var) == trim(varlist1(i))//trim(dirname(j))) then
               lfound = .true.
               exit
            end if
            if (lfound) exit
         end do
      end do
      if (lfound) goto 10
!
!-- Directional layer variables
      do i = 1, nl2
         do j = 0, nd - 1
            do l = nzb_wall, nzt_wall
               write (ls, '(A1,I1)') '_', l
               if (trim(var) == trim(varlist2(i))//trim(ls)//trim(dirname(j))) then
                  lfound = .true.
                  exit
               end if
            end do
            if (lfound) exit
         end do
      end do
      if (.not. lfound) then
         unit = 'illegal'
         return
      end if
10    continue

      if (var(1:9) == 'usm_wshf_' .or. var(1:9) == 'usm_wghf_' .or. &
          var(1:16) == 'usm_wghf_window_' .or. var(1:15) == 'usm_wghf_green_' .or. &
          var(1:10) == 'usm_iwghf_' .or. var(1:17) == 'usm_iwghf_window_' .or. &
          var(1:17) == 'usm_surfwintrans_' .or. &
          var(1:9) == 'usm_qsws_' .or. var(1:13) == 'usm_qsws_veg_' .or. &
          var(1:13) == 'usm_qsws_liq_') &
         then
         unit = 'W/m2'
      elseif (var(1:15) == 'usm_t_surf_wall' .or. var(1:10) == 'usm_t_wall' .or. &
              var(1:12) == 'usm_t_window' .or. var(1:17) == 'usm_t_surf_window' .or. &
              var(1:16) == 'usm_t_surf_green' .or. &
              var(1:11) == 'usm_t_green' .or. var(1:7) == 'usm_swc') &
         then
         unit = 'K'
      elseif (var(1:9) == 'usm_surfz' .or. var(1:11) == 'usm_surfcat') then
         unit = '1'
      else
         unit = 'illegal'
      end if

   end subroutine usm_check_data_output

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for urban surface model
!--------------------------------------------------------------------------------------------------!
   subroutine usm_check_parameters

      use control_parameters, &
         only: bc_pt_b, &
               bc_q_b, &
               constant_flux_layer, &
               large_scale_forcing, &
               lsf_surf, &
               topography

      implicit none

!
!-- Dirichlet boundary conditions are required as the surface fluxes are calculated from the
!-- temperature/humidity gradients in the urban surface model
      if (bc_pt_b == 'neumann' .or. bc_q_b == 'neumann') then
         message_string = 'urban surface model requires setting of bc_pt_b = "dirichlet" and '// &
                          'bc_q_b  = "dirichlet"'
         call message('usm_check_parameters', 'USM0001', 1, 2, 0, 6, 0)
      end if

      if (.not. constant_flux_layer) then
         message_string = 'urban surface model requires constant_flux_layer = .TRUE.'
         call message('usm_check_parameters', 'USM0002', 1, 2, 0, 6, 0)
      end if

      if (.not. radiation) then
         message_string = 'urban surface model requires the radiation model to be switched on'
         call message('usm_check_parameters', 'USM0003', 1, 2, 0, 6, 0)
      end if
!
!-- Surface forcing has to be disabled for LSF in case of enabled urban surface module
      if (large_scale_forcing) then
         lsf_surf = .false.
      end if
!
!-- Topography
      if (topography == 'flat') then
         message_string = 'topography /= "flat" is required when using the urban surface model'
         call message('usm_check_parameters', 'USM0004', 1, 2, 0, 6, 0)
      end if
!
!-- Check if building_type is set within a valid range. First, building_type is set via namelist.
      if (building_type < lbound(building_gen, 2) .and. &
          building_type > ubound(building_gen, 2)) &
         then
         write (message_string, *) 'building_type = ', building_type, ' is out of the valid range'
         call message('usm_check_parameters', 'USM0005', 2, 2, 0, 6, 0)
      end if
!
!-- Check if building_pars is correctly dimensioned.
      if (building_pars_f%from_file) then
         if (building_pars_f%np /= size(building_pars, 1)) then
            write (message_string, *) 'dimension size of static input variable building_pars is ', &
               building_pars_f%np, '&', &
               'dimension size of ', size(building_pars, 1), 'is required'
            call message('usm_check_parameters', 'USM0007', 2, 2, 0, 6, 0)
         end if
      end if

      call usm_check_parameters_input_bpars('building_albedo_type', building_alb_type_f, &
                                            size(building_alb_type, 1))
      call usm_check_parameters_input_bpars('building_emissivity', building_emis_f, &
                                            size(building_emis, 1))
      call usm_check_parameters_input_bpars('building_fraction', building_frac_f, &
                                            size(building_frac, 1))
      call usm_check_parameters_input_bpars('building_general_pars', building_gen_f, &
                                            size(building_gen, 1))
      call usm_check_parameters_input_bpars('building_heat_capacity', building_hcap_f, &
                                            size(building_hcap, 1), size(building_hcap, 2))
      call usm_check_parameters_input_bpars('building_heat_conductivity', building_hcond_f, &
                                            size(building_hcond, 1), size(building_hcond, 2))
      call usm_check_parameters_input_bpars('building_indoor_pars', building_indoor_f, &
                                            size(building_indoor, 1))
      call usm_check_parameters_input_bpars('building_lai', building_lai_f, &
                                            size(building_lai, 1))
      call usm_check_parameters_input_bpars('building_roughness_length', building_z0_f, &
                                            size(building_z0, 1))
      call usm_check_parameters_input_bpars('building_roughness_length_qh', building_z0qh_f, &
                                            size(building_z0qh, 1))
      call usm_check_parameters_input_bpars('building_thickness', building_thick_f, &
                                            size(building_depth, 1), size(building_depth, 2))
      call usm_check_parameters_input_bpars('building_transmissivity', building_trans_f, &
                                            size(building_trans, 1))

   end subroutine usm_check_parameters

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check read-in building parameters from static input file.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_check_parameters_input_bpars(bpars_name, bpars_f, bpars_size, bpars_size_layer)

      character(LEN=*), intent(IN)           ::  bpars_name         !< name of the input data structure

      integer(iwp), intent(IN)               ::  bpars_size         !< required size of the input data structure (type)
      integer(iwp), intent(IN), optional     ::  bpars_size_layer   !< required size of the input data structure (layer)

      type(pars), intent(IN)                 ::  bpars_f            !< input data structure

      if (.not. bpars_f%from_file) return

      if (bpars_f%np /= bpars_size) then
         write (message_string, '(A,I4,A,A,I4,A)') 'type dimension size (', bpars_f%np, &
            ') of static input variable "'//trim(bpars_name), &
            '" does not match the required size of (', bpars_size, ')'
         call message('usm_check_parameters', 'USM0016', 2, 2, 0, 6, 0)
      end if
      if (present(bpars_size_layer)) then
         if (bpars_f%nz /= bpars_size_layer) then
            write (message_string, '(A,I4,A,A,I4,A)') 'layer dimension size (', bpars_f%nz, &
               ') of static input variable "'//trim(bpars_name), &
               '" does not match the required size of (', bpars_size_layer, ')'
            call message('usm_check_parameters', 'USM0016', 2, 2, 0, 6, 0)
         end if
      end if

   end subroutine usm_check_parameters_input_bpars

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format for variables of urban_surface model.
!> It resorts the urban surface module output quantities from surf style indexing into temporary 3D
!> array with indices (i,j,k). It is called from subroutine data_output_3d.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_data_output_3d(av, variable, found, local_pf, nzb_do, nzt_do)

      implicit none

      character(LEN=*), intent(IN)   ::  variable  !< variable name

      character(LEN=varnamelength)   ::  var  !< trimmed variable name

      integer(iwp), intent(IN)       ::  av        !< flag if averaged
      integer(iwp), intent(IN)       ::  nzb_do    !< lower limit of the data output (usually 0)
      integer(iwp), intent(IN)       ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)

      integer(iwp)  ::  ids, idsint, idsidx        !<
      integer(iwp)  ::  i, j, k, iwl, istat, m     !< running indices

      logical ::  downward  !< control flag indicating the output of downward-facing surfaces
      logical ::  eastward  !< control flag indicating the output of east-facing surfaces
      logical ::  northward !< control flag indicating the output of northward-facing surfaces
      logical ::  southward !< control flag indicating the output of southward-facing surfaces
      logical ::  upward    !< control flag indicating the output of upward-facing surfaces
      logical ::  westward  !< control flag indicating the output of westward-facing surfaces

      logical, intent(OUT) ::  found  !<

      real(wp), dimension(nxl:nxr, nys:nyn, nzb_do:nzt_do) ::  local_pf  !< sp - it has to correspond to module data_output_3d
      real(wp), dimension(nzb:nzt + 1, nys:nyn, nxl:nxr)     ::  temp_pf   !< temp array for urban surface output procedure

      found = .true.
      if (.not. data_output_raw) temp_pf = output_fill_value

      ids = -1
      var = trim(variable)
      do i = 0, nd - 1
         k = len(trim(var))
         j = len(trim(dirname(i)))
         if (trim(var(k - j + 1:k)) == trim(dirname(i))) then
            ids = i
            idsint = dirint(ids)
            idsidx = diridx(ids)
            var = var(:k - j)
            exit
         end if
      end do
!
!-- Set direction control flags
      downward = .false.
      eastward = .false.
      northward = .false.
      southward = .false.
      upward = .false.
      westward = .false.
      if (idsint == iup) then
         upward = .true.
      elseif (idsint == idown) then
         downward = .true.
      elseif (idsint == ieast) then
         eastward = .true.
      elseif (idsint == iwest) then
         westward = .true.
      elseif (idsint == inorth) then
         northward = .true.
      elseif (idsint == isouth) then
         southward = .true.
      end if

      if (ids == -1) then
         var = trim(variable)
      end if
      if (var(1:11) == 'usm_t_wall_' .and. len(trim(var)) >= 12) then
!
!--     Wall layers
         read (var(12:12), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:10)
         end if
      end if
      if (var(1:13) == 'usm_t_window_' .and. len(trim(var)) >= 14) then
!
!--     Window layers
         read (var(14:14), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:12)
         end if
      end if
      if (var(1:12) == 'usm_t_green_' .and. len(trim(var)) >= 13) then
!
!--     Green layers
         read (var(13:13), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:11)
         end if
      end if
      if (var(1:8) == 'usm_swc_' .and. len(trim(var)) >= 9) then
!
!--     Green layers soil water content
         read (var(9:9), '(I1)', iostat=istat) iwl
         if (istat == 0 .and. iwl >= nzb_wall .and. iwl <= nzt_wall) then
            var = var(1:7)
         end if
      end if

      select case (trim(var))

      case ('usm_surfz')
!
!--       Array of surface height (z)
!--       Write surface array to temp_pf. Thereby, distinguish between different facings. This is
!--       necessary since the routine for a given quantity can be called several times
!--       (for each facing separately). If a surface element does not belong to the currently
!--       treated facing, do not modify temp_pf's current value.
         if (upward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(max(temp_pf(0, j, i), real(k, KIND=wp)), &
                                        temp_pf(k, j, i), surf_usm%upward(m))
            end do
         elseif (downward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(max(temp_pf(0, j, i), real(k, KIND=wp)), &
                                        temp_pf(k, j, i), surf_usm%downward(m))
            end do
         elseif (eastward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(max(temp_pf(0, j, i), real(k, KIND=wp) + 1.0_sp), &
                                        temp_pf(k, j, i), surf_usm%eastward(m))
            end do
         elseif (westward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(max(temp_pf(0, j, i), real(k, KIND=wp) + 1.0_sp), &
                                        temp_pf(k, j, i), surf_usm%westward(m))
            end do
         elseif (northward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(max(temp_pf(0, j, i), real(k, KIND=wp) + 1.0_sp), &
                                        temp_pf(k, j, i), surf_usm%northward(m))
            end do
         elseif (southward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(max(temp_pf(0, j, i), real(k, KIND=wp) + 1.0_sp), &
                                        temp_pf(k, j, i), surf_usm%southward(m))
            end do
         end if

      case ('usm_surfcat')
!
!--       Surface category
         call write_surface_data_to_temp_pf(real(surf_usm%surface_types, KIND=wp))

      case ('usm_surfwintrans')
!
!--       Transmissivity window tiles
         call write_surface_data_to_temp_pf(surf_usm%transmissivity)

      case ('usm_wshf')
!
!--       Sensible heat flux
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%wshf_eb)
         else
            call write_surface_data_to_temp_pf(wshf_eb_av%val)
         end if

      case ('usm_qsws')
!
!--       Latent heat flux
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%qsws)
         else
            call write_surface_data_to_temp_pf(qsws_av%val)
         end if

      case ('usm_qsws_veg')
!
!--       Latent heat flux from vegetation surfaces
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%qsws_veg)
         else
            call write_surface_data_to_temp_pf(qsws_veg_av%val)
         end if

      case ('usm_qsws_liq')
!
!--       Latent heat flux from surfaces with liquid
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%qsws_liq)
         else
            call write_surface_data_to_temp_pf(qsws_liq_av%val)
         end if

      case ('usm_wghf')
!
!--       Heat flux from ground
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%wghf_eb)
         else
            call write_surface_data_to_temp_pf(wghf_eb_av%val)
         end if

      case ('usm_wghf_window')
!
!--       Heat flux from window ground
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%wghf_eb_window)
         else
            call write_surface_data_to_temp_pf(wghf_eb_window_av%val)
         end if

      case ('usm_wghf_green')
!
!--       Heat flux from green ground
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%wghf_eb_green)
         else
            call write_surface_data_to_temp_pf(wghf_eb_green_av%val)
         end if

      case ('usm_iwghf')
!
!--       Heat flux from indoor ground
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%iwghf_eb)
         else
            call write_surface_data_to_temp_pf(iwghf_eb_av%val)
         end if

      case ('usm_iwghf_window')
!
!--       Heat flux from indoor window ground
         if (av == 0) then
            call write_surface_data_to_temp_pf(surf_usm%iwghf_eb_window)
         else
            call write_surface_data_to_temp_pf(iwghf_eb_window_av%val)
         end if

      case ('usm_t_surf_wall')
!
!--       Surface temperature for wall surfaces
         if (av == 0) then
            call write_surface_data_to_temp_pf(t_surf_wall%val)
         else
            call write_surface_data_to_temp_pf(t_surf_wall_av%val)
         end if

      case ('usm_t_surf_window')
!
!--       Surface temperature for window surfaces
         if (av == 0) then
            call write_surface_data_to_temp_pf(t_surf_window%val)
         else
            call write_surface_data_to_temp_pf(t_surf_window_av%val)
         end if

      case ('usm_t_surf_green')
!
!--       Surface temperature for green surfaces
         if (av == 0) then
            call write_surface_data_to_temp_pf(t_surf_green%val)
         else
            call write_surface_data_to_temp_pf(t_surf_green_av%val)
         end if

      case ('usm_t_wall')
!
!--       Wall temperature for iwl layer of walls
         if (av == 0) then
            call write_surface_data_to_temp_pf(t_wall%val(iwl, :))
         else
            call write_surface_data_to_temp_pf(t_wall_av%val(iwl, :))
         end if

      case ('usm_t_window')
!
!--       Window temperature for iwl layer
         if (av == 0) then
            call write_surface_data_to_temp_pf(t_window%val(iwl, :))
         else
            call write_surface_data_to_temp_pf(t_window_av%val(iwl, :))
         end if

      case ('usm_t_green')
!
!--       Green temperature for iwl layer
         if (av == 0) then
            call write_surface_data_to_temp_pf(t_green%val(iwl, :))
         else
            call write_surface_data_to_temp_pf(t_green_av%val(iwl, :))
         end if

      case ('usm_swc')
!
!--       Soil water content for iwl layer
         if (av == 0) then
            call write_surface_data_to_temp_pf(swc%val(iwl, :))
         else
            call write_surface_data_to_temp_pf(swc_av%val(iwl, :))
         end if

      case DEFAULT
         found = .false.
         return

      end select

!
!-- Rearrange dimensions for NetCDF output
!-- FIXME: this may generate FPE overflow upon conversion from DP to SP
      do j = nys, nyn
         do i = nxl, nxr
            do k = nzb_do, nzt_do
               local_pf(i, j, k) = temp_pf(k, j, i)
            end do
         end do
      end do

   contains

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write surface data onto output array accorrding to its facing.
!--------------------------------------------------------------------------------------------------!
      subroutine write_surface_data_to_temp_pf(surf_array)

         real(wp), dimension(1:surf_usm%ns) ::  surf_array !< treated surface array
!
!--    Write surface array to temp_pf. Thereby, distinguish between different facings. This is
!--    necessary since the routine for a given quantity can be called several times
!--    (for each facing separately). If a surface element does not belong to the currently
!--    treated facing, do not modify temp_pf's current value.
         if (upward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(surf_array(m), temp_pf(k, j, i), surf_usm%upward(m))
            end do
         elseif (downward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(surf_array(m), temp_pf(k, j, i), surf_usm%downward(m))
            end do
         elseif (eastward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(surf_array(m), temp_pf(k, j, i), surf_usm%eastward(m))
            end do
         elseif (westward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(surf_array(m), temp_pf(k, j, i), surf_usm%westward(m))
            end do
         elseif (northward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(surf_array(m), temp_pf(k, j, i), surf_usm%northward(m))
            end do
         elseif (southward) then
            do m = 1, surf_usm%ns
               i = surf_usm%i(m)
               j = surf_usm%j(m)
               k = surf_usm%k(m)
               temp_pf(k, j, i) = merge(surf_array(m), temp_pf(k, j, i), surf_usm%southward(m))
            end do
         end if

      end subroutine write_surface_data_to_temp_pf

   end subroutine usm_data_output_3d

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine defines appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_define_netcdf_grid(variable, found, grid_x, grid_y, grid_z)

      implicit none

      character(LEN=*), intent(IN)  ::  variable  !<
      character(LEN=*), intent(OUT) ::  grid_x    !<
      character(LEN=*), intent(OUT) ::  grid_y    !<
      character(LEN=*), intent(OUT) ::  grid_z    !<

      character(LEN=varnamelength)  ::  var  !<

      logical, intent(OUT)  ::  found  !<

      var = trim(variable)
      if (var(1:9) == 'usm_wshf_' .or. var(1:9) == 'usm_wghf_' .or. &
          var(1:16) == 'usm_wghf_window_' .or. var(1:15) == 'usm_wghf_green_' .or. &
          var(1:10) == 'usm_iwghf_' .or. var(1:17) == 'usm_iwghf_window_' .or. &
          var(1:9) == 'usm_qsws_' .or. var(1:13) == 'usm_qsws_veg_' .or. &
          var(1:13) == 'usm_qsws_liq_' .or. &
          var(1:15) == 'usm_t_surf_wall' .or. var(1:10) == 'usm_t_wall' .or. &
          var(1:17) == 'usm_t_surf_window' .or. var(1:12) == 'usm_t_window' .or. &
          var(1:16) == 'usm_t_surf_green' .or. var(1:11) == 'usm_t_green' .or. &
          var(1:9) == 'usm_surfz' .or. var(1:11) == 'usm_surfcat' .or. &
          var(1:16) == 'usm_surfwintrans' .or. var(1:7) == 'usm_swc') then

         found = .true.
         grid_x = 'x'
         grid_y = 'y'
         grid_z = 'zu'
      else
         found = .false.
         grid_x = 'none'
         grid_y = 'none'
         grid_z = 'none'
      end if

   end subroutine usm_define_netcdf_grid

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wall surface model
!--------------------------------------------------------------------------------------------------!
   subroutine usm_init_wall_heat_model

      implicit none

      integer(iwp) ::  k  !< running index along z-dimension
      integer(iwp) ::  m  !< running index for surface elements

      if (debug_output) call debug_message('usm_init_wall_heat_model', 'start')

!
!-- Calculate wall and window grid spacings. Wall temperature is defined at the center of the
!-- wall layers.
      do m = 1, surf_usm%ns
!
!--    Set-up wall layer discretization
         surf_usm%dz_wall(nzb_wall, m) = surf_usm%zw(nzb_wall, m)
         do k = nzb_wall + 1, nzt_wall
            surf_usm%dz_wall(k, m) = surf_usm%zw(k, m) - surf_usm%zw(k - 1, m)
         end do

         do k = nzb_wall, nzt_wall - 1
            surf_usm%dz_wall_center(k, m) = 0.5_wp* &
                                            (surf_usm%dz_wall(k, m) + surf_usm%dz_wall(k + 1, m))
            if (surf_usm%dz_wall_center(k, m) <= 0.0_wp) then
               write (message_string, '(A,I5,A)') 'invalid wall layer configuration found: '// &
                  'dz_wall_center(k=', k, ') <= 0.0'
               call message('usm_init_wall_heat_model', 'USM0008', 1, 2, 0, 6, 0)
            end if
         end do
         surf_usm%dz_wall_center(nzt_wall, m) = surf_usm%dz_wall(nzt_wall, m)

!
!--    Set-up window layer discretization
         surf_usm%dz_window(nzb_wall, m) = surf_usm%zw_window(nzb_wall, m)
         do k = nzb_wall + 1, nzt_wall
            surf_usm%dz_window(k, m) = surf_usm%zw_window(k, m) - surf_usm%zw_window(k - 1, m)
         end do

         do k = nzb_wall, nzt_wall - 1
            surf_usm%dz_window_center(k, m) = 0.5_wp* &
                                              (surf_usm%dz_window(k, m) + surf_usm%dz_window(k + 1, m))
            if (surf_usm%dz_window_center(k, m) <= 0.0_wp) then
               write (message_string, '(A,I5,A)') 'invalid window layer configuration found: '// &
                  'dz_window_center(k=', k, ') <= 0.0'
               call message('usm_init_wall_heat_model', 'USM0009', 1, 2, 0, 6, 0)
            end if
         end do

         surf_usm%dz_window_center(nzt_wall, m) = surf_usm%dz_window(nzt_wall, m)

!
!--    Set-up green roofs
         if (surf_usm%green_type_roof(m) == 2.0_wp) then
!
!--       Extensive green roof
!--       Set ratio of substrate layer thickness, soil-type and LAI
            soil_type = 3
            surf_usm%lai(m) = 2.0_wp
            surf_usm%zw_green(nzb_wall, m) = 0.05_wp
            surf_usm%zw_green(nzb_wall + 1, m) = 0.10_wp
            surf_usm%zw_green(nzb_wall + 2, m) = 0.15_wp
            surf_usm%zw_green(nzb_wall + 3, m) = 0.20_wp
         else
!
!--       Intensive green roof
!--       Set ratio of substrate layer thickness, soil-type and LAI
            soil_type = 6
            surf_usm%lai(m) = 4.0_wp
            surf_usm%zw_green(nzb_wall, m) = 0.05_wp
            surf_usm%zw_green(nzb_wall + 1, m) = 0.10_wp
            surf_usm%zw_green(nzb_wall + 2, m) = 0.40_wp
            surf_usm%zw_green(nzb_wall + 3, m) = 0.80_wp
         end if

         surf_usm%dz_green(nzb_wall, m) = surf_usm%zw_green(nzb_wall, m)
         do k = nzb_wall + 1, nzt_wall
            surf_usm%dz_green(k, m) = surf_usm%zw_green(k, m) - surf_usm%zw_green(k - 1, m)
         end do

         do k = nzb_wall, nzt_wall - 1
            surf_usm%dz_green_center(k, m) = 0.5_wp* &
                                             (surf_usm%dz_green(k, m) + surf_usm%dz_green(k + 1, m))
            if (surf_usm%dz_green_center(k, m) <= 0.0_wp) then
               write (message_string, '(A,I5,A)') 'invalid green layer configuration found: '// &
                  'dz_green_center(k=', k, ') <= 0.0'
               call message('usm_init_wall_heat_model', 'USM0010', 1, 2, 0, 6, 0)
            end if
         end do

         surf_usm%dz_green_center(nzt_wall, m) = surf_usm%dz_green(nzt_wall, m)

         if (alpha_vangenuchten == 9999999.9_wp) then
            alpha_vangenuchten = soil_pars(0, soil_type)
         end if

         if (l_vangenuchten == 9999999.9_wp) then
            l_vangenuchten = soil_pars(1, soil_type)
         end if

         if (n_vangenuchten == 9999999.9_wp) then
            n_vangenuchten = soil_pars(2, soil_type)
         end if

         if (hydraulic_conductivity == 9999999.9_wp) then
            hydraulic_conductivity = soil_pars(3, soil_type)
         end if

         if (saturation_moisture == 9999999.9_wp) then
            saturation_moisture = m_soil_pars(0, soil_type)
         end if

         if (field_capacity == 9999999.9_wp) then
            field_capacity = m_soil_pars(1, soil_type)
         end if

         if (wilting_point == 9999999.9_wp) then
            wilting_point = m_soil_pars(2, soil_type)
         end if

         if (residual_moisture == 9999999.9_wp) then
            residual_moisture = m_soil_pars(3, soil_type)
         end if

         if (trim(initializing_actions) /= 'read_restart_data' .and. .not. read_spinup_data) &
            then
            do k = nzb_wall, nzt_wall + 1
               swc%val(k, m) = field_capacity
            end do
         end if

         do k = nzb_wall, nzt_wall + 1
            rootfr%val(k, m) = 0.5_wp
            surf_usm%alpha_vg_green(m) = alpha_vangenuchten
            surf_usm%l_vg_green(m) = l_vangenuchten
            surf_usm%n_vg_green(m) = n_vangenuchten
            surf_usm%gamma_w_green_sat(k, m) = hydraulic_conductivity
            swc_sat%val(k, m) = saturation_moisture
            fc%val(k, m) = field_capacity
            wilt%val(k, m) = wilting_point
            swc_res%val(k, m) = residual_moisture
         end do

      end do

      surf_usm%ddz_wall = 1.0_wp/surf_usm%dz_wall
      surf_usm%ddz_wall_center = 1.0_wp/surf_usm%dz_wall_center
      surf_usm%ddz_window = 1.0_wp/surf_usm%dz_window
      surf_usm%ddz_window_center = 1.0_wp/surf_usm%dz_window_center
      surf_usm%ddz_green = 1.0_wp/surf_usm%dz_green
      surf_usm%ddz_green_center = 1.0_wp/surf_usm%dz_green_center

!
!-- Calculate wall heat conductivity (lambda_h) at the _layer level the weighted average
      do m = 1, surf_usm%ns
         do k = nzb_wall, nzt_wall - 1
            surf_usm%lambda_h_layer(k, m) = (surf_usm%lambda_h(k, m)*surf_usm%dz_wall(k, m) + &
                                             surf_usm%lambda_h(k + 1, m)*surf_usm%dz_wall(k + 1, m) &
                                             )*0.5_wp*surf_usm%ddz_wall_center(k, m)
         end do
         surf_usm%lambda_h_layer(nzt_wall, m) = surf_usm%lambda_h(nzt_wall, m)
      end do

      do m = 1, surf_usm%ns
!
!--    Calculate wall heat conductivity (lambda_h) at the _layer level using weighting
         do k = nzb_wall, nzt_wall - 1
            surf_usm%lambda_h_window_layer(k, m) = (surf_usm%lambda_h_window(k, m)* &
                                                    surf_usm%dz_window(k, m) + &
                                                    surf_usm%lambda_h_window(k + 1, m)* &
                                                    surf_usm%dz_window(k + 1, m) &
                                                    )*0.5_wp*surf_usm%ddz_window_center(k, m)
         end do
         surf_usm%lambda_h_window_layer(nzt_wall, m) = surf_usm%lambda_h_window(nzt_wall, m)
      end do

      if (debug_output) call debug_message('usm_init_wall_heat_model', 'end')

   end subroutine usm_init_wall_heat_model

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the urban surface model
!--------------------------------------------------------------------------------------------------!
   subroutine usm_init

      use arrays_3d, &
         only: zw

      implicit none

      integer(iwp) ::  bt                !< short for building type
      integer(iwp) ::  exit_index        !< to store surface element index where z0 limit is exceeded
      integer(iwp) ::  i                 !< loop index x-dirction
      integer(iwp) ::  ind_green         !< index in input list for green
      integer(iwp) ::  ind_wall          !< index in input list for wall
      integer(iwp) ::  ind_win           !< index in input list for window
      integer(iwp) ::  ind_green_frac_w  !< index in input list for green fraction on wall
      integer(iwp) ::  ind_lai_w         !< index in input list for LAI on wall
      integer(iwp) ::  ind_level         !< index in input list for level
      integer(iwp) ::  ilayer            !< loop index input surface element layer
      integer(iwp) ::  is                !< loop index input surface element
      integer(iwp) ::  j                 !< loop index y-dirction
      integer(iwp) ::  k                 !< loop index z-dirction
      integer(iwp) ::  m                 !< loop index surface element

      logical ::  flag_exceed_z0           !< dummy flag to indicate whether roughness length is too high
      logical ::  flag_exceed_z0h          !< dummy flag to indicate whether roughness length for temperature is too high
      logical ::  flag_exceed_z0q          !< dummy flag to indicate whether roughness length for mositure is too high
      logical ::  relative_fraction_error  !< flag indicating if relative surface fractions do not sum up to 1

      real(wp) ::  c                     !<
      real(wp) ::  ground_floor_level_l  !< local height of ground floor level
      real(wp) ::  tin                   !<
      real(wp) ::  twin                  !<
      real(wp) ::  z_agl                 !< height of the surface element above terrain

      if (debug_output) call debug_message('usm_init', 'start')

      call cpu_log(log_point_s(78), 'usm_init', 'start')
!
!-- Surface forcing has to be disabled for LSF in case of enabled urban surface module
      if (large_scale_forcing) then
         lsf_surf = .false.
      end if
!
!-- Store wall indices on surface structure.
      surf_usm%nzb_wall = nzb_wall
      surf_usm%nzt_wall = nzt_wall
!
!-- Calculate constant values
      d_roughness_concrete = 1.0_wp/roughness_concrete
!
!-- Data input from static input file.
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_pars', 'nbuilding_pars', &
                                             parlist=building_pars_f)
!
!-- The deprecated list building_pars has been already read in netcdf_data_input_mod. In case this
!-- list exists in the driver file, give an informative message.
      if (building_pars_f%from_file) then
         message_string = 'building_pars from the static input file is deprecated'
         call message('usm_init', 'USM0015', 0, 0, 0, 6, 0)
      end if

!
!-- Static input of separated building-parameter lists.
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_albedo_type', 'building_surface_type', &
                                             parlist=building_alb_type_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_emissivity', 'building_surface_type', &
                                             parlist=building_emis_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_fraction', 'building_surface_type', &
                                             parlist=building_frac_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_general_pars', 'building_general_par', &
                                             parlist=building_gen_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_heat_capacity', &
                                             'building_surface_type', 'building_surface_layer', &
                                             building_hcap_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_heat_conductivity', &
                                             'building_surface_type', 'building_surface_layer', &
                                             building_hcond_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_indoor_pars', 'building_indoor_par', &
                                             parlist=building_indoor_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_lai', 'building_surface_level', &
                                             parlist=building_lai_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_roughness_length', 'building_surface_level', &
                                             parlist=building_z0_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_roughness_length_qh', &
                                             'building_surface_level', &
                                             parlist=building_z0qh_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_thickness', &
                                             'building_surface_type', 'building_surface_layer', &
                                             building_thick_f)
      call netcdf_data_input_parameter_lists(trim(input_file_static)//trim(coupling_char), &
                                             'building_transmissivity', 'building_surface_level', &
                                             parlist=building_trans_f)
!
!-- Assign building type to internal array.
!-- Level 1 initialization.
      if (surf_usm%ns > 0) surf_usm%btype = building_type
!
!-- Level 2 initialization - building type from static input file. In case of cut-cell topography,
!-- btype has already been set to its final value.
      if (building_type_f%from_file .and. .not. cut_cell_topography) then
         do m = 1, surf_usm%ns
            i = surf_usm%i(m) + surf_usm%ioff(m)
            j = surf_usm%j(m) + surf_usm%joff(m)

            if (building_type_f%var(j, i) /= building_type_f%fill) &
               surf_usm%btype(m) = building_type_f%var(j, i)
         end do
      end if

!
!-- Check for consistent setting of building type. This check needs to be done here due to different
!-- points of initialization of this array. It is not yet availabe in usm_check_paraters in case
!-- of the cut-cell method.
      do m = 1, surf_usm%ns
         if (surf_usm%btype(m) < lbound(building_pars, 2) .or. &
             surf_usm%btype(m) > ubound(building_pars, 2)) then
            i = surf_usm%i(m) + surf_usm%ioff(m)
            j = surf_usm%j(m) + surf_usm%joff(m)
            write (message_string, *) 'building_type = is out of the valid range at (i,j) ', &
               ' = (', i, ',', j, ')'
            call message('usm_init', 'USM0006', 2, 2, myid, 6, 0)
         end if
      end do

!
!-- Flag surface elements belonging to the ground floor level. Therefore, use terrain height array
!-- from file, if available. This flag is later used to control initialization of surface attributes.
      surf_usm%gfl = .false.

      do m = 1, surf_usm%ns
         i = surf_usm%i(m) + surf_usm%ioff(m)
         j = surf_usm%j(m) + surf_usm%joff(m)
         k = surf_usm%k(m)
!
!--    Determine local ground level. Level 1 - default value, level 2 - initialization according
!--    to building type, level 3 - initialization from value in building_gen_f (also from
!--    building_pars_f, but deprecated).
         ground_floor_level_l = ground_floor_level

         ground_floor_level_l = building_gen(ind_gflh, surf_usm%btype(m))

         if (building_pars_f%from_file) then
            if (building_pars_f%pars_xy(ind_gflh, j, i) /= building_pars_f%fill) &
               ground_floor_level_l = building_pars_f%pars_xy(ind_gflh, j, i)
         end if
         if (building_gen_f%from_file) then
            if (building_gen_f%pars_xy(ind_gflh, j, i) /= building_gen_f%fill) &
               ground_floor_level_l = building_gen_f%pars_xy(ind_gflh, j, i)
         end if
!
!--    Determine height of surface element above ground level. Please note, the height of a
!--    surface element is determined with respect to its height above ground of the reference
!--    grid point in the atmosphere. As (j,i) are defined as the surface indices and not
!--    the reference grid point, substract the offset values when assessing the terrain height.
         if (terrain_height_f%from_file) then
            z_agl = zw(k) - terrain_height_f%var(j - surf_usm%joff(m), i - surf_usm%ioff(m))
         else
            z_agl = zw(k)
         end if
!
!--    Set flag for ground level
         if (z_agl <= ground_floor_level_l) surf_usm%gfl(m) = .true.
!
!--    Set ground-floor level attribute to False at horizontally upward-facing surfaces.
!--    This is because ground-floor level properties are not necessarily valid for roofs,
!--    only valid for facades.
         if (surf_usm%upward(m)) surf_usm%gfl(m) = .false.
      end do
!
!-- Initialization of resistances.
      do m = 1, surf_usm%ns
         surf_usm%r_a(m) = 50.0_wp
         surf_usm%r_a_green(m) = 50.0_wp
         surf_usm%r_a_window(m) = 50.0_wp
      end do
!
!-- Initialization of canopy properties
      do m = 1, surf_usm%ns
         surf_usm%r_canopy(m) = 200.0_wp !< canopy_resistance
         surf_usm%r_canopy_min(m) = 200.0_wp !< min_canopy_resistance
         surf_usm%g_d(m) = 0.0_wp   !< canopy_resistance_coefficient
      end do
!
!-- Initialize building_type.
      do m = 1, surf_usm%ns
         surf_usm%building_type(m) = surf_usm%btype(m)
      end do
!
!-- Check if building_type is in a valid range.
      if (surf_usm%ns >= 1) then
         if (any(surf_usm%building_type(:) == 0 .or. &
                 surf_usm%building_type(:) > ubound(building_gen, 2))) then
            message_string = 'building_type is out of the range of possible values'
            call message('usm_init', 'USM0014', 2, 2, 0, 6, 0)
         end if
      end if
!
!-- Initialize urban-type surface attributes. According to initialization in land-surface model,
!-- follow a 3-level approach.
!-- Level 1 and 2: initialization of surfaces via default attributes from the bulk parameter lists,
!-- where the building type is the same for all buildings in the level 1 approach and can vary
!-- for each j,i grid point in the level 2 approach.
!-- The following loop encompasses the level 1 and 2 approach.
      do m = 1, surf_usm%ns

         i = surf_usm%i(m) + surf_usm%ioff(m)
         j = surf_usm%j(m) + surf_usm%joff(m)
!
!--    Set indices and specific variables for horizontal surfaces. Horizontal surfaces are
!--    considered as roof surfaces. If necessary, distinguish between ground-floor and
!--    above-ground floor level.
         if ((.not. cut_cell_topography .and. &
              (surf_usm%upward(m) .or. surf_usm%downward(m))) .or. &
             (cut_cell_topography .and. surf_usm%cut_cell_roof(m))) then

            if (surf_usm%upward(m) .or. surf_usm%cut_cell_roof(m)) then
               surf_usm%isroof_surf(m) = .true.
               surf_usm%surface_types(m) = roof_category
            end if

            ind_green = ind_green_roof
            ind_wall = ind_wall_roof
            ind_win = ind_win_roof
            ind_level = ind_roof
!
!--    Set indices and specific variables for vertical surfaces.
         else
            surf_usm%surface_types(m) = wall_category

            ind_green = merge(ind_green_gfl, ind_green_agfl, surf_usm%gfl(m))
            ind_wall = merge(ind_wall_gfl, ind_wall_agfl, surf_usm%gfl(m))
            ind_win = merge(ind_win_gfl, ind_win_agfl, surf_usm%gfl(m))
            ind_level = merge(ind_gfl, ind_agfl, surf_usm%gfl(m))
         end if

         bt = surf_usm%btype(m)
!
!--    Store name of building type.
         surf_usm%building_type_name(m) = building_type_name(bt)
!
!--    Initialize relative wall- (0), green- (1) and window (2) fractions.
         surf_usm%frac(m, ind_veg_wall) = building_frac(ind_wall, bt)
         surf_usm%frac(m, ind_pav_green) = building_frac(ind_green, bt)
         surf_usm%frac(m, ind_wat_win) = building_frac(ind_win, bt)
!
!--    Intialize LAI.
         surf_usm%lai(m) = building_lai(ind_level, bt)
!
!--    Intialize green roof type (intensive or extensive).
         surf_usm%green_type_roof(m) = building_gen(ind_green_type_roof, bt)

         do ilayer = 0, nzw - 1
!
!--       Intialize heat capacities.
            surf_usm%rho_c_wall(nzb_wall + ilayer, m) = building_hcap(ind_wall, ilayer, bt)
            surf_usm%rho_c_window(nzb_wall + ilayer, m) = building_hcap(ind_win, ilayer, bt)
            surf_usm%rho_c_green(nzb_wall + ilayer, m) = building_hcap(ind_green, ilayer, bt)
!
!--       Intialize heat conductivities.
            surf_usm%lambda_h(nzb_wall + ilayer, m) = building_hcond(ind_wall, ilayer, bt)
            surf_usm%lambda_h_window(nzb_wall + ilayer, m) = building_hcond(ind_win, ilayer, bt)
            surf_usm%lambda_h_green(nzb_wall + ilayer, m) = building_hcond(ind_green, ilayer, bt)
!
!--       Wall layer depths.
            surf_usm%zw(nzb_wall + ilayer, m) = building_depth(ind_wall, ilayer, bt)
            surf_usm%zw_window(nzb_wall + ilayer, m) = building_depth(ind_win, ilayer, bt)
            surf_usm%zw_green(nzb_wall + ilayer, m) = building_depth(ind_green, ilayer, bt)
         end do

         surf_usm%target_temp_summer(m) = building_indoor(ind_theta_int_c_set, bt)
         surf_usm%target_temp_winter(m) = building_indoor(ind_theta_int_h_set, bt)
!
!--    Emissivity of wall-, green- and window fraction.
         surf_usm%emissivity(m, ind_veg_wall) = building_emis(ind_wall, bt)
         surf_usm%emissivity(m, ind_pav_green) = building_emis(ind_green, bt)
         surf_usm%emissivity(m, ind_wat_win) = building_emis(ind_win, bt)
!
!--    Albedo type for wall fraction, green fraction, window fraction.
         surf_usm%albedo_type(m, ind_veg_wall) = building_alb_type(ind_wall, bt)
         surf_usm%albedo_type(m, ind_pav_green) = building_alb_type(ind_green, bt)
         surf_usm%albedo_type(m, ind_wat_win) = building_alb_type(ind_win, bt)
!
!--    Transmissivity at windows.
         surf_usm%transmissivity(m) = building_trans(ind_level, bt)
!
!--    Roughness lengths.
         surf_usm%z0(m) = building_z0(ind_level, bt)
         surf_usm%z0h(m) = building_z0qh(ind_level, bt)
         surf_usm%z0q(m) = building_z0qh(ind_level, bt)
      end do
!
!-- Initialize albedo type via given type from static input file. Please note, even though
!-- the albedo type has been already given by the bulk parameter list, albedo_type overwrites these
!-- values at this moment.
      if (albedo_type_f%from_file) then
         do m = 1, surf_usm%ns
            i = surf_usm%i(m) + surf_usm%ioff(m)
            j = surf_usm%j(m) + surf_usm%joff(m)
            if (albedo_type_f%var(j, i) /= albedo_type_f%fill) &
               surf_usm%albedo_type(m, :) = albedo_type_f%var(j, i)
         end do
      end if
!
!-- Level 3 - initialization via customized parameter lits.
!-- The level 3 initialization via building_pars is deprecated and incomplete (and will be revoved
!-- in future versions). Initialization via deprecated building_pars is done in
!-- usm_init_building_pars_pids.
      if (building_pars_f%from_file) call usm_init_building_pars_pids
!
!-- Initialize from the separated lists.
      do m = 1, surf_usm%ns

         i = surf_usm%i(m) + surf_usm%ioff(m)
         j = surf_usm%j(m) + surf_usm%joff(m)
!
!--    Set indices and specific variables for horizontal surfaces. Horizontal surfaces are
!--    considered as roof surfaces. If necessary, distinguish between ground-floor and
!--    above-ground floor level.
         if ((.not. cut_cell_topography .and. &
              (surf_usm%upward(m) .or. surf_usm%downward(m))) .or. &
             (cut_cell_topography .and. surf_usm%cut_cell_roof(m))) then

            ind_green = ind_green_roof
            ind_wall = ind_wall_roof
            ind_win = ind_win_roof
            ind_level = ind_roof
!
!--    Set indices and specific variables for vertical surfaces or explicitly defined
!--    wall-surfaces.
         else
            ind_green = merge(ind_green_gfl, ind_green_agfl, surf_usm%gfl(m))
            ind_wall = merge(ind_wall_gfl, ind_wall_agfl, surf_usm%gfl(m))
            ind_win = merge(ind_win_gfl, ind_win_agfl, surf_usm%gfl(m))
            ind_level = merge(ind_gfl, ind_agfl, surf_usm%gfl(m))
         end if

!
!--    Initialize relative wall, window and green fractions.
         if (building_frac_f%from_file) then
            if (building_frac_f%pars_xy(ind_wall, j, i) /= building_frac_f%fill) &
               surf_usm%frac(m, ind_veg_wall) = building_frac_f%pars_xy(ind_wall, j, i)
            if (building_frac_f%pars_xy(ind_green, j, i) /= building_frac_f%fill) &
               surf_usm%frac(m, ind_pav_green) = building_frac_f%pars_xy(ind_green, j, i)
            if (building_frac_f%pars_xy(ind_win, j, i) /= building_frac_f%fill) &
               surf_usm%frac(m, ind_wat_win) = building_frac_f%pars_xy(ind_win, j, i)
         end if
!
!--    Initialization green roof type (intensive or extensive).
         if (building_gen_f%from_file) then
            if (building_gen_f%pars_xy(ind_green_type_roof, j, i) /= building_gen_f%fill) &
               surf_usm%green_type_roof(m) = building_gen_f%pars_xy(ind_green_type_roof, j, i)
         end if
!
!--    Intialize LAI.
         if (building_lai_f%from_file) then
            if (building_lai_f%pars_xy(ind_level, j, i) /= building_lai_f%fill) &
               surf_usm%lai(m) = building_lai_f%pars_xy(ind_level, j, i)
         end if
!
!--    Roughness lengths.
         if (building_z0_f%from_file) then
            if (building_z0_f%pars_xy(ind_level, j, i) /= building_z0_f%fill) &
               surf_usm%z0(m) = building_z0_f%pars_xy(ind_level, j, i)
         end if
!
!--    Roughness lengths for moisture and heat.
         if (building_z0qh_f%from_file) then
            if (building_z0qh_f%pars_xy(ind_level, j, i) /= building_z0qh_f%fill) &
               surf_usm%z0h(m) = building_z0qh_f%pars_xy(ind_level, j, i)
            if (building_z0qh_f%pars_xy(ind_level, j, i) /= building_z0qh_f%fill) &
               surf_usm%z0q(m) = building_z0qh_f%pars_xy(ind_level, j, i)
         end if
!
!--    Initialize wall layer thicknesses.
         if (building_thick_f%from_file) then
!
!--       For first layer, depth = thickness.
            if (building_thick_f%pars_xyz(ind_wall, 0, j, i) /= building_thick_f%fill) &
               surf_usm%zw(nzb_wall, m) = building_thick_f%pars_xyz(ind_wall, 0, j, i)
            if (building_thick_f%pars_xyz(ind_win, 0, j, i) /= building_thick_f%fill) &
               surf_usm%zw_window(nzb_wall, m) = building_thick_f%pars_xyz(ind_win, 0, j, i)
            if (building_thick_f%pars_xyz(ind_green, 0, j, i) /= building_thick_f%fill) &
               surf_usm%zw_green(nzb_wall, m) = building_thick_f%pars_xyz(ind_green, 0, j, i)
!
!--       For other layers, depth difference = thickness.
            do ilayer = 1, nzw - 1
               if (building_thick_f%pars_xyz(ind_wall, ilayer, j, i) /= building_thick_f%fill) &
                  surf_usm%zw(nzb_wall + ilayer, m) = surf_usm%zw(nzb_wall + ilayer - 1, m) + &
                                                      building_thick_f%pars_xyz(ind_wall, ilayer, j, i)
               if (building_thick_f%pars_xyz(ind_win, ilayer, j, i) /= building_thick_f%fill) &
                  surf_usm%zw_window(nzb_wall + ilayer, m) = surf_usm%zw_window(nzb_wall + ilayer - 1, m) + &
                                                             building_thick_f%pars_xyz(ind_win, ilayer, j, i)
               if (building_thick_f%pars_xyz(ind_green, ilayer, j, i) /= building_thick_f%fill) &
                  surf_usm%zw_green(nzb_wall + ilayer, m) = surf_usm%zw_green(nzb_wall + ilayer - 1, m) + &
                                                            building_thick_f%pars_xyz(ind_green, ilayer, j, i)
            end do
         end if
!
!--    Intialize heat capacities form buildings_hcap list.
         if (building_hcap_f%from_file) then
            do ilayer = 0, nzw - 1
               if (building_hcap_f%pars_xyz(ind_wall, ilayer, j, i) /= building_hcap_f%fill) &
                  surf_usm%rho_c_wall(nzb_wall + ilayer, m) = &
                  building_hcap_f%pars_xyz(ind_wall, ilayer, j, i)
               if (building_hcap_f%pars_xyz(ind_win, ilayer, j, i) /= building_hcap_f%fill) &
                  surf_usm%rho_c_window(nzb_wall + ilayer, m) = &
                  building_hcap_f%pars_xyz(ind_win, ilayer, j, i)
               if (building_hcap_f%pars_xyz(ind_green, ilayer, j, i) /= building_hcap_f%fill) &
                  surf_usm%rho_c_green(nzb_wall + ilayer, m) = &
                  building_hcap_f%pars_xyz(ind_green, ilayer, j, i)
            end do
         end if
!
!--    Intialize heat conductivities.
         if (building_hcond_f%from_file) then
            do ilayer = 0, nzw - 1
               if (building_hcond_f%pars_xyz(ind_wall, ilayer, j, i) /= building_hcond_f%fill) &
                  surf_usm%lambda_h(nzb_wall + ilayer, m) = &
                  building_hcond_f%pars_xyz(ind_wall, ilayer, j, i)
               if (building_hcond_f%pars_xyz(ind_win, ilayer, j, i) /= building_hcond_f%fill) &
                  surf_usm%lambda_h_window(nzb_wall + ilayer, m) = &
                  building_hcond_f%pars_xyz(ind_win, ilayer, j, i)
               if (building_hcond_f%pars_xyz(ind_green, ilayer, j, i) /= building_hcond_f%fill) &
                  surf_usm%lambda_h_green(nzb_wall + ilayer, m) = &
                  building_hcond_f%pars_xyz(ind_green, ilayer, j, i)
            end do
         end if
!
!--    Intialize emissivities.
         if (building_emis_f%from_file) then
            if (building_emis_f%pars_xy(ind_wall, j, i) /= building_emis_f%fill) &
               surf_usm%emissivity(m, ind_veg_wall) = building_emis_f%pars_xy(ind_wall, j, i)
            if (building_emis_f%pars_xy(ind_win, j, i) /= building_emis_f%fill) &
               surf_usm%emissivity(m, ind_wat_win) = building_emis_f%pars_xy(ind_win, j, i)
            if (building_emis_f%pars_xy(ind_green, j, i) /= building_emis_f%fill) &
               surf_usm%emissivity(m, ind_pav_green) = building_emis_f%pars_xy(ind_green, j, i)
         end if
!
!--    Initialize albedo type.
         if (building_alb_type_f%from_file) then
            if (building_alb_type_f%pars_xy(ind_wall, j, i) /= building_alb_type_f%fill) &
               surf_usm%albedo_type(m, ind_veg_wall) = &
               int(building_alb_type_f%pars_xy(ind_wall, j, i))
            if (building_alb_type_f%pars_xy(ind_win, j, i) /= building_alb_type_f%fill) &
               surf_usm%albedo_type(m, ind_wat_win) = &
               int(building_alb_type_f%pars_xy(ind_win, j, i))
            if (building_alb_type_f%pars_xy(ind_green, j, i) /= building_alb_type_f%fill) &
               surf_usm%albedo_type(m, ind_pav_green) = &
               int(building_alb_type_f%pars_xy(ind_green, j, i))
         end if

!
!--    Initialize transmissivity of windows.
         if (building_trans_f%from_file) then
            if (building_trans_f%pars_xy(ind_level, j, i) /= building_trans_f%fill) &
               surf_usm%transmissivity(m) = building_trans_f%pars_xy(ind_level, j, i)
         end if

      end do
!
!-- Read building surface pars. If present, they override LOD1-LOD3 building pars where applicable.
!-- Same as LOD1-LOD3 initialization, roofs and walls are initialized separately.
!-- This initialization method is not realized for cut-cell topography at the moment.
      if (building_surface_pars_f%from_file .and. .not. cut_cell_topography) then
         do m = 1, surf_usm%ns
            i = surf_usm%i(m)
            j = surf_usm%j(m)
            k = surf_usm%k(m)
!
!--       Iterate over surfaces in column, check height and orientation
            do is = building_surface_pars_f%index_ji(1, j, i), building_surface_pars_f%index_ji(2, j, i)
               if (building_surface_pars_f%coords(4, is) == -surf_usm%koff(m) .and. &
                   building_surface_pars_f%coords(1, is) == k) &
                  then

                  if (building_surface_pars_f%pars(ind_s_wall_frac, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_veg_wall) = &
                     building_surface_pars_f%pars(ind_s_wall_frac, is)

                  if (building_surface_pars_f%pars(ind_s_green_frac_w, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_pav_green) = &
                     building_surface_pars_f%pars(ind_s_green_frac_w, is)

                  if (building_surface_pars_f%pars(ind_s_green_frac_r, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_pav_green) = &
                     building_surface_pars_f%pars(ind_s_green_frac_r, is)
                  !TODO clarify: why should _w and _r be on the same surface?

                  if (building_surface_pars_f%pars(ind_s_win_frac, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_wat_win) = building_surface_pars_f%pars(ind_s_win_frac, is)

                  if (building_surface_pars_f%pars(ind_s_lai_r, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%lai(m) = building_surface_pars_f%pars(ind_s_lai_r, is)

                  if (building_surface_pars_f%pars(ind_s_hc1, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%rho_c_wall(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_hc1, is)
                     surf_usm%rho_c_green(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_hc1, is)
                     surf_usm%rho_c_window(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_hc1, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_hc2, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%rho_c_wall(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_hc2, is)
                     surf_usm%rho_c_green(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_hc2, is)
                     surf_usm%rho_c_window(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_hc2, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_hc3, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%rho_c_wall(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_hc3, is)
                     surf_usm%rho_c_green(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_hc3, is)
                     surf_usm%rho_c_window(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_hc3, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_tc1, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%lambda_h(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_tc1, is)
                     surf_usm%lambda_h_green(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_tc1, is)
                     surf_usm%lambda_h_window(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_tc1, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_tc2, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%lambda_h(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_tc2, is)
                     surf_usm%lambda_h_green(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_tc2, is)
                     surf_usm%lambda_h_window(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_tc2, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_tc3, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%lambda_h(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_tc3, is)
                     surf_usm%lambda_h_green(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_tc3, is)
                     surf_usm%lambda_h_window(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_tc3, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_indoor_target_temp_summer, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%target_temp_summer(m) = &
                     building_surface_pars_f%pars(ind_s_indoor_target_temp_summer, is)

                  if (building_surface_pars_f%pars(ind_s_indoor_target_temp_winter, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%target_temp_winter(m) = &
                     building_surface_pars_f%pars(ind_s_indoor_target_temp_winter, is)

                  if (building_surface_pars_f%pars(ind_s_emis_wall, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%emissivity(m, ind_veg_wall) = &
                     building_surface_pars_f%pars(ind_s_emis_wall, is)

                  if (building_surface_pars_f%pars(ind_s_emis_green, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%emissivity(m, ind_pav_green) = &
                     building_surface_pars_f%pars(ind_s_emis_green, is)

                  if (building_surface_pars_f%pars(ind_s_emis_win, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%emissivity(m, ind_wat_win) = &
                     building_surface_pars_f%pars(ind_s_emis_win, is)

                  if (building_surface_pars_f%pars(ind_s_trans, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%transmissivity(m) = building_surface_pars_f%pars(ind_s_trans, is)

                  if (building_surface_pars_f%pars(ind_s_z0, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%z0(m) = building_surface_pars_f%pars(ind_s_z0, is)

                  if (building_surface_pars_f%pars(ind_s_z0qh, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%z0q(m) = building_surface_pars_f%pars(ind_s_z0qh, is)
                     surf_usm%z0h(m) = building_surface_pars_f%pars(ind_s_z0qh, is)
                  end if

                  exit ! Surface was found and processed
               end if
            end do
         end do

         do m = 1, surf_usm%ns
            i = surf_usm%i(m)
            j = surf_usm%j(m)
            k = surf_usm%k(m)
!
!--       Iterate over surfaces in column, check height and orientation
            do is = building_surface_pars_f%index_ji(1, j, i), building_surface_pars_f%index_ji(2, j, i)
               if (building_surface_pars_f%coords(5, is) == -surf_usm%joff(m) .and. &
                   building_surface_pars_f%coords(6, is) == -surf_usm%ioff(m) .and. &
                   building_surface_pars_f%coords(1, is) == k) &
                  then

                  if (building_surface_pars_f%pars(ind_s_wall_frac, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_veg_wall) = &
                     building_surface_pars_f%pars(ind_s_wall_frac, is)

                  if (building_surface_pars_f%pars(ind_s_green_frac_w, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_pav_green) = &
                     building_surface_pars_f%pars(ind_s_green_frac_w, is)

                  if (building_surface_pars_f%pars(ind_s_green_frac_r, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_pav_green) = &
                     building_surface_pars_f%pars(ind_s_green_frac_r, is)
                  !TODO Clarify: why should _w and _r be on the same surface?

                  if (building_surface_pars_f%pars(ind_s_win_frac, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%frac(m, ind_wat_win) = &
                     building_surface_pars_f%pars(ind_s_win_frac, is)

                  if (building_surface_pars_f%pars(ind_s_lai_r, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%lai(m) = building_surface_pars_f%pars(ind_s_lai_r, is)

                  if (building_surface_pars_f%pars(ind_s_hc1, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%rho_c_wall(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_hc1, is)
                     surf_usm%rho_c_green(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_hc1, is)
                     surf_usm%rho_c_window(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_hc1, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_hc2, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%rho_c_wall(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_hc2, is)
                     surf_usm%rho_c_green(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_hc2, is)
                     surf_usm%rho_c_window(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_hc2, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_hc3, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%rho_c_wall(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_hc3, is)
                     surf_usm%rho_c_green(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_hc3, is)
                     surf_usm%rho_c_window(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_hc3, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_tc1, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%lambda_h(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_tc1, is)
                     surf_usm%lambda_h_green(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_tc1, is)
                     surf_usm%lambda_h_window(nzb_wall:nzb_wall + 1, m) = &
                        building_surface_pars_f%pars(ind_s_tc1, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_tc2, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%lambda_h(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_tc2, is)
                     surf_usm%lambda_h_green(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_tc2, is)
                     surf_usm%lambda_h_window(nzb_wall + 2, m) = &
                        building_surface_pars_f%pars(ind_s_tc2, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_tc3, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%lambda_h(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_tc3, is)
                     surf_usm%lambda_h_green(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_tc3, is)
                     surf_usm%lambda_h_window(nzb_wall + 3, m) = &
                        building_surface_pars_f%pars(ind_s_tc3, is)
                  end if

                  if (building_surface_pars_f%pars(ind_s_indoor_target_temp_summer, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%target_temp_summer(m) = &
                     building_surface_pars_f%pars(ind_s_indoor_target_temp_summer, is)

                  if (building_surface_pars_f%pars(ind_s_indoor_target_temp_winter, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%target_temp_winter(m) = &
                     building_surface_pars_f%pars(ind_s_indoor_target_temp_winter, is)

                  if (building_surface_pars_f%pars(ind_s_emis_wall, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%emissivity(m, ind_veg_wall) = &
                     building_surface_pars_f%pars(ind_s_emis_wall, is)

                  if (building_surface_pars_f%pars(ind_s_emis_green, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%emissivity(m, ind_pav_green) = &
                     building_surface_pars_f%pars(ind_s_emis_green, is)

                  if (building_surface_pars_f%pars(ind_s_emis_win, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%emissivity(m, ind_wat_win) = &
                     building_surface_pars_f%pars(ind_s_emis_win, is)

                  if (building_surface_pars_f%pars(ind_s_trans, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%transmissivity(m) = &
                     building_surface_pars_f%pars(ind_s_trans, is)

                  if (building_surface_pars_f%pars(ind_s_z0, is) /= &
                      building_surface_pars_f%fill) &
                     surf_usm%z0(m) = building_surface_pars_f%pars(ind_s_z0, is)

                  if (building_surface_pars_f%pars(ind_s_z0qh, is) /= &
                      building_surface_pars_f%fill) &
                     then
                     surf_usm%z0q(m) = building_surface_pars_f%pars(ind_s_z0qh, is)
                     surf_usm%z0h(m) = building_surface_pars_f%pars(ind_s_z0qh, is)
                  end if

                  exit ! Surface was found and processed
               end if
            end do
         end do
      end if
!
!-- Check if material fractions of surfaces sum up to 1.
      relative_fraction_error = .false.
      do m = 1, surf_usm%ns
         if ((abs(sum(surf_usm%frac(m, :)) - 1.0_wp)) > 1.0e-5_wp) then
            relative_fraction_error = .true.
         end if
      end do

#if defined( __parallel )
      call MPI_ALLREDUCE(MPI_IN_PLACE, relative_fraction_error, 1, MPI_LOGICAL, MPI_LOR, comm2d, &
                         ierr)
#endif

      if (relative_fraction_error) then
         message_string = 'relative material fractions do not sum-up to one at some surfaces'
         call message('usm_init', 'USM0011', 2, 2, 0, 6, 0)
      end if

!
!-- Initialization of the wall/roof materials.
      call usm_init_wall_heat_model()

!
!-- Init moist green heat capacity with the respective dry value. This is needed for the calculation
!-- of dt_usm.
      surf_usm%rho_c_total_green = surf_usm%rho_c_green

!
!-- Init skin layer properties (can be done after initialization of wall layers).
      do m = 1, surf_usm%ns
         surf_usm%c_surface(m) = surf_usm%rho_c_wall(nzb_wall, m)* &
                                 surf_usm%dz_wall(nzb_wall, m)*0.25_wp
         surf_usm%lambda_surf(m) = surf_usm%lambda_h(nzb_wall, m)* &
                                   surf_usm%ddz_wall(nzb_wall, m)*2.0_wp
         surf_usm%c_surface_green(m) = surf_usm%rho_c_wall(nzb_wall, m)* &
                                       surf_usm%dz_wall(nzb_wall, m)*0.25_wp
         surf_usm%lambda_surf_green(m) = surf_usm%lambda_h_green(nzb_wall, m)* &
                                         surf_usm%ddz_green(nzb_wall, m)*2.0_wp
         surf_usm%c_surface_window(m) = surf_usm%rho_c_window(nzb_wall, m)* &
                                        surf_usm%dz_window(nzb_wall, m)*0.25_wp
         surf_usm%lambda_surf_window(m) = surf_usm%lambda_h_window(nzb_wall, m)* &
                                          surf_usm%ddz_window(nzb_wall, m)*2.0_wp
      end do

!
!-- Check for consistent initialization.
!-- Check if roughness length for momentum, heat, or moisture exceed surface-layer height and
!-- limit local roughness length where necessary (if allowed). If limited, give an informative
!-- message only once in order to avoid the job protocol to be messed-up with messages.
      flag_exceed_z0 = .false.
      flag_exceed_z0h = .false.
      flag_exceed_z0q = .false.
      do m = 1, surf_usm%ns

         if (surf_usm%z0(m) >= 0.5*surf_usm%z_mo(m)) then
            flag_exceed_z0 = .true.
            if (allow_roughness_limitation) then
               surf_usm%z0(m) = 0.5_wp*surf_usm%z_mo(m)
            else
               exit_index = m
               exit
            end if
         end if

         if (surf_usm%z0h(m) >= surf_usm%z_mo(m)) then
            flag_exceed_z0h = .true.
            if (allow_roughness_limitation) then
               surf_usm%z0h(m) = 0.5_wp*surf_usm%z_mo(m)
            else
               exit_index = m
               exit
            end if
         end if

         if (surf_usm%z0q(m) >= surf_usm%z_mo(m)) then
            flag_exceed_z0q = .true.
            if (allow_roughness_limitation) then
               surf_usm%z0q(m) = 0.5_wp*surf_usm%z_mo(m)
            else
               exit_index = m
               exit
            end if
         end if

      end do

      if (flag_exceed_z0 .and. .not. allow_roughness_limitation) then
         write (message_string, '(A,I6,A,I6,A)') &
            'z0 exceeds 0.5 * surface-layer height at building surface grid point (i,j) = (', &
            surf_usm%i(exit_index) + surf_usm%ioff(exit_index), ',', &
            surf_usm%j(exit_index) + surf_usm%joff(exit_index), ')'
         call message('usm_init', 'USM0012', 2, 2, myid, 6, 0)
      end if
#if defined( __parallel )
      call MPI_ALLREDUCE(MPI_IN_PLACE, flag_exceed_z0, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
      if (flag_exceed_z0) then
         write (message_string, *) 'z0 exceeds 0.5 * surface-layer height at building surface(s)'// &
            ' and is limited to that height'
         call message('usm_init', 'USM0013', 0, 0, 0, 6, 0)
      end if

      if (flag_exceed_z0h .and. .not. allow_roughness_limitation) then
         write (message_string, '(A,I6,A,I6,A)') &
            'z0h exceeds 0.5 * surface-layer height at building surface grid point (i,j) = (', &
            surf_usm%i(exit_index) + surf_usm%ioff(exit_index), ',', &
            surf_usm%j(exit_index) + surf_usm%joff(exit_index), ')'
         call message('usm_init', 'USM0012', 2, 2, myid, 6, 0)
      end if
#if defined( __parallel )
      call MPI_ALLREDUCE(MPI_IN_PLACE, flag_exceed_z0h, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
      if (flag_exceed_z0h) then
         write (message_string, *) 'z0h exceeds 0.5 * surface-layer height at building surface(s)'// &
            ' and is limited to that height'
         call message('usm_init', 'USM0013', 0, 0, 0, 6, 0)
      end if

      if (flag_exceed_z0q .and. .not. allow_roughness_limitation) then
         write (message_string, '(A,I6,A,I6,A)') &
            'z0q exceeds 0.5 * surface-layer height at building surface grid point (i,j) = (', &
            surf_usm%i(exit_index) + surf_usm%ioff(exit_index), ',', &
            surf_usm%j(exit_index) + surf_usm%joff(exit_index), ')'
         call message('usm_init', 'USM0012', 2, 2, myid, 6, 0)
      end if
#if defined( __parallel )
      call MPI_ALLREDUCE(MPI_IN_PLACE, flag_exceed_z0q, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
      if (flag_exceed_z0q) then
         write (message_string, *) 'z0q exceeds 0.5 * surface-layer height at building surface(s)'// &
            ' and is limited to that height'
         call message('usm_init', 'USM0013', 0, 0, 0, 6, 0)
      end if

!
!--  Intitialization of the surface and wall/ground/roof temperature. These actions must not be done
!--  in restart runs or when spinup data is read.
      if (trim(initializing_actions) /= 'read_restart_data' .and. &
          .not. read_spinup_data) then
         do m = 1, surf_usm%ns
            i = surf_usm%i(m)
            j = surf_usm%j(m)
            k = surf_usm%k(m)

            t_surf_wall%val(m) = pt(k, j, i)*exner(k)
            t_surf_window%val(m) = pt(k, j, i)*exner(k)
            t_surf_green%val(m) = pt(k, j, i)*exner(k)
            surf_usm%pt_surface(m) = pt(k, j, i)*exner(k)
         end do
!
!--      For the sake of correct initialization, set also q_surface.
!--      Note, at urban surfaces q_surface is initialized with 0.
         if (humidity) then
            do m = 1, surf_usm%ns
               surf_usm%q_surface(m) = 0.0_wp
            end do
         end if
!
!--      Initial values for t_wall.
!--      Outer value is set to surface temperature, inner value is set to wall_inner_temperature
!--      and profile is logaritmic (linear in nz).
!--      Again, initialization is separated between roofs and walls. Start with roofs.
         do m = 1, surf_usm%ns
            if (surf_usm%upward(m) .or. surf_usm%downward(m)) then
               tin = roof_inner_temperature
            else
               tin = wall_inner_temperature
            end if
            twin = window_inner_temperature

            do k = nzb_wall, nzt_wall + 1
               c = real(k - nzb_wall, wp)/real(nzt_wall + 1 - nzb_wall, wp)

               t_wall%val(k, m) = (1.0_wp - c)*t_surf_wall%val(m) + c*tin
               t_window%val(k, m) = (1.0_wp - c)*t_surf_window%val(m) + c*twin
               t_green%val(k, m) = t_surf_wall%val(m)
            end do
         end do
      end if

!
!--  Possibly DO user-defined actions (e.g. define heterogeneous wall surface)
      call user_init_urban_surface

!
!--  Initialize prognostic values for the first timestep
      t_surf_wall_p = t_surf_wall
      t_surf_window_p = t_surf_window
      t_surf_green_p = t_surf_green

      swc_p = swc
      t_wall_p = t_wall
      t_window_p = t_window
      t_green_p = t_green

!
!-- Set initial values for prognostic soil quantities
      if (trim(initializing_actions) /= 'read_restart_data' .and. .not. read_spinup_data) then
         m_liq_usm%val = 0.0_wp
      end if
      m_liq_usm_p%val = m_liq_usm%val
!
!-- Set initial values for diagnostic quantities
      surf_usm%c_liq = 0.0_wp
      surf_usm%qsws_liq = 0.0_wp
      surf_usm%qsws_veg = 0.0_wp

      call cpu_log(log_point_s(78), 'usm_init', 'stop')

      if (debug_output) call debug_message('usm_init', 'end')

   contains

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Initialization of wall and surface properties with building_pars list from PIDS (deprecated).
!--------------------------------------------------------------------------------------------------!
      subroutine usm_init_building_pars_pids

!
!-- Indices of input attributes in building_pars for (above) ground floor level
         integer(iwp) ::  ind_alb_wall_agfl = 38   !< index in input list for albedo_type of wall above ground floor level
         integer(iwp) ::  ind_alb_green_agfl = 39   !< index in input list for albedo_type of green above ground floor level
         integer(iwp) ::  ind_alb_win_agfl = 40   !< index in input list for albedo_type of window fraction above ground floor level
         integer(iwp) ::  ind_emis_wall_agfl = 14   !< index in input list for wall emissivity, above ground floor level
         integer(iwp) ::  ind_emis_wall_gfl = 32   !< index in input list for wall emissivity, ground floor level
         integer(iwp) ::  ind_emis_green_agfl = 15   !< index in input list for green emissivity, above ground floor level
         integer(iwp) ::  ind_emis_green_gfl = 34   !< index in input list for green emissivity, ground floor level
         integer(iwp) ::  ind_emis_win_agfl = 16   !< index in input list for window emissivity, above ground floor level
         integer(iwp) ::  ind_emis_win_gfl = 33   !< index in input list for window emissivity, ground floor level
         integer(iwp) ::  ind_green_frac_w_agfl = 2    !< index in input list for green fraction on wall, above ground floor level
         integer(iwp) ::  ind_green_frac_w_gfl = 23   !< index in input list for green fraction on wall, ground floor level
         integer(iwp) ::  ind_green_frac_r_agfl = 3    !< index in input list for green fraction on roof, above ground floor level
         integer(iwp) ::  ind_hc1_agfl = 6    !< index in input list for heat capacity at first wall layer, above ground floor level
         integer(iwp) ::  ind_hc1_gfl = 26   !< index in input list for heat capacity at first wall layer, ground floor level
         integer(iwp) ::  ind_hc2_agfl = 7
         !< index in input list for heat capacity at second wall layer, above ground floor level
         integer(iwp) ::  ind_hc2_gfl = 27   !< index in input list for heat capacity at second wall layer, ground floor level
         integer(iwp) ::  ind_hc3_agfl = 8    !< index in input list for heat capacity at third wall layer, above ground floor level
         integer(iwp) ::  ind_hc3_gfl = 28   !< index in input list for heat capacity at third wall layer, ground floor level
         integer(iwp) ::  ind_hc4_agfl = 136
         !< index in input list for heat capacity at fourth wall layer, above ground floor level
         integer(iwp) ::  ind_hc4_gfl = 138  !< index in input list for heat capacity at fourth wall layer, ground floor level
         integer(iwp) ::  ind_indoor_target_temp_summer = 12  !<
         integer(iwp) ::  ind_indoor_target_temp_winter = 13  !<
         integer(iwp) ::  ind_lai_r_agfl = 4    !< index in input list for LAI on roof, above ground floor level
         integer(iwp) ::  ind_lai_w_agfl = 5    !< index in input list for LAI on wall, above ground floor level
         integer(iwp) ::  ind_lai_w_gfl = 25   !< index in input list for LAI on wall, ground floor level
         integer(iwp) ::  ind_tc1_agfl = 9
         !< index in input list for thermal conductivity at first wall layer, above ground floor level
         integer(iwp) ::  ind_tc1_gfl = 29   !< index in input list for thermal conductivity at first wall layer, ground floor level
         integer(iwp) ::  ind_tc2_agfl = 10
         !< index in input list for thermal conductivity at second wall layer, above ground floor level
         integer(iwp) ::  ind_tc2_gfl = 30
         !< index in input list for thermal conductivity at second wall layer, ground floor level
         integer(iwp) ::  ind_tc3_agfl = 11
         !< index in input list for thermal conductivity at third wall layer, above ground floor level
         integer(iwp) ::  ind_tc3_gfl = 31   !< index in input list for thermal conductivity at third wall layer, ground floor level
         integer(iwp) ::  ind_tc4_agfl = 137
         !< index in input list for thermal conductivity at fourth wall layer, above ground floor level
         integer(iwp) ::  ind_tc4_gfl = 139
         !< index in input list for thermal conductivity at fourth wall layer, ground floor level
         integer(iwp) ::  ind_thick_1_agfl = 41   !< index for wall layer thickness - 1st layer above ground floor level
         integer(iwp) ::  ind_thick_2_agfl = 42   !< index for wall layer thickness - 2nd layer above ground floor level
         integer(iwp) ::  ind_thick_3_agfl = 43   !< index for wall layer thickness - 3rd layer above ground floor level
         integer(iwp) ::  ind_thick_4_agfl = 44   !< index for wall layer thickness - 4th layer above ground floor level
         integer(iwp) ::  ind_trans_agfl = 17   !< index in input list for window transmissivity, above ground floor level
         integer(iwp) ::  ind_trans_gfl = 35   !< index in input list for window transmissivity, ground floor level
         integer(iwp) ::  ind_wall_frac_agfl = 0    !< index in input list for wall fraction, above ground floor level
         integer(iwp) ::  ind_wall_frac_gfl = 21   !< index in input list for wall fraction, ground floor level
         integer(iwp) ::  ind_wall_frac_r = 89   !< index in input list for wall fraction, roof
         integer(iwp) ::  ind_win_frac_agfl = 1    !< index in input list for window fraction, above ground floor level
         integer(iwp) ::  ind_win_frac_gfl = 22   !< index in input list for window fraction, ground floor level
         integer(iwp) ::  ind_win_frac_r = 102  !< index in input list for window fraction, roof
         integer(iwp) ::  ind_z0_agfl = 18   !< index in input list for z0, above ground floor level
         integer(iwp) ::  ind_z0_gfl = 36   !< index in input list for z0, ground floor level
         integer(iwp) ::  ind_z0qh_agfl = 19   !< index in input list for z0h / z0q, above ground floor level
         integer(iwp) ::  ind_z0qh_gfl = 37   !< index in input list for z0h / z0q, ground floor level

         integer(iwp) ::  ind_wall_frac  !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_green_frac !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_win_frac   !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_lai        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_z0         !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_z0qh       !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_hc1        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_hc2        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_hc3        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_hc4        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_tc1        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_tc2        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_tc3        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_tc4        !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_emis_wall  !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_emis_green !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_emis_win   !< corresponding merged index from gfl, agfl and roof
         integer(iwp) ::  ind_trans      !< corresponding merged index from gfl, agfl and roof

         if (building_pars_f%from_file) then
            do m = 1, surf_usm%ns
               if (surf_usm%upward(m) .or. surf_usm%downward(m) .or. surf_usm%cut_cell_roof(m)) &
                  then
                  i = surf_usm%i(m)
                  j = surf_usm%j(m)

                  ind_wall_frac = ind_wall_frac_r
                  ind_green_frac = ind_green_frac_r_agfl
                  ind_win_frac = ind_win_frac_r
                  ind_lai = ind_lai_r_agfl
                  ind_z0 = ind_z0_agfl
                  ind_z0qh = ind_z0qh_agfl
                  ind_hc1 = ind_hc1_agfl
                  ind_hc2 = ind_hc2_agfl
                  ind_hc3 = ind_hc3_agfl
                  ind_hc4 = ind_hc4_agfl
                  ind_tc1 = ind_tc1_agfl
                  ind_tc2 = ind_tc2_agfl
                  ind_tc3 = ind_tc3_agfl
                  ind_tc4 = ind_tc4_agfl
                  ind_emis_wall = ind_emis_wall_agfl
                  ind_emis_green = ind_emis_green_agfl
                  ind_emis_win = ind_emis_win_agfl
                  ind_trans = ind_trans_agfl
!
!--          Initialize relatvie wall- (0), green- (1) and window (2) fractions
                  if (building_pars_f%pars_xy(ind_wall_frac, j, i) /= building_pars_f%fill) &
                     surf_usm%frac(m, ind_veg_wall) = building_pars_f%pars_xy(ind_wall_frac, j, i)

                  if (building_pars_f%pars_xy(ind_green_frac, j, i) /= building_pars_f%fill) &
                     surf_usm%frac(m, ind_pav_green) = building_pars_f%pars_xy(ind_green_frac, j, i)

                  if (building_pars_f%pars_xy(ind_win_frac, j, i) /= building_pars_f%fill) &
                     surf_usm%frac(m, ind_wat_win) = building_pars_f%pars_xy(ind_win_frac, j, i)

                  if (building_pars_f%pars_xy(ind_lai, j, i) /= building_pars_f%fill) &
                     surf_usm%lai(m) = building_pars_f%pars_xy(ind_lai, j, i)

                  if (building_pars_f%pars_xy(ind_hc1, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall, m) = building_pars_f%pars_xy(ind_hc1, j, i)

                  if (building_pars_f%pars_xy(ind_hc2, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_hc2, j, i)

                  if (building_pars_f%pars_xy(ind_hc3, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_hc3, j, i)

                  if (building_pars_f%pars_xy(ind_hc4, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_hc4, j, i)

                  if (building_pars_f%pars_xy(ind_hc1, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall, m) = building_pars_f%pars_xy(ind_hc1, j, i)

                  if (building_pars_f%pars_xy(ind_hc2, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_hc2, j, i)

                  if (building_pars_f%pars_xy(ind_hc3, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_hc3, j, i)

                  if (building_pars_f%pars_xy(ind_hc4, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_hc4, j, i)

                  if (building_pars_f%pars_xy(ind_hc1, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall, m) = building_pars_f%pars_xy(ind_hc1, j, i)

                  if (building_pars_f%pars_xy(ind_hc2, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_hc2, j, i)

                  if (building_pars_f%pars_xy(ind_hc3, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_hc3, j, i)

                  if (building_pars_f%pars_xy(ind_hc4, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_hc4, j, i)

                  if (building_pars_f%pars_xy(ind_tc1, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall, m) = building_pars_f%pars_xy(ind_tc1, j, i)

                  if (building_pars_f%pars_xy(ind_tc2, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_tc2, j, i)

                  if (building_pars_f%pars_xy(ind_tc3, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_tc3, j, i)

                  if (building_pars_f%pars_xy(ind_tc4, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_tc4, j, i)

                  if (building_pars_f%pars_xy(ind_tc1, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall, m) = building_pars_f%pars_xy(ind_tc1, j, i)

                  if (building_pars_f%pars_xy(ind_tc2, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_tc2, j, i)

                  if (building_pars_f%pars_xy(ind_tc3, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_tc3, j, i)

                  if (building_pars_f%pars_xy(ind_tc4, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_tc4, j, i)

                  if (building_pars_f%pars_xy(ind_tc1, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall, m) = building_pars_f%pars_xy(ind_tc1, j, i)

                  if (building_pars_f%pars_xy(ind_tc2, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_tc2, j, i)

                  if (building_pars_f%pars_xy(ind_tc3, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_tc3, j, i)

                  if (building_pars_f%pars_xy(ind_tc4, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_tc4, j, i)

                  if (building_pars_f%pars_xy(ind_indoor_target_temp_summer, j, i) /= &
                      building_pars_f%fill) &
                     surf_usm%target_temp_summer(m) = building_pars_f%pars_xy(ind_indoor_target_temp_summer, j, i)

                  if (building_pars_f%pars_xy(ind_indoor_target_temp_winter, j, i) /= &
                      building_pars_f%fill) &
                     surf_usm%target_temp_winter(m) = building_pars_f%pars_xy(ind_indoor_target_temp_winter, j, i)

                  if (building_pars_f%pars_xy(ind_emis_wall, j, i) /= building_pars_f%fill) &
                     surf_usm%emissivity(m, ind_veg_wall) = building_pars_f%pars_xy(ind_emis_wall, j, i)

                  if (building_pars_f%pars_xy(ind_emis_green, j, i) /= building_pars_f%fill) &
                     surf_usm%emissivity(m, ind_pav_green) = building_pars_f%pars_xy(ind_emis_green, j, i)

                  if (building_pars_f%pars_xy(ind_emis_win, j, i) /= building_pars_f%fill) &
                     surf_usm%emissivity(m, ind_wat_win) = building_pars_f%pars_xy(ind_emis_win, j, i)

                  if (building_pars_f%pars_xy(ind_trans, j, i) /= building_pars_f%fill) &
                     surf_usm%transmissivity(m) = building_pars_f%pars_xy(ind_trans, j, i)

                  if (building_pars_f%pars_xy(ind_z0, j, i) /= building_pars_f%fill) &
                     surf_usm%z0(m) = building_pars_f%pars_xy(ind_z0, j, i)

                  if (building_pars_f%pars_xy(ind_z0qh, j, i) /= building_pars_f%fill) &
                     surf_usm%z0h(m) = building_pars_f%pars_xy(ind_z0qh, j, i)

                  if (building_pars_f%pars_xy(ind_z0qh, j, i) /= building_pars_f%fill) &
                     surf_usm%z0q(m) = building_pars_f%pars_xy(ind_z0qh, j, i)

                  if (building_pars_f%pars_xy(ind_alb_wall_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%albedo_type(m, ind_veg_wall) = building_pars_f%pars_xy(ind_alb_wall_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_alb_green_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%albedo_type(m, ind_pav_green) = building_pars_f%pars_xy(ind_alb_green_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_alb_win_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%albedo_type(m, ind_wat_win) = building_pars_f%pars_xy(ind_alb_win_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_1_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall, m) = building_pars_f%pars_xy(ind_thick_1_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_2_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_thick_2_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_3_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_thick_3_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_4_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_thick_4_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_1_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall, m) = building_pars_f%pars_xy(ind_thick_1_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_2_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_thick_2_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_3_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_thick_3_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_4_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_thick_4_agfl, j, i)
               end if
            end do

            do m = 1, surf_usm%ns
               if (.not. (surf_usm%upward(m) .or. surf_usm%downward(m)) .or. &
                   surf_usm%cut_cell_wall(m)) then
                  i = surf_usm%i(m) + surf_usm%ioff(m)
                  j = surf_usm%j(m) + surf_usm%joff(m)
!
!--          In order to distinguish between ground floor level and above-ground-floor level
!--          surfaces, set input indices.
                  ind_wall_frac = merge(ind_wall_frac_gfl, ind_wall_frac_agfl, surf_usm%gfl(m))
                  ind_green_frac_w = merge(ind_green_frac_w_gfl, ind_green_frac_w_agfl, surf_usm%gfl(m))
                  ind_win_frac = merge(ind_win_frac_gfl, ind_win_frac_agfl, surf_usm%gfl(m))
                  ind_lai_w = merge(ind_lai_w_gfl, ind_lai_w_agfl, surf_usm%gfl(m))
                  ind_z0 = merge(ind_z0_gfl, ind_z0_agfl, surf_usm%gfl(m))
                  ind_z0qh = merge(ind_z0qh_gfl, ind_z0qh_agfl, surf_usm%gfl(m))
                  ind_hc1 = merge(ind_hc1_gfl, ind_hc1_agfl, surf_usm%gfl(m))
                  ind_hc2 = merge(ind_hc2_gfl, ind_hc2_agfl, surf_usm%gfl(m))
                  ind_hc3 = merge(ind_hc3_gfl, ind_hc3_agfl, surf_usm%gfl(m))
                  ind_hc4 = merge(ind_hc4_gfl, ind_hc4_agfl, surf_usm%gfl(m))
                  ind_tc1 = merge(ind_tc1_gfl, ind_tc1_agfl, surf_usm%gfl(m))
                  ind_tc2 = merge(ind_tc2_gfl, ind_tc2_agfl, surf_usm%gfl(m))
                  ind_tc3 = merge(ind_tc3_gfl, ind_tc3_agfl, surf_usm%gfl(m))
                  ind_tc4 = merge(ind_tc4_gfl, ind_tc4_agfl, surf_usm%gfl(m))
                  ind_emis_wall = merge(ind_emis_wall_gfl, ind_emis_wall_agfl, surf_usm%gfl(m))
                  ind_emis_green = merge(ind_emis_green_gfl, ind_emis_green_agfl, surf_usm%gfl(m))
                  ind_emis_win = merge(ind_emis_win_gfl, ind_emis_win_agfl, surf_usm%gfl(m))
                  ind_trans = merge(ind_trans_gfl, ind_trans_agfl, surf_usm%gfl(m))
!
!--          Initialize relatvie wall- (0), green- (1) and window (2) fractions
                  if (building_pars_f%pars_xy(ind_wall_frac, j, i) /= building_pars_f%fill) &
                     surf_usm%frac(m, ind_veg_wall) = building_pars_f%pars_xy(ind_wall_frac, j, i)

                  if (building_pars_f%pars_xy(ind_green_frac_w, j, i) /= building_pars_f%fill) &
                     surf_usm%frac(m, ind_pav_green) = building_pars_f%pars_xy(ind_green_frac_w, j, i)

                  if (building_pars_f%pars_xy(ind_win_frac, j, i) /= building_pars_f%fill) &
                     surf_usm%frac(m, ind_wat_win) = building_pars_f%pars_xy(ind_win_frac, j, i)

                  if (building_pars_f%pars_xy(ind_lai_w, j, i) /= building_pars_f%fill) &
                     surf_usm%lai(m) = building_pars_f%pars_xy(ind_lai_w, j, i)

                  if (building_pars_f%pars_xy(ind_hc1, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall, m) = building_pars_f%pars_xy(ind_hc1, j, i)

                  if (building_pars_f%pars_xy(ind_hc2, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_hc2, j, i)

                  if (building_pars_f%pars_xy(ind_hc3, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_hc3, j, i)

                  if (building_pars_f%pars_xy(ind_hc4, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_wall(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_hc4, j, i)

                  if (building_pars_f%pars_xy(ind_hc1, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall, m) = building_pars_f%pars_xy(ind_hc1, j, i)

                  if (building_pars_f%pars_xy(ind_hc2, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_hc2, j, i)

                  if (building_pars_f%pars_xy(ind_hc3, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_hc3, j, i)

                  if (building_pars_f%pars_xy(ind_hc4, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_green(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_hc4, j, i)

                  if (building_pars_f%pars_xy(ind_hc1, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall, m) = building_pars_f%pars_xy(ind_hc1, j, i)

                  if (building_pars_f%pars_xy(ind_hc2, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_hc2, j, i)

                  if (building_pars_f%pars_xy(ind_hc3, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_hc3, j, i)

                  if (building_pars_f%pars_xy(ind_hc4, j, i) /= building_pars_f%fill) &
                     surf_usm%rho_c_window(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_hc4, j, i)

                  if (building_pars_f%pars_xy(ind_tc1, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall, m) = building_pars_f%pars_xy(ind_tc1, j, i)

                  if (building_pars_f%pars_xy(ind_tc2, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_tc2, j, i)

                  if (building_pars_f%pars_xy(ind_tc3, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_tc3, j, i)

                  if (building_pars_f%pars_xy(ind_tc4, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_tc4, j, i)

                  if (building_pars_f%pars_xy(ind_tc1, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall, m) = building_pars_f%pars_xy(ind_tc1, j, i)

                  if (building_pars_f%pars_xy(ind_tc2, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_tc2, j, i)

                  if (building_pars_f%pars_xy(ind_tc3, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_tc3, j, i)

                  if (building_pars_f%pars_xy(ind_tc4, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_green(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_tc4, j, i)

                  if (building_pars_f%pars_xy(ind_tc1, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall, m) = building_pars_f%pars_xy(ind_tc1, j, i)

                  if (building_pars_f%pars_xy(ind_tc2, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_tc2, j, i)

                  if (building_pars_f%pars_xy(ind_tc3, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_tc3, j, i)

                  if (building_pars_f%pars_xy(ind_tc4, j, i) /= building_pars_f%fill) &
                     surf_usm%lambda_h_window(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_tc4, j, i)

                  if (building_pars_f%pars_xy(ind_indoor_target_temp_summer, j, i) /= building_pars_f%fill) &
                     surf_usm%target_temp_summer(m) = building_pars_f%pars_xy(ind_indoor_target_temp_summer, j, i)

                  if (building_pars_f%pars_xy(ind_indoor_target_temp_winter, j, i) /= building_pars_f%fill) &
                     surf_usm%target_temp_winter(m) = building_pars_f%pars_xy(ind_indoor_target_temp_winter, j, i)

                  if (building_pars_f%pars_xy(ind_emis_wall, j, i) /= building_pars_f%fill) &
                     surf_usm%emissivity(m, ind_veg_wall) = building_pars_f%pars_xy(ind_emis_wall, j, i)

                  if (building_pars_f%pars_xy(ind_emis_green, j, i) /= building_pars_f%fill) &
                     surf_usm%emissivity(m, ind_pav_green) = building_pars_f%pars_xy(ind_emis_green, j, i)

                  if (building_pars_f%pars_xy(ind_emis_win, j, i) /= building_pars_f%fill) &
                     surf_usm%emissivity(m, ind_wat_win) = building_pars_f%pars_xy(ind_emis_win, j, i)

                  if (building_pars_f%pars_xy(ind_trans, j, i) /= building_pars_f%fill) &
                     surf_usm%transmissivity(m) = building_pars_f%pars_xy(ind_trans, j, i)

                  if (building_pars_f%pars_xy(ind_z0, j, i) /= building_pars_f%fill) &
                     surf_usm%z0(m) = building_pars_f%pars_xy(ind_z0, j, i)

                  if (building_pars_f%pars_xy(ind_z0qh, j, i) /= building_pars_f%fill) &
                     surf_usm%z0h(m) = building_pars_f%pars_xy(ind_z0qh, j, i)

                  if (building_pars_f%pars_xy(ind_z0qh, j, i) /= building_pars_f%fill) &
                     surf_usm%z0q(m) = building_pars_f%pars_xy(ind_z0qh, j, i)

                  if (building_pars_f%pars_xy(ind_alb_wall_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%albedo_type(m, ind_veg_wall) = building_pars_f%pars_xy(ind_alb_wall_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_alb_green_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%albedo_type(m, ind_pav_green) = building_pars_f%pars_xy(ind_alb_green_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_alb_win_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%albedo_type(m, ind_wat_win) = building_pars_f%pars_xy(ind_alb_win_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_1_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall, m) = building_pars_f%pars_xy(ind_thick_1_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_2_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_thick_2_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_3_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_thick_3_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_4_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_thick_4_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_1_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall, m) = building_pars_f%pars_xy(ind_thick_1_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_2_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall + 1, m) = building_pars_f%pars_xy(ind_thick_2_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_3_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall + 2, m) = building_pars_f%pars_xy(ind_thick_3_agfl, j, i)

                  if (building_pars_f%pars_xy(ind_thick_4_agfl, j, i) /= building_pars_f%fill) &
                     surf_usm%zw_green(nzb_wall + 3, m) = building_pars_f%pars_xy(ind_thick_4_agfl, j, i)

               end if
            end do
         end if

      end subroutine usm_init_building_pars_pids

   end subroutine usm_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Wall model as part of the urban surface model. The model predicts vertical and horizontal
!> wall / roof temperatures and window layer temperatures. No window layer temperature calculactions
!> during spinup to increase possible timestep.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_wall_heat_model

      implicit none

      logical             ::  runge_l        !< dummy flag to indicate RK timestepping scheme

      integer(iwp) ::  kw !< grid index - wall depth
      integer(iwp) ::  m  !< running index for surface elements

      real(wp) ::  win_absorp        !< absorption coefficient from transmissivity
      real(wp) ::  win_nonrefl_1side !< non-reflected fraction after outer glass boundary

      real(wp), dimension(1:surf_usm%ns) ::  fac_veg           !< pre-calcualated factor for vegetation fraction
      real(wp), dimension(1:surf_usm%ns) ::  fac_wall          !< pre-calcualated factor for wall fraction

      real(wp), dimension(nzb_wall:nzt_wall, 1:surf_usm%ns) ::  wtend   !< computed tendency for wall fraction
      real(wp), dimension(nzb_wall:nzt_wall, 1:surf_usm%ns) ::  wintend !< computed tendency for window fraction

      type(surf_type), pointer ::  surf !< surface-date type variable

      runge_l = (timestep_scheme(1:5) == 'runge')

      if (debug_output_timestep) then
         write (debug_string, *) 'usm_wall_heat_model: '
         call debug_message(debug_string, 'start')
      end if

      surf => surf_usm
!
!-- Prognostic equation for ground/roof temperature t_wall.
      wtend = 0.0_wp

!
!-- Compute fractions.
      fac_wall = 0.0_wp
      fac_veg = 0.0_wp
!$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
      do m = 1, surf%ns
         if (surf%frac(m, ind_veg_wall) + surf%frac(m, ind_pav_green) > 0.0_wp) then
            fac_wall(m) = surf%frac(m, ind_veg_wall)/ &
                          (surf%frac(m, ind_veg_wall) + surf%frac(m, ind_pav_green))
            fac_veg(m) = surf%frac(m, ind_pav_green)/ &
                         (surf%frac(m, ind_veg_wall) + surf%frac(m, ind_pav_green))
         end if
      end do
!
!-- If indoor model is used inner wall layer is calculated by using iwghf (indoor
!-- wall ground heat flux)
      if (.not. indoor_model) then
!$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
         do m = 1, surf%ns
            surf%iwghf_eb(m) = surf%lambda_h(nzt_wall, m)* &
                               (t_wall%val(nzt_wall + 1, m) - t_wall%val(nzt_wall, m))* &
                               surf%ddz_wall_center(nzt_wall, m)

            surf%iwghf_eb_window(m) = surf%lambda_h_window(nzt_wall, m)* &
                                      (t_window%val(nzt_wall + 1, m) - t_window%val(nzt_wall, m))* &
                                      surf%ddz_window_center(nzt_wall, m)
         end do
      end if
!
!-- Compute tendency terms for wall fraction.
!$OMP PARALLEL DO PRIVATE (kw, m, wtend) SCHEDULE (STATIC)
      do m = 1, surf%ns
         wtend(nzb_wall, m) = (1.0_wp/surf%rho_c_wall(nzb_wall, m))* &
                              (surf%lambda_h_layer(nzb_wall, m)* &
                               (t_wall%val(nzb_wall + 1, m) - t_wall%val(nzb_wall, m))* &
                               surf%ddz_wall_center(nzb_wall, m) &
                               + fac_wall(m)*surf%wghf_eb(m) &
                               - fac_veg(m)*(surf%lambda_h_green(nzt_wall, m)* &
                                             surf%ddz_green_center(nzt_wall, m) &
                                             + surf%lambda_h_layer(nzb_wall, m)* &
                                             surf%ddz_wall_center(nzb_wall, m) &
                                             )/ &
                               (surf%dz_green_center(nzt_wall, m) + &
                                surf%dz_wall_center(nzb_wall, m))*4.0_wp* &
                               (t_wall%val(nzb_wall, m) - t_green%val(nzt_wall, m)) &
                               )*surf%ddz_wall(nzb_wall, m)

         do kw = nzb_wall + 1, nzt_wall - 1
            wtend(kw, m) = (1.0_wp/surf%rho_c_wall(kw, m))* &
                           (surf%lambda_h_layer(kw, m)*(t_wall%val(kw + 1, m) - t_wall%val(kw, m))* &
                            surf%ddz_wall_center(kw, m) - &
                            surf%lambda_h_layer(kw - 1, m)*(t_wall%val(kw, m) - t_wall%val(kw - 1, m))* &
                            surf%ddz_wall_center(kw - 1, m) &
                            )*surf%ddz_wall(kw, m)
         end do
         wtend(nzt_wall, m) = (1.0_wp/surf%rho_c_wall(nzt_wall, m))* &
                              (-surf%lambda_h_layer(nzt_wall - 1, m)* &
                               (t_wall%val(nzt_wall, m) - t_wall%val(nzt_wall - 1, m))* &
                               surf%ddz_wall_center(nzt_wall - 1, m) + surf%iwghf_eb(m) &
                               )*surf%ddz_wall(nzt_wall, m)

         do kw = nzb_wall, nzt_wall
            t_wall_p%val(kw, m) = t_wall%val(kw, m) + dt_3d*(tsc(2)*wtend(kw, m) + &
                                                             tsc(3)*surf%tt_wall_m(kw, m))
         end do
      end do
!
!-- Compute tendency terms for all window fractions.
!-- Skip this during spinup. During the spinup, the tempeature inside window layers is not
!-- calculated to make larger timesteps possible.
      if (.not. spinup_phase) then
         wintend = 0.0_wp

!$OMP PARALLEL DO PRIVATE (kw, m, wintend, win_absorp) SCHEDULE (STATIC)
         do m = 1, surf%ns
!
!--       Reflectivity in glass windows is considered as equal on frontal and rear side of the
!--       glass, which together make total reflectivity (albedo for win fraction).
            win_nonrefl_1side = 1.0_wp - (surf%albedo(m, ind_wat_win) + surf%transmissivity(m) &
                                          + 1.0_wp &
                                          - sqrt((surf%albedo(m, ind_wat_win) &
                                                  + surf%transmissivity(m) + 1.0_wp)**2 &
                                                 - 4.0_wp*surf%albedo(m, ind_wat_win)))/2.0_wp
!
!--       Absorption coefficient is calculated using zw from internal tranmissivity, which only
!--       considers absorption without the effects of reflection.
            win_absorp = -log((surf%transmissivity(m) + surf%albedo(m, ind_wat_win) - 1.0_wp &
                               + win_nonrefl_1side)/win_nonrefl_1side) &
                         /surf%zw_window(nzt_wall, m)
!
!--       Prognostic equation for ground/roof window temperature t_window takes absorption
!--       of shortwave radiation into account
            wintend(nzb_wall, m) = (1.0_wp/surf%rho_c_window(nzb_wall, m)) &
                                   *(surf%lambda_h_window_layer(nzb_wall, m) &
                                     *(t_window%val(nzb_wall + 1, m) - t_window%val(nzb_wall, m)) &
                                     *surf%ddz_window_center(nzb_wall, m) &
                                     + surf%wghf_eb_window(m) &
                                     + surf%rad_sw_in(m)*win_nonrefl_1side &
                                     *(1.0_wp - exp(-win_absorp*surf%zw_window(nzb_wall, m))) &
                                     )*surf%ddz_window(nzb_wall, m)

            do kw = nzb_wall + 1, nzt_wall - 1
               wintend(kw, m) = (1.0_wp/surf%rho_c_window(kw, m)) &
                                *(surf%lambda_h_window_layer(kw, m) &
                                  *(t_window%val(kw + 1, m) - t_window%val(kw, m)) &
                                  *surf%ddz_window_center(kw, m) &
                                  - surf%lambda_h_window_layer(kw - 1, m) &
                                  *(t_window%val(kw, m) - t_window%val(kw - 1, m)) &
                                  *surf%ddz_window_center(kw - 1, m) &
                                  + surf%rad_sw_in(m)*win_nonrefl_1side &
                                  *(exp(-win_absorp*surf%zw_window(kw - 1, m)) &
                                    - exp(-win_absorp*surf%zw_window(kw, m))) &
                                  )*surf%ddz_window(kw, m)

            end do

            wintend(nzt_wall, m) = (1.0_wp/surf%rho_c_window(nzt_wall, m)) &
                                   *(-surf%lambda_h_window_layer(nzt_wall - 1, m) &
                                     *(t_window%val(nzt_wall, m) - t_window%val(nzt_wall - 1, m)) &
                                     *surf%ddz_window_center(nzt_wall - 1, m) &
                                     + surf%iwghf_eb_window(m) &
                                     + surf%rad_sw_in(m)*win_nonrefl_1side &
                                     *(exp(-win_absorp*surf%zw_window(nzt_wall - 1, m)) &
                                       - exp(-win_absorp*surf%zw_window(nzt_wall, m))) &
                                     )*surf%ddz_window(nzt_wall, m)

            do kw = nzb_wall, nzt_wall
               t_window_p%val(kw, m) = t_window%val(kw, m) + dt_3d*(tsc(2)*wintend(kw, m) + &
                                                                    tsc(3)*surf%tt_window_m(kw, m))
            end do
         end do
      end if
!
!-- Calculate weighted Runge-Kutta t_wall tendencies for the next substep.
      if (runge_l) then
         if (intermediate_timestep_count == 1) then
!$OMP PARALLEL DO PRIVATE (kw, m) SCHEDULE (STATIC)
            do m = 1, surf%ns
               do kw = nzb_wall, nzt_wall
                  surf%tt_wall_m(kw, m) = wtend(kw, m)
               end do
            end do
         elseif (intermediate_timestep_count < intermediate_timestep_count_max) then
!$OMP PARALLEL DO PRIVATE (kw, m) SCHEDULE (STATIC)
            do m = 1, surf%ns
               do kw = nzb_wall, nzt_wall
                  surf%tt_wall_m(kw, m) = -9.5625_wp*wtend(kw, m) + &
                                          5.3125_wp*surf%tt_wall_m(kw, m)
               end do
            end do
         end if
      end if
!
!-- Calculate weighted Runge-Kutta t_window tendencies for the next substep. Skip this during
!-- the spinup phase.
      if (.not. spinup_phase) then
         if (runge_l) then
            if (intermediate_timestep_count == 1) then
!$OMP PARALLEL DO PRIVATE (kw, m) SCHEDULE (STATIC)
               do m = 1, surf%ns
                  do kw = nzb_wall, nzt_wall
                     surf%tt_window_m(kw, m) = wintend(kw, m)
                  end do
               end do
            elseif (intermediate_timestep_count < intermediate_timestep_count_max) then
!$OMP PARALLEL DO PRIVATE (kw, m) SCHEDULE (STATIC)
               do m = 1, surf%ns
                  do kw = nzb_wall, nzt_wall
                     surf%tt_window_m(kw, m) = -9.5625_wp*wintend(kw, m) + &
                                               5.3125_wp*surf%tt_window_m(kw, m)
                  end do
               end do
            end if
         end if
      end if

      if (debug_output_timestep) then
         write (debug_string, *) 'usm_wall_heat_model: ', spinup_phase
         call debug_message(debug_string, 'end')
      end if

   end subroutine usm_wall_heat_model

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Green and substrate model as part of the urban surface model. The model predicts ground
!> temperatures.
!>
!> Important: green-heat model crashes due to unknown reason. Green fraction is thus set to zero
!> (in favor of wall fraction).
!> Note, usm_green_heat_model has not been vectorized yet.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_green_heat_model

      implicit none

      integer(iwp)  ::  i, j, k, kw, m      !< running indices

      logical  ::  conserve_water_content = .true.  !<

      real(wp)  ::  drho_l_lv               !< frequently used parameter
      real(wp)  ::  h_vg                    !< Van Genuchten coef. h
      real(wp)  ::  ke, lambda_h_green_sat  !< heat conductivity for saturated soil

      real(wp), dimension(nzb_wall:nzt_wall)   ::  gtend, tend         !< tendency
      real(wp), dimension(nzb_wall:nzt_wall)   ::  root_extr_green    !<

      real(wp), dimension(nzb_wall:nzt_wall + 1) ::  gamma_green_temp   !< temp. gamma
      real(wp), dimension(nzb_wall:nzt_wall + 1) ::  lambda_green_temp  !< temp. lambda

      type(surf_type), pointer                 ::  surf               !< surface-date type variable

      if (debug_output_timestep) then
         write (debug_string, *) 'usm_green_heat_model: '
         call debug_message(debug_string, 'start')
      end if

      drho_l_lv = 1.0_wp/(rho_l*l_v)
!
!-- Set pointer to urban-surface structure. Note, with all surfaces in one array this is
!-- actually not necessary any more. However, for future developments when land- and urban-
!-- surface model will be merged, this will become necessary again.
      surf => surf_usm

!-- Set tendency array for soil moisture to zero
      if (surf%ns > 0) then
         if (intermediate_timestep_count == 1) surf%tswc_m = 0.0_wp
      end if

!$OMP PARALLEL DO PRIVATE (m, i, j, k, kw, lambda_h_green_sat, ke, lambda_green_temp,      &
!$OMP&  gtend, tend, h_vg, gamma_green_temp, m_total, root_extr_green) SCHEDULE (STATIC)
      do m = 1, surf%ns
!
!--    For green fraction at upward-facing walls
         if (surf%frac(m, ind_pav_green) > 0.0_wp .and. surf%upward(m)) then
!
!--      Obtain indices
            i = surf%i(m)
            j = surf%j(m)
            k = surf%k(m)

            do kw = nzb_wall, nzt_wall
!
!--          Calculate volumetric heat capacity of the soil, taking into account water content
               surf%rho_c_total_green(kw, m) = surf%rho_c_green(kw, m)*(1.0_wp - swc_sat%val(kw, m)) &
                                               + rho_c_water*swc%val(kw, m)
!
!--          Calculate soil heat conductivity at the center of the soil layers
               lambda_h_green_sat = lambda_h_green_sm**(1.0_wp - swc_sat%val(kw, m)) &
                                    *lambda_h_water**swc%val(kw, m)

               ke = 1.0_wp + log10(max(0.1_wp, swc%val(kw, m)/swc_sat%val(kw, m)))

               lambda_green_temp(kw) = ke*(lambda_h_green_sat - lambda_h_green_dry) &
                                       + lambda_h_green_dry

            end do
            lambda_green_temp(nzt_wall + 1) = lambda_green_temp(nzt_wall)

!
!--       Calculate soil heat conductivity (lambda_h) at the _center level using weighting
            do kw = nzb_wall, nzt_wall - 1
               surf%lambda_h_green(kw, m) = (lambda_green_temp(kw)*surf%dz_green(kw, m) &
                                             + lambda_green_temp(kw + 1)*surf%dz_green(kw + 1, m) &
                                             )*0.5_wp*surf%ddz_green_center(kw, m)
            end do
            surf%lambda_h_green(nzt_wall, m) = lambda_green_temp(nzt_wall)

            t_green%val(nzt_wall + 1, m) = t_wall%val(nzb_wall, m)
!
!--       Prognostic equation for ground/roof temperature t_green
            gtend(:) = 0.0_wp
            gtend(nzb_wall) = (1.0_wp/surf%rho_c_total_green(nzb_wall, m)) &
                              *(surf%lambda_h_green(nzb_wall, m) &
                                *(t_green%val(nzb_wall + 1, m) - t_green%val(nzb_wall, m)) &
                                *surf%ddz_green_center(nzb_wall, m) &
                                + surf%wghf_eb_green(m) &
                                )*surf%ddz_green(nzb_wall, m)

            do kw = nzb_wall + 1, nzt_wall
               gtend(kw) = (1.0_wp/surf%rho_c_total_green(kw, m)) &
                           *(surf%lambda_h_green(kw, m) &
                             *(t_green%val(kw + 1, m) - t_green%val(kw, m)) &
                             *surf%ddz_green_center(kw, m) &
                             - surf%lambda_h_green(kw - 1, m) &
                             *(t_green%val(kw, m) - t_green%val(kw - 1, m)) &
                             *surf%ddz_green_center(kw - 1, m) &
                             )*surf%ddz_green(kw, m)
            end do

            t_green_p%val(nzb_wall:nzt_wall, m) = t_green%val(nzb_wall:nzt_wall, m) &
                                                  + dt_3d*(tsc(2)*gtend(nzb_wall:nzt_wall) + tsc(3) &
                                                           *surf%tt_green_m(nzb_wall:nzt_wall, m))

!
!--       Calculate t_green tendencies for the next Runge-Kutta step
            if (timestep_scheme(1:5) == 'runge') then
               if (intermediate_timestep_count == 1) then
                  do kw = nzb_wall, nzt_wall
                     surf%tt_green_m(kw, m) = gtend(kw)
                  end do
               elseif (intermediate_timestep_count < intermediate_timestep_count_max) then
                  do kw = nzb_wall, nzt_wall
                     surf%tt_green_m(kw, m) = -9.5625_wp*gtend(kw) + &
                                              5.3125_wp*surf%tt_green_m(kw, m)
                  end do
               end if
            end if

            do kw = nzb_wall, nzt_wall

!
!--          Calculate soil diffusivity at the center of the soil layers
               lambda_green_temp(kw) = (-b_ch*surf%gamma_w_green_sat(kw, m)*psi_sat &
                                        /swc_sat%val(kw, m)) &
                                       *(max(swc%val(kw, m), wilt%val(kw, m)) &
                                         /swc_sat%val(kw, m))**(b_ch + 2.0_wp)

!
!--          Parametrization of Van Genuchten
               if (soil_type /= 7) then
!
!--             Calculate the hydraulic conductivity after Van Genuchten (1980)
                  h_vg = (((swc_res%val(kw, m) - swc_sat%val(kw, m)) &
                           /(swc_res%val(kw, m) - &
                             max(swc%val(kw, m), wilt%val(kw, m))))** &
                          (surf%n_vg_green(m)/(surf%n_vg_green(m) - 1.0_wp)) &
                          - 1.0_wp &
                          )**(1.0_wp/surf%n_vg_green(m))/surf%alpha_vg_green(m)

                  gamma_green_temp(kw) = surf%gamma_w_green_sat(kw, m) &
                                         *(((1.0_wp + (surf%alpha_vg_green(m)*h_vg)** &
                                             surf%n_vg_green(m))** &
                                            (1.0_wp - 1.0_wp/surf%n_vg_green(m)) &
                                            - (surf%alpha_vg_green(m)*h_vg)** &
                                            (surf%n_vg_green(m) - 1.0_wp))**2 &
                                           )/((1.0_wp + (surf%alpha_vg_green(m)*h_vg)** &
                                               surf%n_vg_green(m))** &
                                              ((1.0_wp - 1.0_wp/surf%n_vg_green(m)) &
                                               *(surf%l_vg_green(m) + 2.0_wp)) &
                                              )

!
!--          Parametrization of Clapp & Hornberger
               else
                  gamma_green_temp(kw) = surf%gamma_w_green_sat(kw, m)*(swc%val(kw, m) &
                                                                        /swc_sat%val(kw, m))**(2.0_wp*b_ch + 3.0_wp)
               end if

            end do

!
!--       Prognostic equation for soil moisture content. Only performed, when humidity is enabled in
!--       the atmosphere
            if (humidity) then
!
!--          Calculate soil diffusivity (lambda_w) at the _center level using weighting
!--          To do: replace this with ECMWF-IFS Eq. 8.81
               do kw = nzb_wall, nzt_wall - 1

                  surf%lambda_w_green(kw, m) = (lambda_green_temp(kw)*surf%dz_green(kw, m) &
                                                + lambda_green_temp(kw + 1)*surf%dz_green(kw + 1, m) &
                                                )*0.5_wp*surf%ddz_green_center(kw, m)
                  surf%gamma_w_green(kw, m) = (gamma_green_temp(kw)*surf%dz_green(kw, m) &
                                               + gamma_green_temp(kw + 1)*surf%dz_green(kw + 1, m) &
                                               )*0.5_wp*surf%ddz_green_center(kw, m)

               end do

!
!--          In case of a closed bottom (= water content is conserved), set hydraulic conductivity
!--          to zero so that no water will be lost in the bottom layer.
               if (conserve_water_content) then
                  surf%gamma_w_green(kw, m) = 0.0_wp
               else
                  surf%gamma_w_green(kw, m) = gamma_green_temp(nzt_wall)
               end if

!--          The root extraction (= root_extr * qsws_veg / (rho_l * l_v)) ensures the mass
!--          conservation for water. The transpiration of plants equals the cumulative withdrawals
!--          by the roots in the soil. The scheme takes into account the availability of water in
!--          the soil layers as well as the root fraction in the respective layer. Layer with
!--          moisture below wilting point will not contribute, which reflects the preference of
!--          plants to take water from moister layers.

!
!--          Calculate the root extraction (ECMWF 7.69, the sum of root_extr = 1). The energy
!--          balance solver guarantees a positive transpiration, so that there is no need for an
!--          additional check.
               m_total = 0.0_wp
               do kw = nzb_wall, nzt_wall
                  if (swc%val(kw, m) > wilt%val(kw, m)) then
                     m_total = m_total + rootfr%val(kw, m)*swc%val(kw, m)
                  end if
               end do

               if (m_total > 0.0_wp) then
                  do kw = nzb_wall, nzt_wall
                     if (swc%val(kw, m) > wilt%val(kw, m)) then
                        root_extr_green(kw) = rootfr%val(kw, m)*swc%val(kw, m)/m_total
                     else
                        root_extr_green(kw) = 0.0_wp
                     end if
                  end do
               end if

!
!--          Prognostic equation for soil water content m_soil.
               tend(:) = 0.0_wp

               tend(nzb_wall) = (surf%lambda_w_green(nzb_wall, m) &
                                 *(swc%val(nzb_wall + 1, m) - swc%val(nzb_wall, m)) &
                                 *surf%ddz_green_center(nzb_wall, m) &
                                 - surf%gamma_w_green(nzb_wall, m) &
                                 - (root_extr_green(nzb_wall)*surf%qsws_veg(m) & !+ surf%qsws_soil_green(m)
                                    )*drho_l_lv) &
                                *surf%ddz_green(nzb_wall, m)

               do kw = nzb_wall + 1, nzt_wall - 1
                  tend(kw) = (surf%lambda_w_green(kw, m) &
                              *(swc%val(kw + 1, m) - swc%val(kw, m)) &
                              *surf%ddz_green_center(kw, m) &
                              - surf%gamma_w_green(kw, m) &
                              - surf%lambda_w_green(kw - 1, m) &
                              *(swc%val(kw, m) - swc%val(kw - 1, m)) &
                              *surf%ddz_green_center(kw - 1, m) &
                              + surf%gamma_w_green(kw - 1, m) &
                              - (root_extr_green(kw) &
                                 *surf%qsws_veg(m) &
                                 *drho_l_lv) &
                              )*surf%ddz_green(kw, m)

               end do
               tend(nzt_wall) = (-surf%gamma_w_green(nzt_wall, m) &
                                 - surf%lambda_w_green(nzt_wall - 1, m) &
                                 *(swc%val(nzt_wall, m) - swc%val(nzt_wall - 1, m)) &
                                 *surf%ddz_green_center(nzt_wall - 1, m) &
                                 + surf%gamma_w_green(nzt_wall - 1, m) &
                                 - (root_extr_green(nzt_wall) &
                                    *surf%qsws_veg(m) &
                                    *drho_l_lv) &
                                 )*surf%ddz_green(nzt_wall, m)

               swc_p%val(nzb_wall:nzt_wall, m) = swc%val(nzb_wall:nzt_wall, m) + dt_3d &
                                                 *(tsc(2)*tend(:) + tsc(3) &
                                                   *surf%tswc_m(:, m) &
                                                   )

!
!--          Account for dry soils (find a better solution here!)
               do kw = nzb_wall, nzt_wall
                  if (swc_p%val(kw, m) < 0.0_wp) swc_p%val(kw, m) = 0.0_wp
               end do

!
!--          Calculate m_soil tendencies for the next Runge-Kutta step
               if (timestep_scheme(1:5) == 'runge') then
                  if (intermediate_timestep_count == 1) then
                     do kw = nzb_wall, nzt_wall
                        surf%tswc_m(kw, m) = tend(kw)
                     end do
                  elseif (intermediate_timestep_count < intermediate_timestep_count_max) then
                     do kw = nzb_wall, nzt_wall
                        surf%tswc_m(kw, m) = -9.5625_wp*tend(kw) + 5.3125_wp*surf%tswc_m(kw, m)
                     end do
                  end if
               end if
            end if
!
!--    Vertical walls
         elseif (surf%frac(m, ind_pav_green) > 0.0_wp) then
!
!--            No substrate layer for green walls / only groundbase green walls (ivy i.e.) -> Green layers get
!--            same temperature as first wall layer, therefore no temperature calculations for vertical green
!--            substrate layers now

!
! !
! !--          Obtain indices
!              i = surf%i(m)
!              j = surf%j(m)
!              k = surf%k(m)
!
!              t_green%val(nzt_wall+1,m) = t_wall%val(nzb_wall,m)
! !
! !--          Prognostic equation for green temperature t_green_v
!              gtend(:) = 0.0_wp
!              gtend(nzb_wall) = (1.0_wp / surf%rho_c_green(nzb_wall,m)) *                        &
!                                      ( surf%lambda_h_green(nzb_wall,m) *                        &
!                                        ( t_green%val(nzb_wall+1,m)                              &
!                                        - t_green%val(nzb_wall,m) ) *                            &
!                                        surf%ddz_green(nzb_wall+1,m)                             &
!                                      + surf%wghf_eb(m) ) *                                      &
!                                        surf%ddz_green_stag(nzb_wall,m)
!
!              DO  kw = nzb_wall+1, nzt_wall
!                 gtend(kw) = (1.0_wp / surf%rho_c_green(kw,m))                                   &
!                           * (   surf%lambda_h_green(kw,m)                                       &
!                             * ( t_green%val(kw+1,m) - t_green%val(kw,m) )                       &
!                             * surf%ddz_green(kw+1,m)                                            &
!                           - surf%lambda_h(kw-1,m)                                               &
!                             * ( t_green%val(kw,m) - t_green%val(kw-1,m) )                       &
!                             * surf%ddz_green(kw,m) )                                            &
!                           * surf%ddz_green_stag(kw,m)
!              ENDDO
!
!              t_green_v_p(l)%val(nzb_wall:nzt_wall,m) =                                          &
!                                   t_green%val(nzb_wall:nzt_wall,m)                              &
!                                 + dt_3d * ( tsc(2)                                              &
!                                 * gtend(nzb_wall:nzt_wall) + tsc(3)                             &
!                                 * surf%tt_green_m(nzb_wall:nzt_wall,m) )
!
! !
! !--          Calculate t_green tendencies for the next Runge-Kutta step
!              IF ( timestep_scheme(1:5) == 'runge' )  THEN
!                  IF ( intermediate_timestep_count == 1 )  THEN
!                     DO  kw = nzb_wall, nzt_wall
!                        surf%tt_green_m(kw,m) = gtend(kw)
!                     ENDDO
!                  ELSEIF ( intermediate_timestep_count <                                         &
!                           intermediate_timestep_count_max )  THEN
!                      DO  kw = nzb_wall, nzt_wall
!                         surf%tt_green_m(kw,m) =                                                 &
!                                     - 9.5625_wp * gtend(kw) +                                   &
!                                       5.3125_wp * surf%tt_green_m(kw,m)
!                      ENDDO
!                  ENDIF
!              ENDIF

!
!--       Workaround, set green surface temperature to wall temperature.
            do kw = nzb_wall, nzt_wall + 1
               t_green%val(kw, m) = t_wall%val(nzb_wall, m)
            end do
!
!--       Workaround, rho_c_total_green is used for calculation of the max_dt_green_column
            do kw = nzb_wall, nzt_wall
               surf%rho_c_total_green(kw, m) = surf%rho_c_green(kw, m)
            end do
         end if
      end do

      if (debug_output_timestep) then
         write (debug_string, *) 'usm_green_heat_model: '
         call debug_message(debug_string, 'end')
      end if

   end subroutine usm_green_heat_model

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &usm_par for urban surface model
!--------------------------------------------------------------------------------------------------!
   subroutine usm_parin

      implicit none

      character(LEN=100)  ::  line  !< string containing current line of file PARIN

      integer(iwp) ::  io_status    !< status after reading the namelist file

      logical ::  switch_off_module = .false.  !< local namelist parameter to switch off the module
      !< although the respective module namelist appears in
      !< the namelist file

!> @todo remove this variable with next release (24.10)
      logical ::  usm_wall_mod = .false.       !< urban-surface-parameter, which has no effect any more

      namelist /urban_surface_parameters/ &
         building_type, &
         roof_category, &
         roof_inner_temperature, &
         roughness_concrete, &
         switch_off_module, &
         usm_wall_mod, &
         wall_category, &
         wall_inner_temperature, &
         window_inner_temperature

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
      rewind (11)
      read (11, urban_surface_parameters, IOSTAT=io_status)

!
!-- Action depending on the READ status
      if (io_status == 0) then
!
!--    urban_surface_parameters namelist was found and read correctly. Set flag that indicates that
!--    the urban surface model is switched on.
         if (.not. switch_off_module) urban_surface = .true.

      elseif (io_status > 0) then
!
!--    urban_surface_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
         backspace (11)
         read (11, '(A)') line
         call parin_fail_message('urban_surface_parameters', line)

      end if

!
!-- Issue a warning, if usm_wall_mod (which has no effect any more) has been set.
      if (usm_wall_mod) then
         message_string = 'parameter "usm_wall_mod" has no effect any more, &'// &
                          'please remove it from the urban_surface_parameters namelist, &'// &
                          'an adjustment of conductivity to prevent instability&is not required '// &
                          'any more'
         call message('usm_parin', 'PAC0363', 0, 1, 0, 6, 0)
      end if

   end subroutine usm_parin

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!> Soubroutine reads t_surf and t_wall.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_rrd_local_ftn(k, nxlf, nxlc, nxl_on_file, nxrf, nxr_on_file, nynf, nyn_on_file, &
                                nysf, nysc, nys_on_file, found)

      use control_parameters, &
         only: length, &
               restart_string

      implicit none

      integer(iwp) ::  k                 !< running index over previous input files covering current local domain
      integer(iwp) ::  nxlc              !< index of left boundary on current subdomain
      integer(iwp) ::  nxlf              !< index of left boundary on former subdomain
      integer(iwp) ::  nxl_on_file       !< index of left boundary on former local domain
      integer(iwp) ::  nxrf              !< index of right boundary on former subdomain
      integer(iwp) ::  nxr_on_file       !< index of right boundary on former local domain
      integer(iwp) ::  nynf              !< index of north boundary on former subdomain
      integer(iwp) ::  nyn_on_file       !< index of north boundary on former local domain
      integer(iwp) ::  nysc              !< index of south boundary on current subdomain
      integer(iwp) ::  nysf              !< index of south boundary on former subdomain
      integer(iwp) ::  nys_on_file       !< index of south boundary on former local domain
      integer(iwp) ::  ns_on_file_usm    !< number of surface elements (urban type) on file
!
!-- Note, the save attribute in the following array declaration is necessary, in order to keep the
!-- number of urban surface elements on file during rrd_local calls.
      integer(iwp), dimension(:, :), allocatable, save ::  end_index_on_file    !<
      integer(iwp), dimension(:, :), allocatable, save ::  start_index_on_file  !<

      logical, intent(OUT)  ::  found  !<

      type(surf_type_1d_usm), save ::  tmp_surf   !< temporary variable to read surface data
      type(surf_type_2d_usm), save ::  tmp_wall   !< temporary variable to read wall data

      found = .true.

      select case (restart_string(1:length))

      case ('ns_on_file_usm')
         if (k == 1) then
            read (13) ns_on_file_usm
!
!--          In case of changing mpi topology, this routine could be called more than once.
!--          Hence, arrays need to be deallocated before allocated again.
            if (allocated(tmp_surf%val)) deallocate (tmp_surf%val)
            if (allocated(tmp_wall%val)) deallocate (tmp_wall%val)

!
!--          Allocate temporary arrays for reading data on file. Note, the size of allocated surface
!--          elements do not necessarily need to match the size of present surface elements on
!--          current processor, as the number of processors between restarts can change.
            allocate (tmp_surf%val(1:ns_on_file_usm))
            allocate (tmp_wall%val(nzb_wall:nzt_wall + 1, 1:ns_on_file_usm))
         end if

      case ('usm_start_index')
         if (k == 1) then

            if (allocated(start_index_on_file)) deallocate (start_index_on_file)

            allocate (start_index_on_file(nys_on_file:nyn_on_file, nxl_on_file:nxr_on_file))

            read (13) start_index_on_file

         end if

      case ('usm_end_index')
         if (k == 1) then

            if (allocated(end_index_on_file)) deallocate (end_index_on_file)

            allocate (end_index_on_file(nys_on_file:nyn_on_file, nxl_on_file:nxr_on_file))

            read (13) end_index_on_file

         end if

      case ('t_surf_wall')
         if (k == 1) then
            if (.not. allocated(t_surf_wall%val)) allocate (t_surf_wall%val(1:surf_usm%ns))
            read (13) tmp_surf%val
         end if
         call surface_restore_elements(t_surf_wall%val, tmp_surf%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('t_surf_window')
         if (k == 1) then
            if (.not. allocated(t_surf_window%val)) &
               allocate (t_surf_window%val(1:surf_usm%ns))
            read (13) tmp_surf%val
         end if
         call surface_restore_elements(t_surf_window%val, tmp_surf%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('t_surf_green')
         if (k == 1) then
            if (.not. allocated(t_surf_green%val)) &
               allocate (t_surf_green%val(1:surf_usm%ns))
            read (13) tmp_surf%val
         end if
         call surface_restore_elements(t_surf_green%val, tmp_surf%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('m_liq_usm')
         if (k == 1) then
            if (.not. allocated(m_liq_usm%val)) allocate (m_liq_usm%val(1:surf_usm%ns))
            read (13) tmp_surf%val
         end if
         call surface_restore_elements(m_liq_usm%val, tmp_surf%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('swc')
         if (k == 1) then
            if (.not. allocated(swc%val)) &
               allocate (swc%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            read (13) tmp_wall%val
         end if
         call surface_restore_elements(swc%val, tmp_wall%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('t_wall')
         if (k == 1) then
            if (.not. allocated(t_wall%val)) &
               allocate (t_wall%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            read (13) tmp_wall%val
         end if
         call surface_restore_elements(t_wall%val, tmp_wall%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('t_window')
         if (k == 1) then
            if (.not. allocated(t_window%val)) &
               allocate (t_window%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            read (13) tmp_wall%val
         end if
         call surface_restore_elements(t_window%val, tmp_wall%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case ('t_green')
         if (k == 1) then
            if (.not. allocated(t_green%val)) &
               allocate (t_green%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            read (13) tmp_wall%val
         end if
         call surface_restore_elements(t_green%val, tmp_wall%val, &
                                       surf_usm%start_index, start_index_on_file, &
                                       end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf, &
                                       nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file)

      case DEFAULT

         found = .false.

      end select

   end subroutine usm_rrd_local_ftn

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!> Soubroutine reads t_surf and t_wall.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_rrd_local_mpi

      integer(idp), dimension(nys:nyn, nxl:nxr) ::  global_end_index
      integer(idp), dimension(nys:nyn, nxl:nxr) ::  global_start_index

      logical ::  array_found  !<
      logical ::  data_to_read !< dummy variable

!-- At the moment reading of surface data in combination with cyclic fill is not realized,
!-- so that this is skipped for the moment.
      if (cyclic_fill_initialization) return

      call rd_mpi_io_check_array('usm_global_start', found=array_found)
      if (array_found) call rrd_mpi_io('usm_global_start', global_start_index)

      call rd_mpi_io_check_array('usm_global_end', found=array_found)
      if (array_found) call rrd_mpi_io('usm_global_end', global_end_index)
!
!-- Check if data input for surface-type variables is required. Note, only invoke routine if USM
!-- surface restart data is on file. In case of cyclic fill initialization this is not necessarily
!-- guaranteed. To check this use the array_found control flag.
      if (array_found) then
         call rd_mpi_io_surface_filetypes(surf_usm%start_index, surf_usm%end_index, &
                                          data_to_read, global_start_index, global_end_index)
      else
         data_to_read = .false.
      end if

      if (data_to_read) then
         call rd_mpi_io_check_array('t_surf_wall', found=array_found)
         if (array_found) then
            if (.not. allocated(t_surf_wall%val)) allocate (t_surf_wall%val(1:surf_usm%ns))
            call rrd_mpi_io_surface('t_surf_wall', t_surf_wall%val)
         end if

         call rd_mpi_io_check_array('t_surf_window', found=array_found)
         if (array_found) then
            if (.not. allocated(t_surf_window%val)) allocate (t_surf_window%val(1:surf_usm%ns))
            call rrd_mpi_io_surface('t_surf_window', t_surf_window%val)
         end if

         call rd_mpi_io_check_array('t_surf_green', found=array_found)
         if (array_found) then
            if (.not. allocated(t_surf_green%val)) allocate (t_surf_green%val(1:surf_usm%ns))
            call rrd_mpi_io_surface('t_surf_green', t_surf_green%val)
         end if

         call rd_mpi_io_check_array('m_liq_usm', found=array_found)
         if (array_found) then
            if (.not. allocated(m_liq_usm%val)) allocate (m_liq_usm%val(1:surf_usm%ns))
            call rrd_mpi_io_surface('m_liq_usm', m_liq_usm%val)
         end if

      end if

      call rd_mpi_io_check_array('usm_global_start_2', found=array_found)
      if (array_found) call rrd_mpi_io('usm_global_start_2', global_start_index)

      call rd_mpi_io_check_array('usm_global_end_2', found=array_found)
      if (array_found) call rrd_mpi_io('usm_global_end_2', global_end_index)

      if (array_found) then
         call rd_mpi_io_surface_filetypes(surf_usm%start_index, surf_usm%end_index, &
                                          data_to_read, global_start_index, global_end_index)
      else
         data_to_read = .false.
      end if

      if (data_to_read) then
         call rd_mpi_io_check_array('swc', found=array_found)
         if (array_found) then
            if (.not. allocated(swc%val)) allocate (swc%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            call rrd_mpi_io_surface('swc', swc%val)
         end if

         call rd_mpi_io_check_array('t_wall', found=array_found)
         if (array_found) then
            if (.not. allocated(t_wall%val)) &
               allocate (t_wall%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            call rrd_mpi_io_surface('t_wall', t_wall%val)
         end if

         call rd_mpi_io_check_array('t_window', found=array_found)
         if (array_found) then
            if (.not. allocated(t_window%val)) &
               allocate (t_window%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            call rrd_mpi_io_surface('t_window', t_window%val)
         end if

         call rd_mpi_io_check_array('t_green', found=array_found)
         if (array_found) then
            if (.not. allocated(t_green%val)) &
               allocate (t_green%val(nzb_wall:nzt_wall + 1, 1:surf_usm%ns))
            call rrd_mpi_io_surface('t_green', t_green%val)
         end if
      end if

   end subroutine usm_rrd_local_mpi

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface. It follows the basic ideas and
!> structure of lsm_energy_balance with many simplifications and adjustments.
!> TODO better description
!> No calculation of window surface temperatures during spinup to increase maximum possible timstep
!--------------------------------------------------------------------------------------------------!
   subroutine usm_energy_balance

      call usm_surface_energy_balance

      call usm_green_heat_model                      !< usm_green_heat_model has not been vectorized yet

      call usm_wall_heat_model

   end subroutine usm_energy_balance

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface. It follows the basic ideas and
!> structure of lsm_energy_balance with many simplifications and adjustments.
!> TODO better description
!> No calculation of window surface temperatures during spinup to increase maximum possible timstep
!--------------------------------------------------------------------------------------------------!
   subroutine usm_surface_energy_balance

      implicit none

      logical                           ::  runge_l        !< dummy flag to indicate RK timestepping scheme

      logical, dimension(1:surf_usm%ns) ::  horizontal                !< flag indicating horizontal surfaces
      logical, dimension(1:surf_usm%ns) ::  force_radiation_call_l_v  !< flag to gather information if radiation need to be called

      integer(iwp) ::  i              !< grid index for reference atmosphere cell, x-direction
      integer(iwp) ::  j              !< grid index for reference atmosphere cell, y-direction
      integer(iwp) ::  k              !< grid index for reference atmosphere cell, z-direction
      integer(iwp) ::  kk             !< loop index for depth of vegetation layer
      integer(iwp) ::  m              !< running index for surface elements
      integer(iwp) ::  i_off          !< offset to determine index of surface element, seen from atmospheric grid point, for x
      integer(iwp) ::  j_off          !< offset to determine index of surface element, seen from atmospheric grid point, for y
      integer(iwp) ::  k_off          !< offset to determine index of surface element, seen from atmospheric grid point, for z

      real(wp) ::  coef_1                  !< first coeficient for prognostic equation
      real(wp) ::  coef_2                  !< second  coeficient for prognostic equation
      real(wp) ::  e                       !< water vapour pressure
      real(wp) ::  e_s                     !< water vapour saturation pressure
      real(wp) ::  f1                      !< resistance correction term 1
      real(wp) ::  f3                      !< resistance correction term 3
      real(wp) ::  m_max_depth = 0.0002_wp !< Maximum capacity of the water reservoir (m)
      real(wp) ::  stend_wall              !< tendency for wall surfaces
      real(wp) ::  stend_window            !< tendency for window surfaces
      real(wp) ::  stend_green             !< tendency for green surfaces
      real(wp) ::  tend                    !< tendency
      real(wp) ::  ueff                    !< limited near-surface wind speed - used for calculation of resistance

      real(wp), dimension(1:surf_usm%ns) ::  coef_green_1  !< first coeficient for prognostic green wall equation
      real(wp), dimension(1:surf_usm%ns) ::  coef_green_2  !< second  coeficient for prognostic green wall equation
      real(wp), dimension(1:surf_usm%ns) ::  coef_window_1 !< first coeficient for prognostic window equation
      real(wp), dimension(1:surf_usm%ns) ::  coef_window_2 !< second  coeficient for prognostic window equation
      real(wp), dimension(1:surf_usm%ns) ::  dq_s_dt       !< derivate of q_s with respect to T
      real(wp), dimension(1:surf_usm%ns) ::  drho_l_lv     !< frequently used parameter for green layers
      real(wp), dimension(1:surf_usm%ns) ::  e_s_dt        !< derivate of e_s with respect to T
      real(wp), dimension(1:surf_usm%ns) ::  f_qsws        !< factor for qsws
      real(wp), dimension(1:surf_usm%ns) ::  f_qsws_veg    !< factor for qsws_veg
      real(wp), dimension(1:surf_usm%ns) ::  f_qsws_liq    !< factor for qsws_liq
      real(wp), dimension(1:surf_usm%ns) ::  f2            !< resistance correction term 2
      real(wp), dimension(1:surf_usm%ns) ::  frac_green    !< green fraction, used to restore original values during spinup
      real(wp), dimension(1:surf_usm%ns) ::  frac_win      !< window fraction, used to restore original values during spinup
      real(wp), dimension(1:surf_usm%ns) ::  frac_wall     !< wall fraction, used to restore original values during spinup
      real(wp), dimension(1:surf_usm%ns) ::  f_shf         !< factor for shf_eb
      real(wp), dimension(1:surf_usm%ns) ::  f_shf_green   !< factor for shf_eb green wall
      real(wp), dimension(1:surf_usm%ns) ::  f_shf_window  !< factor for shf_eb window
      real(wp), dimension(1:surf_usm%ns) ::  m_liq_max     !< maxmimum value of the liq. water reservoir
      real(wp), dimension(1:surf_usm%ns) ::  m_total       !< total soil moisture content
      real(wp), dimension(1:surf_usm%ns) ::  qv1           !< specific humidity at first grid level
      real(wp), dimension(1:surf_usm%ns) ::  q_s           !< saturation specific humidity
      real(wp), dimension(1:surf_usm%ns) ::  rho_cp        !< rho_wall_surface * c_p
      real(wp), dimension(1:surf_usm%ns) ::  rho_lv        !< rho_wall_surface * l_v

      type(surf_type), pointer ::  surf              !< surface-date type variable

      if (debug_output_timestep) then
         write (debug_string, *) 'usm_surface_energy_balance_vector:', spinup_phase
         call debug_message(debug_string, 'start')
      end if

      surf => surf_usm

      runge_l = (timestep_scheme(1:5) == 'runge')

      force_radiation_call_l_v = .false.
!
!-- Set control flags
      if (surf%ns > 0) horizontal = (surf%upward(1:surf%ns) .or. surf%downward(1:surf%ns))
!
!-- During spinup set green and window fraction to zero and restore at the end of the loop.
      if (spinup_phase) then
         frac_win = 0.0_wp
         frac_wall = 1.0_wp
         frac_green = 0.0_wp
      else
!$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
         do m = 1, surf%ns
            frac_win(m) = surf%frac(m, ind_wat_win)
            frac_wall(m) = surf%frac(m, ind_veg_wall)
            frac_green(m) = surf%frac(m, ind_pav_green)
         end do
      end if
!
!-- Precalculate frequently used parameters such as rho_cp and qv1.
!$OMP PARALLEL DO PRIVATE (m, i, j, k) SCHEDULE (STATIC)
      do m = 1, surf%ns
!
!--    Get indices of respective grid point.
         i = surf%i(m)
         j = surf%j(m)
         k = surf%k(m)

         rho_cp(m) = c_p*hyp(k)/(r_d*surf%pt1(m)*exner(k))

         if (frac_green(m) > 0.0_wp) then
            rho_lv(m) = rho_cp(m)/c_p*l_v
            drho_l_lv(m) = 1.0_wp/(rho_l*l_v)
         end if

         if (humidity) then
            qv1(m) = q(k, j, i)
         else
            qv1(m) = 0.0_wp
         end if
      end do

!
!-- Calculate aerodyamic resistance.
!$OMP PARALLEL DO PRIVATE (m, i, j, k, ueff) SCHEDULE (STATIC)
      do m = 1, surf%ns
         if (surf%upward(m)) then
!
!--       Calculation for horizontally upward facing surfaces follows LSM formulation.
!--       pt, us, ts are not available for the prognostic time step, data from the
!--       last time step is used here.
            surf%r_a(m) = (surf%pt1(m) - surf%pt_surface(m))/ &
                          (surf%ts(m)*surf%us(m) + 1.0e-20_wp)
         else
!
!--       Get indices of respective grid point.
            i = surf%i(m)
            j = surf%j(m)
            k = surf%k(m)
!
!--       Calculation of r_a for vertical and downward facing horizontal surfaces
!--
!--       Heat transfer coefficient for forced convection along vertical walls follows formulation
!--       in TUF3d model (Krayenhoff & Voogt, 2006)
!--
!--       H = httc (Tsfc - Tair)
!--       httc = rw * (11.8 + 4.2 * Ueff) - 4.0
!--
!--            rw: Wall patch roughness relative to 1.0 for concrete
!--            Ueff: Effective wind speed
!--            - 4.0 is a reduction of Rowley et al (1930) formulation based on
!--            Cole and Sturrock (1977)
!--
!--            Ucan: Canyon wind speed
!--            wstar: Convective velocity
!--            Qs: Surface heat flux
!--            zH: Height of the convective layer
!--            wstar = (g/Tcan*Qs*zH)**(1./3.)
!--       Effective velocity components must always be defined at scalar grid point. The wall
!--       normal component is obtained by simple linear interpolation. (An alternative would be an
!--       logarithmic interpolation.) Parameter roughness_concrete (default value = 0.001) is used
!--       to calculation of roughness relative to concrete. Note, wind velocity is limited
!--       to avoid division by zero. The nominator can become <= 0.0 for values z0 < 3*10E-4.
            ueff = max(sqrt(((u(k, j, i) + u(k, j, i + 1))*0.5_wp)**2 + &
                            ((v(k, j, i) + v(k, j + 1, i))*0.5_wp)**2 + &
                            ((w(k, j, i) + w(k - 1, j, i))*0.5_wp)**2), &
                       ((4.0_wp + 0.1_wp)/(surf%z0(m)*d_roughness_concrete) - 11.8_wp) &
                       /4.2_wp &
                       )

            surf%r_a(m) = rho_cp(m)/(surf%z0(m)*d_roughness_concrete* &
                                     (11.8_wp + 4.2_wp*ueff) - 4.0_wp)
         end if
      end do

      if (surf%ns > 0) then
!
!--    Make sure that the resistance does not drop to zero and does not exceed a maxmium value in
!--    case of zero velocities.
         where (surf%r_a(1:surf%ns) < 1.0_wp) surf%r_a(1:surf%ns) = 1.0_wp
         where (surf%r_a(1:surf%ns) > 300.0_wp) surf%r_a(1:surf%ns) = 300.0_wp
!
!--    Aeorodynamical resistance for the window and green fractions are set to the same value.
         surf%r_a_window(1:surf%ns) = surf%r_a(1:surf%ns)
         surf%r_a_green(1:surf%ns) = surf%r_a(1:surf%ns)
!
!--    Factor for shf_eb.
         f_shf(1:surf%ns) = rho_cp(1:surf%ns)/surf%r_a(1:surf%ns)
         f_shf_window(1:surf%ns) = rho_cp(1:surf%ns)/surf%r_a_window(1:surf%ns)
         f_shf_green(1:surf%ns) = rho_cp(1:surf%ns)/surf%r_a_green(1:surf%ns)

      end if

!$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
      do m = 1, surf%ns
         if (frac_green(m) > 0.0_wp) then
            if (surf%upward(m)) then
!
!--          For upward-facing surfaces, compute soil moisture content. This is required for
!--          correction factor f2.
               m_total(m) = 0.0_wp
               do kk = nzb_wall, nzt_wall + 1
                  m_total(m) = m_total(m) + rootfr%val(nzb_wall, m)* &
                               max(swc%val(nzb_wall, m), wilt%val(nzb_wall, m))
               end do
!
!--          f2: Correction for soil moisture availability to plants (the integrated soil moisture
!--          must thus be considered here). f2 = 0 for very dry soils.
               if (m_total(m) > wilt%val(nzb_wall, m) .and. m_total(m) < fc%val(nzb_wall, m)) then
                  f2(m) = (m_total(m) - wilt%val(nzb_wall, m))/ &
                          (fc%val(nzb_wall, m) - wilt%val(nzb_wall, m))
               elseif (m_total(m) >= fc%val(nzb_wall, m)) then
                  f2(m) = 1.0_wp
               else
                  f2(m) = 1.0e-20_wp
               end if
            else
!
!--          f2 = 1 for vertical surfaces.
               f2(m) = 1.0_wp
            end if
         end if
      end do

!$OMP PARALLEL DO PRIVATE (m, f1, f3, e_s, e, coef_1, coef_2) SCHEDULE (STATIC)
      do m = 1, surf%ns

         if (frac_green(m) > 0.0_wp) then
!
!--       Adapted from LSM:
!--       Second step: calculate canopy resistance r_canopy. f1-f3 here are defined as 1/f1-f3
!--       as in ECMWF documentation.
!--       f1: Correction for incoming shortwave radiation (stomata close at night).
            f1 = min(1.0_wp, (0.004_wp*surf%rad_sw_in(m) + 0.05_wp)/ &
                     (0.81_wp*(0.004_wp*surf%rad_sw_in(m) + 1.0_wp)))
!
!--       Calculate water vapour pressure at saturation.
            e_s = 0.01_wp*magnus_tl(t_surf_green%val(m))
!
!--       f3: Correction for vapour pressure deficit.
            if (surf%g_d(m) /= 0.0_wp) then
!
!--          Calculate vapour pressure.
               e = qv1(m)*surface_pressure/(qv1(m) + 0.622_wp)
               f3 = exp(-surf%g_d(m)*(e_s - e))
            else
               f3 = 1.0_wp
            end if
!
!--       Calculate canopy resistance. In case that c_veg is 0 (bare soils), this calculation is
!--       obsolete, as r_canopy is not used below.
!--       To do: check for very dry soil -> r_canopy goes to infinity.
            surf%r_canopy(m) = surf%r_canopy_min(m)/(surf%lai(m)*f1*f2(m)*f3 + 1.0e-20_wp)
!
!--       Calculate saturation specific humidity.
            q_s(m) = 0.622_wp*e_s/(surface_pressure - e_s)
!
!--       In case of dewfall, set evapotranspiration to zero.
!--       All super-saturated water is then removed from the air.
            if (humidity .and. q_s(m) <= qv1(m)) then
               surf%r_canopy(m) = 0.0_wp
            end if

            if (surf%upward(m)) then
!
!--          Calculate the maximum possible liquid water amount on plants and bare surface. For
!--          vegetated surfaces, a maximum depth of 0.2 mm is assumed, while paved surfaces might
!--          hold up 1 mm of water. The liquid water fraction for paved surfaces is calculated after
!--          Noilhan & Planton (1989), while the ECMWF formulation is used for vegetated surfaces
!--          and bare soils.
               m_liq_max(m) = m_max_depth*(surf%lai(m))
               surf%c_liq(m) = min(1.0_wp, (m_liq_usm%val(m)/m_liq_max(m))**0.67)

!
!--          Calculate coefficients for the total evapotranspiration.
!--          In case of water surface, set vegetation and soil fluxes to zero.
!--          For pavements, only evaporation of liquid water is possible.
               f_qsws_veg(m) = rho_lv(m)*(1.0_wp - surf%c_liq(m))/ &
                               (surf%r_a_green(m) + surf%r_canopy(m))
               f_qsws_liq(m) = rho_lv(m)*surf%c_liq(m)/surf%r_a_green(m)
               f_qsws(m) = f_qsws_veg(m) + f_qsws_liq(m)
            else
               f_qsws_veg(m) = rho_lv(m)*(1.0_wp - 0.0_wp)/ &
                               (surf%r_a_green(m) + surf%r_canopy(m))
               f_qsws_liq(m) = 0.0_wp ! rho_lv(m) * surf%c_liq(m) / surf%r_a_green(m)
               f_qsws(m) = f_qsws_veg(m) + f_qsws_liq(m)
            end if
!
!--       Calculate derivative of q_s for Taylor series expansion.
            e_s_dt(m) = e_s*(17.269_wp/(t_surf_green%val(m) - 35.86_wp) &
                             - 17.269_wp*(t_surf_green%val(m) - degc_to_k)/ &
                             (t_surf_green%val(m) - 35.86_wp)**2 &
                             )
            dq_s_dt(m) = 0.622_wp*e_s_dt(m)/(surface_pressure - e_s_dt(m))
         end if
      end do
!
!-- Add LW up so that it can be removed in prognostic equation.
      if (surf%ns > 0) then
         surf%rad_net_l(1:surf%ns) = surf%rad_sw_in(1:surf%ns) - surf%rad_sw_out(1:surf%ns) + &
                                     surf%rad_lw_in(1:surf%ns) - surf%rad_lw_out(1:surf%ns)
      end if

!
!-- Compute coef_window_1 and coef_window_2.
      coef_window_1 = 0.0_wp
      coef_window_2 = 0.0_wp
      if (.not. spinup_phase) then
!$OMP PARALLEL DO PRIVATE (k, m) SCHEDULE (STATIC)
         do m = 1, surf%ns
            if (frac_win(m) > 0.0_wp) then
!
!--          Get k index of respective grid point.
               k = surf%k(m)
               coef_window_1(m) = surf%rad_net_l(m) + (3.0_wp + 1.0_wp)* &
                                  surf%emissivity(m, ind_wat_win)*sigma_sb* &
                                  t_surf_window%val(m)**4 + f_shf_window(m)*surf%pt1(m) + &
                                  surf%lambda_surf_window(m)*t_window%val(nzb_wall, m)

               coef_window_2(m) = 4.0_wp*surf%emissivity(m, ind_wat_win)*sigma_sb* &
                                  t_surf_window%val(m)**3 + surf%lambda_surf_window(m) + &
                                  f_shf_window(m)/exner(k)
            end if
         end do
      end if
!
!-- Compute coef_green_1 and coef_green_2.
!$OMP PARALLEL DO PRIVATE (k, m) SCHEDULE (STATIC)
      do m = 1, surf%ns
!
!--    Get k index of respective grid point.
         k = surf%k(m)
         if (humidity .and. frac_green(m) > 0.0_wp) then
            coef_green_1(m) = surf%rad_net_l(m) + (3.0_wp + 1.0_wp)* &
                              surf%emissivity(m, ind_pav_green)*sigma_sb* &
                              t_surf_green%val(m)**4 + f_shf_green(m)*surf%pt1(m) + &
                              f_qsws(m)*(qv1(m) - q_s(m) + dq_s_dt(m)*t_surf_green%val(m)) + &
                              surf%lambda_surf_green(m)*t_green%val(nzb_wall, m)

            coef_green_2(m) = 4.0_wp*surf%emissivity(m, ind_pav_green)*sigma_sb* &
                              t_surf_green%val(m)**3 + f_qsws(m)*dq_s_dt(m) + &
                              surf%lambda_surf_green(m) + f_shf_green(m)/exner(k)
         else
            coef_green_1(m) = surf%rad_net_l(m) + (3.0_wp + 1.0_wp)* &
                              surf%emissivity(m, ind_pav_green)*sigma_sb*t_surf_green%val(m)**4 + &
                              f_shf_green(m)*surf%pt1(m) + surf%lambda_surf_green(m)* &
                              t_green%val(nzb_wall, m)
            coef_green_2(m) = 4.0_wp*surf%emissivity(m, ind_pav_green)*sigma_sb* &
                              t_surf_green%val(m)**3 + surf%lambda_surf_green(m) + &
                              f_shf_green(m)/exner(k)
         end if
      end do

!$OMP PARALLEL DO PRIVATE (m, coef_1, coef_2, stend_wall, stend_window,                        &
!$OMP&                     stend_green ) SCHEDULE (STATIC)
      do m = 1, surf%ns
!
!--    Get k index of respective grid point.
         k = surf%k(m)
!
!--    Numerator of the prognostic equation.
!--    Todo: Adjust to tile approach. So far, emissivity for wall (element 0) is used
!--    Rem: Coef +1 corresponds to -lwout included in calculation of radnet_l.
         coef_1 = surf%rad_net_l(m) + (3.0_wp + 1.0_wp)* &
                  surf%emissivity(m, ind_veg_wall)*sigma_sb*t_surf_wall%val(m)**4 + &
                  f_shf(m)*surf%pt1(m) + surf%lambda_surf(m)*t_wall%val(nzb_wall, m)

!
!--    Denominator of the prognostic equation.
         coef_2 = 4.0_wp*surf%emissivity(m, ind_veg_wall)*sigma_sb*t_surf_wall%val(m)**3 + &
                  surf%lambda_surf(m) + f_shf(m)/exner(k)
!
!--    Implicit solution when the surface layer has no heat capacity, otherwise use RK3 scheme.
         t_surf_wall_p%val(m) = (coef_1*dt_3d*tsc(2) + surf%c_surface(m)* &
                                 t_surf_wall%val(m))/(surf%c_surface(m) + coef_2*dt_3d*tsc(2))

         if (.not. spinup_phase .and. frac_win(m) > 0.0_wp) then
            t_surf_window_p%val(m) = (coef_window_1(m)*dt_3d*tsc(2) + &
                                      surf%c_surface_window(m)*t_surf_window%val(m))/ &
                                     (surf%c_surface_window(m) + coef_window_2(m)*dt_3d*tsc(2))
         end if
         t_surf_green_p%val(m) = (coef_green_1(m)*dt_3d*tsc(2) + &
                                  surf%c_surface_green(m)*t_surf_green%val(m))/ &
                                 (surf%c_surface_green(m) + coef_green_2(m)*dt_3d*tsc(2))

!
!--    Add RK3 term.
         t_surf_wall_p%val(m) = t_surf_wall_p%val(m) + dt_3d*tsc(3)* &
                                surf%tt_surface_wall_m(m)
         t_surf_window_p%val(m) = t_surf_window_p%val(m) + dt_3d*tsc(3)* &
                                  surf%tt_surface_window_m(m)
         t_surf_green_p%val(m) = t_surf_green_p%val(m) + dt_3d*tsc(3)* &
                                 surf%tt_surface_green_m(m)

!
!--    Store surface temperature on pt_surface. Further, in case humidity is used, store also
!--    vpt_surface, which is, due to the lack of moisture on roofs, simply assumed to be the
!--    surface temperature.
         surf%pt_surface(m) = (frac_wall(m)*t_surf_wall_p%val(m) &
                               + frac_win(m)*t_surf_window_p%val(m) &
                               + frac_green(m)*t_surf_green_p%val(m) &
                               )/exner(k)

!
!--    Following line is actually not fully correct. In order to overcome this, a q_surface
!--    would be needed, calculated according to the q_surface in the LSM, where it is assumed
!--    that the skin layer is saturated. However, it is not clear whether this makes much sence
!--    in case of walls and windows. Probably only for green surfaces.
         if (humidity) surf%vpt_surface(m) = surf%pt_surface(m)
!
!--    Calculate true tendency
         stend_wall = (t_surf_wall_p%val(m) - t_surf_wall%val(m) - dt_3d*tsc(3)* &
                       surf%tt_surface_wall_m(m))/(dt_3d*tsc(2))
         stend_window = (t_surf_window_p%val(m) - t_surf_window%val(m) - dt_3d*tsc(3)* &
                         surf%tt_surface_window_m(m))/(dt_3d*tsc(2))
         stend_green = (t_surf_green_p%val(m) - t_surf_green%val(m) - dt_3d*tsc(3)* &
                        surf%tt_surface_green_m(m))/(dt_3d*tsc(2))
!
!--    Calculate t_surf tendencies for the next Runge-Kutta step.
         if (runge_l) then
            if (intermediate_timestep_count == 1) then
               surf%tt_surface_wall_m(m) = stend_wall
               surf%tt_surface_window_m(m) = stend_window
               surf%tt_surface_green_m(m) = stend_green
            elseif (intermediate_timestep_count < intermediate_timestep_count_max) then
               surf%tt_surface_wall_m(m) = -9.5625_wp*stend_wall + &
                                           5.3125_wp*surf%tt_surface_wall_m(m)
               surf%tt_surface_window_m(m) = -9.5625_wp*stend_window + &
                                             5.3125_wp*surf%tt_surface_window_m(m)
               surf%tt_surface_green_m(m) = -9.5625_wp*stend_green + &
                                            5.3125_wp*surf%tt_surface_green_m(m)
            end if
         end if
      end do

!$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
      do m = 1, surf%ns
!
!--    In case of fast changes in the skin temperature, it is required to update the radiative
!--    fluxes in order to keep the solution stable.
         if (((abs(t_surf_wall_p%val(m) - t_surf_wall%val(m)) > 1.0_wp) .or. &
              (abs(t_surf_green_p%val(m) - t_surf_green%val(m)) > 1.0_wp) .or. &
              (abs(t_surf_window_p%val(m) - t_surf_window%val(m)) > 1.0_wp)) &
             .and. unscheduled_radiation_calls) &
            then
            force_radiation_call_l_v(m) = .true.
         end if
      end do

!$OMP PARALLEL DO PRIVATE (m, tend, i, i_off, j, j_off, k, k_off ) SCHEDULE (STATIC)
      do m = 1, surf%ns
!
!--    Index offset of surface element point with respect to adjoining atmospheric grid point.
         k_off = surf%koff(m)
         j_off = surf%joff(m)
         i_off = surf%ioff(m)
!
!--    Get indices of respective grid point.
         i = surf%i(m)
         j = surf%j(m)
         k = surf%k(m)
!
!--    Calculate new fluxes
!--    rad_net_l is never used!
         surf%rad_net_l(m) = surf%rad_net_l(m) + &
                             frac_wall(m)*sigma_sb*surf%emissivity(m, ind_veg_wall)* &
                             (t_surf_wall_p%val(m)**4 - t_surf_wall%val(m)**4) + &
                             frac_win(m)*sigma_sb*surf%emissivity(m, ind_wat_win)* &
                             (t_surf_window_p%val(m)**4 - t_surf_window%val(m)**4) + &
                             frac_green(m)*sigma_sb*surf%emissivity(m, ind_pav_green)* &
                             (t_surf_green_p%val(m)**4 - t_surf_green%val(m)**4)

         surf%wghf_eb(m) = surf%lambda_surf(m)* &
                           (t_surf_wall_p%val(m) - t_wall%val(nzb_wall, m))
         surf%wghf_eb_green(m) = surf%lambda_surf_green(m)* &
                                 (t_surf_green_p%val(m) - t_green%val(nzb_wall, m))
         surf%wghf_eb_window(m) = surf%lambda_surf_window(m)* &
                                  (t_surf_window_p%val(m) - t_window%val(nzb_wall, m))
!
!--    Ground/wall/roof surface heat flux.
         surf%wshf_eb(m) = -f_shf(m)*(surf%pt1(m) - t_surf_wall_p%val(m)/exner(k))* &
                           frac_wall(m) &
                           - f_shf_window(m)*(surf%pt1(m) - t_surf_window_p%val(m)/exner(k))* &
                           frac_win(m) &
                           - f_shf_green(m)*(surf%pt1(m) - t_surf_green_p%val(m)/exner(k))* &
                           frac_green(m)
!
!--    Store kinematic surface heat fluxes for utilization in other processes diffusion_s,
!--    surface_layer_fluxes,...
         surf%shf(m) = surf%wshf_eb(m)/c_p
!
!--    If the indoor model is applied, further add waste heat from buildings to the kinematic flux.
         if (indoor_model) then
            surf%shf(m) = surf%shf(m) + surf%waste_heat(m)/c_p
         end if
!
!--    Following line is necessary to remove the density from the flux. For horizontal surfaces
!--    where the heat flux is added to the vertical diffusion term, density cancels out in
!--    diffusion_s.f90. However, for vertical surfaces the density would still be included in the
!--    diffusion terms, meaning that the heat-fluxes at the walls would be overestimated by about
!--    15-20%. Please note, here in the building-surface model density is expressed by
!--    hyp(k) / ( r_d * surf%pt1(m) * exner(k) ).
         if (.not. horizontal(m)) then
            surf%shf(m) = surf%shf(m)*(r_d*surf%pt1(m)*exner(k))/hyp(k)
         end if

         if (humidity .and. frac_green(m) > 0.0_wp) then
!
!--       Calculate true surface resistance.
            if (surf%upward(m)) then
               surf%qsws(m) = -f_qsws(m)*(qv1(m) - q_s(m) + &
                                          dq_s_dt(m)*t_surf_green%val(m) - &
                                          dq_s_dt(m)*t_surf_green_p%val(m))
               surf%qsws(m) = surf%qsws(m)/l_v
               surf%qsws_veg(m) = -f_qsws_veg(m)*(qv1(m) - q_s(m) + &
                                                  dq_s_dt(m)*t_surf_green%val(m) - &
                                                  dq_s_dt(m)*t_surf_green_p%val(m))
               surf%qsws_liq(m) = -f_qsws_liq(m)*(qv1(m) - q_s(m) + &
                                                  dq_s_dt(m)*t_surf_green%val(m) - &
                                                  dq_s_dt(m)*t_surf_green_p%val(m))
               surf%r_s(m) = -rho_lv(m)*(qv1(m) - q_s(m) + &
                                         dq_s_dt(m)*t_surf_green%val(m) - &
                                         dq_s_dt(m)*t_surf_green_p%val(m))/ &
                             (surf%qsws(m) + 1.0e-20) - surf%r_a_green(m)

               if (precipitation) then
!
!--             Calculate change in liquid water reservoir due to dew fall or evaporation of liquid
!--             water. If precipitation is activated, add rain water to qsws_liq and qsws_soil
!--             according to the vegetation coverage. Precipitation_rate is given in mm.
!--             Add precipitation to liquid water reservoir, if possible. Otherwise, add the water
!--             to soil. In case of pavements, the exceeding water amount is implicitely removed as
!--             runoff as qsws_soil is then not used in the soil model
                  if (m_liq_usm%val(m) /= m_liq_max(m)) then
                     surf%qsws_liq(m) = surf%qsws_liq(m) + frac_green(m)*rho_l*l_v*0.001_wp* &
                                        prr(k + k_off, j + j_off, i + i_off)*hyrho(k + k_off)

                  end if
               end if
!
!--          If the air is saturated, check the reservoir water level.
               if (surf%qsws(m) < 0.0_wp) then
!
!--             Check if reservoir is full (avoid values > m_liq_max) In that case, qsws_liq goes to
!--             qsws_soil. In this case qsws_veg is zero anyway (because c_liq = 1), so that tend is
!--             zero and no further check is needed.
                  if (m_liq_usm%val(m) == m_liq_max(m)) then
                     surf%qsws_liq(m) = 0.0_wp
                  end if
!
!--             In case qsws_veg becomes negative (unphysical behavior), let the water enter the
!--             liquid water reservoir as dew on the plant.
                  if (surf%qsws_veg(m) < 0.0_wp) then
                     surf%qsws_liq(m) = surf%qsws_liq(m) + surf%qsws_veg(m)
                     surf%qsws_veg(m) = 0.0_wp
                  end if
               end if

               tend = -surf%qsws_liq(m)*drho_l_lv(m)
               m_liq_usm_p%val(m) = m_liq_usm%val(m) + dt_3d* &
                                    (tsc(2)*tend + tsc(3)*tm_liq_usm_m%val(m))
!
!--          Check if reservoir is overfull -> reduce to maximum
!--          (conservation of water is violated here)
               m_liq_usm_p%val(m) = min(m_liq_usm_p%val(m), m_liq_max(m))
!
!--          Check if reservoir is empty (avoid values < 0.0) (conservation of water is
!--          violated here).
               m_liq_usm_p%val(m) = max(m_liq_usm_p%val(m), 0.0_wp)
!
!--          Calculate m_liq tendencies for the next Runge-Kutta step
               if (runge_l) then
                  if (intermediate_timestep_count == 1) then
                     tm_liq_usm_m%val(m) = tend
                  elseif (intermediate_timestep_count < intermediate_timestep_count_max) then
                     tm_liq_usm_m%val(m) = -9.5625_wp*tend + 5.3125_wp*tm_liq_usm_m%val(m)
                  end if
               end if
            else
!
!--          Downward or vertical surfaces.
               surf%qsws(m) = -f_qsws(m)*(qv1(m) - q_s(m) + dq_s_dt(m)*t_surf_green%val(m) &
                                          - dq_s_dt(m)*t_surf_green_p%val(m))
               surf%qsws(m) = surf%qsws(m)/l_v
!
!--          Following line is necessary to remove the density from the flux. For horizontal
!--          surfaceswhere the heat flux is added to the vertical diffusion term, density cancels
!--          out in diffusion_s.f90. However, for vertical surfaces the density would still be
!--          included in the diffusion terms, meaning that the heat-fluxes at the walls would be
!--          overestimated by about 15-20%. Please note, here in the building-surface model density
!--          is expressed by hyp(k) / ( r_d * surf%pt1(m) * exner(k) )
               if (.not. horizontal(m)) surf%qsws(m) = surf%qsws(m)* &
                                                       (r_d*surf%pt1(m)*exner(k))/hyp(k)

               surf%qsws_veg(m) = -f_qsws_veg(m)*(qv1(m) - q_s(m) + &
                                                  dq_s_dt(m)*t_surf_green%val(m) - &
                                                  dq_s_dt(m)*t_surf_green_p%val(m))
               surf%r_s(m) = -rho_lv(m)*(qv1(m) - q_s(m) + &
                                         dq_s_dt(m)*t_surf_green%val(m) - &
                                         dq_s_dt(m)*t_surf_green_p%val(m))/ &
                             (surf%qsws(m) + 1.0e-20_wp) - surf%r_a_green(m)
               surf%qsws_liq(m) = 0.0_wp  ! - f_qsws_liq(m)  * ( qv1(m) - q_s + dq_s_dt(m) * t_surf_green_h(m)&
!                                                  - dq_s_dt(m) * t_surf_green_h_p(m) )
!
!--          If the air is saturated, check the reservoir water level.
               if (surf%qsws(m) < 0.0_wp) then
!
!--             In case qsws_veg becomes negative (unphysical behavior), let the water enter the
!--             liquid water reservoir as dew on the plant.
                  if (surf%qsws_veg(m) < 0.0_wp) then
                     surf%qsws_veg(m) = 0.0_wp
                  end if
               end if
            end if
         else
            surf%r_s(m) = 1.0e10_wp
         end if

      end do
      force_radiation_call_l = any(force_radiation_call_l_v)
!
!-- Calculation of force_radiation_call:
!-- Make logical OR for all processes.
!-- Force radiation call if at least one processor forces it.
      if (intermediate_timestep_count == intermediate_timestep_count_max - 1) then
#if defined( __parallel )
         if (.not. force_radiation_call) then
            if (collective_wait) call MPI_BARRIER(comm2d, ierr)
            call MPI_ALLREDUCE(force_radiation_call_l, force_radiation_call, &
                               1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
         end if
#else
         force_radiation_call = (force_radiation_call .or. force_radiation_call_l)
#endif
         force_radiation_call_l = .false.
      end if

      if (debug_output_timestep) then
         write (debug_string, *) 'usm_surface_energy_balance: ', spinup_phase
         call debug_message(debug_string, 'end')
      end if

   end subroutine usm_surface_energy_balance

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of time levels for t_surf and t_wall called out from subroutine swap_timelevel
!--------------------------------------------------------------------------------------------------!
   subroutine usm_swap_timelevel(mod_count)

      implicit none

      integer(iwp), intent(IN)  ::  mod_count  !<

      select case (mod_count)

      case (0)
         t_surf_wall => t_surf_wall_1; t_surf_wall_p => t_surf_wall_2
         t_wall => t_wall_1; t_wall_p => t_wall_2
         t_surf_window => t_surf_window_1; t_surf_window_p => t_surf_window_2
         t_window => t_window_1; t_window_p => t_window_2
         t_surf_green => t_surf_green_1; t_surf_green_p => t_surf_green_2
         t_green => t_green_1; t_green_p => t_green_2
      case (1)
         t_surf_wall => t_surf_wall_2; t_surf_wall_p => t_surf_wall_1
         t_wall => t_wall_2; t_wall_p => t_wall_1
         t_surf_window => t_surf_window_2; t_surf_window_p => t_surf_window_1
         t_window => t_window_2; t_window_p => t_window_1
         t_surf_green => t_surf_green_2; t_surf_green_p => t_surf_green_1
         t_green => t_green_2; t_green_p => t_green_1

      end select

   end subroutine usm_swap_timelevel

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate maximum allowed timestep at USM surfaces according to the diffusion criterion.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_timestep

      integer(iwp) ::  k  !< grid index - wall depth
      integer(iwp) ::  m  !< running index for surface elements

      real(wp), dimension(nzb_wall:nzt_wall) ::  max_dt_green_column
      !< allowed timestep for each vertical level of a green surface column
      real(wp), dimension(nzb_wall:nzt_wall) ::  max_dt_wall_column   !< allowed timestep for each vertical level of a wall column
      real(wp), dimension(nzb_wall:nzt_wall) ::  max_dt_win_column    !< allowed timestep for each vertical level of a window column

      real(wp), dimension(:), allocatable ::  max_dt_usm  !< maximum allowed timestep for each usm surface

      allocate (max_dt_usm(1:surf_usm%ns))

!$OMP PARALLEL DO PRIVATE (k, m, max_dt_wall_column, max_dt_win_column, max_dt_green_column) SCHEDULE (STATIC)
      do m = 1, surf_usm%ns
         do k = nzb_wall, nzt_wall
            if (surf_usm%frac(m, ind_veg_wall) > 0.0_wp) then
               max_dt_wall_column(k) = surf_usm%rho_c_wall(k, m)/surf_usm%lambda_h_layer(k, m)* &
                                       (surf_usm%dz_wall(k, m))**2
            else
               max_dt_wall_column(k) = huge(1.0_wp)
            end if
            if (surf_usm%frac(m, ind_wat_win) > 0.0_wp) then
               max_dt_win_column(k) = surf_usm%rho_c_window(k, m)/ &
                                      surf_usm%lambda_h_window_layer(k, m)* &
                                      (surf_usm%dz_window(k, m))**2
            else
               max_dt_win_column(k) = huge(1.0_wp)
            end if
            if (surf_usm%frac(m, ind_pav_green) > 0.0_wp) then
               max_dt_green_column(k) = surf_usm%rho_c_total_green(k, m)/ &
                                        surf_usm%lambda_h_green(k, m)* &
                                        (surf_usm%dz_green(k, m))**2
            else
               max_dt_green_column(k) = huge(1.0_wp)
            end if
         end do
         max_dt_usm(m) = min(minval(max_dt_wall_column(:)), minval(max_dt_win_column(:)), &
                             minval(max_dt_green_column(:)))
      end do
!
!-- Consider a pre-factor (1/8) for the diffusion criterion.
      dt_usm = minval(max_dt_usm)*0.125_wp

#if defined( __parallel )
      if (collective_wait) call MPI_BARRIER(comm2d, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt_usm, 1, MPI_REAL, MPI_MIN, comm2d, ierr)
#endif

      deallocate (max_dt_usm)

   end subroutine usm_timestep

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sampling of USM variables along customized measurement coordinates.
!--------------------------------------------------------------------------------------------------!
   subroutine usm_vm_sampling(variable, var_atmos, i_atmos, j_atmos, k_atmos, ns_atmos, &
                              var_soil, i_soil, j_soil, k_soil, ns_soil, sampled)

      character(LEN=*) ::  variable  !< treated variable

      integer(iwp) ::  i         !< grid index in x-direction
      integer(iwp) ::  j         !< grid index in y-direction
      integer(iwp) ::  k         !< grid index in z-direction
      integer(iwp) ::  m         !< running index over all virtual observation coordinates
      integer(iwp) ::  mm        !< index of surface element which corresponds to the virtual observation coordinate
      integer(iwp) ::  ns_atmos  !< number of sampling points for atmosphere and surface variables
      integer(iwp) ::  ns_soil   !< number of sampling points for soil variables

      integer(iwp), dimension(1:ns_atmos) ::  i_atmos  !< sampling index in x-direction for atmosphere variables
      integer(iwp), dimension(1:ns_atmos) ::  j_atmos  !< sampling index in y-direction for atmosphere variables
      integer(iwp), dimension(1:ns_atmos) ::  k_atmos  !< sampling index in z-direction for atmosphere variables

      integer(iwp), dimension(1:ns_soil) ::  i_soil  !< sampling index in x-direction for soil variables
      integer(iwp), dimension(1:ns_soil) ::  j_soil  !< sampling index in y-direction for soil variables
      integer(iwp), dimension(1:ns_soil) ::  k_soil  !< sampling index in z-direction for soil variables

      logical ::  sampled !< flag indicating whether a variable has been sampled

      real(wp), dimension(1:ns_atmos) ::  var_atmos  !< array to store atmosphere variables

      real(wp), dimension(1:ns_soil) ::  var_soil  !< array to store soil variables

      select case (trim(variable))
!
!--    Soil and wall temperature.
      case ('t_soil')
         do m = 1, ns_soil
            if (j_soil(m) >= nys .and. j_soil(m) <= nyn .and. &
                i_soil(m) >= nxl .and. i_soil(m) <= nxr) &
               then
               k = k_soil(m)
               j = j_soil(m)
               i = i_soil(m)
!
!--             Take only values from horizontally-upward facing surfaces.
               do mm = surf_usm%start_index(j, i), surf_usm%end_index(j, i)
                  var_soil(m) = merge(t_wall%val(k, mm), var_soil(m), surf_usm%upward(mm))
               end do
            end if
         end do
         sampled = .true.

      case DEFAULT

      end select
!
!-- Avoid compiler warning for unused variables by constructing an if condition which is never
!-- fulfilled.
      if (.false. .and. ns_atmos < 0 .and. ns_soil < 0) then
         i_atmos = i_atmos
         j_atmos = j_atmos
         k_atmos = k_atmos
         var_atmos = var_atmos
      end if

   end subroutine usm_vm_sampling

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes t_surf and t_wall data into restart files
!--------------------------------------------------------------------------------------------------!
   subroutine usm_wrd_local

      implicit none

      integer(idp), dimension(nys:nyn, nxl:nxr) ::  global_end_index    !< end index for surface data (MPI-IO)
      integer(idp), dimension(nys:nyn, nxl:nxr) ::  global_start_index  !< start index for surface data (MPI-IO)

      logical  ::  surface_data_to_write  !< switch for MPI-I/O if PE has surface data to write

      if (trim(restart_data_format_output) == 'fortran_binary') then

         call wrd_write_string('ns_on_file_usm')
         write (14) surf_usm%ns

         call wrd_write_string('usm_start_index')
         write (14) surf_usm%start_index

         call wrd_write_string('usm_end_index')
         write (14) surf_usm%end_index

         call wrd_write_string('t_surf_wall')
         write (14) t_surf_wall%val

         call wrd_write_string('t_surf_window')
         write (14) t_surf_window%val

         call wrd_write_string('t_surf_green')
         write (14) t_surf_green%val

         call wrd_write_string('m_liq_usm')
         write (14) m_liq_usm%val

         call wrd_write_string('swc')
         write (14) swc%val

         call wrd_write_string('t_wall')
         write (14) t_wall%val

         call wrd_write_string('t_window')
         write (14) t_window%val

         call wrd_write_string('t_green')
         write (14) t_green%val

      elseif (restart_data_format_output(1:3) == 'mpi') then
!
!--    There is no information about the PE-grid necessary because the restart files consists of the
!--    whole domain. Therefore, ns_on_file_usm are not used with MPI-IO.
         call rd_mpi_io_surface_filetypes(surf_usm%start_index, surf_usm%end_index, &
                                          surface_data_to_write, global_start_index, &
                                          global_end_index)

         call wrd_mpi_io('usm_global_start', global_start_index)
         call wrd_mpi_io('usm_global_end', global_end_index)

         if (surface_data_to_write) then
            call wrd_mpi_io_surface('t_surf_wall', t_surf_wall%val)
            call wrd_mpi_io_surface('t_surf_window', t_surf_window%val)
            call wrd_mpi_io_surface('t_surf_green', t_surf_green%val)

            call wrd_mpi_io_surface('m_liq_usm', m_liq_usm%val)
         end if

         call rd_mpi_io_surface_filetypes(surf_usm%start_index, surf_usm%end_index, &
                                          surface_data_to_write, global_start_index, &
                                          global_end_index)

         call wrd_mpi_io('usm_global_start_2', global_start_index)
         call wrd_mpi_io('usm_global_end_2', global_end_index)

         if (surface_data_to_write) then
            call wrd_mpi_io_surface('swc', swc%val)
            call wrd_mpi_io_surface('t_wall', t_wall%val)
            call wrd_mpi_io_surface('t_window', t_window%val)
            call wrd_mpi_io_surface('t_green', t_green%val)
         end if

      end if

   end subroutine usm_wrd_local

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define building properties
!> Parameters 12, 13, 119 - 135 exclusive used in indoor_model_mod.f90
!> Parameters 0-11, 14-118, 136 - 149 exclusive used in urban_surface_mod.f90
!> Parameters 31, 44 used in indoor_model_mod.f90 and urban_surface_mod.f90
!--------------------------------------------------------------------------------------------------!
   subroutine usm_define_pars

!
!-- Define the building_pars
      building_pars(:, 1) = (/ &
                            0.82_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.18_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1520000.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
                            1512000.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            1512000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.93_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
                            0.81_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            0.81_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            299.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.91_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.7_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            2.9_wp, &  !< parameter 20  - [m] ground floor level height
                            0.82_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.18_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1520000.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
                            1512000.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            1512000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.93_wp, &
                            !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
                            0.81_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            0.81_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.91_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.7_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            36.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            37.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.2_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.38_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.4_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.18_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.36_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.42_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.45_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            1512000.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
                            1512000.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            2112000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            0.52_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
                            0.52_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            2.1_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.02_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.2_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.38_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.4_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            36.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.04_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.06_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.08_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.45_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.45_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.45_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            37.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.04_wp, &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
                            0.06_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.08_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.45_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.45_wp, &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.45_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.02_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.06_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.08_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.1_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            1512000.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            709650.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            709650.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.52_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
                            0.12_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            0.12_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.90_wp, &  !< parameter 100 - [-] wall emissivity roof
                            42.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.02_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.04_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.06_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.08_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.45_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.45_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.45_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.91_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.7_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            37.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.75_wp, &  !< parameter 119 - [-] shading factor
                            0.8_wp, &  !< parameter 120 - [-] g-value windows
                            2.9_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            0.5_wp, &
                            !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5
                            2.0_wp, &
                      !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0
                            0.0_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            3.0_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            260000.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            100.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            0.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            0.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            4.2_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            2.9_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            0.1_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            1.333_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1526000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
                            0.7_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
                            0.7_wp, &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
                            709650.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
                            0.12_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.45_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
                            0.45_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            0.7_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.45_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

      building_pars(:, 2) = (/ &
                            0.75_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.25_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1520000.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
                            79200.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            2112000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.93_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
                            0.046_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            2.1_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            299.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.87_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.65_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            2.5_wp, &  !< parameter 20  - [m] ground floor level height
                            0.75_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.25_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1520000.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
                            79200.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            2112000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.93_wp, &
                            !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
                            0.046_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            2.1_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.87_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.65_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            36.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            37.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.08_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.32_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.34_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.20_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.26_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.32_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.34_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            2112000.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
                            79200.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            2112000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            2.1_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
                            0.05_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            2.1_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.02_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.08_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.32_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.34_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            36.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.04_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.06_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.08_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.19_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.19_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.19_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            37.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.04_wp, &  !< parameter 80  - [m] 2nd cumulative window layer thickness above ground floor level
                            0.06_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.08_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.19_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.19_wp, &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.19_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.02_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.17_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.37_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.39_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            1700000.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            79200.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            2112000.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.16_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
                            0.046_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            2.1_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.93_wp, &  !< parameter 100 - [-] wall emissivity roof
                            42.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.02_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.04_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.06_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.08_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.19_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.19_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.19_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.87_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.65_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            37.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.75_wp, &  !< parameter 119 - [-] shading factor
                            0.7_wp, &  !< parameter 120 - [-] g-value windows
                            1.7_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            0.5_wp, &
                            !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5
                            1.5_wp, &
                      !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0
                            0.0_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            3.5_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            370000.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            80.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            0.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            0.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            4.2_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            2.5_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            0.0_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            2.54_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1526000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
                            0.7_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
                            0.7_wp, &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
                            357200.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
                            0.04_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.19_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
                            0.19_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            0.7_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.19_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

      building_pars(:, 3) = (/ &
                            0.71_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.29_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1520000.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
                            79200.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            1344000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.93_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
                            0.035_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            0.68_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            299.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.8_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.57_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            2.7_wp, &  !< parameter 20  - [m] ground floor level height
                            0.71_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.29_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1520000.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
                            79200.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            1344000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.93_wp, &
                            !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
                            0.035_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            0.68_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.8_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.57_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            36.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            38.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.22_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.58_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.6_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.20_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.32_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.38_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.41_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            2112000.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
                            79200.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            2112000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            2.1_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
                            0.05_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            2.1_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.02_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.22_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.58_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.6_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            36.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.03_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.06_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.09_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.12_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.11_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.11_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.11_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            38.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.03_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.06_wp, &  !< parameter 80  - [m] 2nd cumulative window layer thickness above ground floor level
                            0.09_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.12_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.11_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.11_wp, &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.11_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.02_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.06_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.36_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.38_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            3753600.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            709650.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            79200.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.52_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
                            0.12_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            0.035_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.93_wp, &  !< parameter 100 - [-] wall emissivity roof
                            42.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.03_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.06_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.09_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.12_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.11_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.11_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.11_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.8_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.57_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            38.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.15_wp, &  !< parameter 119 - [-] shading factor
                            0.6_wp, &  !< parameter 120 - [-] g-value windows
                            0.8_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            0.5_wp, &
                            !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5_wp
                            1.5_wp, &
                   !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0_wp
                            0.8_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            2.5_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            165000.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            40.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            0.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            0.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            4.2_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            2.7_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            -2.0_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            1.25_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1526000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
                            0.7_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
                            0.7_wp, &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
                            709650.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
                            0.12_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.11_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
                            0.11_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            0.7_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.11_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

      building_pars(:, 4) = (/ &
                            0.82_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.18_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1520000.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
                            1512000.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            1512000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.93_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
                            0.81_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            0.81_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            299.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.91_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.7_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            2.9_wp, &  !< parameter 20  - [m] ground floor level height
                            0.82_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.18_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1520000.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
                            1512000.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            1512000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.93_wp, &
                            !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
                            0.81_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            0.81_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.91_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.7_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            36.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            37.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.2_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.38_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.4_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.18_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.36_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.42_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.45_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            1512000.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
                            1512000.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            2112000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            0.52_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
                            0.52_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            2.1_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.02_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.2_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.38_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.4_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            36.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.04_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.06_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.08_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.45_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.45_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.45_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            37.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.04_wp, &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
                            0.06_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.08_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.45_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.45_wp, &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.45_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.02_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.06_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.08_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.1_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            1512000.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            709650.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            709650.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.52_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
                            0.12_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            0.12_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.90_wp, &  !< parameter 100 - [-] wall emissivity roof
                            42.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.02_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.04_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.06_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.08_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.45_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.45_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.45_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.91_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.7_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            37.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.75_wp, &  !< parameter 119 - [-] shading factor
                            0.8_wp, &  !< parameter 120 - [-] g-value windows
                            2.9_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            1.0_wp, &
                            !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 1.0_wp, winter 0.2
                            1.0_wp, &
                      !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.0_wp, winter 0.8
                            0.0_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            3.0_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            260000.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            100.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            0.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            7.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            3.0_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            2.9_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            0.1_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            1.333_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1526000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
                            0.7_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
                            0.7_wp, &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
                            709650.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
                            0.12_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.45_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
                            0.45_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            0.7_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.45_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

      building_pars(:, 5) = (/ &
                            0.75_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.25_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1520000.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
                            79200.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            2112000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.93_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
                            0.046_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            2.1_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            299.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.87_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.65_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            2.5_wp, &  !< parameter 20  - [m] ground floor level height
                            0.75_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.25_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1520000.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
                            79200.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            2112000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.93_wp, &
                            !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
                            0.046_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            2.1_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.87_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.65_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            36.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            37.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.08_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.32_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.34_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.20_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.26_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.32_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.34_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            2112000.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
                            79200.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            2112000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            2.1_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
                            0.05_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            2.1_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.02_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.08_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.32_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.34_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            36.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.04_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.06_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.08_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.19_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.19_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.19_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            37.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.04_wp, &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
                            0.06_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.08_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.19_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.19_wp, &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.19_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.02_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.17_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.37_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.39_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            1700000.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            79200.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            2112000.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.16_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
                            0.046_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            2.1_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.93_wp, &  !< parameter 100 - [-] wall emissivity roof
                            42.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.02_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.04_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.06_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.08_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.19_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.19_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.19_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.87_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.65_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            37.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.75_wp, &  !< parameter 119 - [-] shading factor
                            0.7_wp, &  !< parameter 120 - [-] g-value windows
                            1.7_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            1.0_wp, &
                            !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 1.0_wp, winter 0.2
                            1.0_wp, &
                      !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.0_wp, winter 0.8
                            0.0_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            3.5_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            370000.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            80.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            0.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            7.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            3.0_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            2.5_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            0.0_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            2.54_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1526000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
                            0.7_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
                            0.7_wp, &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
                            357200.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
                            0.04_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.19_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
                            0.19_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            0.7_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.19_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

      building_pars(:, 6) = (/ &
                            0.71_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.29_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1520000.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
                            79200.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            1344000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.93_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
                            0.035_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            0.68_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            299.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.8_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.57_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            2.7_wp, &  !< parameter 20  - [m] ground floor level height
                            0.71_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.29_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1520000.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
                            79200.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            1344000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.93_wp, &
                            !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
                            0.035_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            0.68_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.8_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.57_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            36.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            38.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.02_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.22_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.58_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.6_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.20_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.32_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.38_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.41_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            2112000.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
                            79200.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            2112000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            2.1_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
                            0.05_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            2.1_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.02_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.22_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.58_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.6_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            36.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.03_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.06_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.09_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.12_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.11_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.11_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.11_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            38.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.03_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.06_wp, &  !< parameter 80  - [m] 2nd cumulative window layer thickness above ground floor level
                            0.09_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.12_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.11_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.11_wp, &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.11_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.02_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.06_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.36_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.38_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            3753600.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            709650.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            79200.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.52_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
                            0.12_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            0.035_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.93_wp, &  !< parameter 100 - [-] wall emissivity roof
                            42.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.03_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.06_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.09_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.12_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.11_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.11_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.11_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.8_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.57_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            38.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.15_wp, &  !< parameter 119 - [-] shading factor
                            0.6_wp, &  !< parameter 120 - [-] g-value windows
                            0.8_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            1.0_wp, &
                            !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 1.0_wp, winter 0.2
                            1.0_wp, &
                      !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.0_wp, winter 0.8
                            0.8_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            2.5_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            165000.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            40.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            -80.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            7.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            3.0_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            2.7_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            -2.0_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            1.25_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1526000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
                            0.7_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
                            0.7_wp, &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
                            709650.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
                            0.12_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.11_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
                            0.11_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1526000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            0.7_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.11_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

      building_pars(:, 7) = (/ &
                            1.0_wp, &  !< parameter 0   - [-] wall fraction above ground floor level
                            0.0_wp, &  !< parameter 1   - [-] window fraction above ground floor level
                            0.0_wp, &  !< parameter 2   - [-] green fraction above ground floor level
                            0.0_wp, &  !< parameter 3   - [-] green fraction roof above ground floor level
                            1.5_wp, &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
                            1.5_wp, &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
                            1950400.0_wp, &
                            !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (upside) above ground floor level
                            1848000.0_wp, &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
                            1848000.0_wp, &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
                            0.7_wp, &
                            !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (upside) above ground floor level
                            1.0_wp, &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
                            1.0_wp, &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
                            372.15_wp, &  !< parameter 12  - [K] indoor target summer temperature
                            293.15_wp, &  !< parameter 13  - [K] indoor target winter temperature
                            0.93_wp, &  !< parameter 14  - [-] wall emissivity above ground floor level
                            0.86_wp, &  !< parameter 15  - [-] green emissivity above ground floor level
                            0.8_wp, &  !< parameter 16  - [-] window emissivity above ground floor level
                            0.7_wp, &
                            !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
                            0.001_wp, &  !< parameter 18  - [m] z0 roughness above ground floor level
                            0.0001_wp, &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
                            4.0_wp, &  !< parameter 20  - [m] ground floor level height
                            1.0_wp, &  !< parameter 21  - [-] wall fraction ground floor level
                            0.0_wp, &  !< parameter 22  - [-] window fraction ground floor level
                            0.0_wp, &  !< parameter 23  - [-] green fraction ground floor level
                            0.0_wp, &  !< parameter 24  - [-] green fraction roof ground floor level
                            1.5_wp, &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
                            1950400.0_wp, &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (upside) ground floor level
                            1848000.0_wp, &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
                            1848000.0_wp, &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
                            0.7_wp, &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (upside) ground floor level
                            1.0_wp, &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
                            1.0_wp, &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
                            0.93_wp, &  !< parameter 32  - [-] wall emissivity ground floor level
                            0.8_wp, &  !< parameter 33  - [-] window emissivity ground floor level
                            0.86_wp, &  !< parameter 34  - [-] green emissivity ground floor level
                            0.7_wp, &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
                            0.001_wp, &  !< parameter 36  - [m] z0 roughness ground floor level
                            0.0001_wp, &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
                            20.0_wp, &
                        !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
                            5.0_wp, &
                       !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
                            37.0_wp, &
                      !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
                            0.29_wp, &  !< parameter 41  - [m] 1st cumulative wall layer thickness above ground floor level
                            0.4_wp, &  !< parameter 42  - [m] 2nd cumulative wall layer thickness above ground floor level
                            0.695_wp, &  !< parameter 43  - [m] 3rd cumulative wall layer thickness above ground floor level
                            0.985_wp, &  !< parameter 44  - [m] 4th cumulative wall layer thickness above ground floor level
                            20000.0_wp, &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
                            23.0_wp, &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
                            20000.0_wp, &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
                            20000.0_wp, &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
                            23.0_wp, &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
                            10.0_wp, &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
                            1.0_wp, &  !< parameter 51  - [-] wall fraction ground plate
                            0.29_wp, &  !< parameter 52  - [m] 1st cumulative wall layer thickness ground plate
                            0.4_wp, &  !< parameter 53  - [m] 2nd cumulative wall layer thickness ground plate
                            0.695_wp, &  !< parameter 54  - [m] 3rd cumulative wall layer thickness ground plate
                            0.985_wp, &  !< parameter 55  - [m] 4th cumulative wall layer thickness ground plate
                            1950400.0_wp, &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (upside) ground plate
                            1848000.0_wp, &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
                            1848000.0_wp, &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
                            0.7_wp, &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (upside) ground plate
                            1.0_wp, &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
                            1.0_wp, &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
                            0.29_wp, &  !< parameter 62  - [m] 1st cumulative wall layer thickness ground floor level
                            0.4_wp, &  !< parameter 63  - [m] 2nd cumulative wall layer thickness ground floor level
                            0.695_wp, &  !< parameter 64  - [m] 3rd cumulative wall layer thickness ground floor level
                            0.985_wp, &  !< parameter 65  - [m] 4th cumulative wall layer thickness ground floor level
                            20.0_wp, &
                            !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
                            0.003_wp, &  !< parameter 67  - [m] 1st cumulative window layer thickness ground floor level
                            0.006_wp, &  !< parameter 68  - [m] 2nd cumulative window layer thickness ground floor level
                            0.012_wp, &  !< parameter 69  - [m] 3rd cumulative window layer thickness ground floor level
                            0.018_wp, &  !< parameter 70  - [m] 4th cumulative window layer thickness ground floor level
                            1736000.0_wp, &
                            !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
                            1736000.0_wp, &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
                            1736000.0_wp, &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
                            0.57_wp, &
                            !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
                            0.57_wp, &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
                            0.57_wp, &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
                            37.0_wp, &
                            !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
                            5.0_wp, &
                            !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
                            0.003_wp, &  !< parameter 79  - [m] 1st cumulative window layer thickness above ground floor level
                            0.006_wp, &  !< parameter 80  - [m] 2nd cumulative window layer thickness above ground floor level
                            0.012_wp, &  !< parameter 81  - [m] 3rd cumulative window layer thickness above ground floor level
                            0.018_wp, &  !< parameter 82  - [m] 4th cumulative window layer thickness above ground floor level
                            1736000.0_wp, &
                            !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
                            1736000.0_wp, &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
                            1736000.0_wp, &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
                            0.57_wp, &
                            !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
                            0.57_wp, &  !< parameter 87  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
                            0.57_wp, &  !< parameter 88  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
                            1.0_wp, &  !< parameter 89  - [-] wall fraction roof
                            0.29_wp, &  !< parameter 90  - [m] 1st cumulative wall layer thickness roof
                            0.4_wp, &  !< parameter 91  - [m] 2nd cumulative wall layer thickness roof
                            0.695_wp, &  !< parameter 92  - [m] 3rd cumulative wall layer thickness roof
                            0.985_wp, &  !< parameter 93  - [m] 4th cumulative wall layer thickness roof
                            1950400.0_wp, &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
                            1848000.0_wp, &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
                            1848000.0_wp, &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
                            0.7_wp, &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (upside) roof
                            1.0_wp, &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
                            1.0_wp, &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
                            0.93_wp, &  !< parameter 100 - [-] wall emissivity roof
                            19.0_wp, &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 102 - [-] window fraction roof
                            0.003_wp, &  !< parameter 103 - [m] window 1st layer thickness roof
                            0.006_wp, &  !< parameter 104 - [m] window 2nd layer thickness roof
                            0.012_wp, &  !< parameter 105 - [m] window 3rd layer thickness roof
                            0.018_wp, &  !< parameter 106 - [m] window 4th layer thickness roof
                            1736000.0_wp, &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
                            1736000.0_wp, &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
                            1736000.0_wp, &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
                            0.57_wp, &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
                            0.57_wp, &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
                            0.57_wp, &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
                            0.8_wp, &  !< parameter 113 - [-] window emissivity roof
                            0.7_wp, &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
                            37.0_wp, &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
                            0.86_wp, &  !< parameter 116 - [-] green emissivity roof
                            5.0_wp, &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
                            0.0_wp, &  !< parameter 118 - [-] green type roof
                            0.8_wp, &  !< parameter 119 - [-] shading factor
                            100.0_wp, &  !< parameter 120 - [-] g-value windows
                            100.0_wp, &  !< parameter 121 - [W/(m2*K)] u-value windows
                            20.0_wp, &  !< parameter 122 - [1/h] basic airflow without occupancy of the room
                            20.0_wp, &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room
                            0.0_wp, &  !< parameter 124 - [-] heat recovery efficiency
                            1.0_wp, &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
                            1.0_wp, &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heatstorage
                            4.5_wp, &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
                            100000.0_wp, &  !< parameter 128 - [W] maximal heating capacity
                            0.0_wp, &  !< parameter 129 - [W] maximal cooling capacity
                            0.0_wp, &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
                            0.0_wp, &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
                            3.0_wp, &  !< parameter 132 - [m] storey height
                            0.2_wp, &  !< parameter 133 - [m] ceiling construction height
                            0.0_wp, &  !< parameter 134 - [-] anthropogenic heat output for heating
                            0.0_wp, &  !< parameter 135 - [-] anthropogenic heat output for cooling
                            1848000.0_wp, &
                            !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (downside) above ground floor level
                            1.0_wp, &
                            !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (downside) above ground floor level
                            1848000.0_wp, &
                            !< parameter 138 - [J/(m3*K)] heat capacity 4th wall layer (downside) ground floor level
                            1.0_wp, &
                            !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (downside) ground floor level
                            1848000.0_wp, &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (downside) ground plate
                            1.0_wp, &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (downside) ground plate
                            1736000.0_wp, &
                            !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
                            0.57_wp, &
                            !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
                            1736000.0_wp, &
                            !< parameter 144 - [J/(m3*K)] heat capacity 4th window layer (inside) above ground floor level
                            0.57_wp, &
                            !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
                            1848000.0_wp, &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
                            1.0_wp, &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (downside) roof
                            1736000.0_wp, &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
                            0.57_wp &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                            /)

!
!-- Define building properties. gfl: ground floor level, agfl: above ground floor level
!-- Building fractions.
!-- Type:                            1        2         3       4        5        6        7
      building_frac(ind_wall_gfl, 1:7) = &
         (/0.82_wp, 0.75_wp, 0.71_wp, 0.82_wp, 0.75_wp, 0.71_wp, 1.0_wp/)
      building_frac(ind_wall_agfl, 1:7) = &
         (/0.82_wp, 0.75_wp, 0.71_wp, 0.82_wp, 0.75_wp, 0.71_wp, 1.0_wp/)
      building_frac(ind_wall_roof, 1:7) = &
         (/1.00_wp, 1.00_wp, 1.00_wp, 1.00_wp, 1.00_wp, 1.00_wp, 1.0_wp/)
      building_frac(ind_win_gfl, 1:7) = &
         (/0.18_wp, 0.25_wp, 0.29_wp, 0.18_wp, 0.25_wp, 0.29_wp, 0.0_wp/)
      building_frac(ind_win_agfl, 1:7) = &
         (/0.18_wp, 0.25_wp, 0.29_wp, 0.18_wp, 0.25_wp, 0.29_wp, 0.0_wp/)
      building_frac(ind_win_roof, 1:7) = &
         (/0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.0_wp/)
      building_frac(ind_green_gfl, 1:7) = &
         (/0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.0_wp/)
      building_frac(ind_green_agfl, 1:7) = &
         (/0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.0_wp/)
      building_frac(ind_green_roof, 1:7) = &
         (/0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.0_wp/)

!
!-- General building parameters.
!-- Type:             1          2          3          4          5          6          7
      building_gen(ind_gflh, 1:7) = &
         (/2.9000_wp, 2.5000_wp, 2.7000_wp, 2.9000_wp, 2.5000_wp, 2.7000_wp, 4.0000_wp/)
      building_gen(ind_green_type_roof, 1:7) = &
         (/1.0000_wp, 1.0000_wp, 1.0000_wp, 1.0000_wp, 1.0000_wp, 1.0000_wp, 1.0000_wp/)

!
!-- Building LAI.
!-- Type:             1          2          3          4          5          6          7
      building_lai(ind_gfl, 1:7) = &
         (/1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp/)
      building_lai(ind_agfl, 1:7) = &
         (/1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp/)
      building_lai(ind_roof, 1:7) = &
         (/1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp, 1.5000_wp/)

!
!-- Building z0.
!-- Type:             1          2          3          4          5          6          7
      building_z0(ind_gfl, 1:7) = &
         (/0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp/)
      building_z0(ind_agfl, 1:7) = &
         (/0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp/)
      building_z0(ind_roof, 1:7) = &
         (/0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp, 0.0010_wp/)

!
!-- Building z0qh.
!-- Type:             1          2          3          4          5          6          7
      building_z0qh(ind_gfl, 1:7) = &
         (/0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp/)
      building_z0qh(ind_agfl, 1:7) = &
         (/0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp/)
      building_z0qh(ind_roof, 1:7) = &
         (/0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp, 0.0001_wp/)

!
!-- Indoor building parameters.
!-- 1:   indoor target summer temperature
!-- 2:   indoor target winter temperature
!-- 3:   shading factor
!-- 4:   g-value windows
!-- 5:   u-value windows
!-- 6:   basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5
!-- 7:   additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0
!-- 8:   heat recovery efficiency
!-- 9:   dynamic parameter specific effective surface
!-- 10:   dynamic parameter innner heat storage
!-- 11:  ratio internal surface/floor area
!-- 12:  maximal heating capacity
!-- 13:  maximal cooling capacity
!-- 14:  additional internal heat gains dependent on occupancy of the room
!-- 15:  basic internal heat gains without occupancy of the room
!-- 16:  storey height
!-- 17:  ceiling construction height
!-- 18:  anthropogenic heat output for heating
!-- 19:  anthropogenic heat output for cooling

!-- Type: 1            2            3            4            5            6            7
      building_indoor(ind_theta_int_c_set, 1:7) = &
         (/299.1500_wp, 299.1500_wp, 299.1500_wp, 299.1500_wp, 299.1500_wp, 299.1500_wp, 372.1500_wp/)
      building_indoor(ind_theta_int_h_set, 1:7) = &
         (/293.1500_wp, 293.1500_wp, 293.1500_wp, 293.1500_wp, 293.1500_wp, 293.1500_wp, 293.1500_wp/)
      building_indoor(ind_f_c_win, 1:7) = &
         (/0.7500_wp, 0.7500_wp, 0.1500_wp, 0.7500_wp, 0.7500_wp, 0.1500_wp, 0.8000_wp/)
      building_indoor(ind_g_value_win, 1:7) = &
         (/0.8000_wp, 0.7000_wp, 0.6000_wp, 0.8000_wp, 0.7000_wp, 0.6000_wp, 100.0000_wp/)
      building_indoor(ind_u_value_win, 1:7) = &
         (/2.9000_wp, 1.7000_wp, 0.8000_wp, 2.9000_wp, 1.7000_wp, 0.8000_wp, 100.0000_wp/)
      building_indoor(ind_airflow_unocc, 1:7) = &
         (/0.5000_wp, 0.5000_wp, 0.5000_wp, 1.0000_wp, 1.0000_wp, 1.0000_wp, 20.0000_wp/)
      building_indoor(ind_airflow_occ, 1:7) = &
         (/2.0000_wp, 1.5000_wp, 1.5000_wp, 1.0000_wp, 1.0000_wp, 1.0000_wp, 20.0000_wp/)
      building_indoor(ind_eta_ve, 1:7) = &
         (/0.0000_wp, 0.0000_wp, 0.8000_wp, 0.0000_wp, 0.0000_wp, 0.8000_wp, 0.0000_wp/)
      building_indoor(ind_factor_a, 1:7) = &
         (/3.0000_wp, 3.5000_wp, 2.5000_wp, 3.0000_wp, 3.5000_wp, 2.5000_wp, 1.0000_wp/)
      building_indoor(ind_factor_c, 1:7) = &
         (/260000.0_wp, 370000.0_wp, 165000.0_wp, 260000.0_wp, 370000.0_wp, 165000.0_wp, 1.0_wp/)
      building_indoor(ind_lambda_at, 1:7) = &
         (/4.5000_wp, 4.5000_wp, 4.5000_wp, 4.5000_wp, 4.5000_wp, 4.5000_wp, 4.5000_wp/)
      building_indoor(ind_q_h_max, 1:7) = &
         (/100.0000_wp, 80.0000_wp, 40.0000_wp, 100.0000_wp, 80.0000_wp, 40.0000_wp, 0.0000_wp/)
      building_indoor(ind_q_c_max, 1:7) = &
         (/0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, -80.0000_wp, 0.0000_wp/)
      building_indoor(ind_qint_high, 1:7) = &
         (/0.0000_wp, 0.0000_wp, 0.0000_wp, 7.0000_wp, 7.0000_wp, 7.0000_wp, 0.0000_wp/)
      building_indoor(ind_qint_low, 1:7) = &
         (/4.2000_wp, 4.2000_wp, 4.2000_wp, 3.0000_wp, 3.0000_wp, 3.0000_wp, 0.0000_wp/)
      building_indoor(ind_height_storey, 1:7) = &
         (/2.9000_wp, 2.5000_wp, 2.7000_wp, 2.9000_wp, 2.5000_wp, 2.7000_wp, 3.0000_wp/)
      building_indoor(ind_height_cei_con, 1:7) = &
         (/0.2000_wp, 0.2000_wp, 0.2000_wp, 0.2000_wp, 0.2000_wp, 0.2000_wp, 0.2000_wp/)
      building_indoor(ind_params_waste_heat_h, 1:7) = &
         (/0.1000_wp, 0.0000_wp, -2.0000_wp, 0.1000_wp, 0.0000_wp, -2.0000_wp, 0.0000_wp/)
      building_indoor(ind_params_waste_heat_c, 1:7) = &
         (/1.3330_wp, 2.5400_wp, 1.2500_wp, 1.3330_wp, 2.5400_wp, 1.2500_wp, 0.0000_wp/)

!
!-- Window transmissivity.
!-- Type:                     1         2         3         4         5         6         7
      building_trans(ind_gfl, 1:7) = &
         (/0.70_wp, 0.65_wp, 0.57_wp, 0.70_wp, 0.65_wp, 0.57_wp, 0.70_wp/)
      building_trans(ind_agfl, 1:7) = &
         (/0.70_wp, 0.65_wp, 0.57_wp, 0.70_wp, 0.65_wp, 0.57_wp, 0.70_wp/)
      building_trans(ind_roof, 1:7) = &
         (/0.70_wp, 0.65_wp, 0.57_wp, 0.70_wp, 0.65_wp, 0.57_wp, 0.70_wp/)

!
!-- Building albedo types.
!-- Type:                                      1   2   3   4   5   6   7
      building_alb_type(ind_wall_gfl, 1:7) = (/36, 36, 36, 36, 36, 36, 20/)
      building_alb_type(ind_wall_agfl, 1:7) = (/36, 36, 36, 36, 36, 36, 20/)
      building_alb_type(ind_wall_roof, 1:7) = (/42, 42, 42, 42, 42, 42, 19/)
      building_alb_type(ind_win_gfl, 1:7) = (/37, 37, 38, 37, 37, 38, 37/)
      building_alb_type(ind_win_agfl, 1:7) = (/37, 37, 38, 37, 37, 38, 37/)
      building_alb_type(ind_win_roof, 1:7) = (/37, 37, 38, 37, 37, 38, 37/)
      building_alb_type(ind_green_gfl, 1:7) = (/5, 5, 5, 5, 5, 5, 5/)
      building_alb_type(ind_green_agfl, 1:7) = (/5, 5, 5, 5, 5, 5, 5/)
      building_alb_type(ind_green_roof, 1:7) = (/5, 5, 5, 5, 5, 5, 5/)

!
!-- Building emissivities.
!-- Type:                     1         2         3         4         5         6         7
      building_emis(ind_wall_gfl, 1:7) = &
         (/0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp/)
      building_emis(ind_wall_agfl, 1:7) = &
         (/0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp, 0.93_wp/)
      building_emis(ind_wall_roof, 1:7) = &
         (/0.90_wp, 0.93_wp, 0.93_wp, 0.90_wp, 0.93_wp, 0.93_wp, 0.93_wp/)
      building_emis(ind_win_gfl, 1:7) = &
         (/0.91_wp, 0.87_wp, 0.80_wp, 0.91_wp, 0.87_wp, 0.80_wp, 0.80_wp/)
      building_emis(ind_win_agfl, 1:7) = &
         (/0.91_wp, 0.87_wp, 0.80_wp, 0.91_wp, 0.87_wp, 0.80_wp, 0.80_wp/)
      building_emis(ind_win_roof, 1:7) = &
         (/0.91_wp, 0.87_wp, 0.80_wp, 0.91_wp, 0.87_wp, 0.80_wp, 0.80_wp/)
      building_emis(ind_green_gfl, 1:7) = &
         (/0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp/)
      building_emis(ind_green_agfl, 1:7) = &
         (/0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp/)
      building_emis(ind_green_roof, 1:7) = &
         (/0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp, 0.86_wp/)

!
!-- Building heat capacities.
!-- Type: 1           2           3           4           5            6           7
      building_hcap(ind_wall_gfl, 0, 1:7) = &
         (/1520.00_wp, 1520.00_wp, 1520.00_wp, 1520.00_wp, 1520.000_wp, 1520.00_wp, 1950.40_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_gfl, 1, 1:7) = &
         (/1512.00_wp, 79.20_wp, 79.20_wp, 1512.00_wp, 79.200_wp, 79.20_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_gfl, 2, 1:7) = &
         (/1512.00_wp, 2112.00_wp, 1344.00_wp, 1512.00_wp, 2112.000_wp, 1344.00_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_gfl, 3, 1:7) = &
         (/1526.00_wp, 1526.00_wp, 1526.00_wp, 1526.00_wp, 1526.000_wp, 1526.00_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_agfl, 0, 1:7) = &
         (/1520.00_wp, 1520.00_wp, 1520.00_wp, 1520.00_wp, 1520.000_wp, 1520.00_wp, 1950.40_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_agfl, 1, 1:7) = &
         (/1512.00_wp, 79.20_wp, 79.20_wp, 1512.00_wp, 79.200_wp, 79.20_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_agfl, 2, 1:7) = &
         (/1512.00_wp, 2112.00_wp, 1344.00_wp, 1512.00_wp, 2112.000_wp, 1344.00_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_agfl, 3, 1:7) = &
         (/1526.00_wp, 1526.00_wp, 1526.00_wp, 1526.00_wp, 1526.000_wp, 1526.00_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_roof, 0, 1:7) = &
         (/1512.00_wp, 1700.00_wp, 3753.60_wp, 1512.00_wp, 1700.000_wp, 3753.60_wp, 1950.40_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_roof, 1, 1:7) = &
         (/709.65_wp, 79.20_wp, 709.65_wp, 709.65_wp, 79.200_wp, 709.65_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_roof, 2, 1:7) = &
         (/709.65_wp, 2112.00_wp, 79.20_wp, 709.65_wp, 2112.000_wp, 79.20_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_wall_roof, 3, 1:7) = &
         (/1526.00_wp, 1526.00_wp, 1526.00_wp, 1526.00_wp, 1526.000_wp, 1526.00_wp, 1848.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_gfl, 0, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_gfl, 1, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_gfl, 2, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_gfl, 3, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_agfl, 0, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_agfl, 1, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_agfl, 2, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_agfl, 3, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_roof, 0, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_roof, 1, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_roof, 2, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_win_roof, 3, 1:7) = &
         (/1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.00_wp, 1736.000_wp, 1736.00_wp, 1736.00_wp/)* &
         1000.0_wp
      building_hcap(ind_green_gfl, 0, 1:7) = rho_c_soil
      building_hcap(ind_green_gfl, 1, 1:7) = rho_c_soil
      building_hcap(ind_green_gfl, 2, 1:7) = rho_c_soil
      building_hcap(ind_green_gfl, 3, 1:7) = rho_c_soil
      building_hcap(ind_green_agfl, 0, 1:7) = rho_c_soil
      building_hcap(ind_green_agfl, 1, 1:7) = rho_c_soil
      building_hcap(ind_green_agfl, 2, 1:7) = rho_c_soil
      building_hcap(ind_green_agfl, 3, 1:7) = rho_c_soil
      building_hcap(ind_green_roof, 0, 1:7) = rho_c_soil
      building_hcap(ind_green_roof, 1, 1:7) = rho_c_soil
      building_hcap(ind_green_roof, 2, 1:7) = rho_c_soil
      building_hcap(ind_green_roof, 3, 1:7) = rho_c_soil

!
!-- Building heat conductivities.
!-- Type:                       1        2         3         4        5         6         7
      building_hcond(ind_wall_gfl, 0, 1:7) = &
         (/0.93_wp, 0.930_wp, 0.930_wp, 0.93_wp, 0.930_wp, 0.930_wp, 0.70_wp/)
      building_hcond(ind_wall_gfl, 1, 1:7) = &
         (/0.81_wp, 0.046_wp, 0.035_wp, 0.81_wp, 0.046_wp, 0.035_wp, 1.00_wp/)
      building_hcond(ind_wall_gfl, 2, 1:7) = &
         (/0.81_wp, 2.100_wp, 0.680_wp, 0.81_wp, 2.100_wp, 0.680_wp, 1.00_wp/)
      building_hcond(ind_wall_gfl, 3, 1:7) = &
         (/0.70_wp, 0.700_wp, 0.700_wp, 0.70_wp, 0.700_wp, 0.700_wp, 1.00_wp/)
      building_hcond(ind_wall_agfl, 0, 1:7) = &
         (/0.93_wp, 0.930_wp, 0.930_wp, 0.93_wp, 0.930_wp, 0.930_wp, 0.70_wp/)
      building_hcond(ind_wall_agfl, 1, 1:7) = &
         (/0.81_wp, 0.046_wp, 0.035_wp, 0.81_wp, 0.046_wp, 0.035_wp, 1.00_wp/)
      building_hcond(ind_wall_agfl, 2, 1:7) = &
         (/0.81_wp, 2.100_wp, 0.680_wp, 0.81_wp, 2.100_wp, 0.680_wp, 1.00_wp/)
      building_hcond(ind_wall_agfl, 3, 1:7) = &
         (/0.70_wp, 0.700_wp, 0.700_wp, 0.70_wp, 0.700_wp, 0.700_wp, 1.00_wp/)
      building_hcond(ind_wall_roof, 0, 1:7) = &
         (/0.52_wp, 0.160_wp, 0.520_wp, 0.52_wp, 0.160_wp, 0.520_wp, 0.70_wp/)
      building_hcond(ind_wall_roof, 1, 1:7) = &
         (/0.12_wp, 0.046_wp, 0.120_wp, 0.12_wp, 0.046_wp, 0.120_wp, 1.00_wp/)
      building_hcond(ind_wall_roof, 2, 1:7) = &
         (/0.12_wp, 2.100_wp, 0.035_wp, 0.12_wp, 2.100_wp, 0.035_wp, 1.00_wp/)
      building_hcond(ind_wall_roof, 3, 1:7) = &
         (/0.70_wp, 0.700_wp, 0.700_wp, 0.70_wp, 0.700_wp, 0.700_wp, 1.00_wp/)
      building_hcond(ind_win_gfl, 0, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_gfl, 1, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_gfl, 2, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_gfl, 3, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_agfl, 0, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_agfl, 1, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_agfl, 2, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_agfl, 3, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_roof, 0, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_roof, 1, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_roof, 2, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_win_roof, 3, 1:7) = &
         (/0.45_wp, 0.190_wp, 0.110_wp, 0.45_wp, 0.190_wp, 0.110_wp, 0.57_wp/)
      building_hcond(ind_green_gfl, 0, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_gfl, 1, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_gfl, 2, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_gfl, 3, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_agfl, 0, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_agfl, 1, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_agfl, 2, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_agfl, 3, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_roof, 0, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_roof, 1, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_roof, 2, 1:7) = lambda_h_green_sm
      building_hcond(ind_green_roof, 3, 1:7) = lambda_h_green_sm

!
!-- Building layer depths.
!-- Type:             1          2          3          4          5          6          7
      building_depth(ind_wall_gfl, 0, 1:7) = &
         (/0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.2900_wp/)
      building_depth(ind_wall_gfl, 1, 1:7) = &
         (/0.2000_wp, 0.0800_wp, 0.2200_wp, 0.2000_wp, 0.0800_wp, 0.2200_wp, 0.4000_wp/)
      building_depth(ind_wall_gfl, 2, 1:7) = &
         (/0.3800_wp, 0.3200_wp, 0.5800_wp, 0.3800_wp, 0.3200_wp, 0.5800_wp, 0.6950_wp/)
      building_depth(ind_wall_gfl, 3, 1:7) = &
         (/0.4000_wp, 0.3400_wp, 0.6000_wp, 0.4000_wp, 0.3400_wp, 0.6000_wp, 0.9850_wp/)
      building_depth(ind_wall_agfl, 0, 1:7) = &
         (/0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.2900_wp/)
      building_depth(ind_wall_agfl, 1, 1:7) = &
         (/0.2000_wp, 0.0800_wp, 0.2200_wp, 0.2000_wp, 0.0800_wp, 0.2200_wp, 0.4000_wp/)
      building_depth(ind_wall_agfl, 2, 1:7) = &
         (/0.3800_wp, 0.3200_wp, 0.5800_wp, 0.3800_wp, 0.3200_wp, 0.5800_wp, 0.6950_wp/)
      building_depth(ind_wall_agfl, 3, 1:7) = &
         (/0.4000_wp, 0.3400_wp, 0.6000_wp, 0.4000_wp, 0.3400_wp, 0.6000_wp, 0.9850_wp/)
      building_depth(ind_wall_roof, 0, 1:7) = &
         (/0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.0200_wp, 0.2900_wp/)
      building_depth(ind_wall_roof, 1, 1:7) = &
         (/0.0600_wp, 0.1700_wp, 0.0600_wp, 0.0600_wp, 0.1700_wp, 0.0600_wp, 0.4000_wp/)
      building_depth(ind_wall_roof, 2, 1:7) = &
         (/0.0800_wp, 0.3700_wp, 0.3600_wp, 0.0800_wp, 0.3700_wp, 0.3600_wp, 0.6950_wp/)
      building_depth(ind_wall_roof, 3, 1:7) = &
         (/0.1000_wp, 0.3900_wp, 0.3800_wp, 0.1000_wp, 0.3900_wp, 0.3800_wp, 0.9850_wp/)
      building_depth(ind_win_gfl, 0, 1:7) = &
         (/0.0200_wp, 0.0200_wp, 0.0300_wp, 0.0200_wp, 0.0200_wp, 0.0300_wp, 0.0030_wp/)
      building_depth(ind_win_gfl, 1, 1:7) = &
         (/0.0400_wp, 0.0400_wp, 0.0600_wp, 0.0400_wp, 0.0400_wp, 0.0600_wp, 0.0060_wp/)
      building_depth(ind_win_gfl, 2, 1:7) = &
         (/0.0600_wp, 0.0600_wp, 0.0900_wp, 0.0600_wp, 0.0600_wp, 0.0900_wp, 0.0120_wp/)
      building_depth(ind_win_gfl, 3, 1:7) = &
         (/0.0800_wp, 0.0800_wp, 0.1200_wp, 0.0800_wp, 0.0800_wp, 0.1200_wp, 0.0180_wp/)
      building_depth(ind_win_agfl, 0, 1:7) = &
         (/0.0200_wp, 0.0200_wp, 0.0300_wp, 0.0200_wp, 0.0200_wp, 0.0300_wp, 0.0030_wp/)
      building_depth(ind_win_agfl, 1, 1:7) = &
         (/0.0400_wp, 0.0400_wp, 0.0600_wp, 0.0400_wp, 0.0400_wp, 0.0600_wp, 0.0060_wp/)
      building_depth(ind_win_agfl, 2, 1:7) = &
         (/0.0600_wp, 0.0600_wp, 0.0900_wp, 0.0600_wp, 0.0600_wp, 0.0900_wp, 0.0120_wp/)
      building_depth(ind_win_agfl, 3, 1:7) = &
         (/0.0800_wp, 0.0800_wp, 0.1200_wp, 0.0800_wp, 0.0800_wp, 0.1200_wp, 0.0180_wp/)
      building_depth(ind_win_roof, 0, 1:7) = &
         (/0.0200_wp, 0.0200_wp, 0.0300_wp, 0.0200_wp, 0.0200_wp, 0.0300_wp, 0.0030_wp/)
      building_depth(ind_win_roof, 1, 1:7) = &
         (/0.0400_wp, 0.0400_wp, 0.0600_wp, 0.0400_wp, 0.0400_wp, 0.0600_wp, 0.0060_wp/)
      building_depth(ind_win_roof, 2, 1:7) = &
         (/0.0600_wp, 0.0600_wp, 0.0900_wp, 0.0600_wp, 0.0600_wp, 0.0900_wp, 0.0120_wp/)
      building_depth(ind_win_roof, 3, 1:7) = &
         (/0.0800_wp, 0.0800_wp, 0.1200_wp, 0.0800_wp, 0.0800_wp, 0.1200_wp, 0.0180_wp/)
!
!-- Use the same depths for green roofs as for walls; might be a bad choice.
      building_depth(ind_green_gfl, :, :) = building_depth(ind_wall_gfl, :, :)
      building_depth(ind_green_agfl, :, :) = building_depth(ind_wall_agfl, :, :)
      building_depth(ind_green_roof, :, :) = building_depth(ind_wall_roof, :, :)

   end subroutine usm_define_pars

end module urban_surface_mod
