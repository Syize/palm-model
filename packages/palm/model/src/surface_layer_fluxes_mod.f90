!> @file surface_layer_fluxes_mod.f90
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
!> Diagnostic computation of vertical fluxes in the constant flux layer from the values of the
!> variables at grid point k=1 based on Newton iteration.
!>
!> @todo (Re)move large_scale_forcing actions
!> @todo Check/optimize OpenMP directives
!> @todo Simplify if conditions (which flux need to be computed in which case)
!--------------------------------------------------------------------------------------------------!
 MODULE surface_layer_fluxes_mod

    USE arrays_3d,                                                                                 &
        ONLY:  d_exner,                                                                            &
               drho_air_zw,                                                                        &
               e,                                                                                  &
               exner,                                                                              &
               nc,                                                                                 &
               nr,                                                                                 &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               qc,                                                                                 &
               qr,                                                                                 &
               s,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               vpt,                                                                                &
               w,                                                                                  &
               rho_air_zw


    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  c_p,                                                                                &
               g,                                                                                  &
               kappa,                                                                              &
               magnus,                                                                             &
               lv_d_cp,                                                                            &
               pi,                                                                                 &
               rd_d_rv

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    USE chem_modules,                                                                              &
        ONLY:  constant_csflux

    USE cpulog

    USE control_parameters,                                                                        &
        ONLY:  air_chemistry,                                                                      &
               atmosphere_run_coupled_to_ocean,                                                    &
               cloud_droplets,                                                                     &
               constant_heatflux,                                                                  &
               constant_scalarflux,                                                                &
               constant_waterflux,                                                                 &
               cut_cell_topography,                                                                &
               debug_output_timestep,                                                              &
               det_enabled,                                                                        &
               humidity,                                                                           &
               ibc_e_b,                                                                            &
               ibc_pt_b,                                                                           &
               indoor_model,                                                                       &
               land_surface,                                                                       &
               large_scale_forcing,                                                                &
               loop_optimization,                                                                  &
               lsf_surf,                                                                           &
               neutral,                                                                            &
               passive_scalar,                                                                     &
               pt_surface,                                                                         &
               q_surface,                                                                          &
               rho_cp,                                                                             &
               slurb,                                                                              &
               surface_pressure,                                                                   &
               simulated_time,                                                                     &
               time_since_reference_point,                                                         &
               urban_surface,                                                                      &
               use_free_convection_scaling

#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY:  enable_openacc
#endif

    USE kinds

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               microphysics_morrison,                                                              &
               microphysics_seifert

    USE pegrid

    USE land_surface_model_mod,                                                                    &
        ONLY:  aero_resist_kray,                                                                   &
               lsm_start_time

    USE surface_mod,                                                                               &
        ONLY :  surf_type,                                                                         &
                surf_def,                                                                          &
                surf_lsm,                                                                          &
                surf_u,                                                                            &
                surf_usm,                                                                          &
                surf_v,                                                                            &
                surf_w


    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< loop index x direction
    INTEGER(iwp) ::  i_off  !< offset index between surface and reference grid point in x direction
    INTEGER(iwp) ::  j      !< loop index y direction
    INTEGER(iwp) ::  j_off  !< offset index between surface and reference grid point in y direction
    INTEGER(iwp) ::  k      !< loop index z direction
    INTEGER(iwp) ::  k_off  !< offset index between surface and reference grid point in z direction
    INTEGER(iwp) ::  m      !< running index surface elements

    REAL(wp), PARAMETER ::  ol_max   = 1.0E6_wp   !< allowed absolute maximum value Obukhov length
    REAL(wp), PARAMETER ::  ol_min   = 1.0E-6_wp  !< allowed absolute minimum value Obukhov length
    REAL(wp), PARAMETER ::  ol_tol   = 1.0E-4_wp  !< convergence limit for Obukhov length, relative tolerance
    REAL(wp), PARAMETER ::  rib_max  = 1.0E1_wp   !< maximum bulk Richardson number (absolute value)
    REAL(wp), PARAMETER ::  zeta_min = 1.0E-4_wp  !< minimum stability parameter absolute value (neutral limit)

    REAL(wp) ::  e_s     !< saturation water vapor pressure
    REAL(wp) ::  z_mo    !< height of the constant flux layer where MOST is assumed



    TYPE(surf_type), POINTER ::  surf  !< surf-type array, used to generalize subroutines


    SAVE

    PRIVATE

    PUBLIC calc_ol,                                                                                &
           calc_rib,                                                                               &
           init_surface_layer_fluxes,                                                              &
           phi_m,                                                                                  &
           psi_h,                                                                                  &
           psi_m,                                                                                  &
           surface_layer_fluxes

    INTERFACE calc_ol
       MODULE PROCEDURE calc_ol
    END INTERFACE calc_ol

    INTERFACE calc_rib
       MODULE PROCEDURE calc_rib
    END INTERFACE calc_rib

    INTERFACE init_surface_layer_fluxes
       MODULE PROCEDURE init_surface_layer_fluxes
    END INTERFACE init_surface_layer_fluxes

    INTERFACE phi_m
       MODULE PROCEDURE phi_m
    END INTERFACE phi_m

    INTERFACE psi_h
       MODULE PROCEDURE psi_h_dp
       MODULE PROCEDURE psi_h_sp
    END INTERFACE psi_h

    INTERFACE psi_m
       MODULE PROCEDURE psi_m_dp
       MODULE PROCEDURE psi_m_sp
    END INTERFACE psi_m

    INTERFACE surface_layer_fluxes
       MODULE PROCEDURE surface_layer_fluxes
    END INTERFACE surface_layer_fluxes


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Main routine to compute the surface fluxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_layer_fluxes

    IMPLICIT NONE

    IF ( debug_output_timestep )  CALL debug_message( 'surface_layer_fluxes', 'start' )
!
!-- First, calculate the new Obukhov length from precalculated values of log(z/z0) and wind speeds.
!-- As a second step, then calculate new friction velocity, followed by the new scaling
!-- parameters (th*, q*, etc.), and the new surface fluxes, if required. Note, each routine is called
!-- for different surface types. First call for default-type surfaces, then for natural- and
!-- urban-type surfaces.
!-- Start with default-type surfaces
    IF ( surf_def%ns >= 1 )  THEN
       surf => surf_def
!
!--    First, precalculate ln(z/z0) for all surfaces. This is done each timestep, in order
!--    to account for time-dependent roughness or user-modifications.
       CALL calc_ln
!
!--    Calculate temperatures and specific humidity at the first computational grid level
!--    above the surface and store values in the 1d surface data arrays.
       CALL calc_pt_q
!
!--    Store surface temperatures and specific humidity of the Eularian fields in the
!--    1d surface data arrays.
       IF ( .NOT. neutral )  THEN
          CALL store_pt_surface
          IF ( humidity )  THEN
             CALL store_q_surface
             CALL store_vpt_surface
          ENDIF
       ENDIF
!
!--    Calculate surface-parallel absolute velocity on different positions on the staggered grid.
       CALL calc_uvw_abs_s
       CALL calc_uvw_abs_uv
       CALL calc_uvw_abs_w
!
!--    Calculate Richardson flux number and Obukhov length.
       IF ( .NOT. neutral )  THEN
          IF ( ibc_pt_b == 1 )  THEN
             IF ( humidity )  THEN
                CALL calc_rib_with_prescribed_fluxes( surf%ns, surf%downward, surf%k, surf%koff,   &
                                                      surf%pt1, surf%rib, surf%shf, surf%upward,   &
                                                      surf%uvw_abs, surf%z_mo, surf%qsws,          &
                                                      surf%qv1, surf%vpt1 )
             ELSE
                CALL calc_rib_with_prescribed_fluxes( surf%ns, surf%downward, surf%k, surf%koff,   &
                                                      surf%pt1, surf%rib, surf%shf, surf%upward,   &
                                                      surf%uvw_abs, surf%z_mo )
             ENDIF
          ELSE
             IF ( humidity )  THEN
                CALL calc_rib( surf%ns, surf%vpt1, surf%vpt_surface, surf%rib, surf%uvw_abs,       &
                               surf%z_mo )
             ELSE
                CALL calc_rib( surf%ns, surf%pt1, surf%pt_surface, surf%rib, surf%uvw_abs,         &
                               surf%z_mo )
             ENDIF
          ENDIF
          CALL calc_ol( surf%ns, surf%ln_z_z0, surf%ln_z_z0h, surf%ol, surf%rib, surf%z0,          &
                        surf%z0h, surf%z_mo )
       ENDIF
!
!--    Calculate friction velocity representative for different positions on the staggered grid.
       CALL calc_us_s
       CALL calc_us_uv
       CALL calc_us_w
!
!--    Calculate scaling parameters.
       CALL calc_scaling_parameters
!
!--    Calculate surface fluxes for scalars.
       CALL calc_surface_fluxes
!
!--    Calculate surface momentum fluxes. Note, not all fluxes become effective at all surface
!--    orientations. For example, u'w'_0 is zero at vertical walls. Skip this action in case
!--    of cut-cell surfaces. This case, fluxes will be computed separately on dedicated surface
!--    types as the number of surfaces relevant for the staggered grids might be different
!--    compared to the grid-cell center.
       CALL calc_usws( surf%us )
       CALL calc_vsws( surf%us )
       IF ( .NOT. cut_cell_topography )  THEN
          CALL calc_usvs
          CALL calc_vsus
          CALL calc_wsus_wsvs
       ENDIF
!
!--    Calculate surface momentum fluxes on scalar grid. This is required for TKE production.
       CALL calc_usws_vsws_for_tke
!
!--    Calculate aerodynamic resistance.
       IF ( det_enabled )  CALL calc_aerodynamic_resistance
    ENDIF
!
!-- Natural land surfaces.
    IF ( surf_lsm%ns >= 1 )  THEN
       surf => surf_lsm
!
!--    First, precalculate ln(z/z0) for all surfaces. This is done each timestep, in order
!--    to account for time-dependent roughness or user-modifications.
       CALL calc_ln
!
!--    Derive potential temperature and specific humidity at first grid level from the fields
!--    pt and q
       CALL calc_pt_q
!
!--    Calculate surface-parallel absolute velocity on different positions on the staggered grid.
       CALL calc_uvw_abs_s
       CALL calc_uvw_abs_uv
       CALL calc_uvw_abs_w
!
!--    Calculate Richardson flux number and Obukhov length.
       IF ( .NOT. neutral )  THEN
          IF ( ibc_pt_b == 1 )  THEN
             IF ( humidity )  THEN
                CALL calc_rib_with_prescribed_fluxes( surf%ns, surf%downward, surf%k, surf%koff,   &
                                                      surf%pt1, surf%rib, surf%shf, surf%upward,   &
                                                      surf%uvw_abs, surf%z_mo, surf%qsws,          &
                                                      surf%qv1, surf%vpt1 )
             ELSE
                CALL calc_rib_with_prescribed_fluxes( surf%ns, surf%downward, surf%k, surf%koff,   &
                                                      surf%pt1, surf%rib, surf%shf, surf%upward,   &
                                                      surf%uvw_abs, surf%z_mo )
             ENDIF
          ELSE
             IF ( humidity )  THEN
                CALL calc_rib( surf%ns, surf%vpt1, surf%vpt_surface, surf%rib, surf%uvw_abs,       &
                               surf%z_mo )
             ELSE
                CALL calc_rib( surf%ns, surf%pt1, surf%pt_surface, surf%rib, surf%uvw_abs,         &
                               surf%z_mo )
             ENDIF
          ENDIF
          CALL calc_ol( surf%ns, surf%ln_z_z0, surf%ln_z_z0h, surf%ol, surf%rib, surf%z0,          &
                        surf%z0h, surf%z_mo )
       ENDIF
!
!--    Calculate friction velocity representative for different positions on the staggered grid.
       CALL calc_us_s
       CALL calc_us_uv
       CALL calc_us_w
!
!--    Calculate scaling parameters.
       CALL calc_scaling_parameters
!
!--    Calculate surface fluxes for scalars.
       CALL calc_surface_fluxes
!
!--    Calculate surface momentum fluxes. Note, not all fluxes become effective at all surface
!--    orientations. For example, u'w'_0 is zero at vertical walls. Skip this action in case
!--    of cut-cell surfaces. This case, fluxes will be computed separately on dedicated surface
!--    types as the number of surfaces relevant for the staggered grids might be different
!--    compared to the grid-cell center.
       CALL calc_usws( surf%us )
       CALL calc_vsws( surf%us )
       IF ( .NOT. cut_cell_topography )  THEN
          CALL calc_usvs
          CALL calc_vsus
          CALL calc_wsus_wsvs
       ENDIF
!
!--    Calculate surface momentum fluxes on scalar grid. This is required for TKE production.
       CALL calc_usws_vsws_for_tke
    ENDIF
!
!-- Building surfaces.
    IF ( surf_usm%ns >= 1 )  THEN
       surf => surf_usm
!
!--    First, precalculate ln(z/z0) for all surfaces. This is done each timestep, in order
!--    to account for time-dependent roughness or user-modifications.
       CALL calc_ln
!
!--    Derive potential temperature and specific humidity at first grid level from the fields
!--    pt and q
       CALL calc_pt_q
!
!--    Calculate surface-parallel absolute velocity on different positions on the staggered grid.
       CALL calc_uvw_abs_s
       CALL calc_uvw_abs_uv
       CALL calc_uvw_abs_w
!
!--    Calculate Richardson flux number and Obukhov length.
       IF ( .NOT. neutral )  THEN
          IF ( ibc_pt_b == 1 )  THEN
             IF ( humidity )  THEN
                CALL calc_rib_with_prescribed_fluxes( surf%ns, surf%downward, surf%k, surf%koff,   &
                                                      surf%pt1, surf%rib, surf%shf, surf%upward,   &
                                                      surf%uvw_abs, surf%z_mo, surf%qsws,          &
                                                      surf%qv1, surf%vpt1 )
             ELSE
                CALL calc_rib_with_prescribed_fluxes( surf%ns, surf%downward, surf%k, surf%koff,   &
                                                      surf%pt1, surf%rib, surf%shf, surf%upward,   &
                                                      surf%uvw_abs, surf%z_mo )
             ENDIF
          ELSE
             IF ( humidity )  THEN
                CALL calc_rib( surf%ns, surf%vpt1, surf%vpt_surface, surf%rib, surf%uvw_abs,       &
                               surf%z_mo )
             ELSE
                CALL calc_rib( surf%ns, surf%pt1, surf%pt_surface, surf%rib, surf%uvw_abs,         &
                               surf%z_mo )
             ENDIF
          ENDIF
          CALL calc_ol( surf%ns, surf%ln_z_z0, surf%ln_z_z0h, surf%ol, surf%rib, surf%z0,          &
                        surf%z0h, surf%z_mo )
       ENDIF
!
!--    Calculate friction velocity representative for different positions on the staggered grid.
       CALL calc_us_s
       CALL calc_us_uv
       CALL calc_us_w
!
!--    Calculate scaling parameters.
       CALL calc_scaling_parameters
!
!--    Calculate surface fluxes for scalars.
       CALL calc_surface_fluxes
!
!--    Calculate surface momentum fluxes. Note, not all fluxes become effective at all surface
!--    orientations. For example, u'w'_0 is zero at vertical walls. Skip this action in case
!--    of cut-cell surfaces. This case, fluxes will be computed separately on dedicated surface
!--    types as the number of surfaces relevant for the staggered grids might be different
!--    compared to the grid-cell center (see further below).
       CALL calc_usws( surf%us )
       CALL calc_vsws( surf%us )
       IF ( .NOT. cut_cell_topography )  THEN
          CALL calc_usvs
          CALL calc_vsus
          CALL calc_wsus_wsvs
       ENDIF
!
!--    Calculate surface momentum fluxes on scalar grid. This is required for TKE production.
       CALL calc_usws_vsws_for_tke
    ENDIF
!
!-- Compute momentum fluxes on staggered grid in case of cut-cell topography. This case, momentum
!-- fluxes are considered on separate surface arrays for the u-, v- and w-component rather than
!-- the usual surface types which correspond to the grid center. This is necessary as the
!-- number of surfaces relevant for the staggered grids might be different compared to the
!-- grid-cell center in case of the cut-cell topography.
    IF ( cut_cell_topography )  THEN
       IF ( surf_u%ns >= 1 )  THEN
          surf => surf_u

          CALL calc_ln
          CALL calc_uvw_abs_uv
!
!--       Transfer data used for stability correction from the surf_def, surf_lsm and surf_usm
!--       types onto surf_u.
          CALL transfer_ol

          CALL calc_us_uv

          CALL calc_usws( surf%us_uvgrid )
          CALL calc_usvs
       ENDIF

       IF ( surf_v%ns >= 1 )  THEN
          surf => surf_v

          CALL calc_ln
          CALL calc_uvw_abs_uv
!
!--       Transfer data used for stability correction from the surf_def, surf_lsm and surf_usm
!--       types onto surf_v.
          CALL transfer_ol

          CALL calc_us_uv
          CALL calc_vsws( surf%us_uvgrid )
          CALL calc_vsus
       ENDIF

       IF ( surf_w%ns >= 1 )  THEN
          surf => surf_w

          CALL calc_ln
          CALL calc_uvw_abs_w
!
!--       Transfer data used for stability correction from the surf_def, surf_lsm and surf_usm
!--       types onto surf_w.
          CALL transfer_ol

          CALL calc_us_w
          CALL calc_wsus_wsvs
       ENDIF
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'surface_layer_fluxes', 'end' )

 END SUBROUTINE surface_layer_fluxes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializing actions for the surface layer routine.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_surface_layer_fluxes

    IMPLICIT NONE

    CALL location_message( 'initializing surface layer', 'start' )
!
!-- In case of runs with neutral statification, set Obukhov length to a large value
    IF ( neutral )  THEN
       IF ( surf_def%ns >= 1  .AND.  ALLOCATED( surf_def%ol ) )  surf_def%ol = 1.0E10_wp
       IF ( surf_lsm%ns >= 1  .AND.  ALLOCATED( surf_lsm%ol ) )  surf_lsm%ol = 1.0E10_wp
       IF ( surf_usm%ns >= 1  .AND.  ALLOCATED( surf_usm%ol ) )  surf_usm%ol = 1.0E10_wp
    ENDIF

    CALL location_message( 'initializing surface layer', 'finished' )

 END SUBROUTINE init_surface_layer_fluxes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute ln(z/z0).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_ln

!
!-- Note, ln(z/z0h) and ln(z/z0q) is also calculated in neutral situations.
!-- This is because the roughness for scalars are also used for other quantities such as passive
!-- scalar, chemistry and aerosols.
    !$OMP PARALLEL DO PRIVATE( z_mo )
    !$ACC PARALLEL LOOP PRIVATE(z_mo) &
    !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       z_mo = surf%z_mo(m)
       surf%ln_z_z0(m)  = LOG( z_mo / surf%z0(m) )
    ENDDO

    IF ( ALLOCATED( surf%ln_z_z0h ) )  THEN
       !$OMP PARALLEL DO PRIVATE( z_mo )
       !$ACC PARALLEL LOOP PRIVATE(z_mo) &
       !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, surf%ns
          z_mo = surf%z_mo(m)
          surf%ln_z_z0h(m) = LOG( z_mo / surf%z0h(m) )
       ENDDO
    ENDIF

    IF ( ALLOCATED( surf%ln_z_z0q ) )  THEN
       !$OMP PARALLEL DO PRIVATE( z_mo )
       !$ACC PARALLEL LOOP PRIVATE(z_mo) &
       !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, surf%ns
          z_mo = surf%z_mo(m)
          surf%ln_z_z0q(m) = LOG( z_mo / surf%z0q(m) )
       ENDDO
    ENDIF

 END SUBROUTINE calc_ln


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the surface-parallel velocity (relative to the surface)
!> representative for the grid center.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_uvw_abs_s

    IMPLICIT NONE

    INTEGER(iwp) ::  ibit  !< flag to mask computation of relative velocity in case of downward-facing surfaces

    REAL(wp) ::  u_comp   !< u-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  v_comp   !< v-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  vel_dot  !< dot product of velocity and surface normal vector
    REAL(wp) ::  w_comp   !< w-component of the velocity vector interpolated onto respective staggered grid position

    REAL(wp), DIMENSION(1:surf%ns) ::  w_lfc   !< local free convection velocity scale
    !$ACC DECLARE CREATE( w_lfc )
!
!-- Pre-calculate local free-convection scaling parameter if required. This
!-- will maintain a horizontal velocity even for very weak wind convective conditions. SIGN is
!-- used to set w_lfc to zero under stable conditions. Note, free-convection scaling velocity
!-- is only used at horizontal surfaces (multiplication with normal-vector component).
    IF ( use_free_convection_scaling )  THEN
       !$OMP PARALLEL DO
       !$ACC PARALLEL LOOP &
       !$ACC PRESENT(surf) IF(enable_openacc)
       DO  m = 1, surf%ns
          w_lfc(m) = ABS( g / surf%pt1(m) * surf%z_mo(m) * surf%shf(m) * surf%n_s(m,1) )
          w_lfc(m) = ( 0.5_wp * ( w_lfc(m)                                                         &
                                + SIGN( w_lfc(m) , surf%shf(m) * surf%n_s(m,1) ) ) )**(0.33333_wp) &
                     * MERGE( 1.0_wp, 0.0_wp, surf%consider_stability(m) )
       ENDDO
    ELSE
       !$ACC PARALLEL LOOP PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, surf%ns
          w_lfc(m) = 0.0_wp
       ENDDO
    ENDIF

    !$OMP PARALLEL DO PRIVATE(i, ibit, j, k, u_comp, v_comp, w_comp, vel_dot)
    !$ACC PARALLEL LOOP PRIVATE(i, ibit, j, k, u_comp, v_comp, w_comp, vel_dot) &
    !$ACC PRESENT(surf, u, v, w) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)
!
!--    ibit is 1 for upward-facing surfaces, zero for downward-facing surfaces.
       ibit = MERGE( 1, 0, surf%upward(m) )
!
!--    Compute surface-parallel velocity representative for the grid center. Compute dot product
!--    between the grid-center inpolated velocity and surface normal vector components. Note,
!--    horizontal velocity components are considered as relative velocities, which takes coupled
!--    atmosphere ocean surfaces into account (see the k-1 values). Relative velocities, however,
!--    are only considered in case of horizontal upward-facing surfaces (see ibit).
       u_comp = 0.5 * ( u(k,j,i) + u(k,j,i+1) - ( u(k-1,j,i) + u(k-1,j,i+1) ) * ibit )
       v_comp = 0.5 * ( v(k,j,i) + v(k,j+1,i) - ( v(k-1,j,i) + v(k-1,j+1,i) ) * ibit )
       w_comp = 0.5 * ( w(k,j,i) + w(k-1,j,i) )

       vel_dot = u_comp * surf%n_s(m,3) + v_comp * surf%n_s(m,2) + w_comp * surf%n_s(m,1)
       surf%uvw_abs(m) = SQRT( ( u_comp - vel_dot * surf%n_s(m,3) )**2                             &
                             + ( v_comp - vel_dot * surf%n_s(m,2) )**2                             &
                             + ( w_comp - vel_dot * surf%n_s(m,1) )**2 + w_lfc(m)**2 )
    ENDDO

 END SUBROUTINE calc_uvw_abs_s


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the surface-parallel velocity (relative to the surface) for the
!> u- and v-grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_uvw_abs_uv

    IMPLICIT NONE

    LOGICAL  ::  u_grid   !< flag indicating that absolute velocity is required on u-grid needed for north/soutward-facing surfaces

    REAL(wp) ::  u_comp   !< u-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  v_comp   !< v-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  vel_dot  !< dot product of velocity and surface normal vector
    REAL(wp) ::  w_comp   !< w-component of the velocity vector interpolated onto respective staggered grid position

    !$OMP PARALLEL DO PRIVATE(i, j, k, u_comp, u_grid, v_comp, w_comp, vel_dot)
    !$ACC PARALLEL LOOP PRIVATE(i, j, k, u_comp, u_grid, v_comp, w_comp, vel_dot) &
    !$ACC PRESENT(surf, u, v, w) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)
!
!--    Now compute the absolute velocity on the u- or v-grid, depending on the respective
!--    surface grid point orientation.
       u_grid = surf%joff(m) /= 0
!
!--    Compute the surface-parallel absolute velocity on u- and v-grid. Depending on the
!--    considered grid (u- or v-grid), velocities are interpolated to that grid location
!--    beforehand. At north/south facing walls, the respective wall-parallel velocity actually does
!--    not include any contribution from the v-component. Even in case of straight walls this
!--    maintained by the following code (which yields to identical results as before), even though
!--    the respective v-component on the u-grid is not necessarily zero here. However, later on,
!--    this portion cancels out by the normal vector component.
       u_comp = MERGE( u(k,j,i),                                                                   & ! u-grid
                       0.25_wp * ( u(k,j,i) + u(k,j,i+1) + u(k,j-1,i) + u(k,j-1,i+1) ),            & ! v-grid
                       u_grid )
       v_comp = MERGE( 0.25_wp * ( v(k,j,i-1) + v(k,j,i) + v(k,j+1,i-1) + v(k,j+1,i) ),            & ! u-grid
                       v(k,j,i),                                                                   & ! v-grid
                       u_grid )
       w_comp = MERGE( 0.25_wp * ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) + w(k,j,i) ),            & ! u-grid
                       0.25_wp * ( w(k-1,j-1,i) + w(k-1,j,i) + w(k,j-1,i) + w(k,j,i) ),            & ! v-grid
                       u_grid )

       vel_dot = u_comp * surf%n_s(m,3) + v_comp * surf%n_s(m,2) + w_comp * surf%n_s(m,1)
       surf%uvw_abs_uv(m) = SQRT( ( u_comp - vel_dot * surf%n_s(m,3) )**2                          &
                                + ( v_comp - vel_dot * surf%n_s(m,2) )**2                          &
                                + ( w_comp - vel_dot * surf%n_s(m,1) )**2 )
    ENDDO

 END SUBROUTINE calc_uvw_abs_uv


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the surface-parallel velocity (relative to the surface) for the
!> w-grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_uvw_abs_w

    IMPLICIT NONE

    REAL(wp) ::  u_comp   !< u-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  v_comp   !< v-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  vel_dot  !< dot product of velocity and surface normal vector
    REAL(wp) ::  w_comp   !< w-component of the velocity vector interpolated onto respective staggered grid position

    !$OMP PARALLEL DO PRIVATE(i, j, k, u_comp, v_comp, w_comp, vel_dot)
    !$ACC PARALLEL LOOP PRIVATE(i, j, k, u_comp, v_comp, w_comp, vel_dot) &
    !$ACC PRESENT(surf, u, v, w) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)
!
!--    Now compute the surface-parallel absolute velocity on w-grid.
       u_comp = 0.25_wp * ( u(k+1,j,i+1) + u(k+1,j,i) + u(k,j,i+1) + u(k,j,i) )
       v_comp = 0.25_wp * ( v(k+1,j+1,i) + v(k+1,j,i) + v(k,j+1,i) + v(k,j,i) )
       w_comp = w(k,j,i)

       vel_dot = u_comp * surf%n_s(m,3) + v_comp * surf%n_s(m,2) + w_comp * surf%n_s(m,1)
       surf%uvw_abs_w(m) = SQRT( ( u_comp - vel_dot * surf%n_s(m,3) )**2                           &
                               + ( v_comp - vel_dot * surf%n_s(m,2) )**2                           &
                               + ( w_comp - vel_dot * surf%n_s(m,1) )**2 )
    ENDDO

 END SUBROUTINE calc_uvw_abs_w


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the Obukhov length (L). On GPUs, this routine always needs to run with 64-bit precision
!> (in the vector branch), because with 32-bit the algorithm partly does not converge within the
!> maximum of 1000 iterations.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_ol( ns, ln_z_z0, ln_z_z0h, ol, rib, z0, z0h, z_mo )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  ns  !< number of surface points

    REAL(wp), INTENT(IN),    DIMENSION(ns) ::  ln_z_z0   !< logarithm (z/z0)
    REAL(wp), INTENT(IN),    DIMENSION(ns) ::  ln_z_z0h  !< logarithm (z/z0h)
    REAL(wp), INTENT(INOUT), DIMENSION(ns) ::  ol        !< Obukhov length
    REAL(wp), INTENT(IN),    DIMENSION(ns) ::  rib       !< Richardson flux number
    REAL(wp), INTENT(IN),    DIMENSION(ns) ::  z0        !< rougness length for momentum
    REAL(wp), INTENT(IN),    DIMENSION(ns) ::  z0h       !< rougness length for scalar quantities
    REAL(wp), INTENT(IN),    DIMENSION(ns) ::  z_mo      !< constant flux layer height

    INTEGER(iwp) ::  iter  !< Newton iteration step

    LOGICAL, DIMENSION(ns) ::  convergence_reached  !< convergence switch for vectorization

    REAL(wp) ::  f        !< function for Newton iteration: f = Ri - [...]/[...]^2 = 0
    REAL(wp) ::  f_d_ol   !< derivative of f
    REAL(wp) ::  ol_l     !< lower bound of L for Newton iteration
    REAL(wp) ::  ol_m     !< previous value of L for Newton iteration
    REAL(wp) ::  ol_prev  !< previous time step value of L
    REAL(wp) ::  ol_u     !< upper bound of L for Newton iteration
!
!-- Variables for vector branch.
#if defined( _OPENACC )
!-- GPUs require 64 bit precision, because otherwise the algorithm often does not converge.
    REAL(dp) ::  arg_h_zmo_m  !< argument for psi-function
    REAL(dp) ::  arg_h_z0h_m  !< argument for psi-function
    REAL(dp) ::  arg_m_zmo_m  !< argument for psi-function
    REAL(dp) ::  arg_m_z0_m   !< argument for psi-function
    REAL(dp) ::  arg_h_zmo_u  !< argument for psi-function
    REAL(dp) ::  arg_h_z0h_u  !< argument for psi-function
    REAL(dp) ::  arg_m_zmo_u  !< argument for psi-function
    REAL(dp) ::  arg_m_z0_u   !< argument for psi-function
    REAL(dp) ::  arg_h_zmo_l  !< argument for psi-function
    REAL(dp) ::  arg_h_z0h_l  !< argument for psi-function
    REAL(dp) ::  arg_m_zmo_l  !< argument for psi-function
    REAL(dp) ::  arg_m_z0_l   !< argument for psi-function
    REAL(dp) ::  f_vec        !< function for Newton iteration: f = Ri - [...]/[...]^2 = 0
    REAL(dp) ::  f_d_ol_vec   !< derivative of f
    REAL(dp) ::  ol_l_vec     !< lower bound of L for Newton iteration
    REAL(dp) ::  ol_m_vec     !< previous value of L for Newton iteration
    REAL(dp) ::  ol_u_vec     !< upper bound of L for Newton iteration

    REAL(dp), DIMENSION(ns) ::  ol_prev_vec  !< temporary array required for vectorization

    !$ACC DATA CREATE( convergence_reached, ol_prev_vec ) IF(enable_openacc)
#else
    REAL(wp) ::  arg_h_zmo_m  !< argument for psi-function
    REAL(wp) ::  arg_h_z0h_m  !< argument for psi-function
    REAL(wp) ::  arg_m_zmo_m  !< argument for psi-function
    REAL(wp) ::  arg_m_z0_m   !< argument for psi-function
    REAL(wp) ::  arg_h_zmo_u  !< argument for psi-function
    REAL(wp) ::  arg_h_z0h_u  !< argument for psi-function
    REAL(wp) ::  arg_m_zmo_u  !< argument for psi-function
    REAL(wp) ::  arg_m_z0_u   !< argument for psi-function
    REAL(wp) ::  arg_h_zmo_l  !< argument for psi-function
    REAL(wp) ::  arg_h_z0h_l  !< argument for psi-function
    REAL(wp) ::  arg_m_zmo_l  !< argument for psi-function
    REAL(wp) ::  arg_m_z0_l   !< argument for psi-function
    REAL(wp) ::  f_vec        !< function for Newton iteration: f = Ri - [...]/[...]^2 = 0
    REAL(wp) ::  f_d_ol_vec   !< derivative of f
    REAL(wp) ::  ol_l_vec     !< lower bound of L for Newton iteration
    REAL(wp) ::  ol_m_vec     !< previous value of L for Newton iteration
    REAL(wp) ::  ol_u_vec     !< upper bound of L for Newton iteration

    REAL(wp), DIMENSION(ns) ::  ol_prev_vec  !< temporary array required for vectorization
#endif


!
!-- Calculate the Obukhov length using Newton iteration.
    IF ( loop_optimization == 'cache' )  THEN

       !$OMP PARALLEL DO PRIVATE( ol_prev, iter, ol_m, ol_l, ol_u, f, f_d_ol)
       !$ACC PARALLEL LOOP &
       !$ACC PRIVATE(ol_prev, ol_m, ol_l, ol_u, f, f_d_ol) &
       !$ACC PRESENT(ns, ln_z_z0, ln_z_z0h, ol, rib, z0, z0h, z_mo) IF(enable_openacc)
       DO  m = 1, ns
!
!--       Store current value in case the Newton iteration fails.
          ol_prev = ol(m)
!
!--       Flip the sign of the initial Obukhov length if the stability has changed from stable to
!--       unstable or vice versa and set it to a moderate value. A moderate value is also chosen,
!--       if the Obukhov length from the last time step reached the maximum threshold value.
          IF ( rib(m) * ol(m) < 0.0_wp  .OR.  ABS( ol(m) ) == ol_max )  THEN
             IF ( rib(m) > 0.0_wp )  ol(m) =  100.0_wp
             IF ( rib(m) < 0.0_wp )  ol(m) = -100.0_wp
          ENDIF
!
!--       Iteration to find Obukhov length.
          iter = 0
          DO
             iter = iter + 1
!
!--          In case of divergence, use the value of the previous time step.
             IF ( iter > 1000 )  THEN
                ol(m) = ol_prev
                EXIT
             ENDIF

!
!--          Calculate step size for central difference.
             ol_m = ol(m)
             ol_l = ol_m - 0.001_wp * ol_m
             ol_u = ol_m + 0.001_wp * ol_m

             IF ( ibc_pt_b /= 1 )  THEN
!
!--             Calculate f = Ri - [...]/[...]^2 = 0.
                f = rib(m) - ( z_mo(m) / ol_m ) * ( ln_z_z0h(m) - psi_h( z_mo(m) / ol_m )          &
                                                                + psi_h( z0h(m)  / ol_m ) )        &
                                                / ( ln_z_z0(m)  - psi_m( z_mo(m) / ol_m )          &
                                                                + psi_m( z0(m)   / ol_m ) )**2
!
!--             Calculate df/dL.
                f_d_ol = ( - ( z_mo(m) / ol_u ) * ( ln_z_z0h(m) - psi_h( z_mo(m) / ol_u )          &
                                                                + psi_h( z0h(m)  / ol_u ) )        &
                                                / ( ln_z_z0(m)  - psi_m( z_mo(m) / ol_u )          &
                                                                + psi_m( z0(m)   / ol_u ) )**2     &
                           + ( z_mo(m) / ol_l ) * ( ln_z_z0h(m) - psi_h( z_mo(m) / ol_l )          &
                                                                + psi_h( z0h(m)  / ol_l ) )        &
                                                / ( ln_z_z0(m)  - psi_m( z_mo(m) / ol_l )          &
                                                                + psi_m( z0(m)   / ol_l ) )**2     &
                         ) / ( ol_u - ol_l )
             ELSE
!
!--             Calculate f = Ri - 1 /[...]^3 = 0.
                f = rib(m) - ( z_mo(m) / ol_m ) /                                                  &
                             ( ln_z_z0(m) - psi_m( z_mo(m) / ol_m ) + psi_m( z0(m) / ol_m ) )**3

!
!--             Calculate df/dL.
                f_d_ol = ( - ( z_mo(m) / ol_u ) / ( ln_z_z0(m) - psi_m( z_mo(m) / ol_u )           &
                                                               + psi_m( z0(m)   / ol_u ) )**3      &
                           + ( z_mo(m) / ol_l ) / ( ln_z_z0(m) - psi_m( z_mo(m) / ol_l )           &
                                                               + psi_m( z0(m)   / ol_l ) )**3      &
                         ) / ( ol_u - ol_l )
             ENDIF
!
!--          Calculate new L.
             ol(m) = ol_m - f / f_d_ol
!
!--          Ensure that the bulk Richardson number and the Obukhov length have the same sign and
!--          ensure convergence. If the sign is not the same, the above calculated Obukhov length
!--          obviously overshooted to the opposite side, so the next iteration should start with
!--          a smaller value.
             IF ( ol(m) * ol_m < 0.0_wp )  ol(m) = ol_m * 0.5_wp
!
!--          In the deep neutral zone, set L to the maximum allowed value.
             IF ( ABS( ol(m) ) > ol_max )  THEN
                ol(m) = SIGN( ol_max, ol(m) )
                EXIT
             ENDIF
!
!--          Assure that Obukhov length does not become zero.
             IF ( ABS( ol(m) ) < ol_min )  THEN
                ol(m) = SIGN( ol_min, ol(m) )
                EXIT
             ENDIF
!
!--          Check for convergence.
             IF ( ABS( ( ol(m) - ol_m ) /  ol(m) ) < ol_tol )  EXIT

          ENDDO
       ENDDO

    ELSE
!
!--    Calculate the Obukhov length using Newton iteration (vectorized version).
!--    First set arrays required for vectorization.
       !$ACC PARALLEL LOOP &
       !$ACC PRESENT(ol, rib) DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, ns
!
!--       Store current value in case the Newton iteration fails.
          ol_prev_vec(m) = ol(m)
!
!--       Flip the sign of the initial Obukhov length if the stability has changed from stable to
!--       unstable or vice versa and set it to a moderate value. A moderate value is also chosen,
!--       if the Obukhov length from the last time step reached the maximum threshold value.
          IF ( rib(m) * ol(m) < 0.0_wp  .OR.  ABS( ol(m) ) == ol_max )  THEN
             IF ( rib(m) > 0.0_wp )  ol(m) =  100.0_wp
             IF ( rib(m) < 0.0_wp )  ol(m) = -100.0_wp
          ENDIF
!
!--       Initialize convergence flag.
          convergence_reached(m) = .FALSE.
       ENDDO

!
!--    Iteration to find Obukhov length
       iter = 0
       DO
          iter = iter + 1
!
!--       In case of divergence, use the value(s) of the previous time step.
          IF ( iter > 1000 )  THEN
             !$ACC PARALLEL LOOP &
             !$ACC PRESENT(ns, ol) IF(enable_openacc)
             DO  m = 1, ns
                IF ( .NOT. convergence_reached(m) )  ol(m) = ol_prev_vec(m)
             ENDDO
             EXIT
          ENDIF

          !$ACC PARALLEL LOOP PRIVATE(ol_m_vec, ol_l_vec, ol_u_vec, f_vec, f_d_ol_vec) &
          !$ACC PRIVATE(arg_h_zmo_m, arg_h_z0h_m, arg_m_zmo_m, arg_m_z0_m) &
          !$ACC PRIVATE(arg_h_zmo_u, arg_h_z0h_u, arg_m_zmo_u, arg_m_z0_u) &
          !$ACC PRIVATE(arg_h_zmo_l, arg_h_z0h_l, arg_m_zmo_l, arg_m_z0_l) &
          !$ACC PRESENT(ns, ln_z_z0, ln_z_z0h, ol, rib, z0, z0h, z_mo) DEFAULT(NONE) IF(enable_openacc)
          DO  m = 1, ns
             IF ( convergence_reached(m) )  CYCLE

!
!--          Calculate step size for central difference.
             ol_m_vec = ol(m)
             ol_l_vec = ol_m_vec - 0.001_wp * ol_m_vec
             ol_u_vec = ol_m_vec + 0.001_wp * ol_m_vec

             IF ( ibc_pt_b /= 1 )  THEN
!
!--             Calculate f = Ri - [...]/[...]^2 = 0.
                arg_h_zmo_m = z_mo(m) / ol_m_vec
                arg_h_z0h_m = z0h(m)  / ol_m_vec
                arg_m_zmo_m = z_mo(m) / ol_m_vec
                arg_m_z0_m  = z0(m)   / ol_m_vec
                f_vec = rib(m) - ( arg_h_zmo_m ) * ( ln_z_z0h(m) - psi_h( arg_h_zmo_m )            &
                                                                 + psi_h( arg_h_z0h_m ) )          &
                                                 / ( ln_z_z0(m)  - psi_m( arg_m_zmo_m )            &
                                                                 + psi_m( arg_m_z0_m ) )**2
!
!--             Calculate df/dL.
                arg_h_zmo_u = z_mo(m) / ol_u_vec
                arg_h_z0h_u = z0h(m)  / ol_u_vec
                arg_m_zmo_u = z_mo(m) / ol_u_vec
                arg_m_z0_u  = z0(m)   / ol_u_vec
                arg_h_zmo_l = z_mo(m) / ol_l_vec
                arg_h_z0h_l = z0h(m)  / ol_l_vec
                arg_m_zmo_l = z_mo(m) / ol_l_vec
                arg_m_z0_l  = z0(m)   / ol_l_vec
                f_d_ol_vec = ( - ( arg_h_zmo_u ) * ( ln_z_z0h(m) - psi_h( arg_h_zmo_u )            &
                                                                 + psi_h( arg_h_z0h_u ) )          &
                                                 / ( ln_z_z0(m)  - psi_m( arg_m_zmo_u )            &
                                                                 + psi_m( arg_m_z0_u  ) )**2       &
                               + ( arg_h_zmo_l ) * ( ln_z_z0h(m) - psi_h( arg_h_zmo_l )            &
                                                                 + psi_h( arg_h_z0h_l ) )          &
                                                 / ( ln_z_z0(m)  - psi_m( arg_m_zmo_l )            &
                                                                 + psi_m( arg_m_z0_l  ) )**2       &
                         ) / ( ol_u_vec - ol_l_vec )
             ELSE
!
!--             Calculate f = Ri - 1 /[...]^3 = 0.
                arg_h_zmo_m = z_mo(m) / ol_m_vec
                arg_m_z0_m  = z0(m)   / ol_m_vec
                f_vec = rib(m) - ( arg_h_zmo_m ) / ( ln_z_z0(m)  - psi_m( arg_h_zmo_m )            &
                                                                 + psi_m( arg_m_z0_m  ) )**3

!
!--             Calculate df/dL.
                arg_m_zmo_u = z_mo(m) / ol_u_vec
                arg_m_z0_u  = z0(m)   / ol_u_vec
                arg_m_zmo_l = z_mo(m) / ol_l_vec
                arg_m_z0_l  = z0(m)   / ol_l_vec
                f_d_ol_vec = ( - ( arg_m_zmo_u ) / ( ln_z_z0(m)  - psi_m( arg_m_zmo_u )            &
                                                                 + psi_m( arg_m_z0_u  ) )**3       &
                               + ( arg_m_zmo_l ) / ( ln_z_z0(m)  - psi_m( arg_m_zmo_l )            &
                                                                 + psi_m( arg_m_z0_l  ) )**3       &
                             ) / ( ol_u_vec - ol_l_vec )
             ENDIF
!
!--          Calculate new L.
             ol(m) = ol_m_vec - f_vec / f_d_ol_vec

!
!--          Ensure that the bulk Richardson number and the Obukhov length have the same sign and
!--          ensure convergence. If the sign is not the same, the above calculated Obukhov length
!--          obviously overshooted to the opposite side, so the next iteration should start with
!--          a smaller value.
             IF ( ol(m) * ol_m_vec < 0.0_wp )  ol(m) = ol_m_vec * 0.5_wp

!
!--          Check for convergence.
!--          This check does not modify ol, therefore this is done first.
             IF ( ABS( ( ol(m) - ol_m_vec ) /  ol(m) ) < ol_tol )  THEN
                convergence_reached(m) = .TRUE.
             ENDIF
!
!--          In the deep neutral zone, set L to the maximum allowed value.
             IF ( ABS( ol(m) ) > ol_max )  THEN
                ol(m) = SIGN( ol_max, ol(m) )
                convergence_reached(m) = .TRUE.
             ENDIF
          ENDDO
!
!--       Assure that Obukhov length does not become zero.
          !$ACC PARALLEL LOOP &
          !$ACC PRESENT(surf, ns, ol) DEFAULT(NONE) IF(enable_openacc)
          DO  m = 1, surf%ns
             IF ( convergence_reached(m) )  CYCLE
             IF ( ABS( ol(m) ) < ol_min )  THEN
                ol(m) = SIGN( ol_min, ol(m) )
                convergence_reached(m) = .TRUE.
             ENDIF
          ENDDO

          !$ACC UPDATE HOST(convergence_reached) IF(enable_openacc)
          IF ( ALL( convergence_reached ) )  EXIT

       ENDDO  ! End of iteration loop

    ENDIF  ! End of vector branch

    !$ACC END DATA

 END SUBROUTINE calc_ol


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the bulk Richardson number in case of prescribed sensible/latent surface heat fluxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rib_with_prescribed_fluxes( ns, downward, k, koff, pt1, rib, shf, upward,         &
                                             uvw_abs, z_mo, qsws, qv1, vpt1 )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  ns  !< number of surface points

    INTEGER(iwp), INTENT(IN), DIMENSION(ns) ::  k     !< z-index linking surface to prognostic grid point in the PALM 3D-grid
    INTEGER(iwp), INTENT(IN), DIMENSION(ns) ::  koff  !< offset value in z indicating the position
                                                      !< of the surface with respect to the reference grid point

    LOGICAL, INTENT(IN), DIMENSION(ns) ::  downward  !< flag indicating downward-facing surfaces
    LOGICAL, INTENT(IN), DIMENSION(ns) ::  upward    !< flag indicating upward-facing surfaces

    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  pt1      !< potential temperature at first grid level
    REAL(wp), INTENT(OUT), DIMENSION(ns) ::  rib      !< Richardson flux number
    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  shf      !< surface sensible heat flux
    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  uvw_abs  !< absolute surface-parallel velocity on grid center
    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  z_mo     !< constant flux layer height

    REAL(wp), INTENT(IN),  DIMENSION(ns), OPTIONAL ::  qsws  !< surface latent heat flux
    REAL(wp), INTENT(IN),  DIMENSION(ns), OPTIONAL ::  qv1   !< mixing ratio at first grid level
    REAL(wp), INTENT(IN),  DIMENSION(ns), OPTIONAL ::  vpt1  !< virtual potential temperature at first grid level

    IF ( humidity )  THEN
       !$OMP PARALLEL DO
       !$ACC PARALLEL LOOP GANG &
       !$ACC PRESENT(downward, drho_air_zw, k, koff, pt1, qsws, qv1, rib, shf, upward, uvw_abs, vpt1, z_mo) &
       !$ACC DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, ns
          rib(m) = -g * z_mo(m) *                                                                  &
                   ( ( 1.0_wp + 0.61_wp * qv1(m) ) * shf(m) + 0.61_wp * pt1(m) * qsws(m) ) /       &
                   ( uvw_abs(m)**3 * vpt1(m) * kappa**2 + 1.0E-20_wp )
!
!--       Note, at upward or downward-facing surfaces the sensible and latent fluxes include
!--       density. To make Rib dimensionless, multiply with 1/density.
          rib(m) = MERGE( rib(m) * drho_air_zw(k(m)+koff(m)), rib(m), upward(m) .OR. downward(m) )
       ENDDO
       !$ACC END PARALLEL LOOP
    ELSE
       !$OMP PARALLEL DO
       !$ACC PARALLEL LOOP &
       !$ACC PRESENT(downward, drho_air_zw, k, koff, ns, pt1, rib, shf, upward, uvw_abs, z_mo) &
       !$ACC DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, ns
          rib(m) = -g * z_mo(m) * shf(m) / ( uvw_abs(m)**3 * pt1(m) * kappa**2 + 1.0E-20_wp )
!
!--       Note, at upward or downward-facing surfaces the sensible flux includes density.
!--       To make Rib dimensionless, multiply with 1/density.
          rib(m) = MERGE( rib(m) * drho_air_zw(k(m)+koff(m)), rib(m), upward(m) .OR. downward(m) )
       ENDDO
    ENDIF

 END SUBROUTINE calc_rib_with_prescribed_fluxes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the bulk Richardson number for given surface (z0) temperature.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rib( ns, pt1, pt_surface, rib, uvw_abs, z_mo, slurb )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  ns  !< number of surface points

    LOGICAL, OPTIONAL, INTENT(IN) ::  slurb  !< indicates if routine is called from SLUrb model

    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  pt1          !< potential temperature at first grid level
    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  pt_surface   !< skin-surface potential temperature
    REAL(wp), INTENT(OUT), DIMENSION(ns) ::  rib          !< Richardson flux number
    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  uvw_abs      !< absolute surface-parallel velocity on grid center
    REAL(wp), INTENT(IN),  DIMENSION(ns) ::  z_mo         !< constant flux layer height


!
!-- Evaluate bulk Richardson number.
    !$OMP PARALLEL DO
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT(ns, pt1, pt_surface, rib, uvw_abs, z_mo) IF(enable_openacc)
    DO  m = 1, ns
       rib(m) = g * z_mo(m) * ( pt1(m) - pt_surface(m) ) / ( uvw_abs(m)**2 * pt1(m) + 1.0E-20_wp )
    ENDDO

!
!-- For the SLUrb model, limit to |rib| < |rib_max| to dampen possible instabilities during
!-- initialization.
    IF ( PRESENT( slurb ) )  THEN
       IF ( slurb )  THEN
          DO  m = 1, ns
             IF ( ABS( rib(m) ) > rib_max )  rib(m) = SIGN( rib_max, rib(m) )
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE calc_rib


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate friction velocity u* representative for grid center.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_us_s

    IMPLICIT NONE

    !$OMP PARALLEL  DO PRIVATE( z_mo )
    !$ACC PARALLEL LOOP PRIVATE( z_mo ) &
    !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns

       z_mo = surf%z_mo(m)
!
!--    Compute u* at the scalars' grid points. At horizonally upward-facing surfaces or surfaces
!--    with less than 30 degrees slope, use stability correction, else take purely neutral solution.
       surf%us(m) = MERGE( kappa * surf%uvw_abs(m) / ( surf%ln_z_z0(m)                             &
                           - psi_m( z_mo / surf%ol(m) ) + psi_m( surf%z0(m) / surf%ol(m) ) ),      &
                           kappa * surf%uvw_abs(m) / surf%ln_z_z0(m),                              &
                           surf%consider_stability(m) )
    ENDDO

 END SUBROUTINE calc_us_s


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate friction velocity u* representative for the u- or v-grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_us_uv

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( z_mo )
    !$ACC PARALLEL LOOP PRIVATE( z_mo ) &
    !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       z_mo = surf%z_mo(m)
!
!--    Compute u* at the uv grid points. At horizonally upward-facing surfaces or surfaces
!--    with less than 30 degrees slope, use stability correction, else take purely neutral solution.
       surf%us_uvgrid(m) = MERGE( kappa * surf%uvw_abs_uv(m) / ( surf%ln_z_z0(m)                   &
                                - psi_m( z_mo / surf%ol(m) ) + psi_m( surf%z0(m) / surf%ol(m) ) ), &
                                  kappa * surf%uvw_abs_uv(m) / surf%ln_z_z0(m),                    &
                                  surf%consider_stability(m) )

    ENDDO

 END SUBROUTINE calc_us_uv


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate friction velocity u* representative for the w-grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_us_w

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( z_mo )
    !$ACC PARALLEL LOOP PRIVATE( z_mo ) &
    !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       z_mo = surf%z_mo(m)
!
!--    Compute u* at the w grid points. At horizonally upward-facing surfaces or surfaces
!--    with less than 30 degrees slope, use stability correction, else take purely neutral solution.
       surf%us_wgrid(m) = MERGE( kappa * surf%uvw_abs_w(m) / ( surf%ln_z_z0(m)                     &
                                - psi_m( z_mo / surf%ol(m) ) + psi_m( surf%z0(m) / surf%ol(m) ) ), &
                                  kappa * surf%uvw_abs_w(m) / surf%ln_z_z0(m),                     &
                                  surf%consider_stability(m) )
    ENDDO

 END SUBROUTINE calc_us_w

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature, specific humidity, and virtual potential temperature at first
!> grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_pt_q

    IMPLICIT NONE
!
!-- @todo: Following loop need to be split-up.

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
    !$ACC PRESENT(pt, q, surf) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)

#if ! defined( _OPENACC )
       IF ( bulk_cloud_model ) THEN
          surf%pt1(m) = pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i)
          surf%qv1(m) = q(k,j,i) - ql(k,j,i)
       ELSEIF( cloud_droplets ) THEN
          surf%pt1(m) = pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i)
          surf%qv1(m) = q(k,j,i)
       ELSE
          surf%pt1(m) = pt(k,j,i)
          IF ( humidity )  THEN
             surf%qv1(m) = q(k,j,i)
          ELSE
             surf%qv1(m) = 0.0_wp
          ENDIF
       ENDIF

       IF ( humidity )  THEN
          surf%vpt1(m) = pt(k,j,i) * ( 1.0_wp + 0.61_wp * q(k,j,i) )
       ENDIF
#else
       surf%pt1(m) = pt(k,j,i)
       IF ( humidity )  THEN
          surf%qv1(m)  = q(k,j,i)
          surf%vpt1(m) = pt(k,j,i) * ( 1.0_wp + 0.61_wp * q(k,j,i) )
       ELSE
          surf%qv1(m) = 0.0_wp
       ENDIF
#endif
    ENDDO

 END SUBROUTINE calc_pt_q


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store potential temperature at surface grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE store_pt_surface

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k ) &
    !$ACC PRESENT(pt, surf) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%i(m) + surf%ioff(m)
       j = surf%j(m) + surf%joff(m)
       k = surf%k(m) + surf%koff(m)
       surf%pt_surface(m) = pt(k,j,i)
    ENDDO

 END SUBROUTINE store_pt_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store mixing ratio at surface grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE store_q_surface

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k ) &
    !$ACC PRESENT(q, surf, surf%ns, surf%i, surf%ioff, surf%j, surf%joff, surf%k, surf%koff, surf%q_surface) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%i(m) + surf%ioff(m)
       j = surf%j(m) + surf%joff(m)
       k = surf%k(m) + surf%koff(m)
       surf%q_surface(m) = q(k,j,i)
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE store_q_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store virtual potential temperature at surface grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE store_vpt_surface

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k ) &
    !$ACC PRESENT(surf, surf%ns, surf%i, surf%ioff, surf%j, surf%joff, surf%k, surf%koff, surf%vpt_surface, vpt) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%i(m) + surf%ioff(m)
       j = surf%j(m) + surf%joff(m)
       k = surf%k(m) + surf%koff(m)
       surf%vpt_surface(m) = vpt(k,j,i)
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE store_vpt_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the other MOST scaling parameters theta*, q*, (qc*, qr*, nc*, nr*)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_scaling_parameters

    IMPLICIT NONE

    INTEGER(iwp)  ::  lsp   !< running index for chemical species

    LOGICAL  ::  lsm_switch !<

!
!-- Compute theta* at horizontal surfaces
    IF ( constant_heatflux )  THEN
!
!--    For a given heat flux in the surface layer:

       !$OMP PARALLEL DO PRIVATE( k, k_off )
       !$ACC PARALLEL LOOP PRIVATE( k, k_off ) &
       !$ACC PRESENT( surf, drho_air_zw) DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
!
!--       Compute ts for horizontally upward-facing surfaces or surfaces with less than 30 degrees
!--       slope.
          surf%ts(m) = MERGE( -surf%shf(m) * drho_air_zw(k+k_off) / ( surf%us(m) + 1E-30_wp ),     &
                              surf%ts(m),                                                          &
                              surf%consider_stability(m) )
!
!--       ts must be limited, because otherwise overflow may occur in case of us=0 when computing
!--       ol further below.
          IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
          IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp
       ENDDO

    ELSE
!
!--    For a given surface temperature:
       IF ( large_scale_forcing  .AND.  lsf_surf )  THEN

          !$OMP PARALLEL DO PRIVATE( i, i_off, j, j_off, k, k_off )
          DO  m = 1, surf%ns
             i   = surf%i(m)
             j   = surf%j(m)
             k   = surf%k(m)
             i_off = surf%ioff(m)
             j_off = surf%joff(m)
             k_off = surf%koff(m)
!
!--          @todo: This need to be changed to surf%pt_surface later.
             pt(k+k_off,j+j_off,i+i_off) = pt_surface
          ENDDO
       ENDIF

       !$OMP PARALLEL DO PRIVATE( z_mo )
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) &
       !$ACC COPY(surf, surf%ns, surf%z_mo, surf%pt1, surf%pt_surface) &
       !$ACC COPY(surf%ln_z_z0h, surf%ol, surf%z0h, surf%ts, surf%consider_stability) &
       !$ACC PRIVATE(z_mo) IF(enable_openacc)
       DO  m = 1, surf%ns
          z_mo = surf%z_mo(m)
!
!--       Compute ts for horizontally upward-facing surfaces or surfaces with less than 30 degrees
!--       slope.
          surf%ts(m) = MERGE( kappa * ( surf%pt1(m) - surf%pt_surface(m) )                         &
                              / ( surf%ln_z_z0h(m) - psi_h( z_mo / surf%ol(m) )                    &
                                                   + psi_h( surf%z0h(m) / surf%ol(m) ) ),          &
                              surf%ts(m),                                                          &
                              surf%consider_stability(m) )
       ENDDO
       !$ACC END PARALLEL LOOP

    ENDIF
!
!-- Compute theta* again at vertical surfaces. This is only required for natural surfaces when
!-- aerodynamical resistance are computed via MOST relations.
    lsm_switch = land_surface  .AND.  .NOT. aero_resist_kray  .AND.                                &
                 ( ALLOCATED( surf%pavement_surface )  .OR.                                        &
                   ALLOCATED( surf%vegetation_surface )  .OR.                                      &
                   ALLOCATED( surf%water_surface ) )
    !$OMP PARALLEL DO
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT( surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
!
!--    Save already computed values at horizontal surfaces or surfaces with less than 30 degrees
!--    slope (considered by flag consider_stability), otherwise compute t* from the heat-flux and
!--    friction velocity.
       surf%ts(m) = MERGE( -surf%shf(m) / ( surf%us(m) + 1E-30_wp ),                               &
                           surf%ts(m),                                                             &
                           .NOT. ( surf%upward(m)  .OR.  surf%downward(m) )  .AND.  lsm_switch     &
                           .AND.  .NOT. surf%consider_stability(m) )
!
!--    ts must be limited, because otherwise overflow may occur in case of us=0 when computing ol
!--    further below
       IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
       IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp
    ENDDO

!
!-- If required compute q* at horizontal surfaces or surfaces with less than 30 degrees slope.
    IF ( humidity )  THEN
       IF ( constant_waterflux )  THEN
!
!--       For a given water flux in the surface layer
          !$OMP PARALLEL DO PRIVATE( k, k_off )
          !$ACC PARALLEL LOOP GANG &
          !$ACC PRESENT(drho_air_zw, surf, surf%consider_stability, surf%k, surf%koff, surf%ns, surf%qs, surf%qsws, surf%us) &
          !$ACC DEFAULT(NONE) IF(enable_openacc)
          DO  m = 1, surf%ns
             k = surf%k(m)
             k_off = surf%koff(m)
             surf%qs(m) = MERGE( -surf%qsws(m) * drho_air_zw(k+k_off) / ( surf%us(m) + 1E-30_wp ), &
                                 surf%qs(m),                                                       &
                                 surf%consider_stability(m) )
          ENDDO
          !$ACC END PARALLEL LOOP

       ELSE

          IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
             !$OMP PARALLEL DO PRIVATE( i, i_off, j, j_off, k, k_off )
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = surf%k(m)
                i_off = surf%ioff(m)
                j_off = surf%joff(m)
                k_off = surf%koff(m)
                q(k+k_off,j+j_off,i+i_off) = q_surface
             ENDDO
          ENDIF

!
!--       Assume saturation for atmosphere coupled to ocean.
          IF ( atmosphere_run_coupled_to_ocean )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, e_s )
             DO  m = 1, surf%ns
                i   = surf%i(m)
                j   = surf%j(m)
                k   = surf%k(m)
                e_s = magnus( exner(k-1) * pt(k-1,j,i) )
                q(k-1,j,i) = rd_d_rv * e_s / ( 100.0_wp * surface_pressure - e_s )
             ENDDO
          ENDIF

          !$OMP PARALLEL DO PRIVATE( k, k_off, z_mo )
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) &
          !$ACC COPY(surf, surf%ns, surf%k,  surf%koff, surf%z_mo, surf%qs, surf%qv1) &
          !$ACC COPY(surf%q_surface, surf%ln_z_z0q, surf%ol, surf%z0q, surf%consider_stability) &
          !$ACC PRIVATE(k, k_off, z_mo) IF(enable_openacc)
          DO  m = 1, surf%ns
             k = surf%k(m)
             k_off = surf%koff(m)
             z_mo = surf%z_mo(m)
             surf%qs(m) = MERGE( kappa * ( surf%qv1(m) - surf%q_surface(m) )                       &
                                / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                  &
                                                     + psi_h( surf%z0q(m) / surf%ol(m) ) ),        &
                                 surf%qs(m),                                                       &
                                 surf%consider_stability(m) )
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF
!
!--    Compute q* at vertical surfaces or surfaces with more than 30 degrees slope.
       !$OMP PARALLEL DO
       !$ACC PARALLEL LOOP GANG &
       !$ACC PRESENT(surf, surf%consider_stability, surf%downward, surf%ns, surf%qs, surf%qsws, surf%upward, surf%us) &
       !$ACC DEFAULT(NONE) IF(enable_openacc)
       DO  m = 1, surf%ns
!
!--       Save already computed values at horizontal surfaces.
          surf%qs(m) = MERGE( -surf%qsws(m) / ( surf%us(m) + 1E-30_wp ),                           &
                              surf%qs(m),                                                          &
                              .NOT. ( surf%upward(m)  .OR.  surf%downward(m) )   .AND.             &
                              .NOT. surf%consider_stability(m) )
       ENDDO
       !$ACC END PARALLEL LOOP
    ENDIF

!
!-- If required compute s*.
    IF ( passive_scalar )  THEN

       IF ( constant_scalarflux  )  THEN
!
!--       For a given scalar flux in the surface layer and at horizontal surfaces or surfaces with
!--       less than 30 degrees slope.
          !$OMP PARALLEL DO
          DO  m = 1, surf%ns
             surf%ss(m) = MERGE( -surf%ssws(m) / ( surf%us(m) + 1E-30_wp ), surf%ss(m),            &
                                 surf%consider_stability(m) )
          ENDDO
       ELSE

          !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
          DO  m = 1, surf%ns
             i = surf%iref(m)
             j = surf%jref(m)
             k = surf%kref(m)
             k_off = surf%koff(m)
             z_mo = surf%z_mo(m)

             surf%ss(m) = MERGE( kappa * ( s(k,j,i) - s(surf%k(m)+k_off,j,i) )                     &
                                 / ( surf%ln_z_z0h(m) - psi_h( z_mo / surf%ol(m) )                 &
                                               + psi_h( surf%z0h(m) / surf%ol(m) ) ),              &
                                 surf%ss(m),                                                       &
                                 surf%consider_stability(m) )
          ENDDO
       ENDIF
!
!--    Treat vertical surfaces or surfaces with more than 30 degrees slope.
       !$OMP PARALLEL DO
       DO  m = 1, surf%ns
!
!--       Save already computed values at horizontal surfaces.
          surf%ss(m) = MERGE( -surf%ssws(m) / ( surf%us(m) + 1E-30_wp ),                           &
                              surf%ss(m),                                                          &
                              .NOT. ( surf%upward(m)  .OR.  surf%downward(m) )   .AND.             &
                              .NOT. surf%consider_stability(m) )
       ENDDO
    ENDIF

!
!-- If required compute cs* (chemical species).
    IF ( air_chemistry  )  THEN
!
!--    Compute scaling parameters only at horizontal surfaces or surfaces with less than 30
!--    degrees slope, otherwise take the initial values.
       DO  lsp = 1, nvar
          IF ( constant_csflux(lsp) )  THEN
!--          For a given chemical species' flux in the surface layer
             !$OMP PARALLEL DO
             DO  m = 1, surf%ns
                surf%css(lsp,m) = MERGE( -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp ),           &
                                         surf%css(lsp,m),                                          &
                                         surf%consider_stability(m) )
             ENDDO
          ENDIF
       ENDDO
!
!--    Treat vertical surfaces or surfaces with more than 30 degrees slope.
       DO  lsp = 1, nvar
          !$OMP PARALLEL DO
          DO  m = 1, surf%ns
!
!--          Save already computed values at horizontal surfaces.
             surf%css(lsp,m) = MERGE( -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp ),              &
                                      surf%css(lsp,m),                                             &
                                      .NOT. ( surf%upward(m)  .OR.  surf%downward(m) )   .AND.     &
                                      .NOT. surf%consider_stability(m) )
          ENDDO
       ENDDO
    ENDIF

!
!-- If required compute qc* and nc*. Compute scaling parameters only at horizontal surfaces or
!-- surfaces with less than 30 degrees slope, otherwise take the initial values.
    IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
       DO  m = 1, surf%ns
          i = surf%iref(m)
          j = surf%jref(m)
          k = surf%kref(m)
          k_off = surf%koff(m)

          z_mo = surf%z_mo(m)

          surf%qcs(m) = MERGE( kappa * ( qc(k,j,i) - qc(surf%k(m)+k_off,j,i) )                     &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%qcs(m),                                                        &
                               surf%consider_stability(m) )

          surf%ncs(m) = MERGE( kappa * ( nc(k,j,i) - nc(surf%k(m)+k_off,j,i) )                     &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%ncs(m),                                                        &
                               surf%consider_stability(m) )
       ENDDO

    ENDIF

!
!-- If required compute qr* and nr*. Compute scaling parameters only at horizontal surfaces or
!-- surfaces with less than 30 degrees slope, otherwise take the initial values.
    IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
       DO  m = 1, surf%ns
          i = surf%iref(m)
          j = surf%jref(m)
          k = surf%kref(m)
          k_off = surf%koff(m)

          z_mo = surf%z_mo(m)

          surf%qrs(m) = MERGE( kappa * ( qr(k,j,i) - qr(surf%k(m)+k_off,j,i) )                     &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%qrs(m),                                                        &
                               surf%consider_stability(m) )

          surf%nrs(m) = MERGE( kappa * ( nr(k,j,i) - nr(surf%k(m)+k_off,j,i) )                     &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%nrs(m),                                                        &
                               surf%consider_stability(m) )
       ENDDO

    ENDIF

 END SUBROUTINE calc_scaling_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux usws which later enters the u-equation. Note, usws is actually only
!> required at upward-facing surfaces, though the loops runs over all surfaces. Depending on the
!> topography representation, the pre-computed friction velocity is either defined on the grid
!> center (Cartesian topography) or at the u-grid (cut-cell topography).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_usws( us_tmp )

    IMPLICIT NONE

    INTEGER(iwp)  ::  k_surf  !< vertical index of the surface

    REAL(wp)      ::  denom   !< denominator including ln(z/z0) + integrated profile functions to account for stability
    REAL(wp)      ::  u_comp  !< u-component

    REAL(wp), DIMENSION(1:surf%ns) ::  us_tmp !< passed friction velocity corresponding to the s- or u-grid

!
!-- Compute u'w'
    !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, k_surf, z_mo, u_comp, denom )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k, k_off, k_surf, z_mo, u_comp, denom ) &
    !$ACC PRESENT(surf, u, rho_air_zw) PRESENT(us_tmp) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)
       k_surf = surf%k(m) + surf%koff(m)
       k_off = surf%koff(m)
!
!--    Depending on the orientation of the surface the resulting u-component is either defined
!--    as a relative velocity to the surface (in case of coupled atmosphere-ocean runs), or a
!--    absolute value (at downward-facing surfaces). The same is done for the stability
!--    correction which is only considered at upward-facing walls.
       z_mo    = surf%z_mo(m)
       u_comp  = MERGE( u(k,j,i) - u(k_surf,j,i), u(k,j,i), surf%upward(m) )
       denom   = MERGE( surf%ln_z_z0(m) - psi_m( z_mo / surf%ol(m) )                               &
                                        + psi_m( surf%z0(m) / surf%ol(m) ),                        &
                        surf%ln_z_z0(m),                                                           &
                        surf%consider_stability(m) )
!
!--    Please note, the computation of usws is not fully accurate in case of step-like topography.
!--    Actually a further interpolation of ol onto the u-grid, where usws is defined, is required.
!--    However, this is not done since this would require several data transfers between the
!--    surface-data structures. To account for different facings (up/downward), multiply with the
!--    respective normal-vector component for the corresponding facing.
       surf%usws(m) = kappa * u_comp / denom
       surf%usws(m) = -surf%usws(m) * us_tmp(m) * rho_air_zw(k_surf) * surf%n_eff(m)
    ENDDO
!
!-- Mask usws at vertically bounded grid points.
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       surf%usws(m) = MERGE( surf%usws(m), 0.0_wp, surf%upward(m)  .OR.  surf%downward(m) )
    ENDDO

 END SUBROUTINE calc_usws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux vsws which later enters the v-equation. Note, vsws is actually only
!> required at upward-facing surfaces, though the loops runs over all surfaces. Depending on the
!> topography representation, the pre-computed friction velocity is either defined on the grid
!> center (Cartesian topography) or at the v-grid (cut-cell topography).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_vsws( us_tmp )

    IMPLICIT NONE

    INTEGER(iwp)  ::  k_surf  !< vertical index of the surface

    REAL(wp)      ::  denom   !< denominator including ln(z/z0) + integrated profile functions to account for stability
    REAL(wp)      ::  v_comp  !< v-component

    REAL(wp), DIMENSION(1:surf%ns) ::  us_tmp !< passed friction velocity corresponding to the s- or u-grid
!
!-- Compute v'w'
    !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, k_surf, z_mo, v_comp, denom )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k, k_off, k_surf, z_mo, v_comp, denom ) &
    !$ACC PRESENT(surf, v, rho_air_zw) PRESENT(us_tmp) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)
       k_surf = surf%k(m) + surf%koff(m)
       k_off = surf%koff(m)
!
!--    Depending on the orientation of the surface the resulting v-component is either defined
!--    as a relative velocity to the surface (in case of coupled atmosphere-ocean runs), or a
!--    absolute value (at downward-facing surfaces). The same is done for the stability
!--    correction which is only considered at upward-facing walls.
       z_mo    = surf%z_mo(m)
       v_comp  = MERGE( v(k,j,i) - v(k_surf,j,i), v(k,j,i), surf%upward(m) )
       denom   = MERGE( surf%ln_z_z0(m) - psi_m( z_mo / surf%ol(m) )                               &
                                        + psi_m( surf%z0(m) / surf%ol(m) ),                        &
                        surf%ln_z_z0(m),                                                           &
                        surf%consider_stability(m) )
!
!--    Please note, the computation of vsws is not fully accurate in case of step-like topography.
!--    Actually a further interpolation of ol onto the v-grid, where vsws is defined, is required.
!--    However, this is not done since this would require several data transfers between the
!--    surface-data structures. To account for different facings (up/downward), multiply with the
!--    respective normal-vector component for the corresponding facing.
       surf%vsws(m) = kappa * v_comp / denom
       surf%vsws(m) = -surf%vsws(m) * us_tmp(m) * rho_air_zw(k_surf) * surf%n_eff(m)
    ENDDO
!
!-- Mask vsws at vertically bounded grid points.
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT(surf) DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       surf%vsws(m) = MERGE( surf%vsws(m), 0.0_wp, surf%upward(m)  .OR.  surf%downward(m) )
    ENDDO

 END SUBROUTINE calc_vsws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux vsus at vertical surfaces which later enter the u-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_usvs

    IMPLICIT NONE

    REAL(wp) ::  flag_u  !< flag indicating u-grid, used for calculation of horizontal momentum fluxes at vertical surfaces


!
!-- Generalize computation by introducing flags. At north- and south-facing surfaces
!-- u-component is used, at east- and west-facing surfaces v-component is used.
    !$OMP PARALLEL  DO PRIVATE( i, j, k, flag_u )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE(i, j, k) &
    !$ACC PRESENT(surf, surf%iref, surf%jref, surf%kref, surf%ln_z_z0, surf%northward, surf%ns, surf%n_eff, surf%southward, surf%usvs, surf%us_uvgrid, u) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)

       flag_u = MERGE( 1.0_wp, 0.0_wp, surf%northward(m)  .OR.  surf%southward(m) )
!
!--    Compute the fluxes for the horizontal transport of u at vertical surfaces.
!--    Mask fluxes at non-relevant, horizontal surface-orientations by multiplication with 0.
       surf%usvs(m) = kappa * ( flag_u * u(k,j,i) ) / surf%ln_z_z0(m)
       surf%usvs(m) = -surf%usvs(m) * surf%us_uvgrid(m) * surf%n_eff(m)
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE calc_usvs


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux usvs and vsus at vertical surfaces which later enter the v-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_vsus

    IMPLICIT NONE

    REAL(wp) ::  flag_v  !< flag indicating v-grid, used for calculation of horizontal momentum fluxes at vertical surfaces

!
!-- Generalize computation by introducing flags. At north- and south-facing surfaces
!-- u-component is used, at east- and west-facing surfaces v-component is used.
    !$OMP PARALLEL  DO PRIVATE( i, j, k, flag_v )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE(i, j, k) &
    !$ACC PRESENT(surf, surf%eastward, surf%iref, surf%jref, surf%kref, surf%ln_z_z0, surf%ns, surf%n_eff, surf%us_uvgrid, surf%vsus, surf%westward, v) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)

       flag_v = MERGE( 1.0_wp, 0.0_wp, surf%eastward(m)  .OR.  surf%westward(m) )
!
!--    Compute the fluxes for the horizontal transport of v at vertical surfaces.
!--    Mask fluxes at non-relevant, horizontal surface-orientations by multiplication with 0.
       surf%vsus(m) = kappa * ( flag_v * v(k,j,i) ) / surf%ln_z_z0(m)
       surf%vsus(m) = -surf%vsus(m) * surf%us_uvgrid(m) * surf%n_eff(m)
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE calc_vsus


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux wsus and wsvs at vertical surfaces which later enter the w-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_wsus_wsvs

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE(i, j, k) &
    !$ACC PRESENT(surf, surf%eastward, surf%iref, surf%jref, surf%kref, surf%ln_z_z0, surf%northward, surf%ns, surf%southward, surf%n_eff, surf%us_wgrid, surf%westward, surf%wsus_wsvs, w) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)
!
!--    Compute the fluxes for the horizontal transport of w at vertical surfaces. Mask horizontal
!--    surfaces by multiplication with 0.
       surf%wsus_wsvs(m) = kappa * w(k,j,i) / surf%ln_z_z0(m)
       surf%wsus_wsvs(m) = -surf%wsus_wsvs(m) * surf%us_wgrid(m) * surf%n_eff(m) *                 &
                           MERGE( 1.0_wp, 0.0_wp, surf%northward(m)  .OR.  surf%southward(m)  .OR. &
                                                  surf%eastward(m)   .OR.  surf%westward(m) )
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE calc_wsus_wsvs


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux usws and vsws on scalar grid. These fluxes later enter the
!> SGS-TKE-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_usws_vsws_for_tke

    IMPLICIT NONE

    REAL(wp) ::  dum     !< dummy to precalculate logarithm
    REAL(wp) ::  u_comp  !< u-component interpolated onto scalar grid point
    REAL(wp) ::  v_comp  !< v-component interpolated onto scalar grid point
    REAL(wp) ::  w_comp  !< w-component interpolated onto scalar grid point


!
!-- Compute momentum fluxes used for subgrid-scale TKE production at vertical surfaces.
!-- Fluxes required for the TKE prognostic are located at the grid center.
!-- Please note, the note runs over all surfaces but the resulting fluxes will effectively computed
!-- only at vertical surfaces.
    !$OMP PARALLEL DO PRIVATE( i, j, k, dum, u_comp, v_comp, w_comp )
    !$ACC PARALLEL LOOP GANG &
    !$ACC PRIVATE(i, j, k, u_comp, v_comp, w_comp) &
    !$ACC PRESENT(surf, surf%downward, surf%eastward, surf%iref, surf%jref, surf%kref, surf%ln_z_z0) &
    !$ACC PRESENT(surf%mom_flux_tke, surf%ns, surf%northward, surf%southward, surf%tke_production) &
    !$ACC PRESENT(surf%upward, surf%us, surf%westward, u, v, w) &
    !$ACC DEFAULT(NONE) IF(enable_openacc)
    DO  m = 1, surf%ns
       i = surf%iref(m)
       j = surf%jref(m)
       k = surf%kref(m)

       u_comp = MERGE( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ), 0.0_wp,                                 &
                       surf%northward(m)  .OR.  surf%southward(m) )
       v_comp = MERGE( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ), 0.0_wp,                                 &
                       surf%eastward(m)  .OR.  surf%westward(m) )
       w_comp = MERGE( 0.5_wp * ( w(k,j,i) + w(k-1,j,i) ), 0.0_wp,                                 &
                       .NOT. ( surf%upward(m)  .OR.  surf%downward(m) ) )

       dum = kappa / surf%ln_z_z0(m)
!
!--    usvs at north/southward-facing walls (joff/=0) or vsus at
!--    east/westward-facing walls (ioff/=0).
       surf%mom_flux_tke(0,m) = dum * ( u_comp + v_comp )
!
!--    wsvs at north/southward-facing walls (joff/=0) or wsus at
!--    east/westward-facing walls (ioff/=0).
       surf%mom_flux_tke(1,m) = dum * w_comp
!
!--    Finally, compute the momentum fluxes for TKE production. Note, if the distance between
!--    the prognostic grid point and the surface are too close in comparison to the "usual"
!--    half-grid distance, the momentum fluxes are simply set to zero. This is done to avoid
!--    an unrealistic TKE production by shear, leading to a blow-up of TKE values.
       surf%mom_flux_tke(0:1,m) = -surf%mom_flux_tke(0:1,m) * surf%us(m) *                         &
                                   MERGE( 1.0_wp, 0.0_wp, surf%tke_production(m) )
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE calc_usws_vsws_for_tke


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate aerodynamic resistance.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_aerodynamic_resistance

    REAL(wp) ::  ueff  !< limited near-surface wind speed


!
!-- Compute aerodynamic resistance.
    !$OMP PARALLEL DO PRIVATE (m, i, j, k, ueff ) SCHEDULE (STATIC)
    DO  m = 1, surf%ns

       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)
!
!--    Calculate aerodyamic resistance. At upward facing surfaces, use MOST and the vertical
!--    temperature gradient.
       IF ( surf%upward(m)  .OR.  .NOT. aero_resist_kray )  THEN
!
!--       Dimension of r_a is s/m: K / ( K * m / s )
          surf%r_a(m) = ABS( ( surf%pt1(m) - surf%pt_surface(m) ) /                                &
                             ( surf%ts(m) * surf%us(m) + 1.0E-20_wp ) )
!
!--    At surfaces with other orientation, use the formulation used in the TUF3d model.
!--    (Krayenhoff & Voogt, 2006). Note that this formulation is the equivalent to the
!--    ECMWF formulation using drag coefficients. A roughness length of 0.001 is assumed here for
!--    concrete (the inverse, 1000 is used in the nominator for scaling).
       ELSE
!
!--       Limit wind velocity in order to avoid division by zero.
!--       The nominator can become <= 0.0 for values z0 < 3*10E-4.
          ueff = MAX( SQRT(  ( ( u(k,j,i) + u(k,j,i+1) ) * 0.5_wp )**2 +                           &
                             ( ( v(k,j,i) + v(k,j+1,i) ) * 0.5_wp )**2 +                           &
                             ( ( w(k,j,i) + w(k-1,j,i) ) * 0.5_wp )**2 ),                          &
                      1.0_wp / 4.2_wp * ( 4.0_wp / ( surf%z0(m) * 1000.0_wp ) - 11.8_wp ), 0.1_wp )
!
!--       Dimension of r_a is s/m with denominator in W / ( m2 K )
          surf%r_a(m) = rho_cp / ( surf%z0(m) * 1000.0_wp * ( 11.8_wp + 4.2_wp * ueff ) - 4.0_wp )
       ENDIF

!
!--    Make sure that the resistance does not drop to zero for neutral stratification. Also, set a
!--    maximum resistance to avoid the breakdown of MOST for locations with zero wind speed.
       IF ( surf%r_a(m) <   1.0_wp )  surf%r_a(m) =   1.0_wp
       IF ( surf%r_a(m) > 300.0_wp )  surf%r_a(m) = 300.0_wp
    ENDDO

 END SUBROUTINE calc_aerodynamic_resistance


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface fluxes usws, vsws, shf, qsws, (qcsws, qrsws, ncsws, nrsws)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_surface_fluxes

    IMPLICIT NONE

    INTEGER(iwp)  ::  lsp   !< running index for chemical species

    LOGICAL       ::  compute_shf    !< control flag for computation of surface sensible heat flux
    LOGICAL       ::  compute_qsws   !< control flag for computation of surface latent heat flux
    LOGICAL       ::  compute_ssws   !< control flag for computation of surface passive scalar flux
    LOGICAL       ::  compute_qcsws  !< control flag for computation of surface cloud-water content flux (is this necessary?)
    LOGICAL       ::  compute_qrsws  !< control flag for computation of surface cloud-water content flux (is this necessary?)


!
!-- Set control flags to decide whether fluxes need to be computed or not.
    compute_shf = ( .NOT.  constant_heatflux  .AND.  .NOT. neutral  .AND.                          &
                    ( ( time_since_reference_point <=  lsm_start_time  .AND.                       &
                        simulated_time > 0.0_wp )  .OR.  .NOT.  land_surface )  .AND.              &
                    .NOT.  urban_surface )
    compute_qsws = ( .NOT.  constant_waterflux  .AND.                                              &
                     ( ( time_since_reference_point <=  lsm_start_time  .AND.                      &
                         simulated_time > 0.0_wp )  .OR.  .NOT.  land_surface )  .AND.             &
                     .NOT.  urban_surface  .AND.  humidity )
    compute_ssws = ( .NOT.  constant_scalarflux  .AND.  passive_scalar )
    compute_qcsws = ( bulk_cloud_model  .AND.  microphysics_morrison  .AND.  humidity )
    compute_qrsws = ( bulk_cloud_model  .AND.  microphysics_seifert   .AND.  humidity )

    IF ( compute_shf )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
          surf%shf(m) = MERGE( -surf%ts(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),     &
                               0.0_wp,                                                             &
                               surf%upward(m) )
       ENDDO
    ENDIF
    IF ( compute_qsws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
          surf%qsws(m) = MERGE( -surf%qs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),    &
                                0.0_wp,                                                            &
                                surf%upward(m) )
       ENDDO
    ENDIF
    IF ( compute_ssws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
          surf%ssws(m) = MERGE( -surf%ss(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),    &
                                0.0_wp,                                                            &
                                surf%upward(m) )
       ENDDO
    ENDIF
!
!-- Compute (turbulent) fluxes of cloud water content and cloud drop conc.
    IF ( compute_qcsws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
!
!--       Compute (turbulent) fluxes of cloud water content and cloud drop conc.
          surf%qcsws(m) = MERGE( -surf%qcs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),  &
                                 0.0_wp,                                                           &
                                 surf%upward(m) )
          surf%ncsws(m) = MERGE( -surf%ncs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),  &
                                 0.0_wp,                                                           &
                                 surf%upward(m) )
       ENDDO
    ENDIF
    IF ( compute_qrsws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)

          surf%qrsws(m) = MERGE( -surf%qrs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),  &
                                 0.0_wp,                                                           &
                                 surf%upward(m) )
          surf%nrsws(m) = MERGE( -surf%nrs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,1),  &
                                 0.0_wp,                                                           &
                                 surf%upward(m) )
       ENDDO
    ENDIF
!
!-- Compute the vertical chemical species' flux
    DO  lsp = 1, nvar
       IF (  .NOT.  constant_csflux(lsp)  .AND.  air_chemistry )  THEN
          !$OMP PARALLEL DO PRIVATE( k, k_off )
          DO  m = 1, surf%ns
             k = surf%k(m)
             k_off = surf%koff(m)
             surf%cssws(lsp,m) = MERGE( -surf%css(lsp,m) * surf%us(m) * rho_air_zw(k+k_off)        &
                                                                      * surf%n_s(m,1),             &
                                        0.0_wp,                                                    &
                                        surf%upward(m) )
          ENDDO
       ENDIF
    ENDDO

!
!-- Boundary condition for the TKE.
    IF ( ibc_e_b == 2 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j, k, i_off, j_off, k_off )
       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
          i_off = surf%ioff(m)
          j_off = surf%joff(m)
          k_off = surf%koff(m)

          e(k,j,i) = ( surf%us(m) / 0.1_wp )**2
          e(k+k_off,j+j_off,i+i_off) = e(k,j,i)
       ENDDO
    ENDIF

 END SUBROUTINE calc_surface_fluxes

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The subroutine transfers data present on surf_def, surf_lsm and surf_usm onto the surface
!> structures used for the momentum transport. This is e.g. necessary for the Obukhov length, which
!> is not computed at the staggered grid points. It is noted that the number of surfaces on the
!> momentum surface structure is not necessarily the same. In case there is a surface defined
!> on the staggered grid but not at the grid center (only possible with the cut-cell topography),
!> represented by surf_def, surf_lsm or surf_usm, the respective quantity is not transferred but
!> set to a default value. In the special case of the Obukhov length, a neutral limit value is set.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE transfer_ol

    INTEGER(iwp) ::  mm !< running index over surface types on s-grid

!
!-- Loop over surface elements on staggered grid.
    DO  m = 1, surf%ns

       i = surf%i(m)
       j = surf%j(m)
!
!--    Set default value for Obukhov length.
       surf%ol(m) = surf%z_mo(m) / zeta_min

       DO  mm = surf_def%start_index(j,i), surf_def%end_index(j,i)
          IF ( surf%kref(m) == surf_def%kref(mm) )  THEN
             surf%ol(m) = surf_def%ol(mm)
          ENDIF
       ENDDO

       DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
          IF ( surf%kref(m) == surf_lsm%kref(mm) )  THEN
             surf%ol(m) = surf_lsm%ol(mm)
          ENDIF
       ENDDO

       DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
          IF ( surf%kref(m) == surf_usm%kref(mm) )  THEN
             surf%ol(m) = surf_usm%ol(mm)
          ENDIF
       ENDDO

    ENDDO

 END SUBROUTINE transfer_ol


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Integrated stability function for momentum. 64-bit version.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION psi_m_dp( zeta )
    !$ACC ROUTINE SEQ

    USE kinds

    IMPLICIT NONE

    REAL(dp), INTENT(IN) ::  zeta   !< Stability parameter z/L

    REAL(dp) ::  psi_m_dp  !< Integrated similarity function result
    REAL(dp) ::  x         !< dummy variable

    REAL(dp), PARAMETER ::  a = 1.0_dp            !< constant
    REAL(dp), PARAMETER ::  b = 0.66666666666_dp  !< constant
    REAL(dp), PARAMETER ::  c = 5.0_dp            !< constant
    REAL(dp), PARAMETER ::  d = 0.35_dp           !< constant
    REAL(dp), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(dp), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_dp )  THEN
       x = SQRT( SQRT( 1.0_dp  - 16.0_dp * zeta ) )
       psi_m_dp = pi * 0.5_dp - 2.0_dp * ATAN( x ) + LOG( ( 1.0_dp + x )**2                        &
                  * ( 1.0_dp + x**2 ) * 0.125_dp )
    ELSE

       psi_m_dp = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta - bc_d_d
!
!--    Old version for stable conditions (only valid for z/L < 0.5) psi_m = - 5.0_wp * zeta

    ENDIF

 END FUNCTION psi_m_dp


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Integrated stability function for momentum. 32-bit version.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION psi_m_sp( zeta )
    !$ACC ROUTINE SEQ

    USE kinds

    IMPLICIT NONE

    REAL(sp), INTENT(IN) ::  zeta   !< Stability parameter z/L

    REAL(sp) ::  psi_m_sp  !< Integrated similarity function result
    REAL(sp) ::  x         !< dummy variable

    REAL(sp), PARAMETER ::  a = 1.0_sp            !< constant
    REAL(sp), PARAMETER ::  b = 0.66666666666_sp  !< constant
    REAL(sp), PARAMETER ::  c = 5.0_sp            !< constant
    REAL(sp), PARAMETER ::  d = 0.35_sp           !< constant
    REAL(sp), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(sp), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_sp )  THEN
       x = SQRT( SQRT( 1.0_sp  - 16.0_sp * zeta ) )
       psi_m_sp = pi * 0.5_sp - 2.0_sp * ATAN( x ) + LOG( ( 1.0_sp + x )**2                        &
                  * ( 1.0_sp + x**2 ) * 0.125_sp )
    ELSE

       psi_m_sp = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta - bc_d_d
!
!--    Old version for stable conditions (only valid for z/L < 0.5) psi_m = - 5.0_wp * zeta

    ENDIF

 END FUNCTION psi_m_sp


!--------------------------------------------------------------------------------------------------!
! Description:
!------------
!> Integrated stability function for heat and moisture. 64-bit version.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION psi_h_dp( zeta )
    !$ACC ROUTINE SEQ

    USE kinds

    IMPLICIT NONE

    REAL(dp), INTENT(IN) ::  zeta   !< stability parameter z/L

    REAL(dp) ::  psi_h_dp  !< integrated similarity function result
    REAL(dp) ::  x         !< dummy variable

    REAL(dp), PARAMETER ::  a = 1.0_dp            !< constant
    REAL(dp), PARAMETER ::  b = 0.66666666666_dp  !< constant
    REAL(dp), PARAMETER ::  c = 5.0_dp            !< constant
    REAL(dp), PARAMETER ::  d = 0.35_dp           !< constant
    REAL(dp), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(dp), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_dp )  THEN
       x = SQRT( 1.0_dp  - 16.0_dp * zeta )
       psi_h_dp = 2.0_dp * LOG( (1.0_dp + x ) / 2.0_dp )
    ELSE
       psi_h_dp = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_dp                              &
                  + 0.66666666666_dp * a * zeta )**1.5_dp - bc_d_d + 1.0_dp
!
!--    Old version for stable conditions (only valid for z/L < 0.5)
!--    psi_h = - 5.0_wp * zeta
    ENDIF

 END FUNCTION psi_h_dp


!--------------------------------------------------------------------------------------------------!
! Description:
!------------
!> Integrated stability function for heat and moisture. 32-bit version.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION psi_h_sp( zeta )
    !$ACC ROUTINE SEQ

    USE kinds

    IMPLICIT NONE

    REAL(sp), INTENT(IN) ::  zeta   !< stability parameter z/L

    REAL(sp) ::  psi_h_sp  !< integrated similarity function result
    REAL(sp) ::  x         !< dummy variable

    REAL(sp), PARAMETER ::  a = 1.0_sp            !< constant
    REAL(sp), PARAMETER ::  b = 0.66666666666_sp  !< constant
    REAL(sp), PARAMETER ::  c = 5.0_sp            !< constant
    REAL(sp), PARAMETER ::  d = 0.35_sp           !< constant
    REAL(sp), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(sp), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_sp )  THEN
       x = SQRT( 1.0_sp  - 16.0_sp * zeta )
       psi_h_sp = 2.0_sp * LOG( (1.0_sp + x ) / 2.0_sp )
    ELSE
       psi_h_sp = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_sp                              &
                  + 0.66666666666_sp * a * zeta )**1.5_sp - bc_d_d + 1.0_sp
!
!--    Old version for stable conditions (only valid for z/L < 0.5)
!--    psi_h = - 5.0_wp * zeta
    ENDIF

 END FUNCTION psi_h_sp


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates stability function for momentum
!>
!> @author Hauke Wurps
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION phi_m( zeta )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    REAL(wp), INTENT(IN) ::  zeta   !< stability parameter z/L

    REAL(wp) ::  phi_m  !< value of the function

    REAL(wp), PARAMETER ::  a = 16.0_wp  !< constant
    REAL(wp), PARAMETER ::  c = 5.0_wp   !< constant

    IF ( zeta < 0.0_wp )  THEN
       phi_m = 1.0_wp / SQRT( SQRT( 1.0_wp - a * zeta ) )
    ELSE
       phi_m = 1.0_wp + c * zeta
    ENDIF

 END FUNCTION phi_m

 END MODULE surface_layer_fluxes_mod
