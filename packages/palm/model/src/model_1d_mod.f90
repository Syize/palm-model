!> @file model_1d_mod.f90
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
! Description:
! ------------
!> 1D-model to initialize the 3D-arrays.
!> The temperature profile is set as steady and a corresponding steady solution
!> of the wind profile is being computed.
!> All subroutines required can be found within this file.
!>
!> @todo harmonize code with new surface_layer_fluxes module
!> @bug 1D model crashes when using small grid spacings in the order of 1 m
!--------------------------------------------------------------------------------------------------!
 MODULE model_1d_mod

    USE arrays_3d,                                                                                 &
        ONLY:  dd2zu,                                                                              &
               ddzu,                                                                               &
               ddzw,                                                                               &
               dzu,                                                                                &
               dzw,                                                                                &
               pt_init,                                                                            &
               q_init,                                                                             &
               ug,                                                                                 &
               u_init,                                                                             &
               vg,                                                                                 &
               v_init,                                                                             &
               zu

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  g,                                                                                  &
               kappa,                                                                              &
               pi

    USE control_parameters,                                                                        &
        ONLY:  constant_diffusion,                                                                 &
               constant_flux_layer,                                                                &
               dissipation_1d,                                                                     &
               f,                                                                                  &
               humidity,                                                                           &
               ibc_e_b,                                                                            &
               implicit_diffusion_1d,                                                              &
               implicit_timestep_factor,                                                           &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               km_constant,                                                                        &
               message_string,                                                                     &
               mixing_length_1d,                                                                   &
               prandtl_number,                                                                     &
               roughness_length,                                                                   &
               run_description_header,                                                             &
               simulated_time_chr,                                                                 &
               timestep_scheme,                                                                    &
               tsc,                                                                                &
               z0h_factor

    USE indices,                                                                                   &
        ONLY:  nzb,                                                                                &
               nzb_diff,                                                                           &
               nzt

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  myid

    IMPLICIT NONE

    INTEGER(iwp) ::  current_timestep_number_1d = 0  !< current timestep number (1d-model)
    INTEGER(iwp) ::  damp_level_ind_1d               !< lower grid index of damping layer (1d-model)

    LOGICAL ::  run_control_header_1d = .FALSE.  !< flag for output of run control header (1d-model)
    LOGICAL ::  stop_dt_1d = .FALSE.             !< termination flag, used in case of too small timestep (1d-model)

    REAL(wp) ::  alpha_buoyancy                !< model constant according to Koblitz (2013)
    REAL(wp) ::  c_0 = 0.416179145_wp          !< = 0.03^0.25; model constant according to Koblitz (2013)
    REAL(wp) ::  c_1 = 1.52_wp                 !< model constant according to Koblitz (2013)
    REAL(wp) ::  c_2 = 1.83_wp                 !< model constant according to Koblitz (2013)
    REAL(wp) ::  c_3                           !< model constant
    REAL(wp) ::  c_mu                          !< model constant
    REAL(wp) ::  damp_level_1d = -1.0_wp       !< namelist parameter
    REAL(wp) ::  dt_1d = 60.0_wp               !< dynamic timestep (1d-model)
    REAL(wp) ::  dt_max_1d = 300.0_wp          !< timestep limit (1d-model)
    REAL(wp) ::  dt_pr_1d = 9999999.9_wp       !< namelist parameter
    REAL(wp) ::  dt_run_control_1d = 60.0_wp   !< namelist parameter
    REAL(wp) ::  end_time_1d = 864000.0_wp     !< namelist parameter
    REAL(wp) ::  lambda                        !< maximum mixing length
    REAL(wp) ::  qs1d                          !< characteristic humidity scale (1d-model)
    REAL(wp) ::  simulated_time_1d = 0.0_wp    !< updated simulated time (1d-model)
    REAL(wp) ::  sig_diss = 2.95_wp            !< model constant according to Koblitz (2013)
    REAL(wp) ::  sig_e = 2.95_wp               !< model constant according to Koblitz (2013)
    REAL(wp) ::  time_pr_1d = 0.0_wp           !< updated simulated time for profile output (1d-model)
    REAL(wp) ::  time_run_control_1d = 0.0_wp  !< updated simulated time for run-control output (1d-model)
    REAL(wp) ::  ts1d                          !< characteristic temperature scale (1d-model)
    REAL(wp) ::  us1d                          !< friction velocity (1d-model)
    REAL(wp) ::  usws1d                        !< u-component of the momentum flux (1d-model)
    REAL(wp) ::  vsws1d                        !< v-component of the momentum flux (1d-model)
    REAL(wp) ::  z01d                          !< roughness length for momentum (1d-model)
    REAL(wp) ::  z0h1d                         !< roughness length for scalars (1d-model)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  column        !< coefficients for Stone-algorithm
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  diss1d        !< tke dissipation rate (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  diss1d_p      !< prognostic value of tke dissipation rate (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e1d           !< tke (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e1d_p         !< prognostic value of tke (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  kh1d          !< turbulent diffusion coefficient for heat (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  km1d          !< turbulent diffusion coefficient for momentum (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l1d           !< mixing length for turbulent diffusion coefficients (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l1d_init      !< initial mixing length (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l1d_diss      !< mixing length for dissipation (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ri1d          !< gradient Richardson number (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_diss       !< tendency of diss, except diffusion (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_diss_diff  !< diffusion-tendency of diss (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_dissm      !< weighted tendency of diss for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_e          !< tendency of e, except diffusion (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_e_diff     !< diffusion-tendency of e (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_em         !< weighted tendency of e for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_u          !< tendency of u, except diffusion (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_u_diff     !< diffusion-tendency of u (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_um         !< weighted tendency of u for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_v          !< tendency of v, except diffusion (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_v_diff     !< diffusion-tendency of v (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_vm         !< weighted tendency of v for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u1d           !< u-velocity component (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u1d_p         !< prognostic value of u-velocity component (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v1d           !< v-velocity component (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v1d_p         !< prognostic value of v-velocity component (1d-model)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tri  !< matrix-coefficients for Stone-algorithm


!
!-- Initialize 1D model
    INTERFACE init_1d_model
       MODULE PROCEDURE init_1d_model
    END INTERFACE init_1d_model

!
!-- Print profiles
    INTERFACE print_1d_model
       MODULE PROCEDURE print_1d_model
    END INTERFACE print_1d_model

!
!-- Print run control information
    INTERFACE run_control_1d
       MODULE PROCEDURE run_control_1d
    END INTERFACE run_control_1d

!
!-- Main procedure
    INTERFACE time_integration_1d
       MODULE PROCEDURE time_integration_1d
    END INTERFACE time_integration_1d

!
!-- Calculate time step
    INTERFACE timestep_1d
       MODULE PROCEDURE timestep_1d
    END INTERFACE timestep_1d

    SAVE

    PRIVATE
!
!-- Public interfaces
    PUBLIC  init_1d_model

!
!-- Public variables
    PUBLIC  damp_level_1d, damp_level_ind_1d, diss1d, dt_pr_1d, dt_run_control_1d, e1d,            &
            end_time_1d, kh1d, km1d, l1d, ri1d, u1d, us1d, usws1d, v1d, vsws1d


    CONTAINS

 SUBROUTINE init_1d_model

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    IMPLICIT NONE

    CHARACTER(LEN=10) ::  time_to_string  !< function to transform time from real to character string

    INTEGER(iwp) ::  k  !< loop index


!
!-- Allocate required 1D-arrays.
    ALLOCATE( column(nzb+1:nzt+1), diss1d(nzb:nzt+1), diss1d_p(nzb:nzt+1),                         &
              e1d(nzb:nzt+1), e1d_p(nzb:nzt+1), kh1d(nzb:nzt+1),                                   &
              km1d(nzb:nzt+1), l1d(nzb:nzt+1), l1d_init(nzb:nzt+1),                                &
              l1d_diss(nzb:nzt+1), ri1d(nzb:nzt+1), te_diss(nzb:nzt+1), te_diss_diff(nzb:nzt+1),   &
              te_dissm(nzb:nzt+1), te_e(nzb:nzt+1), te_e_diff(nzb:nzt+1),                          &
              te_em(nzb:nzt+1), te_u(nzb:nzt+1), te_u_diff(nzb:nzt+1), te_um(nzb:nzt+1),           &
              te_v(nzb:nzt+1), te_v_diff(nzb:nzt+1), te_vm(nzb:nzt+1), u1d(nzb:nzt+1),             &
              u1d_p(nzb:nzt+1),  v1d(nzb:nzt+1), v1d_p(nzb:nzt+1) )
!
!-- Allocate the tridiagonal matrix.
    ALLOCATE( tri(nzb+1:nzt+1,-1:1) )

!
!-- Initialize arrays.
    IF ( constant_diffusion )  THEN
       km1d = km_constant
       kh1d = km_constant / prandtl_number
    ELSE
       diss1d = 0.0_wp; diss1d_p = 0.0_wp
       e1d = 0.0_wp; e1d_p = 0.0_wp

       kh1d = 0.0_wp; km1d = 0.0_wp
       ri1d = 0.0_wp
!
!--    Compute the mixing length.
       l1d_init(nzb) = 0.0_wp

       IF ( TRIM( mixing_length_1d ) == 'blackadar' )  THEN
!
!--       Blackadar mixing length.
          IF ( f /= 0.0_wp )  THEN
             lambda = 2.7E-4_wp * SQRT( ug(nzt+1)**2 + vg(nzt+1)**2 ) / ABS( f ) + 1E-10_wp
          ELSE
             lambda = 30.0_wp
          ENDIF

          DO  k = nzb+1, nzt+1
             l1d_init(k) = kappa * zu(k) / ( 1.0_wp + kappa * zu(k) / lambda )
          ENDDO

       ELSEIF ( TRIM( mixing_length_1d ) == 'as_in_3d_model' )  THEN
!
!--       Use the same mixing length as in 3D model (LES-mode).
!--       This option has been implementing for testing purposes in order to check, if the
!--       3d-model in LES-mode behaves in the same way as the 1d-model. There is no physical
!--       application of this option.
          DO  k = nzb+1, nzt
             l1d_init(k)  = ( dx * dy * dzw(k) )**0.33333333333333_wp
          ENDDO
          l1d_init(nzt+1) = l1d_init(nzt)

       ENDIF
    ENDIF
    l1d      = l1d_init
    l1d_diss = l1d_init
    u1d      = u_init
    u1d_p    = u_init
    v1d      = v_init
    v1d_p    = v_init

!
!-- Set initial horizontal velocities at the lowest grid levels to a very small value in order to
!-- avoid too small time steps caused by the diffusion limit in the initial phase of a run (at k=1,
!-- dz/2 occurs in the limiting formula!).
    u1d(0:1)   = 0.1_wp
    u1d_p(0:1) = 0.1_wp
    v1d(0:1)   = 0.1_wp
    v1d_p(0:1) = 0.1_wp

!
!-- For u*, theta* and the momentum fluxes plausible values are set.
    IF ( constant_flux_layer )  THEN
!
!--    Without initial friction the flow would not change.
       us1d = 0.1_wp
    ELSE
       diss1d(nzb+1) = 0.001_wp
       e1d(nzb+1)  = 1.0_wp
       km1d(nzb+1) = 1.0_wp
       us1d = 0.0_wp
    ENDIF
    ts1d = 0.0_wp
    usws1d = 0.0_wp
    vsws1d = 0.0_wp
    z01d  = roughness_length
    z0h1d = z0h_factor * z01d
    IF ( humidity )  qs1d = 0.0_wp

!
!-- Tendencies must be preset in order to avoid runtime errors.
    te_diss      = 0.0_wp
    te_diss_diff = 0.0_wp
    te_dissm     = 0.0_wp
    te_e         = 0.0_wp
    te_e_diff    = 0.0_wp
    te_em        = 0.0_wp
    te_um        = 0.0_wp
    te_vm        = 0.0_wp

!
!-- Set model constants.
    IF ( dissipation_1d == 'as_in_3d_model' )  c_0 = 0.1_wp
    c_mu = c_0**4

!
!-- Set start time in hh:mm:ss - format.
    simulated_time_chr = time_to_string( simulated_time_1d )

!
!-- Integrate the 1D-model equations using the Runge-Kutta scheme.
    CALL time_integration_1d

 END SUBROUTINE init_1d_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Runge-Kutta time differencing scheme for the 1D-model, with optional implicit diffusion.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE time_integration_1d

    IMPLICIT NONE

    CHARACTER(LEN=10) ::  time_to_string !< function to transform time from real to character string

    INTEGER(iwp) ::  k !< loop index

    REAL(wp) ::  a         !< auxiliary variable
    REAL(wp) ::  b         !< auxiliary variable
    REAL(wp) ::  dpt_dz    !< vertical temperature gradient
    REAL(wp) ::  flux      !< vertical temperature gradient
    REAL(wp) ::  kmzm      !< Km(z-dz/2)
    REAL(wp) ::  kmzp      !< Km(z+dz/2)
    REAL(wp) ::  l_stable  !< mixing length for stable case
    REAL(wp) ::  pt_0      !< reference temperature
    REAL(wp) ::  uv_total  !< horizontal wind speed


!
!-- Set the time step at beginning.
    dt_1d = 0.01_wp

!
!-- Start of time loop.
    DO  WHILE ( simulated_time_1d < end_time_1d  .AND.  .NOT. stop_dt_1d )

!
!--    In case of implicit scheme for diffusion, calculate the diffusion tendencies first.
!--    For the Runge-Kutta scheme, they are added in the last intermediate timestep.
       IF ( implicit_diffusion_1d )  THEN
!
!--       Diffusion for u-component of velocity.
!--       First set matrix elements for the linear equation system.
!--       The diagonal elements are later used for the v-component, too.
          DO  k = nzb+1, nzt

             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )

             tri(k,1)  = -kmzp * dt_1d * ddzu(k+1) * ddzw(k) * 0.5_wp
             tri(k,-1) = -kmzm * dt_1d * ddzu(k)   * ddzw(k) * 0.5_wp

             IF ( k == 1 )  THEN
                tri(k,0)  = 0.5_wp - tri(k,1)
                column(k) = 0.5_wp * ( u1d(k) + dt_1d * ddzw(k) *                                  &
                                     ( kmzp * ddzu(k+1) * ( u1d(k+1) - u1d(k) ) + 2.0_wp * usws1d) )
             ELSE
                tri(k,0)  = 1.0_wp - tri(k,-1) - tri(k,1)
                column(k) = -tri(k,1) * u1d(k+1) - tri(k,-1) * u1d(k-1) +                        &
                            ( 2.0_wp - tri(k,0) ) * u1d(k)
             ENDIF

          ENDDO

          tri(nzt+1,-1) = 0.0_wp
          tri(nzt+1,0)  = 1.0_wp
          column(nzt+1)  = u1d(nzt+1)

          te_u_diff(nzb) = 0.0_wp
!
!--       Calculate new value via Stone algorithm.
          CALL stone_algorithm( te_u_diff )
!
!--       Re-calculate tendency from new value, to be used later in the prognostic equation.
          te_u_diff(:) = ( te_u_diff(:) - u1d(:) ) / dt_1d
!
!--       Diffusion for v-component of velocity.
!--       First set matrix elements for the linear equation system.
          DO  k = nzb+1, nzt

             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )

             IF ( k == 1 )  THEN
                column(k) = 0.5_wp * ( v1d(k) + dt_1d * ddzw(k) *                                  &
                                     ( kmzp * ddzu(k+1) * ( v1d(k+1) - v1d(k) ) + 2.0_wp * vsws1d) )
             ELSE
                column(k) = -tri(k,1) * v1d(k+1) - tri(k,-1) * v1d(k-1) +                          &
                            ( 2.0_wp - tri(k,0) ) * v1d(k)
             ENDIF

          ENDDO

          column(nzt+1)  = v1d(nzt+1)

          te_v_diff(nzb) = 0.0_wp
!
!--       Calculate new value via Stone algorithm.
          CALL stone_algorithm( te_v_diff )
!
!--       Re-calculate tendency from new value, to be used later in the prognostic equation.
          te_v_diff(:) = ( te_v_diff(:) - v1d(:) ) / dt_1d


!
!--       Diffusion for the TKE.
!--       First set matrix elements for the linear equation system.
          DO  k = nzb+1, nzt

             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )

             tri(k,1)  = tri(k,1)  / sig_e
             tri(k,-1) = tri(k,-1) / sig_e

             IF ( k == 1 )  THEN
                tri(k,0)  = 1.0_wp - tri(k,1)
                column(k) = -tri(k,1) * e1d(k+1) + ( 1.0_wp + tri(k,1) ) * e1d(k)
             ELSE
                tri(k,0)  = 1.0_wp - tri(k,1) - tri(k,-1)
                column(k) = -tri(k,1) * e1d(k+1) - tri(k,-1) * e1d(k-1) +                          &
                            ( 2.0_wp - tri(k,0) ) * e1d(k)
             ENDIF

          ENDDO

          tri(nzt+1,-1) = 0.0_wp
          tri(nzt+1,0)  = 1.0_wp
          column(nzt+1) = e1d(nzt+1)
!
!--       Calculate new value via Stone algorithm.
          CALL stone_algorithm( te_e_diff )
!
!--       Re-calculate tendency from new value, to be used later in the prognostic equation.
          te_e_diff(:) = ( te_e_diff(:) - e1d(:) ) / dt_1d
!
!--       Diffusion for the dissipation.
          IF ( dissipation_1d == 'prognostic' )  THEN
!
!--          First set matrix elements for the linear equation system. They are based on the
!--          matrix elements used for TKE, which have been calculated above.
             DO  k = nzb+1, nzt

                tri(k,1)  = tri(k,1)  * sig_e / sig_diss
                tri(k,-1) = tri(k,-1) * sig_e / sig_diss

                IF ( k == 1 )  THEN
!
!--                Here the equation is different than the case in u, v.
!--                Since Neumann boundary condition is used here, i.e. diss(nzb)=diss(nzb*1),
!--                diss(nzb) should not be included in the prognose, which may result in a value
!--                differeing from diss(nzb+1).
!--                Instead, diss(nzb)=diss(nzb+1) is substituted into the first equation.
                   tri(k,0)  = 1.0_wp - tri(k,1)
                   column(k) = -tri(k,1) * diss1d(k+1) + ( 2.0_wp - tri(k,0) ) * diss1d(k)
                ELSE
                   tri(k,0)  = ( tri(k,0) - 1.0_wp ) * sig_e / sig_diss + 1.0_wp
                   column(k) = -tri(k,1) * diss1d(k+1) - tri(k,-1) * diss1d(k-1) +                 &
                               ( 2.0_wp - tri(k,0) ) * diss1d(k)
                ENDIF

             ENDDO

             tri(nzt+1,-1) = 0.0_wp
             tri(nzt+1,0)  = 1.0_wp
             column(nzt+1)  = diss1d(nzt+1)
!
!--          Calculate new value via Stone algorithm.
             CALL stone_algorithm( te_diss_diff )
!
!--          Re-calculate tendency from new value, to be used later in the prognostic equation.
             te_diss_diff(:) = ( te_diss_diff(:) - diss1d(:) ) / dt_1d

          ENDIF

       ENDIF ! implicit diffusion

!
!--    Depending on the timestep scheme, carry out one or more intermediate timesteps.
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1
          CALL timestep_scheme_steering
!
!--       Compute all tendency terms. If a constant-flux layer is simulated, k starts at nzb+2.
          DO  k = nzb_diff, nzt

             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
!
!--          u-component.
             te_u(k) =  f * ( v1d(k) - vg(k) )
!
!--          Explicit diffusion for u-component.
             IF ( .NOT. implicit_diffusion_1d )  THEN
                te_u_diff(k) = ( kmzp * ( u1d(k+1) - u1d(k) ) * ddzu(k+1)                          &
                               - kmzm * ( u1d(k) - u1d(k-1) ) * ddzu(k)                            &
                               ) * ddzw(k)
             ENDIF
!
!--          v-component.
             te_v(k) = -f * ( u1d(k) - ug(k) )

!
!--          Explicit diffusion for v-component.
             IF ( .NOT. implicit_diffusion_1d )  THEN
                te_v_diff(k) = ( kmzp * ( v1d(k+1) - v1d(k) ) * ddzu(k+1)                          &
                               - kmzm * ( v1d(k) - v1d(k-1) ) * ddzu(k)                            &
                               ) * ddzw(k)
             ENDIF

          ENDDO
          IF ( .NOT. constant_diffusion )  THEN
             DO  k = nzb_diff, nzt
!
!--             TKE and dissipation rate.
                kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
                kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(k)
                   flux =  ( pt_init(k+1)-pt_init(k-1) ) * dd2zu(k)
                ELSE
                   pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                   flux = ( ( pt_init(k+1) - pt_init(k-1) ) +                                      &
                            0.61_wp * ( pt_init(k+1) * q_init(k+1) -                               &
                                        pt_init(k-1) * q_init(k-1)   )                             &
                          ) * dd2zu(k)
                ENDIF

!
!--             Calculate dissipation rate if no prognostic equation is used for dissipation rate.
                IF ( dissipation_1d == 'detering' )  THEN
                   diss1d(k) = c_0**3 * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
                ELSEIF ( dissipation_1d == 'as_in_3d_model' )  THEN
                   diss1d(k) = ( 0.19_wp + 0.74_wp * l1d_diss(k) / l1d_init(k) )                   &
                               * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
                ENDIF
!
!--             TKE
                te_e(k) = km1d(k) * ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2                    &
                                    + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2                    &
                                    )                                                              &
                                    - g / pt_0 * kh1d(k) * flux                                    &
                                    - diss1d(k)
!
!--             Explicit diffusion for TKE.
                IF ( .NOT. implicit_diffusion_1d )  THEN
                   te_e_diff(k) = ( kmzp * ( e1d(k+1) - e1d(k) ) * ddzu(k+1)                       &
                                  - kmzm * ( e1d(k) - e1d(k-1) ) * ddzu(k)                         &
                                  ) * ddzw(k) / sig_e

                ENDIF

                IF ( dissipation_1d == 'prognostic' )  THEN
!
!--                Dissipation rate.
                   IF ( ri1d(k) >= 0.0_wp )  THEN
                      alpha_buoyancy = 1.0_wp - l1d(k) / lambda
                   ELSE
                      alpha_buoyancy = 1.0_wp - ( 1.0_wp + ( c_2 - 1.0_wp )                        &
                                                         / ( c_2 - c_1    ) )                      &
                                              * l1d(k) / lambda
                   ENDIF
                   c_3 = ( c_1 - c_2 ) * alpha_buoyancy + 1.0_wp
                   te_diss(k) = ( km1d(k) *                                                        &
                                  ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2                      &
                                  + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2                      &
                                  ) * ( c_1 + (c_2 - c_1) * l1d(k) / lambda )                      &
                                  - g / pt_0 * kh1d(k) * flux * c_3                                &
                                  - c_2 * diss1d(k)                                                &
                                ) * diss1d(k) / ( e1d(k) + 1.0E-20_wp )

!
!--                Explicit diffusion for dissipation rate.
                   IF ( .NOT. implicit_diffusion_1d )  THEN
                      te_diss_diff(k) = ( kmzp * ( diss1d(k+1) - diss1d(k) ) * ddzu(k+1)           &
                                        - kmzm * ( diss1d(k) - diss1d(k-1) ) * ddzu(k)             &
                                        ) * ddzw(k) / sig_diss
                   ENDIF

                ENDIF

             ENDDO
          ENDIF
!
!--       Tendency terms at the top of the constant-flux layer.
!--       Finite differences of the momentum fluxes are computed using half the normal grid length
!--       (2.0*ddzw(k)) for the sake of enhanced accuracy.
          IF ( constant_flux_layer )  THEN

             k = nzb+1
             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
             IF ( .NOT. humidity )  THEN
                pt_0 = pt_init(k)
                flux =  ( pt_init(k+1)-pt_init(k-1) ) * dd2zu(k)
             ELSE
                pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                flux = ( ( pt_init(k+1) - pt_init(k-1) ) +                                         &
                         0.61_wp * ( pt_init(k+1) * q_init(k+1) -                                  &
                                     pt_init(k-1) * q_init(k-1)   )                                &
                       ) * dd2zu(k)
             ENDIF

!
!--          Calculate dissipation rate if no prognostic equation is used for dissipation rate.
             IF ( dissipation_1d == 'detering' )  THEN
                diss1d(k) = c_0**3 * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
             ELSEIF ( dissipation_1d == 'as_in_3d_model' )  THEN
                diss1d(k) = ( 0.19_wp + 0.74_wp * l1d_diss(k) / l1d_init(k) )                      &
                            * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
             ENDIF

!
!--          u-component.
             te_u(k) = f * ( v1d(k) - vg(k) )
!
!--          Explicit diffusion for u-component.
             IF ( .NOT. implicit_diffusion_1d )  THEN
                te_u_diff(k) = ( kmzp * ( u1d(k+1) - u1d(k) ) * ddzu(k+1) + usws1d                 &
                               ) * 2.0_wp * ddzw(k)
             ENDIF
!
!--          v-component.
             te_v(k) = -f * ( u1d(k) - ug(k) )
!
!--          Explicit diffusion for v-component.
             IF ( .NOT. implicit_diffusion_1d )  THEN
                te_v_diff(k) = ( kmzp * ( v1d(k+1) - v1d(k) ) * ddzu(k+1) + vsws1d                 &
                               ) * 2.0_wp * ddzw(k)
             ENDIF
!
!--          TKE.
             IF ( .NOT. dissipation_1d == 'prognostic' )  THEN
                !> @query why integrate over 2dz
                !>   Why is it allowed to integrate over two delta-z for e
                !>   while for u and v it is not?
                !>   2018-04-23, gronemeier
                te_e(k) = km1d(k) * ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2                    &
                                    + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2                    &
                                    )                                                              &
                                    - g / pt_0 * kh1d(k) * flux                                    &
                                    - diss1d(k)
!
!--             Explicit diffusion for TKE.
                IF ( .NOT. implicit_diffusion_1d )  THEN
                   te_e_diff(k) = ( kmzp * ( e1d(k+1) - e1d(k) ) * ddzu(k+1)                       &
                                  - kmzm * ( e1d(k) - e1d(k-1) ) * ddzu(k)                         &
                                  ) * ddzw(k) / sig_e

                ENDIF
             ENDIF

          ENDIF
!
!--       Prognostic equations for all 1D variables.
          DO  k = nzb+1, nzt

             u1d_p(k) = u1d(k) + dt_1d * ( tsc(2) * te_u(k) + tsc(4) * te_u_diff(k) +              &
                                           tsc(3) * te_um(k) )
             v1d_p(k) = v1d(k) + dt_1d * ( tsc(2) * te_v(k) + tsc(4) * te_v_diff(k) +              &
                                           tsc(3) * te_vm(k) )

          ENDDO
          IF ( .NOT. constant_diffusion )  THEN

             DO  k = nzb+1, nzt
                e1d_p(k) = e1d(k) + dt_1d * ( tsc(2) * te_e(k) + tsc(4) * te_e_diff(k) +           &
                                              tsc(3) * te_em(k) )
             ENDDO
!
!--          Eliminate negative TKE values, which can result from the integration due to numerical
!--          inaccuracies. In such cases the TKE value is reduced to 10 percent of its old value.
             WHERE ( e1d_p < 0.0_wp )  e1d_p = 0.1_wp * e1d

             IF ( dissipation_1d == 'prognostic' )  THEN
                DO  k = nzb+1, nzt
                   diss1d_p(k) = diss1d(k) + dt_1d * ( tsc(2) * te_diss(k) +                       &
                                                       tsc(4) * te_diss_diff(k) +                  &
                                                       tsc(3) * te_dissm(k) )
                ENDDO
                WHERE ( diss1d_p < 0.0_wp )  diss1d_p = 0.1_wp * diss1d
             ENDIF
          ENDIF
!
!--       Calculate tendencies for the next Runge-Kutta step.
          IF ( timestep_scheme(1:5) == 'runge' ) THEN

             IF ( intermediate_timestep_count == 1 )  THEN

                DO  k = nzb+1, nzt
                   te_um(k) = te_u(k) + tsc(1) * te_u_diff(k)
                   te_vm(k) = te_v(k) + tsc(1) * te_v_diff(k)
                ENDDO

                IF ( .NOT. constant_diffusion )  THEN
                   DO  k = nzb+1, nzt
                      te_em(k) = te_e(k) + tsc(1) * te_e_diff(k)
                   ENDDO
                   IF ( dissipation_1d == 'prognostic' )  THEN
                      DO  k = nzb+1, nzt
                         te_dissm(k) = te_diss(k) + tsc(1) * te_diss_diff(k)
                      ENDDO
                   ENDIF
                ENDIF

             ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN

                DO  k = nzb+1, nzt
                   te_um(k) = -9.5625_wp * ( te_u(k) + tsc(1) * te_u_diff(k) ) +                   &
                               5.3125_wp * te_um(k)
                   te_vm(k) = -9.5625_wp * ( te_v(k) + tsc(1) * te_v_diff(k) ) +                   &
                               5.3125_wp * te_vm(k)
                ENDDO

                IF ( .NOT. constant_diffusion )  THEN
                   DO  k = nzb+1, nzt
                      te_em(k) = -9.5625_wp * ( te_e(k) + tsc(1) * te_e_diff(k) ) +                &
                                  5.3125_wp * te_em(k)
                   ENDDO
                   IF ( dissipation_1d == 'prognostic' )  THEN
                      DO  k  = nzb+1, nzt
                         te_dissm(k) = -9.5625_wp * ( te_diss(k) + tsc(1) * te_diss_diff(k) ) +    &
                                        5.3125_wp * te_dissm(k)
                      ENDDO
                   ENDIF
                ENDIF

             ENDIF

          ENDIF
!
!--       Boundary conditions for the prognostic variables.
!--       At the top boundary (nzt+1) u, v, e, and diss keep their initial values (ug(nzt+1),
!--       vg(nzt+1), 0, 0).
!--       At the bottom boundary, Dirichlet condition is used for u and v (0) and Neumann condition
!--       for e and diss (e(nzb)=e(nzb+1)).
          u1d_p(nzb) = 0.0_wp
          v1d_p(nzb) = 0.0_wp
!
!--       Swap the time levels in preparation for the next time step.
          u1d  = u1d_p
          v1d  = v1d_p
          IF ( .NOT. constant_diffusion )  THEN
             e1d  = e1d_p
             IF ( dissipation_1d == 'prognostic' )  THEN
                diss1d = diss1d_p
             ENDIF
          ENDIF
!
!--       Compute diagnostic diffusion quantities.
          IF ( .NOT. constant_diffusion )  THEN

!
!--          First compute the vertical fluxes in the constant-flux layer.
             IF ( constant_flux_layer )  THEN
!
!--             Compute theta* using Ri numbers of the previous time step.
                IF ( ri1d(nzb+1) >= 0.0_wp )  THEN
!
!--                Stable stratification.
                   ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /                              &
                          ( LOG( zu(nzb+1) / z0h1d ) + 5.0_wp * ri1d(nzb+1) *                      &
                                          ( zu(nzb+1) - z0h1d ) / zu(nzb+1)                        &
                          )
                ELSE
!
!--                Unstable stratification.
                   a = SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) )
                   b = SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) / zu(nzb+1) * z0h1d )

                   ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /                              &
                          LOG( (a-1.0_wp) / (a+1.0_wp) * (b+1.0_wp) / (b-1.0_wp) )
                ENDIF
!
!--             Compute the gradient Richardson numbers,
!--             first at the top of the constant-flux layer using u* of the previous time step
!--             (+1E-30, if u* = 0), then in the remaining area.
!--             There, the Ri numbers of the previous time step are used.
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(nzb+1)
                   flux = ts1d
                ELSE
                   pt_0 = pt_init(nzb+1) * ( 1.0_wp + 0.61_wp * q_init(nzb+1) )
                   flux = ts1d + 0.61_wp * pt_init(k) * qs1d
                ENDIF
                ri1d(nzb+1) = zu(nzb+1) * kappa * g * flux / ( pt_0 * ( us1d**2 + 1E-30_wp ) )

             ENDIF

             DO  k = nzb_diff, nzt
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(k)
                   flux = ( pt_init(k+1) - pt_init(k-1) ) * dd2zu(k)
                ELSE
                   pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                   flux = ( ( pt_init(k+1) - pt_init(k-1) )                                        &
                            + 0.61_wp                                                              &
                            * (   pt_init(k+1) * q_init(k+1)                                       &
                                - pt_init(k-1) * q_init(k-1) )                                     &
                          ) * dd2zu(k)
                ENDIF
                IF ( ri1d(k) >= 0.0_wp )  THEN
                   ri1d(k) = g / pt_0 * flux /                                                     &
                              (  ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2                         &
                               + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2                         &
                               + 1E-30_wp                                                          &
                              )
                ELSE
                   ri1d(k) = g / pt_0 * flux /                                                     &
                              (  ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2                         &
                               + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2                         &
                               + 1E-30_wp                                                          &
                              ) * SQRT( SQRT( 1.0_wp - 16.0_wp * ri1d(k) ) )
                ENDIF
             ENDDO
!
!--          Richardson numbers must remain restricted to a realistic value range. It is exceeded
!--          excessively for very small velocities (u,v --> 0).
             WHERE ( ri1d < -5.0_wp )  ri1d = -5.0_wp
             WHERE ( ri1d > 1.0_wp )  ri1d = 1.0_wp
!
!--          Compute u* from the absolute velocity value.
             IF ( constant_flux_layer )  THEN
                uv_total = SQRT( u1d(nzb+1)**2 + v1d(nzb+1)**2 )

                IF ( ri1d(nzb+1) >= 0.0_wp )  THEN
!
!--                Stable stratification.
                   us1d = kappa * uv_total / ( LOG( zu(nzb+1) / z01d )                             &
                                               + 5.0_wp * ri1d(nzb+1) * ( zu(nzb+1) - z01d )       &
                                                 / zu(nzb+1)                                       &
                                             )
                ELSE
!
!--                Unstable stratification.
                   a = 1.0_wp / SQRT( SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) ) )
                   b = 1.0_wp / SQRT( SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) / zu(nzb+1) * z01d ) )
                   us1d = kappa * uv_total / ( LOG( (1.0_wp+b) / (1.0_wp-b) * (1.0_wp-a) /         &
                                                    (1.0_wp+a) ) +                                 &
                                               2.0_wp * ( ATAN( b ) - ATAN( a ) )                  &
                                             )
                ENDIF
!
!--             Compute the momentum fluxes for the diffusion terms.
                usws1d  = - u1d(nzb+1) / uv_total * us1d**2
                vsws1d  = - v1d(nzb+1) / uv_total * us1d**2
!
!--             Boundary condition for the turbulent kinetic energy and dissipation rate at the top
!--             of the constant-flux layer.
!--             Additional Neumann condition de/dz = 0 at nzb is set to ensure compatibility with
!--             the 3D model.
                IF ( ibc_e_b == 2 )  THEN
                   e1d(nzb+1) = ( us1d / c_0 )**2
                ENDIF
                IF ( dissipation_1d == 'prognostic' )  THEN
                   e1d(nzb+1) = ( us1d / c_0 )**2
                   diss1d(nzb+1) = us1d**3 / ( kappa * zu(nzb+1) )
                   diss1d(nzb) = diss1d(nzb+1)
                ENDIF
                e1d(nzb) = e1d(nzb+1)

                IF ( humidity ) THEN
!
!--                Compute q*.
                   IF ( ri1d(nzb+1) >= 0.0_wp )  THEN
!
!--                   Stable stratification.
                      qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /                             &
                             ( LOG( zu(nzb+1) / z0h1d ) + 5.0_wp * ri1d(nzb+1) *                   &
                                             ( zu(nzb+1) - z0h1d ) / zu(nzb+1)                     &
                             )
                   ELSE
!
!--                   Unstable stratification.
                      a = SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) )
                      b = SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) / zu(nzb+1) * z0h1d )
                      qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /                             &
                             LOG( (a-1.0_wp) / (a+1.0_wp) * (b+1.0_wp) / (b-1.0_wp) )
                   ENDIF
                ELSE
                   qs1d = 0.0_wp
                ENDIF

             ENDIF   !  constant_flux_layer
!
!--          Compute the diabatic mixing length. The unstable stratification must not be considered
!--          for l1d (km1d) as it is already considered in the dissipation of TKE via l1d_diss.
!--          Otherwise, km1d would be too large.
             IF ( dissipation_1d /= 'prognostic' )  THEN
                IF ( mixing_length_1d == 'blackadar' )  THEN
                   DO  k = nzb+1, nzt
                      IF ( ri1d(k) >= 0.0_wp )  THEN
                         l1d(k) = l1d_init(k) / ( 1.0_wp + 5.0_wp * ri1d(k) )
                         l1d_diss(k) = l1d(k)
                      ELSE
                         l1d(k) = l1d_init(k)
                         l1d_diss(k) = l1d_init(k) * SQRT( 1.0_wp - 16.0_wp * ri1d(k) )
                      ENDIF
                   ENDDO
                ELSEIF ( mixing_length_1d == 'as_in_3d_model' )  THEN
                   DO  k = nzb+1, nzt
                      dpt_dz = ( pt_init(k+1) - pt_init(k-1) ) * dd2zu(k)
                      IF ( dpt_dz > 0.0_wp )  THEN
                         l_stable = 0.76_wp * SQRT( e1d(k) )                                       &
                                    / SQRT( g / pt_init(k) * dpt_dz ) + 1E-5_wp
                      ELSE
                         l_stable = l1d_init(k)
                      ENDIF
                      l1d(k) = MIN( l1d_init(k), l_stable )
                      l1d_diss(k) = l1d(k)
                   ENDDO
                ENDIF
             ELSE
                DO  k = nzb+1, nzt
                   l1d(k) = c_0**3 * e1d(k) * SQRT( e1d(k) ) / ( diss1d(k) + 1.0E-30_wp )
                ENDDO
             ENDIF
!
!--          Compute the diffusion coefficients for momentum via the corresponding Prandtl-layer
!--          relationship and according to Prandtl-Kolmogorov, respectively.
             IF ( constant_flux_layer )  THEN
                IF ( ri1d(nzb+1) >= 0.0_wp )  THEN
                   km1d(nzb+1) = us1d * kappa * zu(nzb+1) /                                        &
                                 ( 1.0_wp + 5.0_wp * ri1d(nzb+1) )
                ELSE
                   km1d(nzb+1) = us1d * kappa * zu(nzb+1) *                                        &
                                 SQRT( SQRT( 1.0_wp - 16.0_wp * ri1d(nzb+1) ) )
                ENDIF
             ENDIF

             IF ( dissipation_1d == 'prognostic' )  THEN
                DO  k = nzb_diff, nzt
                   km1d(k) = c_mu * e1d(k)**2 / ( diss1d(k) + 1.0E-30_wp )
                ENDDO
             ELSE
                DO  k = nzb_diff, nzt
                   km1d(k) = c_0 * SQRT( e1d(k) ) * l1d(k)
                ENDDO
             ENDIF
!
!--          Add damping layer.
             DO  k = damp_level_ind_1d+1, nzt+1
                km1d(k) = 1.1_wp * km1d(k-1)
                km1d(k) = MIN( km1d(k), 10.0_wp )
             ENDDO
!
!--          Compute the diffusion coefficient for heat via the relationship kh = phim / phih * km.
             DO  k = nzb+1, nzt
                IF ( ri1d(k) >= 0.0_wp )  THEN
                   kh1d(k) = km1d(k)
                ELSE
                   kh1d(k) = km1d(k) * SQRT( SQRT( 1.0_wp - 16.0_wp * ri1d(k) ) )
                ENDIF
             ENDDO

          ENDIF   ! .NOT. constant_diffusion

       ENDDO   ! intermediate step loop

!
!--    Increment simulated time and output times.
       current_timestep_number_1d = current_timestep_number_1d + 1
       simulated_time_1d          = simulated_time_1d + dt_1d
       simulated_time_chr         = time_to_string( simulated_time_1d )
       time_pr_1d                 = time_pr_1d          + dt_1d
       time_run_control_1d        = time_run_control_1d + dt_1d
!
!--    Determine and print out quantities for run control.
       IF ( time_run_control_1d >= dt_run_control_1d )  THEN
          CALL run_control_1d
          time_run_control_1d = time_run_control_1d - dt_run_control_1d
       ENDIF
!
!--    Profile output on file.
       IF ( time_pr_1d >= dt_pr_1d )  THEN
          CALL print_1d_model
          time_pr_1d = time_pr_1d - dt_pr_1d
       ENDIF
!
!--    Determine size of next time step.
       CALL timestep_1d

    ENDDO   ! time loop
!
!-- Set intermediate_timestep_count back to zero. This is required e.g. for initial calls of
!-- calc_mean_profile.
    intermediate_timestep_count = 0

 END SUBROUTINE time_integration_1d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute and print out quantities for run control of the 1D model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE run_control_1d

    IMPLICIT NONE

    INTEGER(iwp) ::  k !< loop index

    REAL(wp) ::  alpha     !< angle of wind vector at top of constant-flux layer
    REAL(wp) ::  energy    !< kinetic energy
    REAL(wp) ::  umax      !< maximum of u
    REAL(wp) ::  uv_total  !< horizontal wind speed
    REAL(wp) ::  vmax      !< maximum of v


    IF ( myid == 0 )  THEN
!
!--    If necessary, write header.
       IF ( .NOT. run_control_header_1d )  THEN
          CALL check_open( 15 )
          WRITE ( 15, 100 )
          run_control_header_1d = .TRUE.
       ENDIF

!
!--    Compute control quantities.
       umax = 0.0_wp; vmax = 0.0_wp; energy = 0.0_wp
       DO  k = nzb+1, nzt+1
          umax = MAX( ABS( umax ), ABS( u1d(k) ) )
          vmax = MAX( ABS( vmax ), ABS( v1d(k) ) )
          energy = energy + 0.5_wp * ( u1d(k)**2 + v1d(k)**2 )
       ENDDO
       energy = energy / REAL( nzt - nzb + 1, KIND=wp )

       uv_total = SQRT( u1d(nzb+1)**2 + v1d(nzb+1)**2 )
       IF ( ABS( v1d(nzb+1) ) < 1.0E-5_wp )  THEN
          alpha = ACOS( SIGN( 1.0_wp , u1d(nzb+1) ) )
       ELSE
          alpha = ACOS( u1d(nzb+1) / uv_total )
          IF ( v1d(nzb+1) <= 0.0_wp )  alpha = 2.0_wp * pi - alpha
       ENDIF
       alpha = alpha / ( 2.0_wp * pi ) * 360.0_wp

       WRITE ( 15, 101 )  current_timestep_number_1d, simulated_time_chr, dt_1d, umax, vmax, us1d, &
                          alpha, energy
!
!--    Write buffer contents to disc immediately.
       FLUSH( 15 )

    ENDIF
!
!-- Formats.
100 FORMAT (///'1D run control output:'/                                                           &
               '------------------------------'//                                                  &
            'ITER.   HHH:MM:SS    DT      UMAX   VMAX    U*   ALPHA   ENERG.'/                     &
            '---------------------------------------------------------------')
101 FORMAT (I7,1X,A10,1X,F6.2,2X,F6.2,1X,F6.2,1X,F6.3,2X,F5.1,2X,F7.2)

 END SUBROUTINE run_control_1d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the time step with respect to the diffusion criterion.
!> The implicit time step scheme allows larger values.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE timestep_1d

    IMPLICIT NONE

    INTEGER(iwp) ::  k !< loop index

    REAL(wp) ::  dt_diff !< time step accorind to diffusion criterion
    REAL(wp) ::  dt_old  !< previous time step
    REAL(wp) ::  fac     !< factor of criterion
    REAL(wp) ::  value   !< auxiliary variable


!
!-- Save previous time step.
    dt_old = dt_1d
!
!-- Compute the currently feasible time step according to the diffusion criterion. At nzb+1 the half
!-- grid length is used.
    fac = 0.125
    IF ( implicit_diffusion_1d )  fac = fac * implicit_timestep_factor

    dt_diff = dt_max_1d
    DO  k = nzb+2, nzt
       value   = fac * dzu(k) * dzu(k) / ( km1d(k) + 1E-20_wp )
       dt_diff = MIN( value, dt_diff )
    ENDDO
    value   = fac * zu(nzb+1) * zu(nzb+1) / ( km1d(nzb+1) + 1E-20_wp )
    dt_1d = MIN( value, dt_diff )
!
!-- Limit the new time step to a maximum of 10 times the previous time step.
    dt_1d = MIN( dt_old * 10.0_wp, dt_1d )
!
!-- Set flag when the time step becomes too small
    IF ( dt_1d < ( 1.0E-15_wp * dt_max_1d ) )  THEN
       stop_dt_1d = .TRUE.

       WRITE( message_string, * ) 'simulaton has been stopped since timestep has fallen below ',   &
                                  'the lower limit dt_1d = ', dt_1d
       CALL message( 'timestep_1d', 'PAC0248', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE timestep_1d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> List output of profiles from the 1D-model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE print_1d_model

    IMPLICIT NONE

    INTEGER(iwp) ::  k !< loop parameter

    LOGICAL, SAVE :: write_first = .TRUE. !< flag for writing header


    IF ( myid == 0 )  THEN
!
!--    Open list output file for profiles from the 1D-model.
       CALL check_open( 17 )
!
!--    Write Header.
       IF ( write_first )  THEN
          WRITE ( 17, 100 )  TRIM( run_description_header )
          write_first = .FALSE.
       ENDIF
!
!--    Write the values.
       WRITE ( 17, 104 )  TRIM( simulated_time_chr )
       WRITE ( 17, 101 )
       WRITE ( 17, 102 )
       WRITE ( 17, 101 )
       DO  k = nzt+1, nzb, -1
          WRITE ( 17, 103)  k, zu(k), u1d(k), v1d(k), pt_init(k), e1d(k), ri1d(k), km1d(k),        &
                            kh1d(k), l1d(k), diss1d(k)
       ENDDO
       WRITE ( 17, 101 )
       WRITE ( 17, 102 )
       WRITE ( 17, 101 )
!
!--    Write buffer contents to disc immediately.
       FLUSH( 17 )

    ENDIF
!
!-- Formats.
100 FORMAT ('# ',A/'#',10('-')/'# 1d-model profiles')
104 FORMAT (//'# Time: ',A)
101 FORMAT ('#',111('-'))
102 FORMAT ('#  k     zu      u          v          pt         e          ',   &
            'Ri         Km         Kh         l          diss')
103 FORMAT (1X,I4,1X,F7.1,9(1X,E10.3))

 END SUBROUTINE print_1d_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solves a tridiagonal equation system via the Stone algorithm.
!> See: Stone, H. (1973): An Efficient Parallel Algorithm for the Solution of a Tridiagonal
!> Linear System of Equations, Journal of the Association for Computing Machinery, 20, 27-38.
!> Only grid points (nzb+1:nzt+1) are calculated. Boundary point (nzb) is not modified.
!> The tridiagonal matrix A in the system Ax=column will be factorised into A=LU, where
!> L has a diagonal of elements 1 and a subdiagonal m,
!> U has a diagonal u and superdiagonal same as that of A.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stone_algorithm( x )

    INTEGER(iwp) ::  k !< loop index

    REAL(wp), DIMENSION(nzb+2:nzt+1) ::  m  !< subdiagonal of L
    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  q  !< auxiliary variable
    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  u  !< diagonal of U
    REAL(wp), DIMENSION(nzb:nzt+1)   ::  x  !< variable to which the solution will be output
    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  y  !< y:=Ux, so Ly=column, used in forward substitution


!
!-- Decomposition.
    CALL decompose

    u(nzb+1) = q(nzb+1)
    DO  k = nzb+2, nzt+1
       u(k) = q(k) / q(k-1)
    ENDDO

    m(nzb+2) = tri(nzb+2,-1) / tri(nzb+1,0)
    DO  k = nzb+2, nzt+1
       m(k) = tri(k,-1) / u(k-1)
    ENDDO
!
!-- Forward substitution.
    CALL linear_forward
!
!-- Backward substitution.
    CALL linear_backward

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Part of Stone algorithm. The variable nomenclature follows the Stone paper.
!> Calculates u1            = d1,
!>            u_(i)         = d_(i) - e_(i)*f_(i-1)/u_(i-1)
!>            q_i           = d_i*q_(i-1) - e_i*f_(i-1)*q_(i-2)
!>            (q_i q(i-1))T = {{d_i, -e_i*f_(i-1)},{1,0}} (q_(i-1), q_(i-2))T
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE decompose

       INTEGER(iwp) ::  i  !< loop index
       INTEGER(iwp) ::  k  !< loop index

       REAL(wp), DIMENSION(nzb+1:nzt+1) ::  ef    !< holds the product -e_i*f_{i-1}
       REAL(wp), DIMENSION(nzb:nzt+1)   ::  qim1  !< holds q(i-1)
       REAL(wp), DIMENSION(nzb-1:nzt+1) ::  qim2  !< holds q(i-2)
       REAL(wp), DIMENSION(nzb+1:nzt+1) ::  tmp   !< temporary variable


       ef(nzb+1) = 0.0_wp
       DO  k = nzb+2, nzt+1
          ef(k) = -tri(k,-1) * tri(k-1,1)
       ENDDO

       qim2(:)   = 1.0_wp
       qim1(nzb) = 1.0_wp
       DO  k = nzb+1, nzt+1
          qim1(k) = tri(k,0)
       ENDDO

       q(nzb+1) = tri(1,0)
       DO  k = nzb+2, nzt+1
          q(k) = tri(k,0) * tri(k-1,0) + ef(k)
       ENDDO

       i = 2
       DO WHILE ( i < nzt-nzb+1 )

          DO  k = nzb+i-1, nzt+1
             tmp(k) = qim1(k) * qim1(k-i+1) + ef(k-i+2) * qim2(k) * qim2(k-i)
          ENDDO
!
!--       Note that Fortran evaluates RHS first. Don't replace the array assignment by a loop, since
!--       it would generate a different result because of recurrence!
          qim1(nzb+i:nzt+1) = q(nzb+i:nzt+1) * qim1(nzb+0:nzt+1-i) +                               &
                              ef(nzb+1:nzt-i+2) * qim1(nzb+i:nzt+1) * qim2(nzb-1:nzt-i)

          DO  k = nzb+i-1, nzt+1
             qim2(k) = tmp(k)
          ENDDO

          DO  k = nzb+i+1, nzt+1
             q(k) = tri(k,0) * qim1(k-1) + ef(k) * qim2(k-2)
         ENDDO

          i = i * 2

       ENDDO

    END SUBROUTINE decompose


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Part of Stone algorithm. The variable nomenclature follows the Stone paper.
!> Calculates y1 = a1, y_i = a_i + b_i+y_(i-1)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE linear_forward

       INTEGER(iwp) ::  i  !< loop index
       INTEGER(iwp) ::  k  !< loop index

       REAL(wp), DIMENSION(nzb+1:nzt+1) ::  mable  !< auxiliary variable to avoid negative sign, holds -m


       y(:) = column(:)
!
!--    Value at nzb+1 is arbitrary, as it does not affect the result.
       mable(nzb+1) = 1.0_wp
       DO  k = nzb+2, nzt+1
          mable(k) = -m(k)
       ENDDO
!
!--    Note that Fortran evaluates RHS first. Don't replace the array assignments by loops, since
!--    they would generate a different result because of recurrence!
       i = 1
       DO WHILE ( i < nzt-nzb+1 )
          y(nzb+i+1:nzt+1)     = y(nzb+i+1:nzt+1) + y(nzb+1:nzt+1-i) * mable(nzb+i+1:nzt+1)
          mable(nzb+i+1:nzt+1) = mable(nzb+i+1:nzt+1) * mable(nzb+1:nzt+1-i)

          i = i * 2
       ENDDO

    END SUBROUTINE linear_forward


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Part of Stone algorithm. The variable nomenclature follows the Stone paper.
!> Calculates xn = an, x_i = a_i + b_i*x_(i+1)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE linear_backward

       INTEGER(iwp) ::  i  !< loop index
       INTEGER(iwp) ::  k  !< loop index

       REAL(wp), DIMENSION(nzb+1:nzt+1) ::  mars  !< auxiliary variables so that the algorithm takes the form of linear_forward
       REAL(wp), DIMENSION(nzb+1:nzt+1) ::  xav   !< auxiliary variables so that the algorithm takes the form of linear_forward


       xav(:) = y(:) / u(:)
!
!--    Value at nzt+1 is arbitrary, as it does not affect the result.
       DO  k = nzb+1, nzt
          mars(k) = -tri(k,1) / u(k)
       ENDDO
       mars(nzt+1) = 1.0_wp
!
!--    Note that Fortran evaluates RHS first. Don't replace the array assignments by loops, since
!--    they would generate a different result because of recurrence!
       i = 1
       DO WHILE ( i < nzt-nzb+1 )
          xav(nzb+1:nzt+1-i)  = xav(nzb+1:nzt+1-i) + xav(nzb+i+1:nzt+1) * mars(nzb+1:nzt+1-i)
          mars(nzb+1:nzt+1-i) = mars(nzb+1:nzt+1-i) * mars(nzb+i+1:nzt+1)

          i = i * 2
       ENDDO

       DO  k = nzb+1, nzt+1
          x(k) = xav(k)
       ENDDO

    END SUBROUTINE linear_backward

 END SUBROUTINE stone_algorithm

 END MODULE
