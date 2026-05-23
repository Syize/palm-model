!> @file poismg.f90
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
! Copyright 2025      pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid V- or W-Cycle scheme.
!>
!> This multigrid method was originally developed for PALM by Joerg Uhlenbrock,
!> September 2000 - July 2001. It has been optimised for speed by Klaus Ketelsen in November 2014.
!>
!> @attention Loop unrolling and cache optimization in SOR-Red/Black method still does not give the
!>             expected speedup!
!--------------------------------------------------------------------------------------------------!
 MODULE poismg_mod

#if defined( __parallel )
    USE MPI
#endif

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               message_string,                                                                     &
               nesting_offline

    USE control_parameters,                                                                        &
        ONLY:  enable_openacc

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz,                                                                     &
               exchange_horiz_int

#if defined( __parallel )
    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_rb
#endif

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nzb,                                                                                &
               nzt

    USE kinds

    USE pegrid


    INTEGER, PRIVATE ::  ind_even_odd  !< border index between even and odd k index

    CHARACTER (LEN=1) ::  cycle_mg = 'w'  !< namelist parameter (see documentation)

    INTEGER(iwp) ::  gamma_mg                     !< switch for steering the multigrid cycle: 1: v-cycle, 2: w-cycle
    INTEGER(iwp) ::  gathered_size                !< number of total domain grid points of the grid level which is gathered on
                                                  !< PE0 (multigrid solver)
    INTEGER(iwp) ::  grid_level                   !< current grid level
    INTEGER(iwp) ::  maximum_grid_level = 0       !< number of grid levels that the multigrid solver is using
    INTEGER(iwp) ::  maximum_grid_level_default   !< maximum number of grid levels that the multigrid solver could use if not restricted by max_mg_grid_levels
    INTEGER(iwp) ::  max_mg_grid_levels = 9999    !< namelist parameter
    INTEGER(iwp) ::  mgcycles = 0                 !< number of cycles that the multigrid solver has actually carried out
    INTEGER(iwp) ::  mg_cycles = 4                !< namelist parameter
    INTEGER(iwp) ::  mg_switch_to_pe0_level = -1  !< namelist parameter
    INTEGER(iwp) ::  ngsrb = 2                    !< namelist parameter
    INTEGER(iwp) ::  ngsrb_initial = 100          !< namelist parameter
    INTEGER(iwp) ::  ngsrb_initial_timesteps = 0  !< namelist parameter
    INTEGER(iwp) ::  subdomain_size               !< number of grid points in (3d) subdomain including ghost points

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  even_odd_level    !< stores ind_even_odd for all MG levels
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  grid_level_count  !< internal switch for steering the multigrid v- and w-cycles
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nxl_mg            !< left-most grid index of subdomain on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nxr_mg            !< right-most grid index of subdomain on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nyn_mg            !< north-most grid index of subdomain on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nys_mg            !< south-most grid index of subdomain on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nzt_mg            !< top-most grid index of subdomain on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  poismg_filtered_holes  !< number of filtered holes on respective grid level


    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  mg_loc_ind  !< array to store index bounds of all PEs of that multigrid level
                                                              !< where data is collected to PE0

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ileft                 !< lower loop index for red/black i loop
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  jsouth                !< lower loop index for red/black j loop
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom               !< lower loop index for red/black k loop
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop                  !< upper loop index for red/black k loop

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ileft_for_nyn_recv    !< lower i loop index for ghost point exchange at north boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ileft_for_nyn_send    !< lower i loop index for ghost point exchange at north boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ileft_for_nys_recv    !< lower i loop index for ghost point exchange at south boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ileft_for_nys_send    !< lower i loop index for ghost point exchange at south boundary

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  jsouth_for_nxl_recv   !< lower j loop index for ghost point exchange at left boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  jsouth_for_nxl_send   !< lower j loop index for ghost point exchange at left boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  jsouth_for_nxr_recv   !< lower j loop index for ghost point exchange at right boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  jsouth_for_nxr_send   !< lower j loop index for ghost point exchange at right boundary

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nxl_recv   !< lower k loop index for ghost point exchange at left boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nxl_send   !< lower k loop index for ghost point exchange at left boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nxr_recv   !< lower k loop index for ghost point exchange at right boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nxr_send   !< lower k loop index for ghost point exchange at right boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nyn_recv   !< lower k loop index for ghost point exchange at north boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nyn_send   !< lower k loop index for ghost point exchange at north boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nys_recv   !< lower k loop index for ghost point exchange at south boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  kbottom_for_nys_send   !< lower k loop index for ghost point exchange at south boundary

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nxl_recv      !< upper k loop index for ghost point exchange at left boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nxl_send      !< upper k loop index for ghost point exchange at left boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nxr_recv      !< upper k loop index for ghost point exchange at right boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nxr_send      !< upper k loop index for ghost point exchange at right boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nyn_recv      !< upper k loop index for ghost point exchange at north boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nyn_send      !< upper k loop index for ghost point exchange at north boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nys_recv      !< upper k loop index for ghost point exchange at south boundary
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  ktop_for_nys_send      !< upper k loop index for ghost point exchange at south boundary

    LOGICAL ::  mg_switch_to_pe0 = .FALSE.  !< internal switch for steering the ghost point exchange
                                            !< in case that data has been collected on PE0
    LOGICAL, DIMENSION(:), ALLOCATABLE ::  unroll   !< flag indicating whether loop unrolling is possible

    REAL(wp) ::  residual_limit = 1.0E-6_wp  !< namelist parameter

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ddx2_mg  !< 1/dx_l**2 (dx_l: grid spacing along x on different multigrid level)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ddy2_mg  !< 1/dy_l**2 (dy_l: grid spacing along y on different multigrid level)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dzu_mg         !< vertical grid spacing (u-grid) for multigrid pressure solver
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dzw_mg         !< vertical grid spacing (w-grid) for multigrid pressure solver
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f1_mg          !< grid factor used in right hand side of Gauss-Seidel equation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f1_mg_b        !< blocked version of f1_mg
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f2_mg          !< grid factor used in right hand side of Gauss-Seidel equation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f2_mg_b        !< blocked version of f2_mg
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f3_mg          !< grid factor used in right hand side of Gauss-Seidel equation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f3_mg_b        !< blocked version of f3_mg
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_air_mg     !< air density profiles on the uv grid for multigrid levels
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_air_mg_b   !< blocked version of rho_air_mg
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_air_zw_mg  !< air density profiles on the w grid for multigrid levels

    TYPE ::  grid_level_flags
       INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  flags  !< topograpyh masking flag on a multigrid level
    END TYPE grid_level_flags

    TYPE(grid_level_flags), DIMENSION(:), ALLOCATABLE ::  gl  !< contains the masking flags for the multigrid levels

    SAVE

    INTERFACE poismg
       MODULE PROCEDURE poismg
    END INTERFACE poismg

    INTERFACE poismg_init
       MODULE PROCEDURE poismg_init
    END INTERFACE poismg_init

    INTERFACE sort_k_to_even_odd_blocks
       MODULE PROCEDURE sort_k_to_even_odd_blocks
       MODULE PROCEDURE sort_k_to_even_odd_blocks_int
       MODULE PROCEDURE sort_k_to_even_odd_blocks_1d
    END INTERFACE sort_k_to_even_odd_blocks

    PUBLIC

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid V- or W-Cycle scheme.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poismg

    USE arrays_3d,                                                                                 &
        ONLY:  d,                                                                                  &
               p_loc

    USE control_parameters,                                                                        &
        ONLY:  current_timestep_number,                                                            &
               ibc_p_t

    USE indices,                                                                                   &
        ONLY:  ngp_3d_inner

 
    IMPLICIT NONE

    INTEGER(iwp) ::  i           !< index variable along x
    INTEGER(iwp) ::  j           !< index variable along y
    INTEGER(iwp) ::  k           !< index variable along z
    INTEGER(iwp) ::  ngsrb_save  !< temporary variable to store the number of red/black iterations

    REAL(wp) ::  maxerror          !<
    REAL(wp) ::  maximum_mgcycles  !<
    REAL(wp) ::  residual_norm     !<

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  p3  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  r   !< residual


    CALL cpu_log( log_point_s(29), 'poismg', 'start' )

!
!-- Initialize arrays and variables used in this subroutine.
!-- If the number of grid points of the gathered grid, which is collected on PE0, is larger than
!-- the number of grid points of a PE, than arrays p3 and r will be enlarged.
    IF ( gathered_size > subdomain_size )  THEN
       ALLOCATE( p3(nzb:nzt_mg(mg_switch_to_pe0_level)+1,nys_mg(                                   &
                 mg_switch_to_pe0_level)-1:nyn_mg(mg_switch_to_pe0_level)+1,                       &
                 nxl_mg(mg_switch_to_pe0_level)-1:nxr_mg(mg_switch_to_pe0_level)+1) )
       ALLOCATE( r(nzb:nzt_mg(mg_switch_to_pe0_level)+1,nys_mg(                                    &
                 mg_switch_to_pe0_level)-1:nyn_mg(mg_switch_to_pe0_level)+1,                       &
                 nxl_mg(mg_switch_to_pe0_level)-1:nxr_mg(mg_switch_to_pe0_level)+1) )
    ELSE
       ALLOCATE( p3(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) )
       ALLOCATE( r(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) )
    ENDIF
    !$ACC DATA CREATE(p3,r) IF(enable_openacc)

    !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
    p3 = 0.0_wp
    !$ACC END KERNELS

!
!-- Ghost boundaries have to be added to divergence array.
!-- Exchange routine needs to know the grid level!
    grid_level = maximum_grid_level
    CALL exchange_horiz( d, 1, grid_level = grid_level )

!
!-- Set bottom and top boundary conditions.
    !$ACC KERNELS PRESENT(d) DEFAULT(NONE) IF(enable_openacc)
    d(nzb,:,:) = d(nzb+1,:,:)
    !$ACC END KERNELS

    IF ( ibc_p_t == 1 ) THEN
       !$ACC KERNELS PRESENT(d) DEFAULT(NONE) IF(enable_openacc)
       d(nzt+1,:,: ) = d(nzt,:,:)
       !$ACC END KERNELS
    ENDIF
!
!-- Set lateral boundary conditions in non-cyclic case.
    IF ( .NOT. bc_lr_cyc )  THEN
       !$ACC KERNELS PRESENT(d) DEFAULT(NONE) IF(enable_openacc)
       IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  d(:,:,nxl-1) = d(:,:,nxl)
       IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  d(:,:,nxr+1) = d(:,:,nxr)
       !$ACC END KERNELS
    ENDIF
    IF ( .NOT. bc_ns_cyc )  THEN
       !$ACC KERNELS PRESENT(d) DEFAULT(NONE) IF(enable_openacc)
       IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  d(:,nyn+1,:) = d(:,nyn,:)
       IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  d(:,nys-1,:) = d(:,nys,:)
       !$ACC END KERNELS
    ENDIF

!
!-- Initiation of the multigrid scheme. Does n cycles until the residual is smaller than the given
!-- limit. The accuracy of the solution of the poisson equation will increase with the number of
!-- cycles. If the number of cycles is preset by the user, this number will be carried out
!-- regardless of the accuracy.
    grid_level_count =  0
    mgcycles         =  0
    IF ( mg_cycles == -1 )  THEN
       maximum_mgcycles = 0
       residual_norm    = 1.0_wp
    ELSE
       maximum_mgcycles = mg_cycles
       residual_norm    = 0.0_wp
    ENDIF

!
!-- At the beginning of a run, large divergence may appear in the vicinity of
!-- topography/buildings. This may significantly increase errors in the advection schemes (for
!-- example generate large peaks in velocities).
!-- The divergence can be reduced by increasing the number of SOR-red-black iterations during the
!-- initial stage of a run.
    IF ( current_timestep_number <= ngsrb_initial_timesteps )  THEN
       ngsrb_save = ngsrb
       ngsrb = ngsrb_initial
    ENDIF

!
!-- Sort input arrays in even/odd blocks along k-dimension.
    CALL sort_k_to_even_odd_blocks( d, grid_level )
    CALL sort_k_to_even_odd_blocks( p_loc, grid_level )

!
!-- The complete multigrid cycles are running in block mode, i.e. over seperate data blocks of even
!-- and odd indices.
    DO WHILE ( residual_norm > residual_limit  .OR.  mgcycles < maximum_mgcycles )

       CALL next_mg_level( d, p_loc, p3, r)

!
!--    Calculate the residual if the user has not preset the number of cycles to be performed.
       IF ( maximum_mgcycles == 0 )  THEN

          CALL resid( d, p_loc, r )

          maxerror = 0.0_wp
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) REDUCTION(+:maxerror) IF(enable_openacc)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   maxerror = maxerror + r(k,j,i)**2
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( maxerror, residual_norm, 1, MPI_REAL, MPI_SUM, comm2d, ierr)
#else
          residual_norm = maxerror
#endif
          residual_norm = SQRT( residual_norm ) / REAL( ngp_3d_inner(0), KIND=wp )

       ENDIF

       mgcycles = mgcycles + 1

!
!--    If the user has not limited the number of cycles, stop the run in case of insufficient
!--    convergence. During the initial phase, do not stop the run but continue after more than 100
!--    cycles have been carried out, since in this phase the residuals sometimes converge very
!--    slowly (e.g. for childs).
       IF ( mg_cycles == -1 )  THEN
          IF ( current_timestep_number <= ngsrb_initial_timesteps  .AND.  mgcycles > 99 )  THEN
             EXIT
          ELSE
             IF ( mgcycles > 999 )  THEN
                WRITE( message_string,'(2A,E30.20)' ) 'no sufficient convergence within 1000 ',    &
                                                      'cycles,&mean residual: ', residual_norm
                CALL message( 'poismg', 'PAC0265', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
       ENDIF

    ENDDO

    !$ACC END DATA

    DEALLOCATE( p3 )

!
!-- For output purposes set pressure to zero at all topography/wall points.
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl-1, nxr+1
       DO  j = nys-1, nyn+1
          DO  k = nzb, nzt+1
             p_loc(k,j,i) = MERGE( 0.0_wp, p_loc(k,j,i), BTEST( gl(grid_level)%flags(k,j,i), 6 ) )
          ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP

!
!-- Result has to be sorted back from even/odd blocks to sequential order.
    CALL sort_k_to_sequential( p_loc )

!
!-- Add Neumann condition at top boundary because this is only implicitly set in redblack, but the
!-- pressure value is required at k=nzt+1 in pres to correct the vertical velocity at k=nzt.
!-- w at nzb is not corrected.
    IF ( ibc_p_t == 1 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       p_loc(nzt+1,:,:) = p_loc(nzt,:,:)
       !$ACC END KERNELS
    ENDIF
!
!-- For the same reason add Neumann conditions at non-cyclic lateral boundaries.
    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
          !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
          p_loc(:,:,nxl-1) = p_loc(:,:,nxl)
          !$ACC END KERNELS
       ENDIF
       IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
          !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
          p_loc(:,:,nxr+1) = p_loc(:,:,nxr)
          !$ACC END KERNELS
       ENDIF
    ENDIF

    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
          !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
          p_loc(:,nyn+1,:) = p_loc(:,nyn,:)
          !$ACC END KERNELS
       ENDIF
       IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
          !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
          p_loc(:,nys-1,:) = p_loc(:,nys,:)
          !$ACC END KERNELS
       ENDIF
    ENDIF

!
!-- Reset the number of SOR-red-black iterations to the value set by user (or default).
    IF ( current_timestep_number <= ngsrb_initial_timesteps )  ngsrb = ngsrb_save

    CALL cpu_log( log_point_s(29), 'poismg', 'stop' )

!    CALL MPI_FINALIZE( ierr )
!    STOP '*** poismg'

 END SUBROUTINE poismg


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the residual of the perturbation pressure.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resid( f_mg, p_mg, r )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< index variable along x
    INTEGER(iwp) ::  j    !< index variable along y
    INTEGER(iwp) ::  k    !< index variable along z
    INTEGER(iwp) ::  l    !< index indicating grid level
    INTEGER(iwp) ::  km1  !< index variable along z dimension (k-1)
    INTEGER(iwp) ::  kp1  !< index variable along z dimension (k+1)

    REAL(wp) ::  pkjim  !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkjip  !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkjmi  !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkjpi  !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkmji  !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkpji  !< pressure at respective grid point (m=-1,p=+1)

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !< velocity divergence
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  p_mg  !< perturbation pressure
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  r     !< residuum of perturbation pressure


    CALL cpu_log( log_point_s(53), 'resid', 'start' )

    l = grid_level

    !$OMP PARALLEL PRIVATE (i,j,k,km1,kp1)
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl_mg(l), nxr_mg(l)
       DO  j = nys_mg(l), nyn_mg(l)
          !DIR$ IVDEP
          !$ACC LOOP VECTOR
          DO  k = ind_even_odd+1, nzt_mg(l)
             km1 = k-ind_even_odd-1
             kp1 = k-ind_even_odd
             pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
             pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
             pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
             pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
             pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
             pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
             r(k,j,i) = f_mg(k,j,i) - rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    - rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    - f2_mg_b(k,l) * pkpji - f3_mg_b(k,l) * pkmji                  &
                                    + f1_mg_b(k,l) * p_mg(k,j,i)
!
!--          Residual within topography should be zero.
             r(k,j,i) = MERGE( 0.0_wp, r(k,j,i), BTEST( gl(l)%flags(k,j,i), 6 ) )
          ENDDO
          !DIR$ IVDEP
          !$ACC LOOP VECTOR
          DO  k = nzb+1, ind_even_odd
             km1 = k+ind_even_odd
             kp1 = k+ind_even_odd+1
             pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
             pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
             pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
             pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
             pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
             pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
             r(k,j,i) = f_mg(k,j,i) - rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    - rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    - f2_mg_b(k,l) * pkpji - f3_mg_b(k,l) * pkmji                  &
                                    + f1_mg_b(k,l) * p_mg(k,j,i)
!
!--          Residual within topography should be zero.
             r(k,j,i) = MERGE( 0.0_wp, r(k,j,i), BTEST( gl(l)%flags(k,j,i), 6 ) )
          ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL
!
!-- Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are implicitly
!-- treated via the flags array.
    CALL exchange_horiz( r, 1, grid_level = grid_level )

!
!-- Dirichlet boundary conditions at bottom and top of the domain. Neumann BCs are implicitly
!-- considered in the calculations above. Points may be within buildings, but that doesn't matter.
    IF ( ibc_p_b == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       r(nzb,:,: ) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    IF ( ibc_p_t == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       r(nzt_mg(l)+1,:,: ) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    CALL cpu_log( log_point_s(53), 'resid', 'stop' )

 END SUBROUTINE resid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the residual on the next coarser grid with "full weighting" scheme
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE restrict( f_mg, r )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< index variable along x on finer grid
    INTEGER(iwp) ::  ic   !< index variable along x on coarser grid
    INTEGER(iwp) ::  j    !< index variable along y on finer grid
    INTEGER(iwp) ::  jc   !< index variable along y on coarser grid
    INTEGER(iwp) ::  k    !< index variable along z on finer grid
    INTEGER(iwp) ::  kc   !< index variable along z on coarser grid
    INTEGER(iwp) ::  km1  !< index variable along z dimension (k-1 on finer level)
    INTEGER(iwp) ::  kp1  !< index variable along z dimension (k+1 on finer level)
    INTEGER(iwp) ::  l    !< index indicating the grid level

    REAL(wp) ::  rkjim    !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjip    !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjmi    !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjmim   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjmip   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjpi    !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjpim   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkjpip   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmji    !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjim   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjip   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjmi   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjmim  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjmip  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjpi   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjpim  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkmjpip  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpji    !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjim   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjip   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjmi   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjmim  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjmip  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjpi   !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjpim  !< residual at respective grid point (m=-1,p=+1)
    REAL(wp) ::  rkpjpip  !< residual at respective grid point (m=-1,p=+1)

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !< residual on coarser grid level

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level+1)+1,nys_mg(grid_level+1)-1:nyn_mg(grid_level+1)+1,  &
                        nxl_mg(grid_level+1)-1:nxr_mg(grid_level+1)+1) ::  r  !< residual on finer grid level


    l   = grid_level

    CALL cpu_log( log_point_s(54), 'restrict', 'start' )

    !$OMP PARALLEL PRIVATE (i,j,k,ic,jc,kc,km1,kp1)
    !$OMP DO SCHEDULE( STATIC )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  ic = nxl_mg(l), nxr_mg(l)
       DO  jc = nys_mg(l), nyn_mg(l)
          i = 2 * ic
!
!--       Calculation for the first point along k.
          j  = 2 * jc
!
!--       Calculation for the other points along k (fine grid index at this point).
!--       kc is the coarse grid index.
          !DIR$ IVDEP
          !$ACC LOOP VECTOR
          DO  k = ind_even_odd+1, nzt_mg(l+1)
             km1 = k-ind_even_odd-1
             kp1 = k-ind_even_odd
             kc  = k-ind_even_odd
!
!--          Use implicit Neumann BCs if the respective gridpoint is inside the building.
             rkjim   = MERGE( r(k,j,i), r(k,j,i-1), BTEST( gl(l+1)%flags(k,j,i-1), 6 ) )
             rkjip   = MERGE( r(k,j,i), r(k,j,i+1), BTEST( gl(l+1)%flags(k,j,i+1), 6 ) )
             rkjpi   = MERGE( r(k,j,i), r(k,j+1,i), BTEST( gl(l+1)%flags(k,j+1,i), 6 ) )
             rkjmi   = MERGE( r(k,j,i), r(k,j-1,i), BTEST( gl(l+1)%flags(k,j-1,i), 6 ) )
             rkjmim  = MERGE( r(k,j,i), r(k,j-1,i-1), BTEST( gl(l+1)%flags(k,j-1,i-1), 6 ) )
             rkjpim  = MERGE( r(k,j,i), r(k,j+1,i-1), BTEST( gl(l+1)%flags(k,j+1,i-1), 6 ) )
             rkjmip  = MERGE( r(k,j,i), r(k,j-1,i+1), BTEST( gl(l+1)%flags(k,j-1,i+1), 6 ) )
             rkjpip  = MERGE( r(k,j,i), r(k,j+1,i+1), BTEST( gl(l+1)%flags(k,j+1,i+1), 6 ) )
             rkmji   = MERGE( r(k,j,i), r(km1,j,i), BTEST( gl(l+1)%flags(km1,j,i), 6 ) )
             rkmjim  = MERGE( r(k,j,i), r(km1,j,i-1), BTEST( gl(l+1)%flags(km1,j,i-1), 6 ) )
             rkmjip  = MERGE( r(k,j,i), r(km1,j,i+1), BTEST( gl(l+1)%flags(km1,j,i+1), 6 ) )
             rkmjpi  = MERGE( r(k,j,i), r(km1,j+1,i), BTEST( gl(l+1)%flags(km1,j+1,i), 6 ) )
             rkmjmi  = MERGE( r(k,j,i), r(km1,j-1,i), BTEST( gl(l+1)%flags(km1,j-1,i), 6 ) )
             rkmjmim = MERGE( r(k,j,i), r(km1,j-1,i-1), BTEST( gl(l+1)%flags(km1,j-1,i-1), 6 ) )
             rkmjpim = MERGE( r(k,j,i), r(km1,j+1,i-1), BTEST( gl(l+1)%flags(km1,j+1,i-1), 6 ) )
             rkmjmip = MERGE( r(k,j,i), r(km1,j-1,i+1), BTEST( gl(l+1)%flags(km1,j-1,i+1), 6 ) )
             rkmjpip = MERGE( r(k,j,i), r(km1,j+1,i+1), BTEST( gl(l+1)%flags(km1,j+1,i+1), 6 ) )
             rkpji   = MERGE( r(k,j,i), r(kp1,j,i), BTEST( gl(l+1)%flags(kp1,j,i), 6 ) )
             rkpjim  = MERGE( r(k,j,i), r(kp1,j,i-1), BTEST( gl(l+1)%flags(kp1,j,i-1), 6 ) )
             rkpjip  = MERGE( r(k,j,i), r(kp1,j,i+1), BTEST( gl(l+1)%flags(kp1,j,i+1), 6 ) )
             rkpjpi  = MERGE( r(k,j,i), r(kp1,j+1,i), BTEST( gl(l+1)%flags(kp1,j+1,i), 6 ) )
             rkpjmi  = MERGE( r(k,j,i), r(kp1,j-1,i), BTEST( gl(l+1)%flags(kp1,j-1,i), 6 ) )
             rkpjmim = MERGE( r(k,j,i), r(kp1,j-1,i-1), BTEST( gl(l+1)%flags(kp1,j-1,i-1), 6 ) )
             rkpjpim = MERGE( r(k,j,i), r(kp1,j+1,i-1), BTEST( gl(l+1)%flags(kp1,j+1,i-1), 6 ) )
             rkpjmip = MERGE( r(k,j,i), r(kp1,j-1,i+1), BTEST( gl(l+1)%flags(kp1,j-1,i+1), 6 ) )
             rkpjpip = MERGE( r(k,j,i), r(kp1,j+1,i+1), BTEST( gl(l+1)%flags(kp1,j+1,i+1), 6 ) )
             f_mg(kc,jc,ic) = 1.0_wp / 64.0_wp *                                                   &
                              ( 8.0_wp * r(k,j,i)                                                  &
                              + 4.0_wp * ( rkjim  + rkjip  + rkjpi   + rkjmi  )                    &
                              + 2.0_wp * ( rkjmim + rkjpim + rkjmip  + rkjpip )                    &
                              + 4.0_wp * rkmji                                                     &
                              + 2.0_wp * ( rkmjim  + rkmjim  + rkmjpi  + rkmjmi  )                 &
                              +          ( rkmjmim + rkmjpim + rkmjmip + rkmjpip )                 &
                              + 4.0_wp * rkpji                                                     &
                              + 2.0_wp * ( rkpjim  + rkpjim  + rkpjpi  + rkpjmi  )                 &
                              +          ( rkpjmim + rkpjpim + rkpjmip + rkpjpip )                 &
                              )
          ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
    !$OMP ENDDO
    !$OMP END PARALLEL

!
!-- Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are implicitly
!-- treated via the flags array.
    CALL exchange_horiz( f_mg, 1, grid_level = grid_level, mg_switch_to_pe0 = mg_switch_to_pe0 )

!
!-- Dirichlet boundary conditions at bottom and top of the domain. Neumann BCs are implicitly
!-- considered in the calculations above. Points may be within buildings, but that doesn't matter.
!-- Note that here f_mg is ordered sequentielly after interpolation on coarse grid. It will be
!-- ordered in odd-even blocks further below.
    IF ( ibc_p_b == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f_mg(nzb,:,:) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    IF ( ibc_p_t == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f_mg(nzt_mg(l)+1,:,:) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

!
!-- Since residual is in sequential order after interpolation, an additional sorting in odd-even
!-- blocks along z dimension is required at this point.
    CALL sort_k_to_even_odd_blocks( f_mg , l)

    CALL cpu_log( log_point_s(54), 'restrict', 'stop' )

 END SUBROUTINE restrict


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the correction of the perturbation pressure to the next finer grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE prolong( p, temp )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !< index variable along x on coarser grid level
    INTEGER(iwp) ::  j   !< index variable along y on coarser grid level
    INTEGER(iwp) ::  k   !< index variable along z on coarser grid level
    INTEGER(iwp) ::  l   !< index indicating finer grid level
    INTEGER(iwp) ::  lm1 !< index for flags indicating coarser grid level (considering the switch to PE0 level)
    INTEGER(iwp) ::  ke  !< index for prolog even
    INTEGER(iwp) ::  ko  !< index for prolog odd
    INTEGER(iwp) ::  kp1 !< index variable along z

    REAL(wp) ::  pkjip    !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkjpi    !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkjpip   !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkpji    !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkpjip   !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkpjpi   !< pressure at respective grid point (m=-1,p=+1)
    REAL(wp) ::  pkpjpip  !< pressure at respective grid point (m=-1,p=+1)

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,  &
                        nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1 ) ::  p  !< perturbation pressure on coarser grid level

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  temp  !< perturbation pressure on finer grid level


    CALL cpu_log( log_point_s(55), 'prolong', 'start' )

    l = grid_level
    ind_even_odd = even_odd_level(grid_level-1)
!
!-- Choose index for the lower level flag array.
    lm1 = grid_level - 1
!
!-- A special index 0 is required when switching from the total domain (switch_to_pe0_level)
!-- to the next finer level, because here the prolongation already calculates on the subdomains
!-- and not on the total domain. The regular flag-array for this level is defined for the total
!-- domain.
    IF ( ( l-1 ) == mg_switch_to_pe0_level )  lm1 = 0

    !$OMP PARALLEL PRIVATE (i,j,k,kp1,ke,ko)
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl_mg(l-1), nxr_mg(l-1)
       DO  j = nys_mg(l-1), nyn_mg(l-1)

          !DIR$ IVDEP
          DO  k = ind_even_odd+1, nzt_mg(l-1)
             kp1 = k - ind_even_odd
             ke  = 2 * ( k-ind_even_odd - 1 ) + 1
             ko  = 2 * k - 1
!
!--          Store pressure at surrounding grid points and apply Neumann boundary conditions in
!--          case of a wall.
             pkjip   = MERGE( p(k,j,i), p(k,j,i+1),     BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkjpi   = MERGE( p(k,j,i), p(k,j+1,i),     BTEST( gl(lm1)%flags(k,j,i), 3 ) )
             pkpji   = MERGE( p(k,j,i), p(kp1,j,i),     BTEST( gl(lm1)%flags(k,j,i), 1 ) )
             pkjpip  = MERGE( p(k,j,i), p(k,j+1,i+1),   BTEST( gl(lm1)%flags(k,j,i), 3 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkpjip  = MERGE( p(k,j,i), p(kp1,j,i+1),   BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkpjpi  = MERGE( p(k,j,i), p(kp1,j+1,i),   BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 3 ) )
             pkpjpip = MERGE( p(k,j,i), p(kp1,j+1,i+1), BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 3 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
!
!--          Points of the coarse grid are directly stored on the next finer grid.
             temp(ko,2*j,2*i) = p(k,j,i)
             temp(ko,2*j,2*i) = MERGE( 0.0_wp, temp(ko,2*j,2*i), BTEST( gl(lm1)%flags(k,j,i), 6 ) )
!
!--          Points between two coarse-grid points.
             temp(ko,2*j,2*i+1) = 0.5_wp * ( p(k,j,i) + pkjip )
             temp(ko,2*j,2*i+1) = MERGE( 0.0_wp, temp(ko,2*j,2*i+1),                               &
                                         BTEST( gl(lm1)%flags(k,j,i), 7 ) )
             temp(ko,2*j+1,2*i) = 0.5_wp * ( p(k,j,i) + pkjpi )
             temp(ko,2*j+1,2*i) = MERGE( 0.0_wp, temp(ko,2*j+1,2*i),                               &
                                         BTEST( gl(lm1)%flags(k,j,i), 8 ) )
             temp(ke,2*j,2*i)   = 0.5_wp * ( p(k,j,i) + pkpji )
             temp(ke,2*j,2*i)   = MERGE( 0.0_wp, temp(ke,2*j,2*i),                                 &
                                         BTEST( gl(lm1)%flags(k,j,i), 9 ) )
!
!--          Points in the center of the planes stretched by four points of the coarse grid cube.
             temp(ko,2*j+1,2*i+1) = 0.25_wp * ( p(k,j,i) + pkjip + pkjpi + pkjpip )
             temp(ko,2*j+1,2*i+1) = MERGE( 0.0_wp, temp(ko,2*j+1,2*i+1),                           &
                                           BTEST( gl(lm1)%flags(k,j,i), 10 ) )
             temp(ke,2*j,2*i+1)   = 0.25_wp * ( p(k,j,i) + pkjip + pkpji + pkpjip )
             temp(ke,2*j,2*i+1)   = MERGE( 0.0_wp, temp(ke,2*j,2*i+1),                             &
                                           BTEST( gl(lm1)%flags(k,j,i), 11 ) )
             temp(ke,2*j+1,2*i)   = 0.25_wp * ( p(k,j,i) + pkjpi + pkpji + pkpjpi )
             temp(ke,2*j+1,2*i)   = MERGE( 0.0_wp, temp(ke,2*j+1,2*i),                             &
                                           BTEST( gl(lm1)%flags(k,j,i), 12 ) )
!
!--          Points in the middle of coarse grid cube.
             temp(ke,2*j+1,2*i+1) = 0.125_wp * ( p(k,j,i) + pkjip + pkjpi  + pkjpip +              &
                                                            pkpji + pkpjip + pkpjpi + pkpjpip )
             temp(ke,2*j+1,2*i+1) = MERGE( 0.0_wp, temp(ke,2*j+1,2*i+1),                           &
                                                   BTEST( gl(lm1)%flags(k,j,i), 13 ) )
          ENDDO

          !DIR$ IVDEP
          DO  k = nzb+1, ind_even_odd
             kp1 = k + ind_even_odd + 1
             ke  = 2 * k
             ko  = 2 * ( k + ind_even_odd )

             pkjip   = MERGE( p(k,j,i), p(k,j,i+1),     BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkjpi   = MERGE( p(k,j,i), p(k,j+1,i),     BTEST( gl(lm1)%flags(k,j,i), 3 ) )
             pkpji   = MERGE( p(k,j,i), p(kp1,j,i),     BTEST( gl(lm1)%flags(k,j,i), 1 ) )
             pkjpip  = MERGE( p(k,j,i), p(k,j+1,i+1),   BTEST( gl(lm1)%flags(k,j,i), 3 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkpjip  = MERGE( p(k,j,i), p(kp1,j,i+1),   BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkpjpi  = MERGE( p(k,j,i), p(kp1,j+1,i),   BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 3 ) )
             pkpjpip = MERGE( p(k,j,i), p(kp1,j+1,i+1), BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 3 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
!
!--          Points of the coarse grid are directly stored on the next finer grid.
             temp(ko,2*j,2*i) = p(k,j,i)
             temp(ko,2*j,2*i) = MERGE( 0.0_wp, temp(ko,2*j,2*i), BTEST( gl(lm1)%flags(k,j,i), 6 ) )
!
!--          Points between two coarse-grid points.
             temp(ko,2*j,2*i+1) = 0.5_wp * ( p(k,j,i) + pkjip )
             temp(ko,2*j,2*i+1) = MERGE( 0.0_wp, temp(ko,2*j,2*i+1),                               &
                                         BTEST( gl(lm1)%flags(k,j,i), 7 ) )
             temp(ko,2*j+1,2*i) = 0.5_wp * ( p(k,j,i) + pkjpi )
             temp(ko,2*j+1,2*i) = MERGE( 0.0_wp, temp(ko,2*j+1,2*i),                               &
                                         BTEST( gl(lm1)%flags(k,j,i), 8 ) )
             temp(ke,2*j,2*i)   = 0.5_wp * ( p(k,j,i) + pkpji )
             temp(ke,2*j,2*i)   = MERGE( 0.0_wp, temp(ke,2*j,2*i),                                 &
                                         BTEST( gl(lm1)%flags(k,j,i), 9 ) )
!
!--          Points in the center of the planes stretched by four points of the coarse grid cube.
             temp(ko,2*j+1,2*i+1) = 0.25_wp * ( p(k,j,i) + pkjip + pkjpi + pkjpip )
             temp(ko,2*j+1,2*i+1) = MERGE( 0.0_wp, temp(ko,2*j+1,2*i+1),                           &
                                           BTEST( gl(lm1)%flags(k,j,i), 10 ) )
             temp(ke,2*j,2*i+1)   = 0.25_wp * ( p(k,j,i) + pkjip + pkpji + pkpjip )
             temp(ke,2*j,2*i+1)   = MERGE( 0.0_wp, temp(ke,2*j,2*i+1),                             &
                                           BTEST( gl(lm1)%flags(k,j,i), 11 ) )
             temp(ke,2*j+1,2*i)   = 0.25_wp * ( p(k,j,i) + pkjpi + pkpji + pkpjpi )
             temp(ke,2*j+1,2*i)   = MERGE( 0.0_wp, temp(ke,2*j+1,2*i),                             &
                                           BTEST( gl(lm1)%flags(k,j,i), 12 ) )
!
!--          Points in the middle of coarse grid cube.
             temp(ke,2*j+1,2*i+1) = 0.125_wp * ( p(k,j,i) + pkjip + pkjpi  + pkjpip +              &
                                                            pkpji + pkpjip + pkpjpi + pkpjpip )
             temp(ke,2*j+1,2*i+1) = MERGE( 0.0_wp, temp(ke,2*j+1,2*i+1),                           &
                                           BTEST( gl(lm1)%flags(k,j,i), 13 ) )
          ENDDO

       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL

    ind_even_odd = even_odd_level(grid_level)
!
!-- Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are implicitly
!-- treated via the flags array.
    CALL exchange_horiz( temp, 1, grid_level = grid_level, mg_switch_to_pe0 = mg_switch_to_pe0 )

!
!-- Dirichlet boundary conditions at bottom and top of the domain. Neumann BCs are implicitly
!-- considered in the calculations above. Points may be within buildings, but that doesn't matter.
    IF ( ibc_p_b == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       temp(nzb,:,: ) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    IF ( ibc_p_t == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       temp(nzt_mg(l)+1,:,: ) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    CALL cpu_log( log_point_s(55), 'prolong', 'stop' )

 END SUBROUTINE prolong


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Relaxation method for the multigrid scheme. A Gauss-Seidel iteration with 3D-Red-Black
!> decomposition (GS-RB) is used.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE redblack( f_mg, p_mg )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  color  !< grid point color, either red or black
    INTEGER(iwp) ::  i      !< index variable along x
    INTEGER(iwp) ::  ic     !< index variable along x
    INTEGER(iwp) ::  j      !< index variable along y
    INTEGER(iwp) ::  jc     !< index variable along y
    INTEGER(iwp) ::  jj     !< index variable along y
    INTEGER(iwp) ::  k      !< index variable along z
    INTEGER(iwp) ::  km1    !< index variable (k-1)
    INTEGER(iwp) ::  kp1    !< index variable (k+1)
    INTEGER(iwp) ::  l      !< grid level
    INTEGER(iwp) ::  n      !< loop variable Gauß-Seidel iterations
    INTEGER(iwp) ::  save_nxl_mg  !< to save nxl_mg on coarsest level 1
    INTEGER(iwp) ::  save_nys_mg  !< to save nys_mg on coarsest level 1


    LOGICAL ::  adjust_lower_i_index  !< adjust lower limit of i loop in case of odd number of grid points
    LOGICAL ::  adjust_lower_j_index  !< adjust lower limit of j loop in case of odd number of grid points

    REAL(wp) ::  pkjim  !< pressure left of i,j,k
    REAL(wp) ::  pkjip  !< pressure right of i,j,k
    REAL(wp) ::  pkjmi  !< pressure south of i,j,k
    REAL(wp) ::  pkjpi  !< pressure north of i,j,k
    REAL(wp) ::  pkmji  !< pressure below i,j,k
    REAL(wp) ::  pkpji  !< pressure above i,j,k

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !< residual of perturbation pressure
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  p_mg  !< perturbation pressure


!    p_mg(:,:,:) = 0.0
    l = grid_level

    unroll(l) = ( MOD( nyn_mg(l)-nys_mg(l)+1, 4 ) == 0  .AND.  MOD( nxr_mg(l)-nxl_mg(l)+1, 2 ) == 0 )

!
!-- The red/black decomposition requires that on the lower i,j indices need to start alternatively
!-- with an even or odd value on the coarsest grid level, depending on the core-id, and if the
!-- subdomain has an uneven number of gridpoints along x/y. The respective steering switches
!-- are set here.
    IF ( l == 1  .AND.  MOD( myidx, 2 ) /= 0  .AND.  MOD( nxl_mg(l) - nxr_mg(l), 2 ) == 0 )  THEN
       adjust_lower_i_index = .TRUE.
       save_nxl_mg = nxl_mg(l)
    ELSE
       adjust_lower_i_index = .FALSE.
    ENDIF
    IF ( l == 1  .AND.  MOD( myidy, 2 ) /= 0  .AND.  MOD( nyn_mg(l) - nys_mg(l), 2 ) == 0 )  THEN
       adjust_lower_j_index = .TRUE.
       save_nys_mg = nys_mg(l)
    ELSE
       adjust_lower_j_index = .FALSE.
    ENDIF


!    WRITE(9,*) ' '
!    WRITE(9,*) '*** level = ', l
    DO  n = 1, ngsrb

       DO  color = 1, 2

          IF ( .NOT. unroll(l) )  THEN

             CALL cpu_log( log_point_s(36), 'redblack_no_unroll_f', 'start' )

             IF ( adjust_lower_i_index )  THEN
                nxl_mg(l) = save_nxl_mg + 1
             ENDIF

             IF ( adjust_lower_j_index )  THEN
                IF ( color == 1 )  THEN
                  nys_mg(l) = save_nys_mg - 1
                ELSE
                  nys_mg(l) = save_nys_mg + 1
                ENDIF
             ENDIF
             !$ACC UPDATE DEVICE (nxl_mg,nys_mg) IF(enable_openacc)
!
!--          Without unrolling of loops, no cache optimization.
!             WRITE(9,*) '*** loop 1:'
!             WRITE(9,'(A,I3,A,I3,A)') 'i = ', nxl_mg(l), ', ', nxr_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'j = ', nys_mg(l) + 2 - color, ', ', nyn_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'k = ', ind_even_odd+1, ', ', nzt_mg(l), ', 1'
             !$OMP PARALLEL PRIVATE (i,j,k,km1,kp1)
             !$OMP DO
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
             !$ACC DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl_mg(l), nxr_mg(l), 2
                DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2
                   !DIR$ IVDEP
!                   !$ACC LOOP VECTOR
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      km1 = k-ind_even_odd-1
                      kp1 = k-ind_even_odd
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP

             IF ( adjust_lower_i_index )  THEN
                nxl_mg(l) = save_nxl_mg - 1
             ENDIF

             IF ( adjust_lower_j_index )  THEN
                IF ( color == 1 )  THEN
                  nys_mg(l) = save_nys_mg + 1
                ELSE
                  nys_mg(l) = save_nys_mg - 1
                ENDIF
             ENDIF
             !$ACC UPDATE DEVICE (nxl_mg,nys_mg) IF(enable_openacc)

!             WRITE(9,*) '*** loop 2:'
!             WRITE(9,'(A,I3,A,I3,A)') 'i = ', nxl_mg(l)+1, ', ', nxr_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'j = ', nys_mg(l) + (color-1), ', ', nyn_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'k = ', ind_even_odd+1, ', ', nzt_mg(l), ', 1'
             !$OMP DO
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
             !$ACC DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl_mg(l)+1, nxr_mg(l), 2
                DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2
                    !DIR$ IVDEP
!                    !$ACC LOOP VECTOR
                    DO  k = ind_even_odd+1, nzt_mg(l)
                      km1 = k-ind_even_odd-1
                      kp1 = k-ind_even_odd
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP

             IF ( adjust_lower_i_index )  THEN
                nxl_mg(l) = save_nxl_mg + 1
             ENDIF

             IF ( adjust_lower_j_index )  THEN
                IF ( color == 1 )  THEN
                  nys_mg(l) = save_nys_mg + 1
                ELSE
                  nys_mg(l) = save_nys_mg - 1
                ENDIF
             ENDIF
             !$ACC UPDATE DEVICE (nxl_mg,nys_mg) IF(enable_openacc)

!             WRITE(9,*) '*** loop 3:'
!             WRITE(9,'(A,I3,A,I3,A)') 'i = ', nxl_mg(l), ', ', nxr_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'j = ', nys_mg(l) + (color-1), ', ', nyn_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'k = ', nzb+1, ', ', ind_even_odd, ', 1'
             !$OMP DO
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
             !$ACC DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl_mg(l), nxr_mg(l), 2
                DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2
                   !DIR$ IVDEP
!                   !$ACC LOOP VECTOR
                   DO  k = nzb+1, ind_even_odd
                      km1 = k+ind_even_odd
                      kp1 = k+ind_even_odd+1
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP

             IF ( adjust_lower_i_index )  THEN
                nxl_mg(l) = save_nxl_mg - 1
             ENDIF

             IF ( adjust_lower_j_index )  THEN
                IF ( color == 1 )  THEN
                  nys_mg(l) = save_nys_mg - 1
                ELSE
                  nys_mg(l) = save_nys_mg + 1
                ENDIF
             ENDIF
             !$ACC UPDATE DEVICE (nxl_mg,nys_mg) IF(enable_openacc)

!             WRITE(9,*) '*** loop 4:'
!             WRITE(9,'(A,I3,A,I3,A)') 'i = ', nxl_mg(l)+1, ', ', nxr_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'j = ', nys_mg(l) + 2 - color, ', ', nyn_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'k = ', nzb+1, ', ', ind_even_odd, ', 1'
             !$OMP DO
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
             !$ACC DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl_mg(l)+1, nxr_mg(l), 2
                DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2
                   !DIR$ IVDEP
!                   !$ACC LOOP VECTOR
                   DO  k = nzb+1, ind_even_odd
                      km1 = k+ind_even_odd
                      kp1 = k+ind_even_odd+1
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
             !$OMP END PARALLEL

             CALL cpu_log( log_point_s(36), 'redblack_no_unroll_f', 'stop' )

          ELSE
!
!--          Loop unrolling along y, only one i loop for better cache use
             CALL cpu_log( log_point_s(38), 'redblack_unroll_f', 'start' )

             !$OMP PARALLEL PRIVATE (i,j,k,ic,jc,km1,kp1,jj)
             !$OMP DO
             !$ACC PARALLEL LOOP GANG COLLAPSE(2) &
             !$ACC DEFAULT(PRESENT) IF(enable_openacc)
             DO  ic = nxl_mg(l), nxr_mg(l), 2
                DO  jc = nys_mg(l), nyn_mg(l), 4
                   i  = ic
                   jj = jc+2-color
                   !DIR$ IVDEP
                   !$ACC LOOP VECTOR
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      km1 = k-ind_even_odd-1
                      kp1 = k-ind_even_odd
                      j   = jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                   i  = ic+1
                   jj = jc+color-1
                   !DIR$ IVDEP
                   !$ACC LOOP VECTOR
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      km1 = k-ind_even_odd-1
                      kp1 = k-ind_even_odd
                      j   = jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                   i  = ic
                   jj = jc+color-1
                   !DIR$ IVDEP
                   !$ACC LOOP VECTOR
                   DO  k = nzb+1, ind_even_odd
                      km1 = k+ind_even_odd
                      kp1 = k+ind_even_odd+1
                      j   = jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                   i  = ic+1
                   jj = jc+2-color
                   !DIR$ IVDEP
                   !$ACC LOOP VECTOR
                   DO  k = nzb+1, ind_even_odd
                      km1 = k+ind_even_odd
                      kp1 = k+ind_even_odd+1
                      j   = jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(kp1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(km1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) *                                        &
                                    ( rho_air_mg_b(k,l) * ddx2_mg(l) * ( pkjip + pkjim )           &
                                    + rho_air_mg_b(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )           &
                                    + f2_mg_b(k,l) * pkpji  + f3_mg_b(k,l) * pkmji                 &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
             !$OMP END PARALLEL

             CALL cpu_log( log_point_s(38), 'redblack_unroll_f', 'stop' )

          ENDIF

!          !$ACC END DATA
!          DO  i = 1, 4
!             WRITE(9,*) '*** exchange loop ', i,' :'
!             WRITE(9,'(A,I3,A,I3,A)') 'i = ', ileft(i,color,l), ', ', nxr_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'j = ', jsouth(i,color,l), ', ', nyn_mg(l), ', 2'
!             WRITE(9,'(A,I3,A,I3,A)') 'k = ', kbottom(i,color,l), ', ', ktop(i,color,l), ', 1'
!          ENDDO
!
!--       Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are
!--       implicitly treated via the flags array.  In case of sufficiently large data,
!--       contiguous buffers are used in exchange_horiz_rb to only exchange data of the respective
!--       color. The threshold of 900 is empirical and may require adjustment to optimize
!--       performance.
!--       Levels where total domain is on PE0 do not require optimized exchange.
#if defined( __parallel )
          IF ( ( ngp_xz(l) >= 900  .OR.  ngp_yz(l) >= 900 )  .AND.  .NOT. mg_switch_to_pe0  .AND.  &
               npex /= 1  .AND.  npey /= 1  )                                                      &
          THEN
             CALL exchange_horiz_rb( p_mg, 1, color = color, kinc = 1,                             &
                                     ileft_for_nys_send = ileft_for_nys_send(:,:,l),               &
                                     ileft_for_nys_recv = ileft_for_nys_recv(:,:,l),               &
                                     ileft_for_nyn_send = ileft_for_nyn_send(:,:,l),               &
                                     ileft_for_nyn_recv = ileft_for_nyn_recv(:,:,l),               &
                                     jsouth_for_nxl_send = jsouth_for_nxl_send(:,:,l),             &
                                     jsouth_for_nxl_recv = jsouth_for_nxl_recv(:,:,l),             &
                                     jsouth_for_nxr_send = jsouth_for_nxr_send(:,:,l),             &
                                     jsouth_for_nxr_recv = jsouth_for_nxr_recv(:,:,l),             &
                                     kbottom_for_nys_send = kbottom_for_nys_send(:,:,l),           &
                                     kbottom_for_nys_recv = kbottom_for_nys_recv(:,:,l),           &
                                     kbottom_for_nyn_send = kbottom_for_nyn_send(:,:,l),           &
                                     kbottom_for_nyn_recv = kbottom_for_nyn_recv(:,:,l),           &
                                     kbottom_for_nxl_send = kbottom_for_nxl_send(:,:,l),           &
                                     kbottom_for_nxl_recv = kbottom_for_nxl_recv(:,:,l),           &
                                     kbottom_for_nxr_send = kbottom_for_nxr_send(:,:,l),           &
                                     kbottom_for_nxr_recv = kbottom_for_nxr_recv(:,:,l),           &
                                     ktop_for_nys_send = ktop_for_nys_send(:,:,l),                 &
                                     ktop_for_nys_recv = ktop_for_nys_recv(:,:,l),                 &
                                     ktop_for_nyn_send = ktop_for_nyn_send(:,:,l),                 &
                                     ktop_for_nyn_recv = ktop_for_nyn_recv(:,:,l),                 &
                                     ktop_for_nxl_send = ktop_for_nxl_send(:,:,l),                 &
                                     ktop_for_nxl_recv = ktop_for_nxl_recv(:,:,l),                 &
                                     ktop_for_nxr_send = ktop_for_nxr_send(:,:,l),                 &
                                     ktop_for_nxr_recv = ktop_for_nxr_recv(:,:,l) )
          ELSE
             CALL exchange_horiz( p_mg, 1, grid_level = grid_level,                                &
                                  mg_switch_to_pe0 = mg_switch_to_pe0 )
          ENDIF
#else
          CALL exchange_horiz( p_mg, 1, grid_level = grid_level,                                   &
                               mg_switch_to_pe0 = mg_switch_to_pe0 )
#endif

!
!--       Dirichlet boundary conditions at bottom and top of the domain. Neumann BCs are implicitly
!--       considered in the calculations above. Points may be within buildings, but that doesn't
!--       matter.
          IF ( ibc_p_b == 0 )  THEN
             !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
             p_mg(nzb,:,: ) = 0.0_wp
             !$ACC END KERNELS
          ENDIF

          IF ( ibc_p_t == 0 )  THEN
             !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
             p_mg(nzt_mg(l)+1,:,: ) = 0.0_wp
             !$ACC END KERNELS
          ENDIF

!          WRITE(9,*) '*** level = ', l, '  color = ', color, '  unroll = ', unroll(l)
!          IF ( .NOT. unroll )  WRITE(9,*)  '    adjust_i = ', adjust_lower_i_index, '  adjust_j = ', adjust_lower_j_index
!          IF ( adjust_lower_i_index )  THEN
!             WRITE(9,*) '    nxl = ', save_nxl_mg, '  nxr = ', nxr_mg(l)
!          ELSE
!             WRITE(9,*) '    nxl = ', nxl_mg(l), '  nxr = ', nxr_mg(l)
!          ENDIF
!          IF ( adjust_lower_j_index )  THEN
!             WRITE(9,*) '    nys = ', save_nys_mg, '  nyn = ', nyn_mg(l)
!          ELSE
!             WRITE(9,*) '    nys = ', nys_mg(l), '  nyn = ', nyn_mg(l)
!          ENDIF
!          DO  k = nzb, nzt_mg(l)+1

!             IF ( adjust_lower_j_index )  THEN
!                DO  j = nyn_mg(grid_level)+1, save_nys_mg-1, -1
!                   WRITE(9,*) '*** j = ', j
!                   DO  k = nzb, nzt_mg(l)+1
!                      IF ( adjust_lower_i_index )  THEN
!                         WRITE(9,'(A,I3,1X,70I1)')  'k=', k, ( NINT( p_mg(k,j,i) ), i = save_nxl_mg-1, nxr_mg(grid_level)+1 )
!                      ELSE
!                         WRITE(9,'(A,I3,1X,70I1)')  'k=', k, ( NINT( p_mg(k,j,i) ), i = nxl_mg(grid_level)-1, nxr_mg(grid_level)+1 )
!                      ENDIF
!                   ENDDO
!                ENDDO
!             ELSE
!                DO  j = nyn_mg(grid_level)+1, nys_mg(grid_level)-1, -1
!                   WRITE(9,*) '*** j = ', j
!                   DO  k = nzb, nzt_mg(l)+1
!                      IF ( adjust_lower_i_index )  THEN
!                         WRITE(9,'(A,I3,1X,70I1)')  'k=', k, ( NINT( p_mg(k,j,i) ), i = save_nxl_mg-1, nxr_mg(grid_level)+1 )
!                      ELSE
!                         WRITE(9,'(A,I3,1X,70I1)')  'k=', k, ( NINT( p_mg(k,j,i) ), i = nxl_mg(grid_level)-1, nxr_mg(grid_level)+1 )
!                      ENDIF
!                   ENDDO
!                ENDDO
!             ENDIF

       ENDDO

    ENDDO

!
!-- Reset lower index limits to their standard values (may happen on coarsest levels only)
    IF ( adjust_lower_i_index )  THEN
       nxl_mg(l) = save_nxl_mg
       !$ACC UPDATE DEVICE (nxl_mg) IF(enable_openacc)
    ENDIF

    IF ( adjust_lower_j_index )  THEN
       nys_mg(l) = save_nys_mg
       !$ACC UPDATE DEVICE (nys_mg) IF(enable_openacc)
    ENDIF

 END SUBROUTINE redblack


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-dimension from sequential into blocks of even and odd. This is required to vectorize the
!> red-black subroutine. Version for 3D-REAL arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sort_k_to_even_odd_blocks( p_mg , glevel )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  glevel  !< grid level

    REAL(wp), DIMENSION(nzb:nzt_mg(glevel)+1,nys_mg(glevel)-1:nyn_mg(glevel)+1,                    &
                        nxl_mg(glevel)-1:nxr_mg(glevel)+1) ::  p_mg  !< array to be sorted

    INTEGER(iwp) ::  i    !< index variable along x
    INTEGER(iwp) ::  ind  !< index variable along z
    INTEGER(iwp) ::  j    !< index variable along y
    INTEGER(iwp) ::  k    !< index variable along z
    INTEGER(iwp) ::  l    !< grid level

    REAL(wp), DIMENSION(nzb:nzt_mg(glevel)+1) ::  tmp  !< odd-even sorted temporary array


    CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'start' )

    l = glevel
    ind_even_odd = even_odd_level(l)

    !$OMP PARALLEL PRIVATE (i,j,k,ind,tmp)
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE(tmp, ind) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl_mg(l)-1, nxr_mg(l)+1
       DO  j = nys_mg(l)-1, nyn_mg(l)+1

!
!--       Sort the data with even k index
          ind = nzb-1
          DO  k = nzb, nzt_mg(l), 2
             ind = ind + 1
             tmp(ind) = p_mg(k,j,i)
          ENDDO
!
!--       Sort the data with odd k index
          DO  k = nzb+1, nzt_mg(l)+1, 2
             ind = ind + 1
             tmp(ind) = p_mg(k,j,i)
          ENDDO

          p_mg(:,j,i) = tmp

       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL

    CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'stop' )

 END SUBROUTINE sort_k_to_even_odd_blocks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-dimension from sequential into blocks of even and odd. This is required to vectorize the
!> red-black subroutine. Version for 1D-REAL arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sort_k_to_even_odd_blocks_1d( f_mg, f_mg_b, glevel )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  glevel  !< grid level

    REAL(wp), DIMENSION(nzb+1:nzt_mg(glevel)) ::  f_mg    !< 1D input array
    REAL(wp), DIMENSION(nzb:nzt_mg(glevel)+1) ::  f_mg_b  !< 1D output array

    INTEGER(iwp) ::  ind  !< index variable along z
    INTEGER(iwp) ::  k    !< index variable along z


    ind = nzb - 1
!
!-- Sort the data with even k index.
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE(ind) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  k = nzb, nzt_mg(glevel), 2
       IF ( k >= nzb+1  .AND.  k <= nzt_mg(glevel) )  THEN
          ind = k / 2
          f_mg_b(ind) = f_mg(k)
       ENDIF
    ENDDO
    !$ACC END PARALLEL LOOP
!
!-- Sort the data with odd k index.
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE(ind) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  k = nzb+1, nzt_mg(glevel)+1, 2
       IF( k >= nzb+1  .AND.  k <= nzt_mg(glevel) )  THEN
          ind = (nzt_mg(glevel) / 2) + (k + 1) / 2
          f_mg_b(ind) = f_mg(k)
       ENDIF
    ENDDO
    !$ACC END PARALLEL LOOP

 END SUBROUTINE sort_k_to_even_odd_blocks_1d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-dimension from sequential into blocks of even and odd. This is required to vectorize the
!> red-black subroutine. Version for 3D-INTEGER arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sort_k_to_even_odd_blocks_int( i_mg , glevel )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  glevel  !< grid level

    INTEGER(iwp), DIMENSION(nzb:nzt_mg(glevel)+1,nys_mg(glevel)-1:nyn_mg(glevel)+1,                &
                            nxl_mg(glevel)-1:nxr_mg(glevel)+1) ::  i_mg   !< array to be sorted

    INTEGER(iwp) :: i        !< index variabel along x
    INTEGER(iwp) :: j        !< index variable along y
    INTEGER(iwp) :: k        !< index variable along z
    INTEGER(iwp) :: l        !< grid level
    INTEGER(iwp) :: ind      !< index variable along z
    INTEGER(iwp),DIMENSION(nzb:nzt_mg(glevel)+1) ::  tmp  !< temporary odd-even sorted array


    CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'start' )

    l = glevel
    ind_even_odd = even_odd_level(l)

    !$OMP PARALLEL PRIVATE (i,j,k,ind,tmp)
    !$OMP DO
    DO  i = nxl_mg(l)-1, nxr_mg(l)+1
       DO  j = nys_mg(l)-1, nyn_mg(l)+1
!
!--       Sort the data with even k index.
          ind = nzb-1
          DO  k = nzb, nzt_mg(l), 2
             ind = ind + 1
             tmp(ind) = i_mg(k,j,i)
          ENDDO
!
!--       Sort the data with odd k index.
          DO  k = nzb+1, nzt_mg(l)+1, 2
             ind = ind + 1
             tmp(ind) = i_mg(k,j,i)
          ENDDO

          i_mg(:,j,i) = tmp

       ENDDO
    ENDDO
    !$OMP END PARALLEL

    CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'stop' )

 END SUBROUTINE sort_k_to_even_odd_blocks_int


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-dimension from blocks of even and odd into sequential
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sort_k_to_sequential( p_mg )

    IMPLICIT NONE

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  p_mg  !< array to be sorted

    INTEGER(iwp) ::  i    !< index variable along x
    INTEGER(iwp) ::  j    !< index variable along y
    INTEGER(iwp) ::  k    !< index variable along z
    INTEGER(iwp) ::  l    !< grid level
    INTEGER(iwp) ::  ind  !< index variable along z

    REAL(wp),DIMENSION(nzb:nzt_mg(grid_level)+1) ::  tmp  !<


    l = grid_level

    !$OMP PARALLEL PRIVATE (i,j,k,ind,tmp)
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE(tmp, ind) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl_mg(l)-1, nxr_mg(l)+1
       DO  j = nys_mg(l)-1, nyn_mg(l)+1

          ind = nzb - 1
          tmp = p_mg(:,j,i)
          DO  k = nzb, nzt_mg(l), 2
             ind = ind + 1
             p_mg(k,j,i) = tmp(ind)
          ENDDO

          DO  k = nzb+1, nzt_mg(l)+1, 2
             ind = ind + 1
             p_mg(k,j,i) = tmp(ind)
          ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL

 END SUBROUTINE sort_k_to_sequential


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gather subdomain data from all PEs.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE mg_gather( f2, f2_sub )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  il      !<
    INTEGER(iwp) ::  ir      !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  jn      !<
    INTEGER(iwp) ::  js      !<
    INTEGER(iwp) ::  k       !<
    INTEGER(iwp) ::  nwords  !<

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f2    !<
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f2_l  !<

    REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,        &
                        mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  f2_sub  !<


    CALL cpu_log( log_point_s(34), 'mg_gather', 'start' )

    f2_l = 0.0_wp

!
!-- Store the local subdomain array on the total array. No ghost boundary values are stored
!-- because the internal ghost boundary values would enter twice the below MPI_ALLREDUCE sum.
    js = mg_loc_ind(3,myid)
    jn = mg_loc_ind(4,myid)
    il = mg_loc_ind(1,myid)
    ir = mg_loc_ind(2,myid)
    DO  i = il, ir
       DO  j = js, jn
          DO  k = nzb+1, nzt_mg(grid_level)
             f2_l(k,j,i) = f2_sub(k,j,i)
          ENDDO
       ENDDO
    ENDDO

!
!-- Bottom and top boundaries have not been set in case of Neumann BCs. Set them to zero here to
!-- avoid problems with MPI_SUM below.
    IF ( ibc_p_b == 1 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f2(nzb,:,:) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    IF ( ibc_p_t == 1 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f2(nzt_mg(grid_level)+1,:,:) = 0.0_wp
       !$ACC END KERNELS
    ENDIF
!
!-- Find out the number of array elements of the total array.
    nwords = SIZE( f2 )

!
!-- Gather subdomain data from all PEs.
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( f2_l(nzb,nys_mg(grid_level)-1,nxl_mg(grid_level)-1),                       &
                        f2(nzb,nys_mg(grid_level)-1,nxl_mg(grid_level)-1), nwords, MPI_REAL,       &
                        MPI_SUM, comm2d, ierr )

!
!-- Bottom and top boundaries must be set in case of Neumann BCs. These values are not used
!-- and set in the remaining parts of the multigrid solver because of the implicit treatment of
!-- BCs via flags.
    IF ( ibc_p_b == 1 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f2(nzb,:,:) = f2(nzb+1,:,:)
       !$ACC END KERNELS
    ENDIF

    IF ( ibc_p_t == 1 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f2(nzt_mg(grid_level)+1,:,:) = f2(nzt_mg(grid_level),:,:)
       !$ACC END KERNELS
    ENDIF

    CALL cpu_log( log_point_s(34), 'mg_gather', 'stop' )

 END SUBROUTINE mg_gather
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Scatter the subdomain data. Since all PEs did the calculations, MPI is not required to scatter
!> them.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE mg_scatter( p2, p2_sub )

    IMPLICIT NONE

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,  &
                        nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

    REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,        &
                        mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  p2_sub  !<


    CALL cpu_log( log_point_s(35), 'mg_scatter', 'start' )

    p2_sub = p2(:,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,                                       &
                  mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1)

    CALL cpu_log( log_point_s(35), 'mg_scatter', 'stop' )

 END SUBROUTINE mg_scatter
#endif

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This is where the multigrid technique takes place. V- and W- Cycle are implemented and steered by
!> the parameter "gamma_mg".
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE next_mg_level( f_mg, p_mg, p3, r )

    IMPLICIT NONE

    INTEGER(iwp) ::  i            !< index variable along x
    INTEGER(iwp) ::  j            !< index variable along y
    INTEGER(iwp) ::  k            !< index variable along z
    INTEGER(iwp) ::  nxl_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nxr_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nyn_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nys_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nzt_mg_save  !< to save index bound of the switch to pe0 level

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: f_mg  !<
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: p_mg  !<
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: p3    !<
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,        &
                        nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: r     !<

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,  &
                        nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  f2  !<
    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,  &
                        nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  f2_sub  !<

#if defined( __parallel )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  p2_sub  !<
#endif

   !$ACC DATA &
   !$ACC CREATE(f2, p2) IF(enable_openacc)


 10 IF ( grid_level == 1 )  THEN

!
!--    Solution on the coarsest grid. Double the number of Gauss-Seidel iterations in order to get a
!--    more accurate solution.
       ngsrb = 2 * ngsrb

       ind_even_odd = even_odd_level(grid_level)


       CALL redblack( f_mg, p_mg )

       ngsrb = ngsrb / 2


    ELSEIF ( grid_level /= 1 )  THEN

       grid_level_count(grid_level) = grid_level_count(grid_level) + 1

!
!--    Solution on the actual grid level.
       ind_even_odd = even_odd_level(grid_level)

       CALL redblack( f_mg, p_mg )

!
!--    Determination of the residual on this level.
       CALL resid( f_mg, p_mg, r )

!--    Restriction of the residual (finer grid values!) to the next coarser grid. Therefore, the
!--    grid level has to be decremented now. nxl..nzt have to be set to the coarse grid values,
!--    because these variables are needed for the exchange of ghost points in routine exchange_horiz.
       grid_level = grid_level - 1

       nxl = nxl_mg(grid_level)
       nys = nys_mg(grid_level)
       nxr = nxr_mg(grid_level)
       nyn = nyn_mg(grid_level)
       nzt = nzt_mg(grid_level)

       IF ( grid_level == mg_switch_to_pe0_level )  THEN
!
!--       From this level on, calculations are done on all PEs redundantly. First, carry out
!--       restriction on the subdomain. Therefore, indices of the level have to be changed to
!--       subdomain values in between (otherwise, the restrict routine would expect the gathered
!--       array).
          nxl_mg_save = nxl_mg(grid_level)
          nxr_mg_save = nxr_mg(grid_level)
          nys_mg_save = nys_mg(grid_level)
          nyn_mg_save = nyn_mg(grid_level)
          nzt_mg_save = nzt_mg(grid_level)
          nxl_mg(grid_level) = mg_loc_ind(1,myid)
          nxr_mg(grid_level) = mg_loc_ind(2,myid)
          nys_mg(grid_level) = mg_loc_ind(3,myid)
          nyn_mg(grid_level) = mg_loc_ind(4,myid)
          nzt_mg(grid_level) = mg_loc_ind(5,myid)
          nxl = mg_loc_ind(1,myid)
          nxr = mg_loc_ind(2,myid)
          nys = mg_loc_ind(3,myid)
          nyn = mg_loc_ind(4,myid)
          nzt = mg_loc_ind(5,myid)

          ALLOCATE( f2_sub(nzb:nzt_mg(grid_level)+1,nys_mg(grid_level)-1:nyn_mg(grid_level)+1,     &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) )

          CALL restrict( f2_sub, r )

!
!--       Restore the correct indices of this level.
          nxl_mg(grid_level) = nxl_mg_save
          nxr_mg(grid_level) = nxr_mg_save
          nys_mg(grid_level) = nys_mg_save
          nyn_mg(grid_level) = nyn_mg_save
          nzt_mg(grid_level) = nzt_mg_save
          nxl = nxl_mg(grid_level)
          nxr = nxr_mg(grid_level)
          nys = nys_mg(grid_level)
          nyn = nyn_mg(grid_level)
          nzt = nzt_mg(grid_level)
!
!--       Gather all arrays from the subdomains of all arrays. They will be redundantly gathered
!--       on all PEs.
#if defined( __parallel )
          CALL mg_gather( f2, f2_sub )
#endif

!
!--       Set switch for routine exchange_horiz, that no ghostpoint exchange has to be carried out
!--       from now on, because PEs contain the total domain.
          mg_switch_to_pe0 = .TRUE.

          DEALLOCATE( f2_sub )

       ELSE

          CALL restrict( f2, r )

       ENDIF

       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       p2 = 0.0_wp
       !$ACC END KERNELS

!
!--    Repeat the same procedure until the coarsest grid is reached.
       CALL next_mg_level( f2, p2, p3, r )

    ENDIF

!
!-- Now follows the prolongation.
    IF ( grid_level >= 2 )  THEN

!
!--    Prolongation of the new residual. The values are transferred from the coarse to the next
!--    finer grid.
       IF ( grid_level == mg_switch_to_pe0_level+1 )  THEN

#if defined( __parallel )
!
!--       At this level, the new residual first has to be scattered from PE0 to the other PEs.
          ALLOCATE( p2_sub(nzb:mg_loc_ind(5,myid)+1,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,     &
                                                    mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) )

          CALL mg_scatter( p2, p2_sub )

!
!--       Therefore, indices of the previous level have to be changed to subdomain values in between
!--       (otherwise, the prolong routine would expect the gathered array).
          nxl_mg_save = nxl_mg(grid_level-1)
          nxr_mg_save = nxr_mg(grid_level-1)
          nys_mg_save = nys_mg(grid_level-1)
          nyn_mg_save = nyn_mg(grid_level-1)
          nzt_mg_save = nzt_mg(grid_level-1)
          nxl_mg(grid_level-1) = mg_loc_ind(1,myid)
          nxr_mg(grid_level-1) = mg_loc_ind(2,myid)
          nys_mg(grid_level-1) = mg_loc_ind(3,myid)
          nyn_mg(grid_level-1) = mg_loc_ind(4,myid)
          nzt_mg(grid_level-1) = mg_loc_ind(5,myid)

!
!--       Set switch for routine exchange_horiz, that ghost point exchange has to be carried out
!--       again from now on
          mg_switch_to_pe0 = .FALSE.

          CALL prolong( p2_sub, p3 )

!
!--       Restore the correct indices of the previous level.
          nxl_mg(grid_level-1) = nxl_mg_save
          nxr_mg(grid_level-1) = nxr_mg_save
          nys_mg(grid_level-1) = nys_mg_save
          nyn_mg(grid_level-1) = nyn_mg_save
          nzt_mg(grid_level-1) = nzt_mg_save

          DEALLOCATE( p2_sub )
#endif

       ELSE

          CALL prolong( p2, p3 )

       ENDIF

!
!--    Computation of the new pressure correction. Therefore, values from prior grids are added up
!--    automatically stage by stage. Don't add the ghost point values, because they are not set in
!--    prolong, and not for p_mg in the coarser grid levels, too, because of implcit treatement of
!--    boundary conditions.
       !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
       !$ACC DEFAULT(PRESENT) IF(enable_openacc)
       DO  i = nxl_mg(grid_level), nxr_mg(grid_level)
          DO  j = nys_mg(grid_level), nyn_mg(grid_level)
             DO  k = nzb+1, nzt_mg(grid_level)
                p_mg(k,j,i) = p_mg(k,j,i) + p3(k,j,i)
             ENDDO
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP
!
!--    Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are implicitly
!--    treated via the flags array.
       CALL exchange_horiz( p_mg, 1, grid_level = grid_level, mg_switch_to_pe0 = mg_switch_to_pe0 )

!
!--    Relaxation of the new solution.
       CALL redblack( f_mg, p_mg )

    ENDIF


!
!-- The following few lines serve the steering of the multigrid scheme.
    IF ( grid_level == maximum_grid_level )  THEN

       GOTO 20

    ELSEIF ( grid_level /= maximum_grid_level  .AND.  grid_level /= 1  .AND.                       &
             grid_level_count(grid_level) /= gamma_mg )  THEN

       GOTO 10

    ENDIF

!
!-- Reset counter for the next call of poismg.
    grid_level_count(grid_level) = 0

!
!-- Continue with the next finer level. nxl..nzt have to be set to the finer grid values, because
!-- these variables are needed for the exchange of ghost points in routine exchange_horiz.
    grid_level = grid_level + 1
    ind_even_odd = even_odd_level(grid_level)

    nxl = nxl_mg(grid_level)
    nxr = nxr_mg(grid_level)
    nys = nys_mg(grid_level)
    nyn = nyn_mg(grid_level)
    nzt = nzt_mg(grid_level)

 20 CONTINUE

    !$ACC END DATA

 END SUBROUTINE next_mg_level


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initial settings for sorting k-dimension from sequential order (alternate even/odd) into blocks
!> of even and odd or vice versa.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_even_odd_blocks

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  l  !<

    LOGICAL, SAVE ::  lfirst = .TRUE.  !<


    IF ( .NOT. lfirst )  RETURN

    ALLOCATE( even_odd_level(maximum_grid_level) )

    ALLOCATE( f1_mg_b(nzb:nzt+1,maximum_grid_level), f2_mg_b(nzb:nzt+1,maximum_grid_level),        &
              f3_mg_b(nzb:nzt+1,maximum_grid_level), rho_air_mg_b(nzb:nzt+1,maximum_grid_level) )

    !$ACC ENTER DATA CREATE(f1_mg_b, f2_mg_b, f3_mg_b, rho_air_mg_b) IF(enable_openacc)

!
!-- Set border index between the even and odd block
    !$ACC PARALLEL LOOP GANG VECTOR COPYOUT(even_odd_level) DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = maximum_grid_level, 1, -1
       even_odd_level(i) = nzt_mg(i) / 2
    ENDDO
    !$ACC END PARALLEL LOOP

!
!-- Sort grid coefficients used in red/black scheme and for calculating the residual to block
!-- (even/odd) structure

    DO  l = maximum_grid_level, 1 , -1
       CALL sort_k_to_even_odd_blocks( f1_mg(nzb+1:nzt_mg(grid_level),l),                          &
                                       f1_mg_b(nzb:nzt_mg(grid_level)+1,l), l )
       CALL sort_k_to_even_odd_blocks( f2_mg(nzb+1:nzt_mg(grid_level),l),                          &
                                       f2_mg_b(nzb:nzt_mg(grid_level)+1,l), l )
       CALL sort_k_to_even_odd_blocks( f3_mg(nzb+1:nzt_mg(grid_level),l),                          &
                                       f3_mg_b(nzb:nzt_mg(grid_level)+1,l), l )
       CALL sort_k_to_even_odd_blocks( rho_air_mg(nzb+1:nzt_mg(grid_level),l),                     &
                                       rho_air_mg_b(nzb:nzt_mg(grid_level)+1,l), l )
    ENDDO

    lfirst = .FALSE.

 END SUBROUTINE init_even_odd_blocks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates wall flags for each grid level of the multigrid-solver.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poismg_init( even_odd_decomposition )

    USE arrays_3d,                                                                                 &
        ONLY:  dzu,                                                                                &
               dzw,                                                                                &
               rho_air,                                                                            &
               rho_air_zw

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t,                                                                            &
               masking_method

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nny,                                                                                &
               nx,                                                                                 &
               ny,                                                                                 &
               nz,                                                                                 &
               topo_flags
#if defined( __parallel )
    USE indices,                                                                                   &
        ONLY:  nbgp
#endif

    IMPLICIT NONE

    LOGICAL, INTENT(IN) ::  even_odd_decomposition  !< switch indicating if called from optimized or non-optimized version

#if defined( __parallel )
    INTEGER(iwp) ::  bufsize         !< size of buffer for sending/receiving contiguous data
    INTEGER(iwp) ::  nzb_l           !< lower index bound along z-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  stored_value    !< temporary variable
#endif
    INTEGER(iwp) ::  color           !< grid point color, either red (1) or black (2)
    INTEGER(iwp) ::  i               !< index variable along x
    INTEGER(iwp) ::  i_topo          !< i index for topo_flags (finest grid)
    INTEGER(iwp) ::  inc             !< incremental parameter for coarsening grid level
    INTEGER(iwp) ::  j               !< index variable along y
    INTEGER(iwp) ::  j_topo          !< j index for topo_flags (finest grid)
    INTEGER(iwp) ::  k               !< index variable along z
    INTEGER(iwp) ::  kbottom_uneven  !< bottom index where the values with (originally) uneven index start
    INTEGER(iwp) ::  k_topo          !< k index for topo_flags (finest grid)
    INTEGER(iwp) ::  l               !< loop variable indication current grid level
    INTEGER(iwp) ::  maximum_grid_level_l  !< maximum number of grid levels without switching to PE 0
    INTEGER(iwp) ::  mg_levels_x              !< maximum number of grid level allowed along x-direction
    INTEGER(iwp) ::  mg_levels_y              !< maximum number of grid level allowed along y-direction
    INTEGER(iwp) ::  mg_levels_z              !< maximum number of grid level allowed along z-direction
    INTEGER(iwp) ::  mg_switch_to_pe0_level_l  !< maximum number of grid level with switching to PE 0
    INTEGER(iwp) ::  ngp             !< number of grid points of topo_tmp array
    INTEGER(iwp) ::  num_wall        !< number of surrounding walls for a single grid point
    INTEGER(iwp) ::  nxl_l           !< index of left PE boundary for multigrid level
    INTEGER(iwp) ::  nxr_l           !< index of right PE boundary for multigrid level
    INTEGER(iwp) ::  nyn_l           !< index of north PE boundary for multigrid level
    INTEGER(iwp) ::  nys_l           !< index of south PE boundary for multigrid level
    INTEGER(iwp) ::  nzt_l           !< index of top PE boundary for multigrid level

    INTEGER(iwp) ::  nxl_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nxr_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nyn_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nys_mg_save  !< to save index bound of the switch to pe0 level
    INTEGER(iwp) ::  nzt_mg_save  !< to save index bound of the switch to pe0 level

!    INTEGER(iwp) ::  nxl_l_f  !< index of left PE boundary for multigrid level for next finer level
!    INTEGER(iwp) ::  nxr_l_f  !< index of right PE boundary for multigrid level for next finer level
!    INTEGER(iwp) ::  nyn_l_f  !< index of north PE boundary for multigrid level for next finer level
!    INTEGER(iwp) ::  nys_l_f  !< index of south PE boundary for multigrid level for next finer level
!    INTEGER(iwp) ::  nzt_l_f  !< index of top PE boundary for multigrid level for next finer level

#if defined( __parallel )
    INTEGER(iwp) ::  ind(5)  !< dummy array containing the subdomain bounds

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ind_all  !< dummy array containing index bounds on subdomain, used for gathering
#endif

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo_tmp         !< temporary array to store topography of the current grid level
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo_tmp_invers  !< invers temporary array to store topography of the current grid level
!    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo_tmp_f  !< temporary array to store topography of the next finer grid level

    REAL(wp) ::  dx_l  !< grid spacing along x on different multigrid levels
    REAL(wp) ::  dy_l  !< grid spacing along y on different multigrid levels


!
!-- Non-uniform subdomains are not allowed, because then subdomains would have different numbers of
!-- coarsening levels.
    IF ( non_uniform_subdomain )  THEN
       message_string = 'multigrid-solver does not allow to use non-uniform subdomains'
       CALL message( 'poismg_init', 'PAC0239', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( cycle_mg == 'w' )  THEN
       gamma_mg = 2
    ELSEIF ( cycle_mg == 'v' )  THEN
       gamma_mg = 1
    ELSE
       message_string = 'unknown multigrid cycle: cycle_mg = "' //  TRIM( cycle_mg ) // '"'
       CALL message( 'poismg_init', 'PAC0031', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Calculate number of allowed/possible grid levels as well as the gridpoint indices on each level.
!-- First calculate number of grid levels possible for the subdomains.
    mg_levels_x = 1
    mg_levels_y = 1
    mg_levels_z = 1

    i = nnx
    DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
       i = i / 2
       mg_levels_x = mg_levels_x + 1
    ENDDO

    j = nny
    DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
       j = j / 2
       mg_levels_y = mg_levels_y + 1
    ENDDO
!
!-- Do not use nnz because it might be > nz due to transposition requirements.
    k = nz
    DO WHILE ( MOD( k, 2 ) == 0  .AND.  k /= 2 )
       k = k / 2
       mg_levels_z = mg_levels_z + 1
    ENDDO
!
!-- The optimized MG-solver does not allow odd values for nz at the coarsest grid level.
    IF ( even_odd_decomposition )  THEN
       IF ( MOD( k, 2 ) /= 0 )  mg_levels_z = mg_levels_z - 1
!
!--    An odd value of nz does not work. The finest level must have an even value.
       IF (  mg_levels_z == 0 )  THEN
          message_string = 'optimized multigrid method requires nz to be even'
          CALL message( 'poisamg_init', 'PAC0241', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    maximum_grid_level = MIN( mg_levels_x, mg_levels_y, mg_levels_z )
    maximum_grid_level_default = maximum_grid_level
!
!-- Restrict the number of multigrid levels if forced by user.
    IF ( max_mg_grid_levels /= 9999 )  THEN
       IF ( maximum_grid_level > max_mg_grid_levels )  THEN
          maximum_grid_level = max_mg_grid_levels
       ENDIF
    ENDIF
!
!-- Check if subdomain sizes prevents any coarsening.
!-- This case, the maximum number of grid levels is 1, i.e. effectively a Gauss-Seidel scheme is
!-- applied rather than a multigrid approach. Give a warning in this case.
    IF ( maximum_grid_level == 1  .AND.  mg_switch_to_pe0_level == -1 )  THEN
       message_string = 'no grid coarsening possible, multigrid ' //                            &
                        'approach effectively reduces to a Gauss-Seidel scheme'
       CALL message( 'poismg_init', 'PAC0242', 0, 1, 0, 6, 0 )
    ENDIF

!
!-- Find out, if the total domain allows more levels. These additional levels are identically
!-- processed on all PEs.
    IF ( numprocs > 1  .AND.  mg_switch_to_pe0_level /= -1 )  THEN

       IF ( mg_levels_z > MIN( mg_levels_x, mg_levels_y ) )  THEN

          mg_switch_to_pe0_level_l = maximum_grid_level

          mg_levels_x = 1
          mg_levels_y = 1

          i = nx+1
          DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
             i = i / 2
             mg_levels_x = mg_levels_x + 1
          ENDDO

          j = ny+1
          DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
             j = j / 2
             mg_levels_y = mg_levels_y + 1
          ENDDO

          maximum_grid_level_l = MIN( mg_levels_x, mg_levels_y, mg_levels_z )
          maximum_grid_level_default = maximum_grid_level_l
!
!--       Restrict the number of multigrid levels if forced by user.
          IF ( max_mg_grid_levels /= 9999 )  THEN
             IF ( maximum_grid_level_l > max_mg_grid_levels )  THEN
                maximum_grid_level_l = max_mg_grid_levels
             ENDIF
          ENDIF

          IF ( maximum_grid_level_l > mg_switch_to_pe0_level_l )  THEN
             mg_switch_to_pe0_level_l = maximum_grid_level_l - mg_switch_to_pe0_level_l + 1
          ELSE
             mg_switch_to_pe0_level_l = 0
          ENDIF

       ELSE

          mg_switch_to_pe0_level_l = 0
          maximum_grid_level_l = maximum_grid_level

       ENDIF

!
!--    Use switch level calculated above only if it is not pre-defined by user.
       IF ( mg_switch_to_pe0_level == 0 )  THEN

          IF ( mg_switch_to_pe0_level_l /= 0 )  THEN
             mg_switch_to_pe0_level = mg_switch_to_pe0_level_l
             maximum_grid_level     = maximum_grid_level_l
          ENDIF

       ELSE
!
!--       Check pre-defined value and reset to default, if neccessary
          IF ( mg_switch_to_pe0_level < mg_switch_to_pe0_level_l  .OR.                             &
               mg_switch_to_pe0_level >= maximum_grid_level_l )  THEN
             message_string = 'mg_switch_to_pe0_level out of range and reset to 0'
             CALL message( 'poismg_init', 'PAC0243', 0, 1, 0, 6, 0 )
             mg_switch_to_pe0_level = 0
          ELSE
!
!--          Use the largest number of possible levels anyway and recalculate the switch level to
!--          this largest number of possible values
             maximum_grid_level = maximum_grid_level_l

          ENDIF

       ENDIF

    ENDIF

    ALLOCATE( grid_level_count(maximum_grid_level),                                             &
              nxl_mg(0:maximum_grid_level), nxr_mg(0:maximum_grid_level),                       &
              nyn_mg(0:maximum_grid_level), nys_mg(0:maximum_grid_level),                       &
              nzt_mg(0:maximum_grid_level) )

    grid_level_count = 0
!
!-- Index zero required as dummy due to definition of arrays f2 and p2 in recursive subroutine
!-- next_mg_level
    nxl_mg(0) = 0; nxr_mg(0) = 0; nyn_mg(0) = 0; nys_mg(0) = 0; nzt_mg(0) = 0

    nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzt_l = nzt

    DO  i = maximum_grid_level, 1 , -1

       IF ( i == mg_switch_to_pe0_level )  THEN
#if defined( __parallel )
!
!--       Save the grid size of the subdomain at the switch level, because it is needed in poismg.
          ind(1) = nxl_l; ind(2) = nxr_l
          ind(3) = nys_l; ind(4) = nyn_l
          ind(5) = nzt_l

          ALLOCATE( ind_all(5*numprocs), mg_loc_ind(5,0:numprocs-1) )
          CALL MPI_ALLGATHER( ind, 5, MPI_INTEGER, ind_all, 5, MPI_INTEGER, comm2d, ierr )

          DO  j = 0, numprocs-1
             DO  k = 1, 5
                mg_loc_ind(k,j) = ind_all(k+j*5)
             ENDDO
          ENDDO

          DEALLOCATE( ind_all )
!
!--       Calculate the grid size of the total domain.
          nxr_l = ( nxr_l-nxl_l+1 ) * npex - 1
          nxl_l = 0
          nyn_l = ( nyn_l-nys_l+1 ) * npey - 1
          nys_l = 0
!
!--       The size of this gathered array must not be larger than the array tend, which is used
!--       in the multigrid scheme as a temporary array. Therefore the subdomain size of an PE is
!--       calculated and the size of the gathered grid. These values are used in routines pres
!--       and poismg.
          subdomain_size = ( nxr - nxl + 2 * nbgp + 1 ) *                                          &
                           ( nyn - nys + 2 * nbgp + 1 ) * ( nzt - nzb + 2 )
          gathered_size  = ( nxr_l - nxl_l + 3 ) * ( nyn_l - nys_l + 3 ) * ( nzt_l - nzb + 2 )

#else
          message_string = 'multigrid gather/scatter impossible in non parallel mode'
          CALL message( 'poismg_init', 'PAC0244', 1, 2, 0, 6, 0 )
#endif
       ENDIF

       nxl_mg(i) = nxl_l
       nxr_mg(i) = nxr_l
       nys_mg(i) = nys_l
       nyn_mg(i) = nyn_l
       nzt_mg(i) = nzt_l

       nxl_l = nxl_l / 2
       nxr_l = nxr_l / 2
       nys_l = nys_l / 2
       nyn_l = nyn_l / 2
       nzt_l = nzt_l / 2

    ENDDO

!
!-- Allocate array for the number of filtered holes on the respective grid level. In case of
!-- masking method, filtering is not required and the number will stay zero..
    ALLOCATE( poismg_filtered_holes(maximum_grid_level) )
    poismg_filtered_holes = 0

!
!-- Temporary problem: Currently calculation of maxerror in routine poismg crashes if grid data
!-- are collected on PE0 already on the finest grid level.
!-- To be solved later.
    IF ( maximum_grid_level == mg_switch_to_pe0_level )  THEN
       message_string = 'grid coarsening on subdomain level cannot be performed'
       CALL message( 'poismg_init', 'PAC0245', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Compute grid spacings s and grid factors for the grid levels with respective density on each
!-- grid.
    ALLOCATE( ddx2_mg(maximum_grid_level) )
    ALLOCATE( ddy2_mg(maximum_grid_level) )
    ALLOCATE( dzu_mg(nzb+1:nzt+1,maximum_grid_level) )
    ALLOCATE( dzw_mg(nzb+1:nzt+1,maximum_grid_level) )
    ALLOCATE( f1_mg(nzb+1:nzt,maximum_grid_level) )
    ALLOCATE( f2_mg(nzb+1:nzt,maximum_grid_level) )
    ALLOCATE( f3_mg(nzb+1:nzt,maximum_grid_level) )
    ALLOCATE( rho_air_mg(nzb:nzt+1,maximum_grid_level) )
    ALLOCATE( rho_air_zw_mg(nzb:nzt+1,maximum_grid_level) )

    dzu_mg(:,maximum_grid_level) = dzu
    rho_air_mg(:,maximum_grid_level) = rho_air
!
!-- Next line to ensure an equally spaced grid.
    dzu_mg(1,maximum_grid_level) = dzu(2)
    rho_air_mg(nzb,maximum_grid_level) = rho_air(nzb) + (rho_air(nzb) - rho_air(nzb+1))

    dzw_mg(:,maximum_grid_level) = dzw
    rho_air_zw_mg(:,maximum_grid_level) = rho_air_zw
    nzt_l = nzt
    DO  l = maximum_grid_level-1, 1, -1
       dzu_mg(nzb+1,l) = 2.0_wp * dzu_mg(nzb+1,l+1)
       dzw_mg(nzb+1,l) = 2.0_wp * dzw_mg(nzb+1,l+1)
       rho_air_mg(nzb,l)    = rho_air_mg(nzb,l+1)    + ( rho_air_mg(nzb,l+1)    -              &
                                                         rho_air_mg(nzb+1,l+1)    )
       rho_air_zw_mg(nzb,l) = rho_air_zw_mg(nzb,l+1) + ( rho_air_zw_mg(nzb,l+1) -              &
                                                         rho_air_zw_mg(nzb+1,l+1) )
       rho_air_mg(nzb+1,l)    = rho_air_mg(nzb+1,l+1)
       rho_air_zw_mg(nzb+1,l) = rho_air_zw_mg(nzb+1,l+1)
       nzt_l = nzt_l / 2
       DO  k = 2, nzt_l+1
          dzu_mg(k,l) = dzu_mg(2*k-2,l+1) + dzu_mg(2*k-1,l+1)
          dzw_mg(k,l) = dzw_mg(2*k-2,l+1) + dzw_mg(2*k-1,l+1)
          rho_air_mg(k,l)    = rho_air_mg(2*k-1,l+1)
          rho_air_zw_mg(k,l) = rho_air_zw_mg(2*k-1,l+1)
       ENDDO
    ENDDO

    nzt_l = nzt
    dx_l  = dx
    dy_l  = dy
    DO  l = maximum_grid_level, 1, -1
       ddx2_mg(l) = 1.0_wp / dx_l**2
       ddy2_mg(l) = 1.0_wp / dy_l**2
       DO  k = nzb+1, nzt_l
          f2_mg(k,l) = rho_air_zw_mg(k,l)   / ( dzu_mg(k+1,l) * dzw_mg(k,l) )
          f3_mg(k,l) = rho_air_zw_mg(k-1,l) / ( dzu_mg(k,l)   * dzw_mg(k,l) )
          f1_mg(k,l) = 2.0_wp * ( ddx2_mg(l) + ddy2_mg(l) ) * rho_air_mg(k,l) +                    &
                       f2_mg(k,l) + f3_mg(k,l)
       ENDDO
       nzt_l = nzt_l / 2
       dx_l  = dx_l * 2.0_wp
       dy_l  = dy_l * 2.0_wp
    ENDDO

!
!-- Copy data to the device, that is only required within the solver.
    !$ACC ENTER DATA &
    !$ACC COPYIN(f1_mg, f2_mg, f3_mg, rho_air_mg) &
    !$ACC COPYIN(nxl_mg, nxr_mg, nys_mg, nyn_mg, nzt_mg) &
    !$ACC COPYIN(ddx2_mg, ddy2_mg) IF(enable_openacc)

#if defined( __parallel )
!
!-- Definition of MPI-derived datatypes coarser level grids.
!-- First re-allocate the respective arrays for the data types. So far, the size is (0:0) and
!-- variables contain only data of the original model grid.
    stored_value = ngp_xz(0)
    DEALLOCATE( ngp_xz )
    ALLOCATE( ngp_xz(0:maximum_grid_level) )
    ngp_xz(0) = stored_value

    stored_value = ngp_xz_int(0)
    DEALLOCATE( ngp_xz_int )
    ALLOCATE( ngp_xz_int(0:maximum_grid_level) )
    ngp_xz_int(0) = stored_value

    stored_value = ngp_yz(0)
    DEALLOCATE( ngp_yz )
    ALLOCATE( ngp_yz(0:maximum_grid_level) )
    ngp_yz(0) = stored_value

    stored_value = ngp_yz_int(0)
    DEALLOCATE( ngp_yz_int )
    ALLOCATE( ngp_yz_int(0:maximum_grid_level) )
    ngp_yz_int(0) = stored_value

    stored_value = type_xz(0)
    DEALLOCATE( type_xz )
    ALLOCATE( type_xz(0:maximum_grid_level) )
    type_xz(0) = stored_value

    stored_value = type_xz_int(0)
    DEALLOCATE( type_xz_int )
    ALLOCATE( type_xz_int(0:maximum_grid_level) )
    type_xz_int(0) = stored_value

    stored_value = type_yz(0)
    DEALLOCATE( type_yz )
    ALLOCATE( type_yz(0:maximum_grid_level) )
    type_yz(0) = stored_value

    stored_value = type_yz_int(0)
    DEALLOCATE( type_yz_int )
    ALLOCATE( type_yz_int(0:maximum_grid_level) )
    type_yz_int(0) = stored_value

!
!-- Definition of MPI-datatyoe using 1 ghost layer only.
    nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzb_l = nzb; nzt_l = nzt

    DO  i = maximum_grid_level, 1 , -1

       ngp_xz(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
       ngp_yz(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)

       ngp_xz_int(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
       ngp_yz_int(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)
!
!--    MPI data type for REAL arrays (xz-layers).
       CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz(i), MPI_REAL, type_xz(i), ierr )
       CALL MPI_TYPE_COMMIT( type_xz(i), ierr )

!
!--    MPI data type for INTEGER arrays (xz-layers).
       CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz_int(i), MPI_INTEGER,             &
                             type_xz_int(i), ierr )
       CALL MPI_TYPE_COMMIT( type_xz_int(i), ierr )

!
!--    MPI data type for REAL arrays (yz-layers).
       CALL MPI_TYPE_VECTOR( 1, ngp_yz(i), ngp_yz(i), MPI_REAL, type_yz(i), ierr )
       CALL MPI_TYPE_COMMIT( type_yz(i), ierr )
!
!--    MPI data type for INTEGER arrays (yz-layers).
       CALL MPI_TYPE_VECTOR( 1, ngp_yz_int(i), ngp_yz_int(i), MPI_INTEGER, type_yz_int(i), ierr )
       CALL MPI_TYPE_COMMIT( type_yz_int(i), ierr )

       nxl_l = nxl_l / 2
       nxr_l = nxr_l / 2
       nys_l = nys_l / 2
       nyn_l = nyn_l / 2
       nzt_l = nzt_l / 2

    ENDDO

#endif

!
!-- Allocate arrays containing the grid level masking flags.
    IF ( mg_switch_to_pe0_level > 0 )  THEN
!
!--    Level 0 contains flags required at the switch to PE0 level
       ALLOCATE( gl(0:maximum_grid_level) )
    ELSE
       ALLOCATE( gl(1:maximum_grid_level) )
    ENDIF

    DO  l = 1, maximum_grid_level
       ALLOCATE( gl(l)%flags(nzb:nzt_mg(l)+1,nys_mg(l)-1:nyn_mg(l)+1,nxl_mg(l)-1:nxr_mg(l)+1) )
       gl(l)%flags = 0
       IF ( l == mg_switch_to_pe0_level )  THEN
          ALLOCATE( gl(0)%flags(nzb:mg_loc_ind(5,myid)+1,                                          &
                                mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,                         &
                                mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) )
          gl(0)%flags = 0
       ENDIF
    ENDDO

!
!-- Initial settings for sorting k-dimension from sequential order (alternate even/odd) into blocks
!-- of even and odd or vice versa (not required for the noopt-version).
    IF ( even_odd_decomposition )  THEN
       grid_level = maximum_grid_level
       CALL init_even_odd_blocks
    ENDIF

!
!-- Grid point increment of the current level.
    inc = 1
    DO  l = maximum_grid_level, 1 , -1
!
!--    Set grid_level as required for exchange_horiz_int.
       grid_level = l

       nxl_l = nxl_mg(l)
       nxr_l = nxr_mg(l)
       nys_l = nys_mg(l)
       nyn_l = nyn_mg(l)
       nzt_l = nzt_mg(l)
!
!--    Set switch for routine exchange_horiz, that no ghostpoint exchange has to be carried out
!--    for this level because PEs contain the total domain.
       IF ( l <= mg_switch_to_pe0_level )  mg_switch_to_pe0 = .TRUE.

!
!--    Depending on the grid level, set the respective bits in case of neighbouring walls.
!--    Bit 0:  wall to the bottom
!--    Bit 1:  wall to the top
!--    Bit 2:  wall to the south
!--    Bit 3:  wall to the north
!--    Bit 4:  wall to the left
!--    Bit 5:  wall to the right
!--    Bit 6:  inside building

!
!--    Allocate temporary array for topography heights on coarser grid level. Initialize all
!--    elements as fluid. This is the setting for the masking method, where the multigrid solver
!--    works like the FFT-solver, it "runs through the topography grid points" and is aware of the
!--    topography only via the divergence, which has been set zero for these points.
       ALLOCATE( topo_tmp(nzb:nzt_l+1,nys_l-1:nyn_l+1,nxl_l-1:nxr_l+1) )
       topo_tmp(:,:,:) = IBSET( topo_tmp(:,:,:), 0 )

!
!--    Set the topography grid points.
       IF ( .NOT. masking_method )  THEN

          DO  i = nxl_l-1, nxr_l+1
             DO  j = nys_l-1, nyn_l+1
                DO  k = nzb, nzt_l+1
                   k_topo = k * inc
                   j_topo = j * inc
                   i_topo = i * inc
!
!--                When levels are equal or below the switch to PE0 level, topo_tmp contains the
!--                total domain, but topo_flags (always) contains only the grid points of the
!--                subdomain, and only they can be stored on topo_tmp.
                   IF ( k_topo >= nzb   .AND.  k_topo <= nzt+1  .AND.                              &
                        j_topo >= nysg  .AND.  j_topo <= nyng   .AND.                              &
                        i_topo >= nxlg  .AND.  i_topo <= nxrg )                                    &
                   THEN
                      topo_tmp(k,j,i) = topo_flags(k_topo,j_topo,i_topo)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
!
!--       For grid levels containing the total domain, collect data from the subdomains.
          IF ( l <= mg_switch_to_pe0_level )  THEN
!
!--          Using MPI_SUM requires that topography points are marked with 1 and atmosphere
!--          with 0. Use a temporary INTEGER array for this
             ALLOCATE( topo_tmp_invers(nzb:nzt_l+1,nys_l-1:nyn_l+1,nxl_l-1:nxr_l+1) )
             DO  i = nxl_l-1, nxr_l+1
                DO  j = nys_l-1, nyn_l+1
                   DO  k = nzb, nzt_l+1
                      IF ( .NOT. BTEST( topo_tmp(k,j,i), 0 ) )  THEN
                         topo_tmp_invers(k,j,i) = 1
                      ELSE
                         topo_tmp_invers(k,j,i) = 0
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ngp = ( nzt_l - nzb + 2 ) * ( nyn_l - nys_l + 3 ) * ( nxr_l - nxl_l + 3 )
#if defined( __parallel )
             CALL MPI_ALLREDUCE( MPI_IN_PLACE, topo_tmp_invers(nzb,nys_l-1,nxl_l-1), ngp,          &
                                 MPI_INTEGER, MPI_SUM, comm2d, ierr )
#endif
             DO  i = nxl_l-1, nxr_l+1
                DO  j = nys_l-1, nyn_l+1
                   DO  k = nzb, nzt_l+1
                      IF ( topo_tmp_invers(k,j,i) == 0 )  THEN
!
!--                      Set atmosphere.
                         topo_tmp(k,j,i) = IBSET( topo_tmp(k,j,i), 0 )
                      ELSE
!
!--                      Set building/wall.
                         topo_tmp(k,j,i) = IBCLR( topo_tmp(k,j,i), 0 )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO

             DEALLOCATE( topo_tmp_invers )

          ENDIF

!
!--       Filter holes that appear in coarser levels.
          IF ( l < maximum_grid_level )  THEN
             DO  i = nxl_l, nxr_l
                DO  j = nys_l, nyn_l
                   DO  k = nzb+1, nzt_l
                      IF ( BTEST( topo_tmp(k,j,i), 0 ) )  THEN
                         num_wall = 0
                         IF ( .NOT. BTEST( topo_tmp(k,j-1,i), 0 ) )  num_wall = num_wall + 1
                         IF ( .NOT. BTEST( topo_tmp(k,j+1,i), 0 ) )  num_wall = num_wall + 1
                         IF ( .NOT. BTEST( topo_tmp(k,j,i-1), 0 ) )  num_wall = num_wall + 1
                         IF ( .NOT. BTEST( topo_tmp(k,j,i+1), 0 ) )  num_wall = num_wall + 1
                         IF ( .NOT. BTEST( topo_tmp(k-1,j,i), 0 ) )  num_wall = num_wall + 1
                         IF ( .NOT. BTEST( topo_tmp(k+1,j,i), 0 ) )  num_wall = num_wall + 1

                         IF ( num_wall >= 4 )  THEN
                            poismg_filtered_holes(l) = poismg_filtered_holes(l) + 1
!
!--                         Set building wall at this point.
                            topo_tmp(k,j,i) = IBCLR( topo_tmp(k,j,i), 0 )
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
#if defined( __parallel )
             CALL MPI_ALLREDUCE( MPI_IN_PLACE, poismg_filtered_holes(l), 1, MPI_INTEGER, MPI_SUM,  &
                                 comm2d, ierr )
#endif
          ENDIF

          topo_tmp(nzt_l+1,:,:) = topo_tmp(nzt_l,:,:)

       ENDIF

!
!--    In case that Neumann-BCs at the bottom (nzb) and/or top (nzt+1) are set, consider them
!--    via setting a wall at these points. This way, the BCs are used via the implicit
!--    Neumann-BCs that the multigrid-solvers assumes at walls.
       IF ( ibc_p_b == 1 )  topo_tmp(nzb,:,:)     = IBCLR( topo_tmp(nzb,:,:), 0 )
       IF ( ibc_p_t == 1 )  topo_tmp(nzt_l+1,:,:) = IBCLR( topo_tmp(nzt_l+1,:,:), 0 )
!
!--    Exchange ghost points on respective multigrid level.
#if defined( __parallel )
       CALL exchange_horiz_int( topo_tmp, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1, type_xz_int(l),    &
                                type_yz_int(l) )
#else
       CALL exchange_horiz_int( topo_tmp, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1 )
#endif
!
!--    Set walls at the total domain boundaries in case of non-cyclic boundary conditions
!--    to consider Neumann-BCs for pressure at these boundaries. Levels equal to or below the
!--    switch to pe0 level contain the total domain, so all cores have to set the non-cylic
!--    conditions.
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s  .OR.  mg_switch_to_pe0 )  THEN
             topo_tmp(:,-1,:) = IBCLR( topo_tmp(:,-1,:), 0 )
          ENDIF
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n  .OR.  mg_switch_to_pe0 )  THEN
             topo_tmp(:,nyn_l+1,:) = IBCLR( topo_tmp(:,nyn_l+1,:), 0 )
          ENDIF
       ENDIF
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l  .OR.  mg_switch_to_pe0 )  THEN
             topo_tmp(:,:,-1) = IBCLR( topo_tmp(:,:,-1), 0 )
          ENDIF
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r  .OR.  mg_switch_to_pe0 )  THEN
             topo_tmp(:,:,nxr_l+1) = IBCLR( topo_tmp(:,:,nxr_l+1), 0 )
          ENDIF
       ENDIF

!
!--    Now set the flags, based on the topography/wall settings above.
       DO  i = nxl_l, nxr_l
          DO  j = nys_l, nyn_l
             DO  k = nzb, nzt_l
!
!--             Inside/outside building (inside building does not need further tests for walls).
                IF ( .NOT. BTEST( topo_tmp(k,j,i), 0 ) )  THEN

                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 6 )

                ELSE
!
!--                Bottom wall.
                   IF ( .NOT. BTEST( topo_tmp(k-1,j,i), 0 ) )  THEN
                      gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 0 )
                   ENDIF
!
!--                Top wall.
                   IF ( .NOT. BTEST( topo_tmp(k+1,j,i), 0 ) )  THEN
                      gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 1 )
                   ENDIF
!
!--                South wall.
                   IF ( .NOT. BTEST( topo_tmp(k,j-1,i), 0 ) )  THEN
                      gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 2 )
                   ENDIF
!
!--                North wall.
                   IF ( .NOT. BTEST( topo_tmp(k,j+1,i), 0 ) )  THEN
                      gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 3 )
                   ENDIF
!
!--                Left wall.
                   IF ( .NOT. BTEST( topo_tmp(k,j,i-1), 0 ) )  THEN
                      gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 4 )
                   ENDIF
!
!--                Right wall.
                   IF ( .NOT. BTEST( topo_tmp(k,j,i+1), 0 ) )  THEN
                      gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 5 )
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO

       gl(l)%flags(nzt_l+1,:,:) = gl(l)%flags(nzt_l,:,:)

#if defined( __parallel )
       CALL exchange_horiz_int( gl(l)%flags, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1, type_xz_int(l), &
                                type_yz_int(l) )
#else
       CALL exchange_horiz_int( gl(l)%flags, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1 )
#endif

!
!--    Set non-cyclic boundary conditions. The ghost layers of the total domain have to be
!--    set as walls (flag 6) to consider Neumann BCs at these boundaries. Flag 6 is e.g. used in
!--    restrict. Levels below or equal the switch to pe0 level contain the total domain, so
!--    if non-cyclic conditions are used flags must always be set for them.
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s  .OR.  mg_switch_to_pe0 )  THEN
             gl(l)%flags(:,-1,:) = gl(l)%flags(:,0,:)
             gl(l)%flags(:,-1,:) = IBSET( gl(l)%flags(:,-1,:), 6 )
          ENDIF
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n  .OR.  mg_switch_to_pe0  )  THEN
             gl(l)%flags(:,nyn_l+1,:) = gl(l)%flags(:,nyn_l,:)
             gl(l)%flags(:,nyn_l+1,:) = IBSET( gl(l)%flags(:,nyn_l+1,:), 6 )
          ENDIF
       ENDIF
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l  .OR.  mg_switch_to_pe0  )  THEN
             gl(l)%flags(:,:,-1) = gl(l)%flags(:,:,0)
             gl(l)%flags(:,:,-1) = IBSET( gl(l)%flags(:,:,-1), 6 )
          ENDIF
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r  .OR.  mg_switch_to_pe0  )  THEN
             gl(l)%flags(:,:,nxr_l+1) = gl(l)%flags(:,:,nxr_l)
             gl(l)%flags(:,:,nxr_l+1) = IBSET( gl(l)%flags(:,:,nxr_l+1), 6 )
          ENDIF
       ENDIF

!
!--    Set special flags to be used in routine prolong. They are not required at ghost layers,
!--    so exchange_horiz does not have to be applied.
!--    Flags indicate for fluid points (k,j,i) in the respective finer level, that the area is
!--    completely covered by buildings/topography in the coarser level in the respective
!--    direction(s), so that interpolation from the adjacent coarse grid points would give wrong
!--    results. Different bits (7-13) are required and set depending on the position of
!--    point (k,j,i) with respect to the coarse grid cube (i.e. on one of the edges, the center of
!--    a plane, or the center of the cube).
       DO  i = nxl_l, nxr_l
          DO  j = nys_l, nyn_l
             DO  k = nzb, nzt_l
!
!--             Walls to the left and right.
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )  .AND.  BTEST( gl(l)%flags(k,j,i+1), 6 ) )  THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 7 )
                ENDIF
!
!--             Walls to the south and north.
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )  .AND.  BTEST( gl(l)%flags(k,j+1,i), 6 ) )  THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 8 )
                ENDIF
!
!--             Walls below and above.
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )  .AND.  BTEST( gl(l)%flags(k+1,j,i), 6 ) )  THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 9 )
                ENDIF
!
!--             Walls around the center of the planes stretched by four points of the
!--             coarse grid cube (xy-plane).
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )  .AND.  BTEST( gl(l)%flags(k,j,i+1), 6 )  .AND.&
                     BTEST( gl(l)%flags(k,j+1,i), 6 ) )                                            &
                THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 10 )
                ENDIF
!
!--             Walls around the center of the planes stretched by four points of the
!--             coarse grid cube (xz-plane).
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )  .AND.  BTEST( gl(l)%flags(k,j,i+1), 6 )  .AND.&
                     BTEST( gl(l)%flags(k+1,j,i), 6 ) )                                            &
                THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 11 )
                ENDIF
!
!--             Walls around the center of the planes stretched by four points of the
!--             coarse grid cube (yz-plane).
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )  .AND.  BTEST( gl(l)%flags(k,j+1,i), 6 )  .AND.&
                     BTEST( gl(l)%flags(k+1,j,i), 6 ) )                                            &
                THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 12 )
                ENDIF
!
!--             Walls around the middle of coarse grid cube.
                IF ( BTEST( gl(l)%flags(k,j,i), 6 )    .AND.  BTEST( gl(l)%flags(k,j,i+1), 6 )  .AND. &
                     BTEST( gl(l)%flags(k,j+1,i), 6 )  .AND.  BTEST( gl(l)%flags(k+1,j,i), 6 ) )   &
                THEN
                   gl(l)%flags(k,j,i) = IBSET( gl(l)%flags(k,j,i), 13 )
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!
!--    Consider walls/topography that may be only one grid point wide on the finer level. On the
!--    current level they may disappear, if these walls are located in between two neighboring
!--    points of the current level.
!       IF ( l < maximum_grid_level )  THEN
!!
!!--       One point wide wall along x:
!          DO  i = nxl_l, nxr_l
!             DO  j = nys_l, nyn_l
!                DO  k = nzb, nzt_l
!!
!!--                Wall to the bottom. There is no point below nzb!
!                   IF ( k > nzb )  THEN
!                      IF ( .NOT. BTEST( topo_tmp_f(k*2-1,j*2,i*2), 0 )  .AND.                      &
!                                 BTEST( topo_tmp_f(k*2,j*2,i*2),   0 ) )                           &
!                      THEN
!                         flags(k,j,i) = IBSET( flags(k,j,i), 0 )
!                      ENDIF
!                   ENDIF
!!
!!--                Wall to the top.
!                   IF ( .NOT. BTEST( topo_tmp_f(k*2+1,j*2,i*2), 0 )  .AND.                         &
!                              BTEST( topo_tmp_f(k*2,j*2,i*2),   0 ) )                              &
!                   THEN
!                      flags(k,j,i) = IBSET( flags(k,j,i), 1 )
!                   ENDIF
!!
!!--                Wall to the south.
!                   IF ( .NOT. BTEST( topo_tmp_f(k*2,j*2-1,i*2), 0 )  .AND.                         &
!                              BTEST( topo_tmp_f(k*2,j*2,i*2),   0 ) )                              &
!                   THEN
!                      flags(k,j,i) = IBSET( flags(k,j,i), 2 )
!                   ENDIF
!!
!!--                Wall to the north.
!                   IF ( .NOT. BTEST( topo_tmp_f(k*2,j*2+1,i*2), 0 )  .AND.                         &
!                              BTEST( topo_tmp_f(k*2,j*2,i*2),   0 ) )                              &
!                   THEN
!                      flags(k,j,i) = IBSET( flags(k,j,i), 3 )
!                   ENDIF
!!
!!--                Wall to the left.
!                   IF ( .NOT. BTEST( topo_tmp_f(k*2,j*2,i*2-1), 0 )  .AND.                         &
!                              BTEST( topo_tmp_f(k*2,j*2,i*2),   0 ) )                              &
!                   THEN
!                      flags(k,j,i) = IBSET( flags(k,j,i), 4 )
!                   ENDIF
!!
!!--                Wall to the right.
!                   IF ( .NOT. BTEST( topo_tmp_f(k*2,j*2,i*2+1), 0 )  .AND.                         &
!                              BTEST( topo_tmp_f(k*2,j*2,i*2),   0 ) )                              &
!                   THEN
!                      flags(k,j,i) = IBSET( flags(k,j,i), 5 )
!                   ENDIF
!
!                ENDDO
!             ENDDO
!          ENDDO
!
!       ENDIF
!
!
!--    Save the topography for the next level, where it is the finer level (f) with respect to
!--    the current level, and save the respective index bounds, too.
!       IF ( ALLOCATED( topo_tmp_f ) )  DEALLOCATE( topo_tmp_f )
!       ALLOCATE( topo_tmp_f(nzb:nzt_l+1,nys_l-1:nyn_l+1,nxl_l-1:nxr_l+1) )
!       topo_tmp_f = topo_tmp
!       nxl_l_f = nxl_l
!       nxr_l_f = nxr_l
!       nys_l_f = nys_l
!       nyn_l_f = nyn_l
!       nzt_l_f = nzt_l

       DEALLOCATE( topo_tmp )
!
!--    Sort flags to even/odd (not required for the noopt-version).
       IF ( even_odd_decomposition )  CALL sort_k_to_even_odd_blocks( gl(l)%flags, l )

       mg_switch_to_pe0 = .FALSE.

!
!--    For the switch to PE0 level a flag array is required in restrict for the subdomains, too.
       IF ( l == mg_switch_to_pe0_level )  THEN
!
!--       Set grid_level as required for exchange_horiz_int.
          grid_level = l
!
!--       Indices of the level have to be changed to subdomain values.
          nxl_mg_save = nxl_mg(l)
          nxr_mg_save = nxr_mg(l)
          nys_mg_save = nys_mg(l)
          nyn_mg_save = nyn_mg(l)
          nzt_mg_save = nzt_mg(l)
          nxl_mg(l) = mg_loc_ind(1,myid)
          nxr_mg(l) = mg_loc_ind(2,myid)
          nys_mg(l) = mg_loc_ind(3,myid)
          nyn_mg(l) = mg_loc_ind(4,myid)
          nzt_mg(l) = mg_loc_ind(5,myid)
          nxl_l = nxl_mg(l)
          nxr_l = nxr_mg(l)
          nys_l = nys_mg(l)
          nyn_l = nyn_mg(l)
          nzt_l = nzt_mg(l)
!
!--       Allocate temporary array for topography heights on coarser grid level. Initialize all
!--       elements as fluid. This is the setting for the masking method, where the multigrid solver
!--       works like the FFT-solver, it "runs through the topography grid points" and is aware of
!--       the topography only via the divergence, which has been set zero for these points.
          ALLOCATE( topo_tmp(nzb:nzt_l+1,nys_l-1:nyn_l+1,nxl_l-1:nxr_l+1) )
          topo_tmp(:,:,:) = IBSET( topo_tmp(:,:,:), 0 )

!
!--       Set the topography grid points.
          IF ( .NOT. masking_method )  THEN

             DO  i = nxl_l, nxr_l
                DO  j = nys_l, nyn_l
                   DO  k = nzb, nzt_l
                      k_topo = k * inc
                      j_topo = j * inc
                      i_topo = i * inc
                      topo_tmp(k,j,i) = topo_flags(k_topo,j_topo,i_topo)
                   ENDDO
                ENDDO
             ENDDO
             topo_tmp(nzt_l+1,:,:) = topo_tmp(nzt_l,:,:)
          ENDIF
!
!--       In case that Neumann-BCs at the bottom (nzb) and/or top (nzt+1) are set, consider them
!--       via setting a wall at these points. This way, the BCs are used via the implicit
!--       Neumann-BCs that the multigrid-solvers assumes at walls.
          IF ( ibc_p_b == 1 )  topo_tmp(nzb,:,:)     = IBCLR( topo_tmp(nzb,:,:), 0 )
          IF ( ibc_p_t == 1 )  topo_tmp(nzt_l+1,:,:) = IBCLR( topo_tmp(nzt_l+1,:,:), 0 )
!
!--       Exchange ghost points on respective multigrid level.
#if defined( __parallel )
          CALL exchange_horiz_int( topo_tmp, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1, type_xz_int(l), &
                                   type_yz_int(l) )
#else
          CALL exchange_horiz_int( topo_tmp, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1 )
#endif
!
!--       Set walls at the total domain boundaries in case of non-cyclic boundary conditions
!--       to consider Neumann-BCs for pressure at these boundaries.
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
                topo_tmp(:,-1,:) = IBCLR( topo_tmp(:,-1,:), 0 )
             ENDIF
             IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
                topo_tmp(:,nyn_l+1,:) = IBCLR( topo_tmp(:,nyn_l+1,:), 0 )
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
                topo_tmp(:,:,-1) = IBCLR( topo_tmp(:,:,-1), 0 )
             ENDIF
             IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
                topo_tmp(:,:,nxr_l+1) = IBCLR( topo_tmp(:,:,nxr_l+1), 0 )
             ENDIF
          ENDIF
!
!--       Now set the flags, based on the topography/wall settings above.
          DO  i = nxl_l, nxr_l
             DO  j = nys_l, nyn_l
                DO  k = nzb, nzt_l
!
!--                Inside/outside building (inside building does not need further tests for walls).
                   IF ( .NOT. BTEST( topo_tmp(k,j,i), 0 ) )  THEN

                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 6 )

                   ELSE
!
!--                   Bottom wall.
                      IF ( .NOT. BTEST( topo_tmp(k-1,j,i), 0 ) )  THEN
                         gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 0 )
                      ENDIF
!
!--                   Top wall.
                      IF ( .NOT. BTEST( topo_tmp(k+1,j,i), 0 ) )  THEN
                         gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 1 )
                      ENDIF
!
!--                   South wall.
                      IF ( .NOT. BTEST( topo_tmp(k,j-1,i), 0 ) )  THEN
                         gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 2 )
                      ENDIF
!
!--                   North wall.
                      IF ( .NOT. BTEST( topo_tmp(k,j+1,i), 0 ) )  THEN
                         gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 3 )
                      ENDIF
!
!--                   Left wall.
                      IF ( .NOT. BTEST( topo_tmp(k,j,i-1), 0 ) )  THEN
                         gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 4 )
                      ENDIF
!
!--                   Right wall.
                      IF ( .NOT. BTEST( topo_tmp(k,j,i+1), 0 ) )  THEN
                         gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 5 )
                      ENDIF

                   ENDIF

                ENDDO
             ENDDO
          ENDDO

          gl(0)%flags(nzt_l+1,:,:) = gl(0)%flags(nzt_l,:,:)

#if defined( __parallel )
          CALL exchange_horiz_int( gl(0)%flags, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1,              &
                                   type_xz_int(l), type_yz_int(l) )
#else
          CALL exchange_horiz_int( gl(0)%flags, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, 1 )
#endif

!
!--       Set non-cyclic boundary conditions.
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
                gl(0)%flags(:,-1,:) = gl(0)%flags(:,0,:)
                gl(0)%flags(:,-1,:) = IBSET( gl(0)%flags(:,-1,:), 6 )
             ENDIF
             IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
                gl(0)%flags(:,nyn_l+1,:) = gl(0)%flags(:,nyn_l,:)
                gl(0)%flags(:,nyn_l+1,:) = IBSET( gl(0)%flags(:,nyn_l+1,:), 6 )
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
                gl(0)%flags(:,:,-1) = gl(0)%flags(:,:,0)
                gl(0)%flags(:,:,-1) = IBSET( gl(0)%flags(:,:,-1), 6 )
             ENDIF
             IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
                gl(0)%flags(:,:,nxr_l+1) = gl(0)%flags(:,:,nxr_l)
                gl(0)%flags(:,:,nxr_l+1) = IBSET( gl(0)%flags(:,:,nxr_l+1), 6 )
             ENDIF
          ENDIF

!
!--       Set special flags to be used in routine prolong. They are not required at ghost layers,
!--       so exchange_horiz does not have to be applied.
!--       Flags indicate for fluid points in the respective finer level, that the area is
!--       completely covered by buildings/topography in the coarser level.
          DO  i = nxl_l, nxr_l
             DO  j = nys_l, nyn_l
                DO  k = nzb, nzt_l
!
!--                Walls to the left and right.
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )  .AND.  BTEST( gl(0)%flags(k,j,i+1), 6 ) )  THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 7 )
                   ENDIF
!
!--                Walls to the south and north.
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )  .AND.  BTEST( gl(0)%flags(k,j+1,i), 6 ) )  THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 8 )
                   ENDIF
!
!--                Walls below and above.
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )  .AND.  BTEST( gl(0)%flags(k+1,j,i), 6 ) )  THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 9 )
                   ENDIF
!
!--                Walls around the center of the planes stretched by four points of the
!--                coarse grid cube (xy-plane).
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )  .AND.  BTEST( gl(0)%flags(k,j,i+1), 6 )  .AND.         &
                        BTEST( gl(0)%flags(k,j+1,i), 6 ) )                                               &
                   THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 10 )
                   ENDIF
!
!--                Walls around the center of the planes stretched by four points of the
!--                coarse grid cube (xz-plane).
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )  .AND.  BTEST( gl(0)%flags(k,j,i+1), 6 )  .AND.         &
                        BTEST( gl(0)%flags(k+1,j,i), 6 ) )                                               &
                   THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 11 )
                   ENDIF
!
!--                Walls around the center of the planes stretched by four points of the
!--                coarse grid cube (yz-plane).
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )  .AND.  BTEST( gl(0)%flags(k,j+1,i), 6 )  .AND.         &
                        BTEST( gl(0)%flags(k+1,j,i), 6 ) )                                               &
                   THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 12 )
                   ENDIF
!
!--                Walls around the middle of coarse grid cube.
                   IF ( BTEST( gl(0)%flags(k,j,i), 6 )    .AND.  BTEST( gl(0)%flags(k,j,i+1), 6 )  .AND.       &
                        BTEST( gl(0)%flags(k,j+1,i), 6 )  .AND.  BTEST( gl(0)%flags(k+1,j,i), 6 ) )            &
                   THEN
                      gl(0)%flags(k,j,i) = IBSET( gl(0)%flags(k,j,i), 13 )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          DEALLOCATE( topo_tmp )
!
!--       Sort flags to even/odd (not required for the noopt-version).
          IF ( even_odd_decomposition )  CALL sort_k_to_even_odd_blocks( gl(0)%flags, l )

!
!--       Restore the correct indices of this level.
          nxl_mg(l) = nxl_mg_save
          nxr_mg(l) = nxr_mg_save
          nys_mg(l) = nys_mg_save
          nyn_mg(l) = nyn_mg_save
          nzt_mg(l) = nzt_mg_save

       ENDIF  ! mg_switch_to_pe0_level

!
!--    Set grid point increment for the next level (use only every 2nd point of the current level).
       inc = inc * 2

    ENDDO

!
!-- Copy the flag arrays to the device via manual deep copy. The TYPE array must be copied first,
!-- and the the flag arrays for each array element.
    IF ( mg_switch_to_pe0_level > 0 )  THEN
       !$ACC ENTER DATA COPYIN(gl(0:maximum_grid_level)) IF(enable_openacc)
    ELSE
       !$ACC ENTER DATA COPYIN(gl(1:maximum_grid_level)) IF(enable_openacc)
    ENDIF
    DO  i = 0, maximum_grid_level
       IF ( i == 0 )  THEN
          IF ( mg_switch_to_pe0_level > 0 )  THEN
             !$ACC ENTER DATA COPYIN(gl(0)%flags(nzb:mg_loc_ind(5,myid)+1,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1)) IF(enable_openacc)
          ENDIF
       ELSE
          !$ACC ENTER DATA COPYIN(gl(i)%flags(nzb:nzt_mg(i)+1,nys_mg(i)-1:nyn_mg(i)+1,nxl_mg(i)-1:nxr_mg(i)+1)) IF(enable_openacc)
       ENDIF
    ENDDO

!
!-- Calculate start indices for the red/black decomposition that are used in the non-unroll loops,
!-- and those required for the ghost point exchange of red/black grid points.
    ALLOCATE( ileft(4,2,maximum_grid_level),                                                       &
              jsouth(4,2,maximum_grid_level),                                                      &
              kbottom(4,2,maximum_grid_level),                                                     &
              ktop(4,2,maximum_grid_level),                                                        &
              unroll(maximum_grid_level) )

    ALLOCATE( ileft_for_nyn_recv(2,2,maximum_grid_level),                                          &
              ileft_for_nyn_send(2,2,maximum_grid_level),                                          &
              ileft_for_nys_recv(2,2,maximum_grid_level),                                          &
              ileft_for_nys_send(2,2,maximum_grid_level),                                          &
              jsouth_for_nxl_recv(2,2,maximum_grid_level),                                         &
              jsouth_for_nxl_send(2,2,maximum_grid_level),                                         &
              jsouth_for_nxr_recv(2,2,maximum_grid_level),                                         &
              jsouth_for_nxr_send(2,2,maximum_grid_level),                                         &
              kbottom_for_nxl_recv(2,2,maximum_grid_level),                                        &
              kbottom_for_nxl_send(2,2,maximum_grid_level),                                        &
              kbottom_for_nxr_recv(2,2,maximum_grid_level),                                        &
              kbottom_for_nxr_send(2,2,maximum_grid_level),                                        &
              kbottom_for_nyn_recv(2,2,maximum_grid_level),                                        &
              kbottom_for_nyn_send(2,2,maximum_grid_level),                                        &
              kbottom_for_nys_recv(2,2,maximum_grid_level),                                        &
              kbottom_for_nys_send(2,2,maximum_grid_level),                                        &
              ktop_for_nxl_recv(2,2,maximum_grid_level),                                           &
              ktop_for_nxl_send(2,2,maximum_grid_level),                                           &
              ktop_for_nxr_recv(2,2,maximum_grid_level),                                           &
              ktop_for_nxr_send(2,2,maximum_grid_level),                                           &
              ktop_for_nyn_recv(2,2,maximum_grid_level),                                           &
              ktop_for_nyn_send(2,2,maximum_grid_level),                                           &
              ktop_for_nys_recv(2,2,maximum_grid_level),                                           &
              ktop_for_nys_send(2,2,maximum_grid_level) )

    DO  l = 1, maximum_grid_level

       unroll(l) = ( MOD( nyn_mg(l)-nys_mg(l)+1, 4 ) == 0  .AND.                                   &
                     MOD( nxr_mg(l)-nxl_mg(l)+1, 2 ) == 0 )

!
!--    Set loop start indices for the red/black decomposition. Four separate loops (1st array index)
!--    are required per color (red or black, 2nd array index).
!--    The optimized solver requires different k indices than the non-optimized version, because
!--    of the even/odd decomposition.
       ileft(1,1:2,l)   = nxl_mg(l)
       jsouth(1,1,l)    = nys_mg(l) + 1
       jsouth(1,2,l)    = nys_mg(l)
       IF ( even_odd_decomposition )  THEN
          kbottom(1,1:2,l) = even_odd_level(l) + 1
          ktop(1,1:2,l)    = nzt_mg(l)
       ELSE
          kbottom(1,1:2,l) = nzb + 1
          ktop(1,1:2,l)    = nzt_mg(l)
       ENDIF

       ileft(2,1:2,l)   = nxl_mg(l) + 1
       jsouth(2,1,l)    = nys_mg(l)
       jsouth(2,2,l)    = nys_mg(l) + 1
       IF ( even_odd_decomposition )  THEN
          kbottom(2,1:2,l) = even_odd_level(l) + 1
          ktop(2,1:2,l)    = nzt_mg(l)
       ELSE
          kbottom(2,1:2,l) = nzb + 1
          ktop(2,1:2,l)    = nzt_mg(l)
       ENDIF

       ileft(3,1:2,l)   = nxl_mg(l)
       jsouth(3,1,l)    = nys_mg(l)
       jsouth(3,2,l)    = nys_mg(l) + 1
       IF ( even_odd_decomposition )  THEN
          kbottom(3,1:2,l) = nzb + 1
          ktop(3,1:2,l)    = even_odd_level(l)
       ELSE
          kbottom(3,1:2,l) = nzb + 2
          ktop(3,1:2,l)    = nzt_mg(l)
       ENDIF

       ileft(4,1:2,l)   = nxl_mg(l) + 1
       jsouth(4,1,l)    = nys_mg(l) + 1
       jsouth(4,2,l)    = nys_mg(l)
       IF ( even_odd_decomposition )  THEN
          kbottom(4,1:2,l) = nzb + 1
          ktop(4,1:2,l)    = even_odd_level(l)
       ELSE
          kbottom(4,1:2,l) = nzb + 2
          ktop(4,1:2,l)    = nzt_mg(l)
       ENDIF

!
!--    The red/black decomposition requires that on the lower i,j indices need to start
!--    alternatively with an even or odd value on the coarsest grid level, depending on the core-id,
!--    and if the subdomain has an uneven number of gridpoints along x/y. The respective index
!--    adjustments are done now.
!--    TODO: This should not be restricted to l=1 in case that non-uniform subdomains will be
!--          allowed.
       IF ( l == 1  .AND.  MOD( myidx, 2 ) /= 0  .AND.  MOD( nxr_mg(l) - nxl_mg(l), 2 ) == 0 )  THEN

          ileft(1,1:2,l) = nxl_mg(l) + 1
          ileft(2,1:2,l) = nxl_mg(l)
          ileft(3,1:2,l) = nxl_mg(l) + 1
          ileft(4,1:2,l) = nxl_mg(l)

       ENDIF

       IF ( l == 1  .AND.  MOD( myidy, 2 ) /= 0  .AND.  MOD( nyn_mg(l) - nys_mg(l), 2 ) == 0 )  THEN

          jsouth(1,1,l) = nys_mg(l)
          jsouth(2,1,l) = nys_mg(l) + 1
          jsouth(3,1,l) = nys_mg(l) + 1
          jsouth(4,1,l) = nys_mg(l)

          jsouth(1,2,l) = nys_mg(l) + 1
          jsouth(2,2,l) = nys_mg(l)
          jsouth(3,2,l) = nys_mg(l)
          jsouth(4,2,l) = nys_mg(l) + 1

       ENDIF

!
!--    Determine j and k start indices for ghost point exchange of lateral boundaries.
!--    The block below is already prepared to work with non-uniform subdomains.
       DO  i = 1, 4
          DO  color = 1, 2
!
!--          Left/right boundaries.
!--          Treat only those cases where the left index starts at the left boundary, because
!--          otherwise the respective loop i does not calculate for boundary points.
!--          There are exactly two cases (two values of i), for which the below if (ileft) is true.
!--          One of the cases contains uneven k values (originally, in the optimized version before
!--          the even/odd decomposition), which starts with even_odd_level(l) + 1 (optimized), or
!--          nzb+1 (non-optimized), the other one starts with with nzb+1 (optimized) or
!--          nzb+2 (non-optmized) and contains the even values.
             IF ( even_odd_decomposition )  THEN
                kbottom_uneven = even_odd_level(l) + 1
             ELSE
                kbottom_uneven = nzb + 1
             ENDIF
             IF ( ileft(i,color,l) == nxl_mg(l) )  THEN

                IF ( kbottom(i,color,l) == kbottom_uneven )  THEN
!
!--                Uneven k (upper half in the optimized version).
                   jsouth_for_nxl_send(1,color,l)  = jsouth(i,color,l)
                   kbottom_for_nxl_send(1,color,l) = kbottom(i,color,l)
                   ktop_for_nxl_send(1,color,l)    = ktop(i,color,l)
                   IF ( MOD( nxr_mg(l) - nxl_mg(l), 2 ) == 0 )  THEN
                      jsouth_for_nxr_send(1,color,l) = jsouth(i,color,l)
                   ELSE
                      jsouth_for_nxr_send(1,color,l) = jsouth(i,color,l) + 1
                   ENDIF
                   kbottom_for_nxr_send(1,color,l) = kbottom(i,color,l)
                   ktop_for_nxr_send(1,color,l)    = ktop(i,color,l)

                ELSE
!
!--                Even k (lower half in the optimized version).
                   jsouth_for_nxl_send(2,color,l)  = jsouth(i,color,l)
                   kbottom_for_nxl_send(2,color,l) = kbottom(i,color,l)
                   ktop_for_nxl_send(2,color,l)    = ktop(i,color,l)
                   IF ( MOD( nxr_mg(l) - nxl_mg(l), 2 ) == 0 )  THEN
                      jsouth_for_nxr_send(2,color,l) = jsouth(i,color,l)
                   ELSE
                      jsouth_for_nxr_send(2,color,l) = jsouth(i,color,l) + 1
                   ENDIF
                   kbottom_for_nxr_send(2,color,l) = kbottom(i,color,l)
                   ktop_for_nxr_send(2,color,l)    = ktop(i,color,l)

                ENDIF

             ENDIF

!
!--          South/north boundaries.
!--          Treat only those cases where the south index starts at the south boundary, because
!--          otherwise the respective loop j does not calculate for boundary points.
             IF ( jsouth(i,color,l) == nys_mg(l) )  THEN

                IF ( kbottom(i,color,l) == kbottom_uneven )  THEN

                   ileft_for_nys_send(1,color,l)   = ileft(i,color,l)
                   kbottom_for_nys_send(1,color,l) = kbottom(i,color,l)
                   ktop_for_nys_send(1,color,l)    = ktop(i,color,l)
                   IF ( MOD( nyn_mg(l) - nys_mg(l), 2 ) == 0 )  THEN
                      ileft_for_nyn_send(1,color,l) = ileft(i,color,l)
                   ELSE
                      ileft_for_nyn_send(1,color,l) = ileft(i,color,l) + 1
                   ENDIF
                   kbottom_for_nyn_send(1,color,l) = kbottom(i,color,l)
                   ktop_for_nyn_send(1,color,l)    = ktop(i,color,l)

                ELSE

                   ileft_for_nys_send(2,color,l) = ileft(i,color,l)
                   kbottom_for_nys_send(2,color,l) = kbottom(i,color,l)
                   ktop_for_nys_send(2,color,l)    = ktop(i,color,l)
                   IF ( MOD( nyn_mg(l) - nys_mg(l), 2 ) == 0 )  THEN
                      ileft_for_nyn_send(2,color,l) = ileft(i,color,l)
                   ELSE
                      ileft_for_nyn_send(2,color,l) = ileft(i,color,l) + 1
                   ENDIF
                   kbottom_for_nyn_send(2,color,l) = kbottom(i,color,l)
                   ktop_for_nyn_send(2,color,l)    = ktop(i,color,l)

                ENDIF

             ENDIF

          ENDDO
       ENDDO

!
!--    Adjust i and j start index for send, so that ghost points along j direction are included in
!--    the exchange.
!--    Two j sweeps required because the above calculations may generate a jsouth index of 3.
!--    The i index also requires adjustment because above calculations sometimes generate an ileft
!--    index of 2.
       DO  j = 1, 2
          DO  i = 1, 2
             DO  color = 1, 2
                IF ( ileft_for_nys_send(i,color,l) > nxl_mg(l)+1 )  THEN
                   ileft_for_nys_send(i,color,l) = ileft_for_nys_send(i,color,l) - 2
                ENDIF
                IF ( ileft_for_nyn_send(i,color,l) > nxl_mg(l)+1 )  THEN
                   ileft_for_nyn_send(i,color,l) = ileft_for_nyn_send(i,color,l) - 2
                ENDIF
                IF ( jsouth_for_nxl_send(i,color,l) > nys_mg(l) )  THEN
                   jsouth_for_nxl_send(i,color,l) = jsouth_for_nxl_send(i,color,l) - 2
                ENDIF
                IF ( jsouth_for_nxr_send(i,color,l) > nys_mg(l) )  THEN
                   jsouth_for_nxr_send(i,color,l) = jsouth_for_nxr_send(i,color,l) - 2
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    ENDDO
!
!-- Exchange the send indices with the respective neighbours, where they are the receive indices.
#if defined( __parallel )
    bufsize = 2 * 2 * maximum_grid_level
!
!-- Send left boundary, receive right one (asynchronous), only red or black points
    req(1:4)  = 0
    req_count = 0
    CALL MPI_ISEND( jsouth_for_nxl_send(1,1,1), bufsize, MPI_INTEGER, pleft, req_count, comm2d,    &
                    req(req_count+1), ierr )
    CALL MPI_IRECV( jsouth_for_nxr_recv(1,1,1), bufsize, MPI_INTEGER, pright, req_count, comm2d,   &
                    req(req_count+2), ierr )
!
!-- Send right boundary, receive left one (asynchronous)
    CALL MPI_ISEND( jsouth_for_nxr_send(1,1,1), bufsize, MPI_INTEGER, pright, req_count+1, comm2d, &
                    req(req_count+3), ierr )
    CALL MPI_IRECV( jsouth_for_nxl_recv(1,1,1), bufsize, MPI_INTEGER, pleft, req_count+1, comm2d,  &
                    req(req_count+4), ierr )

    CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!-- Send south boundary, receive north one (asynchronous)
    req(1:4)  = 0
    req_count = 0
    CALL MPI_ISEND( ileft_for_nys_send(1,1,1), bufsize, MPI_INTEGER, psouth, req_count, comm2d,    &
                    req(req_count+1), ierr )
    CALL MPI_IRECV( ileft_for_nyn_recv(1,1,1), bufsize, MPI_INTEGER, pnorth, req_count, comm2d,    &
                    req(req_count+2), ierr )
!
!-- Send north boundary, receive south one (asynchronous)
    CALL MPI_ISEND( ileft_for_nyn_send(1,1,1), bufsize, MPI_INTEGER, pnorth, req_count+1, comm2d,  &
                    req(req_count+3), ierr )
    CALL MPI_IRECV( ileft_for_nys_recv(1,1,1), bufsize, MPI_INTEGER, psouth, req_count+1, comm2d,  &
                    req(req_count+4), ierr )

    CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!-- Send left boundary, receive right one (asynchronous), only red or black points
    req(1:4)  = 0
    req_count = 0
    CALL MPI_ISEND( kbottom_for_nxl_send(1,1,1), bufsize, MPI_INTEGER, pleft, req_count, comm2d,   &
                    req(req_count+1), ierr )
    CALL MPI_IRECV( kbottom_for_nxr_recv(1,1,1), bufsize, MPI_INTEGER, pright, req_count, comm2d,  &
                    req(req_count+2), ierr )
!
!-- Send right boundary, receive left one (asynchronous)
    CALL MPI_ISEND( kbottom_for_nxr_send(1,1,1), bufsize, MPI_INTEGER, pright, req_count+1,        &
                    comm2d, req(req_count+3), ierr )
    CALL MPI_IRECV( kbottom_for_nxl_recv(1,1,1), bufsize, MPI_INTEGER, pleft, req_count+1,         &
                    comm2d, req(req_count+4), ierr )

    CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!-- Send south boundary, receive north one (asynchronous)
    req(1:4)  = 0
    req_count = 0
    CALL MPI_ISEND( kbottom_for_nys_send(1,1,1), bufsize, MPI_INTEGER, psouth, req_count, comm2d,  &
                    req(req_count+1), ierr )
    CALL MPI_IRECV( kbottom_for_nyn_recv(1,1,1), bufsize, MPI_INTEGER, pnorth, req_count, comm2d,  &
                    req(req_count+2), ierr )
!
!-- Send north boundary, receive south one (asynchronous)
    CALL MPI_ISEND( kbottom_for_nyn_send(1,1,1), bufsize, MPI_INTEGER, pnorth, req_count+1,        &
                    comm2d, req(req_count+3), ierr )
    CALL MPI_IRECV( kbottom_for_nys_recv(1,1,1), bufsize, MPI_INTEGER, psouth, req_count+1,        &
                    comm2d, req(req_count+4), ierr )

    CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!-- Send left boundary, receive right one (asynchronous), only red or black points
    req(1:4)  = 0
    req_count = 0
    CALL MPI_ISEND( ktop_for_nxl_send(1,1,1), bufsize, MPI_INTEGER, pleft, req_count, comm2d,      &
                    req(req_count+1), ierr )
    CALL MPI_IRECV( ktop_for_nxr_recv(1,1,1), bufsize, MPI_INTEGER, pright, req_count, comm2d,     &
                    req(req_count+2), ierr )
!
!-- Send right boundary, receive left one (asynchronous)
    CALL MPI_ISEND( ktop_for_nxr_send(1,1,1), bufsize, MPI_INTEGER, pright, req_count+1, comm2d,   &
                    req(req_count+3), ierr )
    CALL MPI_IRECV( ktop_for_nxl_recv(1,1,1), bufsize, MPI_INTEGER, pleft, req_count+1, comm2d,    &
                    req(req_count+4), ierr )

    CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!-- Send south boundary, receive north one (asynchronous)
    req(1:4)  = 0
    req_count = 0
    CALL MPI_ISEND( ktop_for_nys_send(1,1,1), bufsize, MPI_INTEGER, psouth, req_count, comm2d,     &
                    req(req_count+1), ierr )
    CALL MPI_IRECV( ktop_for_nyn_recv(1,1,1), bufsize, MPI_INTEGER, pnorth, req_count, comm2d,     &
                    req(req_count+2), ierr )
!
!-- Send north boundary, receive south one (asynchronous)
    CALL MPI_ISEND( ktop_for_nyn_send(1,1,1), bufsize, MPI_INTEGER, pnorth, req_count+1, comm2d,   &
                    req(req_count+3), ierr )
    CALL MPI_IRECV( ktop_for_nys_recv(1,1,1), bufsize, MPI_INTEGER, psouth, req_count+1, comm2d,   &
                    req(req_count+4), ierr )

    CALL MPI_WAITALL( 4, req, wait_stat, ierr )
#endif

!
!-- Adjust i and j start index for receive, so that ghost points along j direction are included in
!-- the exchange.
!-- Two j sweeps required because the above calculations may generate a jsouth index of 3.
!-- The i index also requires adjustment because above calculations sometimes generate an ileft
!-- index of 2.
    DO  l = 1, maximum_grid_level
       DO  j = 1, 2
          DO  i = 1, 2
             DO  color = 1, 2

                IF ( ileft_for_nys_recv(i,color,l) > nxl_mg(l)+1 )  THEN
                   ileft_for_nys_recv(i,color,l) = ileft_for_nys_recv(i,color,l) - 2
                ENDIF
                IF ( ileft_for_nyn_recv(i,color,l) > nxl_mg(l)+1 )  THEN
                   ileft_for_nyn_recv(i,color,l) = ileft_for_nyn_recv(i,color,l) - 2
                ENDIF
                IF ( jsouth_for_nxl_recv(i,color,l) > nys_mg(l) )  THEN
                   jsouth_for_nxl_recv(i,color,l) = jsouth_for_nxl_recv(i,color,l) - 2
                ENDIF
                IF ( jsouth_for_nxr_recv(i,color,l) > nys_mg(l) )  THEN
                   jsouth_for_nxr_recv(i,color,l) = jsouth_for_nxr_recv(i,color,l) - 2
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- Unset the grid level.
    grid_level = 0

 END SUBROUTINE poismg_init

 END MODULE poismg_mod
