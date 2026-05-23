!> @file poismg_noopt_mod.f90
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
! Copyright 1997-2025 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid V- or W-Cycle scheme.
!>
!> This multigrid method was originally developed for PALM by Joerg Uhlenbrock,
!> September 2000 - July 2001.
!--------------------------------------------------------------------------------------------------!
 MODULE poismg_noopt_mod

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
               enable_openacc,                                                                     &
               nesting_offline

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

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

    USE poismg_mod

    PRIVATE

    SAVE

    INTERFACE poismg_noopt
       MODULE PROCEDURE poismg_noopt
    END INTERFACE poismg_noopt

    PUBLIC poismg_noopt

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid V- or W-Cycle scheme.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poismg_noopt

    USE arrays_3d,                                                                                 &
        ONLY:  d,                                                                                  &
               p_loc

    USE control_parameters,                                                                        &
        ONLY:  current_timestep_number,                                                            &
               ibc_p_t,                                                                            &
               message_string

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
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  r   !<


    CALL cpu_log( log_point_s(29), 'poismg_noopt', 'start' )
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
       ALLOCATE( r(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1)  )
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

    DO WHILE ( residual_norm > residual_limit  .OR.  mgcycles < maximum_mgcycles )

       CALL next_mg_level_noopt( d, p_loc, p3, r)

!
!--    Calculate the residual if the user has not preset the number of cycles to be performed.
       IF ( maximum_mgcycles == 0 )  THEN

          CALL resid_noopt( d, p_loc, r )

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

    !$ACC END DATA
    DEALLOCATE( p3, r )

 END SUBROUTINE poismg_noopt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the residual of the perturbation pressure.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE resid_noopt( f_mg, p_mg, r )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< index variable along x
    INTEGER(iwp) ::  j    !< index variable along y
    INTEGER(iwp) ::  k    !< index variable along z
    INTEGER(iwp) ::  l    !< index indicating grid level

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


    CALL cpu_log( log_point_s(53), 'resid_noopt', 'start' )

    l = grid_level

    !$OMP PARALLEL PRIVATE (i,j,k)
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl_mg(l), nxr_mg(l)
       DO  j = nys_mg(l), nyn_mg(l)
!          !$ACC LOOP VECTOR
          DO  k = nzb+1, nzt_mg(l)
             pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
             pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
             pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
             pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
             pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
             pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
             r(k,j,i) = f_mg(k,j,i) - rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    - rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    - f2_mg(k,l) * pkpji - f3_mg(k,l) * pkmji                      &
                                    + f1_mg(k,l) * p_mg(k,j,i)
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

    CALL cpu_log( log_point_s(53), 'resid_noopt', 'stop' )

 END SUBROUTINE resid_noopt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the residual on the next coarser grid with "full weighting" scheme.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE restrict_noopt( f_mg, r )

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
    INTEGER(iwp) ::  l    !< index indicating finer grid level

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

    CALL cpu_log( log_point_s(54), 'restrict_noopt', 'start' )

    !$OMP PARALLEL PRIVATE (i,j,k,ic,jc,kc, rkjim,rkjip,rkjpi,rkjmi,rkjmim,rkjpim, &
    !$OMP rkjmip, rkjpip,rkmji,rkmjim,rkmjip,rkmjpi,rkmjmi,rkmjmim,rkmjpim,rkmjmip,&
    !$OMP rkmjpip          )
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  ic = nxl_mg(l), nxr_mg(l)
       DO  jc = nys_mg(l), nyn_mg(l)
          i = 2 * ic
          j = 2 * jc
          !DIR$ IVDEP
!          !$ACC LOOP VECTOR
          DO  kc = nzb+1, nzt_mg(l)
             k = 2 * kc - 1
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
             rkmji   = MERGE( r(k,j,i), r(k-1,j,i), BTEST( gl(l+1)%flags(k-1,j,i), 6 ) )
             rkmjim  = MERGE( r(k,j,i), r(k-1,j,i-1), BTEST( gl(l+1)%flags(k-1,j,i-1), 6 ) )
             rkmjip  = MERGE( r(k,j,i), r(k-1,j,i+1), BTEST( gl(l+1)%flags(k-1,j,i+1), 6 ) )
             rkmjpi  = MERGE( r(k,j,i), r(k-1,j+1,i), BTEST( gl(l+1)%flags(k-1,j+1,i), 6 ) )
             rkmjmi  = MERGE( r(k,j,i), r(k-1,j-1,i), BTEST( gl(l+1)%flags(k-1,j-1,i), 6 ) )
             rkmjmim = MERGE( r(k,j,i), r(k-1,j-1,i-1), BTEST( gl(l+1)%flags(k-1,j-1,i-1), 6 ) )
             rkmjpim = MERGE( r(k,j,i), r(k-1,j+1,i-1), BTEST( gl(l+1)%flags(k-1,j+1,i-1), 6 ) )
             rkmjmip = MERGE( r(k,j,i), r(k-1,j-1,i+1), BTEST( gl(l+1)%flags(k-1,j-1,i+1), 6 ) )
             rkmjpip = MERGE( r(k,j,i), r(k-1,j+1,i+1), BTEST( gl(l+1)%flags(k-1,j+1,i+1), 6 ) )
             rkpji   = MERGE( r(k,j,i), r(k+1,j,i), BTEST( gl(l+1)%flags(k+1,j,i), 6 ) )
             rkpjim  = MERGE( r(k,j,i), r(k+1,j,i-1), BTEST( gl(l+1)%flags(k+1,j,i-1), 6 ) )
             rkpjip  = MERGE( r(k,j,i), r(k+1,j,i+1), BTEST( gl(l+1)%flags(k+1,j,i+1), 6 ) )
             rkpjpi  = MERGE( r(k,j,i), r(k+1,j+1,i), BTEST( gl(l+1)%flags(k+1,j+1,i), 6 ) )
             rkpjmi  = MERGE( r(k,j,i), r(k+1,j-1,i), BTEST( gl(l+1)%flags(k+1,j-1,i), 6 ) )
             rkpjmim = MERGE( r(k,j,i), r(k+1,j-1,i-1), BTEST( gl(l+1)%flags(k+1,j-1,i-1), 6 ) )
             rkpjpim = MERGE( r(k,j,i), r(k+1,j+1,i-1), BTEST( gl(l+1)%flags(k+1,j+1,i-1), 6 ) )
             rkpjmip = MERGE( r(k,j,i), r(k+1,j-1,i+1), BTEST( gl(l+1)%flags(k+1,j-1,i+1), 6 ) )
             rkpjpip = MERGE( r(k,j,i), r(k+1,j+1,i+1), BTEST( gl(l+1)%flags(k+1,j+1,i+1), 6 ) )

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
    !$OMP END PARALLEL

!
!-- Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are implicitly
!-- treated via the flags array.
    CALL exchange_horiz( f_mg, 1, grid_level = grid_level, mg_switch_to_pe0 = mg_switch_to_pe0 )

!
!-- Dirichlet boundary conditions at bottom and top of the domain. Neumann BCs are implicitly
!-- considered in the calculations above. Points may be within buildings, but that doesn't matter.
    IF ( ibc_p_b == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f_mg(nzb,:,: ) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    IF ( ibc_p_t == 0 )  THEN
       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       f_mg(nzt_mg(l)+1,:,: ) = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    CALL cpu_log( log_point_s(54), 'restrict_noopt', 'stop' )

 END SUBROUTINE restrict_noopt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the correction of the perturbation pressure to the next finer grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE prolong_noopt( p, temp )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !< index variable along x on coarser grid level
    INTEGER(iwp) ::  j   !< index variable along y on coarser grid level
    INTEGER(iwp) ::  k   !< index variable along z on coarser grid level
    INTEGER(iwp) ::  l   !< index indicating finer grid level
    INTEGER(iwp) ::  lm1 !< index for flags indicating coarser grid level (considering the switch to PE0 level)

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


    CALL cpu_log( log_point_s(55), 'prolong_noopt', 'start' )

    l = grid_level
!
!-- Choose index for the lower level flag array.
    lm1 = grid_level - 1
!
!-- A special index 0 is required when switching from the total domain (switch_to_pe0_level)
!-- to the next finer level, because here the prolongation already calculates on the subdomains
!-- and not on the total domain. The regular flag-array for this level is defined for the total
!-- domain.
    IF ( ( l-1 ) == mg_switch_to_pe0_level )  lm1 = 0

    !$OMP PARALLEL PRIVATE (i,j,k)
    !$OMP DO
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC DEFAULT(PRESENT) IF(enable_openacc)
    DO  i = nxl_mg(l-1), nxr_mg(l-1)
       DO  j = nys_mg(l-1), nyn_mg(l-1)

          !DIR$ IVDEP
          DO  k = nzb+1, nzt_mg(l-1)
!
!--          Store pressure at surrounding grid points and apply Neumann boundary conditions in
!--          case of a wall.
             pkjip   = MERGE( p(k,j,i), p(k,j,i+1),     BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkjpi   = MERGE( p(k,j,i), p(k,j+1,i),     BTEST( gl(lm1)%flags(k,j,i), 3 ) )
             pkpji   = MERGE( p(k,j,i), p(k+1,j,i),     BTEST( gl(lm1)%flags(k,j,i), 1 ) )
             pkjpip  = MERGE( p(k,j,i), p(k,j+1,i+1),   BTEST( gl(lm1)%flags(k,j,i), 3 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkpjip  = MERGE( p(k,j,i), p(k+1,j,i+1),   BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
             pkpjpi  = MERGE( p(k,j,i), p(k+1,j+1,i),   BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 3 ) )
             pkpjpip = MERGE( p(k,j,i), p(k+1,j+1,i+1), BTEST( gl(lm1)%flags(k,j,i), 1 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 3 )  .OR.     &
                                                        BTEST( gl(lm1)%flags(k,j,i), 5 ) )
!
!--          Points of the coarse grid are directly stored on the next finer grid.
             temp(2*k-1,2*j,2*i) = p(k,j,i)
             temp(2*k-1,2*j,2*i) = MERGE( 0.0_wp, temp(2*k-1,2*j,2*i),                             &
                                          BTEST( gl(lm1)%flags(k,j,i), 6 ) )
!
!--          Points between two coarse-grid points.
             temp(2*k-1,2*j,2*i+1) = 0.5_wp * ( p(k,j,i) + pkjip )
             temp(2*k-1,2*j,2*i+1) = MERGE( 0.0_wp, temp(2*k-1,2*j,2*i+1),                         &
                                            BTEST( gl(lm1)%flags(k,j,i), 7 ))
             temp(2*k-1,2*j+1,2*i) = 0.5_wp * ( p(k,j,i) + pkjpi )
             temp(2*k-1,2*j+1,2*i) = MERGE( 0.0_wp, temp(2*k-1,2*j+1,2*i),                         &
                                            BTEST( gl(lm1)%flags(k,j,i), 8 ))
             temp(2*k,2*j,2*i)     = 0.5_wp * ( p(k,j,i) + pkpji )
             temp(2*k,2*j,2*i)     = MERGE( 0.0_wp, temp(2*k,2*j,2*i),                             &
                                            BTEST( gl(lm1)%flags(k,j,i), 9 ) )
!
!--          Points in the center of the planes stretched by four points of the coarse grid cube.
             temp(2*k-1,2*j+1,2*i+1) = 0.25_wp * ( p(k,j,i) + pkjip + pkjpi + pkjpip )
             temp(2*k-1,2*j+1,2*i+1) = MERGE( 0.0_wp, temp(2*k-1,2*j+1,2*i+1),                     &
                                              BTEST( gl(lm1)%flags(k,j,i), 10 ) )
             temp(2*k,2*j,2*i+1)   = 0.25_wp * ( p(k,j,i) + pkjip + pkpji + pkpjip )
             temp(2*k,2*j,2*i+1)   = MERGE( 0.0_wp, temp(2*k,2*j,2*i+1),                           &
                                            BTEST( gl(lm1)%flags(k,j,i), 11 ) )
             temp(2*k,2*j+1,2*i)   = 0.25_wp * ( p(k,j,i) + pkjpi + pkpji + pkpjpi )
             temp(2*k,2*j+1,2*i)   = MERGE( 0.0_wp, temp(2*k,2*j+1,2*i),                           &
                                            BTEST( gl(lm1)%flags(k,j,i), 12 ) )
!
!--          Points in the middle of coarse grid cube.
             temp(2*k,2*j+1,2*i+1) = 0.125_wp * ( p(k,j,i) + pkjip + pkjpi  + pkjpip +             &
                                                            pkpji + pkpjip + pkpjpi + pkpjpip )
             temp(2*k,2*j+1,2*i+1) = MERGE( 0.0_wp, temp(2*k,2*j+1,2*i+1),                         &
                                            BTEST( gl(lm1)%flags(k,j,i), 13 ) )
          ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL

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

    CALL cpu_log( log_point_s(55), 'prolong_noopt', 'stop' )

 END SUBROUTINE prolong_noopt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Relaxation method for the multigrid scheme. A Gauss-Seidel iteration with 3D-Red-Black
!> decomposition (GS-RB) is used.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE redblack_noopt( f_mg, p_mg )

    USE control_parameters,                                                                        &
        ONLY:  ibc_p_b,                                                                            &
               ibc_p_t

    IMPLICIT NONE

    INTEGER(iwp) ::  color  !< grid point color, either red (1) or black (2)
    INTEGER(iwp) ::  i      !< index variable along x
    INTEGER(iwp) ::  ic     !< index variable along x
    INTEGER(iwp) ::  j      !< index variable along y
    INTEGER(iwp) ::  jc     !< index variable along y
    INTEGER(iwp) ::  jj     !< index variable along y
    INTEGER(iwp) ::  k      !< index variable along z
    INTEGER(iwp) ::  l      !< grid level
    INTEGER(iwp) ::  m      !< loop number
    INTEGER(iwp) ::  n      !< loop variable Gauß-Seidel iterations

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


    l = grid_level

    DO  n = 1, ngsrb

!       p_mg(:,:,:) = 0.0

       DO  color = 1, 2

          IF ( .NOT. unroll(l) )  THEN

             CALL cpu_log( log_point_s(36), 'redblack_no_unroll_noopt', 'start' )

!
!--          Four loops are required to cover all points of one color.
             DO  m = 1,4
!
!--             Without unrolling of loops, no cache optimization.
                !$ACC PARALLEL LOOP INDEPENDENT COLLAPSE(3) GANG &
                !$ACC DEFAULT(PRESENT) IF(enable_openacc)
                DO  i = ileft(m,1,l), nxr_mg(l), 2
                   DO  j = jsouth(m,color,l), nyn_mg(l), 2
                      DO  k = kbottom(m,1,l), nzt_mg(l), 2
                         pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                         pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                         pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                         pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                         pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                         pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                         p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                       &
                                       ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )          &
                                       + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )          &
                                       + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                  &
                                       - f_mg(k,j,i)                                               &
                                       )
!                      p_mg(k,j,i) = color
                      ENDDO
                   ENDDO
                ENDDO
                !$ACC END PARALLEL LOOP

             ENDDO

             CALL cpu_log( log_point_s(36), 'redblack_no_unroll_noopt', 'stop' )

          ELSE

!
!--          Loop unrolling along y, only one i loop for better cache use.
             CALL cpu_log( log_point_s(38), 'redblack_unroll_noopt', 'start' )

             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
             !$ACC DEFAULT(PRESENT) IF(enable_openacc)
             DO  ic = nxl_mg(l), nxr_mg(l), 2
                DO  jc = nys_mg(l), nyn_mg(l), 4
                   i  = ic
                   jj = jc+2-color
!                   !$ACC LOOP VECTOR
                   DO  k = nzb+1, nzt_mg(l), 2
                      j = jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                   i  = ic+1
                   jj = jc+color-1
!                   !$ACC LOOP VECTOR
                   DO  k = nzb+1, nzt_mg(l), 2
                      j =jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                   &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                   i  = ic
                   jj = jc+color-1
!                   !$ACC LOOP VECTOR
                   DO  k = nzb+2, nzt_mg(l), 2
                      j =jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                   i  = ic+1
                   jj = jc+2-color
!                   !$ACC LOOP VECTOR
                   DO  k = nzb+2, nzt_mg(l), 2
                      j =jj
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color

                      j = jj+2
                      pkjip = MERGE( p_mg(k,j,i), p_mg(k,j,i+1), BTEST( gl(l)%flags(k,j,i), 5 ) )
                      pkjim = MERGE( p_mg(k,j,i), p_mg(k,j,i-1), BTEST( gl(l)%flags(k,j,i), 4 ) )
                      pkjpi = MERGE( p_mg(k,j,i), p_mg(k,j+1,i), BTEST( gl(l)%flags(k,j,i), 3 ) )
                      pkjmi = MERGE( p_mg(k,j,i), p_mg(k,j-1,i), BTEST( gl(l)%flags(k,j,i), 2 ) )
                      pkpji = MERGE( p_mg(k,j,i), p_mg(k+1,j,i), BTEST( gl(l)%flags(k,j,i), 1 ) )
                      pkmji = MERGE( p_mg(k,j,i), p_mg(k-1,j,i), BTEST( gl(l)%flags(k,j,i), 0 ) )
                      p_mg(k,j,i) = 1.0_wp / f1_mg(k,l) *                                          &
                                    ( rho_air_mg(k,l) * ddx2_mg(l) * ( pkjip + pkjim )             &
                                    + rho_air_mg(k,l) * ddy2_mg(l) * ( pkjpi + pkjmi )             &
                                    + f2_mg(k,l) * pkpji  + f3_mg(k,l) * pkmji                     &
                                    - f_mg(k,j,i)                                                  &
                                    )
!                      p_mg(k,j,i) = color
                   ENDDO

                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP

             CALL cpu_log( log_point_s(38), 'redblack_unroll_noopt', 'stop' )

          ENDIF

!
!--       Ghost point exchange. Neumann conditions for non-cyclic horizontal boundaries are
!--       implicitly treated via the flags array. In case of sufficiently large data,
!--       contiguous buffers are used in exchange_horiz_rb to only exchange data of the respective
!--       color. The threshold of 900 is empirical and may require adjustment to optimize
!--       performance.
!--       Levels where total domain is on PE0 do not require optimized exchange.
#if defined( __parallel )
          IF ( ( ngp_xz(l) >= 900  .OR.  ngp_yz(l) >= 900 )  .AND.  .NOT. mg_switch_to_pe0  .AND.  &
               npex /= 1  .AND.  npey /= 1  )                                                      &
          THEN
             CALL exchange_horiz_rb( p_mg, 1, color = color, kinc = 2,                             &
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

  !        WRITE(9,*) '*** level = ', l, '  color = ', color, '  unroll = ', unroll(l)
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
  !        DO  k = nzb, nzt_mg(l)+1
  !           WRITE(9,*) '*** k = ', k
!             IF ( adjust_lower_j_index )  THEN
!                DO  j = nyn_mg(grid_level)+1, save_nys_mg-1, -1
!                   IF ( adjust_lower_i_index )  THEN
!                      WRITE(9,'(A,I3,1X,70I1)')  'j=', j, ( NINT( p_mg(k,j,i) ), i = save_nxl_mg-1, nxr_mg(grid_level)+1 )
!                   ELSE
!                      WRITE(9,'(A,I3,1X,70I1)')  'j=', j, ( NINT( p_mg(k,j,i) ), i = nxl_mg(grid_level)-1, nxr_mg(grid_level)+1 )
!                   ENDIF
!                ENDDO
!             ELSE
  !              DO  j = nyn_mg(grid_level)+1, nys_mg(grid_level)-1, -1
!                   IF ( adjust_lower_i_index )  THEN
!                      WRITE(9,'(A,I3,1X,70I1)')  'j=', j, ( NINT( p_mg(k,j,i) ), i = save_nxl_mg-1, nxr_mg(grid_level)+1 )
!                   ELSE
  !                    WRITE(9,'(A,I3,1X,70I1)')  'j=', j, ( NINT( p_mg(k,j,i) ), i = nxl_mg(grid_level)-1, nxr_mg(grid_level)+1 )
!                   ENDIF
  !              ENDDO
!             ENDIF
  !        ENDDO

       ENDDO

    ENDDO

 END SUBROUTINE redblack_noopt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gather subdomain data from all PEs. All PEs will get the same data, so that they can all
!> do the calculations for the levels redundantly (they would have to wait for PE0 anyhow).
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE mg_gather_noopt( f2, f2_sub )

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


    CALL cpu_log( log_point_s(34), 'mg_gather_noopt', 'start' )

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
       f2(nzb,:,:) = 0.0_wp
    ENDIF

    IF ( ibc_p_t == 1 )  THEN
       f2(nzt_mg(grid_level)+1,:,:) = 0.0_wp
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
       f2(nzb,:,:) = f2(nzb+1,:,:)
    ENDIF

    IF ( ibc_p_t == 1 )  THEN
       f2(nzt_mg(grid_level)+1,:,:) = f2(nzt_mg(grid_level),:,:)
    ENDIF

    CALL cpu_log( log_point_s(34), 'mg_gather_noopt', 'stop' )

 END SUBROUTINE mg_gather_noopt
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Scatter the subdomain data. Since all PEs did the calculations, MPI is not required to scatter
!> them.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE mg_scatter_noopt( p2, p2_sub )

    IMPLICIT NONE

    REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,  &
                        nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

    REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,        &
                        mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  p2_sub  !<


    CALL cpu_log( log_point_s(35), 'mg_scatter_noopt', 'start' )

    p2_sub = p2(:,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,                                       &
                  mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1)

    CALL cpu_log( log_point_s(35), 'mg_scatter_noopt', 'stop' )

 END SUBROUTINE mg_scatter_noopt
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This is where the multigrid technique takes place. V- and W- Cycle are implemented and steered by
!> the parameter "gamma_mg".
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE next_mg_level_noopt( f_mg, p_mg, p3, r )

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

       CALL redblack_noopt( f_mg, p_mg )

       ngsrb = ngsrb / 2

!       FLUSH( 9 )
!       CALL MPI_FINALIZE( ierr )
!       STOP 'next_mg_level'


    ELSEIF ( grid_level /= 1 )  THEN

       grid_level_count(grid_level) = grid_level_count(grid_level) + 1

!
!--    Solution on this grid level.
       CALL redblack_noopt( f_mg, p_mg )

!
!--    Determination of the residual on this level.
       CALL resid_noopt( f_mg, p_mg, r )

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

          CALL restrict_noopt( f2_sub, r )
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
          CALL mg_gather_noopt( f2, f2_sub )
#endif
!
!--       Set switch for routine exchange_horiz, that no ghostpoint exchange has to be carried out
!--       from now on, because PEs contain the total domain.
          mg_switch_to_pe0 = .TRUE.

          DEALLOCATE( f2_sub )

       ELSE

          CALL restrict_noopt( f2, r )

       ENDIF

       !$ACC KERNELS DEFAULT(PRESENT) IF(enable_openacc)
       p2 = 0.0_wp
       !$ACC END KERNELS

!
!--    Repeat the same procedure until the coarsest grid is reached.
       CALL next_mg_level_noopt( f2, p2, p3, r )

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

          CALL mg_scatter_noopt( p2, p2_sub )

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

          CALL prolong_noopt( p2_sub, p3 )

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

          CALL prolong_noopt( p2, p3 )

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
       CALL redblack_noopt( f_mg, p_mg )

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
!-- Reset counter for the next call of poismg_noopt.
    grid_level_count(grid_level) = 0

!
!-- Continue with the next finer level. nxl..nzt have to be set to the finer grid values, because
!-- these variables are needed for the exchange of ghost points in routine exchange_horiz.
    grid_level = grid_level + 1
    nxl = nxl_mg(grid_level)
    nxr = nxr_mg(grid_level)
    nys = nys_mg(grid_level)
    nyn = nyn_mg(grid_level)
    nzt = nzt_mg(grid_level)

 20 CONTINUE

    !$ACC END DATA

 END SUBROUTINE next_mg_level_noopt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates wall flags for each grid level of the multigrid-solver.
!--------------------------------------------------------------------------------------------------!

 END MODULE poismg_noopt_mod
