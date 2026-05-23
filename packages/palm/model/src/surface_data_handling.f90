!> @file surface_data_handling.f90
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
!> The module encompasses the allocation and initialization of surface arrays, and handles reading
!> and writing restart data. In addition, a further derived data structure is initialized, to
!> set boundary conditions at surfaces.
!> @todo Revise initialization of surface fluxes (especially for chemistry)
!> @todo Get rid-off deallocation routines in restarts
!--------------------------------------------------------------------------------------------------!
 MODULE surface_data_handling

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  heatflux_input_conversion,                                                          &
               momentumflux_input_conversion,                                                      &
               scalarflux_input_conversion,                                                        &
               waterflux_input_conversion,                                                         &
               x,                                                                                  &
               y,                                                                                  &
               zu,                                                                                 &
               zw

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar,                                                                               &
               spc_names

    USE chem_modules

    USE control_parameters

    USE cut_cell_topography_mod,                                                                   &
        ONLY:  cct

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
               nzt,                                                                                &
               topo_flags

    USE general_utilities,                                                                         &
        ONLY:  gridpoint_id

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE kinds

    USE model_1d_mod,                                                                              &
        ONLY:  ri1d,                                                                               &
               us1d,                                                                               &
               usws1d,                                                                             &
               vsws1d

    USE pegrid

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_surface_filetypes,                                                        &
               rrd_mpi_io,                                                                         &
               rrd_mpi_io_global_array,                                                            &
               rrd_mpi_io_surface,                                                                 &
               total_number_of_surface_elements,                                                   &
               wrd_mpi_io,                                                                         &
               wrd_mpi_io_global_array,                                                            &
               wrd_mpi_io_surface

    USE surface_mod

    IMPLICIT NONE

    INTEGER(iwp), DIMENSION(1) ::  ns_on_file      !< total number of surfaces on subdomain, required for writing restart data

    INTEGER(iwp), DIMENSION(0:6) ::  off_i = (/  0,  0,  0,  0, -1,  1,  0 /) !< offset indicies between reference and surface grid point in x, for
                                                                              !< upward, downward, northward, southward, eastward, westward, model top
    INTEGER(iwp), DIMENSION(0:6) ::  off_j = (/  0,  0, -1,  1,  0,  0,  0 /) !< offset indicies between reference and surface grid point in y
    INTEGER(iwp), DIMENSION(0:6) ::  off_k = (/ -1,  1,  0,  0,  0,  0,  1 /) !< offset indicies between reference and surface grid point in z


    LOGICAL, DIMENSION(:), ALLOCATABLE ::  mapped !< temporary array to check if input surface data read from file has been already mapped onto
                                                  !< the different surface types

    SAVE

    PRIVATE

    INTERFACE init_surface_grid_information
       MODULE PROCEDURE init_surface_grid_information
    END INTERFACE init_surface_grid_information

    INTERFACE init_single_surface_properties
       MODULE PROCEDURE init_single_surface_properties
    END INTERFACE init_single_surface_properties

    INTERFACE init_surfaces
       MODULE PROCEDURE init_surfaces
    END INTERFACE init_surfaces

    INTERFACE init_surface_arrays
       MODULE PROCEDURE init_surface_arrays
    END INTERFACE init_surface_arrays

    INTERFACE surface_rrd_local
       MODULE PROCEDURE surface_rrd_local_ftn
       MODULE PROCEDURE surface_rrd_local_mpi
    END INTERFACE surface_rrd_local

    INTERFACE surface_wrd_local
       MODULE PROCEDURE surface_wrd_local
    END INTERFACE surface_wrd_local

    INTERFACE surface_restore_elements
       MODULE PROCEDURE surface_restore_elements_1d
       MODULE PROCEDURE surface_restore_elements_2d
    END INTERFACE surface_restore_elements

#if defined( _OPENACC )
    INTERFACE enter_surface_arrays
       MODULE PROCEDURE enter_surface_arrays
    END INTERFACE

    INTERFACE exit_surface_arrays
       MODULE PROCEDURE exit_surface_arrays
    END INTERFACE
#endif

!
!-- Public subroutines and functions
    PUBLIC init_single_surface_properties,                                                         &
           init_surfaces,                                                                          &
           init_surface_arrays,                                                                    &
           init_surface_grid_information,                                                          &
           surface_restore_elements,                                                               &
           surface_rrd_local,                                                                      &
           surface_wrd_local

#if defined( _OPENACC )
    PUBLIC enter_surface_arrays,                                                                   &
           exit_surface_arrays
#endif

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data type for containing the boundary conditions at horizontal and vertical surfaces.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_surface_grid_information

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< loop index along x-direction
    INTEGER(iwp) ::  j  !< loop index along y-direction
    INTEGER(iwp) ::  k  !< loop index along y-direction
    INTEGER(iwp) ::  l  !< running index for differently aligned surfaces

    INTEGER(idp) ::  gp_id          !< grid-point ID
    INTEGER(iwp) ::  m              !< running index over all surface-adjacent grid points
    INTEGER(iwp) ::  num_bgp        !< number of boundary grid points on subdomain
    INTEGER(iwp) ::  num_wgp        !< number of surface-adjacent grid points on subdomain
    INTEGER(iwp) ::  num_kji        !< number of surfaces at (j,i)-grid point
    INTEGER(iwp) ::  start_index_hv !< local start index of surface elements at (j,i)-grid point

    INTEGER(idp), DIMENSION(:), ALLOCATABLE ::  bgp_id !< list of boundary-grid point IDs

!
!-- Initialize data structure for horizontal surfaces, i.e. count the number of surface elements,
!-- allocate and initialize the respective index arrays, and set the respective start and end
!-- indices at each (j,i)-location. Note, there are two sets of indices. One which covers all
!-- surface adjacent fluid grid points, and another which covers all boundary grid points.
!-- The first set is generally used to set boundary conditions. However, boundary grid points
!-- can be accessed by different surface adjacent grid points, especially at corner grid points.
!-- This comes with two implications: i) Neumann boundary conditions are overdetermined, and
!-- ii) specific actions applied to the boundary grid points, e.g. applying a cooling rate, can
!-- thus be carried out multiple times. While i) can't be avoided and is actually not relevant
!-- for the numerical solution, ii) would lead to larger discrepancies in the boundary conditions,
!-- which is avoided by using the second index set.
!-- Note, the first set defines the fluid index and an offset value. With this it is possible
!-- to access the boundary grid point (multiple accesses possible). The second set only defines
!-- the boundary grid point.
!-- First, count the number of wall-grid points, independent on surface orientation.
    num_wgp = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Check if current gridpoint belongs to the atmosphere.
             IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                DO  l = 0, 5
                   IF ( .NOT. BTEST( topo_flags(k+off_k(l),j+off_j(l),i+off_i(l)), 0 ) )           &
                      num_wgp = num_wgp + 1
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- Save the number of horizontal surface elements (surface adjacent grid points).
    bc_hv%ns = num_wgp
!
!-- Allocate index arrays and their offset values (indicating their orientation) for the
!-- surface adjacent grid points.
    ALLOCATE( bc_hv%i(1:bc_hv%ns) )
    ALLOCATE( bc_hv%j(1:bc_hv%ns) )
    ALLOCATE( bc_hv%k(1:bc_hv%ns) )
    ALLOCATE( bc_hv%ioff(1:bc_hv%ns) )
    ALLOCATE( bc_hv%joff(1:bc_hv%ns) )
    ALLOCATE( bc_hv%koff(1:bc_hv%ns) )
    ALLOCATE( bc_hv%start_index(nysg:nyng,nxlg:nxrg) )
    ALLOCATE( bc_hv%end_index(nysg:nyng,nxlg:nxrg) )
    bc_hv%start_index = 1
    bc_hv%end_index   = 0

!
!-- Now, create a list of boundary grid point IDs. Therefore, allocate an array with the size of
!-- the surface adjacent grid points. This is the upper limit of all boundary grid points.
!-- Here, this is used to calculate the number of boundary grid points.
    ALLOCATE( bgp_id(1:bc_hv%ns) )
    bgp_id = HUGE( -1_idp )

    num_wgp = 1
    num_bgp = 0
    start_index_hv = 1
    DO  i = nxl, nxr
       DO  j = nys, nyn

          num_kji = 0
          DO  k = nzb+1, nzt
!
!--          Check if current gridpoint belongs to the atmosphere.
             IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                DO  l = 0, 5
                   IF ( .NOT. BTEST( topo_flags(k+off_k(l),j+off_j(l),i+off_i(l)), 0 ) )  THEN
!
!--                   Store indices of surface adjacent grid points and their offset values.
                      bc_hv%i(num_wgp) = i
                      bc_hv%j(num_wgp) = j
                      bc_hv%k(num_wgp) = k
                      bc_hv%ioff(num_wgp) = off_i(l)
                      bc_hv%joff(num_wgp) = off_j(l)
                      bc_hv%koff(num_wgp) = off_k(l)

                      num_kji = num_kji + 1
                      num_wgp = num_wgp + 1
!
!--                   Count the number of boundary grid points using unique grid-point IDs.
                      gp_id = gridpoint_id( k + off_k(l), j + off_j(l), i + off_i(l) )
                      IF ( .NOT. ANY( gp_id == bgp_id(:) ) )  THEN
                         num_bgp = num_bgp + 1
                         bgp_id(num_bgp) = gp_id
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
          bc_hv%start_index(j,i) = start_index_hv
          bc_hv%end_index(j,i)   = bc_hv%start_index(j,i) + num_kji - 1
          start_index_hv         = bc_hv%end_index(j,i) + 1
       ENDDO
    ENDDO

!
!-- Store the number of boundary grid points.
    bc_hv%ns_bgp = num_bgp
!
!-- Allocate index arrays for the boundary grid points.
    ALLOCATE( bc_hv%i_bgp(1:bc_hv%ns_bgp) )
    ALLOCATE( bc_hv%j_bgp(1:bc_hv%ns_bgp) )
    ALLOCATE( bc_hv%k_bgp(1:bc_hv%ns_bgp) )
!
!-- Determine the indices of boundary grid points. Note, the loop runs over all surface adjacent
!-- grid points. First, resize the ID array.
    DEALLOCATE( bgp_id )
    ALLOCATE( bgp_id(1:bc_hv%ns_bgp) )
    bgp_id = HUGE( -1_idp )

    num_bgp = 1
    DO  m = 1, bc_hv%ns
       i = bc_hv%i(m) + bc_hv%ioff(m)
       j = bc_hv%j(m) + bc_hv%joff(m)
       k = bc_hv%k(m) + bc_hv%koff(m)

       gp_id = gridpoint_id( k, j, i )
       IF ( .NOT. ANY( gp_id == bgp_id(:) ) )  THEN
          bgp_id(num_bgp) = gp_id

          bc_hv%i_bgp(num_bgp) = i
          bc_hv%j_bgp(num_bgp) = j
          bc_hv%k_bgp(num_bgp) = k

          num_bgp = num_bgp + 1
       ENDIF
    ENDDO

    DEALLOCATE( bgp_id )
!
!-- Compute the total number of surface adjacent and boundary grid points. This information
!-- is sometimes useful, e.g., to compute domain averages.
    bc_hv%ns_tot     = 0
    bc_hv%ns_bgp_tot = 0
#if defined( __parallel )
    CALL MPI_ALLREDUCE( bc_hv%ns, bc_hv%ns_tot, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( bc_hv%ns_bgp, bc_hv%ns_bgp_tot, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    bc_hv%ns_tot     = bc_hv%ns
    bc_hv%ns_bgp_tot = bc_hv%ns_bgp
#endif

 END SUBROUTINE init_surface_grid_information

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal and vertical surfaces. Counts the number of default-, natural and urban
!> surfaces and allocates memory, respectively.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_surface_arrays

    IMPLICIT NONE

    INTEGER(iwp) ::  i            !< running index x-direction
    INTEGER(iwp) ::  ioff         !< running index x-direction with offset value
    INTEGER(iwp) ::  j            !< running index y-direction
    INTEGER(iwp) ::  joff         !< running index y-direction with offset value
    INTEGER(iwp) ::  k            !< running index z-direction
    INTEGER(iwp) ::  koff         !< running index z-direction with offset value
    INTEGER(iwp) ::  l            !< index variable for surface facing
    INTEGER(iwp) ::  num_def      !< number of default surfaces
    INTEGER(iwp) ::  num_def_up   !< number of upward-facing default surfaces on subdomain
    INTEGER(iwp) ::  num_def_vert !< number of vertical default surfaces on subdomain
    INTEGER(iwp) ::  num_lsm      !< number of natural surfaces
    INTEGER(iwp) ::  num_lsm_up   !< number of upward-facing LSM surfaces on subdomain
    INTEGER(iwp) ::  num_lsm_vert !< number of vertical LSM surfaces on subdomain
    INTEGER(iwp) ::  num_mom_u    !< number of surfaces on u-grid
    INTEGER(iwp) ::  num_mom_v    !< number of surfaces on v-grid
    INTEGER(iwp) ::  num_mom_w    !< number of surfaces on w-grid
    INTEGER(iwp) ::  num_top      !< number of top boundary grid points
    INTEGER(iwp) ::  num_usm      !< number of urban surfaces
    INTEGER(iwp) ::  num_usm_up   !< number of upward-facing USM surfaces on subdomain
    INTEGER(iwp) ::  num_usm_vert !< number of vertical USM surfaces on subdomain

    LOGICAL ::  building             !< flag indicating building grid point
    LOGICAL ::  terrain              !< flag indicating natural terrain grid point
    LOGICAL ::  unresolved_building  !< flag indicating a grid point where actually a building is
                                     !< defined but not resolved by the vertical grid

    num_def = 0
    num_lsm = 0
    num_top = 0
    num_usm = 0

    num_def_up = 0
    num_lsm_up = 0
    num_usm_up = 0

    num_def_vert = 0
    num_lsm_vert = 0
    num_usm_vert = 0

    num_mom_u = 0
    num_mom_v = 0
    num_mom_w = 0
!
!-- Surfaces are classified according to the input data read from static input file. If no input
!-- file is present, all surfaces are classified either as natural, urban, or default, depending on
!-- the setting of land_surface and urban_surface. To control this, use the control flag
!-- topo_no_distinct.
!-- Count number of surfaces (independent on facing) on local domain.
!-- There is no further treatment required in case of the cut-cell topography for the scalar grid,
!-- because the topography flags have been already set before.
!-- However, in contrast to the Cartesian topography, the surface momentum flux is treated
!-- individually in a separate data structure in case of the cut-cell topography.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Check if current scalar grid point belongs to the atmosphere.
             IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                DO  l = 0, 6
                   koff = k + off_k(l)
                   joff = j + off_j(l)
                   ioff = i + off_i(l)

                   IF ( koff == nzt+1  .AND.  l == 6  .AND.  use_top_fluxes )  THEN
                      num_top = num_top + 1
                   ELSEIF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 0 ) ) THEN
!
!--                   Determine flags indicating a terrain surface, a building surface,
                      terrain  = BTEST( topo_flags(koff,joff,ioff), 5 )  .OR.  topo_no_distinct
                      building = BTEST( topo_flags(koff,joff,ioff), 6 )  .OR.  topo_no_distinct
!
!--                   Unresolved_building indicates a surface with equal height as terrain but with a
!--                   non-grid resolved building on top. These surfaces will be flagged as urban
!--                   surfaces.
                      unresolved_building = BTEST( topo_flags(koff,joff,ioff), 5 )  .AND.          &
                                            BTEST( topo_flags(koff,joff,ioff), 6 )  .AND.          &
                                            urban_surface
!
!--                   Land-surface type
                      IF ( land_surface  .AND.  terrain  .AND.  .NOT. unresolved_building )  THEN
                         num_lsm = num_lsm + 1
!
!--                      Count the number of upward-facing LSM surfaces.
                         IF ( off_k(l) == -1 )  THEN
                            num_lsm_up = num_lsm_up + 1
!
!--                      Count the number of vertical LSM surfaces.
                         ELSEIF ( off_k(l) /= 1 )  THEN
                            num_lsm_vert = num_lsm_vert + 1
                         ENDIF
!
!--                   Urban surface tpye
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         num_usm = num_usm + 1
!
!--                      Count the number of upward-facing USM surfaces.
                         IF ( off_k(l) == -1 )  THEN
                            num_usm_up = num_usm_up + 1
!
!--                      Count the number of vertical USM surfaces.
                         ELSEIF ( off_k(l) /= 1 )  THEN
                            num_usm_vert = num_usm_vert + 1
                         ENDIF
!
!--                   Default-surface type
                      ELSEIF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
                         num_def = num_def + 1
!
!--                      Count the number of upward-facing DEF surfaces.
                         IF ( off_k(l) == -1 )  THEN
                            num_def_up = num_def_up + 1
!
!--                      Count the number of vertical DEF surfaces.
                         ELSEIF ( off_k(l) /= 1 )  THEN
                            num_def_vert = num_def_vert + 1
                         ENDIF
!
!--                   Unclassifified surface-grid point. Give error message.
                      ELSE
                         WRITE( message_string, '(A,A,3(I5,A),4(A,L1) )' )                         &
                                  'inconsistent surface element classification &at grid point ',   &
                                  '(i,j,k) = (', i, ',', j, ',', k, ') &', 'LSM = ', land_surface, &
                                  ', USM = ', urban_surface, ', terrain = ', terrain,              &
                                  ', building = ', building
                         CALL message( 'surface_mod', 'PAC0316', 2, 2, myid, 6, 0 )
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
!
!--          Treat surfaces on u- and v-grid individually in case of the cut-cell topography.
             IF ( cut_cell_topography )  THEN
!
!--             Count the number of surfaces relevant for the u-grid.
                IF ( BTEST( topo_flags(k,j,i), 1 ) )  THEN
                   DO  l = 0, 5
                      koff = k + off_k(l)
                      joff = j + off_j(l)
                      ioff = i + off_i(l)
                      IF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 1 ) ) num_mom_u = num_mom_u + 1
                   ENDDO
                ENDIF
!
!--             Count the number of surfaces relevant for the v-grid.
                IF ( BTEST( topo_flags(k,j,i), 2 ) )  THEN
                   DO  l = 0, 5
                      koff = k + off_k(l)
                      joff = j + off_j(l)
                      ioff = i + off_i(l)
                      IF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 2 ) ) num_mom_v = num_mom_v + 1
                   ENDDO
                ENDIF
!
!--             Count the number of surfaces relevant for the w-grid.
                IF ( BTEST( topo_flags(k,j,i), 3 ) )  THEN
                   DO  l = 0, 5
                      koff = k + off_k(l)
                      joff = j + off_j(l)
                      ioff = i + off_i(l)
                      IF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 3 ) ) num_mom_w = num_mom_w + 1
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO

!
!-- Store number of grid points in the data structure.
    surf_def%ns = num_def
    surf_lsm%ns = num_lsm
    surf_top%ns = num_top
    surf_usm%ns = num_usm

    surf_u%ns = num_mom_u
    surf_v%ns = num_mom_v
    surf_w%ns = num_mom_w
!
!-- Allocate required attributes for the different surface types.
    CALL allocate_surface_attributes( surf_def, nys, nyn, nxl, nxr )
    CALL allocate_surface_attributes( surf_lsm, nys, nyn, nxl, nxr )
    CALL allocate_surface_attributes( surf_usm, nys, nyn, nxl, nxr )
!
!-- WARNING!!!!  Following lines are a workaround to avoid compile errors (if array checks are
!--              switched on), that appear because arrays shf_agg and qsws_agg are arguments for
!--              the diffusion_s routine. They need to be allocated in any case, even if the
!--              LSM is not used!
!-- Workaround should be removed!!!!!
    IF ( .NOT. land_surface )  THEN
       ALLOCATE( surf_lsm%shf_agg(1:surf_lsm%ns)  )
       ALLOCATE( surf_lsm%qsws_agg(1:surf_lsm%ns) )
    ENDIF
!
!-- Allocate attributes required to compute the momentum fluxes for the cut-cell topography.
    IF ( cut_cell_topography )  THEN
       CALL allocate_surface_attributes_mom( surf_u, nys, nyn, nxl, nxr )
       CALL allocate_surface_attributes_mom( surf_v, nys, nyn, nxl, nxr )
       CALL allocate_surface_attributes_mom( surf_w, nys, nyn, nxl, nxr )
    ENDIF
!
!-- Allocate required attributes for model top
    CALL allocate_surface_attributes_top( surf_top, nys, nyn, nxl, nxr )
!
!-- Calculate the total number of upward-facing surfaces in the entire model domain of a type.
!-- This is required for statistical evaluation in flow-statistics.
#if defined( __parallel )
    CALL MPI_ALLREDUCE( num_def_up, surf_def%ns_tot_up, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( num_lsm_up, surf_lsm%ns_tot_up, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( num_usm_up, surf_usm%ns_tot_up, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    surf_def%ns_tot_up = num_def_up
    surf_lsm%ns_tot_up = num_lsm_up
    surf_usm%ns_tot_up = num_usm_up
#endif
!
!-- Vertical walls: At the moment this is only required for setting control flags. Note, at the
!-- current stage, the number of vertical surfaces has not been counted yet, nor the offset indices
!-- are set. Hence, determine the number of vertical surfaces by subtracting the number of upward-
!-- simplicity, use offset indices to identify vertical walls instead of northward, southward,
!-- eastward, and westward attribute.
#if defined( __parallel )
    CALL MPI_ALLREDUCE( num_def_vert, surf_def%ns_tot_v, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( num_lsm_vert, surf_lsm%ns_tot_v, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( num_usm_vert, surf_usm%ns_tot_v, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    surf_def%ns_tot_v = num_def_vert
    surf_lsm%ns_tot_v = num_lsm_vert
    surf_usm%ns_tot_v = num_usm_vert
#endif
!
!-- Set the flag for the existence of vertical urban/land surfaces.
    IF ( surf_lsm%ns_tot_v + surf_usm%ns_tot_v > 0 )  vertical_surfaces_exist = .TRUE.

 END SUBROUTINE init_surface_arrays

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter horizontal and vertical surfaces.
!--------------------------------------------------------------------------------------------------!
#if defined( _OPENACC )
 SUBROUTINE enter_surface_arrays

    !$ACC ENTER DATA &
    !$ACC COPYIN(surf_def) &
    !$ACC COPYIN(surf_lsm) &
    !$ACC COPYIN(surf_top) &
    !$ACC COPYIN(surf_usm) &
    !$ACC COPYIN(surf_u) &
    !$ACC COPYIN(surf_v) &
    !$ACC COPYIN(surf_w) IF(enable_openacc)
!
!-- Copy data in surf_def, surf_lsm, surf_usm and surf_top
    CALL enter_surface_attributes( surf_def )
    CALL enter_surface_attributes( surf_lsm )
    CALL enter_surface_attributes( surf_usm )
    CALL enter_surface_attributes_top( surf_top )
    CALL enter_surface_attributes_top( surf_u )
    CALL enter_surface_attributes_top( surf_v )
    CALL enter_surface_attributes_top( surf_w )
!
!-- WARNING!!!!  Following lines are a workaround because _agg arrays only exist for LSM
!--              surfaces and therefore cannot be treated above via enter_surface_attributes.
!-- Workaround should be removed!!!!!
    !$ACC ENTER DATA &
    !$ACC COPYIN(surf_lsm%shf_agg(1:surf_lsm%ns)) &
    !$ACC COPYIN(surf_lsm%qsws_agg(1:surf_lsm%ns)) IF(enable_openacc)

 END SUBROUTINE enter_surface_arrays
#endif

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit horizontal and vertical surfaces.
!--------------------------------------------------------------------------------------------------!
#if defined( _OPENACC )
 SUBROUTINE exit_surface_arrays
!
!-- Delete data in surf_def, surf_lsm, surf_usm and surf_top
    CALL exit_surface_attributes( surf_def )
    CALL exit_surface_attributes( surf_lsm )
    CALL exit_surface_attributes( surf_usm )
    CALL exit_surface_attributes_top( surf_top )

    !$ACC EXIT DATA &
    !$ACC DELETE(surf_def) &
    !$ACC DELETE(surf_lsm) &
    !$ACC DELETE(surf_top) &
    !$ACC DELETE(surf_usm) IF(enable_openacc)

 END SUBROUTINE exit_surface_arrays
#endif

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory surface types, except for top
!> fluxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE deallocate_surface_attributes( surfaces )

    IMPLICIT NONE


    TYPE(surf_type) ::  surfaces  !< respective surface type


    IF ( ALLOCATED( surfaces%start_index ) )  DEALLOCATE ( surfaces%start_index )
    IF ( ALLOCATED( surfaces%end_index   ) )  DEALLOCATE ( surfaces%end_index )
!
!-- Indices to locate surface element
    IF ( ALLOCATED( surfaces%i ) )  DEALLOCATE ( surfaces%i )
    IF ( ALLOCATED( surfaces%j ) )  DEALLOCATE ( surfaces%j )
    IF ( ALLOCATED( surfaces%k ) )  DEALLOCATE ( surfaces%k )

    IF ( ALLOCATED( surfaces%ioff ) )  DEALLOCATE ( surfaces%ioff )
    IF ( ALLOCATED( surfaces%joff ) )  DEALLOCATE ( surfaces%joff )
    IF ( ALLOCATED( surfaces%koff ) )  DEALLOCATE ( surfaces%koff )

    IF ( ALLOCATED( surfaces%iref ) )  DEALLOCATE ( surfaces%iref )
    IF ( ALLOCATED( surfaces%jref ) )  DEALLOCATE ( surfaces%jref )
    IF ( ALLOCATED( surfaces%kref ) )  DEALLOCATE ( surfaces%kref )
!
!-- Arrays indicating surface orientation.
    IF ( ALLOCATED( surfaces%downward  ) )  DEALLOCATE ( surfaces%downward )
    IF ( ALLOCATED( surfaces%eastward  ) )  DEALLOCATE ( surfaces%eastward )
    IF ( ALLOCATED( surfaces%northward ) )  DEALLOCATE ( surfaces%northward )
    IF ( ALLOCATED( surfaces%southward ) )  DEALLOCATE ( surfaces%southward )
    IF ( ALLOCATED( surfaces%upward    ) )  DEALLOCATE ( surfaces%upward )
    IF ( ALLOCATED( surfaces%westward  ) )  DEALLOCATE ( surfaces%westward )

    IF ( ALLOCATED( surfaces%bid ) )              DEALLOCATE ( surfaces%bid )
    IF ( ALLOCATED( surfaces%btype ) )            DEALLOCATE ( surfaces%btype )
    IF ( ALLOCATED( surfaces%ptype ) )            DEALLOCATE ( surfaces%ptype )
    IF ( ALLOCATED( surfaces%vtype ) )            DEALLOCATE ( surfaces%vtype )
    IF ( ALLOCATED( surfaces%wtype ) )            DEALLOCATE ( surfaces%wtype )
    IF ( ALLOCATED( surfaces%cut_cell_wall ) )    DEALLOCATE ( surfaces%cut_cell_wall )
    IF ( ALLOCATED( surfaces%cut_cell_roof ) )    DEALLOCATE ( surfaces%cut_cell_roof )
    IF ( ALLOCATED( surfaces%faces ) )            DEALLOCATE ( surfaces%faces )
    IF ( ALLOCATED( surfaces%num_edges ) )        DEALLOCATE ( surfaces%num_edges )
    IF ( ALLOCATED( surfaces%area ) )             DEALLOCATE ( surfaces%area )
    IF ( ALLOCATED( surfaces%wall_location_c ) )  DEALLOCATE ( surfaces%wall_location_c )

    IF ( ALLOCATED( surfaces%consider_stability ) )  DEALLOCATE ( surfaces%consider_stability )
    IF ( ALLOCATED( surfaces%tke_production     ) )  DEALLOCATE ( surfaces%tke_production )
!
!-- Surface-layer height (distance from surface).
    IF ( ALLOCATED( surfaces%z_mo ) )  DEALLOCATE ( surfaces%z_mo )
!
!-- Surface-parallel wind velocity
    IF ( ALLOCATED( surfaces%uvw_abs    ) )  DEALLOCATE( surfaces%uvw_abs )
    IF ( ALLOCATED( surfaces%uvw_abs_uv ) )  DEALLOCATE( surfaces%uvw_abs_uv )
    IF ( ALLOCATED( surfaces%uvw_abs_w  ) )  DEALLOCATE( surfaces%uvw_abs_w )
!
!-- Pre-calculated ln(z/z0)
    IF ( ALLOCATED( surfaces%ln_z_z0  ) )  DEALLOCATE ( surfaces%ln_z_z0 )
    IF ( ALLOCATED( surfaces%ln_z_z0h ) )  DEALLOCATE ( surfaces%ln_z_z0h )
    IF ( ALLOCATED( surfaces%ln_z_z0q ) )  DEALLOCATE ( surfaces%ln_z_z0q )
!
!-- Roughness
    IF ( ALLOCATED( surfaces%z0  ) )  DEALLOCATE ( surfaces%z0 )
    IF ( ALLOCATED( surfaces%z0h ) )  DEALLOCATE ( surfaces%z0h )
    IF ( ALLOCATED( surfaces%z0q ) )  DEALLOCATE ( surfaces%z0q )
!
!-- Friction velocity
    IF ( ALLOCATED( surfaces%us        ) )  DEALLOCATE ( surfaces%us )
    IF ( ALLOCATED( surfaces%us_uvgrid ) )  DEALLOCATE ( surfaces%us_uvgrid )
    IF ( ALLOCATED( surfaces%us_wgrid  ) )  DEALLOCATE ( surfaces%us_wgrid )
!
!-- Stability parameter
    IF ( ALLOCATED( surfaces%ol ) )  DEALLOCATE ( surfaces%ol )
!
!-- Bulk Richardson number
    IF ( ALLOCATED( surfaces%rib ) )  DEALLOCATE ( surfaces%rib )
!
!-- Vertical momentum fluxes of u and v
    IF ( ALLOCATED( surfaces%usws ) )  DEALLOCATE ( surfaces%usws )
    IF ( ALLOCATED( surfaces%vsws ) )  DEALLOCATE ( surfaces%vsws )
!
!-- Allocate arrays for surface momentum fluxes for u and v. These only apply at at vertical
!-- surfaces.
    IF ( ALLOCATED( surfaces%usvs ) )  DEALLOCATE ( surfaces%usvs )
    IF ( ALLOCATED( surfaces%vsus ) )  DEALLOCATE ( surfaces%vsus )
!
!-- Allocate array for surface momentum flux for w - wsus and wsvs
    IF ( ALLOCATED( surfaces%wsus_wsvs ) )  DEALLOCATE ( surfaces%wsus_wsvs )
!
!-- Allocate array for surface momentum flux for subgrid-scale tke wsus and wsvs; first index usvs
!-- or vsws, second index for wsus or wsvs, depending on surface.
    IF ( ALLOCATED( surfaces%mom_flux_tke ) )  DEALLOCATE ( surfaces%mom_flux_tke )
!
!-- Required in production_e
    IF ( .NOT. constant_diffusion )  THEN
       IF ( ALLOCATED( surfaces%u_0 ) )  DEALLOCATE ( surfaces%u_0 )
       IF ( ALLOCATED( surfaces%v_0 ) )  DEALLOCATE ( surfaces%v_0 )
    ENDIF
!
!-- Characteristic temperature and surface flux of sensible heat
    IF ( ALLOCATED( surfaces%ts  ) )  DEALLOCATE ( surfaces%ts  )
    IF ( ALLOCATED( surfaces%shf ) )  DEALLOCATE ( surfaces%shf )
!
!-- Surface temperature
    IF ( ALLOCATED( surfaces%pt_surface ) )  DEALLOCATE ( surfaces%pt_surface )
!
!-- Characteristic humidity and surface flux of latent heat
    IF ( humidity )  THEN
       IF ( ALLOCATED( surfaces%qs          ) )  DEALLOCATE ( surfaces%qs )
       IF ( ALLOCATED( surfaces%qsws        ) )  DEALLOCATE ( surfaces%qsws )
       IF ( ALLOCATED( surfaces%q_surface   ) )  DEALLOCATE ( surfaces%q_surface )
       IF ( ALLOCATED( surfaces%vpt_surface ) )  DEALLOCATE ( surfaces%vpt_surface )
    ENDIF
!
!-- Aerodynamic resistance.
    IF( ALLOCATED( surfaces%r_a ) )  DEALLOCATE( surfaces%r_a )
!
!-- Characteristic scalar and surface flux of scalar
    IF ( passive_scalar )  THEN
       IF ( ALLOCATED( surfaces%ss   ) )  DEALLOCATE ( surfaces%ss )
       IF ( ALLOCATED( surfaces%ssws ) )  DEALLOCATE ( surfaces%ssws )
    ENDIF
!
!-- Scaling parameter (cs*) and surface flux of chemical species
    IF ( air_chemistry )  THEN
       IF ( ALLOCATED( surfaces%css   ) )  DEALLOCATE ( surfaces%css )
       IF ( ALLOCATED( surfaces%cssws ) )  DEALLOCATE ( surfaces%cssws )
    ENDIF
!
!-- Arrays for storing potential temperature and mixing ratio at first grid level
    IF ( ALLOCATED( surfaces%pt1  ) )  DEALLOCATE ( surfaces%pt1 )
    IF ( ALLOCATED( surfaces%qv1  ) )  DEALLOCATE ( surfaces%qv1 )
    IF ( ALLOCATED( surfaces%vpt1 ) )  DEALLOCATE ( surfaces%vpt1 )

!
!--
    IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
       IF ( ALLOCATED( surfaces%qcs   ) )  DEALLOCATE ( surfaces%qcs )
       IF ( ALLOCATED( surfaces%ncs   ) )  DEALLOCATE ( surfaces%ncs )
       IF ( ALLOCATED( surfaces%qcsws ) )  DEALLOCATE ( surfaces%qcsws )
       IF ( ALLOCATED( surfaces%ncsws ) )  DEALLOCATE ( surfaces%ncsws )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
       IF ( ALLOCATED( surfaces%qrs   ) )  DEALLOCATE ( surfaces%qrs )
       IF ( ALLOCATED( surfaces%nrs   ) )  DEALLOCATE ( surfaces%nrs )
       IF ( ALLOCATED( surfaces%qrsws ) )  DEALLOCATE ( surfaces%qrsws )
       IF ( ALLOCATED( surfaces%nrsws ) )  DEALLOCATE ( surfaces%nrsws )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model .AND. surf_microphysics_ice_phase)  THEN
       IF ( ALLOCATED( surfaces%qis   ) )  DEALLOCATE ( surfaces%qis )
       IF ( ALLOCATED( surfaces%nis   ) )  DEALLOCATE ( surfaces%nis )
       IF ( ALLOCATED( surfaces%qisws ) )  DEALLOCATE ( surfaces%qisws )
       IF ( ALLOCATED( surfaces%nisws ) )  DEALLOCATE ( surfaces%nisws )
    ENDIF
!
!-- Salinity surface flux
    IF ( ALLOCATED( surfaces%sasws ) )  DEALLOCATE ( surfaces%sasws )
!
!-- Arrays for the normalized (to one) surface normal vector
    IF ( ALLOCATED( surfaces%n_eff   ) )  DEALLOCATE ( surfaces%n_eff )

    IF ( ALLOCATED( surfaces%n_s   ) )  DEALLOCATE ( surfaces%n_s )

 END SUBROUTINE deallocate_surface_attributes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for surface types, except for top fluxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE allocate_surface_attributes( surfaces, nys_l, nyn_l, nxl_l, nxr_l,                     &
                                         no_allocate_index_arrays )

    IMPLICIT NONE

    INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

    LOGICAL ::  allocate_index_arrays
    LOGICAL, INTENT(IN), OPTIONAL  :: no_allocate_index_arrays

    TYPE(surf_type) ::  surfaces  !< respective surface type


    IF ( PRESENT( no_allocate_index_arrays ) )  THEN
       allocate_index_arrays = .NOT. no_allocate_index_arrays
    ELSE
       allocate_index_arrays = .TRUE.
    ENDIF
!
!-- Allocate arrays for start and end index of horizontal surface type for each (j,i)-grid point.
!-- This is required e.g. in diffusion_x, which is called for each (j,i). In order to find the
!-- location where the respective flux is store within the surface-type, start- and end-index are
!-- stored for each (j,i). For example, each (j,i) can have several entries where fluxes for
!-- horizontal surfaces might be stored, e.g. for overhanging structures where several upward-facing
!-- surfaces might exist for given (j,i). If no surface of respective type exist at current (j,i),
!-- set indicies such that loop in diffusion routines will not be entered.
    IF ( allocate_index_arrays )  THEN
       ALLOCATE( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
    ENDIF
!
!-- Indices to locate surface element and corresponding atmosphere reference grid points.
    ALLOCATE ( surfaces%i(1:surfaces%ns) )
    ALLOCATE ( surfaces%j(1:surfaces%ns) )
    ALLOCATE ( surfaces%k(1:surfaces%ns) )

    ALLOCATE ( surfaces%ioff(1:surfaces%ns) )
    ALLOCATE ( surfaces%joff(1:surfaces%ns) )
    ALLOCATE ( surfaces%koff(1:surfaces%ns) )

    ALLOCATE ( surfaces%iref(1:surfaces%ns) )
    ALLOCATE ( surfaces%jref(1:surfaces%ns) )
    ALLOCATE ( surfaces%kref(1:surfaces%ns) )
!
!-- Arrays indicating surface orientation
    ALLOCATE ( surfaces%downward(1:surfaces%ns) )
    ALLOCATE ( surfaces%eastward (1:surfaces%ns) )
    ALLOCATE ( surfaces%northward(1:surfaces%ns) )
    ALLOCATE ( surfaces%southward(1:surfaces%ns) )
    ALLOCATE ( surfaces%upward(1:surfaces%ns) )
    ALLOCATE ( surfaces%westward(1:surfaces%ns) )

    ALLOCATE ( surfaces%bid(1:surfaces%ns) )
    ALLOCATE ( surfaces%btype(1:surfaces%ns) )
    ALLOCATE ( surfaces%ptype(1:surfaces%ns) )
    ALLOCATE ( surfaces%vtype(1:surfaces%ns) )
    ALLOCATE ( surfaces%wtype(1:surfaces%ns) )
    ALLOCATE ( surfaces%cut_cell_wall(1:surfaces%ns) )
    ALLOCATE ( surfaces%cut_cell_roof(1:surfaces%ns) )
    ALLOCATE ( surfaces%faces(cct%num_faces_vert,1:surfaces%ns) )
    ALLOCATE ( surfaces%num_edges(1:surfaces%ns) )
    ALLOCATE ( surfaces%area(1:surfaces%ns) )
    ALLOCATE ( surfaces%wall_location_c(cct%dim_3d,1:surfaces%ns) )

    ALLOCATE ( surfaces%consider_stability(1:surfaces%ns) )
    ALLOCATE ( surfaces%tke_production(1:surfaces%ns) )
!
!-- Surface-layer height (distance from surface).
    ALLOCATE ( surfaces%z_mo(1:surfaces%ns) )
!
!-- Surface-parallel wind velocity
    ALLOCATE ( surfaces%uvw_abs(1:surfaces%ns) )
    ALLOCATE ( surfaces%uvw_abs_uv(1:surfaces%ns) )
    ALLOCATE ( surfaces%uvw_abs_w(1:surfaces%ns) )
!
!-- Precalculated ln(z/z0)
    ALLOCATE( surfaces%ln_z_z0(1:surfaces%ns) )
    ALLOCATE( surfaces%ln_z_z0h(1:surfaces%ns) )
    ALLOCATE( surfaces%ln_z_z0q(1:surfaces%ns) )
!
!-- Roughness
    ALLOCATE ( surfaces%z0(1:surfaces%ns) )
    ALLOCATE ( surfaces%z0h(1:surfaces%ns) )
    ALLOCATE ( surfaces%z0q(1:surfaces%ns) )
!
!-- Friction velocity
    ALLOCATE ( surfaces%us(1:surfaces%ns) )
    ALLOCATE ( surfaces%us_uvgrid(1:surfaces%ns) )
    ALLOCATE ( surfaces%us_wgrid(1:surfaces%ns) )
!
!-- Stability parameter
    ALLOCATE ( surfaces%ol(1:surfaces%ns) )
!
!-- Bulk Richardson number
    ALLOCATE ( surfaces%rib(1:surfaces%ns) )
!
!-- Vertical momentum fluxes of u and v
    ALLOCATE ( surfaces%usws(1:surfaces%ns) )
    ALLOCATE ( surfaces%vsws(1:surfaces%ns) )
!
!-- Allocate arrays for surface momentum fluxes for u and v at vertical walls.
    ALLOCATE ( surfaces%usvs(1:surfaces%ns) )
    ALLOCATE ( surfaces%vsus(1:surfaces%ns) )
!
!-- Allocate array for surface momentum flux for w - wsus and wsvs
    ALLOCATE ( surfaces%wsus_wsvs(1:surfaces%ns) )
!
!-- Allocate array for surface momentum flux for subgrid-scale tke wsus and wsvs; first index usvs
!-- or vsws, second index for wsus or wsvs, depending on surface orientation.
    ALLOCATE ( surfaces%mom_flux_tke(0:1,1:surfaces%ns) )
!
!-- Required in production_e
    IF ( .NOT. constant_diffusion )  THEN
       ALLOCATE ( surfaces%u_0(1:surfaces%ns) )
       ALLOCATE ( surfaces%v_0(1:surfaces%ns) )
    ENDIF
!
!-- Characteristic temperature and surface flux of sensible heat
    ALLOCATE ( surfaces%ts(1:surfaces%ns) )
    ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!-- Surface temperature
    ALLOCATE ( surfaces%pt_surface(1:surfaces%ns) )
!
!-- Characteristic humidity, surface flux of latent heat, and surface virtual potential temperature
    IF ( humidity )  THEN
       ALLOCATE ( surfaces%qs(1:surfaces%ns) )
       ALLOCATE ( surfaces%qsws(1:surfaces%ns) )
       ALLOCATE ( surfaces%q_surface(1:surfaces%ns) )
       ALLOCATE ( surfaces%vpt_surface(1:surfaces%ns) )
    ENDIF
!
!-- Aerodynamic resistance.
    ALLOCATE( surfaces%r_a(1:surfaces%ns) )
!
!-- Characteristic scalar and surface flux of scalar
    IF ( passive_scalar )  THEN
       ALLOCATE ( surfaces%ss(1:surfaces%ns) )
       ALLOCATE ( surfaces%ssws(1:surfaces%ns) )
    ENDIF
!
!-- Scaling parameter (cs*) and surface flux of chemical species
    IF ( air_chemistry )  THEN
       ALLOCATE ( surfaces%css(1:nvar,1:surfaces%ns) )
       ALLOCATE ( surfaces%cssws(1:nvar,1:surfaces%ns) )
    ENDIF
!
!-- Arrays for storing potential temperature and mixing ratio at first grid level
    ALLOCATE ( surfaces%pt1(1:surfaces%ns) )
    ALLOCATE ( surfaces%qv1(1:surfaces%ns) )
    ALLOCATE ( surfaces%vpt1(1:surfaces%ns) )
!
!--
    IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
       ALLOCATE ( surfaces%qcs(1:surfaces%ns) )
       ALLOCATE ( surfaces%ncs(1:surfaces%ns) )
       ALLOCATE ( surfaces%qcsws(1:surfaces%ns) )
       ALLOCATE ( surfaces%ncsws(1:surfaces%ns) )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
       ALLOCATE ( surfaces%qrs(1:surfaces%ns) )
       ALLOCATE ( surfaces%nrs(1:surfaces%ns) )
       ALLOCATE ( surfaces%qrsws(1:surfaces%ns) )
       ALLOCATE ( surfaces%nrsws(1:surfaces%ns) )
    ENDIF

!
!--
    IF ( surf_bulk_cloud_model .AND. surf_microphysics_ice_phase)  THEN
       ALLOCATE ( surfaces%qis(1:surfaces%ns) )
       ALLOCATE ( surfaces%nis(1:surfaces%ns) )
       ALLOCATE ( surfaces%qisws(1:surfaces%ns) )
       ALLOCATE ( surfaces%nisws(1:surfaces%ns) )
    ENDIF

!
!-- Salinity surface flux
    IF ( ocean_mode )  ALLOCATE ( surfaces%sasws(1:surfaces%ns) )
!
!-- Allocate arrays for the normalized (to one) surface normal vector
    ALLOCATE ( surfaces%n_eff(1:surfaces%ns)   )
    ALLOCATE ( surfaces%n_s(1:surfaces%ns,1:3) )

 END SUBROUTINE allocate_surface_attributes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for surfaces used to compute momentum fluxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE allocate_surface_attributes_mom( surfaces, nys_l, nyn_l, nxl_l, nxr_l,                 &
                                             no_allocate_index_arrays )

    IMPLICIT NONE

    INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

    LOGICAL ::  allocate_index_arrays
    LOGICAL, INTENT(IN), OPTIONAL  :: no_allocate_index_arrays

    TYPE(surf_type) ::  surfaces  !< respective surface type


    IF ( PRESENT( no_allocate_index_arrays ) )  THEN
       allocate_index_arrays = .NOT. no_allocate_index_arrays
    ELSE
       allocate_index_arrays = .TRUE.
    ENDIF
!
!-- Allocate arrays for start and end index of horizontal surface type for each (j,i)-grid point.
!-- This is required to access surfaces gridpoint-wise, i.e. for each (j,i)-pair individually. For
!-- example, this is required for diffusion_u/v/w, which is called for each (j,i). Each (j,i)-pair
!-- can have several entries for diffenent surfaces located in each column. If no surfaces of
!-- respective type exist at current (j,i), set indicies such that loops will not be entered.
    IF ( allocate_index_arrays )  THEN
       ALLOCATE( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
    ENDIF
!
!-- Indices to locate surface element
    ALLOCATE ( surfaces%i(1:surfaces%ns) )
    ALLOCATE ( surfaces%j(1:surfaces%ns) )
    ALLOCATE ( surfaces%k(1:surfaces%ns) )

    ALLOCATE ( surfaces%ioff(1:surfaces%ns) )
    ALLOCATE ( surfaces%joff(1:surfaces%ns) )
    ALLOCATE ( surfaces%koff(1:surfaces%ns) )

    ALLOCATE ( surfaces%iref(1:surfaces%ns) )
    ALLOCATE ( surfaces%jref(1:surfaces%ns) )
    ALLOCATE ( surfaces%kref(1:surfaces%ns) )
!
!-- Arrays indicating surface orientation
    ALLOCATE ( surfaces%downward(1:surfaces%ns) )
    ALLOCATE ( surfaces%eastward (1:surfaces%ns) )
    ALLOCATE ( surfaces%northward(1:surfaces%ns) )
    ALLOCATE ( surfaces%southward(1:surfaces%ns) )
    ALLOCATE ( surfaces%upward(1:surfaces%ns) )
    ALLOCATE ( surfaces%westward(1:surfaces%ns) )

    ALLOCATE ( surfaces%consider_stability(1:surfaces%ns) )
!
!-- Surface-layer height (distance from surface).
    ALLOCATE ( surfaces%z_mo(1:surfaces%ns) )
!
!-- Surface-parallel wind velocity
    ALLOCATE ( surfaces%uvw_abs(1:surfaces%ns) )
    ALLOCATE ( surfaces%uvw_abs_uv(1:surfaces%ns) )
    ALLOCATE ( surfaces%uvw_abs_w(1:surfaces%ns) )
!
!-- Precalculated ln(z/z0)
    ALLOCATE( surfaces%ln_z_z0(1:surfaces%ns) )
!
!-- Roughness
    ALLOCATE ( surfaces%z0(1:surfaces%ns) )
!
!-- Friction velocity
    ALLOCATE ( surfaces%us(1:surfaces%ns) )
    ALLOCATE ( surfaces%us_uvgrid(1:surfaces%ns) )
    ALLOCATE ( surfaces%us_wgrid(1:surfaces%ns) )
!
!-- Stability parameter
    ALLOCATE ( surfaces%ol(1:surfaces%ns) )
!
!-- Vertical momentum fluxes of u and v
    ALLOCATE ( surfaces%usws(1:surfaces%ns) )
    ALLOCATE ( surfaces%vsws(1:surfaces%ns) )
!
!-- Allocate arrays for surface momentum fluxes for u and v at vertical walls.
    ALLOCATE ( surfaces%usvs(1:surfaces%ns) )
    ALLOCATE ( surfaces%vsus(1:surfaces%ns) )
!
!-- Allocate array for surface momentum flux for w - wsus and wsvs
    ALLOCATE ( surfaces%wsus_wsvs(1:surfaces%ns) )
!
!-- Allocate arrays for the normalized (to one) surface normal vector
    ALLOCATE ( surfaces%n_eff(1:surfaces%ns) )
    ALLOCATE ( surfaces%n_s(1:surfaces%ns,1:3) )

 END SUBROUTINE allocate_surface_attributes_mom


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit memory for upward and downward-facing horizontal surface types, except for top fluxes.
!--------------------------------------------------------------------------------------------------!
#if defined( _OPENACC )
 SUBROUTINE exit_surface_attributes( surfaces )

    IMPLICIT NONE

    TYPE(surf_type) ::  surfaces  !< respective surface type

    !$ACC EXIT DATA &
    !$ACC DELETE(surfaces%start_index(nys:nyn,nxl:nxr)) &
    !$ACC DELETE(surfaces%end_index(nys:nyn,nxl:nxr)) &
    !$ACC DELETE(surfaces%i(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%j(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%k(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%ioff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%joff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%koff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%iref(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%jref(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%kref(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%downward(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%eastward(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%northward(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%southward(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%upward(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%westward(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%consider_stability(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%tke_production(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%z_mo(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%uvw_abs(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%uvw_abs_uv(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%uvw_abs_w(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%ln_z_z0(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%ln_z_z0h(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%ln_z_z0q(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%n_eff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%n_s(1:surfaces%ns,1:3)) &
    !$ACC DELETE(surfaces%z0(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%z0h(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%z0q(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%us(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%us_uvgrid(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%us_wgrid(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%ol(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%rib(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%usws(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%vsws(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%usvs(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%vsus(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%wsus_wsvs(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%mom_flux_tke(0:1,1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%ts(1:surfaces%ns)) &
    !$ACC COPYOUT(surfaces%shf(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%pt_surface(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%pt1(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%qv1(1:surfaces%ns)) IF(enable_openacc)

    IF ( .NOT. constant_diffusion )  THEN
       !$ACC EXIT DATA &
       !$ACC DELETE(surfaces%u_0(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%v_0(1:surfaces%ns)) IF(enable_openacc)
    ENDIF

 END SUBROUTINE exit_surface_attributes
#endif

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter memory for upward and downward-facing horizontal surface types, except for top fluxes.
!--------------------------------------------------------------------------------------------------!
#if defined( _OPENACC )
 SUBROUTINE enter_surface_attributes( surfaces )

    IMPLICIT NONE

    TYPE(surf_type) ::  surfaces  !< respective surface type

    !$ACC ENTER DATA &
    !$ACC COPYIN(surfaces%start_index(nys:nyn,nxl:nxr)) &
    !$ACC COPYIN(surfaces%end_index(nys:nyn,nxl:nxr)) &
    !$ACC COPYIN(surfaces%i(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%j(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%k(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ioff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%joff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%koff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%iref(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%jref(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%kref(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%downward(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%eastward(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%northward(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%southward(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%upward(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%westward(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%consider_stability(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%tke_production(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%z_mo(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%uvw_abs(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%uvw_abs_uv(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%uvw_abs_w(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ln_z_z0(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ln_z_z0h(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ln_z_z0q(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%n_eff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%n_s(1:surfaces%ns,1:3)) &
    !$ACC COPYIN(surfaces%z0(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%z0h(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%z0q(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%us(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%qs(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%us_uvgrid(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%us_wgrid(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ol(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%rib(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%usws(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%vsws(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%qsws(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%usvs(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%vsus(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%wsus_wsvs(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%mom_flux_tke(0:1,1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ts(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%shf(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%pt1(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%qv1(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%vpt1(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%q_surface(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%vpt_surface(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%pt_surface(1:surfaces%ns)) IF(enable_openacc)

    IF ( .NOT. constant_diffusion )  THEN
       !$ACC ENTER DATA &
       !$ACC COPYIN(surfaces%u_0(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%v_0(1:surfaces%ns)) IF(enable_openacc)
    ENDIF

 END SUBROUTINE enter_surface_attributes
#endif

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for model-top fluxes
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE deallocate_surface_attributes_top( surfaces )

    IMPLICIT NONE


    TYPE(surf_type) ::  surfaces !< respective surface type

    DEALLOCATE ( surfaces%start_index )
    DEALLOCATE ( surfaces%end_index )
!
!-- Indices to locate surface (model-top) element
    DEALLOCATE ( surfaces%i )
    DEALLOCATE ( surfaces%j )
    DEALLOCATE ( surfaces%k )

    DEALLOCATE ( surfaces%ioff )
    DEALLOCATE ( surfaces%joff )
    DEALLOCATE ( surfaces%koff )
!
!-- Vertical momentum fluxes of u and v
    DEALLOCATE ( surfaces%usws )
    DEALLOCATE ( surfaces%vsws )
!
!-- Sensible heat flux
    DEALLOCATE ( surfaces%shf )
!
!-- Latent heat flux
    IF ( humidity  .OR.  ocean_run_coupled_to_atmosphere )  THEN
       DEALLOCATE ( surfaces%qsws )
    ENDIF
!
!-- Scalar flux
    IF ( passive_scalar )  THEN
       DEALLOCATE ( surfaces%ssws )
    ENDIF
!
!-- Chemical species flux
    IF ( air_chemistry )  THEN
       DEALLOCATE ( surfaces%cssws )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_morrison )  THEN
       DEALLOCATE ( surfaces%qcsws )
       DEALLOCATE ( surfaces%ncsws )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_seifert )  THEN
       DEALLOCATE ( surfaces%qrsws )
       DEALLOCATE ( surfaces%nrsws )
    ENDIF

!
!--
    IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_ice_phase )  THEN
       DEALLOCATE ( surfaces%qisws )
       DEALLOCATE ( surfaces%nisws )
    ENDIF
!
!-- Salinity flux
    IF ( ocean_mode )  DEALLOCATE ( surfaces%sasws )

 END SUBROUTINE deallocate_surface_attributes_top


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for model-top fluxes
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE allocate_surface_attributes_top( surfaces, nys_l, nyn_l, nxl_l, nxr_l )

    IMPLICIT NONE

    INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
    INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

    TYPE(surf_type) ::  surfaces !< respective surface type

    IF ( .NOT. ALLOCATED( surfaces%start_index ) )  THEN
       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       surfaces%start_index = 0
    ENDIF
    IF ( .NOT. ALLOCATED( surfaces%end_index ) )  THEN
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l) )
       surfaces%end_index = 0
    ENDIF
!
!-- Indices to locate surface (model-top) element
    ALLOCATE ( surfaces%i(1:surfaces%ns) )
    ALLOCATE ( surfaces%j(1:surfaces%ns) )
    ALLOCATE ( surfaces%k(1:surfaces%ns) )

    ALLOCATE ( surfaces%ioff(1:surfaces%ns) )
    ALLOCATE ( surfaces%joff(1:surfaces%ns) )
    ALLOCATE ( surfaces%koff(1:surfaces%ns) )
!
!-- Vertical momentum fluxes of u and v
    ALLOCATE ( surfaces%usws(1:surfaces%ns) )
    ALLOCATE ( surfaces%vsws(1:surfaces%ns) )
!
!-- Sensible heat flux
    ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!-- Latent heat flux
    IF ( humidity  .OR.  ocean_run_coupled_to_atmosphere )  THEN
       ALLOCATE ( surfaces%qsws(1:surfaces%ns) )
    ENDIF
!
!-- Scalar flux
    IF ( passive_scalar )  THEN
       ALLOCATE ( surfaces%ssws(1:surfaces%ns) )
    ENDIF
!
!-- Chemical species flux
    IF ( air_chemistry )  THEN
       ALLOCATE ( surfaces%cssws(1:nvar,1:surfaces%ns) )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_morrison )  THEN
       ALLOCATE ( surfaces%qcsws(1:surfaces%ns) )
       ALLOCATE ( surfaces%ncsws(1:surfaces%ns) )
    ENDIF
!
!--
    IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_seifert )  THEN
       ALLOCATE ( surfaces%qrsws(1:surfaces%ns) )
       ALLOCATE ( surfaces%nrsws(1:surfaces%ns) )
    ENDIF

!
!--
    IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_ice_phase )  THEN
       ALLOCATE ( surfaces%qisws(1:surfaces%ns) )
       ALLOCATE ( surfaces%nisws(1:surfaces%ns) )
    ENDIF
!
!-- Salinity flux
    IF ( ocean_mode )  ALLOCATE ( surfaces%sasws(1:surfaces%ns) )

 END SUBROUTINE allocate_surface_attributes_top


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit memory for model-top fluxes.
!--------------------------------------------------------------------------------------------------!
#if defined( _OPENACC )
 SUBROUTINE exit_surface_attributes_top( surfaces )

    IMPLICIT NONE

    TYPE(surf_type) ::  surfaces  !< respective surface type

    !$ACC EXIT DATA &
    !$ACC DELETE(surfaces%start_index(nys:nyn,nxl:nxr)) &
    !$ACC DELETE(surfaces%end_index(nys:nyn,nxl:nxr)) &
    !$ACC DELETE(surfaces%i(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%j(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%k(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%ioff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%joff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%koff(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%usws(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%vsws(1:surfaces%ns)) &
    !$ACC DELETE(surfaces%shf(1:surfaces%ns)) IF(enable_openacc)

 END SUBROUTINE exit_surface_attributes_top
#endif

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter memory for model-top fluxes.
!--------------------------------------------------------------------------------------------------!
#if defined( _OPENACC )
 SUBROUTINE enter_surface_attributes_top( surfaces )

    IMPLICIT NONE

    TYPE(surf_type) ::  surfaces  !< respective surface type

    !$ACC ENTER DATA &
    !$ACC COPYIN(surfaces%start_index(nys:nyn,nxl:nxr)) &
    !$ACC COPYIN(surfaces%end_index(nys:nyn,nxl:nxr)) &
    !$ACC COPYIN(surfaces%i(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%j(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%k(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%ioff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%joff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%koff(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%usws(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%vsws(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%qsws(1:surfaces%ns)) &
    !$ACC COPYIN(surfaces%shf(1:surfaces%ns)) IF(enable_openacc)

 END SUBROUTINE enter_surface_attributes_top
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize surface elements, i.e. set initial values for surface fluxes, friction velocity,
!> calcuation of start/end indices, etc. Please note, further initialization concerning special
!> surface characteristics, e.g. soil- and vegatation type, building type, etc.,
!> is done in the land-surface, urban-surface, or radiation model, or other modules.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_surfaces

    IMPLICIT NONE

    INTEGER(iwp) ::  i         !< running index x-direction
    INTEGER(iwp) ::  ioff      !< running index x-direction with offset value
    INTEGER(iwp) ::  j         !< running index y-direction
    INTEGER(iwp) ::  joff      !< running index y-direction with offset value
    INTEGER(iwp) ::  k         !< running index z-direction
    INTEGER(iwp) ::  koff      !< running index z-direction with offset value
    INTEGER(iwp) ::  kt_surf   !< uppermost upward-facing surface element across all surface types
    INTEGER(iwp) ::  l         !< index variable for surface facing
    INTEGER(iwp) ::  m         !< running index over surface elements

    INTEGER(iwp) ::  start_index_def !< dummy to determing local start index in surface type for given (j,i), for default surfaces
    INTEGER(iwp) ::  start_index_lsm !< dummy to determing local start index in surface type for given (j,i), for natural surfaces
    INTEGER(iwp) ::  start_index_top !< dummy to determing local start index in surface type for given (j,i), for top surfaces
    INTEGER(iwp) ::  start_index_usm !< dummy to determing local start index in surface type for given (j,i), for urban surfaces
    INTEGER(iwp) ::  start_index_u   !< dummy to determing local start index in surface type for given (j,i), for u-grid surfaces
    INTEGER(iwp) ::  start_index_v   !< dummy to determing local start index in surface type for given (j,i), for v-grid surfaces
    INTEGER(iwp) ::  start_index_w   !< dummy to determing local start index in surface type for given (j,i), for w-grid surfaces
    INTEGER(iwp) ::  num_def         !< current number of surface element, default type
    INTEGER(iwp) ::  num_def_kji     !< dummy to determing local end index in surface type for given (j,i), for default surfaces
    INTEGER(iwp) ::  num_lsm         !< current number of surface element, natural type
    INTEGER(iwp) ::  num_lsm_kji     !< dummy to determing local end index in surface type for given (j,i), for natural surfaces
    INTEGER(iwp) ::  num_top         !< current number of model-top grid point
    INTEGER(iwp) ::  num_top_kji     !< dummy to determing local end index in surface type for given (j,i), for top surfaces
    INTEGER(iwp) ::  num_usm         !< current number of surface element, urban type
    INTEGER(iwp) ::  num_usm_kji     !< dummy to determing local end index in surface type for given (j,i), for urban surfaces
    INTEGER(iwp) ::  num_u           !< current number of surface element, for u-surfaces
    INTEGER(iwp) ::  num_u_kji       !< dummy to determing local end index in surface type for given (j,i), for u-surfaces
    INTEGER(iwp) ::  num_v           !< current number of surface element, for v-surfaces
    INTEGER(iwp) ::  num_v_kji       !< dummy to determing local end index in surface type for given (j,i), for v-surfaces
    INTEGER(iwp) ::  num_w           !< current number of surface element, for w-surfaces
    INTEGER(iwp) ::  num_w_kji       !< dummy to determing local end index in surface type for given (j,i), for w-surfaces


    LOGICAL ::  building             !< flag indicating building grid point
    LOGICAL ::  skip_init_upward_top !< flag to control initialization of upward-top arrays
    LOGICAL ::  terrain              !< flag indicating natural terrain grid point
    LOGICAL ::  unresolved_building  !< flag indicating a grid point where actually a building is defined but not resolved by the grid

!
!-- Initialize surface attributes, store indicies, surfaces orientation, etc.,
    num_def = 1
    num_lsm = 1
    num_top = 1
    num_usm = 1

    IF ( cut_cell_topography )  THEN
       num_u = 1
       num_v = 1
       num_w = 1
    ENDIF

    start_index_def = 1
    start_index_lsm = 1
    start_index_top = 1
    start_index_usm = 1

    IF ( cut_cell_topography )  THEN
       start_index_u = 1
       start_index_v = 1
       start_index_w = 1
    ENDIF

    DO  i = nxl, nxr
       DO  j = nys, nyn
          num_def_kji = 0
          num_lsm_kji = 0
          num_top_kji = 0
          num_usm_kji = 0

          IF ( cut_cell_topography )  THEN
             num_u_kji = 0
             num_v_kji = 0
             num_w_kji = 0
          ENDIF

          DO  k = nzb+1, nzt
!
!--          Check if s-grid point belongs to the atmosphere
             IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                DO  l = 0, 6
                   koff = k + off_k(l)
                   joff = j + off_j(l)
                   ioff = i + off_i(l)
!
!--                Initialize attributes for model top first
                   IF ( koff == nzt+1  .AND.  l == 6  .AND.  use_top_fluxes )  THEN
                      CALL initialize_top( k, j, i, surf_top, num_top, num_top_kji )
                   ELSEIF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( topo_flags(koff,joff,ioff), 5 )  .OR.  topo_no_distinct
                      building = BTEST( topo_flags(koff,joff,ioff), 6 )  .OR.  topo_no_distinct

                      unresolved_building = BTEST( topo_flags(koff,joff,ioff), 5 )  .AND.          &
                                            BTEST( topo_flags(koff,joff,ioff), 6 )  .AND.          &
                                            urban_surface

                      IF ( land_surface  .AND.  terrain  .AND.  .NOT. unresolved_building )  THEN
                         CALL initialize_surfaces( k, j, i, koff, joff, ioff, surf_lsm, num_lsm,   &
                                                   num_lsm_kji )

                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         CALL initialize_surfaces( k, j, i, koff, joff, ioff, surf_usm, num_usm,   &
                                                   num_usm_kji )
                      ELSE
                         CALL initialize_surfaces( k, j, i, koff, joff, ioff, surf_def, num_def,   &
                                                   num_def_kji )
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF

!
!--          Treat surfaces relevant on staggered velocity grid.
             IF ( cut_cell_topography )  THEN
!
!--             u-grid.
                IF ( BTEST( topo_flags(k,j,i), 1 ) )  THEN
                   DO  l = 0, 5
                      koff = k + off_k(l)
                      joff = j + off_j(l)
                      ioff = i + off_i(l)
!
!--                   Initialize attributes for model top first
                      IF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 1 ) )  THEN
                         CALL initialize_surfaces_mom( k, j, i, koff, joff, ioff, surf_u, num_u,   &
                                                       num_u_kji )
                      ENDIF
                   ENDDO
                ENDIF
!
!--             v-grid
                IF ( BTEST( topo_flags(k,j,i), 2 ) )  THEN
                   DO  l = 0, 5
                      koff = k + off_k(l)
                      joff = j + off_j(l)
                      ioff = i + off_i(l)
!
!--                   Initialize attributes for model top first
                      IF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 2 ) )  THEN
                         CALL initialize_surfaces_mom( k, j, i, koff, joff, ioff, surf_v, num_v,   &
                                                       num_v_kji )
                      ENDIF
                   ENDDO
                ENDIF
!
!--             w-grid
                IF ( BTEST( topo_flags(k,j,i), 3 ) )  THEN
                   DO  l = 0, 5
                      koff = k + off_k(l)
                      joff = j + off_j(l)
                      ioff = i + off_i(l)
!
!--                   Initialize attributes for model top first
                      IF ( .NOT. BTEST( topo_flags(koff,joff,ioff), 3 ) )  THEN
                         CALL initialize_surfaces_mom( k, j, i, koff, joff, ioff, surf_w, num_w,   &
                                                       num_w_kji )
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF

          ENDDO
!
!--       Determine start- and end-index at grid point (j,i). Also, for horizontal surfaces more
!--       than 1 horizontal surface element can exist at grid point (j,i) if overhanging structures
!--       are present.
!--       Default surfaces
          surf_def%start_index(j,i) = start_index_def
          surf_def%end_index(j,i)   = surf_def%start_index(j,i) + num_def_kji - 1
          start_index_def           = surf_def%end_index(j,i) + 1
!
!--       Land surfaces
          surf_lsm%start_index(j,i) = start_index_lsm
          surf_lsm%end_index(j,i)   = surf_lsm%start_index(j,i) + num_lsm_kji - 1
          start_index_lsm           = surf_lsm%end_index(j,i) + 1
!
!--       Model top
          surf_top%start_index(j,i) = start_index_top
          surf_top%end_index(j,i)   = surf_top%start_index(j,i) + num_top_kji - 1
          start_index_top           = surf_top%end_index(j,i) + 1
!
!--       Building surfaces
          surf_usm%start_index(j,i) = start_index_usm
          surf_usm%end_index(j,i)   = surf_usm%start_index(j,i) + num_usm_kji - 1
          start_index_usm           = surf_usm%end_index(j,i) + 1

          IF ( cut_cell_topography )  THEN
             surf_u%start_index(j,i) = start_index_u
             surf_u%end_index(j,i)   = surf_u%start_index(j,i) + num_u_kji - 1
             start_index_u           = surf_u%end_index(j,i) + 1

             surf_v%start_index(j,i) = start_index_v
             surf_v%end_index(j,i)   = surf_v%start_index(j,i) + num_v_kji - 1
             start_index_v           = surf_v%end_index(j,i) + 1

             surf_w%start_index(j,i) = start_index_w
             surf_w%end_index(j,i)   = surf_w%start_index(j,i) + num_w_kji - 1
             start_index_w           = surf_w%end_index(j,i) + 1
          ENDIF
!
!--       ATTENTION:
!--       Workaround to prevent vectorization bug on NEC Aurora
          IF ( start_index_def < -99999 )  THEN
             PRINT*, 'i=', i, ' j=',j, ' s=',surf_def%start_index(j,i), ' e=', surf_def%end_index(j,i)
          ENDIF

       ENDDO
    ENDDO

!
!-- Finally, determine the uppermost upward-facing surface element from all types at given
!-- (j,i)-grid point. Later-on, this information is e.g. used to simplify the temporal averaging of
!-- surface data. First, allocate the arrays. Check if arrays have been already allocated.
!-- This is the situation with cyclic-fill initialization, where this routine is called twice.
    skip_init_upward_top = .FALSE.
    IF ( surf_def%ns >= 1 )  THEN
       IF ( .NOT. ALLOCATED( surf_def%upward_top ) )  THEN
          ALLOCATE( surf_def%upward_top(1:surf_def%ns) )
          surf_def%upward_top(:) = .FALSE.
       ELSE
          skip_init_upward_top = .TRUE.
       ENDIF
    ENDIF
    IF ( surf_lsm%ns >= 1 )  THEN
       IF ( .NOT. ALLOCATED( surf_lsm%upward_top ) )  THEN
          ALLOCATE( surf_lsm%upward_top(1:surf_lsm%ns) )
          surf_lsm%upward_top(:) = .FALSE.
       ELSE
          skip_init_upward_top = .TRUE.
       ENDIF
    ENDIF
    IF ( surf_usm%ns >= 1 )  THEN
       IF ( .NOT. ALLOCATED( surf_usm%upward_top ) )  THEN
          ALLOCATE( surf_usm%upward_top(1:surf_usm%ns) )
          surf_usm%upward_top(:) = .FALSE.
       ELSE
          skip_init_upward_top = .TRUE.
       ENDIF
    ENDIF
!
!-- Inititialize upward_top flags. Skip this if the flags have been already set.
    IF ( .NOT. skip_init_upward_top )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          First, determine the vertical index of the uppermost upward-facing surface from
!--          all surface types.
             kt_surf = -HUGE( 1_iwp )
             DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                IF ( surf_def%upward(m)  .AND.  ( surf_def%k(m) + surf_def%koff(m) ) > kt_surf )   &
                THEN
                   kt_surf = surf_def%k(m) + surf_def%koff(m)
                ENDIF
             ENDDO
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m)  .AND.  ( surf_lsm%k(m) + surf_lsm%koff(m) ) > kt_surf )   &
                THEN
                   kt_surf = surf_lsm%k(m) + surf_lsm%koff(m)
                ENDIF
             ENDDO
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                IF ( surf_usm%upward(m)  .AND.  ( surf_usm%k(m) + surf_usm%koff(m) ) > kt_surf )   &
                THEN
                   kt_surf = surf_usm%k(m) + surf_usm%koff(m)
                ENDIF
             ENDDO
!
!--          Based on the previously calculated surface-top index, set the flag accordingly. At a
!--          given (j,i)-index, there is maximum one entry that is .True..
             DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                IF ( surf_def%upward(m)  .AND.  ( surf_def%k(m) + surf_def%koff(m) ) == kt_surf )  &
                THEN
                   surf_def%upward_top(m) = .TRUE.
                ENDIF
             ENDDO
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m)  .AND.  ( surf_lsm%k(m) + surf_lsm%koff(m) ) == kt_surf )  &
                THEN
                   surf_lsm%upward_top(m) = .TRUE.
                ENDIF
             ENDDO
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                IF ( surf_usm%upward(m)  .AND.  ( surf_usm%k(m) + surf_usm%koff(m) ) == kt_surf )  &
                THEN
                   surf_usm%upward_top(m) = .TRUE.
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE init_surfaces

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal surface elements, upward- and downward-facing. Note, horizontal surface
!> type also comprises model-top fluxes, which are, initialized in a different routine.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE initialize_surfaces( k, j, i, k_surf, j_surf, i_surf, surf, num, num_kji )

    IMPLICIT NONE

    INTEGER(iwp) ::  component  !< index of wall_fluxes_ array for respective orientation
    INTEGER(iwp) ::  i          !< x-index of reference grid point
    INTEGER(iwp) ::  i_surf     !< x-index of surface grid point
    INTEGER(iwp) ::  j          !< y-index of reference grid point
    INTEGER(iwp) ::  j_surf     !< y-index of surface grid point
    INTEGER(iwp) ::  k          !< z-index of reference grid point
    INTEGER(iwp) ::  k_surf     !< z-index of surface grid point
    INTEGER(iwp) ::  num        !< current number of surface element
    INTEGER(iwp) ::  num_kji    !< dummy increment
    INTEGER(iwp) ::  lsp        !< running index chemical species
    INTEGER(iwp) ::  lsp_pr     !< running index chemical species??

    TYPE(surf_type) ::  surf     !< respective surface type

!
!-- Store indices of respective surface element and its corresponding reference grid point in
!-- the atmosphere.
    surf%i(num) = i
    surf%j(num) = j
    surf%k(num) = k

    surf%iref(num) = i
    surf%jref(num) = j
    surf%kref(num) = k

    surf%ioff(num) = i_surf - i
    surf%joff(num) = j_surf - j
    surf%koff(num) = k_surf - k
!
!-- Set logical flags indicating the main facing of the surface.
!-- Though for cut cells the orientation might not entirely coincide with
!-- the grid axis, the grid cell is still bounded towards this wall.
    surf%upward(num)    = surf%koff(num) == -1
    surf%downward(num)  = surf%koff(num) ==  1
    surf%northward(num) = surf%joff(num) == -1
    surf%southward(num) = surf%joff(num) ==  1
    surf%eastward(num)  = surf%ioff(num) == -1
    surf%westward(num)  = surf%ioff(num) ==  1

    surf%cut_cell_wall(num) = .FALSE.
    surf%cut_cell_roof(num) = .FALSE.

    surf%bid(num)   = -HUGE( 1_iwp )
    surf%btype(num) = 1
    surf%ptype(num) = -HUGE( 1_iwp )
    surf%vtype(num) = 1
    surf%wtype(num) = -HUGE( 1_iwp )
!
!-- Initialize attributes depending on orientation. In the following, component is set
!-- to properly set wall_heatflux.
    IF ( surf%upward(num) )  THEN
       surf%z_mo(num) = zu(k) - zw(k-1)
       component = 0
    ENDIF
    IF ( surf%downward(num) )  THEN
       surf%z_mo(num) = zw(k) - zu(k)
       component = 5
    ENDIF

    IF ( surf%northward(num)  .OR.  surf%southward(num) )  THEN
       surf%z_mo(num) = 0.5_wp * dy
       IF ( surf%northward(num) )  component = 4
       IF ( surf%southward(num) )  component = 3
    ELSEIF ( surf%eastward(num)  .OR.  surf%westward(num) )  THEN
       surf%z_mo(num) = 0.5_wp * dx
       IF ( surf%eastward(num) )  component = 2
       IF ( surf%westward(num) )  component = 1
    ENDIF

    IF ( surf%upward(num) )  THEN
       surf%consider_stability(num) = .TRUE.
    ELSE
       surf%consider_stability(num) = .FALSE.
    ENDIF

    surf%tke_production(num) = .TRUE.
!
!-- Calculate normal vector. In this first step set normal vector as usual for step-like surfaces.
!-- If cut-cell topography is employed, overwrite this at a later point.
    CALL calculate_surface_orientation( surf, num )
!
!-- Initialize roughness length.
    surf%z0(num)  = roughness_length
    surf%z0h(num) = z0h_factor * roughness_length
    surf%z0q(num) = z0h_factor * roughness_length
!
!-- Initialization in case of 1D pre-cursor run
    IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       IF ( .NOT. constant_diffusion )  THEN
          IF ( constant_flux_layer  .AND.  surf%upward(num) )  THEN
             surf%ol(num)        = surf%z_mo(num) / ( ri1d(nzb+1) + 1.0E-20_wp )
             surf%us(num)        = us1d
             surf%us_uvgrid(num) = us1d
             surf%us_wgrid(num)  = us1d
             surf%usws(num)      = usws1d
             surf%vsws(num)      = vsws1d
          ELSE
             surf%ol(num)        = surf%z_mo(num) / zeta_min
             surf%us(num)        = 0.0_wp
             surf%us_uvgrid(num) = 0.0_wp
             surf%us_wgrid(num)  = 0.0_wp
             surf%usws(num)      = 0.0_wp
             surf%vsws(num)      = 0.0_wp
          ENDIF
       ELSE
          surf%ol(num)        = surf%z_mo(num) / zeta_min
          surf%us(num)        = 0.0_wp
          surf%us_uvgrid(num) = 0.0_wp
          surf%us_wgrid(num)  = 0.0_wp
          surf%usws(num)      = 0.0_wp
          surf%vsws(num)      = 0.0_wp
       ENDIF
!
!-- Initialization in all other cases
    ELSE

       surf%ol(num) = surf%z_mo(num) / zeta_min
!
!--    Very small number is required for calculation of Obukhov length at first timestep
       surf%us(num)        = 1.0E-30_wp
       surf%us_uvgrid(num) = 1.0E-30_wp
       surf%us_wgrid(num)  = 1.0E-30_wp
       surf%usws(num)      = 0.0_wp
       surf%vsws(num)      = 0.0_wp

    ENDIF

    surf%usvs(num) = 0.0_wp
    surf%vsus(num) = 0.0_wp
    surf%wsus_wsvs(num) = 0.0_wp
    surf%mom_flux_tke(0:1,num) = 0.0_wp

    surf%rib(num)        = 0.0_wp
    surf%uvw_abs(num)    = 0.0_wp
    surf%uvw_abs_uv(num) = 0.0_wp
    surf%uvw_abs_w(num)  = 0.0_wp
!
!-- Initialize ln(z/z0)
    surf%ln_z_z0(num)  = LOG( surf%z_mo(num) / surf%z0(num) )
    surf%ln_z_z0h(num) = LOG( surf%z_mo(num) / surf%z0h(num) )
    surf%ln_z_z0q(num) = LOG( surf%z_mo(num) / surf%z0q(num) )

    IF ( .NOT. constant_diffusion )  THEN
       surf%u_0(num) = 0.0_wp
       surf%v_0(num) = 0.0_wp
    ENDIF

    surf%ts(num) = 0.0_wp
!
!-- Set initial value for surface temperature
    surf%pt_surface(num) = pt_surface

    IF ( humidity )  THEN
       surf%qs(num)   = 0.0_wp
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          surf%qcs(num) = 0.0_wp
          surf%ncs(num) = 0.0_wp

          surf%qcsws(num) = 0.0_wp
          surf%ncsws(num) = 0.0_wp

       ENDIF
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          surf%qrs(num) = 0.0_wp
          surf%nrs(num) = 0.0_wp

          surf%qrsws(num) = 0.0_wp
          surf%nrsws(num) = 0.0_wp

          surf%pt1(num)  = 0.0_wp
          surf%qv1(num)  = 0.0_wp
          surf%vpt1(num) = 0.0_wp

       ENDIF

       IF ( surf_bulk_cloud_model .AND. surf_microphysics_ice_phase)  THEN
          surf%qis(num) = 0.0_wp
          surf%nis(num) = 0.0_wp

          surf%qisws(num) = 0.0_wp
          surf%nisws(num) = 0.0_wp
       ENDIF

       surf%q_surface(num)   = q_surface
       surf%vpt_surface(num) = surf%pt_surface(num) * ( 1.0_wp + 0.61_wp * surf%q_surface(num) )
    ENDIF

    IF ( passive_scalar )  surf%ss(num) = 0.0_wp

    IF ( air_chemistry )  THEN
       DO  lsp = 1, nvar
          surf%css(lsp,num)   = 0.0_wp
          surf%cssws(lsp,num) = 0.0_wp
       ENDDO
    ENDIF
!
!-- Inititalize surface fluxes of sensible and latent heat, as well as passive scalar
    surf%shf(num) = 0.0_wp
    IF ( use_surface_fluxes )  THEN

       IF ( surf%upward(num) )  THEN
          IF ( constant_heatflux )  THEN
!
!--          Initialize surface heatflux. However, skip this for now if random_heatflux is set.
!--          This case, shf is initialized later.
             IF ( .NOT. random_heatflux )  THEN
                surf%shf(num) = surface_heatflux * heatflux_input_conversion(k+surf%koff(num))
!
!--             Check if surface heat flux might be replaced by prescribed wall heatflux
                IF ( k-1 /= 0 )  THEN
                   surf%shf(num) = wall_heatflux(component) * heatflux_input_conversion(k+surf%koff(num))
                ENDIF
             ENDIF
          ELSE
             surf%shf(num) = 0.0_wp
          ENDIF
!
!--    Set heat-flux at downward-facing surfaces
       ELSEIF ( surf%downward(num) )  THEN
          surf%shf(num) = wall_heatflux(component) * heatflux_input_conversion(k)
!
!--    Set heat-flux at vertical surfaces, no multiplication with density here.
       ELSE
          surf%shf(num) = wall_heatflux(component)
       ENDIF

       IF ( humidity )  THEN
          surf%qsws(num) = 0.0_wp
          IF ( surf%upward(num) )  THEN
             IF ( constant_waterflux )  THEN
                surf%qsws(num) = surface_waterflux * waterflux_input_conversion(k+surf%koff(num))
                IF ( k-1 /= 0 )  THEN
                   surf%qsws(num) = wall_humidityflux(component) * waterflux_input_conversion(k+surf%koff(num))
                ENDIF
             ELSE
                surf%qsws(num) = 0.0_wp
             ENDIF
          ELSEIF ( surf%downward(num) )  THEN
             surf%qsws(num) = wall_humidityflux(component) * waterflux_input_conversion(k)
!
!--       Set heat-flux at vertical surfaces, no multiplication with density here.
          ELSE
             surf%qsws(num) = wall_humidityflux(component)
          ENDIF
       ENDIF

       IF ( passive_scalar )  THEN
          surf%ssws(num) = 0.0_wp
          IF ( surf%upward(num) )  THEN
             IF ( constant_scalarflux )  THEN
                surf%ssws(num) = surface_scalarflux  * scalarflux_input_conversion(k+surf%koff(num))

                IF ( k-1 /= 0 )  surf%ssws(num) = wall_scalarflux(component)                       &
                                                * scalarflux_input_conversion(k+surf%koff(num))
             ELSE
                surf%ssws(num) = 0.0_wp
             ENDIF
          ELSEIF ( surf%downward(num) )  THEN
             surf%ssws(num) = wall_scalarflux(component) * scalarflux_input_conversion(k)
!
!--       Set flux at vertical surfaces, no multiplication with density here.
          ELSE
             surf%ssws(num) = wall_scalarflux(component)
          ENDIF
       ENDIF

       IF ( air_chemistry )  THEN
          lsp_pr = 1
          DO  WHILE ( TRIM( surface_csflux_name( lsp_pr ) ) /= 'novalue' ) !<'novalue' is the default
             DO  lsp = 1, nvar
!
!--             Assign surface flux for each variable species
                IF ( TRIM( spc_names(lsp) ) == TRIM( surface_csflux_name(lsp_pr) ) )  THEN
                   IF ( surf%upward(num) )  THEN
                      IF ( constant_csflux(lsp_pr) )  THEN
                         surf%cssws(lsp,num) = surface_csflux(lsp_pr) *                            &
                                               scalarflux_input_conversion(k+surf%koff(num))

                         IF ( k-1 /= 0 )  surf%cssws(lsp,num) = wall_csflux(lsp,component) *       &
                                                       scalarflux_input_conversion(k+surf%koff(num))
                      ELSE
                         surf%cssws(lsp,num) = 0.0_wp
                      ENDIF
                   ELSEIF ( surf%downward(num) )  THEN
                      surf%cssws(lsp,num) = wall_csflux(lsp,component) *                           &
                                            scalarflux_input_conversion(k)
!
!--                Set flux at vertical surfaces, no multiplication with density here.
                   ELSE
                      surf%cssws(lsp,num) = wall_csflux(lsp,component)
                   ENDIF
                ENDIF
             ENDDO
             lsp_pr = lsp_pr + 1
          ENDDO
       ENDIF

       IF ( ocean_mode )  THEN
          IF ( surf%upward(num) )  THEN
             surf%sasws(num) = bottom_salinityflux * scalarflux_input_conversion(k+surf%koff(num))
          ELSE
             surf%sasws(num) = 0.0_wp
          ENDIF
       ENDIF

    ENDIF
!
!-- Increment surface indices
    num     = num + 1
    num_kji = num_kji + 1


 END SUBROUTINE initialize_surfaces


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize surfaces considered on the staggered velocity grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE initialize_surfaces_mom( k, j, i, k_surf, j_surf, i_surf, surf, num, num_kji )

    IMPLICIT NONE

    INTEGER(iwp) ::  component  !< index of wall_fluxes_ array for respective orientation
    INTEGER(iwp) ::  i          !< x-index of reference grid point
    INTEGER(iwp) ::  i_surf     !< x-index of surface grid point
    INTEGER(iwp) ::  j          !< y-index of reference grid point
    INTEGER(iwp) ::  j_surf     !< y-index of surface grid point
    INTEGER(iwp) ::  k          !< z-index of reference grid point
    INTEGER(iwp) ::  k_surf     !< z-index of surface grid point
    INTEGER(iwp) ::  num        !< current number of surface element
    INTEGER(iwp) ::  num_kji    !< dummy increment

    TYPE(surf_type) ::  surf     !< respective surface type

!
!-- Store indices of respective surface element and its corresponding reference grid point in
!-- the atmosphere.
    surf%i(num) = i
    surf%j(num) = j
    surf%k(num) = k

    surf%iref(num) = i
    surf%jref(num) = j
    surf%kref(num) = k

    surf%ioff(num) = i_surf - i
    surf%joff(num) = j_surf - j
    surf%koff(num) = k_surf - k
!
!-- Set logical flags indicating the main facing of the surface.
!-- Though for cut cells the orientation might not entirely coincide with
!-- the grid axis, the grid cell is still bounded towards this wall.
    surf%upward(num)    = surf%koff(num) == -1
    surf%downward(num)  = surf%koff(num) ==  1
    surf%northward(num) = surf%joff(num) == -1
    surf%southward(num) = surf%joff(num) ==  1
    surf%eastward(num)  = surf%ioff(num) == -1
    surf%westward(num)  = surf%ioff(num) ==  1
!
!-- Initialize attributes depending on orientation. In the following, component is set
!-- to properly set wall_heatflux.
    IF ( surf%upward(num) )  THEN
       surf%z_mo(num) = zu(k) - zw(k-1)
    ENDIF
    IF ( surf%downward(num) )  THEN
       surf%z_mo(num) = zw(k) - zu(k)
       component = 5
    ENDIF

    IF ( surf%northward(num)  .OR.  surf%southward(num) )  THEN
       surf%z_mo(num) = 0.5_wp * dy
       IF ( surf%northward(num) )  component = 4
       IF ( surf%southward(num) )  component = 3
    ELSEIF ( surf%eastward(num)  .OR.  surf%westward(num) )  THEN
       surf%z_mo(num) = 0.5_wp * dx
       IF ( surf%eastward(num) )  component = 2
       IF ( surf%westward(num) )  component = 1
    ENDIF
!
!-- Calculate normal vector. In this first step set normal vector as usual for step-like surfaces.
!-- If cut-cell topography is employed, overwrite this at a later point.
    CALL calculate_surface_orientation( surf, num )

    IF ( surf%upward(num) )  THEN
       surf%consider_stability(num) = .TRUE.
    ELSE
       surf%consider_stability(num) = .FALSE.
    ENDIF

!
!-- Initialize roughness length.
    surf%z0(num)  = roughness_length
!
!-- Initialization in case of 1D pre-cursor run
    IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       IF ( .NOT. constant_diffusion )  THEN
          IF ( constant_flux_layer  .AND.  surf%upward(num) )  THEN
             surf%ol(num)        = surf%z_mo(num) / ( ri1d(nzb+1) + 1.0E-20_wp )
             surf%us(num)        = us1d
             surf%us_uvgrid(num) = us1d
             surf%us_wgrid(num)  = us1d
             surf%usws(num)      = usws1d
             surf%vsws(num)      = vsws1d
          ELSE
             surf%ol(num)        = surf%z_mo(num) / zeta_min
             surf%us(num)        = 0.0_wp
             surf%us_uvgrid(num) = 0.0_wp
             surf%us_wgrid(num)  = 0.0_wp
             surf%usws(num)      = 0.0_wp
             surf%vsws(num)      = 0.0_wp
          ENDIF
       ELSE
          surf%ol(num)        = surf%z_mo(num) / zeta_min
          surf%us(num)        = 0.0_wp
          surf%us_uvgrid(num) = 0.0_wp
          surf%us_wgrid(num)  = 0.0_wp
          surf%usws(num)      = 0.0_wp
          surf%vsws(num)      = 0.0_wp
       ENDIF
!
!-- Initialization in all other cases
    ELSE

       surf%ol(num) = surf%z_mo(num) / zeta_min
!
!--    Very small number is required for calculation of Obukhov length at first timestep
       surf%us(num)        = 1.0E-30_wp
       surf%us_uvgrid(num) = 1.0E-30_wp
       surf%us_wgrid(num)  = 1.0E-30_wp
       surf%usws(num)      = 0.0_wp
       surf%vsws(num)      = 0.0_wp

    ENDIF

    surf%usvs(num) = 0.0_wp
    surf%vsus(num) = 0.0_wp
    surf%wsus_wsvs(num) = 0.0_wp

    surf%uvw_abs(num)    = 0.0_wp
    surf%uvw_abs_uv(num) = 0.0_wp
    surf%uvw_abs_w(num)  = 0.0_wp
!
!-- Initialize ln(z/z0)
    surf%ln_z_z0(num)  = LOG( surf%z_mo(num) / surf%z0(num) )
!
!-- Increment surface indices
    num     = num + 1
    num_kji = num_kji + 1


 END SUBROUTINE initialize_surfaces_mom


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize model-top fluxes. Currently, only the heatflux and salinity flux can be prescribed,
!> latent flux is zero in this case!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE initialize_top( k, j, i, surf, num_h, num_h_kji )

    IMPLICIT NONE

    INTEGER(iwp)  ::  i          !< running index x-direction
    INTEGER(iwp)  ::  j          !< running index y-direction
    INTEGER(iwp)  ::  k          !< running index z-direction
    INTEGER(iwp)  ::  num_h      !< current number of surface element
    INTEGER(iwp)  ::  num_h_kji  !< dummy increment

    TYPE( surf_type ) ::  surf   !< respective surface type
!
!-- Store indices of respective surface element
    surf%i(num_h) = i
    surf%j(num_h) = j
    surf%k(num_h) = k

    surf%ioff(num_h) = 0
    surf%joff(num_h) = 0
    surf%koff(num_h) = 1
!
!-- Initialize top heat flux
    IF ( constant_top_heatflux )  surf%shf(num_h) = top_heatflux * heatflux_input_conversion(nzt+1)
!
!-- Initialization in case of an ocean run coupled to atmosphere.
    IF ( ocean_run_coupled_to_atmosphere )  THEN
       surf%shf(num_h)  = 0.0_wp
       surf%qsws(num_h) = 0.0_wp
    ENDIF
!
!-- Prescribe latent heat flux at the top
    IF ( humidity )  THEN
       surf%qsws(num_h) = 0.0_wp
       IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_morrison )  THEN
          surf%ncsws(num_h) = 0.0_wp
          surf%qcsws(num_h) = 0.0_wp
       ENDIF
       IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_seifert )  THEN
          surf%nrsws(num_h) = 0.0_wp
          surf%qrsws(num_h) = 0.0_wp
       ENDIF
       IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_ice_phase )  THEN
          surf%nisws(num_h) = 0.0_wp
          surf%qisws(num_h) = 0.0_wp
       ENDIF
    ENDIF
!
!-- Prescribe top scalar flux
    IF ( passive_scalar .AND. constant_top_scalarflux )  surf%ssws(num_h) = top_scalarflux *       &
                                                                  scalarflux_input_conversion(nzt+1)
!
!-- Prescribe top salinity flux
    IF ( ocean_mode .AND. constant_top_salinityflux)  surf%sasws(num_h) = top_salinityflux *       &
                                                                  scalarflux_input_conversion(nzt+1)
!
!-- Top momentum fluxes
    IF ( constant_top_momentumflux )  THEN
       surf%usws(num_h) = top_momentumflux_u * momentumflux_input_conversion(nzt+1)
       surf%vsws(num_h) = top_momentumflux_v * momentumflux_input_conversion(nzt+1)
    ENDIF
!
!-- Increment surface indices
    num_h     = num_h + 1
    num_h_kji = num_h_kji + 1


 END SUBROUTINE initialize_top


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Determine the normalized normal vector of the surface as well as the distance to the wall.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calculate_surface_orientation( surf, num )

    INTEGER(iwp)      ::  num  !< surface index

    TYPE( surf_type ) ::  surf !< treated surface array

    surf%n_s(num,:) = 0.0_wp
    IF ( surf%koff(num) /= 0 )  THEN
       surf%n_s(num,1) = MERGE( 1.0_wp, -1.0_wp, surf%koff(num) == -1 )
       surf%n_eff(num) = surf%n_s(num,1)
    ELSEIF( surf%ioff(num) /= 0 )  THEN
       surf%n_s(num,3) = MERGE( 1.0_wp, -1.0_wp, surf%ioff(num) == -1 )
       surf%n_eff(num) = surf%n_s(num,3)
    ELSEIF( surf%joff(num) /= 0 )  THEN
       surf%n_s(num,2) = MERGE( 1.0_wp, -1.0_wp, surf%joff(num) == -1 )
       surf%n_eff(num) = surf%n_s(num,2)
    ENDIF

 END SUBROUTINE calculate_surface_orientation


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize single surface properties from 2D input arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_single_surface_properties( var_surf, var_2d, ns, fill_value, i_index, i_offset,   &
                                            j_index, j_offset, k_index, k_offset,                  &
                                            input_conversion_factor )

    INTEGER(iwp) ::  m   !< running index over surface elements
    INTEGER(iwp) ::  ns  !< number of surface elements in var_surf

    INTEGER(iwp), DIMENSION(1:ns), INTENT(IN) ::  i_index   !< surface grid index along x
    INTEGER(iwp), DIMENSION(1:ns), INTENT(IN) ::  i_offset  !< surfacd gitd index offset along x
    INTEGER(iwp), DIMENSION(1:ns), INTENT(IN) ::  j_index   !< surface grid index along y
    INTEGER(iwp), DIMENSION(1:ns), INTENT(IN) ::  j_offset  !< surfacd gitd index offset along y
    INTEGER(iwp), DIMENSION(1:ns), INTENT(IN) ::  k_index   !< surface grid index along z
    INTEGER(iwp), DIMENSION(1:ns), INTENT(IN) ::  k_offset  !< surfacd gitd index offset along z

    REAL(wp), INTENT(IN) ::  fill_value !< fill value in var_2d

    REAL(wp), DIMENSION(nzb:nzt+1), OPTIONAL ::  input_conversion_factor !< optional profile of input flux conversion factor
                                                                         !< for an unknow reason the GNU compiler has interface
                                                                         !< problems when INTENT(IN) is added here

    REAL(wp), DIMENSION(1:ns), INTENT(INOUT) ::  var_surf  !< 1D surface variable that should be initialized

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg), INTENT(IN) ::  var_2d  !< input variable

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_space_i  !< grid index along x where surface properties are be defined
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_space_j  !< grid index along y where surface properties are be defined
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  index_space_k  !< grid index along k where surface properties are be defined

    REAL(wp), DIMENSION(nzb:nzt+1) ::  conversion_factor  !< input flux conversion factor to be used


!
!-- Determine surface indices.
    ALLOCATE( index_space_i(1:ns) )
    ALLOCATE( index_space_j(1:ns) )
    ALLOCATE( index_space_k(1:ns) )
    index_space_i(:) = i_index(:) + i_offset(:)
    index_space_j(:) = j_index(:) + j_offset(:)
    index_space_k(:) = k_index(:) + k_offset(:)

    IF ( PRESENT( input_conversion_factor ) )  THEN
       conversion_factor = input_conversion_factor
    ELSE
       conversion_factor = 1.0_wp
    ENDIF

    DO  m = 1, ns
       IF ( var_2d(index_space_j(m),index_space_i(m)) /= fill_value )  THEN
          var_surf(m) = var_2d(index_space_j(m),index_space_i(m))                                  &
                      * conversion_factor(index_space_k(m))
       ENDIF
    ENDDO

    DEALLOCATE( index_space_i, index_space_j, index_space_k )

 END SUBROUTINE init_single_surface_properties


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gathers all surface elements with the same facing (but possibly different type) onto a surface
!> type, and writes binary data into restart files.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_wrd_local


    IMPLICIT NONE

    INTEGER(iwp)                 ::  i              !< running index x-direction
    INTEGER(iwp)                 ::  j              !< running index y-direction
    INTEGER(iwp)                 ::  lsp            !< running index chemical species
    INTEGER(iwp)                 ::  m              !< running index for surface elements on individual surface array
    INTEGER(iwp)                 ::  start_index    !< start index for horizontal surface elements on gathered surface array
    INTEGER(iwp)                 ::  mm             !< running index for surface elements on gathered surface array

    INTEGER(idp),DIMENSION(nys:nyn,nxl:nxr) ::  global_end_index    !< end index for surface data (MPI-IO)
    INTEGER(idp),DIMENSION(nys:nyn,nxl:nxr) ::  global_start_index  !< start index for surface data (MPI-IO)

    LOGICAL ::  surface_data_to_write  !< switch for MPI-I/O if PE has surface data to write

    TYPE(surf_type) ::  surf  !< gathered surfaces, contains all surface types and facings

!
!-- Determine total number of horizontal and vertical surface elements before writing var_list
    CALL surface_last_actions
!
!-- Allocate attributes.
!-- Horizontal upward facing
    surf%ns = ns_on_file(1)
    CALL allocate_surface_attributes( surf, nys, nyn, nxl, nxr )
!
!-- In the following, gather data from surfaces with the same type.
    mm = 1
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
             IF ( ALLOCATED( surf_def%us )  )  surf%us(mm)  = surf_def%us(m)
             IF ( ALLOCATED( surf_def%ts )  )  surf%ts(mm)  = surf_def%ts(m)
             IF ( ALLOCATED( surf_def%qs )  )  surf%qs(mm)  = surf_def%qs(m)
             IF ( ALLOCATED( surf_def%ss )  )  surf%ss(mm)  = surf_def%ss(m)
             IF ( ALLOCATED( surf_def%qcs ) )  surf%qcs(mm) = surf_def%qcs(m)
             IF ( ALLOCATED( surf_def%ncs ) )  surf%ncs(mm) = surf_def%ncs(m)
             IF ( ALLOCATED( surf_def%qis ) )  surf%qis(mm) = surf_def%qis(m)
             IF ( ALLOCATED( surf_def%nis ) )  surf%nis(mm) = surf_def%nis(m)
             IF ( ALLOCATED( surf_def%qrs ) )  surf%qrs(mm) = surf_def%qrs(m)
             IF ( ALLOCATED( surf_def%nrs ) )  surf%nrs(mm) = surf_def%nrs(m)
             IF ( ALLOCATED( surf_def%ol )  )  surf%ol(mm)  = surf_def%ol(m)
             IF ( ALLOCATED( surf_def%rib ) )  surf%rib(mm) = surf_def%rib(m)
             IF ( ALLOCATED( surf_def%pt_surface )  )  surf%pt_surface(mm)  = surf_def%pt_surface(m)
             IF ( ALLOCATED( surf_def%q_surface )   )  surf%q_surface(mm)   = surf_def%q_surface(m)
             IF ( ALLOCATED( surf_def%vpt_surface ) )  surf%vpt_surface(mm) = surf_def%vpt_surface(m)
             IF ( ALLOCATED( surf_def%shf )  )  surf%shf(mm)  = surf_def%shf(m)
             IF ( ALLOCATED( surf_def%qsws ) )  surf%qsws(mm) = surf_def%qsws(m)
             IF ( ALLOCATED( surf_def%ssws ) )  surf%ssws(mm) = surf_def%ssws(m)
             IF ( ALLOCATED( surf_def%css )  )  THEN
                DO  lsp = 1,nvar
                   surf%css(lsp,mm) = surf_def%css(lsp,m)
                ENDDO
             ENDIF
             IF ( ALLOCATED( surf_def%cssws ) )  THEN
                DO  lsp = 1,nvar
                   surf%cssws(lsp,mm) = surf_def%cssws(lsp,m)
                ENDDO
             ENDIF
             IF ( ALLOCATED( surf_def%qcsws ) )  surf%qcsws(mm) = surf_def%qcsws(m)
             IF ( ALLOCATED( surf_def%qrsws ) )  surf%qrsws(mm) = surf_def%qrsws(m)
             IF ( ALLOCATED( surf_def%qisws ) )  surf%qisws(mm) = surf_def%qisws(m)
             IF ( ALLOCATED( surf_def%ncsws ) )  surf%ncsws(mm) = surf_def%ncsws(m)
             IF ( ALLOCATED( surf_def%nisws ) )  surf%nisws(mm) = surf_def%nisws(m)
             IF ( ALLOCATED( surf_def%nrsws ) )  surf%nrsws(mm) = surf_def%nrsws(m)
             IF ( ALLOCATED( surf_def%sasws ) )  surf%sasws(mm) = surf_def%sasws(m)

             IF ( ALLOCATED( surf_def%usws )        )  surf%usws(mm) = surf_def%usws(m)
             IF ( ALLOCATED( surf_def%vsws )        )  surf%vsws(mm) = surf_def%vsws(m)
             IF ( ALLOCATED( surf_def%usvs )        )  surf%usvs(mm) = surf_def%usvs(m)
             IF ( ALLOCATED( surf_def%vsus )        )  surf%vsus(mm) = surf_def%vsus(m)
             IF ( ALLOCATED( surf_def%wsus_wsvs)    )  surf%wsus_wsvs(mm) = surf_def%wsus_wsvs(m)
             IF ( ALLOCATED( surf_def%mom_flux_tke) )                                              &
                surf%mom_flux_tke(0:1,mm) = surf_def%mom_flux_tke(0:1,m)

             mm = mm + 1
          ENDDO

          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             IF ( ALLOCATED( surf_lsm%us )  )  surf%us(mm)  = surf_lsm%us(m)
             IF ( ALLOCATED( surf_lsm%ts )  )  surf%ts(mm)  = surf_lsm%ts(m)
             IF ( ALLOCATED( surf_lsm%qs )  )  surf%qs(mm)  = surf_lsm%qs(m)
             IF ( ALLOCATED( surf_lsm%ss )  )  surf%ss(mm)  = surf_lsm%ss(m)
             IF ( ALLOCATED( surf_lsm%qcs ) )  surf%qcs(mm) = surf_lsm%qcs(m)
             IF ( ALLOCATED( surf_lsm%ncs ) )  surf%ncs(mm) = surf_lsm%ncs(m)
             IF ( ALLOCATED( surf_lsm%qis ) )  surf%qis(mm) = surf_lsm%qis(m)
             IF ( ALLOCATED( surf_lsm%nis ) )  surf%nis(mm) = surf_lsm%nis(m)
             IF ( ALLOCATED( surf_lsm%qrs ) )  surf%qrs(mm) = surf_lsm%qrs(m)
             IF ( ALLOCATED( surf_lsm%nrs ) )  surf%nrs(mm) = surf_lsm%nrs(m)
             IF ( ALLOCATED( surf_lsm%ol )  )  surf%ol(mm)  = surf_lsm%ol(m)
             IF ( ALLOCATED( surf_lsm%rib ) )  surf%rib(mm) = surf_lsm%rib(m)
             IF ( ALLOCATED( surf_lsm%pt_surface ) )   surf%pt_surface(mm)  = surf_lsm%pt_surface(m)
             IF ( ALLOCATED( surf_lsm%q_surface )  )   surf%q_surface(mm)   = surf_lsm%q_surface(m)
             IF ( ALLOCATED( surf_lsm%vpt_surface ) )  surf%vpt_surface(mm) = surf_lsm%vpt_surface(m)
             IF ( ALLOCATED( surf_lsm%shf )  )  surf%shf(mm)  = surf_lsm%shf(m)
             IF ( ALLOCATED( surf_lsm%qsws ) )  surf%qsws(mm) = surf_lsm%qsws(m)
             IF ( ALLOCATED( surf_lsm%ssws ) )  surf%ssws(mm) = surf_lsm%ssws(m)
             IF ( ALLOCATED( surf_lsm%css )  )  THEN
                DO  lsp = 1,nvar
                   surf%css(lsp,mm) = surf_lsm%css(lsp,m)
                ENDDO
             ENDIF
             IF ( ALLOCATED( surf_lsm%cssws ) )  THEN
                DO  lsp = 1, nvar
                   surf%cssws(lsp,mm) = surf_lsm%cssws(lsp,m)
                ENDDO
             ENDIF
             IF ( ALLOCATED( surf_lsm%qcsws ) )  surf%qcsws(mm) = surf_lsm%qcsws(m)
             IF ( ALLOCATED( surf_lsm%qrsws ) )  surf%qrsws(mm) = surf_lsm%qrsws(m)
             IF ( ALLOCATED( surf_lsm%qisws ) )  surf%qisws(mm) = surf_lsm%qisws(m)
             IF ( ALLOCATED( surf_lsm%ncsws ) )  surf%ncsws(mm) = surf_lsm%ncsws(m)
             IF ( ALLOCATED( surf_lsm%nisws ) )  surf%nisws(mm) = surf_lsm%nisws(m)
             IF ( ALLOCATED( surf_lsm%nrsws ) )  surf%nrsws(mm) = surf_lsm%nrsws(m)
             IF ( ALLOCATED( surf_lsm%sasws ) )  surf%sasws(mm) = surf_lsm%sasws(m)

             IF ( ALLOCATED( surf_lsm%usws )        )  surf%usws(mm) = surf_lsm%usws(m)
             IF ( ALLOCATED( surf_lsm%vsws )        )  surf%vsws(mm) = surf_lsm%vsws(m)
             IF ( ALLOCATED( surf_lsm%usvs )        )  surf%usvs(mm) = surf_lsm%usvs(m)
             IF ( ALLOCATED( surf_lsm%vsus )        )  surf%vsus(mm) = surf_lsm%vsus(m)
             IF ( ALLOCATED( surf_lsm%wsus_wsvs)    )  surf%wsus_wsvs(mm) = surf_lsm%wsus_wsvs(m)
             IF ( ALLOCATED( surf_lsm%mom_flux_tke) )                                              &
                surf%mom_flux_tke(0:1,mm) = surf_lsm%mom_flux_tke(0:1,m)

             mm = mm + 1
          ENDDO

          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             IF ( ALLOCATED( surf_usm%us )  )  surf%us(mm)  = surf_usm%us(m)
             IF ( ALLOCATED( surf_usm%ts )  )  surf%ts(mm)  = surf_usm%ts(m)
             IF ( ALLOCATED( surf_usm%qs )  )  surf%qs(mm)  = surf_usm%qs(m)
             IF ( ALLOCATED( surf_usm%ss )  )  surf%ss(mm)  = surf_usm%ss(m)
             IF ( ALLOCATED( surf_usm%qcs ) )  surf%qcs(mm) = surf_usm%qcs(m)
             IF ( ALLOCATED( surf_usm%ncs ) )  surf%ncs(mm) = surf_usm%ncs(m)
             IF ( ALLOCATED( surf_usm%qis ) )  surf%qis(mm) = surf_usm%qis(m)
             IF ( ALLOCATED( surf_usm%nis ) )  surf%nis(mm) = surf_usm%nis(m)
             IF ( ALLOCATED( surf_usm%qrs ) )  surf%qrs(mm) = surf_usm%qrs(m)
             IF ( ALLOCATED( surf_usm%nrs ) )  surf%nrs(mm) = surf_usm%nrs(m)
             IF ( ALLOCATED( surf_usm%ol )  )  surf%ol(mm)  = surf_usm%ol(m)
             IF ( ALLOCATED( surf_usm%rib ) )  surf%rib(mm) = surf_usm%rib(m)
             IF ( ALLOCATED( surf_usm%pt_surface )  )  surf%pt_surface(mm)  = surf_usm%pt_surface(m)
             IF ( ALLOCATED( surf_usm%q_surface )   )  surf%q_surface(mm)   = surf_usm%q_surface(m)
             IF ( ALLOCATED( surf_usm%vpt_surface ) )  surf%vpt_surface(mm) = surf_usm%vpt_surface(m)
             IF ( ALLOCATED( surf_usm%shf )  )  surf%shf(mm)  = surf_usm%shf(m)
             IF ( ALLOCATED( surf_usm%qsws ) )  surf%qsws(mm) = surf_usm%qsws(m)
             IF ( ALLOCATED( surf_usm%ssws ) )  surf%ssws(mm) = surf_usm%ssws(m)
             IF ( ALLOCATED( surf_usm%css )  )  THEN
                DO  lsp = 1,nvar
                   surf%css(lsp,mm) = surf_usm%css(lsp,m)
                ENDDO
             ENDIF
             IF ( ALLOCATED( surf_usm%cssws ) )  THEN
                DO  lsp = 1,nvar
                   surf%cssws(lsp,mm) = surf_usm%cssws(lsp,m)
                ENDDO
             ENDIF
             IF ( ALLOCATED( surf_usm%qcsws ) )  surf%qcsws(mm) = surf_usm%qcsws(m)
             IF ( ALLOCATED( surf_usm%qrsws ) )  surf%qrsws(mm) = surf_usm%qrsws(m)
             IF ( ALLOCATED( surf_usm%qisws ) )  surf%qisws(mm) = surf_usm%qisws(m)
             IF ( ALLOCATED( surf_usm%ncsws ) )  surf%ncsws(mm) = surf_usm%ncsws(m)
             IF ( ALLOCATED( surf_usm%nisws ) )  surf%nisws(mm) = surf_usm%nisws(m)
             IF ( ALLOCATED( surf_usm%nrsws ) )  surf%nrsws(mm) = surf_usm%nrsws(m)
             IF ( ALLOCATED( surf_usm%sasws ) )  surf%sasws(mm) = surf_usm%sasws(m)

             IF ( ALLOCATED( surf_usm%usws )        )  surf%usws(mm) = surf_usm%usws(m)
             IF ( ALLOCATED( surf_usm%vsws )        )  surf%vsws(mm) = surf_usm%vsws(m)
             IF ( ALLOCATED( surf_usm%usvs )        )  surf%usvs(mm) = surf_usm%usvs(m)
             IF ( ALLOCATED( surf_usm%vsus )        )  surf%vsus(mm) = surf_usm%vsus(m)
             IF ( ALLOCATED( surf_usm%wsus_wsvs)    )  surf%wsus_wsvs(mm) = surf_usm%wsus_wsvs(m)
             IF ( ALLOCATED( surf_usm%mom_flux_tke) )                                              &
                surf%mom_flux_tke(0:1,mm) = surf_usm%mom_flux_tke(0:1,m)

             mm = mm + 1
          ENDDO
       ENDDO
    ENDDO
!
!-- Recalculate start- and end indices for gathered surface type.
    start_index = 1
    DO  i = nxl, nxr
       DO  j = nys, nyn
          surf%start_index(j,i) = start_index
          surf%end_index(j,i)   = surf%start_index(j,i) - 1

          DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
             surf%end_index(j,i) = surf%end_index(j,i) + 1
          ENDDO
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             surf%end_index(j,i) = surf%end_index(j,i) + 1
          ENDDO
          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             surf%end_index(j,i) = surf%end_index(j,i) + 1
          ENDDO

          start_index = surf%end_index(j,i) + 1
       ENDDO
    ENDDO
!
!-- Now start writing restart data to file
    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN
!
!--    Output strings for the total number of surfaces on subdomain.
       CALL wrd_write_string( 'ns_on_file' )
       WRITE ( 14 ) ns_on_file(1)

       CALL wrd_write_string( 'surf%start_index' )
       WRITE ( 14 ) surf%start_index

       CALL wrd_write_string( 'surf%end_index' )
       WRITE ( 14 ) surf%end_index

       IF ( ALLOCATED ( surf%us ) )  THEN
          CALL wrd_write_string( 'surf%us' )
          WRITE ( 14 ) surf%us
       ENDIF

       IF ( ALLOCATED ( surf%ts ) )  THEN
          CALL wrd_write_string( 'surf%ts' )
          WRITE ( 14 ) surf%ts
       ENDIF

       IF ( ALLOCATED ( surf%qs ) )  THEN
          CALL wrd_write_string( 'surf%qs' )
          WRITE ( 14 ) surf%qs
       ENDIF

       IF ( ALLOCATED ( surf%ss ) )  THEN
          CALL wrd_write_string( 'surf%ss' )
          WRITE ( 14 ) surf%ss
       ENDIF

       IF ( ALLOCATED ( surf%qcs ) )  THEN
          CALL wrd_write_string( 'surf%qcs' )
          WRITE ( 14 ) surf%qcs
       ENDIF

       IF ( ALLOCATED ( surf%ncs ) )  THEN
          CALL wrd_write_string( 'surf%ncs' )
          WRITE ( 14 ) surf%ncs
       ENDIF

       IF ( ALLOCATED ( surf%qis ) )  THEN
          CALL wrd_write_string( 'surf%qis' )
          WRITE ( 14 ) surf%qis
       ENDIF

       IF ( ALLOCATED ( surf%nis ) )  THEN
          CALL wrd_write_string( 'surf%nis' )
          WRITE ( 14 ) surf%nis
       ENDIF

       IF ( ALLOCATED ( surf%qrs ) )  THEN
          CALL wrd_write_string( 'surf%qrs' )
          WRITE ( 14 ) surf%qrs
       ENDIF

       IF ( ALLOCATED ( surf%nrs ) )  THEN
          CALL wrd_write_string( 'surf%nrs' )
          WRITE ( 14 ) surf%nrs
       ENDIF

       IF ( ALLOCATED ( surf%ol ) )  THEN
          CALL wrd_write_string( 'surf%ol' )
          WRITE ( 14 ) surf%ol
       ENDIF

       IF ( ALLOCATED ( surf%rib ) )  THEN
          CALL wrd_write_string( 'surf%rib' )
          WRITE ( 14 ) surf%rib
       ENDIF

       IF ( ALLOCATED ( surf%pt_surface ) )  THEN
          CALL wrd_write_string( 'surf%pt_surface' )
          WRITE ( 14 ) surf%pt_surface
       ENDIF

       IF ( ALLOCATED ( surf%q_surface ) )  THEN
          CALL wrd_write_string( 'surf%q_surface' )
          WRITE ( 14 ) surf%q_surface
       ENDIF

       IF ( ALLOCATED ( surf%vpt_surface ) )  THEN
          CALL wrd_write_string( 'surf%vpt_surface' )
          WRITE ( 14 ) surf%vpt_surface
       ENDIF

       IF ( ALLOCATED ( surf%shf ) )  THEN
          CALL wrd_write_string( 'surf%shf' )
          WRITE ( 14 ) surf%shf
       ENDIF

       IF ( ALLOCATED ( surf%qsws ) )  THEN
          CALL wrd_write_string( 'surf%qsws' )
          WRITE ( 14 ) surf%qsws
       ENDIF

       IF ( ALLOCATED ( surf%ssws ) )  THEN
          CALL wrd_write_string( 'surf%ssws' )
          WRITE ( 14 ) surf%ssws
       ENDIF

       IF ( ALLOCATED ( surf%css ) )  THEN
          CALL wrd_write_string( 'surf%css' )
          WRITE ( 14 ) surf%css
       ENDIF

       IF ( ALLOCATED ( surf%cssws ) )  THEN
          CALL wrd_write_string( 'surf%cssws' )
          WRITE ( 14 )  surf%cssws
       ENDIF

       IF ( ALLOCATED ( surf%qcsws ) )  THEN
          CALL wrd_write_string( 'surf%qcsws' )
          WRITE ( 14 ) surf%qcsws
       ENDIF

       IF ( ALLOCATED ( surf%ncsws ) )  THEN
          CALL wrd_write_string( 'surf%ncsws' )
          WRITE ( 14 ) surf%ncsws
       ENDIF

       IF ( ALLOCATED ( surf%qisws ) )  THEN
          CALL wrd_write_string( 'surf%qisws' )
          WRITE ( 14 ) surf%qisws
       ENDIF

       IF ( ALLOCATED ( surf%nisws ) )  THEN
          CALL wrd_write_string( 'surf%nisws' )
          WRITE ( 14 ) surf%nisws
       ENDIF

       IF ( ALLOCATED ( surf%qrsws ) )  THEN
          CALL wrd_write_string( 'surf%qrsws' )
          WRITE ( 14 ) surf%qrsws
       ENDIF

       IF ( ALLOCATED ( surf%nrsws ) )  THEN
          CALL wrd_write_string( 'surf%nrsws' )
          WRITE ( 14 ) surf%nrsws
       ENDIF

       IF ( ALLOCATED ( surf%sasws ) )  THEN
          CALL wrd_write_string( 'surf%sasws' )
          WRITE ( 14 ) surf%sasws
       ENDIF

       IF ( ALLOCATED ( surf%usws ) )  THEN
          CALL wrd_write_string( 'surf%usws' )
          WRITE ( 14 ) surf%usws
       ENDIF

       IF ( ALLOCATED ( surf%vsws ) )  THEN
          CALL wrd_write_string( 'surf%vsws' )
          WRITE ( 14 ) surf%vsws
       ENDIF

       IF ( ALLOCATED ( surf%usvs ) )  THEN
          CALL wrd_write_string( 'surf%usvs' )
          WRITE ( 14 ) surf%usvs
       ENDIF

       IF ( ALLOCATED ( surf%vsus ) )  THEN
          CALL wrd_write_string( 'surf%vsus' )
          WRITE ( 14 ) surf%vsus
       ENDIF

       IF ( ALLOCATED ( surf%wsus_wsvs ) )  THEN
          CALL wrd_write_string( 'surf%wsus_wsvs' )
          WRITE ( 14 ) surf%wsus_wsvs
       ENDIF

       IF ( ALLOCATED ( surf%mom_flux_tke ) )  THEN
          CALL wrd_write_string( 'surf%mom_flux_tke' )
          WRITE ( 14 ) surf%mom_flux_tke
       ENDIF
!
!--    Now treat model-top fluxes.
       CALL wrd_write_string( 'ns_on_file_top' )
       WRITE ( 14 ) surf_top%ns

       CALL wrd_write_string( 'surf_top%start_index' )
       WRITE ( 14 ) surf_top%start_index

       CALL wrd_write_string( 'surf_top%end_index' )
       WRITE ( 14 ) surf_top%end_index

       IF ( ALLOCATED ( surf_top%usws ) )  THEN
          CALL wrd_write_string( 'surf_top%usws' )
          WRITE ( 14 ) surf_top%usws
       ENDIF
       IF ( ALLOCATED ( surf_top%vsws ) )  THEN
          CALL wrd_write_string( 'surf_top%vsws' )
          WRITE ( 14 ) surf_top%vsws
       ENDIF
       IF ( ALLOCATED ( surf_top%shf ) )  THEN
          CALL wrd_write_string( 'surf_top%shf' )
          WRITE ( 14 ) surf_top%shf
       ENDIF
       IF ( ALLOCATED ( surf_top%qsws ) )  THEN
          CALL wrd_write_string( 'surf_top%qsws' )
          WRITE ( 14 ) surf_top%qsws
       ENDIF
       IF ( ALLOCATED ( surf_top%ssws ) )  THEN
          CALL wrd_write_string( 'surf_top%ssws' )
          WRITE ( 14 ) surf_top%ssws
       ENDIF
       IF ( ALLOCATED ( surf_top%cssws ) )  THEN
          CALL wrd_write_string( 'surf_top%cssws' )
          WRITE ( 14 ) surf_top%cssws
       ENDIF
       IF ( ALLOCATED ( surf_top%qcsws ) )  THEN
          CALL wrd_write_string( 'surf_top%qcsws' )
          WRITE ( 14 ) surf_top%qcsws
       ENDIF
       IF ( ALLOCATED ( surf_top%ncsws ) )  THEN
          CALL wrd_write_string( 'surf_top%ncsws' )
          WRITE ( 14 ) surf_top%ncsws
       ENDIF
       IF ( ALLOCATED ( surf_top%qrsws ) )  THEN
          CALL wrd_write_string( 'surf_top%qrsws' )
          WRITE ( 14 ) surf_top%qrsws
       ENDIF
       IF ( ALLOCATED ( surf_top%nrsws ) )  THEN
          CALL wrd_write_string( 'surf_top%nrsws' )
          WRITE ( 14 ) surf_top%nrsws
       ENDIF
       IF ( ALLOCATED ( surf_top%qisws ) )  THEN
          CALL wrd_write_string( 'surf_top%qisws' )
          WRITE ( 14 ) surf_top%qisws
       ENDIF
       IF ( ALLOCATED ( surf_top%nisws ) )  THEN
          CALL wrd_write_string( 'surf_top%nisws' )
          WRITE ( 14 ) surf_top%nisws
       ENDIF
       IF ( ALLOCATED ( surf_top%sasws ) )  THEN
          CALL wrd_write_string( 'surf_top%sasws' )
          WRITE ( 14 ) surf_top%sasws
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!
!--    Write-out all surfaces (different types and facings).
!--    All data writen with rd_mpi_io_write_surface are globally indexed 1d-arrays.
       ns_on_file(1) = 0
       CALL rd_mpi_io_surface_filetypes( surf%start_index, surf%end_index, surface_data_to_write,  &
                                         global_start_index, global_end_index )
       ns_on_file(1) = total_number_of_surface_elements

       CALL wrd_mpi_io( 'global_start_index', global_start_index )
       CALL wrd_mpi_io( 'global_end_index', global_end_index )

       IF ( ALLOCATED ( surf%us ) )  THEN
          CALL wrd_mpi_io_surface ( 'surf%us', surf%us )
       ENDIF

       IF ( ALLOCATED ( surf%ts ) )  THEN
          CALL wrd_mpi_io_surface ( 'surf%ts', surf%ts )
       ENDIF

       IF ( ALLOCATED ( surf%qs ) )  THEN
          CALL wrd_mpi_io_surface ( 'surf%qs', surf%qs )
       ENDIF

       IF ( ALLOCATED ( surf%ss ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%ss', surf%ss )
       ENDIF

       IF ( ALLOCATED ( surf%qcs ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qcs', surf%qcs )
       ENDIF

       IF ( ALLOCATED ( surf%ncs ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%ncs', surf%ncs )
       ENDIF

       IF ( ALLOCATED ( surf%qis ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qis', surf%qis )
       ENDIF

       IF ( ALLOCATED ( surf%nis ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%nis', surf%nis )
       ENDIF

       IF ( ALLOCATED ( surf%qrs ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qrs', surf%qrs )
       ENDIF

       IF ( ALLOCATED ( surf%nrs ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%nrs', surf%nrs )
       ENDIF

       IF ( ALLOCATED ( surf%ol ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%ol', surf%ol )
       ENDIF

       IF ( ALLOCATED ( surf%rib ) )  THEN
         CALL wrd_mpi_io_surface( 'surf%rib', surf%rib )
       ENDIF

       IF ( ALLOCATED ( surf%pt_surface ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%pt_surface', surf%pt_surface )
       ENDIF

       IF ( ALLOCATED ( surf%q_surface ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%q_surface', surf%q_surface )
       ENDIF

       IF ( ALLOCATED ( surf%vpt_surface ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%vpt_surface', surf%vpt_surface )
       ENDIF

       IF ( ALLOCATED ( surf%shf ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%shf', surf%shf )
       ENDIF

       IF ( ALLOCATED ( surf%qsws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qsws', surf%qsws )
       ENDIF

       IF ( ALLOCATED ( surf%ssws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%ssws', surf%ssws )
       ENDIF

       IF ( ALLOCATED ( surf%css ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%css', surf%css )
       ENDIF

       IF ( ALLOCATED ( surf%cssws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%cssws', surf%cssws )
       ENDIF

       IF ( ALLOCATED ( surf%qcsws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qcsws', surf%qcsws )
       ENDIF

       IF ( ALLOCATED ( surf%ncsws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%ncsws', surf%ncsws )
       ENDIF

       IF ( ALLOCATED ( surf%qisws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qisws', surf%qisws )
       ENDIF

       IF ( ALLOCATED ( surf%nisws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%nisws', surf%nisws )
       ENDIF

       IF ( ALLOCATED ( surf%qrsws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%qrsws', surf%qrsws )
       ENDIF

       IF ( ALLOCATED ( surf%nrsws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%nrsws', surf%nrsws )
       ENDIF

       IF ( ALLOCATED ( surf%sasws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%sasws', surf%sasws )
       ENDIF

       IF ( ALLOCATED ( surf%usws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%usws', surf%usws )
       ENDIF

       IF ( ALLOCATED ( surf%vsws ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%vsws', surf%vsws )
       ENDIF

       IF ( ALLOCATED ( surf%usvs ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%usvs', surf%usvs )
       ENDIF

       IF ( ALLOCATED ( surf%vsus ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%vsus', surf%vsus )
       ENDIF

       IF ( ALLOCATED ( surf%wsus_wsvs ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%wsus_wsvs', surf%wsus_wsvs )
       ENDIF

       IF ( ALLOCATED ( surf%mom_flux_tke ) )  THEN
          CALL wrd_mpi_io_surface( 'surf%mom_flux_tke', surf%mom_flux_tke )
       ENDIF

       CALL wrd_mpi_io_global_array( 'ns_on_file', ns_on_file )

!
!--    Now treat model-top fluxes. Only proceed when top surfaces are defined, else
!--    rrd_mpi_io_global_array may lead to a segmentation fault.
       IF ( use_top_fluxes )  THEN

          CALL rd_mpi_io_surface_filetypes( surf_top%start_index, surf_top%end_index,              &
                                            surface_data_to_write, global_start_index,             &
                                            global_end_index )

          ns_on_file(1) = total_number_of_surface_elements

          CALL wrd_mpi_io( 'global_start_index_top', global_start_index )
          CALL wrd_mpi_io( 'global_end_index_top', global_end_index )
          CALL wrd_mpi_io_global_array( 'ns_on_file_top', ns_on_file )
!
!--       Check if data is available. In contrast to "normal" surfaces, which are always available,
!--       model-top fluxes are not necessarily present.
          IF ( surface_data_to_write )  THEN
             IF ( ALLOCATED ( surf_top%usws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%usws', surf_top%usws )
             ENDIF
             IF ( ALLOCATED ( surf_top%vsws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%vsws', surf_top%vsws )
             ENDIF
             IF ( ALLOCATED ( surf_top%shf ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%shf', surf_top%shf )
             ENDIF
             IF ( ALLOCATED ( surf_top%qsws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%qsws', surf_top%qsws )
             ENDIF
             IF ( ALLOCATED ( surf_top%ssws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%ssws', surf_top%ssws )
             ENDIF
             IF ( ALLOCATED ( surf_top%cssws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%cssws', surf_top%cssws )
             ENDIF
             IF ( ALLOCATED ( surf_top%qcsws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%qcsws', surf_top%qcsws )
             ENDIF
             IF ( ALLOCATED ( surf_top%ncsws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%ncsws', surf_top%ncsws )
             ENDIF
             IF ( ALLOCATED ( surf_top%qrsws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%qrsws', surf_top%qrsws )
             ENDIF
             IF ( ALLOCATED ( surf_top%nrsws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%nrsws', surf_top%nrsws )
             ENDIF
             IF ( ALLOCATED ( surf_top%qisws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%qisws', surf_top%qisws )
             ENDIF
             IF ( ALLOCATED ( surf_top%nisws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%nisws', surf_top%nisws )
             ENDIF
             IF ( ALLOCATED ( surf_top%sasws ) )  THEN
                CALL wrd_mpi_io_surface ( 'surf_top%sasws', surf_top%sasws )
             ENDIF
          ENDIF

       ENDIF

    ENDIF

 END SUBROUTINE surface_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads surface-related restart data in Fortran binary format. Please note, restart data for a
!> certain surface orientation (e.g. horizontal upward-facing) is stored in one array, even if
!> surface elements may belong to different surface types natural or urban for example). Surface
!> elements are redistributed into its respective surface types within this routine. This allows
!> e.g. changing the surface type after reading the restart data, which might be required in case
!> of cyclic_fill mode.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_rrd_local_ftn( kk, nxlf, nxlc, nxl_on_file, nxrf, nxr_on_file, nynf,           &
                                   nyn_on_file, nysf, nysc, nys_on_file, found )


    IMPLICIT NONE

    INTEGER(iwp) ::  i            !< running index along x-direction, refers to former domain size
    INTEGER(iwp) ::  ic           !< running index along x-direction, refers to current domain size
    INTEGER(iwp) ::  j            !< running index along y-direction, refers to former domain size
    INTEGER(iwp) ::  jc           !< running index along y-direction, refers to former domain size
    INTEGER(iwp) ::  m            !< running index for surface elements, refers to gathered array encompassing all surface types
    INTEGER(iwp) ::  mm           !< running index for surface elements, refers to individual surface types
    INTEGER(iwp) ::  kk           !< running index over previous input files covering current local domain
    INTEGER(iwp) ::  nxlc         !< index of left boundary on current subdomain
    INTEGER(iwp) ::  nxlf         !< index of left boundary on former subdomain
    INTEGER(iwp) ::  nxl_on_file  !< index of left boundary on former local domain
    INTEGER(iwp) ::  nxrf         !< index of right boundary on former subdomain
    INTEGER(iwp) ::  nxr_on_file  !< index of right boundary on former local domain
    INTEGER(iwp) ::  nynf         !< index of north boundary on former subdomain
    INTEGER(iwp) ::  nyn_on_file  !< index of norht boundary on former local domain
    INTEGER(iwp) ::  nysc         !< index of south boundary on current subdomain
    INTEGER(iwp) ::  nysf         !< index of south boundary on former subdomain
    INTEGER(iwp) ::  nys_on_file  !< index of south boundary on former local domain

    LOGICAL ::  surf_match_def  !< flag indicating that surface element is of default type
    LOGICAL ::  surf_match_lsm  !< flag indicating that surface element is of natural type
    LOGICAL ::  surf_match_usm  !< flag indicating that surface element is of urban type
    LOGICAL ::  surfaces_top    !< flag indicating that surface element belongs to surf_top

    LOGICAL, INTENT(OUT) ::  found  !<

    TYPE(surf_type), SAVE ::  surf   !< surface type to handle data on file
    TYPE(surf_type), SAVE ::  surf_t !< surface type to handle model-top data on file


    found = .TRUE.
!
!-- Array "mapped" is allocated further below, after the number of surface elements "ns_on_file"
!-- is known.
    IF ( ALLOCATED( mapped ) )  mapped = .FALSE.

    SELECT CASE ( restart_string(1:length) )
!
!--    Read the number of surface elements and allocate arrays
       CASE ( 'ns_on_file' )
          IF ( kk == 1 )  THEN
             READ ( 13 )  ns_on_file(1)
!
!--          In case of changing mpi topology, this routine could be called more than once.
!--          Hence, arrays need to be deallocated before allocated again.
             IF ( ALLOCATED( surf%start_index ) )  CALL deallocate_surface_attributes( surf )
!
!--          Allocate memory for number of surface elements on file.
!--          Please note, this number is not necessarily the same as the final number of surface
!--          elements on local domain, which is the case if processor topology changes during
!--          restart runs.
             surf%ns = ns_on_file(1)
             CALL allocate_surface_attributes( surf, nys_on_file, nyn_on_file, nxl_on_file,        &
                                               nxr_on_file )
!
!--          Allocate array to indicate already mapped surfaces.
             IF ( .NOT. ALLOCATED( mapped ) )  THEN
                ALLOCATE( mapped(1:surf%ns) )
                mapped = .FALSE.
             ENDIF
          ENDIF
!
!--    Read start and end indices of surface elements at each (ji)-gridpoint
       CASE ( 'surf%start_index' )
          IF ( kk == 1 )  READ ( 13 ) surf%start_index
       CASE ( 'surf%end_index' )
          IF ( kk == 1 )  READ ( 13 ) surf%end_index
!
!--    Read specific attributes
       CASE ( 'surf%us' )
          IF ( ALLOCATED( surf%us )  .AND.  kk == 1 )  READ ( 13 ) surf%us
       CASE ( 'surf%ts' )
          IF ( ALLOCATED( surf%ts )  .AND.  kk == 1 )  READ ( 13 ) surf%ts
       CASE ( 'surf%qs' )
          IF ( ALLOCATED( surf%qs )  .AND.  kk == 1 )  READ ( 13 ) surf%qs
       CASE ( 'surf%ss' )
          IF ( ALLOCATED( surf%ss )  .AND.  kk == 1 )  READ ( 13 ) surf%ss
       CASE ( 'surf%qcs' )
          IF ( ALLOCATED( surf%qcs )  .AND.  kk == 1 )  READ ( 13 ) surf%qcs
       CASE ( 'surf%ncs' )
          IF ( ALLOCATED( surf%ncs )  .AND.  kk == 1 )  READ ( 13 ) surf%ncs
       CASE ( 'surf%qis' )
          IF ( ALLOCATED( surf%qis )  .AND.  kk == 1 )  READ ( 13 ) surf%qis
       CASE ( 'surf%nis' )
          IF ( ALLOCATED( surf%nis )  .AND.  kk == 1 )  READ ( 13 ) surf%nis
       CASE ( 'surf%qrs' )
          IF ( ALLOCATED( surf%qrs )  .AND.  kk == 1 )  READ ( 13 ) surf%qrs
       CASE ( 'surf%nrs' )
          IF ( ALLOCATED( surf%nrs )  .AND.  kk == 1 )  READ ( 13 ) surf%nrs
       CASE ( 'surf%ol' )
          IF ( ALLOCATED( surf%ol )  .AND.  kk == 1 )  READ ( 13 ) surf%ol
       CASE ( 'surf%rib' )
          IF ( ALLOCATED( surf%rib )  .AND.  kk == 1 )  READ ( 13 ) surf%rib
       CASE ( 'surf%pt_surface' )
          IF ( ALLOCATED( surf%pt_surface )  .AND.  kk == 1 )  READ ( 13 ) surf%pt_surface
       CASE ( 'surf%q_surface' )
          IF ( ALLOCATED( surf%q_surface )  .AND.  kk == 1 )  READ ( 13 ) surf%q_surface
       CASE ( 'surf%vpt_surface' )
          IF ( ALLOCATED( surf%vpt_surface )  .AND.  kk == 1 )  READ ( 13 ) surf%vpt_surface
       CASE ( 'surf%shf' )
          IF ( ALLOCATED( surf%shf )  .AND.  kk == 1 )  READ ( 13 ) surf%shf
       CASE ( 'surf%qsws' )
          IF ( ALLOCATED( surf%qsws )  .AND.  kk == 1 )  READ ( 13 ) surf%qsws
       CASE ( 'surf%ssws' )
          IF ( ALLOCATED( surf%ssws )  .AND.  kk == 1 )  READ ( 13 ) surf%ssws
       CASE ( 'surf%css' )
          IF ( ALLOCATED( surf%css )  .AND.  kk == 1 )  READ ( 13 ) surf%css
       CASE ( 'surf%cssws' )
          IF ( ALLOCATED( surf%cssws )  .AND.  kk == 1 )  READ ( 13 ) surf%cssws
       CASE ( 'surf%qcsws' )
          IF ( ALLOCATED( surf%qcsws )  .AND.  kk == 1 )  READ ( 13 ) surf%qcsws
       CASE ( 'surf%ncsws' )
          IF ( ALLOCATED( surf%ncsws )  .AND.  kk == 1 )  READ ( 13 ) surf%ncsws
       CASE ( 'surf%qisws' )
          IF ( ALLOCATED( surf%qisws )  .AND.  kk == 1 )  READ ( 13 ) surf%qisws
       CASE ( 'surf%nisws' )
          IF ( ALLOCATED( surf%nisws )  .AND.  kk == 1 )  READ ( 13 ) surf%nisws
       CASE ( 'surf%qrsws' )
          IF ( ALLOCATED( surf%qrsws )  .AND.  kk == 1 )  READ ( 13 ) surf%qrsws
       CASE ( 'surf%nrsws' )
          IF ( ALLOCATED( surf%nrsws )  .AND.  kk == 1 )  READ ( 13 ) surf%nrsws
       CASE ( 'surf%sasws' )
          IF ( ALLOCATED( surf%sasws )  .AND.  kk == 1 )  READ ( 13 ) surf%sasws
       CASE ( 'surf%usws' )
          IF ( ALLOCATED( surf%usws )  .AND.  kk == 1 )  READ ( 13 ) surf%usws
       CASE ( 'surf%vsws' )
          IF ( ALLOCATED( surf%vsws )  .AND.  kk == 1 )  READ ( 13 ) surf%vsws
       CASE ( 'surf%usvs' )
          IF ( ALLOCATED( surf%usvs )  .AND.  kk == 1 )  READ ( 13 ) surf%usvs
       CASE ( 'surf%vsus' )
          IF ( ALLOCATED( surf%vsus )  .AND.  kk == 1 )  READ ( 13 ) surf%vsus
       CASE ( 'surf%wsus_wsvs' )
          IF ( ALLOCATED( surf%wsus_wsvs )  .AND.  kk == 1 )  READ ( 13 ) surf%wsus_wsvs
       CASE ( 'surf%mom_flux_tke' )
          IF ( ALLOCATED( surf%mom_flux_tke )  .AND.  kk == 1 )  READ ( 13 ) surf%mom_flux_tke
!
!--    Top boundary fluxes
!--    Read the number of surface elements and allocate arrays
       CASE ( 'ns_on_file_top' )
          IF ( kk == 1 )  THEN
             READ ( 13 )  ns_on_file(1)
!
!--          In case of changing mpi topology, this routine could be called more than once.
!--          Hence, arrays need to be deallocated before allocated again.
             IF ( ALLOCATED( surf_t%start_index ) )  THEN
                CALL deallocate_surface_attributes_top( surf_t )
             ENDIF
!
!--          Allocate memory for number of top surface elements on file.
!--          Please note, this number is not necessarily the same as the final number of surface
!--          elements on local domain, which is the case if processor topology changes during
!--          restart runs.
             surf_t%ns = ns_on_file(1)
             CALL allocate_surface_attributes_top( surf_t, nys_on_file, nyn_on_file, nxl_on_file,  &
                                                   nxr_on_file )

          ENDIF
!
!--    Read start and end indices of surface elements at each (ji)-gridpoint
       CASE ( 'surf_top%start_index' )
          IF ( kk == 1 )  READ ( 13 ) surf_t%start_index
       CASE ( 'surf_top%end_index' )
          IF ( kk == 1 )  READ ( 13 ) surf_t%end_index
       CASE ( 'surf_top%usws' )
          IF ( ALLOCATED( surf_t%usws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%usws
       CASE ( 'surf_top%vsws' )
          IF ( ALLOCATED( surf_t%vsws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%vsws
       CASE ( 'surf_top%shf' )
          IF ( ALLOCATED( surf_t%shf )  .AND.  kk == 1 )  READ ( 13 ) surf_t%shf
       CASE ( 'surf_top%qsws' )
          IF ( ALLOCATED( surf_t%qsws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%qsws
       CASE ( 'surf_top%ssws' )
          IF ( ALLOCATED( surf_t%ssws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%ssws
       CASE ( 'surf_top%cssws' )
          IF ( ALLOCATED( surf_t%cssws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%cssws
       CASE ( 'surf_top%qcsws' )
          IF ( ALLOCATED( surf_t%qcsws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%qcsws
       CASE ( 'surf_top%ncsws' )
          IF ( ALLOCATED( surf_t%ncsws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%ncsws
       CASE ( 'surf_top%qrsws' )
          IF ( ALLOCATED( surf_t%qrsws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%qrsws
       CASE ( 'surf_top%nrsws' )
          IF ( ALLOCATED( surf_t%nrsws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%nrsws
       CASE ( 'surf_top%qisws' )
          IF ( ALLOCATED( surf_t%qisws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%qisws
       CASE ( 'surf_top%nisws' )
          IF ( ALLOCATED( surf_t%nisws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%nisws
       CASE ( 'surf_top%sasws' )
          IF ( ALLOCATED( surf_t%sasws )  .AND.  kk == 1 )  READ ( 13 ) surf_t%sasws

       CASE DEFAULT

             found = .FALSE.

    END SELECT

    IF ( .NOT. found )  RETURN

!
!-- Redistribute surface elements on its respective type.
    IF ( .NOT. INDEX( restart_string(1:length), 'ns_on_file'   ) /= 0  .AND.                       &
         .NOT. INDEX( restart_string(1:length), '%start_index' ) /= 0  .AND.                       &
         .NOT. INDEX( restart_string(1:length), '%end_index'   ) /= 0 )                            &
    THEN

       ic = nxlc
       DO  i = nxlf, nxrf
          jc = nysc
          DO  j = nysf, nynf
!
!--          Determine type of surface element, i.e. default, natural, urban, at current grid point.
             surf_match_def  = surf_def%end_index(jc,ic) >= surf_def%start_index(jc,ic)
             surf_match_lsm  = surf_lsm%end_index(jc,ic) >= surf_lsm%start_index(jc,ic)
             surf_match_usm  = surf_usm%end_index(jc,ic) >= surf_usm%start_index(jc,ic)
             surfaces_top    = INDEX( restart_string(1:length), 'surf_top%' ) /= 0
!
!--          Write restart data onto default-type surfaces if required.
             IF ( surf_match_def  .AND.  .NOT. surfaces_top )  THEN
!
!--             Set the start index for the local surface element
                mm = surf_def%start_index(jc,ic)
!
!--             For index pair (j,i) on file loop from start to end index, and in case the local
!--             surface element mm is smaller than the local end index, assign the respective
!--             surface data to this element. Since all surfaces in the restart input array
!--             are gathered into a one-dimensional array and are distributed over different
!--             surface types, it need to be checked that input data is not mapped twice onto
!--             different surface types. Hence, store information of mapped surfaces in a didicated
!--             array and only write input surface data onto the target surface types that has not
!--             been written before.
                DO  m = surf%start_index(j,i), surf%end_index(j,i)
                   IF ( surf_def%end_index(jc,ic) >= mm  .AND.  .NOT. mapped(m) )  THEN
                      CALL restore_surface_elements( surf_def, mm, surf, m )
                      mapped(m) = .TRUE.
                   ENDIF
                   mm = mm + 1
                ENDDO
             ENDIF
!
!--          Same for natural-type surfaces. Please note, it is implicitly assumed that natural
!--          surface elements are below urban surface elements if there are several horizontal
!--          surfaces at (j,i). An example would be bridges or building walls of top of a cliff.
             IF ( surf_match_lsm  .AND.  .NOT. surfaces_top )  THEN
                mm = surf_lsm%start_index(jc,ic)
                DO  m = surf%start_index(j,i), surf%end_index(j,i)
                   IF ( surf_lsm%end_index(jc,ic) >= mm  .AND.  .NOT. mapped(m) )  THEN
                      CALL restore_surface_elements( surf_lsm, mm, surf, m )
                      mapped(m) = .TRUE.
                      mm = mm + 1
                   ENDIF
                ENDDO
             ENDIF
!
!--          Same for urban-type surfaces
             IF ( surf_match_usm  .AND.  .NOT. surfaces_top )  THEN
                mm = surf_usm%start_index(jc,ic)
                DO  m = surf%start_index(j,i), surf%end_index(j,i)
                   IF ( surf_usm%end_index(jc,ic) >= mm  .AND.  .NOT. mapped(m) )  THEN
                      CALL restore_surface_elements( surf_usm, mm, surf, m )
                      mapped(m) = .TRUE.
                      mm = mm + 1
                   ENDIF
                ENDDO
             ENDIF
!
!--          Treat model top surfaces
             IF ( surfaces_top )  THEN
                mm = surf_top%start_index(jc,ic)
                DO  m = surf_t%start_index(j,i), surf_t%end_index(j,i)
                   IF ( surf_top%end_index(jc,ic) >= mm )                                        &
                      CALL restore_surface_elements( surf_top, mm, surf_t, m )
                   mm = mm + 1
                ENDDO
             ENDIF

             jc = jc + 1
          ENDDO
          ic = ic + 1
       ENDDO
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Restores surface elements back on its respective type.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE restore_surface_elements( surf_target, m_target, surf_file, m_file )

    IMPLICIT NONE

    INTEGER(iwp) ::  m_file    !< respective surface-element index of current surface array
    INTEGER(iwp) ::  m_target  !< respecitve surface-element index of surface array on file
    INTEGER(iwp) ::  lsp       !< running index chemical species

    TYPE(surf_type) ::  surf_target  !< target surface type
    TYPE(surf_type) ::  surf_file    !< surface type on file


    IF ( INDEX( restart_string(1:length), '%us' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%us )  .AND.  ALLOCATED( surf_file%us ) )                        &
          surf_target%us(m_target) = surf_file%us(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%ol' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%ol )  .AND.  ALLOCATED( surf_file%ol ) )                        &
          surf_target%ol(m_target) = surf_file%ol(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%pt_surface' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%pt_surface )  .AND.  ALLOCATED( surf_file%pt_surface ) )        &
          surf_target%pt_surface(m_target) = surf_file%pt_surface(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%q_surface' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%q_surface )  .AND.  ALLOCATED( surf_file%q_surface ) )          &
          surf_target%q_surface(m_target) = surf_file%q_surface(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%vpt_surface' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%vpt_surface )  .AND.  ALLOCATED( surf_file%vpt_surface ) )      &
          surf_target%vpt_surface(m_target) = surf_file%vpt_surface(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%ts' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%ts )  .AND.  ALLOCATED( surf_file%ts ) )                        &
          surf_target%ts(m_target) = surf_file%ts(m_file)
    ENDIF
!
!-- Do not restore constant surface fluxes because this would overwrite runtime stettings.
    IF ( .NOT. constant_heatflux )  THEN
       IF ( INDEX( restart_string(1:length), '%shf' ) /= 0 )  THEN
          IF ( ALLOCATED( surf_target%shf )  .AND.  ALLOCATED( surf_file%shf ) )                   &
             surf_target%shf(m_target) = surf_file%shf(m_file)
       ENDIF
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qs )  .AND.  ALLOCATED( surf_file%qs ) )                        &
          surf_target%qs(m_target) = surf_file%qs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qsws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qsws )  .AND.  ALLOCATED( surf_file%qsws ) )                    &
          surf_target%qsws(m_target) = surf_file%qsws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%ss' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%ss )  .AND.  ALLOCATED( surf_file%ss ) )                        &
          surf_target%ss(m_target) = surf_file%ss(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%ssws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%ssws )  .AND.  ALLOCATED( surf_file%ssws ) )                    &
          surf_target%ssws(m_target) = surf_file%ssws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%css' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%css )  .AND.  ALLOCATED( surf_file%css ) )  THEN
          DO  lsp = 1, nvar
             surf_target%css(lsp,m_target) = surf_file%css(lsp,m_file)
          ENDDO
       ENDIF
    ENDIF
    IF ( INDEX( restart_string(1:length), '%cssws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%cssws )  .AND.  ALLOCATED( surf_file%cssws ) )  THEN
          DO  lsp = 1, nvar
             surf_target%cssws(lsp,m_target) = surf_file%cssws(lsp,m_file)
          ENDDO
       ENDIF
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qcs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qcs )  .AND.  ALLOCATED( surf_file%qcs ) )                      &
         surf_target%qcs(m_target) = surf_file%qcs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qcsws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qcsws )  .AND.  ALLOCATED( surf_file%qcsws ) )                  &
          surf_target%qcsws(m_target) = surf_file%qcsws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%ncs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%ncs )  .AND.  ALLOCATED( surf_file%ncs ) )                      &
          surf_target%ncs(m_target) = surf_file%ncs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%ncsws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%ncsws )  .AND.  ALLOCATED( surf_file%ncsws ) )                  &
          surf_target%ncsws(m_target) = surf_file%ncsws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qis' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qis )  .AND.  ALLOCATED( surf_file%qis ) )                      &
         surf_target%qis(m_target) = surf_file%qis(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qisws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qisws )  .AND.  ALLOCATED( surf_file%qisws ) )                  &
          surf_target%qisws(m_target) = surf_file%qisws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%nis' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%nis )  .AND.  ALLOCATED( surf_file%nis ) )                      &
          surf_target%nis(m_target) = surf_file%nis(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%nisws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%nisws )  .AND.  ALLOCATED( surf_file%nisws ) )                  &
          surf_target%nisws(m_target) = surf_file%nisws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qrs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qrs )  .AND.  ALLOCATED( surf_file%qrs ) )                      &
         surf_target%qrs(m_target) = surf_file%qrs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%qrsws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%qrsws )  .AND.  ALLOCATED( surf_file%qrsws ) )                  &
          surf_target%qrsws(m_target) = surf_file%qrsws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%nrs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%nrs )  .AND.  ALLOCATED( surf_file%nrs ) )                      &
          surf_target%nrs(m_target) = surf_file%nrs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%nrsws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%nrsws )  .AND.  ALLOCATED( surf_file%nrsws ) )                  &
          surf_target%nrsws(m_target) = surf_file%nrsws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%sasws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%sasws )  .AND.  ALLOCATED( surf_file%sasws ) )                  &
          surf_target%sasws(m_target) = surf_file%sasws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%usws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%usws )  .AND.  ALLOCATED( surf_file%usws ) )                    &
          surf_target%usws(m_target) = surf_file%usws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%vsws' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%vsws )  .AND.  ALLOCATED( surf_file%vsws ) )                    &
          surf_target%vsws(m_target) = surf_file%vsws(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%usvs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%usvs )  .AND.  ALLOCATED( surf_file%usvs ) )                    &
          surf_target%usvs(m_target) = surf_file%usvs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%vsus' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%vsus )  .AND.  ALLOCATED( surf_file%vsus ) )                    &
          surf_target%vsus(m_target) = surf_file%vsus(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%wsus_wsvs' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%wsus_wsvs )  .AND.  ALLOCATED( surf_file%wsus_wsvs ) )          &
          surf_target%wsus_wsvs(m_target) = surf_file%wsus_wsvs(m_file)
    ENDIF

    IF ( INDEX( restart_string(1:length), '%mom_flux_tke' ) /= 0 )  THEN
       IF ( ALLOCATED( surf_target%mom_flux_tke )  .AND.                                           &
            ALLOCATED( surf_file%mom_flux_tke ) )                                                  &
          surf_target%mom_flux_tke(0:1,m_target) = surf_file%mom_flux_tke(0:1,m_file)
    ENDIF


 END SUBROUTINE restore_surface_elements

 END SUBROUTINE surface_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads surface-related restart data in MPI-IO format. TO_DO: this routine needs to be adjusted for
!> cyclic_fill mode
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_rrd_local_mpi


    IMPLICIT NONE


    INTEGER(iwp) ::  i  !< loop index, x-direction
    INTEGER(iwp) ::  j  !< loop index, y-direction
    INTEGER(iwp) ::  m  !< loop index for surface types - target array
    INTEGER(iwp) ::  mm !< loop index for surface types - file array

    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end_index    !< end index for surface data (MPI-IO)
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start_index  !< start index for surface data (MPI-IO)

    LOGICAL ::  data_to_read    !< cycle in l loop, if no values to read
    LOGICAL ::  surf_match_def  !< flag indicating that surface element is of default type
    LOGICAL ::  surf_match_lsm  !< flag indicating that surface element is of natural type
    LOGICAL ::  surf_match_usm  !< flag indicating that surface element is of urban type

    TYPE(surf_type) ::  surf   !< gathered surfaces, contains all surface types except for the model-top surfaces
    TYPE(surf_type) ::  surf_t !< gathered model-top surfaces

!
!-- Get total number of surface points on the file
    CALL rrd_mpi_io_global_array( 'ns_on_file', ns_on_file )

    IF ( ALLOCATED( surf%start_index ) )  CALL deallocate_surface_attributes( surf )

    ALLOCATE( surf%start_index(nys:nyn,nxl:nxr) )
    ALLOCATE( surf%end_index(nys:nyn,nxl:nxr) )
    surf%start_index = 0
    surf%end_index   = -1

    CALL rrd_mpi_io( 'global_start_index', global_start_index )
    CALL rrd_mpi_io( 'global_end_index' , global_end_index )

    CALL rd_mpi_io_surface_filetypes( surf%start_index, surf%end_index, data_to_read,              &
                                      global_start_index, global_end_index )

    surf%ns = MAX( 2, MAXVAL( surf%end_index ) )

    CALL allocate_surface_attributes( surf, nys, nyn, nxl, nxr, no_allocate_index_arrays = .TRUE. )
!
!-- Allocate array to indicate already mapped surfaces.
    ALLOCATE( mapped(1:surf%ns) )
    mapped = .FALSE.

    IF ( ALLOCATED ( surf%us )  )  CALL rrd_mpi_io_surface( 'surf%us', surf%us )
    IF ( ALLOCATED ( surf%ts )  )  CALL rrd_mpi_io_surface( 'surf%ts', surf%ts )
    IF ( ALLOCATED ( surf%qs )  )  CALL rrd_mpi_io_surface( 'surf%qs', surf%qs )
    IF ( ALLOCATED ( surf%ss )  )  CALL rrd_mpi_io_surface( 'surf%ss', surf%ss )
    IF ( ALLOCATED ( surf%qcs ) )  CALL rrd_mpi_io_surface( 'surf%qcs', surf%qcs )
    IF ( ALLOCATED ( surf%ncs ) )  CALL rrd_mpi_io_surface( 'surf%ncs', surf%ncs )
    IF ( ALLOCATED ( surf%qis ) )  CALL rrd_mpi_io_surface( 'surf%qis', surf%qis )
    IF ( ALLOCATED ( surf%nis ) )  CALL rrd_mpi_io_surface( 'surf%nis', surf%nis )
    IF ( ALLOCATED ( surf%qrs ) )  CALL rrd_mpi_io_surface( 'surf%qrs', surf%qrs )
    IF ( ALLOCATED ( surf%nrs ) )  CALL rrd_mpi_io_surface( 'surf%nrs', surf%nrs )
    IF ( ALLOCATED ( surf%ol )  )  CALL rrd_mpi_io_surface( 'surf%ol',  surf%ol )
    IF ( ALLOCATED ( surf%rib ) )  CALL rrd_mpi_io_surface( 'surf%rib',  surf%rib )
    IF ( ALLOCATED ( surf%pt_surface  ) )                                                          &
       CALL rrd_mpi_io_surface( 'surf%pt_surface', surf%pt_surface )
    IF ( ALLOCATED ( surf%q_surface   ) )                                                          &
       CALL rrd_mpi_io_surface( 'surf%q_surface', surf%q_surface )
    IF ( ALLOCATED ( surf%vpt_surface ) )                                                          &
       CALL rrd_mpi_io_surface( 'surf%vpt_surface', surf%vpt_surface )
    IF ( ALLOCATED ( surf%shf )   )  CALL rrd_mpi_io_surface( 'surf%shf', surf%shf )
    IF ( ALLOCATED ( surf%qsws )  )  CALL rrd_mpi_io_surface( 'surf%qsws', surf%qsws )
    IF ( ALLOCATED ( surf%ssws )  )  CALL rrd_mpi_io_surface( 'surf%ssws', surf%ssws )
    IF ( ALLOCATED ( surf%css )   )  CALL rrd_mpi_io_surface( 'surf%css', surf%css )
    IF ( ALLOCATED ( surf%cssws ) )  CALL rrd_mpi_io_surface( 'surf%cssws', surf%cssws )
    IF ( ALLOCATED ( surf%qcsws ) )  CALL rrd_mpi_io_surface( 'surf%qcsws', surf%qcsws )
    IF ( ALLOCATED ( surf%ncsws ) )  CALL rrd_mpi_io_surface( 'surf%ncsws', surf%ncsws )
    IF ( ALLOCATED ( surf%qisws ) )  CALL rrd_mpi_io_surface( 'surf%qisws', surf%qisws )
    IF ( ALLOCATED ( surf%nisws ) )  CALL rrd_mpi_io_surface( 'surf%nisws', surf%nisws )
    IF ( ALLOCATED ( surf%qrsws ) )  CALL rrd_mpi_io_surface( 'surf%qrsws', surf%qrsws )
    IF ( ALLOCATED ( surf%nrsws ) )  CALL rrd_mpi_io_surface( 'surf%nrsws', surf%nrsws )
    IF ( ALLOCATED ( surf%sasws ) )  CALL rrd_mpi_io_surface( 'surf%sasws', surf%sasws )
    IF ( ALLOCATED ( surf%usws )  )  CALL rrd_mpi_io_surface( 'surf%usws', surf%usws )
    IF ( ALLOCATED ( surf%vsws )  )  CALL rrd_mpi_io_surface( 'surf%vsws', surf%vsws )
    IF ( ALLOCATED ( surf%usvs )  )  CALL rrd_mpi_io_surface( 'surf%usvs', surf%usvs )
    IF ( ALLOCATED ( surf%vsus )  )  CALL rrd_mpi_io_surface( 'surf%vsus', surf%vsus )
    IF ( ALLOCATED ( surf%wsus_wsvs ) )                                                            &
       CALL rrd_mpi_io_surface( 'surf%wsus_wsvs', surf%wsus_wsvs )
    IF ( ALLOCATED ( surf%mom_flux_tke ) )                                                         &
       CALL rrd_mpi_io_surface( 'surf%mom_flux_tke', surf%mom_flux_tke )
!
!-- Redistribute surface elements on its respective type.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          surf_match_def  = surf_def%end_index(j,i) >= surf_def%start_index(j,i)
          surf_match_lsm  = surf_lsm%end_index(j,i) >= surf_lsm%start_index(j,i)
          surf_match_usm  = surf_usm%end_index(j,i) >= surf_usm%start_index(j,i)

          IF ( surf_match_def )  THEN
!
!--          Set the start index for the local surface element
             mm = surf_def%start_index(j,i)
!
!--          For index pair (j,i) on file loop from start to end index, and in case the local
!--          surface element mm is smaller than the local end index, assign the respective
!--          surface data to this element. Since all surfaces in the restart input array
!--          are gathered into a one-dimensional array and are distributed over different surface
!--          types, it need to be checked that input data is not mapped twice onto different
!--          surface types. Hence, store information of mapped surfaces in a didicated array and
!--          only write input surface data onto the target surface types that has not been
!--          written before.
             DO  m = surf%start_index(j,i), surf%end_index(j,i)
                IF ( surf_def%end_index(j,i) >= mm  .AND.  .NOT. mapped(m) )  THEN
                   CALL restore_surface_elements( surf_def, mm, surf, m )
                   mapped(m) = .TRUE.
                   mm = mm + 1
                ENDIF
             ENDDO
          ENDIF
!
!--       Natural- and urban-like horizontal surfaces.
          IF ( surf_match_lsm )  THEN
             mm = surf_lsm%start_index(j,i)
             DO  m = surf%start_index(j,i), surf%end_index(j,i)
                IF ( surf_lsm%end_index(j,i) >= mm  .AND.  .NOT. mapped(m) )  THEN
                   CALL restore_surface_elements( surf_lsm, mm, surf, m )
                   mapped(m) = .TRUE.
                   mm = mm + 1
                ENDIF
             ENDDO
          ENDIF
          IF ( surf_match_usm )  THEN
             mm = surf_usm%start_index(j,i)
             DO  m = surf%start_index(j,i), surf%end_index(j,i)
                IF ( surf_usm%end_index(j,i) >= mm  .AND.  .NOT. mapped(m) )  THEN
                   CALL restore_surface_elements( surf_usm, mm, surf, m )
                   mapped(m) = .TRUE.
                   mm = mm + 1
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO

!
!-- Treat model-top fluxes. Only proceed when top fluxes are defined, else
!-- rrd_mpi_io_global_array may lead to a segmentation fault.
    IF ( use_top_fluxes )  THEN
!
!--    Get total number of model-top points on the file.
       CALL rrd_mpi_io_global_array( 'ns_on_file_top', ns_on_file )

       IF ( ALLOCATED( surf_t%start_index ) )  CALL deallocate_surface_attributes_top( surf_t )

       ALLOCATE( surf_t%start_index(nys:nyn,nxl:nxr) )
       ALLOCATE( surf_t%end_index(nys:nyn,nxl:nxr) )
       surf_t%start_index = 0
       surf_t%end_index   = -1

       CALL rrd_mpi_io( 'global_start_index_top', global_start_index )
       CALL rrd_mpi_io( 'global_end_index_top' , global_end_index )

       CALL rd_mpi_io_surface_filetypes( surf_t%start_index, surf_t%end_index, data_to_read,       &
           global_start_index, global_end_index )
!
!--    Check if data is available. In contrast to "normal" surfaces, which are always available,
!--    model-top fluxes are not necessarily present.
       IF ( data_to_read )  THEN
          surf_t%ns = MAX( 2, MAXVAL( surf_t%end_index ) )

          CALL allocate_surface_attributes_top( surf_t, nys, nyn, nxl, nxr )

          IF ( ALLOCATED ( surf_t%usws )  )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%usws',  surf_t%usws  )
          ENDIF
          IF ( ALLOCATED ( surf_t%vsws )  )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%vsws',  surf_t%vsws  )
          ENDIF
          IF ( ALLOCATED ( surf_t%shf )   )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%shf',   surf_t%shf   )
          ENDIF
          IF ( ALLOCATED ( surf_t%qsws )  )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%qsws',  surf_t%qsws  )
          ENDIF
          IF ( ALLOCATED ( surf_t%ssws )  )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%ssws',  surf_t%ssws  )
          ENDIF
          IF ( ALLOCATED ( surf_t%cssws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%cssws', surf_t%cssws )
          ENDIF
          IF ( ALLOCATED ( surf_t%qcsws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%qcsws', surf_t%qcsws )
          ENDIF
          IF ( ALLOCATED ( surf_t%ncsws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%ncsws', surf_t%ncsws )
          ENDIF
          IF ( ALLOCATED ( surf_t%qrsws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%qrsws', surf_t%qrsws )
          ENDIF
          IF ( ALLOCATED ( surf_t%nrsws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%nrsws', surf_t%nrsws )
          ENDIF
          IF ( ALLOCATED ( surf_t%qisws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%qisws', surf_t%qisws )
          ENDIF
          IF ( ALLOCATED ( surf_t%nisws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%nisws', surf_t%nisws )
          ENDIF
          IF ( ALLOCATED ( surf_t%sasws ) )  THEN
             CALL rrd_mpi_io_surface( 'surf_top%sasws', surf_t%sasws )
          ENDIF
!
!--       Redistribute model-top surface elements.
          DO  i = nxl, nxr
             DO  j = nys, nyn
!
!--             Set the start index for the local surface element.
                mm = surf_top%start_index(j,i)

                DO  m = surf_t%start_index(j,i), surf_t%end_index(j,i)
                   IF ( surf_top%end_index(j,i) >= mm )                                            &
                      CALL restore_surface_elements( surf_top, mm, surf_t, m )
                   mm = mm + 1
                ENDDO
             ENDDO
          ENDDO
       ENDIF

    ENDIF

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Restores surface elements back on its respective type.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE restore_surface_elements( surf_target, m_target, surf_file, m_file )

    INTEGER(iwp) ::  m_file    !< respective surface-element index of current surface array
    INTEGER(iwp) ::  m_target  !< respecitve surface-element index of surface array on file
    INTEGER(iwp) ::  lsp       !< running index chemical species

    TYPE(surf_type) ::  surf_target  !< target surface type
    TYPE(surf_type) ::  surf_file    !< surface type on file

    IF ( ALLOCATED( surf_target%us )  .AND.  ALLOCATED( surf_file%us ) )                           &
       surf_target%us(m_target) = surf_file%us(m_file)

    IF ( ALLOCATED( surf_target%ol )  .AND.  ALLOCATED( surf_file%ol ) )                           &
       surf_target%ol(m_target) = surf_file%ol(m_file)

    IF ( ALLOCATED( surf_target%rib )  .AND.  ALLOCATED( surf_file%rib ) )                         &
       surf_target%rib(m_target) = surf_file%rib(m_file)

    IF ( ALLOCATED( surf_target%ol )  .AND.  ALLOCATED( surf_file%ol ) )                           &
       surf_target%ol(m_target) = surf_file%ol(m_file)

    IF ( ALLOCATED( surf_target%pt_surface )  .AND.  ALLOCATED( surf_file%pt_surface ) )           &
       surf_target%pt_surface(m_target) = surf_file%pt_surface(m_file)

    IF ( ALLOCATED( surf_target%q_surface )  .AND.  ALLOCATED( surf_file%q_surface ) )             &
       surf_target%q_surface(m_target) = surf_file%q_surface(m_file)

    IF ( ALLOCATED( surf_target%vpt_surface )  .AND.  ALLOCATED( surf_file%vpt_surface ) )         &
       surf_target%vpt_surface(m_target) = surf_file%vpt_surface(m_file)

    IF ( ALLOCATED( surf_target%ts )  .AND.  ALLOCATED( surf_file%ts ) )                           &
       surf_target%ts(m_target) = surf_file%ts(m_file)
!
!-- Do not restore constant surface fluxes because this would overwrite runtime stettings.
    IF ( .NOT. constant_heatflux )  THEN
       IF ( ALLOCATED( surf_target%shf )  .AND.  ALLOCATED( surf_file%shf ) )                      &
          surf_target%shf(m_target) = surf_file%shf(m_file)
    ENDIF

    IF ( ALLOCATED( surf_target%qs )  .AND.  ALLOCATED( surf_file%qs ) )                           &
       surf_target%qs(m_target) = surf_file%qs(m_file)

    IF ( ALLOCATED( surf_target%qsws )  .AND.  ALLOCATED( surf_file%qsws ) )                       &
       surf_target%qsws(m_target) = surf_file%qsws(m_file)

    IF ( ALLOCATED( surf_target%ss )  .AND.  ALLOCATED( surf_file%ss ) )                           &
       surf_target%ss(m_target) = surf_file%ss(m_file)

    IF ( ALLOCATED( surf_target%ssws )  .AND.  ALLOCATED( surf_file%ssws ) )                       &
       surf_target%ssws(m_target) = surf_file%ssws(m_file)

    IF ( ALLOCATED( surf_target%css )  .AND.  ALLOCATED( surf_file%css ) )  THEN
       DO  lsp = 1, nvar
          surf_target%css(lsp,m_target) = surf_file%css(lsp,m_file)
       ENDDO
    ENDIF

    IF ( ALLOCATED( surf_target%cssws )  .AND.  ALLOCATED( surf_file%cssws ) )  THEN
       DO  lsp = 1, nvar
          surf_target%cssws(lsp,m_target) = surf_file%cssws(lsp,m_file)
       ENDDO
    ENDIF

    IF ( ALLOCATED( surf_target%qcs )  .AND.  ALLOCATED( surf_file%qcs ) )                         &
      surf_target%qcs(m_target) = surf_file%qcs(m_file)

    IF ( ALLOCATED( surf_target%qcsws )  .AND.  ALLOCATED( surf_file%qcsws ) )                     &
       surf_target%qcsws(m_target) = surf_file%qcsws(m_file)

    IF ( ALLOCATED( surf_target%ncs )  .AND.  ALLOCATED( surf_file%ncs ) )                         &
       surf_target%ncs(m_target) = surf_file%ncs(m_file)

    IF ( ALLOCATED( surf_target%ncsws )  .AND.  ALLOCATED( surf_file%ncsws ) )                     &
       surf_target%ncsws(m_target) = surf_file%ncsws(m_file)

    IF ( ALLOCATED( surf_target%qis )  .AND.  ALLOCATED( surf_file%qis ) )                         &
      surf_target%qis(m_target) = surf_file%qis(m_file)

    IF ( ALLOCATED( surf_target%qisws )  .AND.  ALLOCATED( surf_file%qisws ) )                     &
       surf_target%qisws(m_target) = surf_file%qisws(m_file)

    IF ( ALLOCATED( surf_target%nis )  .AND.  ALLOCATED( surf_file%nis ) )                         &
       surf_target%nis(m_target) = surf_file%nis(m_file)

    IF ( ALLOCATED( surf_target%nisws )  .AND.  ALLOCATED( surf_file%nisws ) )                     &
       surf_target%nisws(m_target) = surf_file%nisws(m_file)

    IF ( ALLOCATED( surf_target%qrs )  .AND.  ALLOCATED( surf_file%qrs ) )                         &
       surf_target%qrs(m_target) = surf_file%qrs(m_file)

    IF ( ALLOCATED( surf_target%qrsws )  .AND.  ALLOCATED( surf_file%qrsws ) )                     &
       surf_target%qrsws(m_target) = surf_file%qrsws(m_file)

    IF ( ALLOCATED( surf_target%nrs )  .AND.  ALLOCATED( surf_file%nrs ) )                         &
       surf_target%nrs(m_target) = surf_file%nrs(m_file)

    IF ( ALLOCATED( surf_target%nrsws )  .AND.  ALLOCATED( surf_file%nrsws ) )                     &
       surf_target%nrsws(m_target) = surf_file%nrsws(m_file)

    IF ( ALLOCATED( surf_target%sasws )  .AND.  ALLOCATED( surf_file%sasws ) )                     &
       surf_target%sasws(m_target) = surf_file%sasws(m_file)

    IF ( ALLOCATED( surf_target%usws )  .AND.  ALLOCATED( surf_file%usws ) )                       &
       surf_target%usws(m_target) = surf_file%usws(m_file)

    IF ( ALLOCATED( surf_target%vsws )  .AND.  ALLOCATED( surf_file%vsws ) )                       &
       surf_target%vsws(m_target) = surf_file%vsws(m_file)

    IF ( ALLOCATED( surf_target%usvs )  .AND.  ALLOCATED( surf_file%usvs ) )                       &
       surf_target%usvs(m_target) = surf_file%usvs(m_file)

    IF ( ALLOCATED( surf_target%vsus )  .AND.  ALLOCATED( surf_file%vsus ) )                       &
       surf_target%vsus(m_target) = surf_file%vsus(m_file)

    IF ( ALLOCATED( surf_target%wsus_wsvs )  .AND.  ALLOCATED( surf_file%wsus_wsvs ) )             &
       surf_target%wsus_wsvs(m_target) = surf_file%wsus_wsvs(m_file)

    IF ( ALLOCATED( surf_target%mom_flux_tke )  .AND.                                              &
         ALLOCATED( surf_file%mom_flux_tke   ) )                                                   &
       surf_target%mom_flux_tke(0:1,m_target) = surf_file%mom_flux_tke(0:1,m_file)

 END SUBROUTINE restore_surface_elements

 END SUBROUTINE surface_rrd_local_mpi

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Counts the number of surface elements- This is required for reading and writing restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_last_actions
!
!-- Sum-up number of surfaces
    ns_on_file(1) = surf_def%ns + surf_lsm%ns + surf_usm%ns

 END SUBROUTINE surface_last_actions

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine maps surface data read from file after restart - 1D arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_restore_elements_1d( surf_target, surf_file, start_index_c,                    &
                                         start_index_on_file, end_index_on_file, nxlc, nysc, nxlf, &
                                         nxrf, nysf, nynf, nys_on_file, nyn_on_file, nxl_on_file,  &
                                         nxr_on_file )

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< running index along x-direction, refers to former domain size
    INTEGER(iwp) ::  ic    !< running index along x-direction, refers to current domain size
    INTEGER(iwp) ::  j     !< running index along y-direction, refers to former domain size
    INTEGER(iwp) ::  jc    !< running index along y-direction, refers to former domain size
    INTEGER(iwp) ::  m     !< surface-element index on file
    INTEGER(iwp) ::  mm    !< surface-element index on current subdomain
    INTEGER(iwp) ::  nxlc  !< index of left boundary on current subdomain
    INTEGER(iwp) ::  nxlf  !< index of left boundary on former subdomain
    INTEGER(iwp) ::  nxrf  !< index of right boundary on former subdomain
    INTEGER(iwp) ::  nysc  !< index of north boundary on current subdomain
    INTEGER(iwp) ::  nynf  !< index of north boundary on former subdomain
    INTEGER(iwp) ::  nysf  !< index of south boundary on former subdomain

    INTEGER(iwp) ::  nxl_on_file  !< leftmost index on file
    INTEGER(iwp) ::  nxr_on_file  !< rightmost index on file
    INTEGER(iwp) ::  nyn_on_file  !< northmost index on file
    INTEGER(iwp) ::  nys_on_file  !< southmost index on file

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index_c  !<
    INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) ::  start_index_on_file  !< start index of surface
                                                                                                      !< elements on file
    INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) ::  end_index_on_file    !< end index of surface
                                                                                                      !< elements on file

    REAL(wp), DIMENSION(:) ::  surf_target  !< target surface type
    REAL(wp), DIMENSION(:) ::  surf_file    !< surface type on file

    ic = nxlc
    DO  i = nxlf, nxrf
       jc = nysc
       DO  j = nysf, nynf
          mm = start_index_c(jc,ic)
          DO  m = start_index_on_file(j,i), end_index_on_file(j,i)
             surf_target(mm) = surf_file(m)
             mm = mm + 1
          ENDDO
          jc = jc + 1
       ENDDO
       ic = ic + 1
    ENDDO


 END SUBROUTINE surface_restore_elements_1d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine maps surface data read from file after restart - 2D arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_restore_elements_2d( surf_target, surf_file, start_index_c,                    &
                                         start_index_on_file, end_index_on_file, nxlc, nysc, nxlf, &
                                         nxrf, nysf, nynf, nys_on_file, nyn_on_file, nxl_on_file,  &
                                         nxr_on_file )

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< running index along x-direction, refers to former domain size
    INTEGER(iwp) ::  ic    !< running index along x-direction, refers to current domain size
    INTEGER(iwp) ::  j     !< running index along y-direction, refers to former domain size
    INTEGER(iwp) ::  jc    !< running index along y-direction, refers to former domain size
    INTEGER(iwp) ::  m     !< surface-element index on file
    INTEGER(iwp) ::  mm    !< surface-element index on current subdomain
    INTEGER(iwp) ::  nxlc  !< index of left boundary on current subdomain
    INTEGER(iwp) ::  nxlf  !< index of left boundary on former subdomain
    INTEGER(iwp) ::  nxrf  !< index of right boundary on former subdomain
    INTEGER(iwp) ::  nysc  !< index of north boundary on current subdomain
    INTEGER(iwp) ::  nynf  !< index of north boundary on former subdomain
    INTEGER(iwp) ::  nysf  !< index of south boundary on former subdomain

    INTEGER(iwp) ::  nxl_on_file  !< leftmost index on file
    INTEGER(iwp) ::  nxr_on_file  !< rightmost index on file
    INTEGER(iwp) ::  nyn_on_file  !< northmost index on file
    INTEGER(iwp) ::  nys_on_file  !< southmost index on file

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index_c !< start index of surface type

    INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) ::  start_index_on_file  !< start index of surface
                                                                                                      !< elements on file
    INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) ::  end_index_on_file    !< end index of surface
                                                                                                      !< elements on file

    REAL(wp), DIMENSION(:,:) ::  surf_target !< target surface type
    REAL(wp), DIMENSION(:,:) ::  surf_file   !< surface type on file

    ic = nxlc
    DO  i = nxlf, nxrf
       jc = nysc
       DO  j = nysf, nynf
          mm = start_index_c(jc,ic)
          DO  m = start_index_on_file(j,i), end_index_on_file(j,i)
             surf_target(:,mm) = surf_file(:,m)
             mm = mm + 1
          ENDDO
          jc = jc + 1
       ENDDO
       ic = ic + 1
    ENDDO

 END SUBROUTINE surface_restore_elements_2d

 END MODULE surface_data_handling
