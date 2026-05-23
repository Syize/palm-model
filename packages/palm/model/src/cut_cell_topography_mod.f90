!> @file cut_cell_topography_mod.f90
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
!> This module contains all routines to set-up and represent cut-cell topography.
!--------------------------------------------------------------------------------------------------!
 MODULE cut_cell_topography_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               dzu,                                                                                &
               x,                                                                                  &
               xu,                                                                                 &
               y,                                                                                  &
               yv,                                                                                 &
               zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  degrees_to_radiants,                                                                &
               pi,                                                                                 &
               radiants_to_degrees

    USE boundary_settings_mod,                                                                     &
        ONLY:  set_lateral_neumann_bc

    USE control_parameters

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz,                                                                     &
               exchange_horiz_2d,                                                                  &
               exchange_horiz_2d_byte,                                                             &
               exchange_horiz_int

    USE general_utilities,                                                                         &
        ONLY:  normalize_vector

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  buildings_f,                                                                        &
               building_id_f,                                                                      &
               char_fill,                                                                          &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               init_model,                                                                         &
               input_file_static,                                                                  &
               input_pids_static,                                                                  &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               int_2d_8bit,                                                                        &
               list_building_ids,                                                                  &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               pids_id,                                                                            &
               terrain_height_f,                                                                   &
               vars_pids

    USE pegrid

    USE surface_mod,                                                                               &
        ONLY:  surf_cct,                                                                           &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_type,                                                                          &
               surf_u,                                                                             &
               surf_usm,                                                                           &
               surf_v,                                                                             &
               surf_w

    IMPLICIT NONE

    TYPE face_data

       INTEGER(iwp) ::  dim_3d = 3               !< number of spatial dimensions
       INTEGER(iwp) ::  dim_vertex_coords  = 4   !< number of vertex coordinates per grid cell
       INTEGER(iwp) ::  dim_vertex_shifts  = 1   !< number of vertex shifts
       INTEGER(iwp) ::  fill_value_int = -1      !< fill value for empty vertexes in polygon
       INTEGER(iwp) ::  num_faces = HUGE(1_iwp)  !< number of cut-cell faces in static input file
       INTEGER(iwp) ::  num_faces_vert = 7       !< maximum number of vertices describing a cut-cell face (set to 7)
       INTEGER(iwp) ::  num_vert = HUGE(1_iwp)   !< number of vertices in domain

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_id     !< corresponding building ID of surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_types  !< corresponding building type of surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_edges    !< list of number of edges per grid cell
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  pav_types    !< corresponding pavement type of surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  types        !< type of surface 1 = lsm, 0 = usm, 3 = def
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  veg_types    !< corresponding vegetation type of surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  wat_types    !< corresponding water type of surface

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  faces          !< list of all cut-cell surfaces, linking vertex number to plane
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  kji            !< (k,j,i)-location of cut-cell surface on finite grid
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  vertex_coords  !< list of all vertex coordinates for radiation

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  area           !< surface area of each cut-cell surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vertex_shifts  !< vertex shift

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  centers   !< mass center coordinate of each cut-cell surfface
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  normals   !< normal vector of each cut-cell surface
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vertices  !< list of all vertex coordinates

    END TYPE face_data

    TYPE face_data_kji

       LOGICAL ::  contain_faces = .FALSE.  !< flag to indicate whether a grid box compasses any cut-cell surface
       LOGICAL ::  def_surface = .FALSE.    !< flag to indicate a default-type cut-cell surface
       LOGICAL ::  face_accessed = .FALSE.  !< flag to indicate that the cut-cell face is accessed and an energy-balance needs to be solved
       LOGICAL ::  lsm_surface = .FALSE.    !< flag to indicate a natural-type cut-cell surface
       LOGICAL ::  roof_surface = .FALSE.   !< flag to indicate a roof cut-cell surface (usm surface)
       LOGICAL ::  wall_surface = .FALSE.   !< flag to indicate a wall cut-cell surface (usm surface)

       INTEGER(iwp) ::  bid             !< building ID
       INTEGER(iwp) ::  btype           !< building type
       INTEGER(iwp) ::  num_access = 0  !< number of LSM and USM surfaces that access the cut-cell face
       INTEGER(iwp) ::  num_edges       !< number of edges per grid cell
       INTEGER(iwp) ::  nv              !< number of vertices per face
       INTEGER(iwp) ::  ptype           !< pavement type
       INTEGER(iwp) ::  vtype           !< vegetation type
       INTEGER(iwp) ::  wtype           !< water type

       INTEGER(iwp), DIMENSION(5) ::  index_access = 0  !< corresponding indices of surf_lsm / surf_usm that access the cut-cell surface

       REAL(wp) ::  area  !< surface area per grid box

       REAL(wp), DIMENSION(3) ::  center      !< mass center coordinates
       REAL(wp), DIMENSION(3) ::  normal      !< normal vector
       REAL(wp), DIMENSION(7) ::  faces = -1  !< vertex numbers

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vertex_coord  !< vertex coordinates for each face

    END TYPE face_data_kji

    TYPE( face_data ) ::  cct           !< static input data but only relevant for local subdomain
    TYPE( face_data ) ::  cct_global_f  !< global input from static driver file, already preprocessed

    TYPE( face_data_kji ), DIMENSION(:,:,:), ALLOCATABLE ::  surf_data  !< cut-cell surface data stored on 3D array (required for processing)

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  cct_vertex_coords  !< indices for vertex coordinates

    INTEGER(iwp), DIMENSION(3,0:5), PARAMETER ::  vertex_shift_nvect = RESHAPE( (/  1, 0, 0,       &
                                                                                   -1, 0, 0,       &
                                                                                    0, 1, 0,       &
                                                                                    0,-1, 0,       &
                                                                                    0, 0, 1,       &
                                                                                    0, 0,-1 /),    &
                                                                                (/ 3, 6 /) )  !< k,j,i of vertex shift

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  cct_vertex_shifts  !< vertex shifts


    SAVE

    PRIVATE
!
!-- Public subroutines.
    PUBLIC cct_check_parameters,                                                                   &
           cct_define_topography,                                                                  &
           cct_init,                                                                               &
           cct_input,                                                                              &
           cct_to_surface_types,                                                                   &
           surface_types_to_cct
!
!-- Public variables.
    PUBLIC cct,                                                                                    &
           cct_vertex_coords,                                                                      &
           cct_vertex_shifts,                                                                      &
           vertex_shift_nvect

    INTERFACE cct_check_parameters
       MODULE PROCEDURE cct_check_parameters
    END INTERFACE cct_check_parameters

    INTERFACE cct_define_topography
       MODULE PROCEDURE cct_define_topography
    END INTERFACE cct_define_topography

    INTERFACE cct_init
       MODULE PROCEDURE cct_init
    END INTERFACE cct_init

    INTERFACE cct_input
       MODULE PROCEDURE cct_input
    END INTERFACE cct_input

    INTERFACE cct_to_surface_types
       MODULE PROCEDURE cct_to_surface_types
    END INTERFACE cct_to_surface_types

    INTERFACE surface_types_to_cct
       MODULE PROCEDURE surface_types_to_cct
    END INTERFACE surface_types_to_cct

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Check parameter settings with respect to the cut-cell topography.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE cct_check_parameters

    IF ( .NOT. allow_roughness_limitation )  THEN
       message_string = 'cut cell topography does not allow_roughness_limitation = .FALSE., &' //  &
                        'will be reset to .TRUE.'
       CALL message( 'cct_check_parameters', 'PAC0362', 0, 1, 0, 6, 0 )
       allow_roughness_limitation = .TRUE.
    ENDIF

 END SUBROUTINE cct_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Define topography and classify grid points into terrain and building based on cut-cell input.
!> Also, prepare surface and vertex data to be processed on a 3d grid and later also for surface
!> classification.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE cct_define_topography( topo )

    INTEGER(iwp) ::  i     !< running index x-direction
    INTEGER(iwp) ::  iref  !< x-index of possibly accessing grid point
    INTEGER(iwp) ::  is    !< x-index of treated cut-cell surface
    INTEGER(iwp) ::  j     !< running index y-direction
    INTEGER(iwp) ::  jref  !< y-index of possibly accessing grid point
    INTEGER(iwp) ::  js    !< y-index of treated cut-cell surface
    INTEGER(iwp) ::  k     !< running index z-direction
    INTEGER(iwp) ::  kref  !< z-index of possibly accessing grid point
    INTEGER(iwp) ::  l     !< running index over different possible access directions
    INTEGER(iwp) ::  m     !< running index over all vertices defining a face
    INTEGER(iwp) ::  n     !< running index over all faces
    INTEGER(iwp) ::  vn    !< vertex number

    INTEGER(iwp), DIMENSION(0:6) ::  off_i = (/  0,  0,  0,  0, -1,  1,  0 /)  !< offset indicies between reference and surface grid point in x, for
                                                                               !< upward, downward, northward, southward, eastward, westward, model top
    INTEGER(iwp), DIMENSION(0:6) ::  off_j = (/  0,  0, -1,  1,  0,  0,  0 /)  !< offset indicies between reference and surface grid point in y
    INTEGER(iwp), DIMENSION(0:6) ::  off_k = (/ -1,  1,  0,  0,  0,  0,  1 /)  !< offset indicies between reference and surface grid point in z

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo  !< passed topography flag array

    LOGICAL ::  found  !< flag to check if grid point has been classified already

    REAL(wp) ::  dd        !< linear equation coefficient
    REAL(wp) ::  dist      !< distance to spanned surface plane approximation
    REAL(wp) ::  dot_norm  !< normalized scalar product
    REAL(wp) ::  normx     !< x-component of normal vector
    REAL(wp) ::  normy     !< x-component of normal vector
    REAL(wp) ::  normz     !< x-component of normal vector
    REAL(wp) ::  xs        !< x-coordinate of face center
    REAL(wp) ::  ys        !< y-coordinate of face center
    REAL(wp) ::  zs        !< z-coordinate of face center

    REAL(wp), DIMENSION(3) ::  dist_max  !< maximum possible distance between prognostic grid point and surface
    REAL(wp), DIMENSION(3) ::  norm      !< direction vector describing the grid line along which a surface is accessed in the diffusion
    REAL(wp), DIMENSION(3) ::  p_0       !< coordinate vector of the prognostic point
    REAL(wp), DIMENSION(3) ::  p_i       !< interception point between access line and cut-cell face


!
!-- Allocate 3d array where all cut-cell information required for topography processing is stored
!-- on. This is usually a sparsely-dense 3D array defined over ghost layers in order to avoid
!-- later exchange of ghost points.
    ALLOCATE( surf_data(nzb:nzt,nysg:nyng,nxlg:nxrg) )

    DO  i = nxl-1, nxr+1
       DO  j = nys-1, nyn+1
          DO  k = nzb, nzt
!
!--          For each (k,j,i) grid point search for the corresponding cut-cell surface. There is
!--          maximum one cut-cell surface per grid cell defined.
             DO  n = 1, cct%num_faces
                IF ( k == cct%kji(1,n)  .AND.  j == cct%kji(2,n)  .AND.  i == cct%kji(3,n) )  THEN
!
!--                Set flag to indicate that the grid cell contains a cut-cell surface.
                   surf_data(k,j,i)%contain_faces = .TRUE.
!
!--                If an energy balance model is employed, distinguish between its general type.
!--                type = 0: natural-type surface (LSM), type = 1 or 2: roof or wall surface,
!--                respectively (USM), type = 3: default surface.
                   IF ( land_surface  .OR.  urban_surface )  THEN
                      surf_data(k,j,i)%lsm_surface  = MERGE( .TRUE., .FALSE., cct%types(n) == 0 )
                      surf_data(k,j,i)%wall_surface = MERGE( .TRUE., .FALSE., cct%types(n) == 1 )
                      surf_data(k,j,i)%roof_surface = MERGE( .TRUE., .FALSE., cct%types(n) == 2 )
                      surf_data(k,j,i)%def_surface  = MERGE( .TRUE., .FALSE., cct%types(n) == 3 )
!
!--                   Furthermore, store the corresponding specific surface type (pavement,
!--                   vegetation, water, building).
                      IF ( land_surface )  THEN
                         surf_data(k,j,i)%ptype = cct%pav_types(n)
                         surf_data(k,j,i)%vtype = cct%veg_types(n)
                         surf_data(k,j,i)%wtype = cct%wat_types(n)
                      ENDIF
                      IF ( urban_surface )  THEN
                         surf_data(k,j,i)%btype = cct%build_types(n)
                      ENDIF
                      surf_data(k,j,i)%bid = cct%build_id(n)

                      IF ( .NOT. surf_data(k,j,i)%lsm_surface  .AND.                               &
                           .NOT. surf_data(k,j,i)%wall_surface  .AND.                              &
                           .NOT. surf_data(k,j,i)%roof_surface )                                   &
                      THEN
                         surf_data(k,j,i)%def_surface = .TRUE.
                      ENDIF
                   ELSE
                      surf_data(k,j,i)%def_surface = .TRUE.
                   ENDIF

                   surf_data(k,j,i)%faces = cct%faces(:,n)
!
!--                Count the number of vertices that span the cut-cell surface.
                   surf_data(k,j,i)%nv = 0
                   DO  m = LBOUND( cct%faces, DIM = 1 ), UBOUND( cct%faces, DIM = 1 )
                      IF ( cct%faces(m,n) /= -1 )  surf_data(k,j,i)%nv = surf_data(k,j,i)%nv + 1
                   ENDDO
!
!--                Store further attributes.
                   surf_data(k,j,i)%num_edges = surf_data(k,j,i)%nv

                   surf_data(k,j,i)%area   = cct%area(n)
                   surf_data(k,j,i)%normal = cct%normals(:,n)
                   surf_data(k,j,i)%center = cct%centers(:,n)

                   ALLOCATE( surf_data(k,j,i)%vertex_coord(1:3,1:surf_data(k,j,i)%nv) )
                   DO  m = 1, surf_data(k,j,i)%nv
                      vn = cct%faces(m,n)
                      surf_data(k,j,i)%vertex_coord(:,m) = cct%vertices(:,vn)
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Clear topography array for the scalar-, u-, v-, and w-grid at the beginning.
    topo(:,:,:) = IBCLR( topo(:,:,:), 0 )
    topo(:,:,:) = IBCLR( topo(:,:,:), 1 )
    topo(:,:,:) = IBCLR( topo(:,:,:), 2 )
    topo(:,:,:) = IBCLR( topo(:,:,:), 3 )
!
!-- Clear the bits for general classification into terrain and building.
    topo(:,:,:) = IBCLR( topo(:,:,:), 11 )
    topo(:,:,:) = IBCLR( topo(:,:,:), 12 )
!
!-- Clear the bits used for the advection scheme.
    topo = IBCLR( topo, 7  )
    topo = IBCLR( topo, 9  )
    topo = IBCLR( topo, 10 )
    topo = IBCLR( topo, 26 )
    topo = IBCLR( topo, 27 )
    topo = IBCLR( topo, 28 )
!
!-- Now check all grid points on the staggerd grid and check if they lie beyond or beneath the
!-- spanned cut-cell surface. Note, the cut-cell surface is usually curved, so that we
!-- need to approximate by its mass center and the average normal vector to define a straight
!-- plane, which is used to check if a grid point lies inside or outside the topography.
!-- Based on this, modify topography flags accordingly.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             DO  n = 1, cct%num_faces
                js = cct%kji(2,n)
                is = cct%kji(3,n)

                IF ( is == i  .AND.  js == j )  THEN
                   zs = cct%centers(1,n)
                   ys = cct%centers(2,n)
                   xs = cct%centers(3,n)

                   normz = cct%normals(1,n)
                   normy = cct%normals(2,n)
                   normx = cct%normals(3,n)
!
!--                Check scalar grid. Therefore, use the general form of a plane equation defined
!--                by its average normal vector and plug-in the grid-point coordinate.
!--                If dist <= 0.0, the grid point lies inside topography.
                   dist = x(i) * normx + y(j) * normy + zu(k) * normz -                            &
                          ( normx * xs + normy * ys + normz * zs )

                   IF ( dist <= 0.0_wp )  THEN
                      topo(k,j,i) = IBSET( topo(k,j,i), 0 )
!
!--                   Classify grid point into terrain or building.
                      IF ( cct%types(n) == 0 )  THEN
                         topo(k,j,i) = IBSET( topo(k,j,i), 11 )
                         topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                      ELSEIF ( cct%types(n) == 1  .OR.  cct%types(n) == 2 )  THEN
                         topo(k,j,i) = IBCLR( topo(k,j,i), 11 )
                         topo(k,j,i) = IBSET( topo(k,j,i), 12 )
                      ELSE
                         topo(k,j,i) = IBSET( topo(k,j,i), 11 )
                         topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                      ENDIF
                   ELSE
                      topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                      topo(k,j,i) = IBCLR( topo(k,j,i), 11 )
                      topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                   ENDIF
!
!--                Check w-grid. If dist <= 0.0, the w-grid point lies beneath the surface.
                   dist = x(i) * normx + y(j) * normy + zw(k) * normz -                            &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  THEN
                      topo(k,j,i) = IBSET( topo(k,j,i), 3 )
                   ELSE
                      topo(k,j,i) = IBCLR( topo(k,j,i), 3 )
                   ENDIF
!
!--                Check u-grid. If dist <= 0.0, the u-grid point lies beneath the surface.
                   dist = xu(i) * normx + y(j) * normy + zu(k) * normz -                           &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  THEN
                      topo(k,j,i) = IBSET( topo(k,j,i), 1 )
                   ELSE
                      topo(k,j,i) = IBCLR( topo(k,j,i), 1 )
                   ENDIF
!
!--                Check v-grid. If dist <= 0.0, the v-grid point lies beneath the surface.
                   dist = x(i) * normx + yv(j) * normy + zu(k) * normz -                           &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  THEN
                      topo(k,j,i) = IBSET( topo(k,j,i), 2 )
                   ELSE
                      topo(k,j,i) = IBCLR( topo(k,j,i), 2 )
                   ENDIF
!
!--                Check the locations where advection fluxes for u in y- and z-direction are
!--                defined.
                   dist = xu(i) * normx + yv(j) * normy + zu(k) * normz -                          &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  topo(k,j,i) = IBSET( topo(k,j,i), 7 )

                   dist = xu(i) * normx + y(j) * normy + zw(k) * normz -                           &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  topo(k,j,i) = IBSET( topo(k,j,i), 9 )
!
!--                Check the locations where advection fluxes of v in x- and z-direction are
!--                defined.
                   dist = xu(i) * normx + yv(j) * normy + zu(k) * normz -                          &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  topo(k,j,i) = IBSET( topo(k,j,i), 10 )

                   dist = x(i) * normx + yv(j) * normy + zw(k) * normz -                           &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  topo(k,j,i) = IBSET( topo(k,j,i), 26 )
!
!--                Check the locations where advection fluxes of w in x- and y-direction are
!--                defined.
                   dist = xu(i) * normx + y(j) * normy + zw(k) * normz -                           &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  topo(k,j,i) = IBSET( topo(k,j,i), 27 )

                   dist = x(i) * normx + yv(j) * normy + zw(k) * normz -                           &
                          ( normx * xs + normy * ys + normz * zs )
                   IF ( dist <= 0.0_wp )  topo(k,j,i) = IBSET( topo(k,j,i), 28 )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Classify topography grid points in terrain an building grid points accordingly. In some cases
!-- the above written algorithm leads to wrong detections so that a building grid point is wrongly
!-- classified as terrain or vice versa. This leads to the situation that corresponding surface
!-- types are not correctly initialized. In order to improve the classification scheme,
!-- run over all topography-classified grid points and check if they will be assessed by a
!-- surface. This case, take the nearby cut-cell surface classification and classify the topography
!-- grid point accordingly.
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Topography grid point at k=nzb is not considered here, but have been already considered
!--       in the above written classification scheme.
          DO  k = nzb+1, nzt
!
!--          Check if current gridpoint belongs to topography.
             IF ( BTEST( topo(k,j,i), 0 ) )  THEN
!
!--             Run loop over all directions.
                DO  l = 0, 6
!
!--                Grid indicies of potential atmosphere grid point.
                   kref = k - off_k(l)
                   jref = j - off_j(l)
                   iref = i - off_i(l)
!
!--                Check if possibly assessing grid point belongs to atmosphere.
                   IF ( .NOT. BTEST( topo(kref,jref,iref), 0 ) )  THEN
                      found = .FALSE.
!
!--                   Define proxy for normal vector representing the grid-line direction.
                      norm = (/ -off_k(l), -off_j(l), -off_i(l) /)
!
!--                   Coordinate vector of assessing atmosphere grid point.
                      p_0 = (/ zu(kref), y(jref), x(iref) /)
!
!--                   Check if there are any cut-cell faces in (kref,jref,iref)-grid box.
                      IF ( surf_data(kref,jref,iref)%contain_faces )  THEN

                         dist_max = (/ 0.5_wp * dx, 0.5_wp * dy, 0.5_wp * dz(1) /)
!
!--                      Compute interception point of grid-line
!--                      (access direction and cut-cell surface).
                         dot_norm = DOT_PRODUCT( norm, surf_data(kref,jref,iref)%normal )

                         dd = 1.0E+8_wp
                         IF ( dot_norm /= 0 )  THEN
                            dd = DOT_PRODUCT( surf_data(kref,jref,iref)%center - p_0,              &
                                              surf_data(kref,jref,iref)%normal ) / dot_norm
                         ENDIF

                         p_i = p_0 + norm * dd

                         dist = SQRT( SUM( ( p_0 - p_i )**2 ) )
!
!--                      If the interception point lies within the same grid box (plus safety range
!--                      because of non-exact representation of normal vector and center
!--                      coordinates).
                         IF ( dist <= 1.5_wp * SQRT( SUM( ( norm * dist_max )**2 ) )  .AND.        &
                              dist > 0.0_wp )                                                      &
                         THEN

                            IF ( surf_data(kref,jref,iref)%lsm_surface )  THEN
                               topo(k,j,i) = IBSET( topo(k,j,i), 11 )
                               topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                            ELSEIF ( surf_data(kref,jref,iref)%wall_surface  .OR.                  &
                                     surf_data(kref,jref,iref)%roof_surface )                      &
                            THEN
                               topo(k,j,i) = IBCLR( topo(k,j,i), 11 )
                               topo(k,j,i) = IBSET( topo(k,j,i), 12 )
                            ELSE
                               topo(k,j,i) = IBSET( topo(k,j,i), 11 )
                               topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                            ENDIF
                            found = .TRUE.

                         ENDIF

                      ENDIF

                      IF ( surf_data(k,j,i)%contain_faces  .AND.  .NOT. found )  THEN

                         dist_max = 5.0_wp * (/ dx, dy, dz(1) /)

                         dot_norm = DOT_PRODUCT( norm, surf_data(k,j,i)%normal )

                         dd = 1.0E+8_wp
                         IF ( dot_norm /= 0 )  THEN
                            dd = DOT_PRODUCT( surf_data(k,j,i)%center - p_0,                       &
                                              surf_data(k,j,i)%normal ) / dot_norm
                         ENDIF

                         p_i = p_0 + norm * dd

                         dist = SQRT( SUM( ( p_0 - p_i )**2 ) )

                         IF ( dist <= SQRT( SUM( ( norm * dist_max )**2 ) )  .AND.  dist > 0.0_wp )&
                         THEN

                            IF ( surf_data(k,j,i)%lsm_surface )  THEN
                               topo(k,j,i) = IBSET( topo(k,j,i), 11 )
                               topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                            ELSEIF ( surf_data(k,j,i)%wall_surface  .OR.                           &
                                     surf_data(k,j,i)%roof_surface )                               &
                            THEN
                               topo(k,j,i) = IBCLR( topo(k,j,i), 11 )
                               topo(k,j,i) = IBSET( topo(k,j,i), 12 )
                            ELSE
                               topo(k,j,i) = IBSET( topo(k,j,i), 11 )
                               topo(k,j,i) = IBCLR( topo(k,j,i), 12 )
                            ENDIF
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- Further conditions for the staggered momemtum grid.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
!
!--          In case scalar grid at i-1 and i are both atmosphere, then u-grid is also
!--          within the atmosphere.
!--          u-grid.
             IF ( .NOT. BTEST( topo(k,j,i), 0 )  .AND.  .NOT. BTEST( topo(k,j,i-1), 0 ) )          &
                topo(k,j,i) = IBCLR( topo(k,j,i), 1 )
!
!--          v grid.
             IF ( .NOT. BTEST( topo(k,j,i), 0 )  .AND.  .NOT. BTEST( topo(k,j-1,i), 0 ) )          &
                topo(k,j,i) = IBCLR( topo(k,j,i), 2 )
!
!--          Vice versa, in case scalar grid at i-1 and i are both topography, then u-grid is also
!--          topography.
!--          u-grid.
             IF ( BTEST( topo(k,j,i), 0 )  .AND.  BTEST( topo(k,j,i-1), 0 ) )                      &
                topo(k,j,i) = IBSET( topo(k,j,i), 1 )
!
!--          v grid.
             IF ( BTEST( topo(k,j,i), 0 )  .AND.  BTEST( topo(k,j-1,i), 0 ) )                      &
                topo(k,j,i) = IBSET( topo(k,j,i), 2 )
          ENDDO

          DO k = nzb, nzt
!
!--         w grid.
            IF ( .NOT. BTEST( topo(k,j,i), 0 )  .AND.  .NOT. BTEST( topo(k+1,j,i), 0 ) )           &
               topo(k,j,i) = IBCLR( topo(k,j,i), 3 )
            IF ( BTEST( topo(k,j,i), 0 )  .AND.  BTEST( topo(k+1,j,i), 0 ) )                       &
               topo(k,j,i) = IBSET( topo(k,j,i), 3 )
          ENDDO
       ENDDO
    ENDDO

    CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )

 END SUBROUTINE cct_define_topography


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Initialize the cut-cell surfaces.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE cct_init

    INTEGER(iwp) ::  i               !< running index x-direction
    INTEGER(iwp) ::  j               !< running index y-direction
    INTEGER(iwp) ::  k               !< running index z-direction
    INTEGER(iwp) ::  m               !< running index over non-active cut-cell surfaces
    INTEGER(iwp) ::  mm              !< running index over active LSM or USM surfaces
    INTEGER(iwp) ::  n_accessed_lsm  !< number of LSM surfaces that access a cut-cell surface
    INTEGER(iwp) ::  n_accessed_usm  !< number of LSM surfaces that access a cut-cell surface
    INTEGER(iwp) ::  num_cc_kji      !< number of non-active cut-cell surfaces on (j,i)-column
    INTEGER(iwp) ::  start_index_cc  !< start index of non-active cut-cell surfaces on (j,i)-column

    LOGICAL ::  check_mismatch  !< error flag

    REAL(wp) ::  dist      !< actual distance between non-active cut-cell surface and active surface
    REAL(wp) ::  dist_min  !< lowest distance between non-active cut-cell surface and active surface

    TYPE( surf_type ), POINTER ::  surf !< surface data type


    IF ( land_surface  .OR.  urban_surface )  THEN
!
!--    Allocate and initialize index arrays required in the LSM + USM structures in order to
!--    link the 3D grid to the surfaces. This is required for exchanging ghost-point data on
!--    surfaces.
       IF ( surf_lsm%ns > 0 )  THEN
          ALLOCATE( surf_lsm%i_cc(1:surf_lsm%ns) )
          ALLOCATE( surf_lsm%j_cc(1:surf_lsm%ns) )
          ALLOCATE( surf_lsm%k_cc(1:surf_lsm%ns) )
          surf_lsm%i_cc = nxl
          surf_lsm%j_cc = nys
          surf_lsm%k_cc = nzb

          ALLOCATE( surf_lsm%accessed(1:surf_lsm%ns) )
          surf_lsm%accessed = .FALSE.
       ENDIF

       IF ( surf_usm%ns > 0 )  THEN
          ALLOCATE( surf_usm%i_cc(1:surf_usm%ns) )
          ALLOCATE( surf_usm%j_cc(1:surf_usm%ns) )
          ALLOCATE( surf_usm%k_cc(1:surf_usm%ns) )
          surf_usm%i_cc = nxl
          surf_usm%j_cc = nys
          surf_usm%k_cc = nzb

          ALLOCATE( surf_usm%accessed(1:surf_usm%ns) )
          surf_usm%accessed = .FALSE.
       ENDIF

    ENDIF
!
!-- First, store normal vector of the cut-cell surfaces on the respective surface elements and
!-- the distance between the surface and the corresponding atmosphere grid point. Check this
!-- surface-element-wise.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          First, for the classical surfaces defined at the s-grid.
             surf => surf_def
             CALL compute_interception( ( i + 0.5_wp ) * dx, ( j + 0.5_wp ) * dy, zu(k),           &
                                        's', .FALSE. )
             surf => surf_lsm
             CALL compute_interception( ( i + 0.5_wp ) * dx, ( j + 0.5_wp ) * dy, zu(k),           &
                                        's', .TRUE., 'lsm' )
             surf => surf_usm
             CALL compute_interception( ( i + 0.5_wp ) * dx, ( j + 0.5_wp ) * dy, zu(k),           &
                                        's', .TRUE., 'usm' )
!
!--          Now for the momemtum surfaces defined at the staggered velocity grid.
             surf => surf_u
             CALL compute_interception( i * dx, ( j + 0.5_wp ) * dy, zu(k), 'u', .FALSE. )

             surf => surf_v
             CALL compute_interception( ( i + 0.5_wp ) * dx, j * dy, zu(k), 'v', .FALSE. )

             surf => surf_w
             CALL compute_interception( ( i + 0.5_wp ) * dx, ( j + 0.5_wp ) * dy, zw(k),           &
                                        'w', .FALSE. )
          ENDDO
!
!--       Set and modify surface properties.
          surf => surf_def
          CALL cct_set_surface_properties

          surf => surf_lsm
          CALL cct_set_surface_properties

          surf => surf_usm
          CALL cct_set_surface_properties

          surf => surf_u
          CALL cct_set_surface_properties

          surf => surf_v
          CALL cct_set_surface_properties

          surf => surf_w
          CALL cct_set_surface_properties
       ENDDO
    ENDDO
!
!-- In case if an energy-balance model is employed, also the radiation model is employed and
!-- requires a surface temperature at each considered surface. However, the number of surf_lsm +
!-- surf_usm surfaces is not necessarily identical with the number of cut-cell surfaces used
!-- in the radiative transfer model to ensure a closed surface (polygon).
!-- Hereafter, lsm and usm surfaces are defined as active surfaces, as the resulting flux is
!-- feed-back to the prognostic equations, while the remaining cut-cell surfaces are not directly
!-- accessed by a prognostic grid point in the diffusion routines. However, as these surfaces
!-- need to be treated by the radiative transfer model, they also need a surface temperature
!-- and basic initialization of radiation properties.
    IF ( land_surface  .OR.  urban_surface )  THEN
!
!--    First, count the number of all cut-cell surfaces.
       surf_cct%ns = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
                IF ( surf_data(k,j,i)%contain_faces )  THEN
                   surf_cct%ns = surf_cct%ns + 1
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Allocate arrays required in surface data structure which contains all cut-cell surfaces,
!--    independent on whether there is an energy balance solved (LSM or USM) or not.
       ALLOCATE( surf_cct%start_index(nys:nyn,nxl:nxr) )
       ALLOCATE( surf_cct%end_index(nys:nyn,nxl:nxr)   )
       surf_cct%start_index = 0
       surf_cct%end_index   = -1

       ALLOCATE( surf_cct%i(1:surf_cct%ns) )
       ALLOCATE( surf_cct%j(1:surf_cct%ns) )
       ALLOCATE( surf_cct%k(1:surf_cct%ns) )

       ALLOCATE( surf_cct%m_index_ref(1:surf_cct%ns) )
       ALLOCATE( surf_cct%num_access(1:surf_cct%ns) )

       ALLOCATE ( surf_cct%n_s(1:surf_cct%ns,1:3) )

       ALLOCATE( surf_cct%face_accessed(1:surf_cct%ns) )
       ALLOCATE( surf_cct%index_access(1:surf_cct%ns,5) )
       ALLOCATE( surf_cct%pt_surface(1:surf_cct%ns) )
       ALLOCATE( surf_cct%albedo(1:surf_cct%ns) )
       ALLOCATE( surf_cct%emissivity(1:surf_cct%ns) )

       ALLOCATE ( surf_cct%faces(cct%num_faces_vert,1:surf_cct%ns) )
       ALLOCATE ( surf_cct%num_edges(1:surf_cct%ns) )
       ALLOCATE ( surf_cct%area(1:surf_cct%ns) )
       ALLOCATE ( surf_cct%wall_location_c(cct%dim_3d,1:surf_cct%ns) )

       ALLOCATE( surf_cct%accessed_lsm(1:surf_cct%ns) )
       ALLOCATE( surf_cct%accessed_usm(1:surf_cct%ns) )
!
!--    Allocate arrays required for the RTM data between cut-cell and LSM/USM surfaces.
       ALLOCATE( surf_cct%rad_lw_dif(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_lw_in(1:surf_cct%ns)  )
       ALLOCATE( surf_cct%rad_lw_out(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_lw_ref(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_lw_res(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_net(1:surf_cct%ns)    )
       ALLOCATE( surf_cct%rad_sw_dif(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_sw_dir(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_sw_in(1:surf_cct%ns)  )
       ALLOCATE( surf_cct%rad_sw_out(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_sw_ref(1:surf_cct%ns) )
       ALLOCATE( surf_cct%rad_sw_res(1:surf_cct%ns) )


       surf_cct%num_access(1:surf_cct%ns) = 0
!
!--    Transfer 3D cut-cell information onto dedicated surf-type array used in the RTM.
       start_index_cc = 1
       m = 1
       DO  i = nxl, nxr
          DO  j = nys, nyn

             num_cc_kji = 0
             DO  k = nzb, nzt
                IF ( surf_data(k,j,i)%contain_faces )  THEN
                   num_cc_kji = num_cc_kji + 1
!
!--                Store properties.
                   surf_cct%n_s(m,:) = surf_data(k,j,i)%normal(:)

                   surf_cct%num_edges(m)         = surf_data(k,j,i)%num_edges
                   surf_cct%faces(:,m)           = surf_data(k,j,i)%faces(:)
                   surf_cct%wall_location_c(:,m) = surf_data(k,j,i)%center(:)
                   surf_cct%area(m)              = surf_data(k,j,i)%area

                   surf_cct%i(m) = i
                   surf_cct%j(m) = j
                   surf_cct%k(m) = k
!
!--                Check if the cut-cell surface is accessed by a nearby grid point, i.e. if it is
!--                active or not.
                   IF ( surf_data(k,j,i)%face_accessed )  THEN
                      surf_cct%face_accessed(m) = .TRUE.
                      surf_cct%num_access(m) = surf_data(k,j,i)%num_access

                      surf_cct%index_access(m,:) = surf_data(k,j,i)%index_access(:)
                      IF ( surf_cct%num_access(m) == 0 )  surf_cct%face_accessed(m) = .FALSE.
                   ELSE
                      surf_cct%face_accessed(m) = .FALSE.
                   ENDIF

                   m = m + 1
                ENDIF
             ENDDO
!
!--          Store start- and end-indices at (j,i)-column.
             surf_cct%start_index(j,i) = start_index_cc
             surf_cct%end_index(j,i)   = surf_cct%start_index(j,i) + num_cc_kji - 1
             start_index_cc            = surf_cct%end_index(j,i) + 1

          ENDDO
       ENDDO
!
!--    Build a connection between the active LSM and USM surfaces and the non-active cut-cell
!--    surfacesm, where at the moment no energy balance is solved but the surface data is taken
!--    from the closest nearby LSM/USM surface.
!--    In order to distinguish between LSM and USM surfaces, store different sign: LSM (+), USM (-).
       DO  m = 1, surf_cct%ns

          IF ( .NOT. surf_cct%face_accessed(m) )  THEN

             dist_min = HUGE( 1.0_wp )

             DO  mm = 1, surf_lsm%ns
                dist = SQRT( ( surf_cct%wall_location_c(1,m) - surf_lsm%wall_location_c(1,mm) )**2 &
                           + ( surf_cct%wall_location_c(2,m) - surf_lsm%wall_location_c(2,mm) )**2 &
                           + ( surf_cct%wall_location_c(3,m) - surf_lsm%wall_location_c(3,mm) )**2 &
                           )
                IF ( dist <= dist_min )  surf_cct%m_index_ref(m) = mm
             ENDDO

             DO  mm = 1, surf_usm%ns
                dist = SQRT( ( surf_cct%wall_location_c(1,m) - surf_usm%wall_location_c(1,mm) )**2 &
                           + ( surf_cct%wall_location_c(2,m) - surf_usm%wall_location_c(2,mm) )**2 &
                           + ( surf_cct%wall_location_c(3,m) - surf_usm%wall_location_c(3,mm) )**2 &
                           )
                IF ( dist <= dist_min )  surf_cct%m_index_ref(m) = -mm
             ENDDO

          ENDIF
       ENDDO
!
!--    Now, build a connection between the the cut-cell surface considered in the RTM and the LSM
!--    and USM surfaces. As a cut-cell surface in the RTM can be composed of multiple LSM/USM
!--    surfaces (the surface is always the same, but the atmosphere information for each surf_lsm/
!--    surf_usm surface might be different), create a list of all active surfaces that belong to it.
!--    This will be used to aggregate the respective surface information for the RTM and later on,
!--    distribute the resulting radiation fluxes onto the corresponding LSM and USM surfaces.
       DO  m = 1, surf_cct%ns

          IF ( surf_cct%face_accessed(m) )  THEN
!
!--          First, count the number of connections between the cut-cell surface and the LSM and USM
!--          surfaces.
             n_accessed_lsm = 0
             n_accessed_usm = 0
             DO  mm = 1, surf_cct%num_access(m)
                IF ( surf_cct%index_access(m,mm) > 0 )  THEN
                   n_accessed_lsm = n_accessed_lsm + 1
                ELSEIF ( surf_cct%index_access(m,mm) < 0 )  THEN
                   n_accessed_usm = n_accessed_usm + 1
                ENDIF
             ENDDO
!
!--          Allocate memory for the conenctions and initialize the arrays.
             surf_cct%accessed_lsm(m)%n_accessed = n_accessed_lsm
             surf_cct%accessed_usm(m)%n_accessed = n_accessed_usm

             ALLOCATE( surf_cct%accessed_lsm(m)%m_list(1:surf_cct%accessed_lsm(m)%n_accessed) )
             ALLOCATE( surf_cct%accessed_usm(m)%m_list(1:surf_cct%accessed_usm(m)%n_accessed) )
!
!--          Generate list of
             n_accessed_lsm = 0
             n_accessed_usm = 0
             DO  mm = 1, surf_cct%num_access(m)
                IF ( surf_cct%index_access(m,mm) > 0 )  THEN
                   n_accessed_lsm = n_accessed_lsm + 1
                   surf_cct%accessed_lsm(m)%m_list(n_accessed_lsm) = surf_cct%index_access(m,mm)
                ELSEIF ( surf_cct%index_access(m,mm) < 0 )  THEN
                   n_accessed_usm = n_accessed_usm + 1
                   surf_cct%accessed_usm(m)%m_list(n_accessed_usm) = -surf_cct%index_access(m,mm)
                ENDIF
             ENDDO
          ENDIF

       ENDDO

    ENDIF
!
!-- Check for consistent type classification.
    IF ( land_surface )  THEN
       check_mismatch = .FALSE.
       DO  m = 1, surf_lsm%ns
          IF ( surf_lsm%ptype(m) <= 0  .AND.  surf_lsm%vtype(m) <= 0  .AND.                        &
               surf_lsm%wtype(m) <= 0 )  check_mismatch = .TRUE.
       ENDDO

#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, check_mismatch, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
       IF ( check_mismatch )  THEN
          message_string = 'no land-surface type defined at cut-cell surface'
          CALL message( 'cct_init', 'PAC0361', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( urban_surface )  THEN
       check_mismatch = .FALSE.
       DO  m = 1, surf_usm%ns
          IF ( surf_usm%btype(m) <= 0 )  check_mismatch = .TRUE.
       ENDDO

#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, check_mismatch, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr)
#endif
       IF ( check_mismatch )  THEN
          message_string = 'no building type defined at cut-cell surface'
          CALL message( 'cct_init', 'PAC0361', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    DEALLOCATE( surf_data )


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Check if the already defined surf_def, surf_lsm and surf_usm surfaces (defined based on
!> topography flags) correspond to a cut-cell surface. This case, re-compute the orientation of the
!> surface and its distance.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE compute_interception( p_0x, p_0y, p_0z, grid, eb_active, type_accessed )

    CHARACTER(LEN=*), INTENT(IN) :: grid                    !< type of staggered grid
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: type_accessed !< type of surface that accesses

    INTEGER(iwp) ::  grid_bit  !< bit position of corresponding grid in topo_flags
    INTEGER(iwp) ::  ia        !< temporary variable to hold value of num_access
    INTEGER(iwp) ::  i_off     !< grid index in x-direction plus offset
    INTEGER(iwp) ::  j_off     !< grid index in y-direction plus offset
    INTEGER(iwp) ::  k_off     !< grid index in z-direction plus offset
    INTEGER(iwp) ::  m         !< running index over surfaces

    LOGICAL ::  eb_active !< flag indicating energy-balance active surfaces (LSM and USM)
    LOGICAL ::  found     !< flag to indicate if the cut-cell surface has been already accessed

    REAL(wp) ::  dd                !< linear equation coefficient
    REAL(wp) ::  dist              !< distance between p_0 and interception point p_i
    REAL(wp) ::  dot_norm          !< dot product between grid line and face normal vector
    REAL(wp) ::  fac_correct = 5.0 !< correction factor
    REAL(wp) ::  p_0x              !< x-coordinate of the prognostic point
    REAL(wp) ::  p_0y              !< y-coordinate of the prognostic point
    REAL(wp) ::  p_0z              !< z-coordinate of the prognostic point

    REAL(wp), DIMENSION(3) ::  dist_max1 !< maximum possible distance between prognostic grid point
                                         !< and surface when surface is located in the same grid cell
    REAL(wp), DIMENSION(3) ::  dist_max2 !< maximum possible distance between prognostic grid point
                                         !< and surface when surface is located in the adjacent grid cell
    REAL(wp), DIMENSION(3) ::  norm_line !< direction vector describing the grid line along which a surface is accessed in the diffusion
    REAL(wp), DIMENSION(3) ::  p_0       !< coordinate vector of the prognostic point
    REAL(wp), DIMENSION(3) ::  p_i       !< interception point between access line and cut-cell face

!
!-- Depending on the location of the staggered grid, determine the maximum allowed
!-- distances between the prognostic grid point and the surface. dist_max1 defines the distance
!-- between the grid point and the surface if the surface lies within the grid box, while
!-- dist_max2 defines the distance when the cut-cell surface lies in the adjacent grid box.
    SELECT CASE ( TRIM( grid ) )

       CASE ( 's' )
          dist_max1 = (/ 0.5_wp * dx, 0.5_wp * dy, 0.5_wp * dz(1) /)
          dist_max2 = (/          dx,          dy,          dz(1) /)
          grid_bit  = 0
       CASE ( 'u' )
          dist_max1 = (/ dx, 0.5_wp * dy, 0.5_wp * dz(1) /)
          dist_max2 = (/ dx,          dy,          dz(1) /)
          grid_bit  = 1
       CASE ( 'v' )
          dist_max1 = (/ 0.5_wp * dx, dy, 0.5_wp * dz(1) /)
          dist_max2 = (/          dx, dy,          dz(1) /)
          grid_bit  = 2
       CASE ( 'w' )
          dist_max1 = (/ 0.5_wp * dx, 0.5_wp * dy, dz(1) /)
          dist_max2 = (/          dx,          dy, dz(1) /)
          grid_bit  = 3

    END SELECT
!
!-- Define coordinate vector of prognostic grid point.
    p_0 = (/ p_0z, p_0y, p_0x /)

    DO  m = surf%start_index(j,i), surf%end_index(j,i)

       IF ( .NOT. BTEST( topo_flags(k,j,i), grid_bit ) )  CYCLE

       i_off = surf%i(m) + surf%ioff(m)
       j_off = surf%j(m) + surf%joff(m)
       k_off = surf%k(m) + surf%koff(m)
!
!--    The PALM surface element and the given cut-cell surface must refer to the same or the
!--    adjacent grid cell in the direction of access. Also, there is only one offset possible,
!--    either along x-, y- or z-direction.
       IF ( TRIM( grid ) == 's' )  THEN
          IF ( ABS( k - k_off ) > 1  .OR.  ABS( j - j_off ) > 1  .OR.  ABS( i - i_off ) > 1 )      &
             CYCLE
          IF ( ABS( k - k_off ) + ABS( j - j_off ) + ABS( i - i_off ) > 1 )                        &
             CYCLE
       ENDIF
!
!--    Check if the actual grid box or the adjacent grid box contain any cut-cell surfaces.
       IF ( surf_data(k,j,i)%contain_faces  .OR.                                                   &
            surf_data(k_off,j_off,i_off)%contain_faces )  THEN
!
!--       Check if the line, defined by the vector (surf%ioff, surf%ioff, surf%joff)
!--       intercepts with the cut-cell face.
!--       First, check the cut-cell face in the grid box (i,j,k).
          norm_line = REAL( (/ surf%koff(m), surf%joff(m), surf%ioff(m) /), KIND = wp )

          dd = 1E+8_wp
          found = .FALSE.
          IF ( surf_data(k,j,i)%contain_faces )  THEN
             dot_norm = DOT_PRODUCT( norm_line, surf_data(k,j,i)%normal )

             IF ( dot_norm /= 0 )  THEN
                dd = DOT_PRODUCT( surf_data(k,j,i)%center - p_0, surf_data(k,j,i)%normal )         &
                   / dot_norm
             ENDIF
!
!--          Compute interception point.
             p_i = p_0 + norm_line * dd
!
!--          If the interception point lies between grid point (k,j,i) and the grid cell face,
!--          this is the right surface and we can store the normal vector and the distance to
!--          the cut-cell face. This is the case when the absolute value of
!--          p_0 - p_i (i.e. the distance) is smaller/equal than half a grid spacing in that
!--          direction. Moreover, store information that the face is accessed by a surface. This,
!--          however, only need to be done for the surface types that are involved for the
!--          surface-energy balance.
             dist = SQRT( SUM( ( p_0 - p_i )**2 ) )

             IF ( dist <= SQRT( SUM( ( norm_line * dist_max1 )**2 ) )  .AND.  dist > 0.0_wp )  THEN

                surf%n_s(m,:) = surf_data(k,j,i)%normal(:)
                surf%z_mo(m) = dist
!
!--             For energy-balance surfaces also store the corresponding type, later on used to
!--             initialize the surface properties in the LSM/USM accordingly.
                IF ( PRESENT( type_accessed ) )  THEN
                   IF ( type_accessed == 'lsm' )  THEN
                      surf%ptype(m) = surf_data(k,j,i)%ptype
                      surf%vtype(m) = surf_data(k,j,i)%vtype
                      surf%wtype(m) = surf_data(k,j,i)%wtype
                   ENDIF
                   IF ( type_accessed == 'usm' )  THEN
                      surf%btype(m) = surf_data(k,j,i)%btype
                      surf%bid(m)   = surf_data(k,j,i)%bid
                   ENDIF
                ENDIF
!
!--             Cut-cell surface found.
                found = .TRUE.
!
!--             Store data particularly required for radiation-effective surfaces.
                IF ( eb_active )  THEN
                   surf_data(k,j,i)%face_accessed = .TRUE.

                   IF ( .NOT. ANY( surf_data(k,j,i)%index_access == m ) )  THEN
                      surf_data(k,j,i)%num_access = surf_data(k,j,i)%num_access + 1
                      surf_data(k,j,i)%index_access(surf_data(k,j,i)%num_access) = m

                      IF ( type_accessed == 'lsm' )  THEN
                         surf_data(k,j,i)%index_access(surf_data(k,j,i)%num_access) = m
                      ELSE
                         surf_data(k,j,i)%index_access(surf_data(k,j,i)%num_access) = -m
                      ENDIF
                   ENDIF

                   IF ( surf_data(k,j,i)%wall_surface )  surf%cut_cell_wall(m) = .TRUE.
                   IF ( surf_data(k,j,i)%roof_surface )  surf%cut_cell_roof(m) = .TRUE.
!
!--                Store number of vertices, edges and the wall location.
                   surf%num_edges(m)         = surf_data(k,j,i)%num_edges
                   surf%faces(:,m)           = surf_data(k,j,i)%faces(:)
                   surf%wall_location_c(:,m) = surf_data(k,j,i)%center(:)
                   surf%area(m)              = surf_data(k,j,i)%area
!
!--                Store indices that link the surface to a 3D field.
                   surf%i_cc(m) = i
                   surf%j_cc(m) = j
                   surf%k_cc(m) = k

                   surf%accessed(m) = .TRUE.
                ENDIF
             ENDIF

          ENDIF
!
!--       Repeat this action for the possible cut-cell surface in the adjacent grid box, if no
!--       valid interception point has been found yet. Please see comments before.
          IF ( surf_data(k_off,j_off,i_off)%contain_faces  .AND.  .NOT. found )  THEN

             dot_norm = DOT_PRODUCT( norm_line, surf_data(k_off,j_off,i_off)%normal )

             IF ( dot_norm /= 0 )  THEN
                dd = DOT_PRODUCT( surf_data(k_off,j_off,i_off)%center - p_0,                       &
                                  surf_data(k_off,j_off,i_off)%normal ) / dot_norm
             ENDIF

             p_i = p_0 + norm_line * dd

             dist = MIN( SQRT( SUM( ( p_0 - p_i )**2 ) ),                                          &
                         SQRT( SUM( ( p_0 - surf_data(k_off,j_off,i_off)%center )**2 ) ) )
!
!--          Here, consider a correction factor. In rare cases it can happen that the plane
!--          approximation by the mass center and the average normal vector yield to unrealistic
!--          distances. In order to assure that all surfaces are initialized, treat those surfaces
!--          nevertheless and limit the distance between the surface and the prognostic grid point.
             IF ( dist <= fac_correct * SQRT( SUM( ( norm_line * dist_max2 )**2 ) )  .AND.         &
                  dist > 0.0_wp )  THEN

                surf%n_s(m,:) = surf_data(k_off,j_off,i_off)%normal(:)
                surf%z_mo(m) = MIN( dist, SQRT( SUM( ( norm_line * dist_max2 )**2 ) ) )

                IF ( PRESENT( type_accessed ) )  THEN
                   IF ( type_accessed == 'lsm' )  THEN
                      surf%ptype(m) = surf_data(k_off,j_off,i_off)%ptype
                      surf%vtype(m) = surf_data(k_off,j_off,i_off)%vtype
                      surf%wtype(m) = surf_data(k_off,j_off,i_off)%wtype
                   ENDIF
                   IF ( type_accessed == 'usm' )  THEN
                      surf%btype(m) = surf_data(k_off,j_off,i_off)%btype
                      surf%bid(m)   = surf_data(k_off,j_off,i_off)%bid
                   ENDIF
                ENDIF

                found = .TRUE.

                IF ( eb_active )  THEN
                   surf_data(k_off,j_off,i_off)%face_accessed = .TRUE.

                   IF ( .NOT. ANY( surf_data(k_off,j_off,i_off)%index_access == m ) )  THEN
                      ia = surf_data(k_off,j_off,i_off)%num_access + 1
                      surf_data(k_off,j_off,i_off)%num_access = ia

                      IF ( type_accessed == 'lsm' )  THEN
                         surf_data(k_off,j_off,i_off)%index_access(ia) = m
                      ELSE
                         surf_data(k_off,j_off,i_off)%index_access(ia) = -m
                      ENDIF
                   ENDIF

                   IF ( surf_data(k_off,j_off,i_off)%wall_surface )  surf%cut_cell_wall(m) = .TRUE.
                   IF ( surf_data(k_off,j_off,i_off)%roof_surface )  surf%cut_cell_roof(m) = .TRUE.

                   surf%num_edges(m)         = surf_data(k_off,j_off,i_off)%num_edges
                   surf%faces(:,m)           = surf_data(k_off,j_off,i_off)%faces(:)
                   surf%wall_location_c(:,m) = surf_data(k_off,j_off,i_off)%center(:)
                   surf%area(m)              = surf_data(k_off,j_off,i_off)%area
!
!--                Store indices that link the surface to a 3D field.
                   surf%i_cc(m) = i_off
                   surf%j_cc(m) = j_off
                   surf%k_cc(m) = k_off

                   surf%accessed(m) = .TRUE.
                ENDIF
             ENDIF

          ENDIF

       ENDIF

    ENDDO

 END SUBROUTINE compute_interception


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Check and modify general surface settings.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE cct_set_surface_properties

    INTEGER(iwp) ::  iref  !< reference index in x-direction
    INTEGER(iwp) ::  jref  !< reference index in x-direction
    INTEGER(iwp) ::  kref  !< reference index in x-direction
    INTEGER(iwp) ::  m     !< running index over surfaces

    REAL(wp) ::  dxyz  !< distance from grid-cell center to face


!
!-- Normalize normal vector (if not done already in the pre-processing) and set effective normal
!-- vector component.
    DO  m = surf%start_index(j,i), surf%end_index(j,i)

       CALL normalize_vector( surf%n_s(m,:) )

       IF ( surf%upward(m)     .OR.  surf%downward(m)  )  surf%n_eff(m) = surf%n_s(m,1)
       IF ( surf%northward(m)  .OR.  surf%southward(m) )  surf%n_eff(m) = surf%n_s(m,2)
       IF ( surf%eastward(m)   .OR.  surf%westward(m)  )  surf%n_eff(m) = surf%n_s(m,3)

    ENDDO
!
!-- If the distance to the wall is too close, Monin-Obukhov relations can become critical.
!-- To prevent this, change the atmosphere reference grid point so that the distance between the
!-- reference grid point in the atmosphere and the surface is sufficiently large.
    DO  m = surf%start_index(j,i), surf%end_index(j,i)

       IF ( surf%upward(m)     .OR.  surf%downward(m)  )  dxyz = 0.5_wp * dzu(surf%k(m))
       IF ( surf%northward(m)  .OR.  surf%southward(m) )  dxyz = 0.5_wp * dy
       IF ( surf%eastward(m)   .OR.  surf%westward(m)  )  dxyz = 0.5_wp * dx

       IF ( surf%z_mo(m) < dxyz )  THEN
!
!--       Compute new reference grid point.
          iref = surf%i(m) - surf%ioff(m)
          jref = surf%j(m) - surf%joff(m)
          kref = surf%k(m) - surf%koff(m)
!
!--       Check if this grid point lies within the atmosphere. If this is the case, update the
!--       reference indices and ajdust z_mo accordingly, which is the current distance + one
!--       grid spacing. If no other reference grid point can be used, artificially increase
!--       the surface-layer height. Since this can happen only in narrow canyons that are
!--       poorly resolved anywhere, the expected error is small.
          IF ( BTEST( topo_flags(kref,jref,iref), 0 ) )  THEN
             surf%iref(m) = iref
             surf%jref(m) = jref
             surf%kref(m) = kref

             surf%z_mo(m) = surf%z_mo(m) + 2.0_wp * dxyz
          ELSE
             surf%z_mo(m) = dxyz
          ENDIF
       ENDIF

    ENDDO
!
!-- If the distance to the wall is still too close (e.g. because it was not possible to increase
!-- it), set flag tke_production to false, in order to disable TKE production at near surface grid
!-- points. This is necessary to avoid extremely large diffusion coefficients (only happens in
!-- rare cases of narrow canyons).
    IF ( ALLOCATED( surf%tke_production ) )  THEN

       DO  m = surf%start_index(j,i), surf%end_index(j,i)
          IF ( surf%upward(m)     .OR.  surf%downward(m)  )  dxyz = 0.5_wp * dzu(surf%k(m))
          IF ( surf%northward(m)  .OR.  surf%southward(m) )  dxyz = 0.5_wp * dy
          IF ( surf%eastward(m)   .OR.  surf%westward(m)  )  dxyz = 0.5_wp * dx

          IF ( surf%z_mo(m) <= 0.5_wp * dxyz )  THEN
             surf%tke_production(m) = .FALSE.
          ENDIF
       ENDDO

    ENDIF
!
!-- Set flag at surfaces to indicate whether stability needs to be considered or not.
!-- In the step-like topography, stability is only considered at upward-facing surfaces.
!-- In the cut-cell topography, however, even east/west/nort/southward-facing flagged surfaces can
!-- be approximately "flat", i.e. they face upward. In order to evaluate, if stability is
!-- considered or not, compute the slope of the surface. Stability is considered at slope
!-- angles <= 30 degrees. (The exact angle needs to be tested somehow!)
    IF ( ALLOCATED( surf%consider_stability ) )  THEN
       DO  m = surf%start_index(j,i), surf%end_index(j,i)
          surf%consider_stability(m) =                                                             &
                   MERGE( .TRUE., .FALSE., ACOS( surf%n_s(m,1) ) * radiants_to_degrees <= 30.0_wp )
       ENDDO
    ENDIF

 END SUBROUTINE cct_set_surface_properties

 END SUBROUTINE cct_init


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Input cut-cell data from static input file. Already pre-processed polygon data is read.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE cct_input

    INTEGER(iwp) ::  n          !< running index over the number of total cut-cell surfaces in domain
    INTEGER(iwp) ::  num_faces  !< counter and running variable for the cut-cell surfaces

!
!-- Check for input from static file.
    IF ( .NOT. input_pids_static )  THEN
       message_string = 'cut_cell_topopgraphy = .T. but static driver file is missing'
       CALL message( 'cct_input', 'PAC0364', 1, 2, 0, 6, 0 )
    ENDIF

#if defined( __netcdf )
!
!-- Open file and inquire a list of variables.
    CALL open_read_file( TRIM( input_file_static ) // TRIM( coupling_char ), pids_id )

    CALL inquire_num_variables( pids_id, num_var_pids )

    ALLOCATE( vars_pids(1:num_var_pids) )
    CALL inquire_variable_names( pids_id, vars_pids )
!
!-- Read dimensions.
    CALL get_dimension_length( pids_id, cct_global_f%dim_3d, 'dim_3d' )
    CALL get_dimension_length( pids_id, cct_global_f%num_faces, 'cct_num_faces' )
    CALL get_dimension_length( pids_id, cct_global_f%num_faces_vert,                               &
                               'cct_max_num_vertices_per_face' )
    CALL get_dimension_length( pids_id, cct_global_f%num_vert, 'cct_num_vert' )
    CALL get_dimension_length( pids_id, cct_global_f%dim_vertex_coords, 'cct_dim_vertex_coords' )
    CALL get_dimension_length( pids_id, cct_global_f%dim_vertex_shifts, 'cct_dim_vertex_shifts' )
    CALL get_attribute( pids_id, 'empty_vert', cct_global_f%fill_value_int, .TRUE. )
!
!-- Read actual vertex data.
    IF ( check_existence( vars_pids, 'cct_vertices' ) )  THEN
       ALLOCATE( cct_global_f%vertices(1:cct_global_f%dim_3d,1:cct_global_f%num_vert ) )
       CALL get_variable( pids_id, 'cct_vertices', cct_global_f%vertices, cct_global_f%dim_3d,     &
                          cct_global_f%num_vert )
    ELSE
       message_string = 'variable "cct_vertices" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read adjusted vertex data for radiation.
    IF ( check_existence( vars_pids, 'cct_vertex_coords' ) )  THEN
       ALLOCATE( cct_global_f%vertex_coords(1:cct_global_f%dim_vertex_coords,                      &
                                            1:cct_global_f%num_vert ) )
       CALL get_variable( pids_id, 'cct_vertex_coords', cct_global_f%vertex_coords,                &
                          cct_global_f%dim_vertex_coords, cct_global_f%num_vert )
    ELSE
       message_string = 'variable "cct_vertex_coords" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read adjusted vertex data for radiation.
    IF ( check_existence( vars_pids, 'cct_vertex_shifts' ) )  THEN
       ALLOCATE( cct_global_f%vertex_shifts(1:cct_global_f%num_vert ) )
       CALL get_variable( pids_id, 'cct_vertex_shifts', cct_global_f%vertex_shifts,                &
                          1, cct_global_f%num_vert )
    ELSE
       message_string = 'variable "cct_vertex_shifts" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read polygon data.
    IF ( check_existence( vars_pids, 'cct_vertices_per_face' ) )  THEN
       ALLOCATE( cct_global_f%faces(1:cct_global_f%num_faces_vert,1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_vertices_per_face', cct_global_f%faces,                    &
                          cct_global_f%num_faces_vert, cct_global_f%num_faces )
    ELSE
       message_string = 'variable "cct_vertices_per_face" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read number of edges of face.
    IF ( check_existence( vars_pids, 'cct_num_vertices_per_face' ) )  THEN
       ALLOCATE( cct_global_f%num_edges(1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_num_vertices_per_face', cct_global_f%num_edges )
    ELSE
       message_string = 'variable "cct_num_vertices_per_face" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read normal vector.
    IF ( check_existence( vars_pids, 'cct_face_normal_vector' ) )  THEN
       ALLOCATE( cct_global_f%normals(1:cct_global_f%dim_3d,1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_face_normal_vector', cct_global_f%normals,                 &
                          cct_global_f%dim_3d, cct_global_f%num_faces )
    ELSE
       message_string = 'variable "cct_face_normal_vector" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read surface area.
    IF ( check_existence( vars_pids, 'cct_face_area' ) )  THEN
       ALLOCATE( cct_global_f%area(1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_face_area', cct_global_f%area,                             &
                          1, cct_global_f%num_faces )
    ELSE
       message_string = 'variable "cct_face_area" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read mass center.
    IF ( check_existence( vars_pids, 'cct_face_center' ) )  THEN
       ALLOCATE( cct_global_f%centers(1:cct_global_f%dim_3d,1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_face_center', cct_global_f%centers,                        &
                          cct_global_f%dim_3d, cct_global_f%num_faces )
    ELSE
       message_string = 'variable "cct_face_center" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read index tuples.
    IF ( check_existence( vars_pids, 'cct_3d_grid_indices' ) )  THEN
       ALLOCATE( cct_global_f%kji(1:cct_global_f%dim_3d,1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_3d_grid_indices', cct_global_f%kji,                        &
                          cct_global_f%dim_3d, cct_global_f%num_faces )
    ELSE
       message_string = 'variable "cct_3d_grid_indices" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read general type-classification.
    IF ( check_existence( vars_pids, 'cct_surface_type_classification' ) )  THEN
       ALLOCATE( cct_global_f%types(1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_surface_type_classification',cct_global_f%types )
    ELSE
       message_string = 'variable "cct_surface_type_classification" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Read vegetation-, pavement-, and water-types associated with cut-cell surfaces.
    IF ( land_surface )  THEN
       IF ( check_existence( vars_pids, 'cct_vegetation_type_classification' ) )  THEN
          ALLOCATE( cct_global_f%veg_types(1:cct_global_f%num_faces ) )
          CALL get_variable( pids_id, 'cct_vegetation_type_classification', cct_global_f%veg_types )
       ELSE
          message_string = 'variable "cct_vegetation_type_classification" missing in static driver'
          CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( check_existence( vars_pids, 'cct_pavement_type_classification' ) )  THEN
          ALLOCATE( cct_global_f%pav_types(1:cct_global_f%num_faces ) )
          CALL get_variable( pids_id, 'cct_pavement_type_classification', cct_global_f%pav_types )
       ELSE
          message_string = 'variable "cct_pavement_type_classification" missing in static driver'
          CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( check_existence( vars_pids, 'cct_water_type_classification' ) )  THEN
          ALLOCATE( cct_global_f%wat_types(1:cct_global_f%num_faces ) )
          CALL get_variable( pids_id, 'cct_water_type_classification', cct_global_f%wat_types )
       ELSE
          message_string = 'variable "cct_water_type_classification" missing in static driver'
          CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Read building-type associated with cut-cell surfaces.
    IF ( urban_surface )  THEN
       IF ( check_existence( vars_pids, 'cct_building_type_classification' ) )  THEN
          ALLOCATE( cct_global_f%build_types(1:cct_global_f%num_faces ) )
          CALL get_variable( pids_id, 'cct_building_type_classification', cct_global_f%build_types )
       ELSE
          message_string = 'variable "cct_building_type_classification" missing in static driver'
          CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Read building ID associated with cut-cell surfaces.
    IF ( check_existence( vars_pids, 'cct_building_id_classification' ) )  THEN
       ALLOCATE( cct_global_f%build_id(1:cct_global_f%num_faces ) )
       CALL get_variable( pids_id, 'cct_building_id_classification', cct_global_f%build_id )
    ELSE
       message_string = 'variable "cct_building_id_classification" missing in static driver'
       CALL message( 'cct_input', 'PAC0360', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Close topography input file and clean-up memory.
    CALL close_input_file( pids_id )
    DEALLOCATE( vars_pids )
#endif

!
!-- Reduce the array size for faster topography processing. First, count the number of faces on
!-- local subdomain (including ghost layers).
    num_faces = 0
    DO  n = 1, cct_global_f%num_faces
       IF ( nysg <= cct_global_f%kji(2,n)  .AND.  cct_global_f%kji(2,n) <= nyng  .AND.             &
            nxlg <= cct_global_f%kji(3,n)  .AND.  cct_global_f%kji(3,n) <= nxrg )                  &
       THEN
          num_faces = num_faces + 1
       ENDIF
    ENDDO
!
!-- Allocate arrays for cropped data structure.
    cct%num_faces         = num_faces
    cct%dim_3d            = cct_global_f%dim_3d
    cct%num_faces_vert    = cct_global_f%num_faces_vert
    cct%num_vert          = cct_global_f%num_vert
    cct%dim_vertex_coords = cct_global_f%dim_vertex_coords

    ALLOCATE( cct%area(1:cct%num_faces ) )
    ALLOCATE( cct%types(1:cct%num_faces ) )
    ALLOCATE( cct%num_edges(1:cct%num_faces ) )
    ALLOCATE( cct%build_id(1:cct%num_faces ) )
    IF ( land_surface )  THEN
       ALLOCATE( cct%veg_types(1:cct%num_faces ) )
       ALLOCATE( cct%pav_types(1:cct%num_faces ) )
       ALLOCATE( cct%wat_types(1:cct%num_faces ) )
    ENDIF
    IF ( urban_surface )  ALLOCATE( cct%build_types(1:cct%num_faces ) )

    ALLOCATE( cct%centers(1:cct%dim_3d,1:cct%num_faces ) )
    ALLOCATE( cct%kji(1:cct%dim_3d,1:cct%num_faces ) )
    ALLOCATE( cct%normals(1:cct%dim_3d,1:cct%num_faces ) )

    ALLOCATE( cct%faces(1:cct%num_faces_vert,1:cct%num_faces ) )

    ALLOCATE( cct%vertices(1:cct%dim_3d,1:cct%num_vert ) )
    ALLOCATE( cct%vertex_coords(1:cct%dim_vertex_coords,1:cct%num_vert ) )
    ALLOCATE( cct%vertex_shifts(1:cct%num_vert ) )
!
!-- Arrays which do not depend on the number of faces (vertex data), is simply copied.
    cct%vertices      = cct_global_f%vertices
    cct%vertex_coords = cct_global_f%vertex_coords
    cct%vertex_shifts = cct_global_f%vertex_shifts
!
!-- Now, all data depending on the number of faces is filled for the local subdomain.
    num_faces = 0
    DO  n = 1, cct_global_f%num_faces
       IF ( nysg <= cct_global_f%kji(2,n)  .AND.  cct_global_f%kji(2,n) <= nyng  .AND.             &
            nxlg <= cct_global_f%kji(3,n)  .AND.  cct_global_f%kji(3,n) <= nxrg )                  &
       THEN

          num_faces = num_faces + 1

          cct%area(num_faces)      = cct_global_f%area(n)
          cct%types(num_faces)     = cct_global_f%types(n)
          cct%num_edges(num_faces) = cct_global_f%num_edges(n)
          cct%build_id(num_faces)  = cct_global_f%build_id(n)
          IF ( land_surface )  THEN
             cct%veg_types(num_faces) = cct_global_f%veg_types(n)
             cct%pav_types(num_faces) = cct_global_f%pav_types(n)
             cct%wat_types(num_faces) = cct_global_f%wat_types(n)
          ENDIF
          IF ( land_surface )  THEN
             cct%build_types(num_faces) = cct_global_f%build_types(n)
          ENDIF

          cct%centers(:,num_faces) = cct_global_f%centers(:,n)
          cct%kji(:,num_faces)     = cct_global_f%kji(:,n)
          cct%normals(:,num_faces) = cct_global_f%normals(:,n)
          cct%faces(:,num_faces)   = cct_global_f%faces(:,n)

       ENDIF
    ENDDO
!
!-- Deallocate all arrays in cct_global_f (global data).
    IF ( ALLOCATED( cct_global_f%area          ) )  DEALLOCATE( cct_global_f%area )
    IF ( ALLOCATED( cct_global_f%types         ) )  DEALLOCATE( cct_global_f%types )
    IF ( ALLOCATED( cct_global_f%num_edges     ) )  DEALLOCATE( cct_global_f%num_edges )
    IF ( ALLOCATED( cct_global_f%build_id      ) )  DEALLOCATE( cct_global_f%build_id )
    IF ( ALLOCATED( cct_global_f%veg_types     ) )  DEALLOCATE( cct_global_f%veg_types )
    IF ( ALLOCATED( cct_global_f%pav_types     ) )  DEALLOCATE( cct_global_f%pav_types )
    IF ( ALLOCATED( cct_global_f%wat_types     ) )  DEALLOCATE( cct_global_f%wat_types )
    IF ( ALLOCATED( cct_global_f%build_types   ) )  DEALLOCATE( cct_global_f%build_types )
    IF ( ALLOCATED( cct_global_f%centers       ) )  DEALLOCATE( cct_global_f%centers )
    IF ( ALLOCATED( cct_global_f%kji           ) )  DEALLOCATE( cct_global_f%kji )
    IF ( ALLOCATED( cct_global_f%normals       ) )  DEALLOCATE( cct_global_f%normals )
    IF ( ALLOCATED( cct_global_f%faces         ) )  DEALLOCATE( cct_global_f%faces )
    IF ( ALLOCATED( cct_global_f%vertices      ) )  DEALLOCATE( cct_global_f%vertices )
    IF ( ALLOCATED( cct_global_f%vertex_coords ) )  DEALLOCATE( cct_global_f%vertex_coords )
    IF ( ALLOCATED( cct_global_f%vertex_shifts ) )  DEALLOCATE( cct_global_f%vertex_shifts )

 END SUBROUTINE cct_input


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Distributes data, in particular RTM-relevant arrays, onto LSM and USM surfaces.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE cct_to_surface_types

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var_tmp !< temporary 3d array used to exchange ghost-point data


!
!-- Set temporary variable to a huge number.
    var_tmp = HUGE( 1.0_wp )

    CALL exchange_cut_cell_data( surf_cct%rad_lw_dif )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_lw_dif,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_lw_dif,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_lw_in )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_lw_in,       &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_lw_in,       &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_lw_out )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_lw_out,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_lw_out,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_lw_ref )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_lw_ref,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_lw_ref,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_lw_res )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_lw_res,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_lw_res,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_net )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_net,         &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_net,         &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_sw_dif )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_sw_dif,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_sw_dif,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_sw_dir )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_sw_dir,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_sw_dir,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_sw_in )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_sw_in,       &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_sw_in,       &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_sw_out )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_sw_out,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_sw_out,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_sw_ref )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_sw_ref,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_sw_ref,      &
                            surf_usm%accessed )

    CALL exchange_cut_cell_data( surf_cct%rad_sw_res )
    IF ( surf_lsm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_lsm%i_cc, surf_lsm%j_cc, surf_lsm%k_cc, surf_lsm%rad_sw_res,      &
                            surf_lsm%accessed )
    IF ( surf_usm%ns > 0 )                                                                         &
       CALL map_on_surface( surf_usm%i_cc, surf_usm%j_cc, surf_usm%k_cc, surf_usm%rad_sw_res,      &
                            surf_usm%accessed )

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Writes cut-cell surface data onto 3D array and exchanges ghost-point data.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE exchange_cut_cell_data( var_surf_cc )

       INTEGER(iwp) ::  i  !< grid index x-direction
       INTEGER(iwp) ::  j  !< grid index y-direction
       INTEGER(iwp) ::  k  !< grid index z-direction
       INTEGER(iwp) ::  m  !< running index over cut-cell surface

       REAL(wp), DIMENSION(:) ::  var_surf_cc  !< passed variable on surface data structure for cut-cells


       DO  m = 1, surf_cct%ns
          i = surf_cct%i(m)
          j = surf_cct%j(m)
          k = surf_cct%k(m)
          var_tmp(k,j,i) = var_surf_cc(m)
       ENDDO

       CALL exchange_horiz( var_tmp, nbgp )

    END SUBROUTINE exchange_cut_cell_data


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Writes 3D data onto surfaces used in the prognostic equations of the energy-balance models.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE map_on_surface( i_cc, j_cc, k_cc, var_surf, accessed )

       INTEGER(iwp) ::  m  !< running index over surfaces

       INTEGER(iwp), DIMENSION(:) ::  i_cc  !< passed x-indices
       INTEGER(iwp), DIMENSION(:) ::  j_cc  !< passed y-indices
       INTEGER(iwp), DIMENSION(:) ::  k_cc  !< passed z-indices

       LOGICAL, DIMENSION(:) ::  accessed  !< flag indicating if surface elemement accesses a cut-cell

       REAL(wp), DIMENSION(:) ::  var_surf  !< passed variable on surface data structure


       DO  m = 1, SIZE( var_surf )
          IF ( accessed(m) )  var_surf(m) = var_tmp(k_cc(m),j_cc(m),i_cc(m))
       ENDDO

    END SUBROUTINE map_on_surface

 END SUBROUTINE cct_to_surface_types


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Updates information on non-active energy-balance surfaces. For now, this information is taken
!> from the nearest active surface. In the future, this might be realized by solving the
!> energy balance explicitly.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_types_to_cct

    INTEGER(iwp) ::  m      !< running index over cut-cell surfaces
    INTEGER(iwp) ::  m_ref  !< reference index
    INTEGER(iwp) ::  mm     !< reference index

    REAL(wp) ::  var_sum  !< dummy used to sum-up data


    DO  m = 1, surf_cct%ns

       IF ( .NOT. surf_cct%face_accessed(m) )  THEN
!
!--       Take data from nearest LSM or USM surface.
          m_ref = ABS(surf_cct%m_index_ref(m))

          IF ( surf_cct%m_index_ref(m) > 0 )  THEN
             surf_cct%pt_surface(m) = surf_lsm%pt_surface(m_ref)
             surf_cct%albedo(m)     = SUM( surf_lsm%frac(m_ref,:) * surf_lsm%albedo(m_ref,:) )
             surf_cct%emissivity(m) = SUM( surf_lsm%frac(m_ref,:) * surf_lsm%emissivity(m_ref,:) )
          ELSE
             surf_cct%pt_surface(m) = surf_usm%pt_surface(m_ref)
             surf_cct%albedo(m)     = SUM( surf_usm%frac(m_ref,:) * surf_usm%albedo(m_ref,:) )
             surf_cct%emissivity(m) = SUM( surf_usm%frac(m_ref,:) * surf_usm%emissivity(m_ref,:) )
          ENDIF

       ELSE
!
!--       Aggregate surface temperature, albedo and emissivity from accessing LSM / USM surfaces.
          var_sum = 0.0_wp
          DO  mm = 1, surf_cct%accessed_lsm(m)%n_accessed
             m_ref = surf_cct%accessed_lsm(m)%m_list(mm)
             var_sum = var_sum + surf_lsm%pt_surface(m_ref)
          ENDDO
          DO  mm = 1, surf_cct%accessed_usm(m)%n_accessed
             m_ref = surf_cct%accessed_usm(m)%m_list(mm)
             var_sum = var_sum + surf_usm%pt_surface(m_ref)
          ENDDO
          surf_cct%pt_surface(m) = var_sum / ( surf_cct%accessed_lsm(m)%n_accessed         &
                                             + surf_cct%accessed_usm(m)%n_accessed )
!
!--       Aggregate surface albedo.
          var_sum = 0.0_wp
          DO  mm = 1, surf_cct%accessed_lsm(m)%n_accessed
             m_ref = surf_cct%accessed_lsm(m)%m_list(mm)
             var_sum = var_sum + SUM( surf_lsm%frac(m_ref,:) * surf_lsm%albedo(m_ref,:) )
          ENDDO
          DO  mm = 1, surf_cct%accessed_usm(m)%n_accessed
             m_ref = surf_cct%accessed_usm(m)%m_list(mm)
             var_sum = var_sum + SUM( surf_usm%frac(m_ref,:) * surf_usm%albedo(m_ref,:) )
          ENDDO
          surf_cct%albedo(m) = var_sum / ( surf_cct%accessed_lsm(m)%n_accessed             &
                                         + surf_cct%accessed_usm(m)%n_accessed )
!
!--       Aggregate emissivity.
          var_sum = 0.0_wp
          DO  mm = 1, surf_cct%accessed_lsm(m)%n_accessed
             m_ref = surf_cct%accessed_lsm(m)%m_list(mm)
             var_sum = var_sum + SUM( surf_lsm%frac(m_ref,:) * surf_lsm%emissivity(m_ref,:) )
          ENDDO
          DO  mm = 1, surf_cct%accessed_usm(m)%n_accessed
             m_ref = surf_cct%accessed_usm(m)%m_list(mm)
             var_sum = var_sum + SUM( surf_usm%frac(m_ref,:) * surf_usm%emissivity(m_ref,:) )
          ENDDO
          surf_cct%emissivity(m) = var_sum / ( surf_cct%accessed_lsm(m)%n_accessed         &
                                             + surf_cct%accessed_usm(m)%n_accessed )
       ENDIF

    ENDDO

 END SUBROUTINE surface_types_to_cct


 END MODULE cut_cell_topography_mod
