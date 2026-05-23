!> @file surface_mod.f90
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
!> Surface module defines derived data structures to treat surface-adjacent grid cells. Three
!> different types of surfaces are defined: default surfaces, natural surfaces, and urban surfaces.
!> @todo Clean up urban-surface variables (some of them are not used any more)
!--------------------------------------------------------------------------------------------------!
 MODULE surface_mod

    USE kinds

    IMPLICIT NONE

!
!-- Data type used to identify grid-points where horizontal boundary conditions are applied
    TYPE bc_type

       INTEGER(iwp) ::  ns     !< number of wall-adjacent grid points (on s-grid) in the subdomain
       INTEGER(iwp) ::  ns_tot !< number of wall-adjacent grid points (on s-grid) in the total domain

       INTEGER(iwp) ::  ns_bgp     !< number of boundary grid points (on s-grid) in the subdomain
       INTEGER(iwp) ::  ns_bgp_tot !< number of boundary grid points (on s-grid) in the total domain

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i      !< x-index of wall-adjacent grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_bgp  !< x-index of boundary grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ioff   !< offset value in x indicating the position
                                                          !< of the surface with respect to the reference grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j      !< y-index of wall-adjacent grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_bgp  !< y-index of boundary grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  joff   !< offset value in y indicating the position
                                                          !< of the surface with respect to the reference grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k      !< z-index of wall-adjacent grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_bgp  !< z-index of boundary grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  koff   !< offset value in z indicating the position
                                                          !< of the surface with respect to the reference grid point

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  end_index    !< end index within surface data type for given (j,i)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  start_index  !< start index within surface data type for given (j,i)

    END TYPE bc_type
!
!-- Data structure which gathers information from all surface elements of all types on subdomain for
!-- output in surface_data_output_mod
    TYPE surf_out_type

       INTEGER(iwp) ::  ns             !< number of surface elements on subdomain
       INTEGER(iwp) ::  ns_total       !< total number of surface elements
       INTEGER(iwp) ::  npoints        !< number of points / vertices which define a surface element (on subdomain)
       INTEGER(iwp) ::  npoints_total  !< total number of points / vertices which define a surface element
       INTEGER(iwp) ::  num_vert       !< local number of vertices
       INTEGER(iwp) ::  num_faces_vert !< number of vertices that define polygons


       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  s  !< coordinate for NetCDF output, number of the surface element

       REAL(wp) ::  fillvalue = -9999.0_wp  !< fillvalue for surface elements which are not defined

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  area      !< surface element area for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  azimuth   !< azimuth orientation coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  es_utm    !< E-UTM coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  ns_utm    !< E-UTM coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  xs        !< x-coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  ys        !< y-coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  zs        !< z-coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  zenith    !< zenith orientation coordinate for NetCDF output
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  var_out   !< output variable
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_av    !< variable used for averaging
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  points    !< points  / vertices of a surface element
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  polygons  !< polygon data of a surface element
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vertices  !< polygon data of a surface element
    END TYPE surf_out_type
!
!-- Data type used to identify and treat surface-adjacent grid points
    TYPE surf_type

       INTEGER(iwp) ::  ns             !< number of surface elements on subdomain on s-grid
       INTEGER(iwp) ::  ns_tot_up = 0  !< number of upward-facing surface elements within the entire model domain
       INTEGER(iwp) ::  ns_tot_v  = 0  !< number of vertical surface elements within the entire model domain
       INTEGER(iwp) ::  ns_u = 0       !< number of surface elements on subdomain on u-grid
       INTEGER(iwp) ::  ns_v = 0       !< number of surface elements on subdomain on v-grid
       INTEGER(iwp) ::  ns_w = 0       !< number of surface elements on subdomain on v-grid

       INTEGER(iwp) ::  nzb_soil       !< lower index of soil grid in LSM
       INTEGER(iwp) ::  nzb_wall       !< lower index of wall grid in USM
       INTEGER(iwp) ::  nzt_soil       !< upper index of soil grid in LSM
       INTEGER(iwp) ::  nzt_wall       !< upper index of wall grid in USM

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i     !< x-index linking surface to prognostic grid point in the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_cc  !< x-index linking surface to location of a cut-cell surface in 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ioff  !< offset value in x indicating the position
                                                         !< of the surface with respect to the reference grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  iref  !< x-index linking surface to atmosphere grid point which enters MOST relation (this
                                                         !< is usually equal to i index, but might be differ from it in case of cut-cell surfaces)
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j     !< y-index linking surface to prognostic grid point in the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_cc  !< x-index linking surface to location of a cut-cell surface in 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  joff  !< offset value in y indicating the position
                                                         !< of the surface with respect to the reference grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  jref  !< y-index linking surface to atmosphere grid point which enters MOST relation (this
                                                         !< is usually equal to i index, but might be differ from it in case of cut-cell surfaces)
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k     !< z-index linking surface to prognostic grid point in the PALM 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_cc  !< x-index linking surface to location of a cut-cell surface in 3D-grid
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  koff  !< offset value in z indicating the position
                                                         !< of the surface with respect to the reference grid point
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  kref  !< k-index linking surface to atmosphere grid point which enters MOST relation (this
                                                         !< is usually equal to i index, but might be differ from it in case of cut-cell surfaces)

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  start_index  !< Start index within surface data type for given (j,i)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  end_index    !< End index within surface data type for given (j,i)

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  consider_stability !< flag indicating surface where stability correction in MOST is considered
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  downward           !< flag indicating downward-facing surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  eastward           !< flag indicating eastward-facing surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  northward          !< flag indicating northward-facing surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  southward          !< flag indicating southward-facing surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  tke_production     !< flag to enable/disable near-surface TKE production
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  upward             !< flag indicating upward-facing surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  upward_top         !< flag indicating the uppermost upward-facing surface at (j,i)
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  westward           !< flag indicating westward-facing surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_mo         !< distance to surface (equals surface-layer height for MOST assumptions)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uvw_abs      !< absolute surface-parallel velocity on grid center
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uvw_abs_uv   !< absolute surface-parallel velocity on u- or v-grid
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uvw_abs_w    !< absolute surface-parallel velocity on w-grid
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us           !< friction velocity valid for grid center
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_uvgrid    !< friction velocity valid for u- or v-grid
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_wgrid     !< friction velocity valid for w-grid
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts           !< scaling parameter temerature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs           !< scaling parameter humidity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ss           !< scaling parameter passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcs          !< scaling parameter qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncs          !< scaling parameter nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qis          !< scaling parameter qi
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nis          !< scaling parameter ni
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrs          !< scaling parameter qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrs          !< scaling parameter nr

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol           !< Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib          !< Richardson bulk number

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0           !< roughness length for momentum
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0h          !< roughness length for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0q          !< roughness length for humidity

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt1          !< potential temperature at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qv1          !< mixing ratio at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt1         !< virtual potential temperature at first grid level
!
!--    Pre-defined arrays for ln(z/z0)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ln_z_z0      !< ln(z/z0)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ln_z_z0h     !< ln(z/z0h)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ln_z_z0q     !< ln(z/z0q)
!
!--    Define arrays for surface fluxes
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  usws      !< vertical momentum flux usws for u-component at horizontal surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsws      !< vertical momentum flux vsws for v-component at horizontal surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wsus_wsvs !< vertical momentum flux wsus and wsvs for w-component at vertical surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  usvs      !< momentum flux usvs  for the u-component at north/south-adjacent vertical surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsus      !< momentum flux vsuss for the v-component at east/west-adjacent vertical surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf       !< surface flux sensible heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws      !< surface flux latent heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ssws      !< surface flux passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcsws     !< surface flux qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncsws     !< surface flux nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qisws     !< surface flux qi
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nisws     !< surface flux ni
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrsws     !< surface flux qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrsws     !< surface flux nr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  sasws     !< surface flux salinity

!
!--    For surface-atmosphere coupling through diffusion_s, special arrays are introduced in order
!--    to allow for aggregated surfaces (mixed tiles, currently in the case of SLUrb). In the time
!--    integration scheme, the surface_layer_fluxes, which computes the scaling parameters for the
!--    LSM surfaces, is called before LSM itself. Thus modifying the fluxes after LSM is called
!--    will affect the computation of the scaling parameters in surface_layer_fluxes, which is
!--    against the tile approach, where the tiles should be independent from each other.
!--    Therefore, additional arrays to store the aggregated values (shf_agg, qsws_agg) are
!--    introduced here. These are used for surface-atmosphere coupling by the diffusion_s routine.
!--    LSM and surface_layer_fluxes continue to use the unaggregated values. For momentum fluxes
!--    these are not needed as surface_layer_fluxes has no dependency on them if the surface
!--    models are used.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_agg  !< surface flux aggregated latent heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_agg   !< surface flux aggregated sensible heat, needed in case of urban modification
                                                         !< of LSM flux by another surface model, e.g. SLUrb. Otherwise equal to shf.

!
!--    Arrays to represent the normal vector and the wall location.
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  n_eff    !< normal vector component on grid center in the respective facing direction

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  n_s             !< full normal vector on grid center
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wall_location_c !< coordinates for the center point of the surface
!
!--    Data required for the cut-cell method.
!--    Polygon vertices locations, slanted faces.
       INTEGER(iwp) ::  num_faces_vert !< maximum number of vertices in polygon
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  bid       !< building ID associated with cut-cell surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  btype     !< building type associated with cut-cell surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_edges !< list of number of edges
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ptype     !< pavement type associated with cut-cell surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  vtype     !< vegetation type associated with cut-cell surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  wtype     !< water type associated with cut-cell surface

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  faces   !< vertex data defining the cut-cell face

       LOGICAL, DIMENSION(:), ALLOCATABLE ::  accessed           !< flag to identify that cut-cell data has been set on surface
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  cut_cell_wall      !< flag indicating cut-cell wall surface
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  cut_cell_roof      !< flag indicating cut-cell roof surface

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  area      !< list of all normal vector coordinates
!
!--    Surface fluxes for chemistry
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  css    !< scaling parameter chemical species
!
!--    Surface fluxes for SALSA
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  answs  !< surface flux aerosol number: dim 1: flux, dim 2: bin
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  amsws  !< surface flux aerosol mass: dim 1: flux, dim 2: bin
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gtsws  !< surface flux gaseous tracers: dim 1: flux, dim 2: gas
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  cssws  !< surface flux chemical species
!
!--    Surface-related variables for DET
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dmsws         !< surface total flux dust mass: dim 1: flux, dim 2: bin
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dm_depo_flux  !< surface deposition flux: dim 1: bin, dim 2: flux
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dm_emis_flux  !< surface emission flux: dim 1: bin, dim 2: flux
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  us_t          !< surface threshold friction velocity: dim 1: velocity, dim 2: bin
!
!--    Required for horizontal walls in production_e
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_0  !< virtual velocity component (see production_e_init for further explanation)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_0  !< virtual velocity component (see production_e_init for further explanation)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  mom_flux_tke  !< momentum flux usvs, vsus, wsus, wsvs at vertical surfaces at grid
                                                               !< center (used in production_e)
!
!--    Variables required for LSM as well as for USM
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  building_type_name    !< building type name at surface element
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  pavement_type_name    !< pavement type name at surface element
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  vegetation_type_name  !< water type at name surface element
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  water_type_name       !< water type at name surface element

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nzt_pavement     !< top index for pavement in soil
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  building_type    !< building type at surface element
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  pavement_type    !< pavement type at surface element
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  vegetation_type  !< vegetation type at surface element
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  water_type       !< water type at surface element

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  albedo_type  !< albedo type, for each fraction
                                                                  !< (wall,green,window or vegetation,pavement water)

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  building_surface  !< flag parameter indicating that the surface element is covered
                                                                 !< by buildings (no LSM actions, not implemented yet)
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  building_covered  !< flag indicating that buildings are on top of orography,
                                                                 !< only used for vertical surfaces in LSM
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  pavement_surface    !< flag parameter for pavements
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  water_surface       !< flag parameter for water surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  vegetation_surface  !< flag parameter for natural land surfaces

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  albedo    !< broadband albedo for each surface fraction
                                                           !< (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  albedo_uv !< broadband albedo in the UV-spectral range for each surface fraction
                                                           !< (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  emissivity  !< emissivity of the surface, for each fraction
                                                             !< (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  frac  !< relative surface fraction
                                                       !< (LSM: vegetation, water, pavement; USM: wall, green, window)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  aldif       !< albedo for longwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  aldir       !< albedo for longwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  asdif       !< albedo for shortwave diffusive radiation, solar angle of 60 deg.
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  asdir       !< albedo for shortwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_aldif  !< albedo for longwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_aldir  !< albedo for longwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_asdif  !< albedo for shortwave diffusive radiation, solar angle of 60 deg.
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_asdir  !< albedo for shortwave direct radiation, solar angle of 60 degrees
!
!--    Define arrays for soil and wall-layer depth. At the moment soil layer depth can only be
!--    one-dimensional.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zs  !< soil layer depths (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw  !< wall layer depths (m)

#if defined( __tenstream )
!
!--    Declare TenStream related variables.
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  ts_albedo !< TS: albedo for each surface fraction, which will be used finally by TS
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  ts_aldif  !< TS: albedo for longwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  ts_aldir  !< TS: albedo for longwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  ts_asdif  !< TS: albedo for shortwave diffusive radiation, solar angle of 60 deg.
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  ts_asdir  !< TS: albedo for shortwave direct radiation, solar angle of 60 degrees

#endif
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  q_surface        !< skin-surface mixing ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  pt_surface       !< skin-surface temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  vpt_surface      !< skin-surface virtual temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  rad_net          !< net radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  rad_net_l        !< net radiation, used in USM
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h         !< heat conductivity of soil/ wall (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_layer   !< heat conductivity of soil/ wall interpolated to layer edge (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_green   !< heat conductivity of green soil (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_window  !< heat conductivity of windows (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_window_layer  !< heat conductivity of windows interpolated to layer edge(W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_def     !< default heat conductivity of soil (W/m/K)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_in   !< incoming longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out  !< emitted longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_dif  !< incoming longwave radiation from sky
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_ref  !< incoming longwave radiation from reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_res  !< residual longwave radiation in surface after last reflection step
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in   !< incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_out  !< emitted shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_dir  !< direct incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_dif  !< diffuse incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_ref  !< incoming shortwave radiation from reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_res  !< residual shortwave radiation in surface after last reflection step

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq             !< liquid water coverage (of vegetated area)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_veg             !< vegetation coverage
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  f_sw_in           !< fraction of absorbed shortwave radiation by the surface layer
                                                                 !< (not implemented yet)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf               !< ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  g_d               !< coefficient for dependence of r_canopy
                                                                 !< on water vapour pressure deficit
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lai               !< leaf area index
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surface_u  !< coupling between surface and soil (depends on vegetation type)
                                                                 !< (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surface_s  !< coupling between surface and soil (depends on vegetation type)
                                                                 !< (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_liq          !< surface flux of latent heat (liquid water portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_soil         !< surface flux of latent heat (soil portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_veg          !< surface flux of latent heat (vegetation portion)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a           !< aerodynamic resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_green     !< aerodynamic resistance at green fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_window    !< aerodynamic resistance at window fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy      !< canopy resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil        !< soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil_min    !< minimum soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_s           !< total surface resistance (combination of r_soil and r_canopy)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy_min  !< minimum canopy (stomatal) resistance

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_10cm  !< near surface air potential temperature at distance 10 cm from
                                                        !< the surface (K)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  alpha_vg         !< coef. of Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_w         !< hydraulic diffusivity of soil (?)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w          !< hydraulic conductivity of soil (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_sat      !< hydraulic conductivity at saturation
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  l_vg             !< coef. of Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_fc             !< soil moisture at field capacity (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_res            !< residual soil moisture
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_sat            !< saturation soil moisture (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_wilt           !< soil moisture at permanent wilting point (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  n_vg             !< coef. Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total_def  !< default volumetric heat capacity of the (soil) layer (J/m3/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total      !< volumetric heat capacity of the actual soil matrix (J/m3/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  root_fr          !< root fraction within the soil layers

!--    Indoor model variables
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_prev      !< indoor temperature for facade element
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  waste_heat  !< waste heat
!
!--    Urban surface variables
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  surface_types  !< array of types of wall parameters

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  isroof_surf   !< flag indicating roof surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  gfl           !< flag indicating ground floor level surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  target_temp_summer  !< indoor target temperature summer
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  target_temp_winter  !< indoor target temperature summer

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface           !< heat capacity of the wall surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface_green     !< heat capacity of the green surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface_window    !< heat capacity of the window surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  green_type_roof     !< type of the green roof
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf         !< heat conductivity between air and surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf_green   !< heat conductivity between air and green surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf_window  !< heat conductivity between air and window surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_wall      !< thickness of the wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_green     !< thickness of the green wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_window    !< thickness of the window wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  transmissivity      !< transmissivity of windows

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsl  !< reflected shortwave radiation for local surface in i-th reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutll  !< reflected + emitted longwave radiation for local surface
                                                          !< in i-th reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfhf     !< total radiation flux incoming to minus outgoing from local surface

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_wall_m    !< surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_window_m  !< window surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_green_m   !< green surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf_eb              !< wall heat flux of sensible heat in wall normal direction

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb          !< wall ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_window   !< window ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_green    !< green ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb         !< indoor wall ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_window  !< indoor window ground heat flux

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out_change_0  !<

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw   !< shortwave radiation falling to local surface including radiation
                                                          !< from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw  !< total shortwave radiation outgoing from nonvirtual surfaces surfaces
                                                          !< after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw   !< longwave radiation falling to local surface including radiation from
                                                          !< reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw  !< total longwave radiation outgoing from nonvirtual surfaces surfaces
                                                          !< after all reflection

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  n_vg_green      !< vangenuchten parameters
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  alpha_vg_green  !< vangenuchten parameters
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_vg_green      !< vangenuchten parameters


       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_wall         !< volumetric heat capacity of the material ( J m-3 K-1 )
                                                                    !< (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall            !< wall grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_wall           !< 1/dz_wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall_center     !< wall grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_wall_center    !< 1/dz_wall_center
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_wall_m          !< t_wall prognostic array
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_window       !< volumetric heat capacity of the window material ( J m-3 K-1 )
                                                                    !< (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_window          !< window grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_window         !< 1/dz_window
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_window_center   !< window grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_window_center  !< 1/dz_window_center
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_window_m        !< t_window prognostic array
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_window          !< window layer depths (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_green        !< volumetric heat capacity of the green material ( J m-3 K-1 )
                                                                    !< (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total_green  !< volumetric heat capacity of the moist green material
                                                                    !< ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_green           !< green grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_green          !< 1/dz_green
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_green_center    !< green grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_green_center   !< 1/dz_green_center
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_green_m         !< t_green prognostic array
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_green           !< green layer depths (m)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_green_sat    !< hydraulic conductivity
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_w_green       !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_w_green_layer !< lambda_w_green at center of layer
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_green        !< hydraulic conductivity
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_green_layer  !< gamma_w_green at center of layer
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tswc_m               !<

    END TYPE surf_type

!
!-- Derived type for the SLUrb model.
    TYPE surf_type_slurb
!
!--    Grid definition.
       INTEGER(iwp) ::  ns  !< total number of SLUrb tiles

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: i  !< i index of the tile
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: j  !< j index of the tile

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: m  !< tile index m matching the tile j,i

       REAL(wp), DIMENSION(:), ALLOCATABLE :: dt_max  !< time step limit for model physical processes (s)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_road  !< road layer thickness (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_roof  !< roof layer thickness (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall  !< wall layer thickness (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_win   !< window layer thickness (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_win   !< cumulative window thickness (m)
!
!--    Tile-aggregated quantities.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  albedo_urb       !< effective urban albedo
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  emiss_urb        !< effective urban emissivity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_urb           !< urban Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_urb         !< total urban latent heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_in_urb   !< incoming longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out_urb  !< outgoing longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in_urb   !< incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_out_urb  !< outgoing shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ram_urb         !< urban aerodynamic resistance for momentum
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_urb         !< urban bulk-Richardson number
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_urb         !< total urban sensible heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_2m_urb        !< urban 2-metre temperature (extrapolated) (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_c_urb         !< complete (area-weighted) urban surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_h_urb         !< effective urban surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_rad_urb       !< urban radiative surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  usws_urb        !< urban momentum flux (u-component)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsws_urb        !< urban momentum flux (v-component)
!
!--    Model prognostic variables.
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  m_liq_road    !< liquid water reservoir on roads
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  m_liq_road_p  !< prog. liquid water reservoir on roads
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  m_liq_roof    !< liquid water reservoir on roofs
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  m_liq_roof_p  !< prog. liquid water reservoir on roofs
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  q_can         !< canyon mixing ratio (kg/kg)
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  q_can_p       !< prognostic canyon mixing ratio (kg/kg)
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  t_can         !< canyon air temperature (K)
       REAL(wp), DIMENSION(:), POINTER, CONTIGUOUS   ::  t_can_p       !< prog. canyon temperature (K)

       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_road      !< road temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_road_p    !< prog. road temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_roof      !< roof temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_roof_p    !< prog. roof temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_wall_a    !< wall A temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_wall_a_p  !< prog. wall A temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_wall_b    !< wall B temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_wall_b_p  !< prog. wall B temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_win_a     !< window A temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_win_a_p   !< prog. window A temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_win_b     !< window B temperature (K)
       REAL(wp), DIMENSION(:,:), POINTER, CONTIGUOUS ::  t_win_b_p   !< prog. window B temperature (K)
!
!--    Tendencies of the prognostic variables.
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  tm_liq_road  !< road liquid water reservoir tendency
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  tm_liq_roof  !< roof liquid water reservoir tendency
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  tq_can       !< canyon mixing ratio tendency (kg/kg/s)
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  tt_can       !< canyon temperature tendency (K/s)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_road    !< road temperature tendency (K/s)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_roof    !< roof temperature tendency (K/s)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_wall_a  !< wall A temperature tendency (K/s)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_wall_b  !< wall B temperature tendency (K/s)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_win_a   !< window A temperature tendency (K/s)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_win_b   !< window B temperature tendency (K/s)
!
!--    Diagnostic surface thermodynamic variables.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_road    !< road surface potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_roof    !< roof surface potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_wall_a  !< wall A surface potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_wall_b  !< wall B surface potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_win_a   !< window A surface potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_win_b   !< window A surface potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_road     !< road surface mixing ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_roof     !< roof surface mixing ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs_road    !< road surface saturation mixing ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs_roof    !< roof surface saturation mixing ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt_road   !< road surface virtual potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt_roof   !< roof surface virtual potential temperature
!
!--    Diagnostic internal sensible heat fluxes.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_can       !< sensible heat flux between the street canyon and the atmosphere
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_external  !< sensible heat flux external to the model (e.g. industry)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_road      !< road surface sensible heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_roof      !< roof surface sensible heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_traffic   !< traffic sensible heat flux (input-only)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_wall_a    !< wall A sensible heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_wall_b    !< wall B sensible heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_win_a     !< window A sensible heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_win_b     !< window B sensible heat flux
!
!--    Diagnostic internal latent heat fluxes.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_can       !< latent heat flux between the street canyon and the atmosphere
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_external  !< latent heat flux external to the model (e.g. industry)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_liq_road  !< roof latent heat flux (liquid incl. precipitation)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_liq_roof  !< roof latent heat flux (liquid incl. precipitation)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_road      !< road latent heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_roof      !< roof latent heat flux
!
!--    Liquid water coverages (storages).
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq_road  !< liquid water coverage on road
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq_roof  !< liquid water coverage on roof
!
!--    Diagnostic ground heat fluxes.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_road    !< road ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_roof    !< roof indoor heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_wall_a  !< wall A indoor heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_wall_b  !< wall B indoor heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_win_a   !< window A indoor heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf_win_b   !< window B indoor heat flux
!
!--    Model internal radiation fluxes.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_can     !< net longwave radiative at canyon top (downwards)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_road    !< net longtwave radiative flux on road
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_roof    !< net longwave radiative flux on roof
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_urb     !< urban aggegated net longwave radiative flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_wall_a  !< net longwave radiative flux on wall A
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_wall_b  !< net longwave radiative flux wall B
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_win_a   !< net longwave radiative flux on wall A
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_net_win_b   !< net longwave radiative flux window B
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in_road     !< incoming shortwave radiative flux on road
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in_win_a    !< incoming shortwave radiative flux on window A
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in_win_b    !< incoming shortwave radiative flux on window B
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_road    !< net shortwave radiative flux on road
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_roof    !< net shortwave radiative flux on roof
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_urb     !< urban aggegated net shortwave radiative flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_wall_a  !< net shortwave radiative flux on wall A
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_wall_b  !< net shortwave radiative flux on wall B
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_win_a   !< net shortwave radiative flux on window A
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_net_win_b   !< net shortwave radiative flux on wall B
!
!--    Surface layer model diagnostic variables.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_can      !< canyon top Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_road     !< road Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol_roof     !< rroof Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_can      !< street canyon virtual potential temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_can     !< canyon top bulk Richardson number
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_road    !< road bulk Richardson number
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib_roof    !< roof bulk Richardson number
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_can      !< friction velocity for canyon resistance calculation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uv_abs_can  !< horizontal wind speed in street caynon at half-height
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uv_eff_can  !< effective horizontal wind speed in street canyon at half-height
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt_can     !< street canyon virtual potential temperature (K)
!
!--    Aerodynamic resistances for heat.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_can     !< street canyon air aerodynamic resistance for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_facade  !< wall and window aerodynamic resistance for heat (combined)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_road    !< road aerodynamic resistance for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_roof    !< roof aerodynamic resistance for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_wall_a  !< wall A aerodynamic resistance for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_wall_b  !< wall B aerodynamic resistance for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_win_a   !< wall A aerodynamic resistance for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rah_win_b   !< wall B aerodynamic resistance for heat
!
!--    Local friction velocities for roofs and roads.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_road  !< friction velocity for roads
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_roof  !< friction velocity for roofs
!
!--    Diagnostic variables, defined at the first atmospheric grid level.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt1      !< potential temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  q1       !< specific humidity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us_urb   !< friction velocity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uv_abs1  !< horizontal wind speed
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uv_eff1  !< effective horizontal wind speed
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt1     !< virtual potential temperature
!
!--    Parameters for the whole urban tile.
       LOGICAL,  DIMENSION(:), ALLOCATABLE ::  anisotropic_canyon  !< boolean flag to mark anisotropic canyon
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  f_bld               !< fractional area occupied by buldings (plan area fraction)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  f_bld_frn           !< frontal area fraction of buildings
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  f_win               !< window fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  h_bld               !< building height
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  hw_can              !< canyon aspect ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  svf_road            !< sky-view factor for road
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  svf_wall            !< sky-view-factor for walls
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  theta_can           !< canyon orientation / road direction in radians
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0_urb              !< aerodynamic roughness length of the urban surface
!
!--    Material properties.
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  albedo_road         !< albedo of the road
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  albedo_roof         !< albedo of the roof
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  albedo_wall         !< albedo of the wall
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  albedo_wall_win     !< weighted average of wall and window albedos for reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  albedo_win          !< albedo of the window
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  emiss_road          !< emissivity of the road
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  emiss_roof          !< emissivity of the roof
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  emiss_wall          !< emissivity of the wall
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  emiss_win           !< emissivity of the window
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  transmissivity_win  !< transmissivity of the window layers
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z0_road             !< aerodynamic roughness length for momentum for roads
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z0_roof             !< aerodynamic roughness length for momentum of roofs
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z0_wall             !< aerodynamic roughness length for walls and windows
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z0h_road            !< aerodynamic roughness length for heat for roads
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z0h_roof            !< aerodynamic roughness length for heat for roofs

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  absorption_win  !< fraction of absorbed shortwave radiation over glass sheet
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  c_road          !< total heat capacity of the road
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  c_roof          !< total (specific c * layer depth) heat capacity of the roof
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  c_wall          !< total heat capacity of the wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  c_win           !< total heat heat capacity of the window
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_road     !< thermal conductivity of the road
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_roof     !< thermal conductivity of the roof
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_wall     !< thermal conductivity of the wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_win      !< effective thermal conductivity of the window
!
!--    Building indoor parameters.
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  t_indoor  !< building indoor temperature (K)
!
!--    Soil parameters.
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  t_soil  !< fixed soil top temperature (K)
!
!--    Pre-computed total layer conductivities.
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  conductivity_road  !< total conductivity bewtween road layers (lambda_h / dz)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  conductivity_roof  !< total conductivity between roof layers (lambda_h / dz)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  conductivity_wall  !< total conductivity between wall layers (lambda_h / dz)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  conductivity_win   !< total conductivity between window layers (lambda_h / dz)
!
!--    Pre-computed variables and coefficients.
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  sw_ref_denom      !< SW radiation reflection denominator
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  uv_abs_can_coef   !< coefficient for the canyon wind speed
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  wall_hor_a_ratio  !< wall-to-horizontal area ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z_mo              !< reference height for MOST for the atmosphere
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  z_mo_can          !< canyon reference height for MOST (canyon half-height)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lw_road_coef  !< LW radiation coefficients for roads
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lw_roof_coef  !< LW radiation coefficients for roofs
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lw_wall_coef  !< LW radiation coefficients for walls
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lw_win_coef   !< LW radiation coefficients for walls

    END TYPE surf_type_slurb
!
!-- Data structure to link grid-line surfaces with energy-balance surfaces, where an
!-- energy-balance surface may interact with more than one atmosphere cell.
    TYPE grid_line_to_eb_surface
       INTEGER(iwp) :: n_gls !< number of grid-line surfaces (or atmosphere cells affected) represented by an EB-surface

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_a   !< x-indices of adjacent atmosphere cells
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_a   !< y-indices of adjacent atmosphere cells
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_a   !< z-indices of adjacent atmosphere cells

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  m_ind !< reference indices in surface types

       LOGICAL, DIMENSION(:), ALLOCATABLE ::  def_surface !< flag indicating a default surface
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  lsm_surface !< flag indicating a natural surface
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  usm_surface !< flag indicating a building surface

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_a !< potential temperature in reference atmosphere cell
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_a  !< mixing ratio in reference atmosphere cell
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf  !< surface sensible flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws !< surface latent flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us   !< friction velocity in the m-th surface

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  nv !< normal vector
    END TYPE grid_line_to_eb_surface
!
!-- Data structure used to build a connection between cut-cell surfaces and LSM / USM surfaces.
    TYPE index_list

       INTEGER(iwp) ::  n_accessed = 0 !< number of LSM or USM surfaces that access the cut-cell surface

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  m_list !< list of indices

    END TYPE index_list

!
!-- Data structure that encompasses all surfaces necessary to solve the energy balance
    TYPE surf_type_cct

       CHARACTER(LEN=3), DIMENSION(:,:), ALLOCATABLE ::  face_accessed_by !< type of surface that accesses the cut-cell face

       LOGICAL, DIMENSION(:), ALLOCATABLE ::  face_accessed !< flag indicating if the cut-cell surface is accessed via a prognostic grid point

       INTEGER(iwp) ::  ns = 0          !< number of energy-balance surfaces

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i           !< x-index where the surface is located
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j           !< y-index where the surface is located
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k           !< z-index where the surface is located
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  m_index_ref !< index with respect to nearest reference energy-balance active surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_access  !< number of LSM and USM surfaces that access the cut-cell face

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  end_index    !< End index within surface data type for given (j,i)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  index_access !< corresponding indices that access the cut-cell face
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  start_index  !< Start index within surface data type for given (j,i)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  albedo              !< broadband albedo
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  emissivity          !< emissivity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_surface          !< skin-surface temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_dif          !< incoming longwave radiation (from sky)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_in           !< incoming longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out          !< outgoing longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out_change_0 !< derivative of outgoing longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_ref          !< incoming longwave radiation from reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_res          !< residual longwave radiation in surface after last reflection step
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_net             !< net radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_dif          !< incoming shortwave radiation (diffuse)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_dir          !< incoming shortwave radiation (direct)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in           !< incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_out          !< outgoing shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_ref          !< incoming shortwave radiation from reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_res          !< residual shortwavewave radiation in surface after last reflection step
!
!--    Geometric data.
!--    Polygon vertices locations, slanted faces.
       INTEGER(iwp)                              ::  num_faces_vert !< maximum number of vertices in polygon
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  faces
       INTEGER(iwp), DIMENSION(:),   ALLOCATABLE ::  num_edges !< list of number of edges

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  area            !< list of all normal vector coordinates
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  normals         !< list of all normal vector coordinates
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  n_s             !< normal vector on grid center
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wall_location_c !< coordinates for the center point of the surface
!
!--    Data type required for tranferring data between cut-cell surfaces as treated in the RTM and the
!--    actual LSM and USM surfaces.
       TYPE( index_list ), DIMENSION(:), ALLOCATABLE ::  accessed_lsm !< accessing LSM surfaces
       TYPE( index_list ), DIMENSION(:), ALLOCATABLE ::  accessed_usm !< accessing USM surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us   !< friction velocity in the m-th surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts           !< scaling parameter temerature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs           !< scaling parameter humidity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ss           !< scaling parameter passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcs          !< scaling parameter qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncs          !< scaling parameter nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qis          !< scaling parameter qi
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nis          !< scaling parameter ni
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrs          !< scaling parameter qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrs          !< scaling parameter nr

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol           !< Obukhov length

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0           !< roughness length for momentum
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0h          !< roughness length for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0q          !< roughness length for humidity

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt1          !< potential temperature at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qv1          !< mixing ratio at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt1         !< virtual potential temperature at first grid level
!
!--    Define arrays for surface fluxes
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  usws      !< vertical momentum flux usws for u-component at horizontal surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsws      !< vertical momentum flux vsws for v-component at horizontal surfaces

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf      !< surface flux sensible heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws     !< surface flux latent heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ssws     !< surface flux passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcsws    !< surface flux qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncsws    !< surface flux nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qisws    !< surface flux qi
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nisws    !< surface flux ni
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrsws    !< surface flux qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrsws    !< surface flux nr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  sasws    !< surface flux salinity

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  q_surface        !< skin-surface mixing ratio

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  vpt_surface      !< skin-surface virtual temperature

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_10cm  !< near surface air potential temperature at distance 10 cm from


       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf               !< ground heat flux

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a           !< aerodynamic resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_green     !< aerodynamic resistance at green fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy      !< canopy resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil        !< soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil_min    !< minimum soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_s

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uvw_abs      !< absolute surface-parallel velocity on grid center
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  waste_heat  !< waste heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb         !< indoor wall ground heat flux

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  frac  !< relative surface fraction

    END TYPE surf_type_cct

    TYPE (bc_type) ::  bc_hv  !< data structure that includes all surface-adjacent grid points, used to set boundary conditions

    TYPE(surf_out_type) ::  surf_out  !< variable which contains all surface output information

    TYPE(surf_type), TARGET ::  surf_def  !< default surfaces
    TYPE(surf_type), TARGET ::  surf_lsm  !< land surfaces
    TYPE(surf_type), TARGET ::  surf_top  !< model top surfaces (e.g. for ocean simulations or ocean-atmosphere coupling)
    TYPE(surf_type), TARGET ::  surf_usm  !< building surfaces

    TYPE(surf_type), TARGET ::  surf_u !< surface on u-grid in case of cut-cell approach
    TYPE(surf_type), TARGET ::  surf_v !< surface on v-grid in case of cut-cell approach
    TYPE(surf_type), TARGET ::  surf_w !< surface on w-grid in case of cut-cell approach

    TYPE( surf_type_cct ), TARGET ::  surf_cct !< cut-cell surfaces that won't be accessed

    TYPE (surf_type_slurb) ::  surf_slurb   !< SLUrb surfaces

    INTEGER(iwp), PARAMETER ::  ind_veg_wall  = 0  !< index for vegetation / wall-surface fraction, used for access of albedo,
                                                   !< emissivity, etc., for each surface type
    INTEGER(iwp), PARAMETER ::  ind_pav_green = 1  !< index for pavement / green-wall surface fraction, used for access of albedo,
                                                   !< emissivity, etc., for each surface type
    INTEGER(iwp), PARAMETER ::  ind_wat_win   = 2  !< index for water / window-surface fraction, used for access of albedo,
                                                   !< emissivity, etc., for each surface type

    LOGICAL ::  vertical_surfaces_exist     = .FALSE.  !< flag indicating that there are vertical urban/land surfaces
                                                       !< in the domain (required to activiate RTM)

    LOGICAL ::  surf_bulk_cloud_model       = .FALSE.  !< use cloud microphysics
    LOGICAL ::  surf_microphysics_morrison  = .FALSE.  !< use 2-moment Morrison (add. prog. eq. for nc and qc)
    LOGICAL ::  surf_microphysics_seifert   = .FALSE.  !< use 2-moment Seifert and Beheng scheme
    LOGICAL ::  surf_microphysics_ice_phase = .FALSE.  !< use 2-moment Seifert and Beheng scheme

    REAL(wp), DIMENSION(0:20) ::  soil_moisture = -9999.0  !< NAMELIST soil moisture content (m3/m3)

!
!-- Urban fraction for SLUrb/DCEP.
    REAL(wp), ALLOCATABLE ::  fr_urb(:,:)  !< fraction of urban parts in a grid element

!
!-- DCEP related variables.
    REAL(wp), ALLOCATABLE ::  albedo_urb(:,:)   !< effective urban albedo
    REAL(wp), ALLOCATABLE ::  albedop_urb(:,:)  !< effective urban albedo
    REAL(wp), ALLOCATABLE ::  emiss_urb(:,:)    !< urban emissivity
    REAL(wp), ALLOCATABLE ::  t_grad_urb(:,:)   !< effective urban radiation temperature


    SAVE

    PRIVATE
!
!-- Public variables
    PUBLIC albedo_urb,                                                                             &
           albedop_urb,                                                                            &
           bc_hv,                                                                                  &
           emiss_urb,                                                                              &
           fr_urb,                                                                                 &
           ind_pav_green,                                                                          &
           ind_veg_wall,                                                                           &
           ind_wat_win,                                                                            &
           soil_moisture,                                                                          &
           surf_bulk_cloud_model,                                                                  &
           surf_cct,                                                                               &
           surf_def,                                                                               &
           surf_lsm,                                                                               &
           surf_microphysics_ice_phase,                                                            &
           surf_microphysics_morrison,                                                             &
           surf_microphysics_seifert,                                                              &
           surf_out_type,                                                                          &
           surf_out,                                                                               &
           surf_slurb,                                                                             &
           surf_top,                                                                               &
           surf_type,                                                                              &
           surf_type_cct,                                                                          &
           surf_u,                                                                                 &
           surf_usm,                                                                               &
           surf_v,                                                                                 &
           surf_w,                                                                                 &
           t_grad_urb,                                                                             &
           vertical_surfaces_exist

 END MODULE surface_mod
