!> @file init_vertical_profiles.f90
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
! Copyright 2017-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Siegfried Raasch
!
! Description:
! ------------
!> Inititalizes the vertical profiles of scalar quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_vertical_profiles( vertical_gradient_level_ind, vertical_gradient_level,          &
                                    vertical_gradient, initial_profile, surface_value,             &
                                    bc_top_gradient, quantity )
    USE arrays_3d,                                                                                 &
        ONLY:  dzu,                                                                                &
               zu

    USE control_parameters,                                                                        &
        ONLY:  neutral,                                                                            &
               ocean_mode

    USE indices,                                                                                   &
        ONLY:  nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE kinds

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  quantity  !< name of scalar quantity to be treated

    INTEGER(iwp), DIMENSION(1:12), INTENT(INOUT) ::  vertical_gradient_level_ind  !< vertical grid indices for gradient levels

    REAL(wp), INTENT(OUT) ::  bc_top_gradient  !< model top gradient
    REAL(wp), INTENT(IN)  ::  surface_value    !< surface value of the respecitve quantity

    REAL(wp), DIMENSION(0:nz+1), INTENT(OUT)   ::  initial_profile          !< initialisation profile
    REAL(wp), DIMENSION(1:12),   INTENT(INOUT) ::  vertical_gradient        !< given vertical gradient
    REAL(wp), DIMENSION(1:12),   INTENT(INOUT) ::  vertical_gradient_level  !< given vertical gradient level

    INTEGER(iwp) ::  i     !< loop counter
    INTEGER(iwp) ::  imax  !< maximum number of scalar levels where gradients are changing, including model bottom and top
    INTEGER(iwp) ::  k     !< loop counter

    REAL(wp), DIMENSION(1:12) ::  scalar_gradient  !< vertical gradient to be used above the respective index
    REAL(wp), DIMENSION(1:12) ::  scalar_level     !< heights at which gradients are changing
    REAL(wp), DIMENSION(1:12) ::  scalar_value     !< scalar values at heights where gradients are changing


    scalar_level = HUGE( 1.0_wp )
!
!-- If no gradients are given, set initial profiles to surface value and gradient to zero gradient.
    IF ( vertical_gradient_level(1) == HUGE( 1.0_wp )  .OR.                                        &
         ( TRIM( quantity ) == 'pt'  .AND.  neutral ) )                                            &
    THEN
       initial_profile(:)             = surface_value
       vertical_gradient_level(1)     = 0.0_wp
       vertical_gradient_level_ind(1) = 0
       bc_top_gradient                = 0.0_wp
       RETURN
    ENDIF

    IF ( .NOT. ocean_mode )  THEN

!
!--    First calculated and store scalar values at heights where the gradients are changing.
!--    Add bottom level, if first gradient is given for z > 0.
       IF ( vertical_gradient_level(1) > 0.0_wp )  THEN
          scalar_level(1)       = 0.0_wp
          scalar_gradient(1)    = 0.0_wp
          scalar_level(2:11)    = vertical_gradient_level(1:10)
          scalar_gradient(2:11) = vertical_gradient(1:10)
          imax = 11
       ELSE
          scalar_level(1:10)    = vertical_gradient_level(1:10)
          scalar_gradient(1:10) = vertical_gradient(1:10)
          imax = 10
       ENDIF
!
!--    Remove undefined levels or levels above domain top.
       DO WHILE ( scalar_level(imax) == HUGE( 1.0_wp )  .OR.  scalar_level(imax) > zu(nzt+1) )
          imax = imax - 1
       ENDDO
!
!--    Add point at model top, if uppermost defined level is below the top.
       IF ( scalar_level(imax) < zu(nzt+1) )  THEN
          imax = imax + 1
          scalar_level(imax)    = zu(nzt+1)
          scalar_gradient(imax) = 0.0_wp
       ENDIF
!
!--    Now calculate scalar values for all heights that have been determined above.
       scalar_value(1) = surface_value
       DO  i = 2, imax
          scalar_value(i) = scalar_value(i-1) + scalar_gradient(i-1) / 100.0_wp *                  &
                                                ( scalar_level(i) - scalar_level(i-1) )
       ENDDO
!
!--    Finally, calculate scalar values at points of the vertical model grid by linearly
!--    interpolating between the values calculated above for those heights where gradients are
!--    changing.
       i = 1
       initial_profile(nzb) = scalar_value(1)
       DO  k = nzb+1, nzt
          DO WHILE ( scalar_level(i+1) <= zu(k)  .AND.  i < imax )
             i = i + 1
          ENDDO
          initial_profile(k) = scalar_value(i) + ( scalar_value(i+1) - scalar_value(i) ) *         &
                                                 ( zu(k) - scalar_level(i) ) /                     &
                                                 ( scalar_level(i+1) - scalar_level(i) )
       ENDDO
       initial_profile(nzt+1) = scalar_value(imax)
!
!--    Store bottom/top levels that may have been additionally added above, and determine grid
!--    indices at (or below) gradients are applied. These quantities (gradient/gradient_level/
!--    gradient_level_ind) are only used for output purpose in the _rc file.
       vertical_gradient(:)       = scalar_gradient(:)
       vertical_gradient_level(:) = scalar_level(:)
       DO  i = 1, 12
          IF ( vertical_gradient_level(i) == HUGE( 1.0_wp ) )  EXIT
          DO  k = nzb, nzt+1
             IF ( vertical_gradient_level(i) >= zu(k) )  THEN
                vertical_gradient_level_ind(i) = k
             ELSE
                EXIT
             ENDIF
          ENDDO
       ENDDO
!
!--    Store gradient at the top boundary for possible Neumann boundary condition.
       bc_top_gradient  = ( initial_profile(nzt+1) - initial_profile(nzt) ) / dzu(nzt+1)

    ELSE

!
!--    First calculated and store scalar values at heights where the gradients are changing.
!--    Also store values at bottom and top of the domain.
       IF ( vertical_gradient_level(1) < 0.0_wp )  THEN
!
!--       Add ocean surface level, if first gradient is given for z < 0.
          scalar_level(1)       = 0.0_wp
          scalar_gradient(1)    = 0.0_wp
          scalar_level(2:11)    = vertical_gradient_level(1:10)
          scalar_gradient(2:11) = vertical_gradient(1:10)
          imax = 11
       ELSE
          scalar_level(1:10)    = vertical_gradient_level(1:10)
          scalar_gradient(1:10) = vertical_gradient(1:10)
          imax = 10
       ENDIF
!
!--    Remove undefined levels or levels below the bottom of the domain.
       DO WHILE ( scalar_level(imax) == HUGE( 1.0_wp )  .OR.  scalar_level(imax) < zu(0) )
          imax = imax - 1
       ENDDO
!
!--    Add point at model bottom, if the lowest defined level is above the bottom.
       IF ( scalar_level(imax) > zu(0) )  THEN
          scalar_level(imax+1) = zu(0)
          imax = imax + 1
       ENDIF
!
!--    Now calculate scalar values for all heights that have been determined above.
       scalar_value(1) = surface_value
       DO  i = 2, imax
          scalar_value(i) = scalar_value(i-1) - scalar_gradient(i-1) / 100.0_wp *                  &
                                                ( scalar_level(i-1) - scalar_level(i) )
       ENDDO
!
!--    Finally, calculate scalar values at vertical grid points by linearly interpolating between
!--    the values calculated above for those heights where gradients are changing.
!--    In ocean mode, profiles are constructed starting from the ocean surface, which is at the top
!--    of the model domain.
!--    Note that due to the staggered grid there is no scalar point directly defined at the
!--    sea surface (nzt+1 is half a grid spacing "above" sea level).
       i = 1
       initial_profile(nzt+1) = scalar_value(1) + 0.5_wp * dzu(nzt+1) *                            &
                                                  scalar_gradient(1) / 100.0_wp
       vertical_gradient_level_ind(1) = 0
       DO  k = nzt, nzb+1, -1
          DO WHILE ( scalar_level(i+1) >= zu(k)  .AND.  i < imax )
             i = i + 1
             IF ( i < 11 )  vertical_gradient_level_ind(i) = k + 1
          ENDDO
          initial_profile(k) = scalar_value(i) - ( scalar_value(i) - scalar_value(i+1) ) *         &
                                                 ( scalar_level(i) - zu(k) ) /                     &
                                                 ( scalar_level(i) - scalar_level(i+1) )
       ENDDO
       initial_profile(nzb) = scalar_value(imax)

!
!--    Store bottom/top levels that may have been additionally added above, and determine grid
!--    indices at (or below) gradients are applied. These quantities (gradient/gradient_level/
!--    gradient_level_ind) are only used for output purpose in the _rc file.
       vertical_gradient(:)       = scalar_gradient(:)
       vertical_gradient_level(:) = scalar_level(:)
       DO  i = 1, 12
          IF ( vertical_gradient_level(i) == HUGE( 1.0_wp ) )  EXIT
          DO  k = nzt, nzb+1, -1
             IF ( vertical_gradient_level(i) <= zu(k) )  THEN
                vertical_gradient_level_ind(i) = k
             ELSE
                EXIT
             ENDIF
          ENDDO
       ENDDO

    ENDIF
!
!-- Avoid negative values of scalars.
    IF ( TRIM( quantity ) /= 'ug'  .AND.  TRIM( quantity ) /= 'vg'  .AND.                          &
         TRIM( quantity ) /= 'ws' )                                                                &
    THEN
       DO  k = nzb, nzt+1
          IF ( initial_profile(k) < 0.0_wp )  initial_profile(k) = 0.0_wp
       ENDDO
    ENDIF

 END SUBROUTINE init_vertical_profiles
