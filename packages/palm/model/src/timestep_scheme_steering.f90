!> @file timestep_scheme_steering.f90
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
!> Depending on the timestep scheme set the steering factors for the prognostic equations.
!> tsc(1): switch for adding the explicit vertical diffusion tendency to the total tendency
!>         in intermediate time steps
!> tsc(2), tsc(3): factors for Runge-Kutta substep tendencies
!> tsc(4): factor for implicit vertical diffusion tendency
!> tsc(5): switch for adding damping tendencies, or not
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE timestep_scheme_steering

    USE control_parameters,                                                                        &
        ONLY:  intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               implicit_diffusion_1d,                                                              &
               timestep_scheme,                                                                    &
               tsc

#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY: enable_openacc
#endif

    USE kinds

    IMPLICIT NONE


    IF ( timestep_scheme(1:5) == 'runge' )  THEN
!
!--    Runge-Kutta schemes (here the factors depend on the respective intermediate step)
       IF ( timestep_scheme == 'runge-kutta-2' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             tsc(1:5) = (/ 1.0_wp, 1.0_wp,  0.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             tsc(1:5) = (/ 1.0_wp, 0.5_wp, -0.5_wp, 0.0_wp, 1.0_wp /)
          ENDIF
       ELSE
          IF ( intermediate_timestep_count == 1 )  THEN
             tsc(1:5) = (/ 1.0_wp,  1.0_wp /  3.0_wp,          0.0_wp, 0.0_wp, 0.0_wp /)
          ELSEIF ( intermediate_timestep_count == 2 )  THEN
             tsc(1:5) = (/ 1.0_wp, 15.0_wp / 16.0_wp, -25.0_wp/48.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             tsc(1:5) = (/ 1.0_wp,  8.0_wp / 15.0_wp,   1.0_wp/15.0_wp, 0.0_wp, 1.0_wp /)
          ENDIF
       ENDIF
!
!--    For implicit diffusion, the diffusion tendency is only added for the last intermediate step.
!--    For explicit diffusion, the factor is the same as for the other tendency terms.
       IF ( implicit_diffusion_1d )  THEN
          tsc(1) = 0.0_wp
          IF ( intermediate_timestep_count == intermediate_timestep_count_max )  tsc(4) = 1.0_wp
       ELSE
          tsc(4) = tsc(2)
       ENDIF

       !$ACC UPDATE DEVICE(tsc(1:5)) IF(enable_openacc)

    ELSEIF ( timestep_scheme == 'euler' )  THEN
!
!--    Euler scheme
       tsc(1:5) = (/ 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)

    ENDIF


 END SUBROUTINE timestep_scheme_steering
