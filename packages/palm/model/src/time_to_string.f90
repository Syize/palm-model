!> @file time_to_string.f90
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
!> Transforming the time from real to character-string hh:mm:ss.
!--------------------------------------------------------------------------------------------------!
 FUNCTION time_to_string( time )

    USE kinds

    IMPLICIT NONE

    CHARACTER(LEN=10) ::  time_to_string  !<

    INTEGER(iwp) ::  hours    !< number of hours
    INTEGER(iwp) ::  minutes  !< number of minutes
    INTEGER(iwp) ::  seconds  !< number of seconds

    REAL(wp) ::  rest_time  !< fraction of second

    REAL(wp), INTENT(IN) ::  time  !< model time (usually time_since_reference_point)


!
!-- Calculate the number of hours, minutes, and seconds. Use magnitude because negative times
!-- occur in spinup-runs with coupling_start_time > 0.
    hours     = INT( ABS( time ) / 3600.0_wp )
    rest_time = ABS( time ) - hours * 3600_wp
    minutes   = INT( rest_time / 60.0_wp )
    seconds   = rest_time - minutes * 60

!
!-- Build the string.
    IF ( time < 0.0 )  THEN
       WRITE (time_to_string,'(''-'',I3.3,'':'',I2.2,'':'',I2.2)')  hours, minutes, seconds
    ELSE
       WRITE (time_to_string,'(1X,I3.3,'':'',I2.2,'':'',I2.2)')  hours, minutes, seconds
    ENDIF

!    IF ( hours < 100 )  THEN
!       WRITE (time_to_string,'(I2.2,'':'',I2.2,'':'',I2.2)')  hours, minutes, seconds
!    ELSE
!       WRITE (time_to_string,'(I3.3,'':'',I2.2,'':'',I2.2)')  hours, minutes, seconds
!    ENDIF

!
!-- Negative times appear in case of coupling_start_time > 0.
    IF ( time < 0 )  time_to_string(1:1) = '-'

 END FUNCTION time_to_string
