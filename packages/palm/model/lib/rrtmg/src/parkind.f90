      module parkind

      use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

      implicit none
      save

!------------------------------------------------------------------
! rrtmg kinds
! Define integer and real kinds for various types.
!
! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!
!     integer kinds
!     -------------
!
      integer, parameter :: kind_ib = int64       ! 8 byte integer
      integer, parameter :: kind_im = int32       ! 4 byte integer
      integer, parameter :: kind_in = kind(1)     ! native integer

!
!     real kinds
!     ----------
!
#if defined( __single_precision )
      integer, parameter :: kind_rb = real32      ! 4 byte real
#else
      integer, parameter :: kind_rb = real64      ! 8 byte real
#endif
      integer, parameter :: kind_rm = real32      ! 4 byte real
      integer, parameter :: kind_rn = kind(1.0)   ! native real

      end module parkind
