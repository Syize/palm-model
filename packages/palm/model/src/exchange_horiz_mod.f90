!> @file exchange_horiz.f90
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
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting of cyclic lateral
!> boundary conditions for the total domain .
!--------------------------------------------------------------------------------------------------!
 MODULE exchange_horiz_mod

#if defined( __parallel )
    USE MPI
#endif

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTERFACE exchange_horiz
       MODULE PROCEDURE exchange_horiz
    END INTERFACE exchange_horiz

#if defined( __parallel )
    INTERFACE exchange_horiz_rb
       MODULE PROCEDURE exchange_horiz_rb
    END INTERFACE exchange_horiz_rb
#endif

    INTERFACE exchange_horiz_int
       MODULE PROCEDURE exchange_horiz_int
    END INTERFACE exchange_horiz_int

    INTERFACE exchange_horiz_2d
       MODULE PROCEDURE exchange_horiz_2d
    END INTERFACE exchange_horiz_2d

    INTERFACE exchange_horiz_2d_byte
       MODULE PROCEDURE exchange_horiz_2d_byte
    END INTERFACE exchange_horiz_2d_byte

    INTERFACE exchange_horiz_2d_int
       MODULE PROCEDURE exchange_horiz_2d_int
    END INTERFACE exchange_horiz_2d_int

    PRIVATE

    PUBLIC exchange_horiz,                                                                         &
           exchange_horiz_int,                                                                     &
           exchange_horiz_2d,                                                                      &
           exchange_horiz_2d_byte,                                                                 &
           exchange_horiz_2d_int

#if defined( __parallel )
    PUBLIC exchange_horiz_rb
#endif


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting of cyclic lateral
!> boundary conditions for the total domain.
!> This routine is for REAL 3d-arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz( ar, nbgp_local, alternative_communicator, grid_level, mg_switch_to_pe0 )

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  synchronous_exchange,                                                               &
               use_contiguous_buffer
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt

#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY:  enable_openacc,                                                                     &
               message_string

    USE openacc,                                                                                   &
        ONLY:  acc_is_present
#endif

    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  alternative_communicator  !< alternative MPI communicator to be used
    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  grid_level                !< multigrid grid level to be used

#if defined( __parallel )
    INTEGER(iwp) ::  bufsize       !< size of buffer for sending/receiving contiguous data
#endif
    INTEGER(iwp) ::  communicator  !< communicator that is used as argument in MPI calls
    INTEGER(iwp) ::  i             !< loop index
    INTEGER(iwp) ::  j             !< loop index
    INTEGER(iwp) ::  k             !< loop index
    INTEGER(iwp) ::  l             !< grid level to choose the dreived datatypes
    INTEGER(iwp) ::  left_pe       !< id of left pe that is used as argument in MPI calls
    INTEGER(iwp) ::  nbgp_local    !< number of ghost point layers
    INTEGER(iwp) ::  north_pe      !< id of north pe that is used as argument in MPI calls
    INTEGER(iwp) ::  right_pe      !< id of right pe that is used as argument in MPI calls
    INTEGER(iwp) ::  south_pe      !< id of south pe that is used as argument in MPI calls

    LOGICAL, OPTIONAL, INTENT(IN) ::  mg_switch_to_pe0

    LOGICAL ::  switch_to_pe0  !< local switch telling if total domain is gathered on one PE


    REAL(wp), DIMENSION(nzb:nzt+1,nys-nbgp_local:nyn+nbgp_local,                                   &
                        nxl-nbgp_local:nxr+nbgp_local) ::  ar !< 3d-array for which exchange is done

#if defined( __parallel )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  recvbuf1  !< buffer for sending contiguous ghost layer data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  recvbuf2  !< buffer for sending contiguous ghost layer data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  sendbuf1  !< buffer for sending contiguous ghost layer data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  sendbuf2  !< buffer for sending contiguous ghost layer data
#endif


    CALL cpu_log( log_point_s(2), 'exchange_horiz', 'start' )

#if defined( _OPENACC )
!
!-- Security check for arrays from modules that are still not ported.
    IF ( enable_openacc  .AND.  .NOT. acc_is_present( ar ) )  THEN
       message_string = 'array not on device'
       CALL message( 'exchange_horiz', 'PAC0358', 1, 2, 0, 6, 0 )
    ENDIF
#endif

!
!-- Set the communicator to be used
    IF ( PRESENT( alternative_communicator ) )  THEN
!
!--    Alternative communicator is to be used
       communicator = communicator_configurations(alternative_communicator)%mpi_communicator
       left_pe  = communicator_configurations(alternative_communicator)%pleft
       right_pe = communicator_configurations(alternative_communicator)%pright
       south_pe = communicator_configurations(alternative_communicator)%psouth
       north_pe = communicator_configurations(alternative_communicator)%pnorth

    ELSE
!
!--    Main communicator is to be used
       communicator = comm2d
       left_pe  = pleft
       right_pe = pright
       south_pe = psouth
       north_pe = pnorth

    ENDIF

!
!-- Set the grid level.
    IF ( PRESENT( grid_level ) )  THEN
       l = grid_level
    ELSE
       l = 0
    ENDIF

!
!-- Set the switch that tells if data of the toal domain is gathered on the PE
    IF ( PRESENT( mg_switch_to_pe0 ) )  THEN
       switch_to_pe0 = mg_switch_to_pe0
    ELSE
       switch_to_pe0 = .FALSE.
    ENDIF

#if defined( __parallel )
!
!-- Exchange in y-direction of lateral boundaries. Must be done before exchange in x, because
!-- in case of buffering, only computational grid points (nxl:nxr) will be exchanged.
    IF ( npey == 1  .OR.  switch_to_pe0 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged within the PE memory
       IF ( PRESENT( alternative_communicator ) )  THEN
          IF ( alternative_communicator == 1  .OR.  alternative_communicator == 3 )  THEN
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  k = nzb, nzt+1
                DO  j = 1, nbgp_local
                   DO  i = nxl-nbgp_local, nxr+nbgp_local
                      ar(k,nys-nbgp_local-1+j,i) = ar(k,nyn-nbgp_local+j,i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  k = nzb, nzt+1
                DO  j = 1, nbgp_local
                   DO  i = nxl-nbgp_local, nxr+nbgp_local
                      ar(k,nyn+j,i) = ar(k,nys+j-1,i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
          ENDIF
       ELSE
          IF ( bc_ns_cyc )  THEN
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  k = nzb, nzt+1
                DO  j = 1, nbgp_local
                   DO  i = nxl-nbgp_local, nxr+nbgp_local
                      ar(k,nys-nbgp_local-1+j,i) = ar(k,nyn-nbgp_local+j,i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  k = nzb, nzt+1
                DO  j = 1, nbgp_local
                   DO  i = nxl-nbgp_local, nxr+nbgp_local
                      ar(k,nyn+j,i) = ar(k,nys+j-1,i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
          ENDIF
       ENDIF

    ELSE

       IF ( synchronous_exchange )  THEN
!
!--       Send front boundary, receive rear one (synchronous)
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_SENDRECV( ar(nzb,nys,nxl-nbgp_local),   1, type_xz(l), south_pe, 0,             &
                             ar(nzb,nyn+1,nxl-nbgp_local), 1, type_xz(l), north_pe, 0,             &
                             communicator, status, ierr )
          !$ACC END HOST_DATA
!
!--       Send rear boundary, receive front one (synchronous)
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_SENDRECV( ar(nzb,nyn-nbgp_local+1,nxl-nbgp_local), 1, type_xz(l), north_pe, 1,  &
                             ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1, type_xz(l), south_pe, 1,  &
                             communicator, status, ierr )
          !$ACC END HOST_DATA

       ELSE

!
!--       Asynchroneous exchange
          req(1:4)  = 0
          req_count = 0

          IF ( use_contiguous_buffer )  THEN

!
!--          Allocate buffers for sending / receiving ghost layer data
             ALLOCATE( recvbuf1(nzb:nzt+1,nxl:nxr,nbgp_local) )
             ALLOCATE( recvbuf2(nzb:nzt+1,nxl:nxr,nbgp_local) )
             ALLOCATE( sendbuf1(nzb:nzt+1,nxl:nxr,nbgp_local) )
             ALLOCATE( sendbuf2(nzb:nzt+1,nxl:nxr,nbgp_local) )
             !$ACC DATA CREATE(recvbuf1,recvbuf2,sendbuf1,sendbuf2) IF(enable_openacc)

             bufsize = ( nzt - nzb + 2 ) * ( nxr - nxl + 1 ) * nbgp_local

!
!--          Pack data for sending south boundary.
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl, nxr
                DO  j = 1, nbgp_local
                   DO  k = nzb, nzt+1
                      sendbuf1(k,i,j) = ar(k,nys+j-1,i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
!
!--          Pack data for sending north boundary.
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl, nxr
                DO  j = 1, nbgp_local
                   DO  k = nzb, nzt+1
                      sendbuf2(k,i,j) = ar(k,nyn-nbgp_local+j,i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
!
!--          Send south boundary, receive north one (asynchronous)
             !$ACC HOST_DATA USE_DEVICE(sendbuf1) IF(enable_openacc)
             CALL MPI_ISEND( sendbuf1(nzb,nxl,1), bufsize, MPI_REAL, south_pe, req_count,          &
                             communicator, req(req_count+1), ierr )
             !$ACC END HOST_DATA

             !$ACC HOST_DATA USE_DEVICE(recvbuf1) IF(enable_openacc)
             CALL MPI_IRECV( recvbuf1(nzb,nxl,1), bufsize, MPI_REAL, north_pe, req_count,          &
                             communicator, req(req_count+2), ierr )
             !$ACC END HOST_DATA
!
!--          Send north boundary, receive south one (asynchronous)
             !$ACC HOST_DATA USE_DEVICE(sendbuf2) IF(enable_openacc)
             CALL MPI_ISEND( sendbuf2(nzb,nxl,1), bufsize, MPI_REAL, north_pe, req_count+1,        &
                             communicator, req(req_count+3), ierr )
             !$ACC END HOST_DATA

             !$ACC HOST_DATA USE_DEVICE(recvbuf2) IF(enable_openacc)
             CALL MPI_IRECV( recvbuf2(nzb,nxl,1), bufsize, MPI_REAL, south_pe, req_count+1,        &
                             communicator, req(req_count+4), ierr )
             !$ACC END HOST_DATA

             CALL MPI_WAITALL( 4, req, wait_stat, ierr )
!
!--          Unpack received north boundary data
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl, nxr
                DO  j = 1, nbgp_local
                   DO  k = nzb, nzt+1
                      ar(k,nyn+j,i) = recvbuf1(k,i,j)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
!
!--          Unpack received south boundary data
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = nxl, nxr
                DO  j = 1, nbgp_local
                   DO  k = nzb, nzt+1
                      ar(k,nys-nbgp_local-1+j,i) = recvbuf2(k,i,j)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP

             !$ACC END DATA
             DEALLOCATE( recvbuf1, recvbuf2, sendbuf1, sendbuf2 )

          ELSE
!
!--          Send front boundary, receive rear one (asynchronous)
             !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
             CALL MPI_ISEND( ar(nzb,nys,nxl-nbgp_local),   1, type_xz(l), south_pe,                &
                             req_count, communicator, req(req_count+1), ierr )
             !$ACC END HOST_DATA
             !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
             CALL MPI_IRECV( ar(nzb,nyn+1,nxl-nbgp_local), 1, type_xz(l), north_pe,                &
                             req_count, communicator, req(req_count+2), ierr )
             !$ACC END HOST_DATA
!
!--          Send rear boundary, receive front one (asynchronous)
             !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
             CALL MPI_ISEND( ar(nzb,nyn-nbgp_local+1,nxl-nbgp_local), 1, type_xz(l), north_pe,     &
                             req_count+1, communicator, req(req_count+3), ierr )
             !$ACC END HOST_DATA
             !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1, type_xz(l), south_pe,     &
                             req_count+1, communicator, req(req_count+4), ierr )
             !$ACC END HOST_DATA

             CALL MPI_WAITALL( 4, req, wait_stat, ierr )

          ENDIF   ! contiguous buffer or MPI derived type

       ENDIF   ! synchronous or asynchronous

    ENDIF   ! npey == 1  / switch_to_pe0 or not

!
!-- Exchange in x-direction of lateral boundaries. All grid points including ghost points
!-- (nys-nbgp_local:nyn+nbgp_local) are exchanged.
    IF ( npex == 1  .OR.  switch_to_pe0 )  THEN
!
!--    One-dimensional decomposition along y, boundary values can be exchanged within the PE memory.
       IF ( PRESENT( alternative_communicator ) )  THEN
          IF ( alternative_communicator <= 2 )  THEN
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = 1, nbgp_local
                DO  j = nys-nbgp_local, nyn+nbgp_local
                   DO  k = nzb, nzt+1
                      ar(k,j,nxl-nbgp_local-1+i) = ar(k,j,nxr-nbgp_local+i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = 1, nbgp_local
                DO  j = nys-nbgp_local, nyn+nbgp_local
                   DO  k = nzb, nzt+1
                      ar(k,j,nxr+i) = ar(k,j,nxl+i-1)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
          ENDIF
       ELSE
          IF ( bc_lr_cyc )  THEN
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = 1, nbgp_local
                DO  j = nys-nbgp_local, nyn+nbgp_local
                   DO  k = nzb, nzt+1
                      ar(k,j,nxl-nbgp_local-1+i) = ar(k,j,nxr-nbgp_local+i)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
             !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
             DO  i = 1, nbgp_local
                DO  j = nys-nbgp_local, nyn+nbgp_local
                   DO  k = nzb, nzt+1
                      ar(k,j,nxr+i) = ar(k,j,nxl+i-1)
                   ENDDO
                ENDDO
             ENDDO
             !$ACC END PARALLEL LOOP
          ENDIF
       ENDIF

    ELSE

       IF ( synchronous_exchange )  THEN
!
!--       Send left boundary, receive right one (synchronous)
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_SENDRECV( ar(nzb,nys-nbgp_local,nxl),   1, type_yz(l), left_pe,  0,             &
                             ar(nzb,nys-nbgp_local,nxr+1), 1, type_yz(l), right_pe, 0,             &
                             communicator, status, ierr )
          !$ACC END HOST_DATA
!
!--       Send right boundary, receive left one (synchronous)
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_SENDRECV( ar(nzb,nys-nbgp_local,nxr+1-nbgp_local), 1, type_yz(l), right_pe, 1,  &
                             ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1, type_yz(l), left_pe,  1,  &
                             communicator, status, ierr )
          !$ACC END HOST_DATA

       ELSE

!
!--       Asynchroneous exchange
          req(1:4)  = 0
          req_count = 0

!
!--       Send left boundary, receive right one (asynchronous)
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_ISEND( ar(nzb,nys-nbgp_local,nxl),   1, type_yz(l), left_pe,                    &
                          req_count, communicator, req(req_count+1), ierr )
          !$ACC END HOST_DATA
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxr+1), 1, type_yz(l), right_pe,                   &
                          req_count, communicator, req(req_count+2), ierr )
          !$ACC END HOST_DATA
!
!--       Send right boundary, receive left one (asynchronous)
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_ISEND( ar(nzb,nys-nbgp_local,nxr+1-nbgp_local), 1, type_yz(l), right_pe,        &
                          req_count+1, communicator, req(req_count+3), ierr )
          !$ACC END HOST_DATA
          !$ACC HOST_DATA USE_DEVICE(ar) IF(enable_openacc)
          CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1, type_yz(l), left_pe,         &
                          req_count+1, communicator, req(req_count+4), ierr )
          !$ACC END HOST_DATA

          CALL MPI_WAITALL( 4, req, wait_stat, ierr )

       ENDIF   ! synchronous or asynchronous

    ENDIF   ! npex == 1  / switch_to_pe0 or not


#else

!
!-- Lateral boundary conditions in the non-parallel case.
!-- Case dependent, because in GPU mode still not all arrays are on device. This workaround has to
!-- be removed later. Also, since PGI compiler 12.5 has problems with array syntax, explicit loops
!-- are used.
    IF ( PRESENT( alternative_communicator ) )  THEN
       IF ( alternative_communicator <= 2 )  THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = 1, nbgp_local
             DO  j = nys-nbgp_local, nyn+nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,j,nxl-nbgp_local-1+i) = ar(k,j,nxr-nbgp_local+i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = 1, nbgp_local
             DO  j = nys-nbgp_local, nyn+nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,j,nxr+i) = ar(k,j,nxl-1+i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF
    ELSE
       IF ( bc_lr_cyc )  THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = 1, nbgp_local
             DO  j = nys-nbgp_local, nyn+nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,j,nxl-nbgp_local-1+i) = ar(k,j,nxr-nbgp_local+i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = 1, nbgp_local
             DO  j = nys-nbgp_local, nyn+nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,j,nxr+i) = ar(k,j,nxl-1+i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF
    ENDIF

    IF ( PRESENT( alternative_communicator ) )  THEN
       IF ( alternative_communicator == 1  .OR.  alternative_communicator == 3 )  THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = nxl-nbgp_local, nxr+nbgp_local
             DO  j = 1, nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,nys-nbgp_local-1+j,i) = ar(k,nyn-nbgp_local+j,i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = nxl-nbgp_local, nxr+nbgp_local
             DO  j = 1, nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,nyn+j,i) = ar(k,nys-1+j,i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF
    ELSE
       IF ( bc_ns_cyc )  THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = nxl-nbgp_local, nxr+nbgp_local
             DO  j = 1, nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,nys-nbgp_local-1+j,i) = ar(k,nyn-nbgp_local+j,i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lacc) IF(enable_openacc)
          DO  i = nxl-nbgp_local, nxr+nbgp_local
             DO  j = 1, nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,nyn+j,i) = ar(k,nys-1+j,i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF
    ENDIF

#endif

    CALL cpu_log( log_point_s(2), 'exchange_horiz', 'stop' )

 END SUBROUTINE exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting of cyclic lateral
!> boundary conditions for the total domain.
!> This routine for the special exchange of red/black grid points used in the multigrid solver.
!> Only the respective red or black points will be exchanged. For the exchange, data is always
!> stored in contiguous buffer arrays.
!> This routine can be (and is) used in parallel mode only.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE exchange_horiz_rb( ar, nbgp_local, color, kinc, ileft_for_nys_send,                    &
                               ileft_for_nys_recv, ileft_for_nyn_send, ileft_for_nyn_recv,         &
                               jsouth_for_nxl_send, jsouth_for_nxl_recv, jsouth_for_nxr_send,      &
                               jsouth_for_nxr_recv, kbottom_for_nys_send, kbottom_for_nys_recv,    &
                               kbottom_for_nyn_send, kbottom_for_nyn_recv, kbottom_for_nxl_send,   &
                               kbottom_for_nxl_recv, kbottom_for_nxr_send, kbottom_for_nxr_recv,   &
                               ktop_for_nys_send, ktop_for_nys_recv, ktop_for_nyn_send,            &
                               ktop_for_nyn_recv, ktop_for_nxl_send, ktop_for_nxl_recv,            &
                               ktop_for_nxr_send, ktop_for_nxr_recv )

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s

#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY:  enable_openacc,                                                                     &
               message_string

    USE openacc,                                                                                   &
        ONLY:  acc_is_present
#endif

    INTEGER(iwp) ::  color  !< red or black colored pressure points in multigrid solver
    INTEGER(iwp) ::  kinc   !< increment of k loop (= 2, for even/odd sorting = 1)

    INTEGER(iwp), DIMENSION(2,2) ::  ileft_for_nyn_recv    !< lower loop index for i loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ileft_for_nyn_send    !< lower loop index for i loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ileft_for_nys_recv    !< lower loop index for i loop (ghost point exchange of right boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ileft_for_nys_send    !< lower loop index for i loop (ghost point exchange of right boundary)

    INTEGER(iwp), DIMENSION(2,2) ::  jsouth_for_nxl_recv   !< lower loop index for j loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  jsouth_for_nxl_send   !< lower loop index for j loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  jsouth_for_nxr_recv   !< lower loop index for j loop (ghost point exchange of right boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  jsouth_for_nxr_send   !< lower loop index for j loop (ghost point exchange of right boundary)

    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nxl_recv  !< lower loop index for k loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nxl_send  !< lower loop index for k loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nxr_recv  !< lower loop index for k loop (ghost point exchange of right boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nxr_send  !< lower loop index for k loop (ghost point exchange of right boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nyn_recv  !< lower loop index for k loop (ghost point exchange of north boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nyn_send  !< lower loop index for k loop (ghost point exchange of north boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nys_recv  !< lower loop index for k loop (ghost point exchange of south boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  kbottom_for_nys_send  !< lower loop index for k loop (ghost point exchange of south boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nxl_recv     !< upper loop index for k loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nxl_send     !< upper loop index for k loop (ghost point exchange of left boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nxr_recv     !< upper loop index for k loop (ghost point exchange of right boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nxr_send     !< upper loop index for k loop (ghost point exchange of right boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nyn_recv     !< upper loop index for k loop (ghost point exchange of north boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nyn_send     !< upper loop index for k loop (ghost point exchange of north boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nys_recv     !< upper loop index for k loop (ghost point exchange of south boundary)
    INTEGER(iwp), DIMENSION(2,2) ::  ktop_for_nys_send     !< upper loop index for k loop (ghost point exchange of south boundary)

    INTEGER(iwp) ::  bufsize       !< size of buffer for sending/receiving contiguous data
    INTEGER(iwp) ::  communicator  !< communicator that is used as argument in MPI calls
    INTEGER(iwp) ::  i             !< loop index
    INTEGER(iwp) ::  ii            !< loop index
    INTEGER(iwp) ::  il            !< loop index
    INTEGER(iwp) ::  j             !< loop index
    INTEGER(iwp) ::  jj            !< loop index
    INTEGER(iwp) ::  js            !< loop index
    INTEGER(iwp) ::  k             !< loop index
    INTEGER(iwp) ::  kb            !< bottom loop index
    INTEGER(iwp) ::  ki            !< local increment of k loop
    INTEGER(iwp) ::  kt            !< top loop index
    INTEGER(iwp) ::  nxr_rb        !< upper index for red/black decomposition (pressure)
    INTEGER(iwp) ::  nyn_rb        !< upper index for red/black decomposition (pressure)
    INTEGER(iwp) ::  left_pe       !< id of left pe that is used as argument in MPI calls
    INTEGER(iwp) ::  nbgp_local    !< number of ghost point layers
    INTEGER(iwp) ::  north_pe      !< id of north pe that is used as argument in MPI calls
    INTEGER(iwp) ::  right_pe      !< id of right pe that is used as argument in MPI calls
    INTEGER(iwp) ::  south_pe      !< id of south pe that is used as argument in MPI calls

    REAL(wp), DIMENSION(nzb:nzt+1,nys-nbgp_local:nyn+nbgp_local,                                   &
                        nxl-nbgp_local:nxr+nbgp_local) ::  ar !< 3d-array for which exchange is done

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  recvbuf1  !< buffer for sending contiguous ghost layer data
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  recvbuf2  !< buffer for sending contiguous ghost layer data
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sendbuf1  !< buffer for sending contiguous ghost layer data
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sendbuf2  !< buffer for sending contiguous ghost layer data


    CALL cpu_log( log_point_s(47), 'exchange_horiz_rb', 'start' )

#if defined( _OPENACC )
!
!-- Security check for arrays from modules that are still not ported.
    IF ( enable_openacc  .AND.  .NOT. acc_is_present( ar ) )  THEN
       message_string = 'array not on device'
       CALL message( 'exchange_horiz_rb', 'PAC0358', 1, 2, 0, 6, 0 )
    ENDIF
#endif

!
!-- Main communicator is to be used.
    communicator = comm2d
    left_pe  = pleft
    right_pe = pright
    south_pe = psouth
    north_pe = pnorth

!
!-- Exchange in y-direction of lateral boundaries. Must be done before exchange in x, because
!-- in case of buffering, only computational grid points (nxl:nxr) will be exchanged.
    IF ( npey == 1 )  THEN

       IF ( bc_ns_cyc )  THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
          DO  k = nzb, nzt+1
             DO  j = 1, nbgp_local
                DO  i = nxl-nbgp_local, nxr+nbgp_local
                   ar(k,nys-nbgp_local-1+j,i) = ar(k,nyn-nbgp_local+j,i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
          DO  k = nzb, nzt+1
             DO  j = 1, nbgp_local
                DO  i = nxl-nbgp_local, nxr+nbgp_local
                   ar(k,nyn+j,i) = ar(k,nys+j-1,i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF

    ELSE

       req(1:4)  = 0
       req_count = 0

       ki = kinc
!
!--    Special case for exchanging only red/black colored pressure points.
!--    Allocate buffers for sending / receiving ghost layer data.
!--    Only 2d-arrays are required since number of ghost layers is limited to one.
!--    Third index is used here, since the arrays are declared as 3d.
       nxr_rb = ( nxr - nxl ) / 2

       ALLOCATE( recvbuf1(nzb+1:nzt,0:nxr_rb) )
       ALLOCATE( recvbuf2(nzb+1:nzt,0:nxr_rb) )
       ALLOCATE( sendbuf1(nzb+1:nzt,0:nxr_rb) )
       ALLOCATE( sendbuf2(nzb+1:nzt,0:nxr_rb) )
       !$ACC DATA CREATE(recvbuf1,recvbuf2,sendbuf1,sendbuf2) IF(enable_openacc)

       bufsize = ( nzt - nzb ) * ( nxr_rb + 1 )

!
!--    Pack data for sending south/north boundary.
       ii = -1
       il = ileft_for_nys_send(1,color)
       kb = kbottom_for_nys_send(1,color)
       kt = ktop_for_nys_send(1,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  i = il, nxr, 2
          ii = ii + 1
          DO  k = kb, kt, ki
             sendbuf1(k,ii) = ar(k,nys,i)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

       ii = -1
       il = ileft_for_nyn_send(1,color)
       kb = kbottom_for_nyn_send(1,color)
       kt = ktop_for_nyn_send(1,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  i = il, nxr, 2
          ii = ii + 1
          DO  k = kb, kt, ki
             sendbuf2(k,ii) = ar(k,nyn,i)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

       ii = -1
       il = ileft_for_nys_send(2,color)
       kb = kbottom_for_nys_send(2,color)
       kt = ktop_for_nys_send(2,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  i = il, nxr, 2
          ii = ii + 1
          DO  k = kb, kt, ki
             sendbuf1(k,ii) = ar(k,nys,i)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

       ii = -1
       il = ileft_for_nyn_send(2,color)
       kb = kbottom_for_nyn_send(2,color)
       kt = ktop_for_nyn_send(2,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  i = il, nxr, 2
          ii = ii + 1
          DO  k = kb, kt, ki
             sendbuf2(k,ii) = ar(k,nyn,i)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

!
!--    Send south boundary, receive north one (asynchronous), only red or black points.
       !$ACC HOST_DATA USE_DEVICE(sendbuf1) IF(enable_openacc)
       CALL MPI_ISEND( sendbuf1(nzb+1,0), bufsize, MPI_REAL, south_pe, req_count, communicator,    &
                       req(req_count+1), ierr )
       !$ACC END HOST_DATA

       !$ACC HOST_DATA USE_DEVICE(recvbuf1) IF(enable_openacc)
       CALL MPI_IRECV( recvbuf1(nzb+1,0), bufsize, MPI_REAL, north_pe, req_count, communicator,    &
                       req(req_count+2), ierr )
       !$ACC END HOST_DATA
!
!--    Send north boundary, receive south one (asynchronous)
       !$ACC HOST_DATA USE_DEVICE(sendbuf2) IF(enable_openacc)
       CALL MPI_ISEND( sendbuf2(nzb+1,0), bufsize, MPI_REAL, north_pe, req_count+1, communicator,  &
                       req(req_count+3), ierr )
       !$ACC END HOST_DATA

       !$ACC HOST_DATA USE_DEVICE(recvbuf2) IF(enable_openacc)
       CALL MPI_IRECV( recvbuf2(nzb+1,0), bufsize, MPI_REAL, south_pe, req_count+1, communicator,  &
                       req(req_count+4), ierr )
       !$ACC END HOST_DATA

       CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!--    Unpack received south/north boundary data.
       IF ( .NOT. bc_dirichlet_n  .AND.  .NOT. bc_radiation_n )  THEN

          ii = -1
          il = ileft_for_nyn_recv(1,color)
          kb = kbottom_for_nyn_recv(1,color)
          kt = ktop_for_nyn_recv(1,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  i = il, nxr, 2
             ii = ii + 1
             DO  k = kb, kt, ki
                ar(k,nyn+1,i) = recvbuf1(k,ii)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

          ii = -1
          il = ileft_for_nyn_recv(2,color)
          kb = kbottom_for_nyn_recv(2,color)
          kt = ktop_for_nyn_recv(2,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  i = il, nxr, 2
             ii = ii + 1
             DO  k = kb, kt, ki
                ar(k,nyn+1,i) = recvbuf1(k,ii)
                ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

       ENDIF

       IF ( .NOT. bc_dirichlet_s  .AND.  .NOT. bc_radiation_s )  THEN

          ii = -1
          il = ileft_for_nys_recv(1,color)
          kb = kbottom_for_nys_recv(1,color)
          kt = ktop_for_nys_recv(1,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  i = il, nxr, 2
             ii = ii + 1
             DO  k = kb, kt, ki
                ar(k,nys-1,i) = recvbuf2(k,ii)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

          ii = -1
          il = ileft_for_nys_recv(2,color)
          kb = kbottom_for_nys_recv(2,color)
          kt = ktop_for_nys_recv(2,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  i = il, nxr, 2
             ii = ii + 1
             DO  k = kb, kt, ki
                ar(k,nys-1,i) = recvbuf2(k,ii)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

       ENDIF

       !$ACC END DATA
       DEALLOCATE( recvbuf1, recvbuf2, sendbuf1, sendbuf2 )

    ENDIF

!
!-- Exchange in x-direction of lateral boundaries. All grid points including ghost points
!-- (nys-nbgp_local:nyn+nbgp_local) are exchanged.
    IF ( npex == 1 )  THEN

       IF ( bc_lr_cyc )  THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
          DO  i = 1, nbgp_local
             DO  j = nys-nbgp_local, nyn+nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,j,nxl-nbgp_local-1+i) = ar(k,j,nxr-nbgp_local+i)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(enable_openacc)
          DO  i = 1, nbgp_local
             DO  j = nys-nbgp_local, nyn+nbgp_local
                DO  k = nzb, nzt+1
                   ar(k,j,nxr+i) = ar(k,j,nxl+i-1)
                ENDDO
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
       ENDIF

    ELSE

       req(1:4)  = 0
       req_count = 0

       ki = kinc
!
!--    Special case for exchanging red/black colored pressure points.
!--    Always using contiguous buffers!
!--    Allocate buffers for sending / receiving ghost layer data.
       nyn_rb = ( nyn - nys ) / 2

       ALLOCATE( recvbuf1(nzb+1:nzt,-1:nyn_rb+1) )
       ALLOCATE( recvbuf2(nzb+1:nzt,-1:nyn_rb+1) )
       ALLOCATE( sendbuf1(nzb+1:nzt,-1:nyn_rb+1) )
       ALLOCATE( sendbuf2(nzb+1:nzt,-1:nyn_rb+1) )
       !$ACC DATA CREATE(recvbuf1,recvbuf2,sendbuf1,sendbuf2) IF(enable_openacc)

       bufsize = ( nzt - nzb ) * ( nyn_rb + 3 )

!
!--    Pack data for sending left/right boundary.
       jj = -2
       js = jsouth_for_nxl_send(1,color)
       kb = kbottom_for_nxl_send(1,color)
       kt = ktop_for_nxl_send(1,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  j = js, nyn+1, 2
          jj = jj + 1
          DO  k = kb, kt, ki
             sendbuf1(k,jj) = ar(k,j,nxl)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

       jj = -2
       js = jsouth_for_nxr_send(1,color)
       kb = kbottom_for_nxr_send(1,color)
       kt = ktop_for_nxr_send(1,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  j = js, nyn+1, 2
          jj = jj + 1
          DO  k = kb, kt, ki
             sendbuf2(k,jj) = ar(k,j,nxr)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

       jj = -2
       js = jsouth_for_nxl_send(2,color)
       kb = kbottom_for_nxl_send(2,color)
       kt = ktop_for_nxl_send(2,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  j = js, nyn+1, 2
          jj = jj + 1
          DO  k = kb, kt, ki
             sendbuf1(k,jj) = ar(k,j,nxl)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

       jj = -2
       js = jsouth_for_nxr_send(2,color)
       kb = kbottom_for_nxr_send(2,color)
       kt = ktop_for_nxr_send(2,color)
       !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
       DO  j = js, nyn+1, 2
          jj = jj + 1
          DO  k = kb, kt, ki
             sendbuf2(k,jj) = ar(k,j,nxr)
          ENDDO
       ENDDO
       !$ACC END PARALLEL LOOP

!
!--    Send left boundary, receive right one (asynchronous), only red or black points
       !$ACC HOST_DATA USE_DEVICE(sendbuf1) IF(enable_openacc)
       CALL MPI_ISEND( sendbuf1(nzb+1,-1), bufsize, MPI_REAL, left_pe, req_count, communicator,    &
                       req(req_count+1), ierr )
       !$ACC END HOST_DATA

       !$ACC HOST_DATA USE_DEVICE(recvbuf1) IF(enable_openacc)
       CALL MPI_IRECV( recvbuf1(nzb+1,-1), bufsize, MPI_REAL, right_pe, req_count, communicator,   &
                       req(req_count+2), ierr )
       !$ACC END HOST_DATA
!
!--    Send right boundary, receive left one (asynchronous)
       !$ACC HOST_DATA USE_DEVICE(sendbuf2) IF(enable_openacc)
       CALL MPI_ISEND( sendbuf2(nzb+1,-1), bufsize, MPI_REAL, right_pe, req_count+1, communicator, &
                       req(req_count+3), ierr )
       !$ACC END HOST_DATA

       !$ACC HOST_DATA USE_DEVICE(recvbuf2) IF(enable_openacc)
       CALL MPI_IRECV( recvbuf2(nzb+1,-1), bufsize, MPI_REAL, left_pe, req_count+1, communicator,  &
                       req(req_count+4), ierr )
       !$ACC END HOST_DATA

       CALL MPI_WAITALL( 4, req, wait_stat, ierr )

!
!--    Unpack received left/right boundary data.
       IF ( .NOT. bc_dirichlet_r  .AND.  .NOT. bc_radiation_r )  THEN

          jj = -2
          js = jsouth_for_nxr_recv(1,color)
          kb = kbottom_for_nxr_recv(1,color)
         kt = ktop_for_nxr_recv(1,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  j = js, nyn+1, 2
             jj = jj + 1
             DO  k = kb, kt, ki
               ar(k,j,nxr+1) = recvbuf1(k,jj)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

          jj = -2
          js = jsouth_for_nxr_recv(2,color)
          kb = kbottom_for_nxr_recv(2,color)
          kt = ktop_for_nxr_recv(2,color)
         !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  j = js, nyn+1, 2
             jj = jj + 1
             DO  k = kb, kt, ki
                ar(k,j,nxr+1) = recvbuf1(k,jj)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

       ENDIF

       IF ( .NOT. bc_dirichlet_l  .AND.  .NOT. bc_radiation_l )  THEN

          jj = -2
          js = jsouth_for_nxl_recv(1,color)
          kb = kbottom_for_nxl_recv(1,color)
          kt = ktop_for_nxl_recv(1,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  j = js, nyn+1, 2
             jj = jj + 1
             DO  k = kb, kt, ki
                ar(k,j,nxl-1) = recvbuf2(k,jj)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

          jj = -2
          js = jsouth_for_nxl_recv(2,color)
          kb = kbottom_for_nxl_recv(2,color)
          kt = ktop_for_nxl_recv(2,color)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(enable_openacc)
          DO  j = js, nyn+1, 2
             jj = jj + 1
             DO  k = kb, kt, ki
                ar(k,j,nxl-1) = recvbuf2(k,jj)
             ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

       ENDIF

       !$ACC END DATA
       DEALLOCATE( recvbuf1, recvbuf2, sendbuf1, sendbuf2 )

    ENDIF   ! npex == 1

    CALL cpu_log( log_point_s(47), 'exchange_horiz_rb', 'stop' )

 END SUBROUTINE exchange_horiz_rb
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_int( ar, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, nbgp_local, type_xz_in,     &
                                type_yz_in, alternative_communicator )

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc

    USE indices,                                                                                   &
        ONLY:  nzb

    INTEGER(iwp), OPTIONAL ::  alternative_communicator  !< alternative MPI communicator to be used

    INTEGER(iwp) ::  communicator  !< communicator that is used as argument in MPI calls
    INTEGER(iwp) ::  left_pe       !< id of left pe that is used as argument in MPI calls
    INTEGER(iwp) ::  nbgp_local    !< number of ghost points
    INTEGER(iwp) ::  north_pe      !< id of north pe that is used as argument in MPI calls
    INTEGER(iwp) ::  nxl_l         !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l         !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l         !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l         !< local index bound at current grid level, south side
    INTEGER(iwp) ::  nzt_l         !< local index bound at current grid level, top
    INTEGER(iwp) ::  right_pe      !< id of right pe that is used as argument in MPI calls
    INTEGER(iwp) ::  south_pe      !< id of south pe that is used as argument in MPI calls
    INTEGER(iwp) ::  type_xz       !< MPI datatype of exchange 3D data - left/right
    INTEGER(iwp) ::  type_yz       !< MPI datatype of exchange 3D data - north/south

    INTEGER(iwp), OPTIONAL ::  type_xz_in !< passed MPI datatype to exchange 3D data between left/right MPI ranks
    INTEGER(iwp), OPTIONAL ::  type_yz_in !< passed MPI datatype to exchange 3D data between north/south MPI ranks

    INTEGER(iwp), DIMENSION(nzb:nzt_l+1,nys_l-nbgp_local:nyn_l+nbgp_local,                         &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array


!
!-- Set MPI datatype depending on the requested task.
    IF ( PRESENT( type_xz_in ) )  THEN
       type_xz = type_xz_in
    ELSE
#if defined( __parallel )
       type_xz = type_xz_int(0)
#endif
    ENDIF

    IF ( PRESENT( type_yz_in ) )  THEN
       type_yz = type_yz_in
    ELSE
#if defined( __parallel )
       type_yz = type_yz_int(0)
#endif
    ENDIF

!
!-- Set the communicator to be used.
    IF ( PRESENT( alternative_communicator ) )  THEN
!
!--    Alternative communicator is to be used.
       communicator = communicator_configurations(alternative_communicator)%mpi_communicator
       left_pe  = communicator_configurations(alternative_communicator)%pleft
       right_pe = communicator_configurations(alternative_communicator)%pright
       south_pe = communicator_configurations(alternative_communicator)%psouth
       north_pe = communicator_configurations(alternative_communicator)%pnorth

    ELSE
!
!--    Main communicator is to be used.
       communicator = comm2d
       left_pe  = pleft
       right_pe = pright
       south_pe = psouth
       north_pe = pnorth

    ENDIF


#if defined( __parallel )
    IF ( npex == 1 )  THEN
!
!--    One-dimensional decomposition along y, boundary values can be exchanged within the PE memory.
       IF ( PRESENT( alternative_communicator ) )  THEN
          IF ( alternative_communicator <= 2 )  THEN
             ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
             ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
          ENDIF
       ELSE
          IF ( bc_lr_cyc )  THEN
             ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
             ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
          ENDIF
       ENDIF
    ELSE
!
!--    Send left boundary, receive right one (synchronous).
       CALL MPI_SENDRECV( ar(nzb,nys_l-nbgp_local,nxl_l),   1, type_yz, left_pe,  0,               &
                          ar(nzb,nys_l-nbgp_local,nxr_l+1), 1, type_yz, right_pe, 0,               &
                          communicator, status, ierr )
!
!--    Send right boundary, receive left one (synchronous).
       CALL MPI_SENDRECV( ar(nzb,nys_l-nbgp_local,nxr_l+1-nbgp_local), 1, type_yz, right_pe, 1,    &
                          ar(nzb,nys_l-nbgp_local,nxl_l-nbgp_local),   1, type_yz, left_pe,  1,    &
                          communicator, status, ierr )
    ENDIF


    IF ( npey == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged within the PE memory.
       IF ( PRESENT( alternative_communicator ) )  THEN
          IF ( alternative_communicator == 1  .OR.  alternative_communicator == 3 )  THEN
             ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
             ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
          ENDIF
       ELSE
          IF ( bc_ns_cyc )  THEN
             ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
             ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
          ENDIF
       ENDIF

    ELSE

!
!--    Send front boundary, receive rear one (synchronous).
       CALL MPI_SENDRECV( ar(nzb,nys_l,nxl_l-nbgp_local),   1, type_xz, south_pe, 0,               &
                          ar(nzb,nyn_l+1,nxl_l-nbgp_local), 1, type_xz, north_pe, 0,               &
                          communicator, status, ierr )
!
!--    Send rear boundary, receive front one (synchronous).
       CALL MPI_SENDRECV( ar(nzb,nyn_l-nbgp_local+1,nxl_l-nbgp_local), 1, type_xz, north_pe, 1,    &
                          ar(nzb,nys_l-nbgp_local,nxl_l-nbgp_local),   1, type_xz, south_pe, 1,    &
                          communicator, status, ierr )

    ENDIF

#else
!
!-- Lateral boundary conditions in the non-parallel case.
    IF ( PRESENT( alternative_communicator ) )  THEN
       IF ( alternative_communicator <= 2 )  THEN
          ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
          ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
       ENDIF
    ELSE
       IF ( bc_lr_cyc )  THEN
          ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
          ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
       ENDIF
    ENDIF

    IF ( PRESENT( alternative_communicator ) )  THEN
       IF ( alternative_communicator == 1  .OR.  alternative_communicator == 3 )  THEN
          ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
          ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
       ENDIF
    ELSE
       IF ( bc_ns_cyc )  THEN
          ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
          ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
       ENDIF
    ENDIF

#endif

 END SUBROUTINE exchange_horiz_int

! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic boundary conditions,
!> respectively, for 2D-arrays.
!>
!> TODO: This routine requires adjustments to be used for GPUs. The simple check further below
!>       is commented, because the routine is called for output purposes in average_3d_data, so that
!>       the check would stop the run since many of the output quantities do not yet exist on the
!>       device.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d( ar )

#if ! defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc
#endif
#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY:  message_string
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg

    REAL(wp) ::  ar(nysg:nyng,nxlg:nxrg)  !<


!#if defined( _OPENACC )
!    message_string = 'routine not prepared for using openacc'
!    CALL message( 'exchange_horiz_2d', 'PAC0359', 1, 2, 0, 6, 0 )
!#endif


    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( npex == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged within the PE memory
       ar(:,nxlg:nxl-1) = ar(:,nxr-nbgp+1:nxr)
       ar(:,nxr+1:nxrg) = ar(:,nxl:nxl+nbgp-1)

    ELSE
!
!--    Send left boundary, receive right one
       CALL MPI_SENDRECV( ar(nysg,nxl), 1, type_y, pleft,  0,                                      &
                          ar(nysg,nxr+1), 1, type_y, pright, 0,                                    &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nysg,nxr+1-nbgp), 1, type_y, pright,  1,                              &
                          ar(nysg,nxlg), 1, type_y, pleft,   1,                                    &
                          comm2d, status, ierr )


    ENDIF

    IF ( npey == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged within the PE memory
       ar(nysg:nys-1,:) = ar(nyn-nbgp+1:nyn,:)
       ar(nyn+1:nyng,:) = ar(nys:nys+nbgp-1,:)

    ELSE
!
!--    Send front boundary, receive rear one
       CALL MPI_SENDRECV( ar(nys,nxlg), 1, type_x, psouth, 0,                                      &
                          ar(nyn+1,nxlg), 1, type_x, pnorth, 0,                                    &
                          comm2d, status, ierr )
!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn+1-nbgp,nxlg), 1, type_x, pnorth, 1,                               &
                          ar(nysg,nxlg), 1, type_x, psouth, 1,                                     &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxlg:nxl-1) = ar(:,nxr-nbgp+1:nxr)
       ar(:,nxr+1:nxrg) = ar(:,nxl:nxl+nbgp-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nysg:nys-1,:) = ar(nyn-nbgp+1:nyn,:)
       ar(nyn+1:nyng,:) = ar(nys:nys+nbgp-1,:)
    ENDIF

#endif

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic boundary conditions,
!> respectively, for 2D 8-bit integer arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d_byte( ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local )

#if ! defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc
#endif
#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY:  message_string
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    INTEGER(iwp) ::  nbgp_local  !< number of ghost layers to be exchanged
    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side

    INTEGER(ibp), DIMENSION(nys_l-nbgp_local:nyn_l+nbgp_local,                                     &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array


#if defined( _OPENACC )
    message_string = 'routine not prepared for using openacc'
    CALL message( 'exchange_horiz_2d_byte', 'PAC0359', 1, 2, 0, 6, 0 )
#endif

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( npex == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged within the PE memory
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)

    ELSE
!
!--    Send left boundary, receive right one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxl_l),   1,                                         &
                          type_y_byte, pleft,  0,                                                  &
                          ar(nys_l-nbgp_local,nxr_l+1), 1,                                         &
                          type_y_byte, pright, 0,                                                  &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxr_l+1-nbgp_local), 1,                              &
                          type_y_byte, pright, 1,                                                  &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,                              &
                          type_y_byte, pleft,  1,                                                  &
                          comm2d, status, ierr )

    ENDIF

    IF ( npey == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged within the PE memory
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)


    ELSE
!
!--    Send front boundary, receive rear one
       CALL MPI_SENDRECV( ar(nys_l,nxl_l-nbgp_local),   1,                                         &
                          type_x_byte, psouth, 0,                                                  &
                          ar(nyn_l+1,nxl_l-nbgp_local), 1,                                         &
                          type_x_byte, pnorth, 0,                                                  &
                          comm2d, status, ierr )

!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn_l+1-nbgp_local,nxl_l-nbgp_local), 1,                              &
                          type_x_byte, pnorth, 1,                                                  &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,                              &
                          type_x_byte, psouth, 1,                                                  &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)
    ENDIF

#endif

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d_byte


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic boundary conditions,
!> respectively, for 2D 32-bit integer arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d_int( ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local )

#if ! defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc
#endif
#if defined( _OPENACC )
    USE control_parameters,                                                                        &
        ONLY:  message_string
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    INTEGER(iwp) ::  nbgp_local  !< number of ghost layers to be exchanged
    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side

    INTEGER(iwp), DIMENSION(nys_l-nbgp_local:nyn_l+nbgp_local,                                     &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array


#if defined( _OPENACC )
    message_string = 'routine not prepared for using openacc'
    CALL message( 'exchange_horiz_2d_int', 'PAC0359', 1, 2, 0, 6, 0 )
#endif

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( npex == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged within the PE memory
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)

    ELSE
!
!--    Send left boundary, receive right one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxl_l),   1, type_y_int, pleft,  0,                  &
                          ar(nys_l-nbgp_local,nxr_l+1), 1, type_y_int, pright, 0,                  &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxr_l+1-nbgp_local), 1, type_y_int, pright, 1,       &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1, type_y_int, pleft,  1,       &
                          comm2d, status, ierr )

    ENDIF

    IF ( npey == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged within the PE memory
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)


    ELSE
!
!--    Send front boundary, receive rear one
       CALL MPI_SENDRECV( ar(nys_l,nxl_l-nbgp_local),   1, type_x_int, psouth, 0,                  &
                          ar(nyn_l+1,nxl_l-nbgp_local), 1, type_x_int, pnorth, 0,                  &
                          comm2d, status, ierr )

!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn_l+1-nbgp_local,nxl_l-nbgp_local), 1, type_x_int, pnorth, 1,       &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1, type_x_int, psouth, 1,       &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)
    ENDIF

#endif

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d_int

 END MODULE exchange_horiz_mod
