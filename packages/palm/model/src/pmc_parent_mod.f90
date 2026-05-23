!> @file pmc_parent_mod.f90
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
! Authors:
! --------
!> @author Klaus Ketelsen (no affiliation)
!
! Description:
! ------------
!> Parent part of palm model coupler.
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_parent

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

    USE MPI

    USE kinds

    USE pmc_general,                                                                               &
        ONLY:  arraydef,                                                                           &
               childdef,                                                                           &
               da_namedef,                                                                         &
               da_namelen,                                                                         &
               pedef,                                                                              &
               pmc_g_setname,                                                                      &
               pmc_max_array,                                                                      &
               pmc_max_models

    USE pmc_handle_communicator,                                                                   &
        ONLY:  m_model_comm,                                                                       &
               m_model_rank,                                                                       &
               m_model_npes,                                                                       &
               m_to_child_comm,                                                                    &
               m_world_rank,                                                                       &
               pmc_parent_for_child,                                                               &
               use_one_sided_communication

    USE pmc_mpi_wrapper,                                                                           &
        ONLY:  pmc_alloc_mem,                                                                      &
               pmc_bcast,                                                                          &
               pmc_send_to_child,                                                                  &
               pmc_recv_from_child,                                                                &
               pmc_waitall

    IMPLICIT NONE

    INTEGER ::  next_array_in_list = 0  !<

    REAL(wp),DIMENSION(:), POINTER ::  base_array_cp  !< base array for child to parent transfer
    REAL(wp),DIMENSION(:), POINTER ::  base_array_pc  !< base array for parent to child transfer

    TYPE childindexdef
       INTEGER ::  nrpoints  !<

       INTEGER, DIMENSION(:,:), ALLOCATABLE ::  index_list_2d  !<
    END TYPE childindexdef

    TYPE(childdef), DIMENSION(pmc_max_models) ::  children  !<

    TYPE(childindexdef), DIMENSION(pmc_max_models) ::  indchildren  !<

    SAVE

    PRIVATE

!
!-- Public functions.
    PUBLIC pmc_parent_for_child

!
!-- Public variables, constants and types.
    PUBLIC children,                                                                               &
           pmc_p_clear_next_array_list,                                                            &
           pmc_p_fillbuffer,                                                                       &
           pmc_p_finalize,                                                                         &
           pmc_p_getbuffer,                                                                        &
           pmc_p_getnextarray,                                                                     &
           pmc_p_get_child_npes,                                                                   &
           pmc_p_parentinit,                                                                       &
           pmc_p_setind_and_allocmem,                                                              &
           pmc_p_set_active_data_array,                                                            &
           pmc_p_set_dataarray,                                                                    &
           pmc_p_split_and_broadcast_index_list

    INTERFACE pmc_p_clear_next_array_list
       MODULE PROCEDURE pmc_p_clear_next_array_list
    END INTERFACE pmc_p_clear_next_array_list

    INTERFACE pmc_p_fillbuffer
       MODULE PROCEDURE pmc_p_fillbuffer
    END INTERFACE pmc_p_fillbuffer

    INTERFACE pmc_p_finalize
       MODULE PROCEDURE  pmc_p_finalize
    END INTERFACE pmc_p_finalize

    INTERFACE pmc_p_getbuffer
       MODULE PROCEDURE pmc_p_getbuffer
    END INTERFACE pmc_p_getbuffer

    INTERFACE pmc_p_getnextarray
       MODULE PROCEDURE pmc_p_getnextarray
    END INTERFACE pmc_p_getnextarray

    INTERFACE pmc_p_get_child_npes
       MODULE PROCEDURE pmc_p_get_child_npes
    END INTERFACE pmc_p_get_child_npes

    INTERFACE pmc_p_parentinit
       MODULE PROCEDURE  pmc_p_parentinit
    END INTERFACE pmc_p_parentinit

    INTERFACE pmc_p_setind_and_allocmem
       MODULE PROCEDURE pmc_p_setind_and_allocmem
    END INTERFACE pmc_p_setind_and_allocmem

    INTERFACE pmc_p_set_active_data_array
       MODULE PROCEDURE pmc_p_set_active_data_array
    END INTERFACE pmc_p_set_active_data_array

    INTERFACE pmc_p_set_dataarray
       MODULE PROCEDURE pmc_p_set_dataarray_2d
       MODULE PROCEDURE pmc_p_set_dataarray_3d
       MODULE PROCEDURE pmc_p_set_dataarray_ip2d
    END INTERFACE pmc_p_set_dataarray

    INTERFACE pmc_p_split_and_broadcast_index_list
       MODULE PROCEDURE pmc_p_split_and_broadcast_index_list
    END INTERFACE pmc_p_split_and_broadcast_index_list

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> If this model is intended to be a parent, initialize parent part of parent-child data transfer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_parentinit

    INTEGER(iwp) ::  childid  !<
    INTEGER(iwp) ::  i        !<
    INTEGER(iwp) ::  istat    !<
    INTEGER(iwp) ::  j        !<


    DO  i = 1, SIZE( pmc_parent_for_child ) - 1

       childid = pmc_parent_for_child( i )

       children(childid)%model_comm = m_model_comm
       children(childid)%inter_comm = m_to_child_comm(childid)

!
!--    Get rank and size.
       CALL MPI_COMM_RANK( children(childid)%model_comm, children(childid)%model_rank, istat )
       CALL MPI_COMM_SIZE( children(childid)%model_comm, children(childid)%model_npes, istat )
       CALL MPI_COMM_REMOTE_SIZE( children(childid)%inter_comm, children(childid)%inter_npes,      &
                                  istat )

!
!--    Intra communicator is used for MPI_GET. .FALSE. means lower core numbers.
       CALL MPI_INTERCOMM_MERGE( children(childid)%inter_comm, .FALSE.,                            &
                                 children(childid)%intra_comm, istat )
       CALL MPI_COMM_RANK( children(childid)%intra_comm, children(childid)%intra_rank, istat )
!
!--    Allocate the PEs structure (TYPE pe_def) for all children of this parent.
       ALLOCATE( children(childid)%pes(children(childid)%inter_npes) )
!
!--    Allocate array of TYPE arraydef for all child PEs to store information of the transfer array.
       DO  j = 1, children(childid)%inter_npes
         ALLOCATE( children(childid)%pes(j)%array_list(pmc_max_array) )
       ENDDO
!
!--    Although array names on parent and child are the same, they are only set on the child side to
!--    avoid redundancy in the code. Therefore, get them from the child.
       CALL pmc_p_get_da_names_from_child( childid )

    ENDDO

 END SUBROUTINE pmc_p_parentinit


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> PE 0 transfers the index list, which contains all the parent grid points involved in
!> parent-child data transfer to that PE on which this grid cell is located.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_split_and_broadcast_index_list( childid, index_list )

     INTEGER(iwp) ::  i                             !<
     INTEGER(iwp) ::  ian                           !<
     INTEGER(iwp) ::  ip                            !<
     INTEGER(iwp) ::  istat                         !<
     INTEGER(iwp) ::  maximum_number_of_gridpoints  !<

     INTEGER(iwp), INTENT(IN) ::  childid  !<

     INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  index_list  !<

     INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  gridpoints_on_pe  !<

     INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  local_index_list  !<


!
!--  PE0 transfers (IF) and the other PEs receive (ELSE).
     IF ( m_model_rank == 0 )  THEN
!
!--     Compute maximum number of grid points located on one parent PE.
        ALLOCATE( gridpoints_on_pe(0:m_model_npes-1) )
        gridpoints_on_pe = 0

        DO  i = 1, SIZE( index_list, 2 )
           gridpoints_on_pe(index_list(6,i)) = gridpoints_on_pe(index_list(6,i)) + 1
        ENDDO

        maximum_number_of_gridpoints = MAXVAL( gridpoints_on_pe )

!
!--     Allocate temp array for PE dependent transfer of index_list.
        ALLOCATE( local_index_list(SIZE( index_list, 1 ),maximum_number_of_gridpoints) )

        DO  ip = 0, m_model_npes-1
!
!--        Split into parent processes.
           ian = 0

           DO  i = 1, SIZE( index_list, 2 )
              IF ( index_list(6,i) == ip )  THEN
                 ian = ian + 1
                 local_index_list(:,ian) = index_list(:,i)
              ENDIF
           ENDDO
!
!--        Send data to other parent processes.
           IF ( ip == 0 )  THEN
              indchildren(childid)%nrpoints = ian
!
!--           Allocate array for index_list_2d. Note, the array will also be allocated in case
!--           ian = 0, in order to avoid errors when array bounds are checked.
              ALLOCATE( indchildren(childid)%index_list_2d(6,1:ian) )
              IF ( ian > 0 )  THEN
                  indchildren(childid)%index_list_2d(:,1:ian) = local_index_list(:,1:ian)
              ENDIF
           ELSE
              CALL MPI_SEND( ian, 1, MPI_INTEGER, ip, 1000, m_model_comm, istat )
              IF ( ian > 0 )  THEN
                  CALL MPI_SEND( local_index_list, 6*ian, MPI_INTEGER, ip, 1001, m_model_comm,     &
                                 istat )
              ENDIF
           ENDIF
        ENDDO

        DEALLOCATE( local_index_list )
        DEALLOCATE( gridpoints_on_pe )

     ELSE

        CALL MPI_RECV( indchildren(childid)%nrpoints, 1, MPI_INTEGER, 0, 1000, m_model_comm,       &
                       MPI_STATUS_IGNORE, istat )
        ian = indchildren(childid)%nrpoints
!
!--     Allocate array for index_list_2d. Note, the array will also be allocated in case ian=0, in
!--     order to avoid errors when array bounds are checked.
        ALLOCATE( indchildren(childid)%index_list_2d(6,1:ian) )
        IF ( ian > 0 )  THEN
           CALL MPI_RECV( indchildren(childid)%index_list_2d, 6*ian, MPI_INTEGER, 0, 1001,         &
                          m_model_comm, MPI_STATUS_IGNORE, istat)
        ENDIF
     ENDIF

     CALL pmc_p_prepare_index_list_for_child( childid, children(childid),                          &
                                              indchildren(childid)%index_list_2d,                  &
                                              indchildren(childid)%nrpoints )

 END SUBROUTINE pmc_p_split_and_broadcast_index_list


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Before creating an array list with arrays scheduled for parent to child transfer,
!> make sure that the list is empty.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_clear_next_array_list

    IMPLICIT NONE


!
!-- next_array_in_list is a global variable in pmc_parent_mod.
    next_array_in_list = 0

 END SUBROUTINE pmc_p_clear_next_array_list


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gets the next array in the list and returns its name, or just returns .FALSE. if there are
!> no more arrays in the list.
!--------------------------------------------------------------------------------------------------!
 LOGICAL FUNCTION pmc_p_getnextarray( childid, myname )

    CHARACTER(LEN=*), INTENT(OUT) ::  myname  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    TYPE(pedef),    POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar  !<


!
!-- next_array_in_list is a global variable in pmc_parent_mod.
    next_array_in_list = next_array_in_list + 1

!
!-- Array names are the same on all children PEs, so take first PE to get the name.
    ape => children(childid)%pes(1)

    IF ( next_array_in_list > ape%nr_arrays )  THEN
!
!--    All arrays are done.
       pmc_p_getnextarray = .FALSE.

    ELSE

       ar => ape%array_list(next_array_in_list)
       myname = ar%name

!
!--    Return .TRUE. if there is still an array in the list.
       pmc_p_getnextarray = .TRUE.

    ENDIF

 END FUNCTION pmc_p_getnextarray


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add 2d-REAL array to the list of arrays scheduled for the parent-child transfer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_set_dataarray_2d( childid, array, array_2 )

    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER ::  array  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER, OPTIONAL ::  array_2  !<

    TYPE(C_PTR) ::  array_adr   !<
    TYPE(C_PTR) ::  second_adr  !<


    dims      = 1
    nrdims    = 2
    dims(1)   = SIZE( array, 1 )
    dims(2)   = SIZE( array, 2 )
    array_adr = C_LOC( array )

    IF ( PRESENT( array_2 ) )  THEN
       second_adr = C_LOC( array_2 )
       CALL pmc_p_setarray( childid, nrdims, dims, array_adr, second_adr = second_adr )
    ELSE
       CALL pmc_p_setarray( childid, nrdims, dims, array_adr )
    ENDIF

 END SUBROUTINE pmc_p_set_dataarray_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add 2d-INTEGER array to the list of arrays scheduled for the parent-child transfer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_set_dataarray_ip2d( childid, array )

    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(idp), INTENT(IN), DIMENSION(:,:), POINTER ::  array  !<

    TYPE(C_PTR) ::  array_adr  !<


    dims      = 1
    nrdims    = 2
    dims(1)   = SIZE( array, 1 )
    dims(2)   = SIZE( array, 2 )
    array_adr = C_LOC( array )

    CALL pmc_p_setarray( childid, nrdims, dims, array_adr , dimkey = 22 )

 END SUBROUTINE pmc_p_set_dataarray_ip2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add 3d-REAL array to the list of arrays scheduled for the parent-child transfer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_set_dataarray_3d( childid, array, nz_cl, nz, array_2, ks_cl, ke_cl )

    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp), INTENT(IN) ::  nz       !<
    INTEGER(iwp), INTENT(IN) ::  nz_cl    !<
    INTEGER(iwp), INTENT(IN) ::  ks_cl    !<
    INTEGER(iwp), INTENT(IN) ::  ke_cl    !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER ::  array  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER, OPTIONAL ::  array_2  !<

    TYPE(C_PTR) ::  array_adr   !<
    TYPE(C_PTR) ::  second_adr  !<


    nrdims  = 3
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )
    dims(3) = SIZE( array, 3 )
    dims(4) = nz_cl + dims(1) - nz  ! works for first dimension 1:nz as well as for 0:nz+1

    array_adr = C_LOC( array )
!
!-- In PALM's pointer version, two indices have to be stored internally for the two prognostic
!-- time levels. The active address of the data array is set in swap_timelevel.
    IF ( PRESENT( array_2 ) )  THEN
       second_adr = C_LOC( array_2 )
       CALL pmc_p_setarray( childid, nrdims, dims, array_adr, second_adr = second_adr, ks=ks_cl,   &
                            ke=ke_cl )
    ELSE
       CALL pmc_p_setarray( childid, nrdims, dims, array_adr, ks=ks_cl, ke=ke_cl )
    ENDIF

 END SUBROUTINE pmc_p_set_dataarray_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines the transfer buffer size and the buffer indices which are the starting indices of each
!> array (variable) and sends the buffer indices to all child PEs using MPI_ALLTOALL and receives
!> the corresponding indices from the child. After these operations, the transfer buffer is
!> allocated (pmc_alloc_mem) and the RMA-window for parent to child transfer is created and its
!> base-array pointer is assigned to the buffer. Finally, preparations are made for receiving
!> data from the child. This includes determination of the receive buffer size and indices
!> (array starting indices), allocation of the receive buffer (pmc_alloc_mem) and setting the
!> receive buffer pointer to the structure children%pes%array_list%recvbuf.
!>
!> Naming convention for appendices:   _pc  -> parent to child transfer
!>                                     _cp  -> child to parent transfer
!>                                     send -> parent to child transfer
!>                                     recv -> child to parent transfer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_setind_and_allocmem( childid )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    INTEGER(iwp) ::  arlen            !<
    INTEGER(iwp) ::  i                !<
    INTEGER(iwp) ::  ierr             !<
    INTEGER(iwp) ::  j                !<
    INTEGER(iwp) ::  local_nr_arrays  !< store number of arrays in  local variiab le
    INTEGER(iwp) ::  myindex          !<
    INTEGER(iwp) ::  total_npes       !< total number of PEs parent and child

    INTEGER(idp) ::  bufsize          !< size of MPI data window

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_s  !< array of send indices to be sent to children
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_r  !< array of receive indices to be received from children

    TYPE(C_PTR) ::  base_ptr  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<


    CALL MPI_COMM_SIZE( children(childid)%intra_comm, total_npes, ierr )
!
!-- Parent to child direction.
    myindex = 1
!    
!-- An initial (base) buffer size is required because inactive processes need a non-zero size in 
!-- MPI_WIN_CREATE. Eight is just an arbitrary choice.
    bufsize = 8

!
!-- All child PEs get the same number of arrays. Therefore, the number of arrays from the first
!-- child PE can be used as the dimension.
    local_nr_arrays = children(childid)%pes(1)%nr_arrays

    ALLOCATE( myindex_s(local_nr_arrays,0:total_npes-1) )
    ALLOCATE( myindex_r(local_nr_arrays,0:total_npes-1) )

    myindex_s = 0

!
!-- First stride: compute size and set index.
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)

       DO  j = 1, ape%nr_arrays

          ar  => ape%array_list(j)
          IF ( ar%nrdims == 2 )  THEN
             arlen = ape%nrele
          ELSEIF ( ar%nrdims == 3 )  THEN
             arlen = ape%nrele * ar%a_dim(4)
          ELSE
             arlen = -1
          ENDIF
          ar%sendindex = myindex
!
!--       Using intra communicator for MPI_ALLTOALL, the child PE numbers are after the parent ones.
          myindex_s(j,i-1+children(childid)%model_npes) = myindex

          myindex = myindex + arlen
          bufsize = bufsize + arlen
          ar%sendsize = arlen

       ENDDO

    ENDDO

!
!-- Using MPI_ALLTOALL to send indices from parent to children. The data coming back from the
!-- child PEs are ignored. The respective call on the child side is in pmc_c_setind_and_allocmem.
    CALL MPI_ALLTOALL( myindex_s, local_nr_arrays, MPI_INTEGER, myindex_r, local_nr_arrays,        &
                       MPI_INTEGER, children(childid)%intra_comm, ierr )

!
!-- Using MPI_Alltoall to receive indices from children.
    myindex_s = 0
    myindex_r = 0
    CALL MPI_ALLTOALL( myindex_s, local_nr_arrays, MPI_INTEGER, myindex_r, local_nr_arrays,        &
                       MPI_INTEGER, children(childid)%intra_comm, ierr )
!
!-- Create RMA (One Sided Communication) window for data buffer parent to child transfer.
!-- The buffer of MPI_GET (counterpart of transfer) can be PE-local, i.e. it can but must not be
!-- part of the MPI RMA window. Only one RMA window is required to prepare the data for:
!--          parent -> child transfer on the parent side
!-- and for:
!--          child -> parent transfer on the child side
!-- ATTENTION: In the following it is assumed that all data type variables have the size of wp.
!--            This does not work with wp=sp, since particle INTEGER arrays are always INTEGER8.
!--            Error PMC0039 is issued in such a case.
    IF ( use_one_sided_communication )  THEN

       CALL pmc_alloc_mem( base_array_pc, bufsize )
       children(childid)%totalbuffersize = bufsize * wp

       winsize = bufsize * wp
       CALL MPI_WIN_CREATE( base_array_pc, winsize, wp, MPI_INFO_NULL,                             &
                            children(childid)%intra_comm, children(childid)%win_parent_child, ierr )
    ELSE

       ALLOCATE( base_array_pc( bufsize ) )

    ENDIF

!
!-- Second stride: set buffer pointer.
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%sendbuf = C_LOC( base_array_pc(ar%sendindex) )
          IF ( ar%sendindex + ar%sendsize > bufsize )  THEN
             WRITE( message_string, '(A,I4,4I7,1X,A)' ) 'parent buffer too small ',i ,             &
                    ar%sendindex, ar%sendsize, ar%sendindex + ar%sendsize, bufsize, TRIM( ar%name )
             CALL message( 'pmc_p_setind_and_allocmem', 'PMC0004', 3, 2, 0, 6, 0 )
          ENDIF
       ENDDO

    ENDDO

!
!-- Child to parent direction.
!-- An initial (base) buffer size is required because inactive processes need a non-zero size in 
!-- MPI_WIN_CREATE. Eight is just an arbitrary choice.
    bufsize = 8

!
!-- First stride: compute size and set index.
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
!
!--       Receive index from child.
          IF ( ar%nrdims == 3 )  THEN
             bufsize = MAX( bufsize, INT( ape%nrele * ar%a_dim(4), MPI_ADDRESS_KIND ) )
          ELSE
             bufsize = MAX( bufsize, INT( ape%nrele, MPI_ADDRESS_KIND ) )
          ENDIF
          ar%recvindex = myindex_r(j,i-1+children(childid)%model_npes)
       ENDDO

    ENDDO

    DEALLOCATE( myindex_s )
    DEALLOCATE( myindex_r )

!
!-- Create RMA data buffer. The buffer for MPI_GET can be PE local, i.e. it can but must not be
!-- part of the MPI RMA window.
    CALL pmc_alloc_mem( base_array_cp, bufsize, base_ptr )
    children(childid)%totalbuffersize = bufsize * wp

    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

!
!-- Second stride: set buffer pointer.
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvbuf = base_ptr
       ENDDO

    ENDDO

 END SUBROUTINE pmc_p_setind_and_allocmem


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fill buffer in the RMA window to enable the child to fetch the data with MPI_GET. Or send the
!> data to child using pmc_send_to_child in case of two-sided communication.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_fillbuffer( childid, particle_transfer )

    INTEGER(iwp) ::  ierr     !<
    INTEGER(iwp) ::  ij       !<
    INTEGER(iwp) ::  ip       !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  myindex  !<
    INTEGER(iwp) ::  nr_ar    !< number of arrays
    INTEGER(iwp) ::  nr_req   !< number of non blocking requests


    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(iwp), DIMENSION(1) ::  buf_shape  !<

    INTEGER(idp), POINTER, DIMENSION(:) ::  ibuf  !<

    INTEGER(iwp), DIMENSION(children(childid)%inter_npes) :: req  !< array of requests

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  idata_2d  !<

    LOGICAL ::  particle_transfer_l  !< local variable with default .FALSE. to hold value of optional argument particle_transfer

    LOGICAL, INTENT(IN), OPTIONAL ::  particle_transfer  !<

    REAL(wp), POINTER, DIMENSION(:)     ::  buf      !<
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(arraydef), POINTER ::  ar   !<
    TYPE(pedef),    POINTER ::  ape  !<

    TYPE :: bufdef
       INTEGER(iwp)                            ::  nr    !<
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  buf   !<
       INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  ibuf  !<
    END TYPE

    TYPE(bufdef), DIMENSION(children(childid)%inter_npes) ::  bufall  !< data buffer for all requests


!
!-- Set local switch to determine if particle transfer is active.
    particle_transfer_l = .FALSE.
    IF ( PRESENT( particle_transfer ) )    particle_transfer_l = particle_transfer

    IF ( use_one_sided_communication )  THEN
!
!--    The RMA window is using passive target synchronization, therefore no MPI_WIN... calls here.
       DO  ip = 1, children(childid)%inter_npes

          ape => children(childid)%pes(ip)

          DO  j = 1, ape%nr_arrays

             ar => ape%array_list(j)
             myindex = 1
             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN
!
!--             2d-REAL*8 array.
                buf_shape(1) = ape%nrele
                CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
                CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
!
!--             Copy from PALM 2-d array to sendbuf in RMA window.
                DO  ij = 1, ape%nrele
                   buf(myindex) = data_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                   myindex = myindex + 1
                ENDDO

             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
!
!--             3d-REAL*8 array.
                buf_shape(1) = ape%nrele*ar%a_dim(4)
                CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
                CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
!
!--             Copy from PALM 3d-array (e.g. u,v,...) to sendbuf in RMA window.
                DO  ij = 1, ape%nrele
                   buf(myindex:myindex+ar%a_dim(4)-1) =                                            &
                                           data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i)
                   myindex = myindex + ar%a_dim(4)
                ENDDO

             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
!
!--             2d-INTEGER*8 array for particle transfer.
                buf_shape(1) = ape%nrele
                CALL C_F_POINTER( ar%sendbuf, ibuf, buf_shape )
                CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   ibuf(myindex) = idata_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                   myindex = myindex + 1
                ENDDO

             ENDIF

          ENDDO

       ENDDO

!
!--    Parent has filled the RMA Window. The child uses MPI_LOCK/MPI_UNLOCK for synchronization.
!--    The next barrier ensures that the child is not starting to get the data.
       CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
!
!--    Wait until child has received all data. There is a correponding barrier in pmc_c_getbuffer
!--    after the child received all data.
!--    The barrier is just for security reasons. It could be removed.
       CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

    ELSE
!
!--    Two-sided communication.

       nr_ar = 0
       DO  ip = 1, children(childid)%inter_npes
          ape => children(childid)%pes(ip)
          nr_ar = MAX( nr_ar, ape%nr_arrays )
       ENDDO

       DO  j = 1, nr_ar

          req    = 0
          nr_req = 0
!
!--       Using one buffer for every ip value is the safe version. It may not be necessary on send
!--       side.
          DO  ip = 1, children(childid)%inter_npes
             IF ( ALLOCATED(bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
             IF ( ALLOCATED(bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
          ENDDO

          DO  ip = 1, children(childid)%inter_npes

             ape => children(childid)%pes(ip)
             ar  => ape%array_list(j)

             IF ( j > ape%nr_arrays )  CYCLE

             myindex = 1
             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN
!
!--             2d-REAL*8 array.
                IF ( ape%nrele > 0 )  THEN
                   ALLOCATE( bufall(ip)%buf(ape%nrele) )
                   CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
                   DO  ij = 1, ape%nrele
                      bufall(ip)%buf(myindex) = data_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                      myindex = myindex + 1
                   ENDDO
                   nr_req = nr_req + 1
                   CALL pmc_send_to_child( childid, bufall(ip)%buf, ape%nrele, ip-1, 2100+j,       &
                                           req=req(nr_req), ierr=ierr )
                ENDIF

             ELSEIF ( ar%dimkey == 3  .AND. .NOT. particle_transfer_l )  THEN
!
!--             3d-REAL*8 array.
                IF ( ape%nrele > 0 )  THEN
                   ALLOCATE( bufall(ip)%buf(ape%nrele*ar%a_dim(4)) )
                   CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
                   DO  ij = 1, ape%nrele
                      bufall(ip)%buf(myindex:myindex+ar%a_dim(4)-1) =                              &
                                           data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i)
                      myindex = myindex + ar%a_dim(4)
                   ENDDO
                   nr_req = nr_req + 1
                   CALL pmc_send_to_child( childid, bufall(ip)%buf, ape%nrele*ar%a_dim(4), ip-1,   &
                                           2100+j, req=req(nr_req), ierr=ierr )
                ENDIF

             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
!
!--             2d-INTEGER*8 array for particle transfer.
                IF ( ape%nrele > 0 )  THEN
                   ALLOCATE( bufall(ip)%ibuf(ape%nrele) )
                   CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )
                   DO  ij = 1, ape%nrele
                      bufall(ip)%ibuf(myindex) = idata_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                      myindex = myindex + 1
                   ENDDO
                   nr_req = nr_req + 1
                   CALL pmc_send_to_child( childid, bufall(ip)%ibuf, ape%nrele, ip-1, 2200+j,      &
                                           req=req(nr_req), ierr=ierr )
                ENDIF

             ENDIF

          ENDDO

          IF ( nr_req > 0 )  CALL pmc_waitall( nr_req, req, ierr )

       ENDDO

       DO  ip = 1, children(childid)%inter_npes
          IF ( ALLOCATED( bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
          IF ( ALLOCATED( bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
       ENDDO

    ENDIF

 END SUBROUTINE pmc_p_fillbuffer


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Receive data from child by MPI_GET. Or receive the data from child using pmc_recv_from_child
!> in case of two-sided communication.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_getbuffer( childid, particle_transfer, child_process_nr )

    INTEGER(iwp) ::  ierr       !<
    INTEGER(iwp) ::  ij         !<
    INTEGER(iwp) ::  ip         !<
    INTEGER(iwp) ::  ip_start   !<
    INTEGER(iwp) ::  ip_end     !<
    INTEGER(iwp) ::  j          !<
    INTEGER(iwp) ::  myindex    !<
    INTEGER(iwp) ::  nr         !<
    INTEGER(iwp) ::  nr_ar      !<
    INTEGER(iwp) ::  nr_req     !<
    INTEGER(iwp) ::  target_pe  !<


    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  child_process_nr  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_disp  !<

    INTEGER(idp), POINTER, DIMENSION(:) ::  ibuf  !<

    INTEGER(iwp), DIMENSION(children(childid)%inter_npes) ::  req  !<

    INTEGER(iwp), DIMENSION(1) ::  buf_shape  !<

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  idata_2d  !<

    LOGICAL ::  particle_transfer_l  !<

    LOGICAL, INTENT(IN), OPTIONAL ::  particle_transfer  !<

    REAL(wp), POINTER, DIMENSION(:)     ::  buf      !<
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(arraydef), POINTER ::  ar   !<
    TYPE(pedef),    POINTER ::  ape  !<

    TYPE :: bufdef
       INTEGER(iwp)                            ::  nr    !<
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  buf   !<
       INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  ibuf  !<
    END TYPE

    TYPE(bufdef), DIMENSION(children(childid)%inter_npes) ::  bufall  !<


!
!-- Set local switch to determine if particle transfer is active.
    particle_transfer_l = .FALSE.
    IF ( PRESENT( particle_transfer ) )  particle_transfer_l = particle_transfer

    IF( PRESENT( child_process_nr ) )  THEN
       ip_start = child_process_nr
       ip_end   = child_process_nr
    ELSE
       ip_start = 1
       ip_end   = children(childid)%inter_npes
    ENDIF

    IF ( ip_start == 1  .AND.   use_one_sided_communication )  THEN
!
!--    Wait for child to fill buffer.
       CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
    ENDIF

    IF ( use_one_sided_communication )  THEN

       DO  ip = ip_start, ip_end

          ape => children(childid)%pes(ip)

          DO  j = 1, ape%nr_arrays

             ar => ape%array_list(j)

             IF ( ar%recvindex < 0 )  CYCLE

             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele
             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele * ar%a_dim(4)
             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
                nr = ape%nrele
             ELSE
!
!--             The value of particle_transfer_l is .T., if pmc_p_getbuffer is called from
!--             pmc_particle_interface, and .F. if called from pmc_interface_mod. Depending on this,
!--             arrays that do or that do not belong to the LPM are skipped (not transferred).
                CYCLE
             ENDIF

             buf_shape(1) = nr

             IF ( particle_transfer_l )  THEN
!
!--             Here only arrays with dimkey=22 are treated.
                CALL C_F_POINTER( ar%recvbuf, ibuf, buf_shape )
             ELSE
!
!--             Here only arrays with dimkey=2/3 are treated.
                CALL C_F_POINTER( ar%recvbuf, buf, buf_shape )
             ENDIF

!
!--          MPI passive target RMA.
             IF ( nr > 0 )  THEN

                target_disp = ar%recvindex - 1
!
!--             Child processes are located behind parent process.
                target_pe = ip - 1 + m_model_npes
                CALL MPI_WIN_LOCK( MPI_LOCK_SHARED, target_pe, 0,                                  &
                                   children(childid)%win_parent_child, ierr )
                IF ( particle_transfer_l )  THEN
!
!--                Here only arrays with dimkey=22 are treated.
                   CALL MPI_GET( ibuf, nr, MPI_INTEGER8, target_pe, target_disp, nr, MPI_INTEGER8, &
                                 children(childid)%win_parent_child, ierr )
                ELSE
!
!--                Here only arrays with dimkey=2 or 3 are treated.
                   CALL MPI_GET( buf, nr, MPI_REAL, target_pe, target_disp, nr, MPI_REAL,          &
                                 children(childid)%win_parent_child, ierr )
                ENDIF
                CALL MPI_WIN_UNLOCK( target_pe, children(childid)%win_parent_child, ierr )

             ENDIF

             myindex = 1

             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   data_2d(ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex)
                   myindex = myindex + 1
                ENDDO

             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
                DO  ij = 1, ape%nrele
                   data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i) =                      &
                                                                 buf(myindex:myindex+ar%a_dim(4)-1)
                   myindex = myindex + ar%a_dim(4)
                ENDDO

             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   idata_2d(ape%locind(ij)%j,ape%locind(ij)%i) = ibuf(myindex)
                   myindex = myindex + 1
                ENDDO

             ENDIF

          ENDDO

       ENDDO

       IF ( ip_start == 1 )  CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

    ELSE
!
!--    Two-sided communication. 
!--    For non blocking receive (MPI_IRECV) additional buffers are required to supply space for data
!--    of all receives. The data is completly available after pmc_waitall.
       nr_ar = 0

       DO  ip = 1, ip_start, ip_end
          ape => children(childid)%pes(ip)
          nr_ar = MAX( nr_ar, ape%nr_arrays )
       ENDDO

       DO  j = 1, nr_ar

          req    = 0
          nr_req = 0
          DO  ip = ip_start, ip_end
             IF ( ALLOCATED( bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
             IF ( ALLOCATED( bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
          ENDDO

          DO  ip = ip_start, ip_end

             ape => children(childid)%pes(ip)
             ar  => ape%array_list(j)

             IF ( j > ape%nr_arrays )  CYCLE

             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele
             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele * ar%a_dim(4)
             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
                nr = ape%nrele
             ELSE
!
!--             The value of particle_transfer_l is .T., if pmc_p_getbuffer is called from
!--             pmc_particle_interface, and .F. if called from pmc_interface_mod. Depending on this,
!--             arrays that do or that do not belong to the LPM are skipped (not transferred).
                CYCLE
             ENDIF

             bufall(ip)%nr = nr
             buf_shape(1)  = nr

             IF ( nr > 0 )  THEN
                IF ( particle_transfer_l )  THEN
!
!--                Here only arrays with dimkey=22 are treated.
                   ALLOCATE( bufall(ip)%ibuf(nr) )
                   nr_req = nr_req + 1
                   CALL pmc_recv_from_child( childid, bufall(ip)%ibuf ,nr , ip-1, 3200+j,          &
                                             req=req(nr_req), ierr=ierr )
                ELSE
!
!--                Here only arrays with dimkey=2/3 are treated.
                   ALLOCATE( bufall(ip)%buf(nr) )
                   bufall(ip)%buf = 0.0
                   nr_req = nr_req + 1
                   CALL pmc_recv_from_child( childid, bufall(ip)%buf ,nr , ip-1, 3100+j,           &
                                             req=req(nr_req), ierr=ierr )
                ENDIF
             ENDIF

          ENDDO

          IF ( nr_req > 0)   CALL pmc_waitall( nr_req, req, ierr )

          DO  ip = ip_start, ip_end

             ape => children(childid)%pes(ip)
             ar  => ape%array_list(j)

             IF ( j > ape%nr_arrays )  CYCLE

             nr = bufall(ip)%nr

             IF ( nr <= 0 )  CYCLE

             myindex = 1
             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   data_2d(ape%locind(ij)%j,ape%locind(ij)%i) = bufall(ip)%buf(myindex)
                   myindex = myindex + 1
                ENDDO

             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN

                NULLIFY( data_3d )
                CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
                DO  ij = 1, ape%nrele
                   data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i) =                      &
                                                      bufall(ip)%buf(myindex:myindex+ar%a_dim(4)-1)
                   myindex = myindex + ar%a_dim(4)
                ENDDO

             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   idata_2d(ape%locind(ij)%j,ape%locind(ij)%i) = bufall(ip)%ibuf(myindex)
                   myindex = myindex + 1
                ENDDO

             ENDIF

          ENDDO

       ENDDO

       DO  ip = ip_start, ip_end
          IF ( ALLOCATED(bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
          IF ( ALLOCATED(bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
       ENDDO

    ENDIF

 END SUBROUTINE pmc_p_getbuffer


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Broadcast the name and couple_index of arrays to be transferred from child PE 0 to parent PEs.
!> Then call pmc_g_setname to set the name and couple index in children%pes%array_list. 
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_get_da_names_from_child( childid )

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    TYPE(da_namedef) ::  myname  !<


    DO
!
!--    This loop is over the arrays to be transferred. Note the exit mechanism.
!--    myname%nameonparent and myname%couple_index come from the child side for 
!--    each array in a sequence until couple_index = -1 is found.    
       CALL pmc_bcast( myname%couple_index, 0, comm=m_to_child_comm(childid) )

       IF ( myname%couple_index == -1 )  EXIT

       CALL pmc_bcast( myname%nameonparent, 0, comm=m_to_child_comm(childid) )
!
!--    This call sets the array name and couple id into hildren%pes%array_list. 
       CALL pmc_g_setname( children(childid), myname%couple_index, myname%nameonparent )

    ENDDO

 END SUBROUTINE pmc_p_get_da_names_from_child


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the next array in the list and its dimensions into children%npes%array_list for all
!> child PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_setarray( childid, nrdims, dims, array_adr, second_adr, dimkey, ks, ke )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp), INTENT(IN) ::  nrdims   !<
    INTEGER(iwp), INTENT(IN), OPTIONAL :: ks  !<
    INTEGER(iwp), INTENT(IN), OPTIONAL :: ke  !<

    INTEGER(iwp), INTENT(IN), OPTIONAL :: dimkey  !<

    INTEGER(iwp), INTENT(IN), DIMENSION(:) :: dims  !<

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  local_dimkey  !<

    TYPE(C_PTR), INTENT(IN) :: array_adr  !<

    TYPE(C_PTR), INTENT(IN), OPTIONAL ::  second_adr  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<


!
!-- Check for valid dimkey.
    IF ( PRESENT( dimkey ) )  THEN
       local_dimkey = dimkey
    ELSE
       local_dimkey = nrdims
    ENDIF

    IF ( local_dimkey /= 2  .AND.  local_dimkey /= 3  .AND.  local_dimkey /= 22 )  THEN
       WRITE( message_string, '(A,I4)' )  'invalid dimkey = ', local_dimkey
       CALL message( 'pmc_p_setarray', 'PMC0041', 3, 2, 0, 6, 0 )
    ENDIF

    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       IF ( PRESENT( dimkey ) )  THEN
          ar%dimkey = dimkey
       ELSE
          ar%dimkey = nrdims
       ENDIF

       ar%a_dim  = dims
       ar%data   = array_adr
       IF ( PRESENT( second_adr ) )  THEN
          ar%po_data(1) = array_adr
          ar%po_data(2) = second_adr
       ELSE
          ar%po_data(1) = C_NULL_PTR
          ar%po_data(2) = C_NULL_PTR
       ENDIF

       IF( PRESENT (ks) )  THEN
          ar%ks = ks
       ENDIF
       IF( PRESENT (ke) )  THEN
          ar%ke = ke
       ENDIF

    ENDDO

 END SUBROUTINE pmc_p_setarray


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the data pointer to the current time level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_set_active_data_array( childid, iactive )

    INTEGER(iwp) :: ip  !<
    INTEGER(iwp) :: j   !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp), INTENT(IN) ::  iactive  !<

    TYPE(pedef), POINTER ::  ape  !<

    TYPE(arraydef), POINTER ::  ar   !<


    DO  ip = 1, children(childid)%inter_npes
       ape => children(childid)%pes(ip)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
!
!--       2d-array variables (dimkey = 22) are not treated because they do not have two time levels.
          IF ( MOD( ar%dimkey, 10 ) == 2 )  CYCLE
          ar%data = ar%po_data(iactive)
       ENDDO
    ENDDO

 END SUBROUTINE pmc_p_set_active_data_array


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gets the number of child PEs.
!--------------------------------------------------------------------------------------------------!
 INTEGER FUNCTION pmc_p_get_child_npes( child_id )

   INTEGER(iwp), INTENT(IN) ::  child_id  !<


   pmc_p_get_child_npes = children(child_id)%inter_npes

   RETURN

 END FUNCTION pmc_p_get_child_npes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prepare and transfer the local index list of parent PEs to the current child PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_prepare_index_list_for_child( childid, mychild, index_list, nrp )

    IMPLICIT NONE

    INTEGER(iwp) ::  arraysize  !< size of array to hold the index list
    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  ierr       !<
    INTEGER(iwp) ::  ind        !<
    INTEGER(iwp) ::  i2         !<
    INTEGER(iwp) ::  j          !<
    INTEGER(iwp) ::  rempe      !<

    INTEGER(iwp), INTENT(IN) ::  childid  !<
    INTEGER(iwp), INTENT(IN) ::  nrp      !<

    INTEGER(iwp), INTENT(IN), DIMENSION(:,:) ::  index_list  !<

    TYPE(childdef), INTENT(INOUT) ::  mychild  !<

    INTEGER(iwp), DIMENSION(mychild%inter_npes) ::  remind  !<

    INTEGER(iwp), DIMENSION(:), POINTER ::  remindw  !<

    INTEGER(iwp), DIMENSION(mychild%inter_npes)   ::  child_nr_val       !<
    INTEGER(iwp), DIMENSION(mychild%inter_npes)   ::  child_start_index  !<
    INTEGER(iwp), DIMENSION(mychild%inter_npes*2) ::  rldef              !<

    TYPE(pedef), POINTER ::  ape  !<


!
!-- First, clear the number of elements for every remote child PE.
    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ape%nrele = 0
    ENDDO

!
!-- Count the number of parent grid cells to be sent to the individual child PEs.
    DO  j = 1, nrp
!
!--    PE number on remote PE
       rempe = index_list(5,j) + 1
       ape => mychild%pes(rempe)
!
!--    Increment the number of grid points in the index_list for this child PE.
       ape%nrele = ape%nrele + 1
    ENDDO

    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ALLOCATE( ape%locind(ape%nrele) )
    ENDDO

    remind = 0

!
!-- Second, create lists. Loop over number of parent grid points.
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1
       ape => mychild%pes(rempe)
       remind(rempe)     = remind(rempe) + 1
       ind               = remind(rempe)
       ape%locind(ind)%i = index_list(1,j)
       ape%locind(ind)%j = index_list(2,j)
    ENDDO

    rldef(1) = 0          ! index on remote PE 0
    rldef(2) = remind(1)  ! number of elements on remote PE 0

    CALL pmc_send_to_child( childid, rldef, 2, 0, 1001, ierr )
!
!-- Reserve buffer for index array.
    DO  i = 2, mychild%inter_npes
       i2          = ( i - 1 ) * 2 + 1
       rldef(i2)   = rldef(i2-2) + rldef(i2-1) * 2  ! index on remote PE
       rldef(i2+1) = remind(i)                      ! number of grid points on remote PE
       CALL pmc_send_to_child( childid, rldef(i2:i2+1), 2, i-1, 1001, ierr )
    ENDDO

    i2 = 2 * mychild%inter_npes - 1
    arraysize = ( rldef(i2) + rldef(i2+1) ) * 2
    arraysize = MAX( arraysize, 1 )

!
!-- Create the 2D index list and send to child.
    ALLOCATE( remindw(arraysize) )
    child_start_index = -1
    child_nr_val      = 0

    DO  j = 1, nrp
!
!--    PE number on remote PE.
       rempe = index_list(5,j) + 1

       ape => mychild%pes(rempe)
       i2    = rempe * 2 - 1
       ind   = rldef(i2) + 1
       remindw(ind)   = index_list(3,j)
       remindw(ind+1) = index_list(4,j)
       rldef(i2)      = rldef(i2) + 2
       IF ( child_start_index(rempe) == -1 )  THEN
          child_start_index(rempe) = ind
       ENDIF
       child_nr_val(rempe) = child_nr_val(rempe) + 2

    ENDDO

    DO  i = 1, mychild%inter_npes
       IF ( child_nr_val(i) > 0 )  THEN
          CALL pmc_send_to_child( childid, remindw(child_start_index(i):), child_nr_val(i), i-1,   &
                                  1002, ierr )
       ENDIF
    ENDDO

 END SUBROUTINE pmc_p_prepare_index_list_for_child


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Free window and memory of pmc data buffer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_p_finalize( childid )

    INTEGER(iwp), INTENT(IN) ::  childid  !<

    INTEGER(iwp) :: ierr  !<


    IF ( use_one_sided_communication )  THEN
       CALL MPI_FREE_MEM( base_array_pc, ierr )
       CALL MPI_WIN_FREE( children(childid)%win_parent_child, ierr )
    ENDIF

 END SUBROUTINE pmc_p_finalize
#endif

 END MODULE pmc_parent
