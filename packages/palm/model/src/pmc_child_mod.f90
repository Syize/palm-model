MODULE pmc_child

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
!> Child part of Palm Model Coupler.
!--------------------------------------------------------------------------------------------------!

#if defined( __parallel )

    USE, INTRINSIC ::  iso_c_binding

    USE MPI

    USE kinds

    USE pmc_general,                                                                               &
        ONLY:  arraydef,                                                                           &
               childdef,                                                                           &
               da_desclen,                                                                         &
               da_namedef,                                                                         &
               da_namelen,                                                                         &
               pedef,                                                                              &
               pmc_da_name_err,                                                                    &
               pmc_g_setname,                                                                      &
               pmc_max_array

    USE pmc_handle_communicator,                                                                   &
        ONLY:  m_model_comm,                                                                       &
               m_model_npes,                                                                       &
               m_model_rank,                                                                       &
               m_to_parent_comm,                                                                   &
               use_one_sided_communication

    USE pmc_mpi_wrapper,                                                                           &
        ONLY:  pmc_alloc_mem,                                                                      &
               pmc_bcast,                                                                          &
               pmc_inter_bcast,                                                                    &
               pmc_recv_from_parent,                                                               &
               pmc_send_to_parent,                                                                 &
               pmc_waitall

    IMPLICIT NONE

    PRIVATE

    INTEGER(iwp) ::  myindex = 0             !< counter and unique number for data arrays
    INTEGER(iwp) ::  next_array_in_list = 0  !< index of array in the list of arrays to be coupled

    REAL(wp), DIMENSION(:), POINTER ::  base_array_cp  !< 1d-array used as basis to assign (distribute the data buffers)
    REAL(wp), DIMENSION(:), POINTER ::  base_array_pc  !< 1d-array used as basis to assign (distribute the data buffers)

    TYPE(childdef), PUBLIC ::  me  !< child data structure on the child side 

    SAVE

    INTERFACE pmc_c_childinit
        MODULE PROCEDURE pmc_c_childinit
    END INTERFACE pmc_c_childinit

    INTERFACE pmc_c_clear_next_array_list
        MODULE PROCEDURE pmc_c_clear_next_array_list
    END INTERFACE pmc_c_clear_next_array_list

    INTERFACE pmc_c_finalize
        MODULE PROCEDURE pmc_c_finalize
    END INTERFACE pmc_c_finalize

    INTERFACE pmc_c_getbuffer
        MODULE PROCEDURE pmc_c_getbuffer
    END INTERFACE pmc_c_getbuffer

    INTERFACE pmc_c_getnextarray
        MODULE PROCEDURE pmc_c_getnextarray
    END INTERFACE pmc_c_getnextarray

    INTERFACE pmc_c_get_2d_index_list
        MODULE PROCEDURE pmc_c_get_2d_index_list
    END INTERFACE pmc_c_get_2d_index_list

    INTERFACE pmc_c_putbuffer
        MODULE PROCEDURE pmc_c_putbuffer
    END INTERFACE pmc_c_putbuffer

    INTERFACE pmc_c_setind_and_allocmem
        MODULE PROCEDURE pmc_c_setind_and_allocmem
    END INTERFACE pmc_c_setind_and_allocmem

    INTERFACE pmc_c_set_dataarray
        MODULE PROCEDURE pmc_c_set_dataarray_2d
        MODULE PROCEDURE pmc_c_set_dataarray_3d
        MODULE PROCEDURE pmc_c_set_dataarray_ip2d
    END INTERFACE pmc_c_set_dataarray

    INTERFACE pmc_c_set_dataarray_name
        MODULE PROCEDURE pmc_c_set_dataarray_name
        MODULE PROCEDURE pmc_c_set_dataarray_name_lastentry
    END INTERFACE pmc_c_set_dataarray_name

    PUBLIC pmc_c_childinit,                                                                        &
           pmc_c_clear_next_array_list,                                                            &
           pmc_c_finalize,                                                                         &
           pmc_c_getbuffer,                                                                        &
           pmc_c_getnextarray,                                                                     &
           pmc_c_putbuffer,                                                                        &
           pmc_c_setind_and_allocmem,                                                              &
           pmc_c_set_dataarray,                                                                    &
           pmc_c_set_dataarray_name,                                                               &
           pmc_c_get_2d_index_list

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> If this model is intended to be a child, i.e. is not the root model, initialize child part
!> of parent-child data transfer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_childinit

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  istat  !<


!
!-- Get/define the MPI environment.
    me%model_comm = m_model_comm
    me%inter_comm = m_to_parent_comm

    CALL MPI_COMM_RANK( me%model_comm, me%model_rank, istat )
    CALL MPI_COMM_SIZE( me%model_comm, me%model_npes, istat )
    CALL MPI_COMM_REMOTE_SIZE( me%inter_comm, me%inter_npes, istat )
!
!-- Intra-communicator is used for MPI_GET. .TRUE. means high core numbers
    CALL MPI_INTERCOMM_MERGE( me%inter_comm, .TRUE., me%intra_comm, istat )
    CALL MPI_COMM_RANK( me%intra_comm, me%intra_rank, istat )

    ALLOCATE( me%pes(me%inter_npes) )
!
!-- Allocate an array of type arraydef for all parent processes to store information of the
!-- transfer array.
    DO  i = 1, me%inter_npes
       ALLOCATE( me%pes(i)%array_list(pmc_max_array) )
    ENDDO

 END SUBROUTINE pmc_c_childinit


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store the name of arrays which shall be coupled in a list.
!> The names on parent and child side are always identical, but stored independently for historical
!> reasons. Maybe simplified later. 
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_set_dataarray_name( arrayname )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  arrayname  !< name of the current array in the list

    INTEGER(iwp) ::  mype  !<

    TYPE(da_namedef) ::  myname  !<


!
!-- Check length of array names.
    IF ( LEN( TRIM( arrayname) ) > da_namelen )  THEN
       message_string = 'PMC array name "' // TRIM( arrayname ) // '" too long'
       CALL message( 'pmc_c_set_dataarray_name', 'PMC0040', 3, 2, 0, 6, 0 )
    ENDIF

    IF ( m_model_rank == 0 )  THEN
       myindex = myindex + 1
       myname%couple_index = myindex
       myname%parentdesc   = 'parent'
       myname%nameonparent = TRIM( arrayname )
       myname%childdesc    = 'child'
       myname%nameonchild  = TRIM( arrayname )
    ENDIF

!
!-- Broadcast the complete description of a transfer array to all child PEs.
    CALL pmc_bcast( myname%couple_index, 0, comm=m_model_comm )
    CALL pmc_bcast( myname%parentdesc,   0, comm=m_model_comm )
    CALL pmc_bcast( myname%nameonparent, 0, comm=m_model_comm )
    CALL pmc_bcast( myname%childdesc,    0, comm=m_model_comm )
    CALL pmc_bcast( myname%nameonchild,  0, comm=m_model_comm )

!
!-- Broadcast the complete description of a transfer array to all parent PEs.
!-- Only the root PE of these broadcasts to parent PEs is using intra communicator.
    IF ( m_model_rank == 0 )  THEN
        mype = MPI_ROOT
    ELSE
        mype = MPI_PROC_NULL
    ENDIF

    CALL pmc_bcast( myname%couple_index, mype, comm=m_to_parent_comm )
    CALL pmc_bcast( myname%nameonparent, mype, comm=m_to_parent_comm )

    CALL pmc_g_setname( me, myname%couple_index, myname%nameonchild )

 END SUBROUTINE pmc_c_set_dataarray_name


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Puts the value of -1 to the couple_index to mark the end of the array list.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_set_dataarray_name_lastentry( lastentry )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) ::  lastentry  !< marker for the end of the list of arrays 
!
!-- Local variables
    INTEGER ::  mype  !<

    TYPE(da_namedef) ::  myname  !<


    IF ( .NOT. lastentry )  RETURN

    myname%couple_index = -1

    IF ( m_model_rank == 0 )  THEN
       mype = MPI_ROOT
    ELSE
       mype = MPI_PROC_NULL
    ENDIF

    CALL pmc_bcast( myname%couple_index, mype, comm=m_to_parent_comm )

 END SUBROUTINE pmc_c_set_dataarray_name_lastentry


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Get the index list of parent PEs on the child PE.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_get_2d_index_list

    IMPLICIT NONE

    INTEGER(iwp) ::  arraysize  !< size of array to hold the index list
    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  i2         !<
    INTEGER(iwp) ::  ierr       !<
    INTEGER(iwp) ::  j          !<
    INTEGER(iwp) ::  nr         !<

    INTEGER, DIMENSION(me%inter_npes*2) ::  nrele  !< number of grid points in a horizontal slice of a 3d-array

    INTEGER, DIMENSION(:), POINTER ::  myind  !<

    TYPE(pedef), POINTER ::  ape  !> pointer to pedef structure


    DO  i = 1, me%inter_npes
       CALL pmc_recv_from_parent( nrele((i-1)*2+1:(i-1)*2+2), 2, i-1, 1001, ierr )
    ENDDO

!
!-- Allocate memory for index array.
    arraysize = 0
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       i2 = ( i-1 ) * 2 + 1
       nr = nrele(i2+1)
       IF ( nr > 0 )  THEN
          ALLOCATE( ape%locind(nr) )
       ELSE
          NULLIFY( ape%locind )
       ENDIF
       arraysize = MAX( nr, arraysize )
    ENDDO

    ALLOCATE( myind(2*arraysize) )

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       nr = nrele(i*2)
       IF ( nr > 0 )  THEN
          CALL pmc_recv_from_parent( myind ,2*nr, i-1, 1002, ierr )
          DO  j = 1, nr
             ape%locind(j)%i = myind(2*j-1)
             ape%locind(j)%j = myind(2*j)
          ENDDO
          ape%nrele = nr
       ELSE
          ape%nrele = -1
       ENDIF
    ENDDO

    DEALLOCATE( myind )

 END SUBROUTINE pmc_c_get_2d_index_list


 SUBROUTINE pmc_c_clear_next_array_list

    IMPLICIT NONE


!
!-- next_array_in_list is a global variable in pmc_child_mod.
    next_array_in_list = 0

 END SUBROUTINE pmc_c_clear_next_array_list



 LOGICAL FUNCTION pmc_c_getnextarray( myname )

    CHARACTER(LEN=*), INTENT(OUT) ::  myname  !<

    TYPE(pedef), POINTER    :: ape  !<
    TYPE(arraydef), POINTER :: ar   !<


!
!-- next_array_in_list is a global variable in pmc_child_mod.
    next_array_in_list = next_array_in_list + 1
!
!-- Array names are the same on all child PEs, so take first process to get the name.
    ape => me%pes(1)
!
!-- Check if all arrays have been processed.
    IF ( next_array_in_list > ape%nr_arrays )  THEN
       pmc_c_getnextarray = .FALSE.
       RETURN
    ENDIF

    ar => ape%array_list( next_array_in_list )

    myname = ar%name
!
!-- Return .TRUE. if another array follows.
!-- If all array have been processed, the RETURN statement a couple of lines above is executed.
    pmc_c_getnextarray = .TRUE.

 END FUNCTION pmc_c_getnextarray


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Set the array names of the partial copies of the parent arrays on actual child (2d-real arrays)
! Generic name: pmc_c_set_dataarray_2d.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_set_dataarray_2d( array )

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER ::  array  !<

    TYPE(C_PTR)             ::  array_adr  !<
    TYPE(arraydef), POINTER ::  ar         !<
    TYPE(pedef), POINTER    ::  ape        !<


    dims    = 1
    nrdims  = 2
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )
!
!-- Using C-pointer for storing the address of the array has among others the advantage that only
!-- one array_adr variable is necessary in the array list structure to handle 2d-real, 3d-real, and
!-- 2d-integer arrays.
    array_adr = C_LOC( array )
!
!-- Fill the array_list structure for every parent PE that is communicating with this child.
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims              ! number of dimensions of this array
       ar%dimkey = nrdims              ! key of dimension (used in pmc_c_getbuffer)
       ar%a_dim  = dims                ! array with sizes of the dimensions
       ar%data   = array_adr
    ENDDO

 END SUBROUTINE pmc_c_set_dataarray_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Set the array names of the partial copies of the parent arrays on actual child (2d-integer arrays)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_set_dataarray_ip2d( array )

    IMPLICIT NONE

    INTEGER(idp), INTENT(IN) , DIMENSION(:,:), POINTER ::  array  !<

    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), DIMENSION(4) ::  dims  !<

    TYPE(C_PTR)             ::  array_adr  !<
    TYPE(arraydef), POINTER ::  ar         !<
    TYPE(pedef), POINTER    ::  ape        !<


    dims    = 1
    nrdims  = 2
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )

    array_adr = C_LOC( array )
!
!-- Fill the array_list structure for every parent PE that is communicating with this child
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%dimkey = 22
       ar%a_dim  = dims
       ar%data   = array_adr
    ENDDO

 END SUBROUTINE pmc_c_set_dataarray_ip2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Set the array names of the partial copies of the parent arrays on actual child (3d-real arrays).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_set_dataarray_3d (array)

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  nrdims  !<

    INTEGER(iwp), DIMENSION (4) ::  dims  !<

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER ::  array  !<

    TYPE(C_PTR)             ::  array_adr  !<
    TYPE(pedef), POINTER    ::  ape        !<
    TYPE(arraydef), POINTER ::  ar         !<


    dims    = 1
    nrdims  = 3
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )
    dims(3) = SIZE( array, 3 )

    array_adr = C_LOC( array )
!
!-- Fill the array_list structure for every parent PE that is communicating with this child.
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%dimkey = nrdims
       ar%a_dim  = dims
       ar%data   = array_adr
    ENDDO

 END SUBROUTINE pmc_c_set_dataarray_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Receives the array starting indices from all the parent PEs (MPI_ALLTOALL) and determines the
!> transfer-buffer size for parent to child data transfer. Allocates the receive buffer and sets the
!> receive-buffer pointer to me%pes%array_list%recvbuf. Determines the array starting indices of the
!> of the child and sends them by MPI_ALLTOALL to parent PEs for child to parent data transfer.
!> Determines the buffer size for child to parent transfer and then allocates the transfer buffer
!> (pmc_alloc_mem) and creates the RMA-window for child to parent transfer and its base-array
!> pointer is assigned to the sendbuffer me%pes%array_list%sendbuf.
!>
!> Naming convention for appendices:  _pc  -> parent to child transfer
!>                                    _cp  -> child to parent transfer
!>                                    recv -> parent to child transfer
!>                                    send -> child to parent transfer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_setind_and_allocmem

    USE control_parameters,                                                                        &
        ONLY:  message_string

    IMPLICIT NONE

    INTEGER(iwp), PARAMETER ::  noindex = -1  !<

    INTEGER(iwp) ::  arlen           !<
    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  ierr            !<
    INTEGER(iwp) ::  j               !<
    INTEGER(iwp) ::  local_nr_arrays !<
    INTEGER(iwp) ::  myindex         !<
    INTEGER(iwp) ::  total_npes      !<

    INTEGER(idp) ::  bufsize  !< size of MPI data window

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_s  !< array of parent send indices to be received from parent
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_r  !< array of parent receive indices to be sent to parent


    TYPE(pedef), POINTER    ::  ape       !<
    TYPE(arraydef), POINTER ::  ar        !<

    Type(C_PTR)             ::  base_ptr  !<


    CALL MPI_COMM_SIZE( me%intra_comm, total_npes, ierr )

    local_nr_arrays = me%pes(1)%nr_arrays

    ALLOCATE( myindex_s(local_nr_arrays,0:total_npes-1) )
    ALLOCATE( myindex_r(local_nr_arrays,0:total_npes-1) )

    myindex_s = 0

!
!-- Receive indices from parent. The respective call on the parent side is in
!-- pmc_p_setind_and_allocmem.
    CALL MPI_ALLTOALL( myindex_s, local_nr_arrays, MPI_INTEGER, myindex_r, local_nr_arrays,        &
                       MPI_INTEGER, me%intra_comm, ierr )
!    
!-- An initial (base) buffer size is required because inactive processes need a non-zero size in 
!-- MPI_WIN_CREATE. Eight is just an arbitrary choice.
    bufsize = 8
!
!-- Parent to child direction.
!-- First stride: compute size and set index.
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvindex = myindex_r(j,i-1)
!
!--       Determine max, because child buffer is allocated only once.
!--       All 2d- and 3d-arrays use the same buffer.
          IF ( ar%nrdims == 3 )  THEN
             bufsize = MAX( bufsize, INT( ar%a_dim(1)*ar%a_dim(2)*ar%a_dim(3), MPI_ADDRESS_KIND ) )
          ELSE
             bufsize = MAX( bufsize, INT( ar%a_dim(1)*ar%a_dim(2), MPI_ADDRESS_KIND ) )
          ENDIF
       ENDDO
    ENDDO
!
!-- Create RMA (one sided communication) data buffer.
!-- The buffer for MPI_GET can be PE local, i.e. it can but must not be part of the MPI RMA window.
    CALL pmc_alloc_mem( base_array_pc, bufsize, base_ptr )
!
!-- Total buffer size in byte.
    me%totalbuffersize = bufsize * wp
!
!-- Second stride: set buffer pointer.
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvbuf = base_ptr
       ENDDO
    ENDDO

!
!-- Child to parent direction.
    myindex = 1
!    
!-- An initial (base) buffer size is required because inactive processes need a non-zero size in 
!-- MPI_WIN_CREATE. Eight is just an arbitrary choice.    
    bufsize = 8

    myindex_s = 0
    myindex_r = 0

    DO  i = 1, me%inter_npes

       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays

          ar => ape%array_list(j)
          IF ( ar%nrdims == 2 )  THEN
             arlen = ape%nrele
          ELSEIF( ar%nrdims == 3 )  THEN
             arlen = ape%nrele * ar%a_dim(1)
          ENDIF

          IF ( ape%nrele > 0 )  THEN
             ar%sendindex = myindex
          ELSE
             ar%sendindex = noindex
          ENDIF

          myindex_s(j,i-1) = ar%sendindex

          IF ( ape%nrele > 0 )  THEN
             ar%sendsize = arlen
             myindex     = myindex + arlen
             bufsize     = bufsize + arlen
          ENDIF

       ENDDO

    ENDDO
!
!-- Send indices to parent.
    CALL MPI_ALLTOALL( myindex_s, local_nr_arrays, MPI_INTEGER, myindex_r, local_nr_arrays,        &
                       MPI_INTEGER, me%intra_comm, ierr)

    DEALLOCATE( myindex_s )
    DEALLOCATE( myindex_r )

!
!-- Create RMA (one sided communication) window for data buffer child to parent transfer.
!-- The buffer of MPI_GET (counter part of transfer) can be PE-local, i.e. it can but must not be
!-- part of the MPI RMA window. Only one RMA window is required to prepare the data:
!--        for parent -> child transfer on the parent side
!-- and
!--        for child -> parent transfer on the child side.
    IF ( use_one_sided_communication )  THEN

       CALL pmc_alloc_mem( base_array_cp, bufsize )
!
!--    Total buffer size in byte.
       me%totalbuffersize = bufsize * wp

       winsize = me%totalbuffersize

       CALL MPI_WIN_CREATE( base_array_cp, winsize, wp, MPI_INFO_NULL, me%intra_comm,              &
                            me%win_parent_child, ierr )
    ELSE

       ALLOCATE( base_array_cp(bufsize) )

    ENDIF

    CALL MPI_BARRIER( me%intra_comm, ierr )
!
!-- Second stride: set buffer pointer.
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          IF ( ape%nrele > 0 )  THEN
             ar%sendbuf = C_LOC( base_array_cp(ar%sendindex) )
!
!--          Check the child buffer size.
             IF ( ar%sendindex+ar%sendsize > bufsize )  THEN
                WRITE( message_string, '(A,I4,4I7,1X,A)' ) 'child buffer too small ', i,           &
                                                           ar%sendindex, ar%sendsize,              &
                                                           ar%sendindex+ar%sendsize, bufsize,      &
                                                           TRIM( ar%name )
                CALL message( 'pmc_c_setind_and_allocmem', 'PMC0042', 3, 2, 0, 6, 0 )
             ENDIF
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE pmc_c_setind_and_allocmem


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Receive data from parent by MPI_GET. Or receive the data from parent using pmc_recv_from_parent
!> in case of two-sided communication.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_getbuffer( particle_transfer )

    IMPLICIT NONE

    INTEGER(iwp) ::  ierr     !<
    INTEGER(iwp) ::  ij       !<
    INTEGER(iwp) ::  ip       !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  myindex  !<
    INTEGER(iwp) ::  nr       !< number of elements to get from parent
    INTEGER(iwp) ::  nr_ar    !<
    INTEGER(iwp) ::  nr_req   !<


    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_disp  !<

    INTEGER(idp), POINTER, DIMENSION(:) ::  ibuf  !<

    INTEGER, DIMENSION(1) ::  buf_shape  !<

    INTEGER(iwp), DIMENSION(me%inter_npes) ::  req  !<

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  idata_2d  !<

    LOGICAL ::  particle_transfer_l  !<  local variable with default .FALSE. to hold value of optional argument particle_transfer

    LOGICAL, INTENT(IN), OPTIONAL ::  particle_transfer  !<

    REAL(wp), POINTER, DIMENSION(:)     ::  buf      !<
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(pedef), POINTER    ::  ape  !<
    TYPE(arraydef), POINTER ::  ar   !<

    TYPE ::  bufdef
       INTEGER(iwp)                            ::  nr    !<
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  buf   !<
       INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  ibuf  !<
    END TYPE

    TYPE(bufdef), DIMENSION(me%inter_npes) ::  bufall  !<


!
!-- Set local switch to determine if particle transfer is active.
    particle_transfer_l = .FALSE.
    IF ( PRESENT( particle_transfer) )  particle_transfer_l = particle_transfer

    IF ( use_one_sided_communication )  THEN
!
!--    Wait for buffer to be filled.
!--    The parent side (in pmc_p_fillbuffer) is filling the buffer in the MPI RMA window. When the
!--    filling is complete, an MPI_BARRIER is called. The child is not allowd to access the parent-
!--    buffer before it is completely filled. Therefore the following barrier is required.
       CALL MPI_BARRIER( me%intra_comm, ierr )

       DO  ip = 1, me%inter_npes
          ape => me%pes(ip)
          DO  j = 1, ape%nr_arrays
             ar => ape%array_list(j)

             IF ( ar%dimkey == 2  .AND.  .NOT.  particle_transfer_l )  THEN
                nr = ape%nrele
             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele * ar%a_dim(1)
             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
                nr = ape%nrele
             ELSE
!
!--             The value of particle_transfer_l is .T., if pmc_c_getbuffer is called from
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
!--          MPI passive target RMA. One data array is fetchted from MPI RMA window on parent.
             IF ( nr > 0 )  THEN
                target_disp = ar%recvindex - 1
                CALL MPI_WIN_LOCK( MPI_LOCK_SHARED , ip-1, 0, me%win_parent_child, ierr )
                IF ( particle_transfer_l )  THEN
!
!--                Here only arrays with dimkey=22 are treated.
                   CALL MPI_GET( ibuf, nr*8, MPI_BYTE, ip-1, target_disp, nr*8, MPI_BYTE,          &
                                 me%win_parent_child, ierr )
                ELSE
!
!--                Here only arrays with dimkey=2/3 are treated.
                   CALL MPI_GET( buf, nr, MPI_REAL, ip-1, target_disp, nr, MPI_REAL,               &
                                 me%win_parent_child, ierr )
                ENDIF
                CALL MPI_WIN_UNLOCK( ip-1, me%win_parent_child, ierr )
             ENDIF
             myindex = 1

             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   data_2d(ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex)
                   myindex = myindex + 1
                ENDDO

             ELSEIF ( ar%dimkey == 3  .AND.  .NOT.  particle_transfer_l )  THEN

                CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
                DO  ij = 1, ape%nrele
                   data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex:myindex+ar%a_dim(1)-1)
                   myindex = myindex+ar%a_dim(1)
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
!
!--    This barrier is because there is a corresponding barrier in pmc_p_fillbuffer where the parent
!--    waits until the child received all data.
!--    The barrier is just for security reasons. It could be removed.
       CALL MPI_BARRIER( me%intra_comm, ierr )

    ELSE
!
!--    For non blocking receive (MPI_IRECV) additional buffers are required to supply space for data
!--    of all receives. Data is completly available after pcm_waitall.
!--    Compute maximum number of arrays for all target processes (should be the same for all ip).
       nr_ar = 0
       DO  ip = 1, me%inter_npes
          ape => me%pes(ip)
          nr_ar = MAX(nr_ar,ape%nr_arrays)
       ENDDO
!
!--    For nonblocking receive ip and j loops have been swapped.
       DO  j = 1, nr_ar

          req    = 0
          nr_req = 0
          DO  ip = 1, me%inter_npes
             IF ( ALLOCATED( bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
             IF ( ALLOCATED( bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
          ENDDO
!
!--       First loop, initial non blocking receive.
          DO  ip = 1, me%inter_npes

             ape => me%pes(ip)
             ar  => ape%array_list(j)

             IF ( j > ape%nr_arrays )  CYCLE

             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele
             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
                nr = ape%nrele * ar%a_dim(1)
             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
                nr = ape%nrele
             ELSE
!
!--             The value of particle_transfer_l is .T., if pmc_c_getbuffer is called from
!--             pmc_particle_interface, and .F. if called from pmc_interface_mod. Depending on this,
!--             arrays that do or that do not belong to the LPM are skipped (not transferred).
                CYCLE
             ENDIF

             bufall(ip)%nr  = nr
             buf_shape(1) = nr

             IF ( nr > 0 )  THEN
                IF ( particle_transfer_l )  THEN
!
!--                Here only arrays with dimkey=22 are treated.
                   ALLOCATE( bufall(ip)%ibuf(nr) )
                   nr_req = nr_req+1
                   CALL pmc_recv_from_parent( bufall(ip)%ibuf, nr, ip-1, 2200+j, req=req(nr_req),  &
                                              ierr=ierr)
                ELSE
!
!--                Here only arrays with dimkey=2/3 are treated.
                   ALLOCATE( bufall(ip)%buf(nr) )
                   bufall(ip)%buf = 0.0
                   nr_req = nr_req+1
                   CALL pmc_recv_from_parent( bufall(ip)%buf, nr, ip-1, 2100+j, req=req(nr_req),   &
                                              ierr=ierr)
                ENDIF
             ENDIF

          ENDDO

          IF ( nr_req > 0 )  CALL pmc_waitall( nr_req, req, ierr )

!
!--       Second loop, get data from buffers.
          DO  ip = 1, me%inter_npes

             ape => me%pes(ip)
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
                   data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i) =                                  &
                                                      bufall(ip)%buf(myindex:myindex+ar%a_dim(1)-1)
                   myindex = myindex+ar%a_dim(1)
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

       DO  ip = 1, me%inter_npes
          IF( ALLOCATED( bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
          IF( ALLOCATED( bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
       ENDDO

    ENDIF

 END SUBROUTINE pmc_c_getbuffer


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fill buffer in RMA window to enable the parent to fetch the data with MPI_GET. Or send the data
!> to parent using pmc_send_to_parent in case of two-sided communication.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_putbuffer( particle_transfer )

    IMPLICIT NONE

    INTEGER(iwp) ::  ierr     !<
    INTEGER(iwp) ::  ij       !<
    INTEGER(iwp) ::  ip       !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  myindex  !<
    INTEGER(iwp) ::  nr_ar    !<
    INTEGER(iwp) ::  nr_req   !<


    INTEGER(iwp), DIMENSION(1) ::  buf_shape  !<

    INTEGER(idp), POINTER, DIMENSION(:) ::  ibuf  !<

    INTEGER(iwp), DIMENSION(me%inter_npes) :: req  !<

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  idata_2d  !<

    LOGICAL, INTENT(IN), OPTIONAL ::  particle_transfer  !<

    LOGICAL ::  particle_transfer_l  !<

    REAL(wp), POINTER, DIMENSION(:)     ::  buf  !<
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE  :: bufdef
       INTEGER(iwp)                            ::  nr    !<
       REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  buf   !<
       INTEGER(idp), ALLOCATABLE, DIMENSION(:) ::  ibuf  !<
    END TYPE

    TYPE(bufdef), DIMENSION(me%inter_npes) ::  bufall  !<

    TYPE(pedef),    POINTER ::  ape  !<
    TYPE(arraydef), POINTER ::  ar   !<


!
!-- Set local switch to determine if particle transfer is active.
    particle_transfer_l = .FALSE.
    IF ( PRESENT( particle_transfer) )  particle_transfer_l = particle_transfer

!
!-- Wait for empty buffer. Switch RMA epoche.
    IF ( use_one_sided_communication )  THEN

       DO  ip = 1, me%inter_npes

          ape => me%pes(ip)

          DO  j = 1, ape%nr_arrays
             ar => aPE%array_list(j)
             myindex = 1

             IF ( ar%dimkey == 2  .AND.  .NOT. particle_transfer_l )  THEN
!
!--             2d-REAL*8 array.
                buf_shape(1) = ape%nrele
                CALL C_F_POINTER( ar%sendbuf, buf,     buf_shape     )
                CALL C_F_POINTER( ar%data,    data_2d, ar%a_dim(1:2) )
                DO  ij = 1, ape%nrele
                   buf(myindex) = data_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                   myindex = myindex + 1
                ENDDO

             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
!
!--             3d-REAL*8 array.
                buf_shape(1) = ape%nrele*ar%a_dim(1)
                CALL C_F_POINTER( ar%sendbuf, buf,     buf_shape     )
                CALL C_F_POINTER( ar%data,    data_3d, ar%a_dim(1:3) )
                DO  ij = 1, ape%nrele
                   buf(myindex:myindex+ar%a_dim(1)-1) = data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i)
                   myindex = myindex + ar%a_dim(1)
                ENDDO

             ELSEIF ( ar%dimkey == 22  .AND.  particle_transfer_l )  THEN
!
!--             2d-INTEGER*8 array for particle transfer.
                buf_shape(1) = ape%nrele
                CALL C_F_POINTER( ar%sendbuf, ibuf,     buf_shape     )
                CALL C_F_POINTER( ar%data,    idata_2d, ar%a_dim(1:2) )

                DO  ij = 1, ape%nrele
                   ibuf(myindex) = idata_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                   myindex = myindex + 1
                ENDDO

             ENDIF

          ENDDO

       ENDDO

       CALL MPI_BARRIER( me%intra_comm, ierr )

       CALL MPI_BARRIER( me%intra_comm, ierr )

    ELSE
!
!--    Two-sided communication.       
       nr_ar = 0
       DO  ip = 1, me%inter_npes
          ape => me%pes(ip)
          nr_ar = MAX( nr_ar, ape%nr_arrays )
       ENDDO
!
!--    For non blocking send/receive swap ip and j loop.
       DO  j = 1, nr_ar

          req    = 0
          nr_req = 0
          DO  ip = 1, me%inter_npes
             IF ( ALLOCATED(bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
             IF ( ALLOCATED(bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
          ENDDO

          DO  ip = 1, me%inter_npes

             ape => me%pes(ip)
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
                   nr_req = nr_req+1
                   CALL pmc_send_to_parent( bufall(ip)%buf, ape%nrele, ip-1, 3100+j,               &
                                            req=req(nr_req), ierr = ierr)
                ENDIF

             ELSEIF ( ar%dimkey == 3  .AND.  .NOT. particle_transfer_l )  THEN
!
!--             3d-REAL*8 array.
                IF ( ape%nrele > 0 )  THEN
                   ALLOCATE( bufall(ip)%buf(ape%nrele*ar%a_dim(1)) )
                   CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
                   DO  ij = 1, ape%nrele
                      bufall(ip)%buf(myindex:myindex+ar%a_dim(1)-1) =                              &
                                                       data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i)
                      myindex = myindex + ar%a_dim(1)
                   ENDDO
                   nr_req = nr_req+1
                   CALL pmc_send_to_parent( bufall(ip)%buf, ape%nrele*ar%a_dim(1), ip-1, 3100+j,   &
                                            req=req(nr_req), ierr = ierr)
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
                   nr_req = nr_req+1
                   CALL pmc_send_to_parent( bufall(ip)%ibuf, ape%nrele, ip-1, 3200+j,              &
                                            req=req(nr_req), ierr = ierr)
                ENDIF

             ENDIF

          ENDDO

          IF ( nr_req > 0 )  CALL pmc_waitall( nr_req, req, ierr )

       ENDDO

       DO  ip = 1, me%inter_npes
          IF ( ALLOCATED(bufall(ip)%buf)  )  DEALLOCATE( bufall(ip)%buf  )
          IF ( ALLOCATED(bufall(ip)%ibuf) )  DEALLOCATE( bufall(ip)%ibuf )
       ENDDO

    ENDIF

 END SUBROUTINE pmc_c_putbuffer


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Free MPI window and respective memory of pmc data buffer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_c_finalize

    INTEGER(iwp) ::  ierr


    IF ( use_one_sided_communication )  THEN
       CALL MPI_FREE_MEM( base_array_cp, ierr )
       CALL MPI_WIN_FREE( me%win_parent_child, ierr )
    ENDIF

 END SUBROUTINE pmc_c_finalize
#endif

 END MODULE pmc_child
