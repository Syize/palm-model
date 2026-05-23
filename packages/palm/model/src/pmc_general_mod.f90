!> @file pmc_general_mod.f90
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
! Authors:
! --------
!> @author Klaus Ketelsen (no affiliation)
!
! Description:
! ------------
!> Structure definition and utilities of the Palm Model Coupler.
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_general

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

    USE kinds

    USE MPI

    IMPLICIT NONE

    INTEGER(iwp) ::  pmc_max_array  !< max # of arrays which can be coupled
                                    !< - will be determined dynamically in pmc_interface

    INTEGER(iwp), PARAMETER ::  da_desclen       =  8  !<
    INTEGER(iwp), PARAMETER ::  da_namelen       = 16  !<
    INTEGER(iwp), PARAMETER ::  pmc_da_name_err  = 10  !<
    INTEGER(iwp), PARAMETER ::  pmc_max_models   = 64  !<
    INTEGER(iwp), PARAMETER ::  pmc_status_ok    =  0  !<

    TYPE ::  ij_index  !< pair of indices in horizontal plane
       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
    END TYPE

    TYPE ::  arraydef
       CHARACTER(LEN=da_namelen) ::  name  !< name of array

       INTEGER(iwp) ::  coupleindex  !< ID of array
       INTEGER(iwp) ::  dimkey       !< key for NR dimensions and array type (2 = 2d-real, 3 = 3d-real, or 22  = 2d-integer*8))
       INTEGER(iwp) ::  nrdims       !< number of dimensions
       INTEGER(iwp) ::  recvsize     !< size in receive buffer
       INTEGER(iwp) ::  sendsize     !< size in send buffer
       INTEGER(idp) ::  recvindex    !< index in receive buffer
       INTEGER(idp) ::  sendindex    !< index in send buffer
       INTEGER(iwp) ::  ks           !< start index in z direction (3d arrays on parent only)
       INTEGER(iwp) ::  ke           !< end   index in z direction (3d arrays on parent only)

       INTEGER(iwp), DIMENSION(4) ::  a_dim  !< size of dimensions

       TYPE(C_PTR) ::  data     !< pointer of data in parent space
                                !< c pointers are used because they can point to different data types
                                !< e.g. REAL(2d), REAL(3d), INTEGER, etc.
                                !< otherwise, separate arraydef types would be required for each of them
       TYPE(C_PTR) ::  sendbuf  !< data pointer in send buffer
       TYPE(C_PTR) ::  recvbuf  !< data pointer in receive buffer

       TYPE(arraydef), POINTER ::  next  !<

       TYPE(C_PTR), DIMENSION(2) ::  po_data  !< base pointers, pmc_p_set_active_data_array
                                              !< sets active pointer to respective time level (e.g. u_1 or u_2)
    END TYPE arraydef

    TYPE ::  pedef
       INTEGER(iwp) ::  nr_arrays = 0  !< number of arrays which will be transfered
       INTEGER(iwp) ::  nrele          !< number of elements along x and y, same for all arrays

       TYPE(arraydef), POINTER, DIMENSION(:) ::  array_list  !< list of data arrays to be transfered

       TYPE(ij_index), POINTER, DIMENSION(:) ::  locind  !< i,j index local array for remote PE
    END TYPE pedef

    TYPE ::  childdef
       INTEGER(iwp) ::  inter_comm        !< inter communicator model and child
       INTEGER(iwp) ::  inter_npes        !< number of PEs child model
       INTEGER(iwp) ::  intra_comm        !< intra communicator model and child
       INTEGER(iwp) ::  intra_rank        !< rank within intra_comm
       INTEGER(iwp) ::  model_comm        !< communicator of this model
       INTEGER(iwp) ::  model_npes        !< number of PEs this model
       INTEGER(iwp) ::  model_rank        !< rank of this model
       INTEGER(idp) ::  totalbuffersize   !< size of RMA window
       INTEGER(iwp) ::  win_parent_child  !< MPI RMA window object for preparing data on parent AND child side

       TYPE(pedef), DIMENSION(:), POINTER ::  pes  !< list of all child PEs on parent or list of all parents on child
    END TYPE childdef

!
!-- In the following data type, only the couple_index and one array name would be sufficient,
!-- because arrays always have the same name on parent and child. childdesc and parentdesc could
!-- be removed.
    TYPE ::  da_namedef  !< data array name definition
       CHARACTER(LEN=da_desclen) ::  childdesc     !< child array description
       CHARACTER(LEN=da_namelen) ::  nameonchild   !< name of array within child
       CHARACTER(LEN=da_namelen) ::  nameonparent  !< name of array within parent
       CHARACTER(LEN=da_desclen) ::  parentdesc    !< parent array description

       INTEGER(iwp) ::  couple_index  !< unique number of array
    END TYPE da_namedef

    SAVE

    PRIVATE

!
!-- Public functions.
    PUBLIC pmc_g_setname

!
!-- Public variables, constants and types.
    PUBLIC arraydef,                                                                               &
           childdef,                                                                               &
           da_desclen,                                                                             &
           da_namedef,                                                                             &
           da_namelen,                                                                             &
           pedef,                                                                                  &
           pmc_da_name_err,                                                                        &
           pmc_max_array,                                                                          &
           pmc_max_models,                                                                         &
           pmc_status_ok

    INTERFACE pmc_g_setname
       MODULE PROCEDURE pmc_g_setname
    END INTERFACE pmc_g_setname


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Increment the number of array, set array name and couple_id to the array_list in the
!> structure PEs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_g_setname( child, couple_index, aname )

    CHARACTER(LEN=*), INTENT(IN) ::  aname  !<

    INTEGER(iwp) ::  i  !<

    INTEGER(iwp), INTENT(IN) ::  couple_index  !< ID of array

    TYPE(childdef), INTENT(INOUT) ::  child  !< an element of the child data structure "children"
                                             !< this subroutine is called from a children loop

    TYPE(pedef), POINTER ::  ape  !<


!
!-- Assign the new array to next free index in the array list. Set name of array in arraydef
!-- structure.
    DO  i = 1, child%inter_npes
       ape => child%pes(i)
       ape%nr_arrays = ape%nr_arrays + 1
       ape%array_list(ape%nr_arrays)%name        = aname
       ape%array_list(ape%nr_arrays)%coupleindex = couple_index
    ENDDO

 END SUBROUTINE pmc_g_setname
#endif


 END MODULE pmc_general
