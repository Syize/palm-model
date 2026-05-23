---
title: Self nesting
---
# Self nesting

---

!!! warning
    This site is  Work in Progress.

    ToDo:

    - [ ] add description of particle transfer
    - [ ] add ocean-coupling related descriptions if necessary 

## General Features

See the [guide
section](../../../../../Guide/LES_Model/Modules/Nesting/self_nesting)
for general features and usage of the self-nesting.

## Nesting modes

Nesting data transfer between parent and child models can be made in two alternative modes:

- **One-way** means that children get solution data from their parents for setting
  boundary conditions on their nested boundaries but parents do not use
  any solution data from their children

- **Two-way** means that also parents use solution data received from their
  children to modify the parent solution in the area occupied by a
  child domain using so-called post-insertion procedure (Clark and Hall, 1991)

In case of two-way coupling and cascading arrangement of more than one
nested domain, the order of data-transfer operations (from parent to
child and from child to parent) can be made in three
alternative orders:

- **cascade**: first from parent to child consequently all the way down to the
  lowest level child and then from child to parent consequently all the way up to the root level
- **overlap**: simultaneous parent to child operations and child to parent operations
- **mixed** (default): simultaneous operations from parent to child and consequent operations from child to parent

Cascade is following the formally correct order of operations and it is
the slowest data-transfer mode because it involves more waiting time
than the other modes. Overlap is the fastest one, but it ignores the
information transferred on previous stages of the cascade. Mixed is
the default and the recommended data-transfer mode as it is faster
than cascade and does not ignore the data anterpolated on the previous
stages of the cascade. In practice, this is only important in the
child to parent transfer (anterpolation) in which 3-D data is
substituted to the parent domain while in the parent to child transfer
this is not critical because the data is only substituted on the child
boundaries which are located in different places on different-level nested
domains within the cascade.

Child boundary conditions on the nested boundaries are set by
interpolation which is described and discussed by Hellsten et
al. (2021). Also the child-data post-insertion to parent (called
anterpolation), and the phases of the time-step algorithm on which
these actions take place, as well as some other features are reported
in detail in Hellsten et al. (2021). It should be noted that the
interpolation and anterpolation operations are performed for the
intermediate velocity field, i.e. before the pressure-correction
step. The reason for this is that otherwise, a second
pressure-correction step would have been needed for the anterpolated
fields and this would be computationally very costly.


## Overall code arrangement

The nesting data-transfer pattern is complicated since both child and
parent models are internally parallellized and their subdomain
divisions are generally different and non-matching. Therefore, a
sending process typically has to send data to more than one receiving
processes, a receiving process has to receive data from more than one
sending processes, and subdomains may be only partially involved in
the transfer. This is illustrated in Figure 1 by a simple
example. This complicated data-transfer pattern is implemented in the
Palm Model Coupler (PMC). It is largely based on a hierarchy of
defined data types and creative use of pointers.

![parent_child_datatransfer_pattern](../../Figures/P_C_Mapping.png){width=80%}
<br> **Figure 1:** A simple example of the data transfer patterns
between parent and child on horizontal plane. Blue color represents
parent and red child. The designations like P0P and P0C stand for
process zero of parent and process zero of child, respectively, and so
on. The list on the left side shows which child processes each parent
process sends data, and the list on the right side indicates which
parent processes each child process sends data. Gridlines are drawn as
thin dotted lines and the thin solid lines are subdomain borders. The
red dashed rectangle shows the full horizontal extent of the child's
parent-grid arrays.

The data transfer is based on the parent grid. A parent only sends and
receives data while interpolation from parent to child grid and
anterpolation from child to parent grid are always made by child. This
significantly reduces the amount of data to be transferred. Children
have so-called child's parent-grid arrays which are defined in the
parent grid and span over the child-occupied subvolume of the parent
domain with two layers of ghost nodes beyond the outer boundaries but
not beyond the internal boundaries.

PMC consists of the PMC core including five modules, the PMC-interface
module, and the PMC-particle interface module. The PMC core is
connected to the rest of the code through the PMC interface only. The
PMC core modules do not access any other modules than other PMC
modules with the exception of module `kinds`, which includes the data
type kind-definitions. The PMC module source-code files are:

``` 
pmc_interface_mod.f90
pmc_particle_interface.f90 
    pmc_general_mod.f90
    pmc_handle_communicator_mod.f90 
    pmc_parent_mod.f90
    pmc_child_mod.f90 
    pmc_mpi_wrapper_mod.f90 
```

The PMC core has not been originally developed for PALM but for global
to regional climate- and weather forecast model coupling, and it has
been originally written in the C-language, and only later, during the
PALM-nesting development it was translated to Fortran. Therefore it
looks very different compared with other PALM source
code. Therefore the PMC source code is described here in more detail
than other parts of the PALM code. The interface modules `pmc_interface`
and `pmc_particle_interface` are developed specifically for PALM to
provide interface to the PMC core.

The `pmc_interface` is a relatively large module containing
subroutines for the upper-level nesting initiating and preparation
operations, all child-solution initialization operations,
interpolation and anterpolation operations as well as calling-routines
to the actual data transfer routines in the PMC core. All the
subroutine names in the PMC interface begin with `pmci_`. Similarly
all subroutines in `pmc_general` start with `pmc_g_`, in `pmc_parent`
with `pmc_p_` and in `pmc_child` with `pmc_c_`.


## Nesting initiation and initial preparations

**Initiation of a nested run.** In the beginning of the PALM main
program subroutine `pmci_init` in the `pmc_interface` is called.  This
subroutine calls `pmc_init_model` which resides in the
`pmc_handle_communicator` module, one of the five PMC core modules. By
calling `pmc_parin` this subroutine first detects if the
[&nesting_parameters](../../../Namelists/#nesting-parameters) namelist
is given in the root namelist file (suffix `_p3d`) or not. If it is not found,
`pmc_init_model` returns such that the run will continue as a
non-nested run. If it is found, `pmc_parin` reads the
`nesting_parameters` namelist and then `pmc_init_model` determines the
processing-element space for each model and creates the model-specific
communicators for each models' process groups according to which model
a process belongs. These process groups are said to have different
"colors".  This is done by splitting the global world communicator
`MPI_COMM_WORLD` by calling `MPI_COMM_SPLIT` according to couple id
`m_my_cpl_id` as the color and the global rank (the rank in `MPI_COMM_WORLD`)
`m_my_cpl_rank` of the current process.

Next, `pmc_init_model` creates inter-communicators for each
parent-child pair by combining the local communicators of different
colors by calling `MPI_INTERCOMM_CREATE`. This is done by using a copy
of `MPI_COMM_WORLD` as a peer communicator meaning that
`MPI_INTERCOMM_CREATE` is given the local communicator, the peer
communicator and the global rank of
the first process on the remote process group. This way it is able to
identify the remote group. `pmc_init_model` counts the number of
children of the current model and calls `MPI_INTERCOMM_CREATE` to
create the inter-communicators between the current model and all of
its children `m_to_child_comm(i)`, and between the current model and
its parent `m_to_parent_comm`.

The next call from the main program to PMC interface is to
`pmci_modelconfiguration`. It makes the initial preparations for all
models by calling other PMC-interface subroutines which e.g. set up
the coordinates of each model domain (`pmci_setup_coordinates`),
determines the number of variables to be coupled (`pmci_num_arrays`),
sets up children (`pmci_setup_child`) and parents
(`pmci_setup_parent`). The last two are large subroutines making calls
to several lower-level subroutines.


**Child setup.** The subroutine `pmci_setup_child` performs the
following operations. First it calls `pmc_c_childinit` which is in
module `pmc_child`.  It collects the information on the local
model-communicator `m_model_comm` and on the inter-communicator
`m_to_parent_comm` created by `pmc_init_model` into a `childdef`-type
variable `me`. The `childdef`-type is defined in the module
`pmc_general`. This information includes the local model communicator,
the inter-communicator and the rank of the current process and size of
the model communicator and the remote size of the
inter-communicator. Then it creates an intra-communicator from the
inter-communicator by calling `MPI_INTERCOMM_MERGE`. The reason why
the inter-communicators have to be merged to intra-communicators is
that the Remote Memory Access windows (RMA) used in the one-sided MPI
communication cannot employ inter-communicators, they must be based on
intra-communicators.  The intra-communicator `me%intra_comm`
includes first all the parent processes active for `me` (the current
child), and then the current child processes. Then it allocates a
pointer array `pes` for storing the child processes.  This array is an
element of the `childdef`-type variable `me` and is itself of type
`pedef`. The array `pes` in turn includes the type `arraydef` array
`array_list` which is also allocated here for each parent process
connected with the current child.  To understand the data structures,
it is important to get familiar with the data types `arraydef`,
`pedef`, and `childdef` defined in module `pmc_general`. Their
definitions are listed below.

``` Fortran

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

```

After returning from `pmc_c_childinit`, `pmci_setup_child` sets names
for the child's parent grid arrays by calling
`pmc_c_set_dataarray_name` in module `pmc_child`. It sets the name
information of each variable to be coupled in type `da_namedef`
structure `myname`. In the end of `pmc_c_set_dataarray_name`
subroutine `pmc_g_setname` of module `pmc_general` is called. It sets
`myname%couple_index` and `myname%nameonchild` to the `childdef`-type
data structure `me`.

Next, `pmci_setup_child` exchanges child- and parent-grid information
with its parent using the transfer subroutines `pmc_send_to_parent`
and `pmc_recv_from_parent` defined in the module
`pmc_mpi_wrapper_mod`.  After this
`pmci_map_child_grid_to_parent_grid` is called to determine the local
index bounds `ipl`, `ipr`, `jps`, and `jpn` of the data-exchange area
in the parent-grid horizontal index space. This is the area shown by
red dashed line in Figure 1. Next, the child receives the so-called
index list from its parent by calling `pmc_c_get_2d_index_list`. The
index list is described below in the context of parent setup. The next
task is to create the child's parent grid arrays by calling
`pmci_create_childs_parent_grid_arrays` sequently for each variable to
be coupled. The list of variables to be coupled is a linked-list type
construction controlled by function `pmc_c_getnextarray` and based on
the name information already set by
`pmc_c_set_dataarray_name`. Subroutine
`pmci_create_childs_parent_grid_arrays` allocates the requested
child's parent-grid array, for example `uc` for `u` and associates a
temporary pointer `p_3d` (or `p_2d` if the variable in question is a
2-D array or i_2d if it is an integer 2-D array). Then it passes this
temporary pointer to the subroutine `pmc_c_set_dataarray` which sets
this array pointer and the array dimensions inside the structure `me`
for every process on the parent side `i = 1, me%inter_npes`. The core of
this subroutine is given here as an example on how the data structures
are typically operated in the PMC.

```Fortran
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
```

Note that the type of `ar%data` is C-pointer. This is because the
one-sided MPI-communication operates on raw pointers (bare starting
memory addresses) and cannot thus handle Fortran-pointers. Therefore
the PMC code includes quite a lot of pointer conversions and C-pointer
returning memory allocations using `MPI_ALLOC_MEM` actually embedded
in the wrapper routine `pmc_alloc_mem` in the module
`pmc_mpi_wrapper_mod`. For this reason there is the USE-statement

```Fortran
    USE, INTRINSIC :: ISO_C_BINDING
```

in the declaration section of most of the PMC-modules. A C-pointer to
an array is simply the starting memory address of the array. This is
determined here by the function `C_LOC()`. The `childdef`-type
datastructure `me` includes the `pedef`-type structure `pes` which
contains elements for every parent process communicating with the
present child. It in turn contains the `arraydef` type `array_list`,
here shortnamed as `ar`, which contains all the arrays to be
coupled. This is the structure from which the final transfer buffer
will be formed in the actual child to parent data-transfer will be
extracted in `pmc_c_putbuffer`. This will be explained in section Data
transfer. Because `ar%data` is just a starting address, also the
dimensions have to be stored in the data structure.

The next phase is to prepare for the data receiving from parent that
will be done in `pmc_c_getbuffer` and data passing to parent that will
be done in `pmc_c_putbuffer`, see section Data transfer for more
details. These preparations are made by
`pmc_c_setind_and_allocmem`. It is important to understand that the
size of the RMA window, which is created by the parent, covers all the
processes on the receiving side (`i = 1, me%inter_npes`) and all the
coupled variables (`j = 1, me%pes(i)%nr_arrays`). However, the child
gets the data one by one in `pmc_c_getbuffer`. Therefore start indices
to the window are needed. These start indices are received from the
parent side by `MPI_ALLTOALL` in a two-dimensional array
`myindex_r`. The elements of `myindex_r` are then stored in
`me%pes(i)%array_list(j)%recvindex` in process- (`i`) and array- (`j`)
loops. In the same nested loops the receive buffer size is
determined. The receive buffer is not covering the whole RMA window,
instead it is only of the size of single child's parent grid array,
i.e. `ar%a_dim(1)*ar%a_dim(2)*ar%a_dim(3)`. Maximum size of all `i`-
and `j`-values is used for all of them. Next `base_array_pc`is
allocated by calling `pmc_alloc_mem`, which is a wrapper routine in
`module_pmc_mpi_wrapper_mod` module which calls `MPI_ALLOC_MEM` and
assigns the resulting C-pointer to a Fortran pointer by
`C_F_POINTER`. Next the C-pointer `base_ptr`to the receive buffer is
set for all `i` and `j` although it is the same for all of them
because only one buffer is allocated and it is reused for all transfer
iterations `i` and `j` in `pmc_c_getbuffer`. This completes the
preparations for data receiving from the parent on the child side.

Next the preparations for child to parent data transfer are carried
out in the latter part of `pmc_c_setind_and_allocmem`. First the
send-buffer size, `sendsize` (size of one transfer unit in the
`i-j`-loop) and the start indices `myindex_s` are determined in loops
as follows.

``` Fortran

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

```

The initial value of `bufsize` must not be zero, because RMA window
cannot be of zero size. Therefore an arbitrary small size of 8 is
initially given. The size of a horizontal layer of an individual array
to be transferred is given by `nrele` and its vertical dimension is
`ar%a_dim(1)`. Note that `nrele` may be smaller than the horizontal
size of a child's parent grid array on the sending process (which is
`(jpn-jps+1)*(ipr-ipl+1)`), because the receiving parent process
subdomain may overlap with the sending child subdomain only partially,
see Figure 1. To indicate a passive parent process (no overlap at all
and `nrele=0`) `noindex` is set to -1. The index `i-1`in
`myindex_s(j,i-1)` is because MPI process id-numbering always starts
from zero, but the process-loop `i` here starts from 1. After this,
the start indices are sent to parent by `MPI_ALLTOALL`. Next,
`base_array_cp` is allocated using `pmc_alloc_mem` and the RMA window
is created on the intra-communicator `me%intra_comm` by calling
`MPI_WIN_CREATE`. Finally, the buffer pointers are set as follows


``` Fortran

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          IF ( ape%nrele > 0 )  THEN
             ar%sendbuf = C_LOC( base_array_cp(ar%sendindex) )
             ...
          ENDIF
       ENDDO
    ENDDO

```

`C_LOC` returns C-pointer to the `ar%sendindex`-indicated entry in
`base_array_cp`. Here the three dots just indicate seven code lines
containing a buffer-size check and error-messaging. These code lines
are excluded from here for clarity and conciseness. This completes the
child-side preparations for child to parent transfer and the
subroutine `pmc_c_setind_and_allocmem`.

Next `pmci_setup_child` call its internal subroutine
`pmci_define_index_mapping` which defines the mappings from parent
grid to child grid. For example `iflo(ip)` and `ifuo(ip)`indicate the
lower and upper child grid points in the x-direction that correspond
to the parent grid point `ip`. Such mapping arrays are defined for all
three directions and separately for scalar nodes and staggered
velocity-component nodes in each direction.  These mappings are needed in the
interpolation and anterpolation. It also determines and saves the
numbers of child-grid nodes used for anterpolation for each
parent-grid node.

In the end of `pmci_setup_child` there are still few operations
left. First, the consistency of the nested-grid setup is checked by
`pmci_check_grid_matching`. Then the child domain boundary-face areas
are computed by `pmci_compute_face_areas` to be used in the
mass-conservation control in
`pmci_ensure_nest_mass_conservation`. Lastly, the anterpolation
starting height as a function of x and y is determined by
`pmci_compute_kpb_anterp`. This is needed for canopy-restricted
anterpolation.


**Parent setup.** The subroutine `pmci_setup_parent` performs largely
corresponding operations as `pmci_child_setup`, but there are
differences, too. One difference is that a parent may have several
children while a child can only have one parent. Therefore
`pmci_setup_parent` handles the parent-child pairs (called couplers)
in `childid`-loops over the number of children of the present parent.
Another difference is that on the parent side in
`pmci_set_array_pointer` the temporary pointer `p_3d` is set to the
whole array, for instance `u`, instead of a coupling array like for
instance `uc` on the child side. Such coupling arrays do not exist on
the parent side. Moreover, because the time level swapping of the
Runge-Kutta time-integration scheme, also another temporary pointer
`p_3d_sec` has to be set for the second time level variable, for
instance `u_2`. Also, on the parent side `dims` has four elements
instead of three like on the child side. The first three are the
dimensions of the actual prognostic array, for instance `u`, and the
fourth is the vertical dimension of the corresponding child's
parent-grid array, for example `uc`. Perhaps the largest difference to
`pmci_setup_child` is that parent creates the index list, prepares the
complete parent-child mapping and finally transfers the mapping
information to child.

The index list is an integer array which maps the child index- and
subdomain space on the parent index- and subdomain space. Its
dimension is `(6,index_list_size)` where the six elements of the first
dimension are `i`- and `j`-index on the parent grid, `i`- and
`j`-index on the child's parent grid arrays, child process id and
parent process id. The second dimension `index_list_size` is the
number of grid points on the parent-grid x-y-plane subspace covered by
the child domain including the ghost-node layers, see Figure 1. Note
that the index list contains only horizontal index-information. This
is a consequence of the two-dimensional domain decomposition. The only
vertical index information needed is the vertical dimension of the
transfer domains and this is provided by `dims(4)`.  Index list is
created by parent process 0 by subroutine
`pmci_create_index_list`. After creating the index list it calls
module `pmc_parent` subroutine `pmc_p_split_and_broadcast_index_list`
which splits the index list into local index lists according to the
parent process id of the index-list elements and sends each local
index list to the parent process it belongs to. Then it calls
`pmc_p_prepare_index_list_for_child`, which arranges the data
according to the child process id it will be sent to, defines the
start indices and number of elements to be sent to each child process
and the `i` and `j` indices of the child's parent grid array. Finally
it sends these information to child processes and `pmc_child` module
subroutine `pmc_c_get_2d_index_list` receives this information as
mentioned above in the description of `pmc_setup_child`, and
substitutes the number of elements in `me%pes(i)%nrele` and the `i`
and `j` grid indices in `me%pes(i)%locind(j)%i` and
`me%pes(i)%locind(j)%j`, respectively. Here the loop index `i` runs
over all the remote processes of the intercommunicator `me%inter_comm`
(usually many of them are not involved with the current transfer and
in those loop-iterations nothing is done). The loop index `j` runs
over the array elements on each child process (note that these `i` and
`j` are not the grid-indices). These actions complete the mapping
between the child index- and subdomain space on the parent index- and
subdomain space.

After setting up the parent-child mapping, the variables to be coupled
are set up by calling `pmci_set_array_pointer` to associate the
variables to be coupled to the temporary pointer as already mentioned
above. After that step, `pmc_p_setind_and_allocmem` is called for each
child of the current parent. It makes the preparations for the actual
data transfer which is made through the RMA-window that is created
here. These parent-side preparations are made largely in a similar
fashion as the corresponding operations on the child side in
`pmc_c_setind_and_allocmem` already described above.  Like in
`pmc_c_setind_and_allocmem`, here are two base arrays: `base_array_pc`
for parent to child transfer and `base_array_cp` for child to parent
transfer.  `pmc_p_setind_and_allocmem` first determines the start
indices `myindex_s` and the necessary size of `base_array_pc`. It
stores the start indices in
`children(childid)%pes(i)%array_list(j)%sendindex` according to the
process and variable. Then it sends the `myindex_s` to the child side
using `MPI_ALLTOALL` over the intra communicator
`children(childid)%intra_comm` which is the same communicator as
`me%intra_comm` on the child side. The corresponding receive has
already been discussed above in the description of `pmci_setup_child`.
After this, the base array is allocated and the RMA window is created
on `children(childid)%intra_comm` by calling `MPI_WIN_CREATE`. Finally
it sets the buffer pointer to the base array in
`children(childid)%pes(i)%array_list(j)%sendbuffer` which is of type
`C_PTR`. Note that `base_array_pc` and `base_array_cp` themselves are
Fortran pointers.  Below is shown the part of
`pmc_p_setind_and_allocmem` doing this. The three dots just indicate
five code lines containing a buffer-size check and
error-messaging. These code lines are excluded from here for clarity
and conciseness.


``` Fortran
    DO  i = 1, children(childid)%inter_npes
       ape => children(childid)%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%sendbuf = C_LOC( base_array_pc(ar%sendindex) )
          ...
       ENDDO
    ENDDO
```

After this point the corresponding operations are made for child to
parent transfer in a similar fashion as the parent to child transfer
preparations on the child side in `pmc_c_setind_and_allocmem`.


## Child initialization

Nested-domain solutions can be initialized by interpolating from an
existing root solution. The first-level child initial conditions can
be interpolated from the root solution (which is either from a precursor run or just the initial root setting). Then
the second-level child initial conditions can be interpolated from
their parents which are first level children, and so on until all
children have got their initial conditions. There are two alternative
methods for this initial-condition interpolation

- three-dimensional initialization (default)
- homogeneous initialization (see nesting parameter [homogeneous_initialization_child](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--homogeneous_initialization_child) )

In the `pmc_interface_mod`, there are subroutines:
`pmci_child_initialize` and `pmci_parent_initialize` called from the
main program before the beginning of the time
integration. `pmci_child_initialize` is called first to make sure that
models being both child and parent have first obtained their own initial
condition before transferring data further to their own children.

In case of three-dimensional initialization, `pmci_parent_initialize`
simply fills the buffers of the current parent by calling
`pmc_p_fillbuffer` for its children to receive the data by calling
`pmc_c_getbuffer` in `pmci_child_initialize`. The subroutines
`pmc_p_fillbuffer` and `pmc_p_fillbuffer` are described in section
[Data transfer](#data-transfer).  After this, `pmci_child_initialize` calls
`pmci_interp_all` and `exchange_horiz` for all variables to be
coupled.  The subroutine `pmci_interp_all` employs zeroth-order
interpolation which simply means that the parent-grid value in the
child's parent-grid array, for instance `uc`, is copied to all
child-grid points inside the current parent-grid cell.

In case of homogeneous initialization, `pmci_parent_initialize` first
calls `pmci_send_domain_averaged_profiles` in a loop for all children
of the current parent. The subroutine
`pmci_send_domain_averaged_profiles` calls `pmci_compute_average` to
compute horizontally averaged z-profile of the parent solution and
then sends it to child by calling `pmc_send_to_child`. The horizontal
averaging spans only the horizontal area occupied by the current
child. This process is repeated for all variables to be coupled except
`w` which is simply set to zero.  In `pmci_child_initialize`, the
vertically zeroth-order interpolated to the child grid.


## Data transfer

**Upper-level data transfer subroutines.** The upper-level data
transfer subroutines: `pmci_datatrans`, `pmci_child_datatrans`, and
`pmci_parent_datatrans` are found in the pmc interface module while
the lower-level data transfer subroutines are found in
`pmc_parent_mod` and `pmc_child_mod`. The hierarchy of these
subroutines, including also the upper-level interpolation and
anterpolation subroutines, is shown below.

```
pmci_datatrans
    pmci_parent_datatrans
        pmc_p_fillbuffer
        pmc_p_getbuffer
    pmci_child_datatrans
        pmc_c_getbuffer
        pmci_interpolation
        pmci_anterpolation
        pmc_c_putbuffer

``` 

The subroutine `pmci_datatrans` is called from the subroutine
`time_integration`. It branches according to the control parameters
`nesting_mode` and `nesting_datatransfer_mode` and calls
`pmci_parent_datatrans` and `pmci_child_datatrans` in a corresponding
order according to these control parameters. Subroutine
`pmci_parent_datatrans` calls the actual data transfer routine
`pmc_p_fillbuffer` or `pmc_p_getbuffer` depending on the transfer
direction (`parent_to_child` or `child_to parent`). After calling
`pmc_p_getbuffer`, it ensures as a safety measure that the velocity
components are zero inside buildings and under the terrain surface
after the anterpolation. In `parent_to_child`-direction subroutine
`pmci_child_datatrans` first calls `pmc_c_getbuffer` and then
`pmci_interpolation`. In `child_to_parent`-direction it first calls
`pmci_anterpolation` and then `pmc_c_putbuffer`.


**Parent send.** The subroutine `pmc_p_fillbuffer` in module
`child_mod` takes care of sending data to its children. It is called
separately for each child of the current parent and `childid` comes in
as a parameter. Then loops over the processes (`ip`) and arrays (`j`)
are run, and within them the transfer data buffer size is first
determined and stored in `buf_shape`. In case of normal
three-dmensional data of type `REAL(wp)` it is 
`children(childid)%pes(ip)%nrele * children(childid)%pes(ip)%array_list(j)%a_dim(4)`. 
The latter factor is the vertical dimension of the transfer domain, not the whole
vertical size of the parent domain. Next, the sendbuffer
`children(childid)%pes(ip)%array_list(j)%sendbuf` which is a C-pointer
is associated to the Fortran-pointer array `buf` by calling
`C_F_POINTER`. Similarly
`children(childid)%pes(ip)%array_list(j)%data` is associated with
`data_3d`.  Then the sendbuffer is filled through `buf`. See the truncated source
code insert below.

``` Fortran

    DO  ip = 1, children(childid)%inter_npes
        ape => children(childid)%pes(ip)
        DO  j = 1, ape%nr_arrays
            ar => ape%array_list(j)
            myindex = 1
            ...
            
            buf_shape(1) = ape%nrele*ar%a_dim(4)
            CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
            CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
            ...
!
!--         Copy from PALM 3d-array (e.g. u,v,...) to sendbuf in RMA window.
            DO  ij = 1, ape%nrele
                buf(myindex:myindex+ar%a_dim(4)-1) =                                            &
                              data_3d(ar%ks:ar%ke+2,ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + ar%a_dim(4)
            ENDDO
            ...
        ENDDO
    ENDDO
```

Note that the `ij`-loop runs over the horizontal dimensions of the
transfer domain (`nrele`) substituting the whole vertical stride of
length `ar%a_dim(1)` at each loop-iteration to `buf`. For those
`ip`-loop iterations for which `nrele` = 0, nothing is substituted to
`buf` and `myindex` is not incremented. After completing these loops
shown above, the sendbuffer is completely filled and there is a call
to `MPI_BARRIER` to make sure that the child side will not proceed to
receive the data until the buffer is completely filled. There is a
corresponding call to `MPI_BARRIER` on the child side in the beginning
of `pmc_c_getbuffer`. There is also second call to call to
`MPI_BARRIER` right after the first one to let the parent know that
the current child has completed its receive and the parent may go on
to the next call to send to its next child, if any. The counterpart of
this second `MPI_BARRIER` is in the end of `pmc_c_getbuffer`.

For setting `buf_shape`and for the `C_F_POINTER`-calls, there are also
corresponding branches for transfer of two-dimensional real-number
data and for two-dimensional integer-valued data. These are not shown
above but replaced by three dots for clarity and conciseness.

PMC also includes an option for using two-sided communication, but it
is not recommended in normal use as it is slower than one-sided
communication. However, note that one-sided communication is causing problems with some MPI-libraries, which can only be circumvented via setting of [use_one_sided_communication](../../../Namelists/#nesting_parameters--use_one_sided_communication) = *.F.*.  If two-sided communication is activated, the data is
substituted one variable (array) at time to `bufall(ip)` which is an
array of local defined type `bufdef` and `ip` is the child
target-process index as in the code inserts above. It is statically
allocated for all `children(childid)%inter_npes` child processes
involved and its elements contain transfer buffer `buf` for
real-valued data, and `ibuf` for integer-valued data, and the buffers
are dynamically re-allocated for each transfer of each variable. Here the
order of process and array loops is reversed, i.e. the array loop runs
outer, and the process loop inner. Data for each child target
processes is substituted in `buf` (or `ibuf`), and sent separately to
each receiving child process by calling `pmc_send_to_child`.


**Child receive.** The subroutine `pmc_c_getbuffer` in module
`child_mod` takes care of receiving data from the parent through the
RMA window which is filled by `pmc_p_fillbuffer`as described above. It
begins with a call to `MPI_BARRIER` which has its counterpart in
`pmc_p_fillbuffer` after the filling operations as described above to
make sure that the so-called epoch is not started before the parent
has completed its buffer filling. Then again loops over the processes
(`ip`) and arrays (`j`) are run, and within them the transfer data
buffer size is first determined and stored in `nr` which in case of
normal three-dimensional data of type `REAL(wp)` is `me%pes(ip)%nrele
* me%pes(ip)%array_list(j)%a_dim(1)`.  Then the actual transfer buffer
of type `C_PTR` is associated to the Fortran pointer `buf` by calling
`C_F_POINTER`. After this, if `nr` > 0, the subroutine starts the
epoch by calling `MPI_WIN_LOCK` and to get the data to `buf` by
calling `MPI_GET` followed by ending the epoch by calling
`MPI_WIN_UNLOCK`. This method of synchronization is known as passive
target synchronization. Next the data is substituted into the
three-dimensional array `data_3d`. See the truncated source code
insert below.

``` Fortran

    DO  ip = 1, me%inter_npes
        ape => me%pes(ip)
        DO  j = 1, ape%nr_arrays
            ar => ape%array_list(j)
            ...
            myindex = 1
            ...
            CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
            DO  ij = 1, ape%nrele
                data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex:myindex+ar%a_dim(1)-1)
                myindex = myindex+ar%a_dim(1)
            ENDDO
	    ...
        ENDDO
    ENDDO

```

Note that the `ij`-loop runs over the horizontal dimensions of the transfer
domain (`nrele`) substituting the whole vertical stride of length
`ar%a_dim(1)` at each iteration. This completes the array- and
process-loops and the whole receive sequence and finally
`MPI_BARRIER`is called again to inform the parent side that the
receive is completed. There is a corresponding call to `MPI_BARRIER`
on the parent side.

For the optional two-sided communication mode, `pmc_c_getbuffer` has a
receive buffer `bufall` similar to that on the parent side and
statically allocated for all involved `me%inter_npes` parent
processes. Its receive buffers `buf`and `ibuf`are dynamically
re-allocated for each transfer of each variable. Naturally the array and
process loops run in the same order as on the parent side, i.e. the
array loop runs outer. Data is received to `bufall(ip)%buf` or `ibuf` by
calling `pmc_recv_from_parent`. All the data for the current variable
is received to `bufall` within this process loop. After this, but
within the same array loop, there is a second process loop to
substitute the data from buffer to the three-dimensional array. This
procedure is relatively similar to that of the one-sided communication
branch shown above.

In child to parent transfer `pmc_c_putbuffer` sends data from child to
parent. It works very much the same way as `pmc_p_fillbuffer` and it
can be understood just by studying how `pmc_p_fillbuffer`
works. Therefore `pmc_c_putbuffer` is not further explained here.

Also `pmc_p_getbuffer` works similarly to `pmc_c_getbuffer` and can be
almost completely understood just by understanding
`pmc_c_putbuffer`. However, there is one exception worth to be
explained here. The range of child processes in the scope can be
limited to one by an optional argument `child_process_nr`. If this
optional argument is present, the process loop is restricted to only
this particular process indicated by `child_process_nr`. If this
argument is not present, the process loops run normally from 1 to
`children(childid)%inter_npes`. The reason for this arrangement is in
particle coupling. It is never used in coupling the Eulerian field
variables.


## Interpolation of child's boundary conditions

As soon as a child has received all the coupled variables in its
child's parent-grid arrays, interpolation of child's boundary
conditions on its nested boundaries starts as already explained in
section [Data transfer](#data-transfer). In the `pmc_interface_mod`, there is an
intermediate-level calling routine `pmci_interpolation` for calling
the actual interpolation routines. There are three interpolation
subroutines: `pmci_interp_lr` for handling the left and right
boundaries, `pmci_interp_sn` for handling the south and north
boundaries, and `pmci_interp_bt` for handling the bottom and top
boundaries. Originally, it was assumed that all child domains are set
on the bottom boundary of the root domain and it was not possible to
set the bottom boundary as nested, but an elevated-nest feature has
been implemented later and the bottom-boundary handling was added into
`pmci_interp_t` then renamed as `pmci_interp_bt`. The calling routine
`pmci_interpolation` first calls `pmci_interp_lr` for all necessary
variables if `nesting_bounds` is not set to 'vertical_only' or
'cyclic_along_x'. It is called first for all necessary variables to
set the boundary conditions on the left boundary and then for the
right boundary. Then, similar sequence of calls to `pmci_interp_sn` is
conducted to set the boundary conditions on the south boundary and
then for the north boundary if [nesting_bounds](../../../Namelists/#nesting_parameters--nesting_bounds) is not set to
*'vertical_only'* or *'cyclic_along_y'*. Finally, the boundary conditions
are set to top and bottom (if the nested domain is elevated)
boundaries by calling `pmci_interp_bt`.

The actual interpolation routines, for instance `pmci_interp_lr` start
by setting the necessary indices which depend on the variable in
question, staggered or non-staggered as projected on the boundary
plane, and on the boundary itself, in this case left or right.

Next, a so-called boundary work-array is updated. Such a boundary
work-array, in this case `work_arr_lr`, which extends by one
parent-grid cell past internal subdomain boundaries is needed because
the child's parent-grid arrays here designated as `parent_array` must
definitely not reach beyond any internal boundaries. Otherwise
anterpolation would be messed up around the internal boundaries. The
work arrays are first filled up from the `parent_array` and then its
ghost-nodes are exchanged by calling `MPI_SENDRECV`. These work
arrays: `work_arr_lr`, `work_arr_sn`, and `work_arr_bt` are allocated
by `pmci_allocate_workarrays` and their exchange data types to
facilitate swift exchange are defined by
`pmci_create_workarray_exchange_datatypes`. These two subroutines are
called from `pmci_define_index_mapping` which in turn is called by
`pmci_setup_child` in the initial preparation phase.

Then, the actual interpolation operations follow and there are
separate branches for all velocity components and for any scalar
variable. All the algorithms employed here are introduced and
explained by Hellsten et al. (2021) in detail. Therefore they are not
explained here.


## Anterpolation 

Anterpolation means to use the child solution for replacing the parent
solution within the subvolume occupied by the child domain (or a
subset of it). The child solution must first be restricted to the
parent grid. This is made by integrating the child solution such that
the velocity components are anterpolated onto two-dimensional faces
while scalars are anterpolated into three dimensional volumes (grid
cells), see (Clark and Hall, 1991) and (Hellsten at al., 2021). The
resulting fields are stored in child's parent-grid arrays. This
operation is done by subroutine `pmci_anterp_var` as described below.
After this phase is completed `pmc_c_child_datatrans` calls
`pmc_c_putbuffer` and the parent receives all the anterpolated
data directly in its variables by calling `pmc_p_getbuffer`.

There is a calling routine `pmci_anterpolation` called from
`pmci_child_datatrans`. This calling routine calls the actual
anterpolation routine `pmci_anterp_var` for all necessary variables
with the index-mapping arrays passed as arguments. Note that the
index-mappings are different for the velocity components (staggered)
and scalar variables (non-staggered). The number of
child-grid faces for each velocity component and child-grid cells for
scalars for all parent grid nodes are precomputed, stored in `ijkfc`,
and passed here as an argument. The anterpolation is just integration
(summation) over the child-grid nodes on the parent-grid face for
velocity components or in the parent-grid cell for scalars. The result
is stored in the corresponding child's parent-grid array. 

In practice, anterpolation cannot span over the whole child
domain. Narrow buffer-zones are defined just inside the nested
boundaries of the child domain. The default width of these buffers is
two parent-grid spacings. This can be increased by setting the
parameter [anterpolation_buffer_width](../../../Namelists/#nesting_parameters--anterpolation_buffer_width) in the
[&nesting_parameters](../../../Namelists/#nesting-parameters) namelist
to a value higher than two. The reason for these buffers is that the
nodes directly involved with the interpolation should not be
anterpolated because it would typically lead to an unstable feedback
loop and blow up of the solution.

In some cases, anterpolation may produce unphysical features in the
solution. For example if a canopy drag (plant or building canopy)
becomes remarkably different in child and parent due to the better
resolution of child, the mean-flow velocities in the parent may
increase or decrease in the area occupied by the child. This usually
induces unphysical secondary flow phenomena in the parent solution. In
an attempt to avoid such problems, anterpolation can be switched off
in the lowest part of the child domain by setting the parameter
[anterpolation_starting_height](../../../Namelists/#nesting_parameters--anterpolation_starting_height) in the
[&nesting_parameters](../../../Namelists/#nesting-parameters) namelist
to a suitable value. Note that this parameter is a grid index, not height
in metres, and it is relative to the local terrain surface. Note also,
that it influences all the domains the same way, i.e. it cannot be
specified separately for different domains. This feature has been
tested and studied by Hellsten et al. (2021).


## Mass conservation correction

The subroutine `time_integration` calls the
`pmc_interface_mod`-subroutine `pmci_ensure_nest_mass_conservation`
right before the pressure-correction step to make sure that overall
mass-conservation is not violated in any of the child domains.  This
is carried out by first integrating the total mass flux through all
the nested boundaries (including the bottom boundary in case of an elevated child). Then these are summed over the domain to
find out the possible mass imbalance. Then the correction velocity is
defined by dividing the imbalance by total boundary area of nested
boundaries. These areas are precalcuted within the initial
preparations. The correction velocity is then added onto the nested
boundaries with correct sign on each boundary.

## PMC finalization

In the end of a nested run, the main program calls `pmci_finalize`
which calls `pmc_p_finalize` and `pmc_c_finalize`. These subroutines call
`MPI_FREE_MEM` to release the memory allocated for the base arrays of
parent and child, respectively. Then they call `MPI_WIN_FREE` to free
the RMA windows.


## Instructions for adding a new variable to the coupling

If a new variable needs to be added to the set of variables to be
coupled or an existing one removed, it needs intervention in three
parts of the code. These are:

- in `pmci_setup_child`
- in `pmci_num_arrays`
- within the child initialization routines `pmci_child_initialize` and `pmci_send_domain_averaged_profiles`

First, in the part of `pmci_setup_child` where
`pmc_c_set_dataarray_name` is called for each coupled variable, one has to
add (or remove) `CALL pmc_c_set_dataarray_name('varname')`, where
`'varname'`is the name of the variable in question

In `pmci_num_arrays` one must ensure that `pmc_max_array`is set
equal to the actual number of variables to be coupled
including the added new variables as `pmc_max_array` is used to allocate
`me%pes(i)%array_list()`and `children(childid)%pes(j)%array_list()`.

Within `pmci_child_initialize`, one has to add (or remove) the
corresponding call to `pmci_interp_all` for three-dimensional
initialization. For homogeneous initialization one has to add (or
remove) the corresponding calls and `pmci_interp_1d` in
`pmci_child_initialize` and to `pmci_compute_average` and 
`pmc_send_to_child` in `pmci_send_domain_averaged_profiles`.



## References

**Clark, T. and Hall, W. (1991):** Multi-domain simulations of the time dependent
Navier-Stokes equations: benchmark error analysis of some nesting
procedures, J. Comput. Phys., 92, 456–481.

**Hellsten, A., Ketelsen, K., Suehring, M., Auvinen, M., Maronga, B., Knigge, C., Barmpas, F., Tsegas, G., and Moussiopoulos, N., Raasch, S. (2021):** A Nested Multi-Scale System Implemented in the Large-Eddy Simulation Model PALM model system 6.0. Geosci. Model Dev., 14(6), 3185-3214, [doi.org/10.5194/gmd-14-3185-2021](https://gmd.copernicus.org/articles/14/3185/2021/).
