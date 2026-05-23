---
title: User-defined output quantities
---
## User-defined output quantities

A very typical request is the calculation and output of quantities which are not provided via PALM's standard output. The basic user interface includes a number of subroutines which allow the calculation and output of user-defined quantities as:

- vertical profiles (horizontally averaged)
- time series
- 2d cross section, 3d volume data or masked data
- spectra
- virtual flights.

Examples for the generation of such output quantities are given and explained below, step by step. The respective user-interface subroutines already contain sample code lines (written as comment lines) for defining, calculating and output of quantities. Output times, averaging intervals, etc. are steered by the same variables as used for the standard PALM output quantities, e.g. [dt_data_output](../../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_data_output).


## Vertical profiles

This example explains the output of the quantity "turbulent resolved-scale horizontal momentum flux (*u\*v\**). If more than one user-defined quantity shall be output, the following steps have to be carried out in the same way for each of the quantities.

1. The quantity has to be given a unique string identifier, e.g. *u\*v\**. This identifier must be different from the identifiers used for the PALM standard output (see the [table of profile quantities](../../../../../Reference/LES_Model/Output_quantities/#vertical-profile-quantities)). The identifier must not contain more than **20** characters. In case that 2d cross section output is defined on one single level only (see the 2d section chapter, paragraph 9, further below), the identifier string must contain an asterisk (`"*"`). <br><br>To switch on output of the quantity, the string identifier has to be added to [data_output_pr_user](../../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_pr_user), eg.:

    `data_output_pr_user = 'u*v*'`,

2. For the quantity, an identification number, a physical unit, and the vertical grid on which it is defined (u- or w-grid), has to be assigned (subroutine [`user_check_data_output_pr`](../../../../../Reference/LES_Model/Modules/User_Interface#user_check_data_output_pr)):

        CASE ( 'u*v*' )

            user_pr_index = pr_palm + 1 ! identification number
            dopr_index(var_count) = user_pr_index
            dopr_unit = 'm2s2' ! physical unit
            unit = dopr_unit
            hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! vertical grid

    The identification number (user_pr_index) must be within the range [ `pr_palm+1 , pr_palm+max_pr_user` ], where `max_pr_user` is the number of user-defined profiles as given by parameter [data_output_pr_user](../../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_pr_user) in the respective PALM run. The physical unit has to be given with respect to the netCDF conventions. If no unit is given, PALM will abort. The vertical grid has to be either `zu` (`u`-grid) or `zw` (`w`-grid). 

3. The quantity has to be calculated for all grid points (subroutine [`user_statistics`](../../../../../Reference/LES_Model/Modules/User_Interface#user_statistics)): 
         
        !$OMP DO
        DO  i = nxl, nxr
            DO  j = nys, nyn
                DO  k = nzb+1, nzt

                    sums_l(k,pr_palm+1,tn) = sums_l(k,pr_palm+1,tn) + &
                    ( 0.5*(u(k, j, i)+u(k, j, i+1))-hom(k, 1, 1, sr))*&
                    ( 0.5*(v(k, j, i)+v(k, j+1, i))-hom(k, 1, 2, sr)) &
                    * rmask(j,i,sr)

                ENDDO
            ENDDO
        ENDDO

    The turbulent resolved-scale momentum flux *u\*v\** is defined as the product of the deviations of the horizontal velocities from their respective horizontally averaged mean values. These mean values are stored in array `hom(..,1,1,sr)` and `hom(..,1,2,sr)` for the u- and v-component, respectively. Since due to the staggered grid, *u* and *v* are not defined at the same grid points, they have to be interpolated appropriately (here to the center of the grid box). The result of the calculation is stored in `array sums_l`. The second index of this array is the identification number of the profile which must match the one given in the previous step 2.

    **Attention:** All quantities that are calculated here, need to be switched on via parameter `data_output_pr_user`, otherwise the run may crash, because the array `sums_l` is not allocated with its correct size.


## Time series

This example shows the output of two time series for the maxima of the absolute values of the horizontal velocities *u* and *v*. If more than one user-defined quantity shall be output, the following steps have to be carried out in the same way for each of the quantities.

1. For each time series quantity a label and a unit has to be given in subroutine [`user_check_data_output_ts`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_check_data_output_ts), which will be used for the netCDF file. They must not contain more than **13** characters. The value of `dots_num` and `dots_num_user` have to be increased by the number of new time series quantities. The old value of `dots_num` has to be stored in `dots_num_palm`:

        dots_num_palm = dots_num
        
        dots_num = dots_num + 1
        dots_num_user = dots_num_user + 1
        dots_label(dots_num) = 'abs_umax''
        dots_unit(dots_num) = 'm/s'
        
        dots_num = dots_num + 1
        dots_num_user = dots_num_user + 1
        dots_label(dots_num) = 'abs_vmax'
        dots_unit(dots_num) = 'm/s'
    
2. These quantities are calculated and output in subroutine [`user_statistics`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_statistics) for every defined statistic region `sr`, but at least for the "total domain" region: 

        ts_value_l(1) = ABS( u_max )
        ts_value_l(2) = ABS( v_max )
    
    Subsequently, values are collected by PE0, because only PE0 outputs the time series. However, collection includes taking the sum over all processors. The sum may have to normalized, depending on the calculated quantity. For serial runs, no action is required:

        #if defined( __parallel )

           IF ( collective_wait ) CALL MPI_BARRIER( comm2d, ierr )
           CALL MPI_ALLREDUCE( ts_value_l(1),ts_value(dots_num_palm+1,sr),dots_num_user, MPI_REAL, MPI_MAX, comm2d, ierr )

        #else

           ts_value(dots_num_palm+1:dots_num_palm+dots_num_user,sr) = ts_value_l

        #endif

    **HINT:** If the time series value that is calculated has the same value on all PEs, the `MPI_ALLREDUCE` call can be replaced by `ts_value(dots_num_palm+1:dots_num_palm+dots_num_user,sr) = ts_value_l`.


## 2d cross sections, 3d volume data or masked data

This example explains how to output the quantity "square of the u-component" (Note: this quantity could of course easily be calculated from the u-component by postprocessing the PALM output so that calculation within PALM is not necessarily required). If more than one user-defined quantity shall be output, the following steps have to be carried out in the same way for each of the quantities.

1. For output of **2d cross sections** and **3d volume data**, the quantity has to be given a unique string identifier, e.g. *'u2'*. This identifier must be different from the identifiers used for the PALM standard output (see the [table of 2d output quantities](../../../../../Reference/LES_Model/Output_quantities/#2d-array-quantities) or the [table of 3d output quantities](../../../../../Reference/LES_Model/Output_quantities/#3d-array-quantities)). The identifier must not contain more than **20** characters. To switch on output of this quantity, the string identifier has to be assigned to the parameter [data_output_user](../../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_user), eg.:

        data_output_user = 'u2', 'u2_xy_av'

    The pure string *'u2'* switches on the output of instantaneous 3d volume data. Output of cross section data and time averaged data is switched on by additionally appending the strings *'_xy', '_xz', '_yz'*, and/or *'_av'* (for a detailed explanation see parameter [data_output](../../../../../Reference/LES_Model/Namelists/#runtime_parameters--data_output)).

2. For output of **masked data**, the quantity has to be given a unique string identifier, e.g. *'u2'*. This identifier must be different from the identifiers used for the PALM standard output (see the [table of masked output quantities](../../../../../Reference/LES_Model/Output_quantities/#masked-array-quantities)). To switch on output of this quantity, the string identifier has to be assigned to the parameter [data_output_masks_user](../../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_masks_user), eg.:

        data_output_masks_user(1,:) = 'u2'
    
3. In order to store the quantities' data within PALM, a 3d data array has to be declared in module [`user`](../../../../../Reference/LES_Model/Modules/User_Interface#user): 

        REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: u2, u2_av

    The second array `u2_av` is required in case that output of time averaged data is requested. It is used to store the sum of the data of the respective time levels over which averaging is requested. For the output of **masked data**, the arrays must be declared as public. Otherwise, they are unknown quantities in [`user_data_output_mask`](../../../../../Reference/LES_Model/Modules/User_Interface#user_data_output_mask) and an error is issued.

4. The data array has to be allocated in subroutine [`user_init`](../../../../../Reference/LES_Model/Modules/User_Interface#user_init):

        ALLOCATE( u2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

5. The quantity has to be given a unit in subroutine [`user_check_data_output`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_check_data_output):

        CASE ( 'u2' )
            unit = 'm2/s2' 

    Otherwise, PALM will abort.

6. The vertical grid on which the quantity is defined (given by the levels `zu` or `zw`, on which the u- or w-component of the velocity are defined) has to be specified for the netCDF output files in subroutine [`user_define_netcdf_grid`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_define_netcdf_grid): 

        CASE ( 'u2', 'u2_xy', 'u2_xz', 'u2_yz' )
            grid = 'zu'

    As the example shows, this grid has to be defined for the 3d volume data as well as for all of the three cross sections.

7. After each time step, the quantity has to be calculated at all grid points in subroutine [`user_actions`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_actions) at location `'after_integration'`:

        CASE ( 'after_integration' )
        !
        !-- Enter actions to be done after every time integration (before data output).
        !-- Sample for user-defined output:
        DO  i = nxlg, nxrg
           DO  j = nysg, nyng
              DO  k = nzb, nzt+1
                 u2(k,j,i) = u(k,j,i)**2
              ENDDO
           ENDDO
        ENDDO

8. In case that output of time-averaged data is requested, the sum- and average-operations as well as the allocation of the sum-array have to be carried out in subroutine [`user_3d_data_averaging`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_3d_data_averaging):

        IF ( mode == 'allocate' ) THEN
           ...
            CASE ( 'u2' )
               IF ( .NOT. ALLOCATED( u2_av ) )  THEN
                  ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               ENDIF
               u2_av = 0.0
               ...

        ELSEIF ( mode == 'sum' ) THEN
           ...
           CASE ( 'u2' )
              IF ( ALLOCATED( u2_av ) )  THEN
                 DO  i = nxlg, nxrg
                    DO  j = nysg, nyng
                       DO  k = nzb, nzt+1
                          u2_av(k,j,i) = u2_av(k,j,i) + u2(k,j,i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              ...

        ELSEIF ( mode == 'average' ) THEN
           ...
           CASE ( 'u2' )
              IF ( ALLOCATED( u2_av ) )  THEN
                 DO  i = nxlg, nxrg
                    DO  j = nysg, nyng
                       DO  k = nzb, nzt+1
                          u2_av(k,j,i) = u2_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              ...

9. For output of **2d cross sections**, the data of the quantity has to be resorted to array `local_pf` in subroutine [`user_data_output_2d`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_data_output_2d). Also the vertical grid, on which the quantity is defined, has to be set again:

        CASE ( 'u2_xy', 'u2_xz', 'u2_yz' )
           IF ( av == 0 )  THEN

              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    DO  k = nzb_do, nzt_do
                       local_pf(i,j,k) = u2(k,j,i)
                    ENDDO
                 ENDDO
              ENDDO

            ELSE

               IF ( .NOT. ALLOCATED( u2_av ) )  THEN
                  ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                  u2_av = REAL( fill_value, KIND = wp )
               ENDIF
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb_do, nzt_do
                        local_pf(i,j,k) = u2_av(k,j,i)
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF
            grid = 'zu'

    The `ELSE` case is only required in case that output of time-averaged data is requested.

    As a special case, xy cross section output can also be defined at one single level at height `k=nzb+1` on the u-grid. This features is useful for output of surface data (e.g. heat fluxes). In this case, the corresponding 2d data has to be resorted to the array `local_pf(i,j,nzb+1)`. In addition to this, the grid defined in [`user_define_netcdf_grid`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_define_netcdf_grid) as well as in [`user_data_output_2d`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_data_output_2d) must be set to `grid = 'zu1'`. Furthermore, the identifier string must contain an asterisk (`'*'`).

        CASE ( 'u2*_xy' )
           IF ( av == 0 )  THEN
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    local_pf(i,j,nzb+1) = u2(j,i)
                 ENDDO
              ENDDO
           ELSE
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    local_pf(i,j,nzb+1) = u2_av(j,i)
                 ENDDO
              ENDDO
           ENDIF

           grid = 'zu1'
           two_d = .TRUE.

    Note that `two_d = .TRUE.` is required for output of just a 2d data section.

10. For output of 3d volume data, the data of the quantity has to be resorted to array local_pf in subroutine [`user_data_output_3d`](../../../../../Reference/LES_Model/Modules/User_Interface/#user_data_output_3d): 

        CASE ( 'u2' )
           IF ( av == 0 )  THEN
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    DO  k = nzb, nz_do
                       local_pf(i,j,k) = u2(k,j,i)
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              DO  i = nxlg, nxrg
                 DO  j = nysg, nyng
                    DO  k = nzb, nz_do
                       local_pf(i,j,k) = u2_av(k,j,i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF 

    The `ELSE` case is required only in case that output of time-averaged data is requested.

11. For output of masked data, the data of the quantity has to be resorted to array `local_pf` in subroutine [`user_data_output_mask`](../../../../../Reference/LES_Model/Modules/User_Interface#user_data_output_mask): 

        CASE ( 'u2' )
           IF ( av == 0 )  THEN
              DO  i = 1, mask_size_l(mid,1)
                 DO  j = 1, mask_size_l(mid,2)
                    DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = u2(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              DO  i = 1, mask_size_l(mid,1)
                 DO  j = 1, mask_size_l(mid,2)
                    DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = u2_av(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF

    The `ELSE` case is required only in case that output of time-averaged data is requested. 

12. In case of job chains, the sum array has to be output to the (binary) restart file (local filename [BINOUT](../../../../../Reference/LES_Model/IO-Files/#restart-output)) in subroutine [`user_last_actions`](../../../../../Reference/LES_Model/Modules/User_Interface#user_last_actions):

        IF ( ALLOCATED( u2_av ) )  THEN
           WRITE ( 14 ) 'u2_av'
           WRITE ( 14 ) u2_av
        ENDIF
    
    Otherwise, the time-average calculated in the next restart-run may be wrong.<br><br>In the restart run, this quantity has to be read from the restart file by the following code in subroutine [`user_read_restart_data`](../../../../../Reference/LES_Model/Modules/User_Interface#user_read_restart_data):

        IF ( initializing_actions == 'read_restart_data' )  THEN

           READ ( 13 ) field_char

           DO WHILE ( TRIM( field_char ) /= '*** end user ***' )

              DO  k = 1, overlap_count

                 nxlf = nxlfa(i,k)
                 nxlc = nxlfa(i,k) + offset_xa(i,k)
                 nxrf = nxrfa(i,k)
                 nxrc = nxrfa(i,k) + offset_xa(i,k)
                 nysf = nysfa(i,k)
                 nysc = nysfa(i,k) + offset_ya(i,k)
                 nynf = nynfa(i,k)
                 nync = nynfa(i,k) + offset_ya(i,k)

                 SELECT CASE ( TRIM( field_char ) )

                     CASE ( 'u2_av' )
                        IF ( .NOT. ALLOCATED( u2_av ) )  THEN
                           ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                        ENDIF
                        IF ( k == 1 )  READ ( 13 ) tmp_3d
                        u2_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                     CASE DEFAULT
                        WRITE( message_string, * ) 'unknown variable named "', &
                                                   TRIM( field_char ), '" found in', &
                                                   '&data from prior run on PE ', myid
                        CALL message( 'user_read_restart_data', 'UI0012', 1, 2, 0, 6, 0 )

                 END SELECT

              ENDDO

              READ ( 13 )  field_char

           ENDDO

        ENDIF


## Spectra

This example explains the output of spectra of the quantity "turbulent resolved-scale horizontal momentum flux" (*u\*v\**). If more than one user-defined quantity shall be output, the following steps have to be carried out in the same way for each of the quantities.

1. The calculation of user-defined spectra is closely linked with the calculation of [user-defined vertical profiles](#vertical-profiles) and [user-defined 3d volume data](#2d-cross-sections-3d-volume-data-or-masked-data). Therefore, the following prerequisites apply for each user-defined spectrum quantity:

    - From [user-defined vertical profiles](#vertical-profiles) steps 2 and 3. See the sample code (as comment lines) for *'u\*v\*'* and `ustvst`, respectively. (Actual output of vertical profiles - step 1 - is not required.)
    
    - From [user-defined 3d volume data](#2d-cross-sections-3d-volume-data-or-masked-data) steps 2, 3, 4, 5, and 6. See the sample code (as comment lines) for *'u\*v\*'* and `ustvst`, respectively. (Actual output of 3d volume data - step 1 - is not required.)
    
    - The quantity has to be given a unique string identifier, e.g. *'u\*v\*'*. This identifier must be different from the identifiers used for the PALM standard output (see [data_output_sp](../../../../../Reference/LES_Model/Namelists/#spectra_parameters--data_output_sp)). To switch on output of this quantity, the string identifier has to be added to parameter [data_output_sp](../../../../../Reference/LES_Model/Namelists/#spectra_parameters--data_output_sp), eg.: 

            data_output_sp = 'u*v*'

    All three items require a naming convention of identical identifiers, e.g. [data_output_pr_user](../../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_pr_user) = *'u\*v\*'*, [data_output_user](../../../../../Reference/LES_Model/Namelists/#user_parameters--data_output_user) = *'u\*v\*'*, and [data_output_sp](../../../../../Reference/LES_Model/Namelists/#spectra_parameters--data_output_sp) = *'u\*v\*'*. This naming convention applies only in case of user-defined spectra.

2. Edit and modify subroutine [`user_spectra`](../../../../../Reference/LES_Model/Modules/User_Interface#user_spectra), contained in file `user_spectra.f90`, as follows:

        IF ( mode == 'preprocess' )  THEN

           SELECT CASE ( TRIM( data_output_sp(m) ) )

              CASE ( 'u', 'v', 'w', 'pt', 'q' )
                 !-- Not allowed here since these are the standard quantities used in
                 !-- preprocess_spectra.

              CASE ( 'u*v*' )
                  pr = pr_palm+1
                  d(nzb+1:nzt, nys:nyn, nxl:nxr) = ustvst(nzb+1:nzt, nys:nyn, nxl:nxr)

              CASE DEFAULT
                  message_string = 'Spectra of ' // &
                                   TRIM( data_output_sp(m) ) // ' can not be calculated'
                  CALL message( 'user_spectra', 'USR0006', 0, 1, 0, 6, 0 )

           END SELECT

        ELSEIF ( mode == 'data_output' )  THEN

           SELECT CASE ( TRIM( data_output_sp(m) ) )

              CASE ( 'u', 'v', 'w', 'pt', 'q' )
                 !-- Not allowed here since these are the standard quantities used in
                 !-- data_output_spectra.

              CASE ( 'u*v*' )
                 pr = 6

              CASE DEFAULT
                 message_string = 'Spectra of ' // &
                                  TRIM( data_output_sp(m) ) // ' are not defined'
                 CALL message( 'user_spectra', 'USR0007', 0, 0, 0, 6, 0 )

           END SELECT

        ENDIF 

    Note that spectra output requires the additional namelist [&spectra_parameters](../../../../../Reference/LES_Model/Namelists/#spectra-parameters).


## Flight measurements

This example shows the output of two user-defined quantities for the absolute values of the horizontal velocities `u` and `v` captured by virtual flight measurements. The given quantities will be calculated and output for each leg.

1. At first, the number of user-defined quantities has to be given in the subroutine [`user_init_flight`](../../../../../Reference/LES_Model/Modules/User_Interface#user_init_flight) by 

        num_var_fl_user = num_var_fl_user + 2
    
    The subroutine is contained in file `virtual_flight_mod.f90`.

2. In the following, give a label and a unit for each user-defined quantity in subroutine [`user_init_flight`](../../../../../Reference/LES_Model/Modules/User_Interface#user_init_flight), which will be used for the netCDF file. They must not contain more than **13** characters. 

        CASE ( 1 )
           dofl_label(k) = TRIM(label_leg) // '_' // 'abs_u'
           dofl_unit(k) = 'm/s'
           k = k + 1

        CASE ( 2 )
           dofl_label(k) = TRIM(label_leg) // '_' // 'abs_v'
           dofl_unit(k) = 'm/s'
           k = k + 1

3. The user-defined quantities are calculated in subroutine [`user_flight`](../../../../../Reference/LES_Model/Modules/User_Interface#user_flight) at every timestep. Note, the aregument in the `CASE()` statement must be set accordingly to the settings used in [`user_init_flight`](../../../../../Reference/LES_Model/Modules/User_Interface#user_init_flight).

        CASE ( 1 )
           DO  i = nxl-1, nxr+1
              DO  j = nys-1, nyn+1
                 DO  k = nzb, nzt
                    var(k,j,i) = ABS( u(k,j,i) )
                 ENDDO
              ENDDO
           ENDDO

        CASE ( 2 )
           DO  i = nxl-1, nxr+1
              DO  j = nys-1, nyn+1
                 DO  k = nzb, nzt
                    var(k,j,i) = ABS( v(k,j,i) )
                 ENDDO
              ENDDO
           ENDDO

Flight measurements as well as data output of the respective user-defined quantities are done automatically.


## User-defined domains

By default, the values of the time series quantities and the horizontally averaged vertical profiles always refer to the total model domain. Independently, time series or profiles can be computed and output additionally for up to 9 different user-defined domains. In principle, steering is done via the initialization parameter [statistic_regions](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--statistic_regions).

These domains have to be defined within the user-defined routine [`user_init`](../../../../../Reference/LES_Model/Modules/User_Interface#user_init). The domains are defined via a mask array named `rmask`, which has to be set to value *1.0* for all horizontal grid points belonging to the user-defined domain, and to value *0.0* for those points not belonging to the user-defined domain. In the PALM code, `rmask` is declared as:

    REAL :: rmask(nysg:nyng,nxlg:nxrg,0:9)

The first two indices denote the array bounds (including the ghost points) in y and x-direction for each subdomain (don't confuse the subdomain with the user-defined domain!). The third index determines the user-defined domain, where *0* indicates the total model domain and *1* to *9* the user-defined domains.

The following example should illustrate the settings for two user-defined domains. The first domain is determined by all grid points which lie within a circle with center in the center of the model domain, and whose diameter is equal to half of the total horizontal domain size (a square-shaped total domain is assumed). The second domain is defined by all points outside of this domain. This requires the following lines of code in routine [`user_init`](../../../../../Reference/LES_Model/Modules/User_Interface#user_init):

    USE grid_variables
    USE indices
    USE statistics
    ...

    disc_center_x = dx * (nx + 1)/2
    disc_center_y = dy * (ny + 1)/2
    disc_radius = 0.5 * disc_center_x
    DO  i = nxlg, nxrg
       x = i * dx
       DO  j = nysg, nyng
          y = j * dy
          radial_distance = SQRT( ( x - disc_center_x )**2 + ( y - disc_center_y )**2 )
          IF ( radial_distance > disc_radius )  THEN
             rmask(j,i,1) = 0.0
             rmask(j,i,2) = 1.0
          ELSE
             rmask(j,i,1) = 1.0
             rmask(j,i,2) = 0.0
          ENDIF
       ENDDO
    ENDDO

The module `statistics` must be used here because it contains `rmask`. Likewise, the modules `grid_variables` and `indices` are required in this example because grid spacings and indices are used. All array elements of `rmask` (`rmask(:,:,:)`) are preset by the model with *1.0*. In no case this assignment must be changed for the total domain (`rmask(:,:,0)`)!

Computations and output for the user-defined domains only take place if [statistic_regions](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--statistic_regions) is set ≥ *1*.

Names for user-defined domains can be assigned via the initialization parameter [region](../../../../../Reference/LES_Model/Namelists/#user_parameters--region). Names of the selected user-defined domains are output to the _header and _rc files within the user-defined routine [`user_header`](../../../../../Reference/LES_Model/Modules/User_Interface#user_header).
