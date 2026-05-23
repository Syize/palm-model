---
title: User Interfaces
---
## Interfaces for user-defined code

The following table lists and describes the available interfaces for user-defined code in the model, the names and possible arguments of the subroutines involved, as well as actions which can be accomplished in these subroutines. The respective name of a source code file is the subroutine name followed by `.f90`. Exception: the file containing the module `user` has the name `user_module.f90`.

| Call location      | Subroutine name and argument(s) | Purpose |
|--------------------|---------------------------------|------------------|
| ----------------------------------------------| <a id="user">`Module user`</a> | All user-defined variables which are used outside the respective local scopes of the user-defined subroutines have to be declared here. |
| `module_interface.f90`              | <a id="user_parin">`user_parin`</a>                                                                        |                  |
| `module_interface.f90`              | `user_check_parameters`                                                             |                  |
| `module_interface.f90`              | <a id="user_check_data_output_ts">`user_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )`</a>            |                  |
| `module_interface.f90`              | <a id="user_check_data_output_pr">`user_check_data_output_pr( variable, var_count, unit, dopr_unit )`</a>                |                  |
| `module_interface.f90`              | <a id="user_check_data_output">`user_check_data_output( variable, unit ) `</a>                                          |                  |
| `module_interface.f90`              | `user_init_arrays`                                                                  |                  |
| `module_interface.f90`              | <a id="user_init">`user_init`</a>                                                                         |                  |
| `module_interface.f90`              | <a id="user_header">`user_header ( io )`</a>                                                                 |                  |
| `module_interface.f90`              | <a id="user_actions">`user_actions( location )`</a><br> or<br> `user_actions( i, j, location )`                        |                  |
| `module_interface.f90`              | `user_exchange_horiz( location )`                                                   |                  |
| `module_interface.f90`              | `user_prognostic_equations`<br> or<br> `user_prognostic_equations( i, j, i_omp_start, tn )`   |                  |
| `module_interface.f90`              | `user_boundary_conditions`                                                          |                  |
| `module_interface.f90`              | `user_swap_timelevel( swap_mode )`                                                  |                  |
| `module_interface.f90`              | <a id="user_3d_data_averaging">`user_3d_data_averaging( mode, variable )`</a>                                          |                  |
| `module_interface.f90`              | <a id="user_data_output_2d">`user_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )`</a> |                  |
| `module_interface.f90`              | <a id="user_data_output_3d">`user_data_output_3d( av, variable, found, local_pf, resorted, nzb_do, nzt_do )`</a>    |                  |
| `module_interface.f90`              | <a id="user_statistics">`user_statistics( mode, sr, tn )`</a>                                                  |                  |
| `module_interface.f90`              | `user_rrd_global( found )`                                                          |                  |
| `module_interface.f90`              | `user_rrd_global`                                                                   |                  |
| `module_interface.f90`              | `user_wrd_global`                                                                   |                  |
| `module_interface.f90`              | `user_rrd_local( map_index, nxlf, nxlc, nxl_on_file ,`<br> `nxrf, nxrc, nxr_on_file,`<br> `nynf, nync, nyn_on_file,`<br> `nysf, nysc, nys_on_file,`<br> `tmp_3d, found )`,<br> or<br> `user_rrd_local `                                     |                  |
| `module_interface.f90`              | `user_wrd_local`                                                                    |                  |
| `module_interface.f90`              | <a id="user_last_actions">`user_last_actions`</a>                                                                 |                  |
| `data_output_mask.f90`              | <a id="user_data_output_mask">`user_data_output_mask( av, domask(mid,av,ivar), found, local_pf, mid )`</a>            |                  |
| `data_output_spectra.f90` <br>and<br> `spectra_mod.f90`          | <a id="user_spectra">`user_spectra( 'data_output', m, pr )`</a>                                              |                  |
| `init_3d_model.f90`                 | `user_init_3d_model`                                                                |                  |
| `lagrangian_particle_model_mod.f90` | `user_lpm_init`                                                                     |                  |
| `lagrangian_particle_model_mod.f90` | `user_lpm_advec( i, j, k )`                                                         |                  |
| `land_surface_model_mod.f90`        | `user_init_land_surface`                                                           |                  |
| `netcdf_interface_mod.f90`          | <a id="user_define_netcdf_grid">`user_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )`</a><br> or<br> `user_define_netcdf_grid( do3d(av,i), found, grid_x, grid_y, grid_z )`<br> or<br> `user_define_netcdf_grid( data_output_sp(i), found, grid_x, grid_y, grid_z )`          |                  |
| `plant_canopy_model_mod.f90`        | `user_init_plant_canopy`                                                            |                  |
| `radiation_model_mod.f90`           | `user_init_radiation`                                                               |                  |
| `topography_mod.f90`                | `user_init_grid( topo )`                                                            |                  |
| `urban_surface_mod.f90`             | `user_init_urban_surface`                                                           |                  |
| `virtual_flight_mod.f90`            | <a id="user_init_flight">`user_init_flight( init )`</a><br> or<br> `user_init_flight( init, k, i, label_leg )`             |                  |
| `virtual_flight_mod.f90`            | <a id="user_flight">`user_flight( var_u, n )`</a>                                                           |                  |
