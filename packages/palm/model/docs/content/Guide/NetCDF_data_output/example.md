---
title: Example of a PALM netCDF data set
---

# Example of a PALM netCDF data set
---

The netCDF dataset described here contains data of instantaneous horizontal cross sections and has been created using the settings of the example parameter file (see the ... (**link to (new) turorial "My first LES (convective boundary layer) needs to be added here**) ), i.e. it contains section data of the w-velocity-component and of the potential temperature for vertical grid levels with index `k = 2` and `k = 10`, selected by the respective parameter settings [data_output](../../../../Reference/LES_Model/Namelists/#runtime_parameters--data_output) = *'w_xy'*, *'pt_xy'*, and [section_xy](../../../../Reference/LES_Model/Namelists/#runtime_parameters--section_xy) = *2, 10*. Output has been created after every 900 s ([dt_data_output](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_data_output) = *900.0*). Because of [end_time](../../../../Reference/LES_Model/Namelists/#runtime_parameters--end_time) = *3600.0*, the file contains data of 4 time levels (*t* = *900, 1800, 2700, 3600* s). Supposed that the name of the netCDF dataset is `example_cbl_xy.nc`, the general content of this file (file header) as generated via command 
```
ncdump -h example_cbl_xy.nc
```
is as follows.
The original `ncdump` output is displayed using monospace, additional explanations are given behind exclamation marks. 

``` fortran
netcdf example_cbl_xy {                          ! filename
dimensions:                                      ! 41 gridpoints along x and y, 4 timelevels
        time = UNLIMITED ; // (4 currently)      ! unlimited means that additional time levels can be added (e.g. by
                                                 ! restart jobs)
        zu_xy = 2 ;                              ! vertical dimension (2, because two cross sections are selected);
        zw_xy = 2 ;                              ! there are two different vertical dimensions zu and zw because due
        zu1_xy = 1 ;                             ! to the staggered grid the z-levels of variables are those of the
        x = 41 ;                                 ! u- or the w-component of the velocity
        xu = 41 ;
        y = 41 ;
        yv = 41 ;
variables:                                       ! precision, dimensions, and units of the variables
        double time(time) ;                      ! the variables containing the time levels and grid point co-
                time:units = "seconds" ;         ! ordinates have the same names as the respective dimensions
        double zu_xy(zu_xy) ;
                zu_xy:units = "meters" ;
        double zw_xy(zw_xy) ;
                zw_xy:units = "meters" ;
        double zu1_xy(zu1_xy) ;
                zu1_xy:units = "meters" ;
        double ind_z_xy(zu_xy) ;
                ind_z_xy:units = "gridpoints" ;
        double x(x) ;
                x:units = "meters" ;
        double xu(xu) ;
                xu:units = "meters" ;
        double y(y) ;
                y:units = "meters" ;
        double yv(yv) ;
                yv:units = "meters" ;
        float w_xy(time, zw_xy, y, x) ;          ! array of the vertical velocity; it has 4 dimensions: x and y,
                w_xy:long_name = "w_xy" ;        ! because it is a horizontal cross section, zw_xy, which defines
                w_xy:units = "m/s" ;             ! the vertical levels of the sections, and time, for the time levels
        float pt_xy(time, zu_xy, y, x) ;         ! array of the potential temperature, which is defined on the u-grid
                pt_xy:long_name = "pt_xy" ;    
                pt_xy:units = "K" ;

// global attributes:
                :Conventions = "COARDS" ;
                :title = "PALM 3.7a  Rev: 584  run: example_cbl.00  host: lcsgih  12-10-10 08:52:28" ; ! PALM run-
                                                 ! identifier
                :VAR_LIST = ";w_xy;pt_xy;" ;     ! the list of output quantities contained in this dataset;
                                                 ! this global attribute can be used by FORTRAN programs to identify
                                                 ! and read the quantities contained in the file

data:

 time = 905.7, 1800.46, 2711.96, 3610.86 ;       ! values of the four time levels

 zu_xy = 75, 475 ;                               ! heights of the two selected cross sections (u-grid)

 zw_xy = 100, 500 ;

 zu1_xy = 25 ;

 x = 25, 75, 125, 175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, ! x-coordinates of the grid points
    775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375,
    1425, 1475, 1525, 1575, 1625, 1675, 1725, 1775, 1825, 1875, 1925, 1975,
    2025 ;

 xu = 0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
    750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
    1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950,
    2000 ;

 y = 25, 75, 125, 175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725,  ! y-coordinates of the grid points
    775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375,
    1425, 1475, 1525, 1575, 1625, 1675, 1725, 1775, 1825, 1875, 1925, 1975,
    2025 ;

 yv = 0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
    750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
    1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950,
    2000 ;
}
```

If the option `-h` is omitted in the `ncdump` call, then all values of all variables (here: grid point data of `w_xy` and `pt_xy`) are output to the terminal. 

```
ncdump -v pt_xy example_cbl_xy.nc
```
will display only the grid point data of the quantity specified via option `-v`.
