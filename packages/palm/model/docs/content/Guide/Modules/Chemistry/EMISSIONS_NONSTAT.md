# Moving Emission Mode Input Data

For specific cases, non-stationary, i.e. moving, point-wise volume-source emissions can be explicitly provided. The spatially and temporally varying emissions are stored in a separate netCDF file in the `INPUT` directory, and are read automatically. The respective input file is indicated by the suffix `_emis_instat`. Point-wise, time-dependent volume-source emissions are read during the simulation and ascribed to the corresponding grid cells. Each point-source acts individually in space and time and can emit different species.

Please note, in order to avoid numerical implications caused by numerical dispersion and dissipation errors, point-sources should be resolved by at least three grid points along the x-, y- and z-axis (Ardeshiri et al., 2020).



Please note that all emissions are to be expressed in SI units, i.e. kg, mol, m, and s. Further, unless otherwise indicated, reactive gas phase emission species, such as NO<sub>2</sub> or O<sub>3</sub> are to be provided in moles. On the other hand, tracer or inert species, such as PM<sub>10</sub> or pollen (in the pollen module), are to be specified in kilograms.

## Temporal Emission Profiles

Temporal profiles of emission are stored using the following format:

`YYYY-MM-DD HH:mm:ss +ZZ`

where `YYYY-MM-DD` represents the date in full numeric format, `HH:mm:ss` represents the local time in 24-hour format, `ZZ` is the time zone. It is generally recommended that the profiles are specified in coordinated universal time (UTC) to maximize reusability. With this time format, the temporal resolution of the prescribed emissions is limited to 1s at the moment.

The user can specify timestamps with arbitrary start and end times, as well as intervals, in chronological order. The emission data is assigned onto the numerical grid according to the actual model time, while the simulation start is defined as by `origin_date_time` in the `_p3d` file). If the simulation time is in between two emission timesteps, both the spatial position and emission strength are linearly interpolated. With this approach, point-wise, moving sources such as from ships, airplanes, or cars, can be explicitly considered.

## netCDF File Format

Although the netCDF 4 API can be used to create the emission input file, storage of the emission data follows the netCDF 3 convention to maintain compatibility with other netCDF files used in the model. In particular, strings are to be stored as arrays of characters, as opposed to the `STRING` variable type. User-defined data structures are not used.

Each point-source acts individually in space and time and can emit different species. In order to consider different moving point sources, each point source is distinguished by an ID, being an integer number starting at one and added to the corresponding dimensions and variables. The total number of point sources must be defined by the global attribute `num_emission_path`, with `ID` running from 1 to `num_emission_path`.

### Dimensions

The definition of moving emission variables are predicated on the following mandatory dimensions:

| Dimension name | Value      | Description                           |
|----------------|------------|---------------------------------------|
| `ntime<ID>`    | at least 1 | Number of emissions timesteps         |
| `field_length` | 23         | Fixed length of all string variables  |
| `nspecies<ID>` | at least 1 | Number of emissions species           |
| `nvsrc<ID>`    | at least 1 | Number of volumetric emission sources |

Note that the `field_length` dimension must be set to 23. Other dimensions must be greater than zero.

### Variables

The mandatory variables for storage of emissions can then be defined using the above dimensions:

| Variable name    | netCDF data type | Dimension(s)                     | Description                                          |
|------------------|------------------|----------------------------------|------------------------------------------------------|
| `timestamp<ID>` | `NC_CHAR`        | (`ntime<ID>`, `field_length`)    | Individual time stamps for each set of emission data |
| `species<ID>`    | `NC_CHAR`        | (`nspecies<ID>`, `field_length`) | Names of individual chemical species                 |
| `vsrc<ID>_eutm`  | `NC_FLOAT`       | (`ntime<ID>`, `nvsrc<ID>`)       | EUTM coordinate                                      |
| `vsrc<ID>_nutm`  | `NC_FLOAT`       | (`ntime<ID>`, `nvsrc<ID>`)       | NUTM coordinate                                      |
| `vsrc<ID>_zag`   | `NC_FLOAT`       | (`ntime<ID>`, `nvsrc<ID>`)       | height above surface level                           |
| `vsrc_[species]` | `NC_FLOAT`       | (`ntime<ID>`, `nvsrc<ID>`)       | Volumetric emission of chemical species              |

The species indicated in the variable `species<ID>` can be of any name, but in general they should appear in the active kinetic mechanism used in the chemistry model. Species names not appeared in said mechanism and corresponding emissions will be ignored. All volumetric emission sources are expressed in terms of mol/(m<sup>3</sup>s) for gas-phase species and kg/(m<sup>3</sup>s) for particulate matter (PM).

### Example

The following is the header of an example netCDF emission input file for non-stationary moving emissions.

```
global attributes:
        num_emission_path = 2
```

```
dimensions:
        field_length = 23 ;
        ntime1 = 10 ;
        nvsrc1 = 27 ;
        nspecies1 = 2 ;

        ntime2 = 100 ;
        nvsrc2 = 9 ;
        nspecies2 = 1 ;

variables:
        char timestamp1(ntime1, field_length) ;
                timestamp1:description = "Time stamps" ;
        char species1(nspecies1, field_length) ;
                species1:description = "Emission species" ;
        int nspecies1(nspecies1) ;
                nspecies1:units = "-" ;
                nspecies1:standard_name = "number of species" ;
        int ntime1(ntime1) ;
                ntime1:units = "-" ;
                ntime1:standard_name = "number of timestamps" ;
        int nvsrc1(nvsrc1) ;
                nvsrc1:units = "-" ;
                nvsrc1:standard_name = "number of sources" ;
        float vsrc1_eutm(ntime1, nvsrc1) ;
                vsrc1_eutm:_FillValue = -9999.9f ;
        float vsrc1_nutm(ntime1, nvsrc1) ;
                vsrc1_nutm:_FillValue = -9999.9f ;
        float vsrc1_zag(ntime1, nvsrc1) ;
                vsrc1_zag:_FillValue = -9999.9f ;
        float vsrc1_PM2.5(ntime1, nvsrc1) ;
                vsrc1_PM2.5:_FillValue = -9999.9f ;
                vsrc1_PM2.5:unit = "kg/(m3 s)" ;
        float vsrc1_PM10(ntime1, nvsrc1) ;
                vsrc1_PM10:_FillValue = -9999.9f ;
                vsrc1_PM10:unit = "kg/(m3 s)" ;
        char timestamp2(ntime2, field_length) ;
                timestamp2:description = "Time stamps" ;
        char species2(nspecies2, field_length) ;
                species2:description = "Emission species" ;
        int nspecies2(nspecies2) ;
                nspecies2:units = "-" ;
                nspecies2:standard_name = "number of species" ;
        int ntime2(ntime2) ;
                ntime2:units = "-" ;
                ntime2:standard_name = "number of timestamps" ;
        int nvsrc2(nvsrc2) ;
                nvsrc2:units = "-" ;
                nvsrc2:standard_name = "number of sources" ;
        float vsrc2_eutm(ntime2, nvsrc2) ;
                vsrc2_eutm:_FillValue = -9999.9f ;
        float vsrc2_nutm(ntime2, nvsrc2) ;
                vsrc2_nutm:_FillValue = -9999.9f ;
        float vsrc2_zag(ntime2, nvsrc2) ;
                vsrc2_zag:_FillValue = -9999.9f ;
        float vsrc2_PM10(ntime2, nvsrc2) ;
                vsrc2_PM10:_FillValue = -9999.9f ;
                vsrc2_PM10:unit = "kg/(m3 s)" ;
```


## References

- Ardeshiri, H., M. Cassiani, S. Y. Park, A. Stohl, I. Pisso, A. S. Dinger (2020): On the convergence and capability of the large-eddy simulation of concentration fluctuations in passive plumes for a neutral boundary layer at infinite reynolds number, Boundary-Layer Meteol., 176, 291-327, DOI: https://doi.org/10.1007/s10546-020-00537-6
