# SLUrb user guide

## When to use SLUrb?
SLUrb is a single-layer urban canopy representation aimed for simulations where grid resolution is not sufficient to resolve individual obstacles (e.g. buildings) and the flow around them. It can be used in e.g. mesoscale studies, in coarse resolution domains of nested setups, or in any other context where modelling the processes within the urban canopy itself with LES is not necessary. The required resolution to resolve the urban canopies with LES depends heavily on meteorological conditions and the urban form itself. In general, SLUrb is intended to be used with grid resolutions in the order of ten meters ($\geq 10~\mathrm{m}$).

## Job preparation
SLUrb is enabled by adding `&slurb_parameters` into p3d namelist, which is also used to define the general configuration for the model. In addition, a minimal set of model surface parameters can be set for testing purposes using only the namelist parameters (homogeneous initialization). Spatially heterogeneous initialization of surface parameters is realized through a netCDF driver, see [SLUrb's model driver documentation](../../../../../Reference/LES_Model/IO-Files/PIDS_SLURB) for details. The value set using netCDF input takes priviledge over the namelist initialization. Thus, the namelist-based initialization may be amended by the netCDF input only for some grid cells.

If all the mandatory inputs are provided using the netCDF driver and the default configuration is otherwise suitable for your case, [`slurb_parameters`](../../../../../Reference/LES_Model/Namelists/#slurb-parameters) must still be included in the namelist to enable the module, but can be left empty:

```
&slurb_parameters
/
```

### Namelist example

A minimum setup with homogeneous surface and anisotropic south-north oriented street canyons could look like the following:

```
&slurb_parameters
    urban_fraction = 0.5,
    urban_roughness_length = 0.5,
    building_area_fraction = 0.3,
    building_frontal_fraction = 0.11,
    building_height = 20.0,
    window_fraction = 0.2,
    street_canyon_aspect_ratio = 0.5,
    building_type = 2,
    pavement_type = 2,
    anisotropic_street_canyons = .T.,
    street_canyon_orientation = 0.,
    soil_temperature = 288.0,
/
```

As SLUrb requires the use of the land surface model (PALM-LSM) to model the natural fluxes, [`land_surface_parameters`](../../../../../Reference/LES_Model/Namelists/#land-surface-parameters) has to be defined in the namelist as well. For heterogeneous surfaces, the netCDF driver must be used.


### Input preparation
The shape of the urban form in SLUrb is defined using morphological parameters, which can be computed from e.g. urban plans or high-resolution maps of urban topography. The morphological parameters SLUrb uses are standard parameters used in urban climate, with references on how these can be computed being readily available (e.g. [Lipson et al., 2022](https://doi.org/10.3389/fenvs.2022.866398)). Where high-resolution surface data is not readily available, for example maps of Local Climate Zones with corresponding look-up table values can be utilized.

These parameters, however, are well-defined only for reference areas with a size of at least $100\times100~\mathrm{m}^2$. As the grid cell size in PALM simulations is typically smaller than this, some upsampling will be needed, with e.g. nearest neighbour sampling from a coarser data. An alternative to direct upsampling is to use a sliding (rolling) window to compute them, with a window of at least $100\times100~\mathrm{m}^2$ in size centered around the PALM grid cell. This method has the benefit of providing smoother spatial gradients for the morphological parameters compared to direct upsampling with the nearest neighbour method.

As SLUrb doesn't currently implement a vegetation model on its own, the non-urban fraction is modelled as vegetation or water by the PALM's land surface model (LSM). Therefore, [the static driver](../../../../../Reference/LES_Model/IO-Files/Drivers/static) is still needed to define vegetation and water surfaces. It is recommended to set up the land surface in the static driver so that it corresponds to the dominant non-urban surface type in the given urban cell (e.g. short grass or water). For urban cells and low vegetation, it is recommended to set the rougness length for momentum (value of [vegetation_pars(4,y,x)](../../../../../Reference/LES_Model/IO-Files/Drivers/static/#static--variable--vegetation_pars) in the static driver) to match the urban roughness length in order to avoid underestimation of the total surface rougness. This affects only the momentum flux from the vegetation and is likely to be revised in the future versions of SLUrb.

Note that SLUrb doesn't currently model heat fluxes for other anthropogenic sources other than heat diffusion through building walls and roofs. To include emissions from e.g. traffic, industry and HVAC, dynamic inputs `shf_traffic`, `shf_external` and `qsws_external` need to be prescribed by the user.

## Outputs
SLUrb provides a range of output quantities specific to the urban environment, with both instantaneous and temporally averaged outputs available. All possible outputs are listed and described as a part of [the output documentation](../../../../../Reference/LES_Model/Output_quantities#slurb-quantities). Activating SLUrb has some side effects to the PALM's core surface outputs: `shf*`, `qsws*` and surface radiation balance outputs will become aggregated fluxes. If the user wishes to output the urban and natural surfaces separately, SLUrb provides additional output quantities to access the non-aggregated surface fluxes (e.g. `slurb_shf_urban*` and `slurb_shf_lsm*` for urban and natural sensible heat flux).

# Available namelist parameters

{{ include_palm_namelist('slurb_parameters', as_table=True) }}
{{ include_palm_namelist('slurb_parameters') }}
