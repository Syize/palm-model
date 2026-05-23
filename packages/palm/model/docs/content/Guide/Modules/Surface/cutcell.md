---

!!! warning
    This site is Work in Progress.

    ToDo:

    - [ ] give more detailed explanations about usage and features


# Cut-Cell Topography

## Cut-Cell Topography related Parameters and Inputs

To enable a cut-cell representation in PALM, the flag parameter [cut_cell_topography](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--cut_cell_topography) needs to be set.
In order to have information on the location and orientation of the cut-cell surface, further cut-cell specific input is required, which is described in the [static input file documentation](../../../../../Reference/LES_Model/IO-Files/Drivers/static) and in the [technical description of the cut-cell method](../../../../../Reference/LES_Model/Modules/Surface/cutcell).
The [PALM-GeM tool](https://gitlab.palm-model.org/static_driver/palm_gem) has the capability to generate the cut-cell surface inputs.
To visualize data on 3D topology, the postprocessing tool `slanted_surface_output.py` to be found in the repository folder `packages/palm/model/share/utils/` can be used.

**A mixture of Cartesian step-like surfaces and cut-cell surfaces is currently not possible, e.g. vertically oriented building walls can not be represented via cut-cell surfaces.**
