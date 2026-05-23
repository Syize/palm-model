---
title: Self Nesting
---
# Self Nesting

---

!!! warning
    This site is  Work in Progress.


## Purpose
The motivation of self nesting is to allow local resolution refinements within the total modelling domain. This makes it possible to model a larger domain and still have locally high resolution in principal areas of interest or other areas where the solution requires higher resolution than elsewhere. Having the highest required resolution all over the domain would lead to computationally very heavy setup or alternatively too small domain with boundary conditions too close to the principal areas of interest. Using the self nesting this problem can largely be avoided by concentrating the highest resolution to certain areas only and using somewhat lower resolution elsewhere.

## General Information
PALM allows nested configurations consisting of a number of model domains nested inside the largest domain called the **root** domain or root model. Domains nested into the root model are called **child** domains. Up to 63 child domains are allowed. PALM is run for all these domains actually as they were individual model runs but they are all run together as one parallel run in which also the inter-domain (inter-model) communication is handled by MPI. This means that there are two levels of parallelization: the usual intra-model parallelization (required by the domain decomposition) and the inter-model communication. Only one run is launched by the user, but the input data must be prepared for each model domain separately. Similarly output is done separately from each model (see further below).

Children can be recursively nested within each other, so a model domain can be **parent** and **child** at the same time (see child 1 in Fig. 1). Child domains can also be parallel to each other sharing the same parent domain (see child 1 and child 3 in Fig. 1). In other words, cascades of nested domains can be set up as well as parallel nested domains or parallel cascades. So, one parent may have more than one child and these children can again have
one or more children. The parent-child relations are defined using domain ids and parent ids (see further below).

![parent_child_nesting](../../Figures/self-nesting_parent_child_nesting.png){width=50%} <br>
**Figure 1:** 2d-sketch of parent/child domain layout. The root domain (which is always a parent domain) is displayed in light blue. Child 1 and child 3 are independent children embedded in the root domain, while child 2 is embedded into child 1. Therefore, child 1 acts as a parent as well as a child. The arrow symbol indicates that there must be a distance of minimum 4 parent grid spacings between the parent boundary and the child boundary.

By default, data exchange between parents and children is carried out at each PALM time step (but see chapter [Synchroneous / Asynchroneous Coupling](#synchroneous-asynchroneous-coupling)). Each model first calculates its time step depending on the time step requirements (CFL-criterion, etc.), and the minimum time step from all models is then used for all models.

Nesting data transfer can be made in two alternative modes:

- **One-way** means that children get solution data from their parents for setting
  boundary conditions on their nested boundaries but parents do not use
  any solution data from their children.

- **Two-way** means that also parents use solution data received from their
  children to modify the parent solution in the area occupied by a
  child domain using so-called post-insertion procedure. This is also called anterpolation in general.

Data of all prognostic variables is exchanged, except for the SGS-TKE. This exchange would have no real benefit and a coupling of SGS-TKE is everything else than straightforward since it strongly depends on the grid resolution. In case of two-way coupling and cascading arrangement of more than one nested domain, the order of data-transfer operations (from parent to
child and from child to parent) can be selected in alternative orders. For further details about this, as well as boundary-conditions, interpolation, and anterpolation see the [reference section](../../../../Reference/LES_Model/Modules/Nesting/self_nesting.md#nesting-modes-interpolation-and-anterpolation).


## Setup Requirements
- All child domains must be completely embedded within their parent domain. Overlapping of parallel child domains is allowed, but only for one-way nesting mode (see nesting parameter [nesting_mode](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--nesting_mode)). A child must always have only one parent domain (see Fig. 2). By default, lower boundaries of children match the lower boundary (surface) of the parent, but children can also be elevated with respect to their parent.

- It is important that the parent-grid lines always match child-grid lines meaning e.g. that the grid-spacing ratio (GSR) must be integer valued in each direction (see Fig. 3). This also means that the child-grid dimensions must be divisible by the corresponding GSR. Furthermore, this is required also on the subdomain level, i.e. even the subdomain dimensions must be divisible by the corresponding GSR. If e.g. the child grid spacing is 5 m, the parent grid spacing must be 10 m, 15 m, 20 m, etc.. It is recommended to avoid grid spacing ratios larger than *5*, because the zones at the child boundaries where turbulence adjusts to the finer grid spacing may become too large (and maybe even larger than the horizontal extent of the child).

- Vertical grid stretching is only allowed in the root domain above the top level of the highest located child domain.

- The 2d domain decomposition of child domains (as forced by the parallelization) must guarantee that the sub-domain size is always larger than the parent grid spacing in the respective direction.

- Anterpolation must not be done right next to the nested boundaries because this would create an unstable feedback loop usually leading to blow up of the solution. Therefore there are buffer zones right inside the nest boundaries where anterpolation is not made, as indicated by the arrow symbol in Fig. 1. The default width of these zones is two parent-grid cells, and the width can be set to different values via nesting parameter [anterpolation_buffer_width](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--anterpolation_buffer_width). Anterpolation can also be switched off within canopy layers in order to avoid some unphysical
secondary flow phenomena (see [anterpolation_starting_height](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--anterpolation_starting_height) and Hellsten et al. (2021).

![parent_child_nesting](../../Figures/self-nesting_parent_child_nesting_wrong.png){width=50%} <br>
**Figure 2:** 2d-sketch of non-allowed child domain position. Overlapping of children positions as marked by the red cross is only forbidden in case of two-way nesting mode. 

![parent_child_nesting](../../Figures/self-nesting_parent_child_grids.png){width=50%} <br>
**Figure 3:** Sketch of the alignment of parent and child grid lines. All parent grid lines must exactly match a child grid line.




## Basic Usage / Settings

A nested setup requires additional input files and generates additional output files. They all have a suffix that gives the respective child id (see further below). In the following, `<ri>` stands for the run-identifier of a respective run.

- Each model/domain has its own parameter file: `<ri>_p3d` (root domain),  `<ri>_p3d_N02` (first child), `<ri>_p3d_N03` (second child), etc..
- Additional input files (e.g. static and dynamic drivers) have to be provided for each child domain using the domain tags, e.g. `<ri>_static_N02`, `<ri>_static_N03`, etc..
- The data output is done for each domain in a separate file, using the domain tags, e.g. for the run-control output in file `<ri>_rc` (root domain), `<ri>_rc_N02` (first child), `<ri>_rc_N03` (second child), etc..

The nesting is switched on by adding the additional namelist [&nesting_parameters](../../../../../Reference/LES_Model/Namelists/#nesting-parameters) to the namelist file of the root domain only (`<ri>_p3d`)! The list below shows typical settings for a nested run with one child.
```Fortran
&nesting_parameters 
         domain_layouts = 'coarse', 1, -1, 64, 0.0,   0.0,   0.0,
                          'fine',   2,  1, 64, 320.0, 160.0, 0.0,
            
         nesting_mode         = 'two-way',  
         nesting_bounds       = '3d_nested’,
         synchronize_timestep = .FALSE.,

/
```
The most important parameter to define the number and position of children to be used is [domain_layouts](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--domain_layouts). For each model (root as well as children) 7 values have to be given.

1. The first one is an arbitrary string to give the model a name.
2. It is followed by an id that identifies the model. The root domain must have the id *1*, the first child the id *2*, and so on. Among others, the id is used as a suffix for the respective I/O files.
3. The third column gives the id of the parent, to which the child belongs. The root domain has no parent as indicated by the parent id *-1*. The parent of the first child is the root domain, which has the id *1*.
4. The fourth column defines the number of cores to be assigned to the respective model. The sum of cores over all models must match the total number of cores that is given via palmrun option `-X`.
5. The fifth column gives the x-coordinate of the leftmost grid point of the child in m relative to the leftmost grid point of the root domain. The root domain x-coordinate must always be *0.0*.
6. The sixth column gives the y-coordinate of the southernmost grid point of the child in m relative to the southernmost grid point of the root domain. The root domain y-coordinate must always be *0.0*.
7. The seventh column gives the vertical shift (in meters) of a child with respect to the root domain. For more details see section [Elevated children](#elevated-children).

Parameter [nesting_bounds](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--nesting_bounds) is used to set principal nesting setup types. *'3d_nested'* means a full 3d-nested simulation, where the child domains along all directions have a smaller size than their respective parents, as shown in Figs. 1-3. This is the default value. Other possible configurations are described in section [Nesting Setup Types*](#nesting-setup-types).

## One-Way / Two-Way Coupling
The use of one-way or two-way coupling via parameter [nesting_mode](../../../../../Reference/LES_Model/Namelists/#nesting_parameters-nesting_mode) has certain advantages and disadvanteges.

**One-way coupling advanteges:**

- Saves computational costs since anterpolation and child to parent data transfer are omitted.
- It is safer to use than two-way coupling since there is no risk of producing unphysical secondary-flow features. 

**One-way coupling disadvanteges:**

- No feedback from child to parent. If e.g. a scalar is released only in the child domain, it will not appear in the parent. Also wakes of wind turbines that are placed in a child would not be seen in the parent.
- The decoupling of turbulence may cause some flow discontinuities at the interface between parent and child.
- Buildings or steep orography very near the lateral boundaries of children may generate very large velocities (see section [Notes, Shortcomings, and Open Issues](#notes-shortcomings-and-open-issues))

**Two-way coupling advanteges:**

- Effects of highly resolved child domains are well represented in the parent.

**Two-way coupling disadvanteges:**

- Higher computational costs compared to one-way coupling.
- Simulations may become numerically unstable because of feedbacks between interpolation and anterpolation. This can often be avoided by adjusting the nesting parameter [anterpolation_buffer_width](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--anterpolation_buffer_width).
- Because the child solution often differs from the parent solution (which is to be expected because of the higher spatial resolution), artificial horizontal temperature gradients may appear at the boundaries between children and parents, causing secondary (thermal) circulations. They may be well visible in temporally averaged velocity fields, especially for horizontally homogeneous setups. This is an inherent feature of the two-way coupling and can't be avoided.

If the focus of a study is on the child-domain data only, one-way coupling is recommended in general. Results of two-way coupling setups should, if possible, always be compared with a reference run without nesting, which uses the spatial resolution of the child for the complete root domain, in order to identify possible numerical artifacts.

Regardless of using one-way or two-way coupling, the 3d child domains are always initialized with the 3d data received from its parent domain. Instead of initialization with the 3d data, a horizontally homogeneous initialization using the parent data can be chosen. For more details see nesting parameter [homogeneous_initialization_child](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--homogeneous_initialization_child).


## Synchroneous / Asynchroneous Coupling
Nesting parameter [synchronize_timestep](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--synchronize_timestep) enables or disables timestep synchronization between parent and child. If [synchronize_timestep](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--synchronize_timestep) = *.T.*, the time-integration for each model is done using the smallest required time step among all models (usually the innermost child domain). [synchronize_timestep](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--synchronize_timestep) = *.F.* allows models with larger time steps to use them, i.e. no time step synchronization between parent and child domains is performed. The model with the largest time step (usually the outermost parent domain) then determines the time of the next data exchange between parent and child domains. This helps to relax computational requirements, since models with coarser grid spacings may carry out less time steps.

Note that [synchronize_timestep](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--synchronize_timestep) = *.T.* is only allowed in conjunction with [nesting_mode](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--nesting_mode) = *'one-way'* and [particle_coupling](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--particle_coupling) = *.F.*.

Setting [synchronize_timestep](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--synchronize_timestep) = *.F.* usually requires to adjust the number of cores assigned to each model, in order to avoid load imbalances. See section [How to Optimize Usage of Computational Resources](#how-to-optimize-usage-of-computational-resources).


## Nesting Setup Types

The setup type can be chosen via parameter [nesting_bounds](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--nesting_bounds). Allowed values are:

- **'3d_nested'**<br> It is the most general nesting setup type. This means that child domains do not share any other boundaries with their parent except their bottom boundaries. However, even the child bottom boundary may be detached from the parent's bottom boundary (see section [Elevated Children](#elevated-children)). In 3-d nesting, there must be a clearance of at least four parent-grid spacings (two if the PW advection scheme is used) between parent- and child boundaries. However, it is recommended to use a much larger clearence especially on the upstream side of the child because adaptation of the solution to finer resolution of child may
take quite long distance, the longer the higher grid spacing ratio is (see Hellsten
et al., 2021).

- **'vertical_only'**<br> is another type of nesting setup. In this case only top (and bottom if the nest is elevated) are nested and the lateral boundaries are collocated with the parent lateral boundaries and cyclic conditions ([bc_lr](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_lr) = *'cyclic'* and [bc_ns](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_ns) = *'cyclic'*) are applied on these boundaries in both parent and child domains. In this case two-way coupling is recommended. Of course, grid spacings and number of grid points for root and child must be chosen to exactly match the horizontal domain size along x or y.<br><br> For example, it is a way to better resolve the complete surface layer. Recursive nesting is possible, i.e. the child can contain further children that either use 3d-nesting, 2d-nesting, or pure vertical nesting. See [Giersch and Raasch (2023)](https://doi.org/10.1007/s10546-023-00792-3) for an application of vertical nesting with PALM. It is highly recommended to not use the y-shift setting (see parameter [y_shift](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--y_shift)) for the cyclic conditions in vertically only nested runs, because it easily leads to different amount of shift in parent and child and thus conflicting solutions.

- **'cyclic_along_x'**<br> is a combination of the vertical only and 3-d nesting types of setup. In this case the left and right boundaries are set as in the vertical only case but the south and north boundaries are nested and these child boundaries must be located inside the
parent domain as in the 3-D nesting case (see Fig. 4). [bc_lr](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_lr) = *'cyclic'* must be set for both root and child domain. Also in this case, it is not recommended to use the y-shift setting for the cyclic conditions.

- **'cyclic_along_y'**<br> is similar to cyclic along x but with the cyclic boundaries being south and north (see Fig. 5). [bc_ns](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--bc_ns) = *'cyclic'* must be set for both root and child domain.

![2d_nesting_x](../../Figures/self-nesting_2d_nesting_x.png){width=50%} <br>
**Figure 4:** Sketch of a 2d nested setup with cyclic boundary conditions along x. Child and parent size along x exactly match.

![2d_nesting_y](../../Figures/self-nesting_2d_nesting_y.png){width=50%} <br>
**Figure 5:** Sketch of a 2d nested setup with cyclic boundary conditions along y. Child and parent size along y exactly match.


## Elevated Children
Many nested setups assume that the height of the bottom boundary in all models is the same (z = 0.0 m). In contrast to that, the bottom child boundary can be elevated with respect to the root model. The vertical shift is given in the seventh column of [domain_layouts](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--domain_layouts). A value of *0.0* means that the bottom boundary of the child matches the bottom boundary of the root domain. The root domain always requires to set the vertical shift to *0.0*. The bottom coordinate of a child (as given by the vertical shift) must match a vertical grid level of the root domain. If e.g. the vertical grid spacing of the root domain is [dz](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--dz) = *20.0*, only values of *20.0*, *40.0*, *60.0*, etc. are allowed for the vertical shift.

Elevated children may be used for various reasons:

- to better resolve specific layers, e.g. the entrainment zone at the top of the convective boundary layer ([nesting_bounds](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--nesting_bounds) = *'vertical_only'* shoud be set here).
- to save computational resources, if the child is put in an area with high orography. Here the vertical shift should be given the value of orography at the child grid point with lowest orography. With the usual vertical shift of *0.0*, a lot of child grid points may lie within the orography, which might be a waste of memory and computational resources.
- to better resolve the flow around the top of a tall building or in the rotor blade area of a wind turbine.


## RANS-RANS / RANS-LES Nesting
PALM can run either in LES or in RANS mode. Different turbulence closures (three for LES mode, two for RANS mode) are available (see [turbulence_closure](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--turbulence_closure)).

Nesting can be applied for both modes:

- RANS – RANS nesting ( 1-way or 2-way coupling )
- LES – LES nesting ( 1-way or 2-way coupling )
- RANS – LES nesting ( 1-way coupling only ) <br> Here the root (parent) model must be run in RANS mode, while only the child is running in LES mode (see Fig. 6). It requires to switch on a mechanism to initiate turbulence at the lateral child boundaries (see [Synthetic Turbulence Generator](../../Synthetic_turbulence)). Buffer zones at the boundaries of the LES child are usually required to allow for full development of turbulence. The width of the buffer zones have to be determined by test runs. There are no rules of thumb for that.

![RANS LES nesting](../../Figures/self-nesting_RANS_LES_nesting.png){width=50%} <br>
**Figure 6:** Sketch of a RANS-LES nesting setup. Inside the root model running in RANS mode there is a child model running in RANS mode, too. It contains a child running in LES mode, which in turn contains another child running in LES mode.

![adjustment zone](../../Figures/self-nesting_adjustment_zone.png){width=50%} <br>
**Figure 7:** Sketch of a RANS-LES nesting setup where the turbulence adjustment zones in the child are marked by a dashed area. The turbulence generator adds turbulence at all grid points of the child boundary (green line). If the wind is coming from south west (as indicated by the arrows) throughout the whole simulation, no adjustment zones would be required at the right and north boundaries. Adjustment zones mean that no child data from the adjustment zones should be used for analysis.

A good way to find out the required width of adjustment zones is to generated turbulence spectra at different distances from the child boundary and look for convergence of spectra. **Attention:** The flow adjustment zones significantly increase with increasing parent/child grid spacing ratio. For more details see [Hellsten et al. (2021)](https://doi.org/10.5194/gmd-14-3185-2021).


## Nesting with the Lagrangian Particle Model

The Lagrangian Particle Model allows to use nesting in general. Nesting for Lagrangian particles is activated via nesting parameter [particle_coupling](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--particle_coupling) = *.T.*. Lagrangian particles are then transferred between domains. A particle only exists once, either in the parent or in the child domain, i.e. if it enters a child domain, it is removed from the respective parent, and if it leaves the child, it is removed from that child and added to the parent domain again. The [&particle_parameters](../../../../../Reference/LES_Model/Namelists/#particle-parameters) namelist is only allowed to be given for the root domain. Particles that are released at the position of children are transferred to the respective child before the first time step.

The change of stochastic subgrid-scale particle speeds (see particle parameter [use_sgs_for_particles](../../../../../Reference/LES_Model/Namelists/#particle_parameters--use_sgs_for_particles)) that should appear at the partticle transsition from parent to child or vice versa is not considered in the code. This may have effects near the surface layer, where SGS velocity contributions are usually high. Results should be carefully analyzed if [use_sgs_for_particles](../../../../../Reference/LES_Model/Namelists/#particle_parameters--use_sgs_for_particles)) = *.T.* has been set.


## First Activation of Nesting in a Restart Run

In order to save computational resources, the nesting can be activated after the root model has advanced for a specific time interval. For example, the root model usually requires a spinup time to fully develop turbulence. During this time it makes no real sense to run the child model. **Note: This method allows to activate only one child per restart.**

First prepare the namelist files for the root and child model. The root model requires a file for the initial and the restart run (`_p3d`, `_p3dr`) while the child requires a file only for the restart run (`_p3dr_N02`).

Omit the [nesting_parameters](/Reference/LES_Model/Namelists/#nesting-parameters) namelist in the `_p3d` file of the root model only, and choose [end_time](/Reference/LES_Model/Namelists/#runtime_parameters--end_time) as the time when the nesting shall be activated.

The restart file for the root model (`_p3dr`) must contain the [nesting_parameters](/Reference/LES_Model/Namelists/#nesting-parameters) namelist. Start the initial run for the root domain via command 
```bash
palmrun .... -a "d3# restart"
```
Set [end_time](/Reference/LES_Model/Namelists/#runtime_parameters--end_time) in the `_p3dr` file of the root domain and the `_p3dr` file of the child to the desired value. In `_p3dr` set the nesting namelist parameter [init_child_id](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--init_child_id) = *2*. In `_p3dr_N02` you should set [initialization_actions](../../../../../Reference/LES_Model/Namelists/#initialitation_parameters--initialization_actions) = *'read_restart_data'* to allow for possible further restarts. If the child has been activated via [init_child_id](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--init_child_id), it will not read data from the restart file but will be initialized with data from the root model instead.

Start a restart run via command
```bash
palmrun .... -a "d3r"
```
The run will use the `_p3dr` file of the root domain, and file `_p3d3_N02` file of the child.

If option `-a "d3r restart"` has been given, further restart runs may be carried out via command

```bash
palmrun .... -a "d3r"
```
but the parameter [init_child_id](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--init_child_id) has to be removed from the `_p3dr` file before, because otherwise the child would be initialized again with data from the root model and not use the restart data.

In case a second child shall be added to the second restart, [init_child_id](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--init_child_id) = *3* has to be set and and additional information for defining the second child has to be added to [domain_layouts](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--domain_layouts) in file `_p3dr`. Also, a namelist file `_p3dr_N03` has to be provided. In the same way, children can be added for further restarts, but just one at a time. 

Keep in mind to always use activation string `restart` to store restart data for the next restart run, if necessary.

Another simple way of delayed activation of children is to run the root domain first (without children) as a precursor run for a specific length of time, and feed the data from the end of this run to the main run using [initializing_actions](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--initializing_actions) = *'cyclic_fill'* for the root domain, with the same size as the precursor-run domain. All children models can then be activated at the same time by adding the [&nesting_parameters](../../../../../Reference/LES_Model/Namelists/#nesting-parameters) namelist to the `_p3d` file.


## How to Optimize Usage of Computational Resources

Take care of a reasonable load balance between children and parents. The subdomain sizes of parents and children in terms of grid points should be comparable, i.e. they should cover about the same number of grid points. You can find information about the subdomain sizes that a setup is using in the run-control files (files with suffixes `_rc`, `_rc_N02`, etc. in the `MONITORING` folder). Otherwise the parent always waits for the child or vice versa after a time step. Also non-synchronized timesteps need to be considered. Here the parent model may use a smaller number of time steps than a child, and so less cores should be assigned to the parent. Any load imbalance can be checked via the CPU-measures files (files with suffixes `_cpu`, `_cpu_N02`, etc. in the `MONITORING` folder). CPU times indicated by lines starting with strings `nesting`, `pmc parent send`, `pmc parent recv`, `pmc child send`, or `pmc child recv` should not dominate the total CPU time required for the run.

A possible load imbalance may also have other reasons, e.g. because chemistry is only switched on in a child, or if excessive data output is done by a specific child only, but there may be many other reasons. The CPU-measures files usually helps to identify the reasons for load imbalance.

Another reason for load imbalance is that data transfer between parents and children is not treated well by the cluster network. In such cases you should contact the IT-support of your cluster.


## A Simple Complete Nesting Setup

For a simple complete nesting setup see the example files [cbl_particle_nested_p3d](https://gitlab.palm-model.org/palm/model/-/blob/master/tests/cases/cbl_particle_nested/INPUT/cbl_particle_nested_p3d) and [cbl_particle_nested_p3d_N02](https://gitlab.palm-model.org/palm/model/-/blob/master/tests/cases/cbl_particle_nested/INPUT/cbl_particle_nested_p3d_N02). It is a setup for a simple convective boundary layer with a weak mean wind and a constant surface sensible heat flux. The root domain uses a grid spacing of 50 m, and the child domain a grid spacing of 25 m. The Lagrangian particle model is activated. Note that the respective [&particle_parameters](../../../../../Reference/LES_Model/Namelists/#particle-parameters) namelist is only contain in the `_p3d` file. The [&nesting_parameters](../../../../../Reference/LES_Model/Namelists/#nesting-parameters) section in the `_p3d` shows that two cores each are used for the root and the child model, so the `palmrun` option `-X4` has to be used to run this setup.

These files may be used to develop and run more complex setups.



## References

For more detailed scientific and technical information about the self-nesting see the [reference section](../../../../Reference/LES_Model/Modules/Nesting/self_nesting.md) and [Hellsten et al. (2021)](https://doi.org/10.5194/gmd-14-3185-2021).

## Notes, Shortcomings and Open Issues

1.  In general, the two-way nesting mode is more challenging for users because it involves
a risk of unphysical features in the solution. One must be very careful when applying it. The one-way mode is more safe in this sense.
2. In case of one-way nesting, avoid to place buildings or steep orography at the lateral boundaries of children. At inflow boundaries this may cause very high velocities because of flow adjustment. For example, if the horizontal velocity at the inflow is large, a building immediately placed behind the inflow boundary may cause very high vertical velocities. Buildings or any changes in orography should keep a horizontal distance of at least three grid points from the lateral boundaries. If required, appropriately filter the static driver data.
3. A y-shift of data at the inflow (see parameter [y-shift](../../../../../Reference/LES_Model/Namelists/#initialization_parameters--y_shift)) does not work for pure vertical nesting or nesting with cyclic conditions along x in root and children models.
4. For using nesting in the ocean mode see [nested ocean runs](../../Ocean/#nested-ocean-runs).
5. So far, using a spin-up for child domains is not allowed (see [spinup_time](/Reference/LES_Model/Namelists/#initialization_parameters--spinup_time)).
6. Specific OpenMPI libraries, networks (`Omnipath`), and related system libraries (e.g. `libfabric`) have problems with the so-called one-sided communication, which is an MPI-3 feature that is used for the nesting communication. If you see unexpected and sometimes non-reproducable crashes of nested PALM runs (especially when using more than 2 children), try to switch to two-sided communication by setting the nesting parameter [use_one_sided_communication](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--use_one_sided_communication) = *.F.*. Note that this will decrease the performance of a run, but usually by not more than 10%. [use_one_sided_communication](../../../../../Reference/LES_Model/Namelists/#nesting_parameters--use_one_sided_communication) = *.F.* does not work for nesting with the Lagrangian particle model switched on.
7. The change of stochastic subgrid-scale particle speeds (see particle parameter [use_sgs_for_particles](../../../../../Reference/LES_Model/Namelists/#particle_parameters--use_sgs_for_particles)) that should appear at the particle transition from parent to child or vice versa is not considered in the code. This may have effects near the surface layer, where SGS velocity contributions are usually high. Results should be carefully analyzed.

