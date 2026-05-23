---
title: Overview
---
## User-interface Overview

Programming user-defined code extensions usually requires exact knowledge of the internal structure of PALM. The latest PALM publications and the technical documentation are usually not sufficient and must be supplemented by a rigorous study of the model's source code. Programming experiences with Fortran and MPI are absolutely required. Please note that support for developing user-interfaces via the trouble ticket system can not be given.

Changes of the default model code by the user should be avoided whenever possible and are reserved to the developer group of PALM. The corrections, revisions and extensions of the model accomplished by this group are published in the technical/numerical documentation and respective new releases are accessible to users.

However, users frequently may feel the need to extend the model code for his/her own simulations. For this purpose, a set of interfaces is available, which can be used to add user-defined code to the model. This chapter describes the programming of such user-defined code extensions.

By default, user-defined code embedded via subroutine calls at various places in the model code. These subroutines have pre-defined names, which must not be changed. Their basic versions are a contained (called) in the default code and can be found in the source code files `user_***.f90`. The basic versions accomplish nearly no executable code. They are pure templates, which can be extended by the user as required. Any executable code already contained in these basic versions should not be changed. An example of such a basic version (file `user_init.f90`) is given below:

```Fortran
SUBROUTINE user_init

!------------------------------------------------------------------------------
!
!
! Description:
! -----------
! Execution of user-defined initializing actions
!------------------------------------------------------------------------------
!

 USE control_parameters
 USE user

 IMPLICIT NONE

!
!-- Here the user defined initializing actions follow:


END SUBROUTINE user_init
```

Access of variables defined within PALM's default code is realized via individual `FORTRAN` - modules defined inside `PALM`. The appropriate modules (many of them can be found in source code file `modules.f90`) must be declared in the user-defined routines by means of `USE` statements. In the example above, this is realized for PALM module `control_parameters`. This yields access to most of the existing parameters for steering the model.

In addition, the module `user` appears in the example above. This is a user-defined module (it is part of file `user_module.f90`) and can be used for communication between the user-defined routines. Additional variables can be declared in this module as desired. It is not used outside of the user code, and should never be used in or added to PALM's default code.

A very typical request of users is the calculation and output of quantities which are not part of PALM's standard output. Several routines in the basic user interface are already designed and prepared for calculating and output of such quantities.

As already explained above, the contents of the files `user_***.f90` can be used as a basis for extensions. However, these files should not be manipulated directly. Instead, copies of these files should be created and modified instead.


## Compiling and Linking

To use user-interfaces within a PALM-run requires the following steps:

1. **Copy** the appropriate default (empty) user-interface files `user_***.f90` to the folder defined in the configuration file via variable `user_source_path`. The default path for this folder is `$HOME/palm/current_version/JOBS/$run_identifier/USER_CODE`, so in case the run identifier is `example`:  
        
        cd ~/palm/current_version
        mkdir -p JOBS/example/USER_CODE
        cp ..../packages/palm/model/src/user_***.f90  JOBS/example/USER_CODE  

    **Attention**: Don't change names of the user-interface files! 

2. **Modify** the interface routines.

3. **Add** the namelist `&user_parameters` to your namelist file (suffix `_p3d`). This is required to activate the user-interface during the run. If the namelist is not given, the user-interface will be compiled but not used. The namelist can be empty (`&user_parameters /`), if no specific user-parameters have been defined.

4. **Start** a PALM run by entering  

        palmrun -r example ...

The files `user_***.f90` will be **automatically compiled** as part of this run, before PALM is executed, and will **replace** PALM‘s respective default user-interface routines in `..../packages/palm/model/src/`. The compiled binaries are put into folder `SOURCES_FOR_RUN_$configuration_identifier_$run_identifier`, where they replace the default binaries that have been created via `palmbuild` before. The `SOURCES_FOR_RUN` folder is created in the directory defined by variable `fast_io_catalog` in the configuration file. It is newly generated for each manual call of `palmrun`. Since the run identifier is part of the user-interface's SOURCE path, you may defined and use different interfaces for different runs.

It is also possible to add **additional routines** that are not part of the default user-interface (see the [list of available interface routines](../../../../Reference/LES_Model/Modules/User_Interface.md)) by

- **appending** the code of such routines to the default user-interface file `user_additional_routines.f90` or

- **creating** new files with new names in the `USER_CODE` folder.

These files must have names different from the default PALM source code files. In case of new files, the default `Makefile` must be copied too, e.g.:

```
cp ..../packages/palm/model/src/Makefile  JOBS/example/USER_CODE
```

Edit the `Makefile` and add new files and possible dependencies, following the standard "make"-rules.

### How to avoid re-compilation

It is important to know, that in general the modified user-interface files cannot be pre-compiled with `palmbuild`. This would not make sense because the user-interface may differ among different PALM runs. That would required to use different `MAKE_DEPOSITORIES` depending on the run-identifier (option `-r`). Therefore, for each manual call `palmrun` newly compiles the user-interface. Strictly speaking, `palmrun` internally calls `palmbuild` with special options, and `palmbuild` then compiles the user-interface and generates the `SOURCES_FOR_RUN` folder (see above).

Compiling user-interfaces may require long time. Re-compilation with each call of `palmrun`  can be avoided by adding `palmrun` option `-V`. Then, `palmrun`/`palmbuild` tries to re-use the compiled interface from the last call of `palmrun` with same configuration identifier and run identifier (given by options `-c` and `-r`), by re-using the respective `SOURCES_FOR_RUN` folder that has been generated by the previous run. Automatically generated calls of `palmrun` (restart runs in job chains) use the same method, so they do not need to compile the user-interface again. Before using option `-V` be aware that the user-interface code hasn't changed.

## User-defined parameters

As for the model in general, also the user-defined code may require steering by parameters. For each run, the model should be able to read in current values of such parameters. The declaration of user-defined parameters is done in the user-defined module [`user`](../../../../Reference/LES_Model/Modules/User_Interface#user) (file `user_module.f90`). This module must be declared in all relevant user-defined routines via a `USE user` statement, in order to make the parameters available.

Values can be assigned to user-defined parameters via namelist group `&user_parameters` in the `_p3d` file. Parameters must be declared within the respective namelist statement in the user-defined subroutine [`user_parin`](../../../../Reference/LES_Model/Modules/User_Interface/#user_parin).

The following example illustrates the procedure. The example assumes that a `LOGICAL` variable named `abcd` has been declared for steering the user-defined code. This declaration must be done in file `user_module.f90`:

```Fortran
LOGICAL ::  abcd = .FALSE.
```

Pay attention that in this example a default value (*.FALSE.*) is assigned to the variable and that it will keep this value if nothing else has been assigned in the `_p3d` file. In routine [user_parin](../../../../Reference/LES_Model/Modules/User_Interface/#user_parin) the namelist must be extended by the name of the new variable: 

```Fortran
NAMELIST /user_parameters/  abcd, data_output_masks_user, data_output_pr_user, data_output_user, region
```

The listed parameters following `abcd` are defined in [user_parin](../../../../Reference/LES_Model/Modules/User_Interface/#user_parin) by default. A complete list of pre-defined parameters can be found [here](../../../../Reference/LES_Model/Namelists/#user-parameters). For `abcd`, a value can be assigned in the `_p3d` file:: 

```Fortran
 &user_parameters abcd = .T., ... ,/
```

User defined parameters in the `_p3d` file are considered as runtime parameters, i.e. they must be specified again for each restart run.

A control output of user-parameter values is recommended, to later check values used during the respective model run in the `_rc` file.. For this purpose, the user-defined subroutine [user_header](../../../../Reference/LES_Model/Modules/User_Interface/#user_header) is available, which writes into the `_header` and `_rc` files.

The user-module is automatically activated when the `&user_parameters` namelist appears in the `_p3d` file.

## User-defined output

See [here](output).

## Further information

For a complete list of available user-interface routines see the [reference section](../../../../Reference/LES_Model/Modules/User_Interface).


