---
title: Palmrun Script
---
# The PALM Run Script

---

!!! warning
    This site is Work in Progress.

    ToDo:

    - [ ] Several links need to be fixed
    - [ ] Text need to be split into Guide / Reference part


The main script to execute PALM is called `palmrun`. This chapter describes the actions carried out by `palmrun` and gives a complete list of its available [options](http://127.0.0.1:8000/Reference/LES_Model/Scripts/options/). PALM can run in different modes:

- [**Interactive mode**](#interactive-mode)<br>
 PALM executes (almost) immediately within your terminal session after entering the    `palmrun` command.
 
- [**Batch mode**](#batch-mode)<br>
PALM job is submitted by `palmrun` to a queuing/batch system (e.g. PBS, SLURM, ...), where it is scheduled for execution.

The handling of PALM differs between interactive and batch mode, and it slightly varies, depending if a job is submitted to the:
 
- [**Local host**](#running-palm-in-batch-on-a-local-computer)<br>
The system that you are currently logged in.

- [**Remote host**](#running-palm-in-batch-on-a-remote-computer)<br>
Any remote computer with a batch system, to which you have `ssh` access, but are not logged in. The remote host becomes the local host after login to the remote host via `ssh`.

## Interactive Mode

The following instructions assume, that PALM has been correctly installed. Switch to your working directory (the one that contains the configuration file `.palm.config.default`) and enter 
```
palmrun -r example_cbl -c default -a "d3#" -X 4
```

The progress of execution can be followed on the terminal, where informative messages will be output (the execution can be stopped at any time by typing `Ctrl+C`). Some general settings will be listed on the terminal and are prompted for confirmation:
```
*** palmrun
    will be executed.     Please wait ...

    Reading the configuration file... 
    Reading the I/O files... 

  *** INFORMATIVE: additional source code directory
      "/home/<local_username>/palm/current_version/JOBS/example_cbl/USER_CODE" 
      does not exist or is not a directory.
      No source code will be used from this directory!

#------------------------------------------------------------------------# 
| palmrun                               Thu Jun 26 08:49:32 AM CEST 2025 | 
| Version: PALM release ...                                    | 
|                                                                        | 
| called on:               <local host name>                             | 
| config. identifier:      imuk (execute on IP: 111.11.111.111)          | 
| running in:              interactive run mode                          | 
| number of cores:         4                                             | 
| tasks per node:          4 (number of nodes: 1)                        | 
|                                                                        | 
| cpp directives:          -cpp -D__parallel ...                         | 
| compiler options:        -O3 -g -fbacktrace -ffpe-trap=invalid,zero ...| 
| linker options:          -O3 -g -fbacktrace -L ... -l...               | 
|                                                                        | 
| run identifier:          example_cbl                                   | 
| activation string list:  d3#                                           | 
#------------------------------------------------------------------------#

 >>> everything o.k. (y/n) ?
```
Listed settings are determined by the given `palmrun` options and settings in the [configuration file](../config_file.md) (here `.palm.config.default`).
**Entering** `n` aborts `palmrun`, **entering** `y` starts execution of PALM, and some more informative messages will appear on the terminal. 
```
***  PALMRUN will now continue to execute on this machine

  *** creating executable and other sources for the local host

make: Nothing to be done for 'all'.
  *** executable and other sources created

  *** changed to temporary directory: /localdata/......./example_cbl.23751

  *** providing INPUT-files:
  ----------------------------------------------------------------------------
  >>> INPUT: /home/....../palm/current_version/JOBS/example_cbl/INPUT/example_cbl_p3d  to  PARIN
  *** INFORMATIVE: some optional INPUT-files are not present
  ----------------------------------------------------------------------------
  *** all INPUT-files provided 


  *** execution of INPUT-commands:
  ----------------------------------------------------------------------------
  >>> ulimit -s unlimited
  ----------------------------------------------------------------------------


  *** execution starts in directory
      "/localdata/....../example_cbl.23751"
  ----------------------------------------------------------------------------

  *** running on: hostname hostname hostname hostname
  *** execute command:
      "mpiexec -machinefile hostfile -n 4 ./palm" 

   08:55:30   -finished-   reading environment parameters from ENVPAR
   08:55:30   -start----   reading NAMELIST parameters from PARIN
   08:55:30   -finished-   reading NAMELIST parameters from PARIN
   08:55:30   -start----   creating virtual PE grids + MPI derived data types
   08:55:30   -finished-   creating virtual PE grids + MPI derived data types
   08:55:30   -start----   Setup topography
   08:55:31   -finished-   Setup topography
   08:55:31   -start----   checking parameters
   08:55:31   -finished-   checking parameters
   08:55:31   -start----   model initialization
   08:55:31   -start----   initializing surface layer
   08:55:31   -finished-   initializing surface layer
   08:55:32   -finished-   model initialization
   08:55:32   -start----   time-stepping

      [XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX]   0.0 left
   09:01:29   -finished-   time-stepping
   09:01:29   -start----   calculating cpu statistics
   09:01:29   -finished-   calculating cpu statistics

  ----------------------------------------------------------------------------
  *** execution finished 
```
In case that `palmrun` has proceeded to this point (`-finished-  time-stepping` and `*** execution finished`) without giving warning- or error-messages, the PALM simulation has finished successfully. The displayed progress bar (`xxxxx`) roughly allows to estimate how long the run still needs to finish.

Subsequent messages give information about post processing and copying of output data:
```
*** post-processing: now executing "mpiexec -machinefile hostfile -n 1 combine_plot_fields.x" ... 
 
 *** combine_plot_fields ***
     uncoupled run
 
 
     no XY-section data available
 
     no XZ-section data available
 
     no YZ-section data available
  
     no 3D-data file available

     
  *** execution of OUTPUT-commands:
  ----------------------------------------------------------------------------
  >>> [[ -f LIST_PROFIL_1D ]] && cat LIST_PROFIL_1D >> LIST_PROFILE
  >>> [[ -f LIST_PROFIL ]] && cat LIST_PROFIL >> LIST_PROFILE
  ----------------------------------------------------------------------------


  *** saving OUTPUT-files:       local time: ..:..:..
  ----------------------------------------------------------------------------
  >>> OUTPUT: RUN_CONTROL  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/MONITORING/example_cbl_rc.000

  >>> OUTPUT: HEADER  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/MONITORING/example_cbl_header.000

  >>> OUTPUT: CPU_MEASURES  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/MONITORING/example_cbl_cpu.000

  >>> OUTPUT: DATA_1D_PR_NETCDF  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/OUTPUT/example_cbl_pr.000.nc

  >>> OUTPUT: DATA_1D_TS_NETCDF  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/OUTPUT/example_cbl_ts.000.nc

  >>> OUTPUT: DATA_2D_XY_NETCDF  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/OUTPUT/example_cbl_xy.000.nc

  >>> OUTPUT: DATA_2D_XZ_NETCDF  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/OUTPUT/example_cbl_xz.000.nc

  >>> OUTPUT: DATA_2D_XZ_AV_NETCDF  to
              /home/<local_username>/palm/current_version/JOBS/example_cbl/OUTPUT/example_cbl_xz_av.000.nc

  ----------------------------------------------------------------------------
  *** all OUTPUT-files saved       local time: ..:..:..

 --> palmrun finished
```
Output files can be found at their respective positions as listed in the terminal output. Most of PALM's output files are written in NetCDF format and are copied to subdirectory `OUTPUT`. Some general information files are written in ASCII format and are copied to folder `MONITORING`. All available output files of PALM are listed [here](../IO-Files/).


## Batch Mode

Larger simulation setups usually cannot be run interactively, since larger amounts of required resources (memory as well as cpu-time) are only provided through batch environments. `palmrun` supports two different ways to run PALM in batch mode. In both cases it automatically generates a batch job, i.e. a file containing directives for the batch system plus commands to execute PALM, which is then either submitted to the local computer or to a remote machine. Running PALM in batch mode usually requires to manually modify and extend the [configuration file](../config_file.md), and that a batch system (e.g. slurm, pbs, ...) is installed on the respective machine. Because of the large variety of required batch directives, the installer is not able to automatically configure the configuration file. Please carefully check the batch system documentation or contact the system-support about the required batch-directives.

### Running PALM in Batch on a Local Computer

The local computer is the one where commands are enter in a terminal and executed. This might a local PC/workstation, or a login-node of a cluster-system / computer center which is accessed via ssh. Regardless of the computer, it is assumed that PALM has been successfully installed on that machine, either using the installer or via manual installation.

For running PALM in batch mode, additional options are required for the `palmrun` command to specify the system resources requested by the job. Additional modifications of the configuration file may be required, too (see further below). A minimum set of additional `palmrun` options is
```
palmrun  ....-b -c <configuration identifier>  -m <memory> -t <cputime>
             -X <total number of cores>  -T <MPI tasks per node>  -q <queue>
```
**Note**: The first option `-b` causes `palmrun` to create a batch job running on the local computer!

Before entering the above command, additional information needs to be added to the configuration file. **Best practice** would be to create a new file, e.g. by copying the default file `.palm.config.default` to e.g. `.palm.config.batch`, and then editing the new file. On a system that allows both batch and interactive mode in the same software environment, the same configuration file may be used for `palmrun` in either of the modes. For any newly created configuration file, the PALM source code has to be compiled again. Assuming the above name of the configuration file, compilation is done via
```
palmbuild -c batch
```
The compiled code will be put into folder `MAKE_DEPOSITORY_batch`.

More details and a complete description of the configuration file can be found [here](../config_file.md).

Based on the `palmrun` arguments, environment variables (for a description of available variables see [here](../config_file.md#variables-used-in-batch-job-directives)) will be set by `palmrun`as described below. The following list of automatically set environment variables assumes a `palmrun` call 
```
palmrun .... -t 5400  -X 48  -T 12  -q medium
```
- `{{ '{{run_id}}' }}` = example_cbl.#####  
  where ##### is a five digit random number which is unique for each job. The `run_id` is used for different purposes, e.g. it defines the name under which the job can be found in the queuing system.

- `{{ '{{cpu_hours}}' }}` = 1, `{{ '{{cpu_minutes}}' }}` = 30 and `{{ '{{cpu_seconds}}' }}` = 0
calculated from option `-t`.

- `{{ '{{mpi_tasks}}' }}` = 48
  as given by option `-X`.

- `{{ '{{tasks_per_node}}' }}` = 12
  as given by option  `-T`.

- `{{ '{{nodes}}' }}` = 4
  calculated from `-X`/`-T` . If `-X` is not a multiple of `-T`, `nodes` is incremented by one, e.g. `-X 49 -T 12` gives `nodes` = 5.

- `{{ '{{queue}}' }}` = medium
  as given by option `-q`.
  
After confirming the `palmrun` settings with `y`, following information (in extracts) will be output to the terminal:


```
 >>> everything o.k. (y/n) ?  y

  ***  batch-job will be created and submitted

   *** creating executable and other sources for the local host
make: Nothing to be done for 'all'.
   *** executable and other sources created

  *** submit the job
 <<<submit message from batch system>>>

  --> palmrun finished
```

Before the batch job is finally submitted, `palmrun` generates a folder named `SOURCES_FOR_RUN_<run_identifier>` which is located in the `fast_io_catalog` and which contains various files required for the run (e.g. the PALM executable, PALM's source code and object files, copies of the configuration files, etc.). Messages `*** executable and other sources created` and `*** input files have been copied` tell you that this folder has been created. `make: Nothing to be done for 'all'.` means that no user interface needs to be compiled. After the job submission, the batch system usually prompts a message (`<<<submit message from batch system>>>`) which contains the batch system id under which the job is registered in the queueing system (can be used e.g. to cancel that job). The job is now queued and will be finished depending on how many other jobs are waiting in the queue. The main task of the job is to execute the `palmrun` command again, that has been entered in the terminal, but now on the compute nodes of the system. A job protocol file with name `<configuration identifier>_<run identifier>` as given with `palmrun` options `-c` and `-r` (here it will be `batch_example_cbl`) will be put in the folder that has been set via variable `local_jobcatalog` in the configuration file (`.palm.config.batch`). Check contents of this file carefully. Beside some additional information, it mainly contains the output of the `palmrun` command as like during interactive execution, e.g. information is given to where the output files have been copied, etc..

Typically, batch systems allow to run jobs only for a limited time, e.g. 12 hours. See chapter [job chains and restart jobs](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/runs) on how `palmrun` is used to create so-called job chains in order to carry out simulations which exceed the time limit for single jobs.

### Running PALM in Batch on a Remote Computer

The `palmrun`command can be used on the local computer (e.g. local PC or workstation) to submit a batch job to a remote computer. `palmrun` copies required input files from the local computer to the remote machine and transfers output files back to the local machine, depending on the settings in the `.palm.iofiles` file. The job protocol file will also be automatically copied back to the local computer.

To use this `palmrun` feature, additional settings in the configuration file are required. Furthermore, the PALM-code needs to be pre-compile for the remote machine using `palmbuild`. **The installer can not be used to install PALM on that machine.** Most of the settings must be done manually.

Furthermore, passwordless ssh/scp access is required from the local computer to the remote computer, as well as from the remote to the local computer. In remote mode, `palmrun` and `palmbuild` are heavily using `ssh` and `scp` commands. If there is no established passwordless access, passwords would need to be entered several times before the batch job is finally submitted. Moreover, the job protocol file and any output files cannot be transferred back to the local computer because there is no connection to the job which could be used to provide passwords for these transfers.

The configuration file requires specific settings for remote batch jobs. For this it would be convenient to create a new configuration file based on an already existing one, e.g. via
```
cp  .palm.config.batch  .palm.config.batch_remote
```
where `batch_remote` can be any string to identify the remote host. Edit this file as described [here](../config_file.md#additional-directives-for-batch-jobs-on-remote-hosts).

After setting up the configuration file and before calling `palmrun`, `palmbuild` needs to be called to compile the source code and generate the PALM executable for the remote host:
```
palmbuild -c batch_remote
```
Keep in mind that the configuration file `.palm.config.batch_remote` requires correct settings valid for the remote computer (compiler name, compiler options, include and library paths, etc.). If `palmbuild` has not been called, `palmrun` will automatically do it.

If `palmbuild` succeeded, the `palmrun` command can be entered like
```
palmrun -r example_cbl -c batch_remote ......
```
After confirming the `palmrun` settings with `y`, similar information as for local batch jobs will be output to the terminal. `palmrun` finally terminates with messsage `--> palmrun finished`. The batch job is now queued on the remote system. After the job has been finished, the job protocol will be transferred back to the local computer and put into the folder defined by `local_jobcatalog`. If this file does not appear, because e.g. the transfer failed, the protocol file can be found on the remote host in the folder defined by `remote_jobcatalog`. Like in case of batch jobs running on local computers, check the contents of this file carefully. Beside some additional information, it mainly contains the output of the `palmrun` command as for interactive execution, and especially information about where to find the output files on the local computer.

**Note**: Large PALM-setups (those using large number of grid points) can produce extremely large output files which would require long time for transferring them to your local system and which might have sizes that exceed the capacity of your local discs. See the description of the [I/O file connection configuration](../../../IO-Files/IO_Configuration_File) on how to control copying of INPUT/OUTPUT files.

