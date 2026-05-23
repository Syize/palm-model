---
title: Configuration File
---
# The PALM configuration file

---

!!! warning
    This site is  Work in Progress.

Running PALM with the `palmrun` script or compiling PALM with `palmbuild` requires a **configuration file** in the working directory from where the scripts are called. The configuration file contains information about 

* the compiler and compiler options to be used 
* paths to required libraries (e.g. MPI, NetCDF, or fftw) 
* pre-processor directives to switch on/off special features of PALM 
* paths to be used for storing PALM's input and output files 
* paths where user-interface routines and job protocol files can be found 
* the execute command to be used for starting the PALM executable 
* settings of other UNIX environment variable 
* batch job directives 
* unix commands to be executed before / after the PALM code is started, or that should be carried out in case of errors that appeared during the run 

If PALM has been successfully installed by the automatic installer, the installer creates a configuration file named `.palm.config.default` in the working directory of the user, by default `$HOME/palm/current_version`. 
This file should contain all relevant information to run PALM on the respective computer in [interactive mode](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/palmrun#interactive). In case of any problems, please check the [installation guide](../../../../Get_Started/#installation) or [FAQ](https://palm.muk.uni-hannover.de/trac/wiki/Help/FAQ) page. If these pages do not help, contact us via the [ticket system](https://palm.muk.uni-hannover.de/trac/newticket).

For [batch jobs](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/palmrun#batch) and in [remote mode](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/palmrun#batch_remote), the configuration file has to be modified. 

## How to generate a configuration file manually
**This section needs to be moved to the Guide. Furthermore, the `.palm.config.default` in the repository requires to be updated!**

 As a first step, copy the default template of the configuration file (which is part of the PALM download) to your working directory: 
```
  cd $HOME/palm/current_version
  cp <palm release>/packages/palm/model/share/config/.palm.config.default .
```

## Configuration file format
The configuration file is an ASCII file which may be opened with any editor. Lines need to start with either
`#`, `%`, `IC:`, `OC:`, `EC:`, `BD:`, or `BDT:`. Empty lines are allowed.

* `#` starts comment lines

    `# this is a comment line`

* `%` in the first column defines unix environment variables used in `palmrun` and `palmbuild`:

    `%var value of var`

    A variable named `var` with value `value of var` is created. There must be at least one blank between the variable name and its value. The value may contain an arbitrary number of blanks. The value can contain already defined variables: 

    `%default_folder  /work/abcd`

    `%subfolder1      $default_folder/efgh`

    which means that `subfolder1` has the value `/work/abcd/efgh`. You can also use any variables that are already defined within the `palmrun`/`palmbuild` scripts. The most important one is `run_identifier`, given by `palmrun` option `-r`, and which defines the so-called *run identifier*. This variable is mentioned here because it is used in the default configuration files for naming I/O files and to better sort I/O files from different PALM runs in an organized folder structure. Beside the value replacement using the UNIX shell syntax (i.e. `$abcd` gives the value of variable `abcd`), another way to insert values of environment variables is to write them in double curly brackets, i.e. `{ {abcd} }` will be replaced by the value of `abce`. This way of variable replacement is required for some special variables and in the job directives. 


* `IC:` defines any unix command that is executed by `palmrun` just before the PALM model is started. A typical input command might be

    `IC:ulimit -s unlimited`

    to increase the system's stack size (if this is not unlimited by default), which is required for PALM runs with larger memory demands.

* `OC:` defines unix commands that are executed by `palmrun` just after the PALM model has stopped. For example, you may inform yourself about termination of the program by sending an email:

    `OC:echo` "PALM simulation $run_identifier has finished" | mailx  username@email-address

* `EC:` defines unix commands that shall be executed in case that the PALM model or the `palmrun` script terminated because of any kind of error. You can restrict execution of error commands to specific kinds of error, e.g. errors that appear during PALM execution:

    `EC:[[ \$locat = execution ]]  &&  error-command`
    
    See the `palmrun` script (.../packages/palm/model/bin/palmrun) for other specific locations that are used in this script. 

* `BD:` defines directives that are required for batch jobs, i.e. if PALM shall be run in batch mode. Explanations for batch directives are given further below. 

* `BDT:` defines directives for an additional batch job that is required in case of running PALM in batch mode on a remote host. This additional job transfers the job's log file of the main PALM job back to the local host. 


## List of environment variables and their meaning

The default template contains settings for the minimum set of variables that always need to be defined, plus suggestions for other variables that you may need to uncomment, in case that you like/need to use them. The following table lists all relevant variables and their meaning. No default values are assigned for these variables in scripts `palmrun` and `palmbuild`. Values given in the template file are for a Linux system with Intel Fortran compiler, mpich-, fftw-, and NetCDF4-library that is locally used by the PALM group. **last sentence requires update**


| Variable name     | Meaning    |
|-------------------|------------|
| base_data         | Directory where the PALM I/O-files are stored by default. This variable is used in the [file-connection file](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/palm_iofiles).<br> **Attention:** Since this variable is also used to determine file locations on remote hosts, it must use the `~` instead of `$HOME`. |
| base_directory    | Working directory from where `palmrun` or `palmbuild` is called. **Attention:** The configuration file(s) `.palm.config....` must reside here.|
| compiler_name <a name="compiler_name"></a> | Name of the Fortran compiler to be used to create the PALM executable. Typically, this is the name of a wrapper script like `mpif90` or e.g. `mpiifx`, which automatically invokes the required MPI library and MPI include file. If you don't have a wrapper script, you may need to explicitly give compiler options (see [compiler_options](#compiler_options)) to provide paths to the library / include file. If you like to run PALM without MPI (serial mode, or OpenMP parallelization), you should not use a wrapper script and give the original compiler name instead. |
| compiler_name_ser <a name="compiler_name_ser"></a> | Fortran compiler name to create non-MPI executables. This name is required, because `palmbuild` generates several helper programs for pre-/post-processing, which run in serial mode on just one core. Here you give the original compiler name, like `gfortran`, `ifx`, etc..  |
| compiler_options <a name="compiler_options"></a>  | Options to be used by the compiler that has been specified by [compiler_name](#compiler_name) and [compiler_name_ser](#compiler_name_ser) in order to compile the PALM and utilities source code. See the list of [recommended compiler options](../../../../Guide/LES_Model/Recommended_compiler_options) for specific compilers. Library paths do not have to be given here (although you can do that), but paths to INCLUDE files may need to be specified. |
| cpp_options       | Preprocessor directives to be used for compiling the PALM code. They allow for conditional compilation using the `-D` compiler option. Compiling PALM with MPI support requires options `-D__parallel -DMPI_REAL=MPI_DOUBLE_PRECISION -DMPI_2REAL=MPI_2DOUBLE_PRECISION`. Many compilers require to set an additional option to run the Fortran preprocessor on source files before compilation (e.g. `-fpp` for the Intel compiler). This option has to be given here, too. Alternatively, you can provide it as part of the [compiler_options](#compiler_options). See [cpp options](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/cpp_options) for a complete list of preprocessor define strings that are used in the PALM code.
| defaultqueue   <a name="defaultqueue"></a>   |  Batch job queue to be used if no queue is explicitly given with `palmrun` option `-q`.  |
| execute_command  <a name="execute_command"></a> | MPI command to start the PALM executable. Please see your local MPI documentation about which command needs to be used on your system. The name of the PALM executable, usually the last argument of the execute command, must be `palm`. Typically, the command requires to give several further options like the number of MPI threads to be started, or the number of compute nodes to be used. Values of these options may change from run to run. Don't give specific values here and use variables (written in double curly brackets) instead, which will be automatically replaced by `palmrun` with values that you have specified via respective `palmrun` options. As an example `mpirun -n {{ '{{mpi_tasks}}' }} -N {{ '{{tasks_per_node}}' }} palm` will be interpreted as `mpirun -n 240 -N 24 palm` if called via `palmrun ... -X240 -T24 ....` See the batch job section below about further variables that are recognized by `palmrun`. |
| execute_command_for_combine | Command to start the post processing tool `combine_plot_fields`. By default, the execute command given by `execute_command` will be used, with string *"palm"* replaced by string *"combine_plot_fields.x"*. This might not work, especially if `execute_command` contains options for number of cores or number of cores per node to be used. Since `combine_plot_fields` is not parallelized, it must be executed on one core only. In such cases, you need to add explicit setting of `execute_command_for_combine` to your configuration file. For a SLURM batch system the additional line may read `%execute_command_for_combine srun --propagate=STACK -n 1 --ntasks-per-node=1 combine_plot_fields.x`. |
| fast_io_catalog   | Path to a file system with fast discs (if available). This folder is used so store the temporary catalog generated by `palmrun` during each run. It should also be used to store large I/O files (e.g. restart data or 3D-output) in order to reduce I/O time. This variable is used in the default `.palm.iofiles` for the restart data files. The folder must be accessible from all compute nodes, i.e. it must reside in a global file system.<br> **WARNING:** `/tmp` will only work on single node systems! In case of batch jobs on remote hosts, the variable refers to a folder on the remote host. The variable has no default value and must be set by the user.  |
| hostfile          | Name of the hostfile that is used by MPI to determine the nodes on which the MPI processes are started. `palmrun` automatically generates the hostfile if you set `%hostfile auto`. All MPI processes will then be started on the node on which `palmrun` is executed. The real name of the hostfile will then be set to `hostfile` (instead of `auto`) and, depending on your local MPI implementation, you may have to give this name in the [execute_command](#execute_command). MPI implementations on large computer centers usually do not require to explicitly specify a hostfile (in such a case you can remove this line from the configuration file), or the batch systems provides a hostfile which name you may access via environment variables (e.g. `$PBS_NODEFILE`) and which needs to be given in the [execute_command](#execute_command). Please see your local system / batch system documentation about the hostfile policy on your system.  |
| linker_options    | Compiler options to be used to link the PALM executable. Typically, these are paths to libraries used by PALM, e.g. NetCDF, FFTW, MPI, etc. You may repeat the options that you have given with [compiler_options](#compiler_options) here. See your local system documentation / software manuals for required path settings. Requirements differ from system to system and also depend on the respective libraries that you are using.  |
| local_ip          | IP-address of your local computer / the computer on which you call the `palmrun`/`palmbuild` command. You may use `127.0.0.0` if you are running PALM in interactive mode or in batch mode on your local computer. The address is only used to identify where to send the output data in case of batch jobs on a remote host.  |
| local_jobcatalog <a name="local_jobcatalog"></a> | Folder on the local host to store the batch job protocols. In case of batch jobs running on remote hosts, the job protocol will be created on the [remote_jobcatalog](#remote_jobcatalog), and after completion of the job, it is sent via `scp` to the `local_jobcatalog`.  |
| local_username    | Your username on the local computer / the computer on which you call the `palmrun`/`palmbuild` command. The local username is required for running batch jobs on a remote host in order to allow the batch job to access your local system (e.g. for sending back output data or for automatically starting restart runs).  |
| login_init_cmd    | Special commands to be carried out at login or start of batch jobs on the remote host. You may specify here a command, e.g. for setting up special system environments in batch jobs. It is carried out as first command in the batch job.  |
| make_options      | Options for the UNIX `make`-command, which is used by `palmbuild` to compile the PALM code. In order to speed up compilation, you may use the `-j` option, which specifies the number of jobs to run simultaneously. If you have e.g. 4 cores on your local computer system, then `-j 4` starts 4 instances of the Fortran compiler, i.e. 4 Fortran-files are compiled simultaneously (if the dependencies allow for that). Do not try to start more instances than the number of available cores, because this will decrease the compiler performance significantly.  |
| memory            | Memory request per MPI process (or CPU core) in MByte.<br> **Attention:** `palmrun` option`-m` overwrites this setting.  |
| module_commands   | Module command(s) for loading required software / libraries. In case that you have a `modules` package on your system, you can specify here the command(s) to load the specific software / libraries that your PALM run requires, e.g. the compiler, the NetCDF software, the MPI library, etc. Alternatively, you can load the modules from your shell profile (e.g. `.bashrc`), but then all your PALM runs will use the same settings. An example for a Cray system to use fftw and parallel NetCDF is `module load fftw cray-hdf5-parallel cray-netcdf-hdf5parallel`. The commands are carried out at the beginning of a batch job, or before PALM is compiled with `palmbuild`.  |
| remote_ip <a name="remote_ip"></a> | IP-address of the remote system where the batch job shall be started. On large cluster systems this will usually be the address of a login node. Setting this variable in the configuration file will cause `palmrun` to run in remote batch job mode, i.e. a batch job will be created and send to the remote system automatically without giving `palmrun`-option `-b`.  |
| remote_jobcatalog <a name="remote_jobcatalog"></a> |  In case of batch jobs running on remote hosts, the job protocol will be put in this folder, and then automatically transferred back via scp to the [local_jobcatalog](#local_jobcatalog). The transfer is done by a separate small batch job, which directives are defined by the `BDT:` lines. The variable has no default value and must be set by the user. Absolute paths need to be given.<br> **Attention:** Using `$HOME` is not allowed / does not work. |
| remote_loginnode  | Name of the login node of the remote computer. Nodes on big compute clusters are separated into compute nodes and login nodes (and sometimes I/O nodes). Some computer centers only allow the login nodes to establish ssh/scp connections to addresses outside the computing center. In such cases, since `palmrun` is executed on the compute nodes, it first has to send the output data to the login node, from where it is then forwarded to your local computer. If the compute nodes on your remote host do not allow direct ssh/scp connections to your local computer, you need to provide the name of the login node of the remote host. Typically, this is a mnemonic name like *loginnode1* and not an IP-address (like *111.111.11.11*). Several login nodes often exist. You just have to give one of them. If you do not provide a name, you probably will not receive data from the PALM run on your local host.  |
| remote_username   |  Your username on the remote computer that is given by [remote_ip](#remote_ip).  |
| source_path       |  Path to PALM's Fortran source files. This is the place where the installer has installed PALM.  |
| ssh_key           |  Name of the file from which the identity (private key) for public key authentication is read. This file is assumed to be in folder `$HOME/.ssh`. By default (if you omit this variable), file `id_dsa` or `id_rsa` is used.  |
| submit_command  <a name="submit_command"></a>  |  Command that has to be used to submit batch jobs on your system (either on the local, or on the remote host), including required option. See documentation of your batch system / computing center to find out which command has to be used. An example for a `slurm` batch system could be `sbatch`.  |
| user_source_path  |  Path to the [user interface routines](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/userint). The variable `run_identifier` that may be used in the default path is replaced by the argument given with `palmrun`-option `-r`.  |

 You may add further variables to this list, which might e.g. be required for batch directives (see below). 


## Batch job directives
Running PALM in batch mode requires to add those batch directives to the configuration file that are needed by the specific batch system. Add the string `BD:` at the beginning of each directive. Because of a large variety of batch systems with different syntax, and because many computer centers further modify the directives, we can only give a general example here, which is for a `slurm` batch system om a cluster used by the PALM developers.

Batch directives required for this system read 
```
BD:#!/bin/bash
BD:#SBATCH --job-name={{ '{{job_id}}' }}
BD:#SBATCH --time={{ '{{cpu_hours}}' }}:{{ '{{cpu_minutes}}' }}:{{ '{{cpu_seconds}}' }}
BD:#SBATCH --ntasks={{ '{{mpi_tasks}}' }} 
BD:#SBATCH --nodes={{ '{{nodes}}' }}
BD:#SBATCH --ntasks-per-node={{ '{{tasks_per_node}}' }}
BD:#SBATCH --partition={{ '{{queue}}' }}
BD:#SBATCH --output={{ '{{job_protocol_file}}' }}
BD:#SBATCH --error={{ '{{job_protocol_file}}' }}
```

Strings in double curly brackets are interpreted as variables and are replaced by `palmrun` based on settings via specific `palmrun` options or settings in the environment variable section of the configuration file. From the given batch directives, `palmrun` generates a batch script (file), also called batch job, which is then submitted to the batch system using the submit command that has been defined via the variable [submit_command](#submit_command). If you like to check the generated batch script, then run `palmrun` with additional option `-F`, which will write the batch script to file `jobfile.#####` in your current working directory, where `#####` is a 5-digit random number (which is part of the so-called `run_id`). A batch job will not be submitted. 

In addition to the batch directives, the configuration file requires further information to be set for using the batch system, which is done by adding / modifying variable assignments. A minimum set of variables to be added / modified: 

```
# to be added
%submit_command      sbatch
%defaultqueue        small
%memory              1500

# to be modified
%local_jobcatalog    /home/username/local_jobcatalog
%fast_io_catalog     /gfs2/work/username
%execute_command     srun  -n {{ '{{mpi_tasks}}' }}  -N {{ '{{tasks_per_node}}' }}  ./palm
```

**Given values are just examples!** The installer may have already included these variable settings as comment lines (starting with `#`). Then just remove the `#` and provide proper values. 

## Variables used in batch job directives
The following variables are frequently used in batch directives and recognized by `palmrun` by default: 


| Variable name     | Meaning      | Value     |
|-------------------|--------------|-----------|
| cores             | total number of cores requested by the job | as given by `palmrun` option `-X`| 
| cpu_hours<br>cpu_minutes<br>cpu_seconds<br>cputime | cpu time requested by the job split in hours, minutes and seconds. `cputime` is the requested time in seconds. | calculated from `palmrun` option `-t`, e.g. in the above example, `palmrun -t 3662 ...` will generate `BD:#SBATCH --time=1:1:2`|
| job_protocol_file | name of the file (including path) to which the job protocol is written | generated from `palmrun` options `-c` and `-r` and the path set by environment variable `local_jobcatalog`. As an example, if `local_jobcatalog = /home/user/job_queue`, the call of `palmrun -r testrun -c mycluster ....` gives a job protocol file `home/user/job_queue/mycluster_testrun`.  |
| memory | requested memory in MByte | as given by `palmrun` option `-m` or as set in the configuration file via `%memory`. Option overwrites the setting in the configuration file. |
| mpi_tasks | total number of MPI tasks to be started | calculated as `cores / threads_per_task` |
| nodes | number of compute nodes requested by the job | calculated from `palmrun` options `-X` (total number of cores to be used), `-T` (number of MPI tasks to be started on each node), and `-O` (number of OpenMP threads to be started by each MPI task). `nodes` is calculated as `cores / ( tasks_per_node * threads_per_task)`. `threads_per_task` is one in pure MPI applications. If `tasks_per_node * threads_per_task` is not an integral divisor of the total number of cores, less tasks/threads will run on the last node. |
| previous_job | run_id of a previous job | as given with `palmrun` option `-W`. Can be used to define job dependencies. The run_id should be the one that has been assigned by the batch system to the (previous) job. | 
| project_account | account number under which the batch job shall run | argument of `palmrun` option `-A` |
| queue | batch queue to which the job is submitted | as given by `palmrun` option `-q`. If the option is omitted, a default queue defined by variable [default_queue](#default_queue) is used. |
| run_id | job name under which you can find the job in the respective job queue | argument of `palmrun` option `-r` plus a 5-digit random number, separated by a dot, e.g. `palmrun -r abcde ...` may generate `abcde.12345`. | 
| tasks_per_node | number of MPI tasks to be started on each requested node | as given by `palmrun` option `-T` |
| threads_per_task | number of OpenMP threads to be started by each MPI task | as given by `palmrun` option `-O` | 
| timestring | requested CPU time in format hh:mm:ss | calculated from `palmrun` option `-t` |

Instead of using variables in the batch directives, you may set fixed values, but with a loss of flexibility. 



## Additional directives for batch jobs on remote hosts

If `palmrun` is used in remote batch mode, i.e. the batch job is submitted from a local computer to a remote computer, additional batch job directives are required to guarantee that the job protocol file is sent back to the local computer after the batch job has finished on the remote system. Since the job protocol file is often only available after the job has finished, a small additional job is started at the end of the batch job, which only purpose is to transfer the job protocol from the remote to the local system. Batch directives for this job are given in the configuration file too. Add the string `BDT:` at the beginning of each directive. As for the main job directives (that start with `BD:`), only a general example is given here, which is again for a `slurm` batch system. 

```
BDT:#!/bin/bash
BDT:#SBATCH --job-name=job_transfer
BDT:#SBATCH --time=00:30:00
BDT:#SBATCH --ntasks=1 
BDT:#SBATCH --nodes=1
BDT:#SBATCH --ntasks-per-node=1
BDT:#SBATCH --partition={{ '{{queue}}' }}
BDT:#SBATCH --output={{ '{{job_transfer_protocol_file}}' }}
BDT:#SBATCH --error={{ '{{job_transfer_protocol_file}}' }}
```

Keep in mind to request only few resources because this job just carries out a file transfer via `scp`. Computing centers often offer a special queue for these kind of jobs. If so, `{{ '{{queue}}' }}` should be replaced by the specific name of that queue. The variable `job_transfer_protocol_file` is determined by `palmrun`. In case that you did not receive the job protocol, you may look into the protocol file of this transfer job. You can/should find this file under the name `last_job_transfer_protocol` on the remote host in the directory defined by [remote_jobcatalog](#remote_jobcatalog). A new job overwrites the transfer protocol of a previous job. 


## UNIX commands

The configuration file can be used to specify any Unix commands to be executed before / after the PALM code is started, or that should be carried out in case of errors that appeared during the run. Add string `IC:`, `OC:`, or `EC:` at the beginning of the respective command.

Input commands are often used for specific environment settings that are required to run PALM, e.g. 

```
IC:ulimit -s unlimited
IC:export FI_PROVIDER=psm2
IC:export FI_PSM2_CONN_TIMEOUT=20
IC:export PSM2_MEMORY=large
IC:export PSM2_MQ_RECVREQS_MAX=268435456
IC:export I_MPI_HYDRA_BRANCH_COUNT=128
IC:export I_MPI_ADJUST_ALLTOALL=3
```

On many systems, the stack size is limited by default and needs to be extended to allow for runs with larger memory requests (i.e. larger number of grid points). Settings for MPI optimization or debugging should be done here, too. Commands are carried out immediately before the PALM executable is started via the command given by [execute_command](#execute_command). 

Output commands are carried out immediately after the PALM executable has finished, and before the output files are copied. A typical example would be to send a message as a reminder that execution has finished: 

```
OC:echo "PALM simulation $run_identifier has finished" | mailx  username@email-address
```

 Error commands may be used for debugging. They are carried out if `palmrun` terminates with messages starting with`+++`, like `+++ runtime error occured`. The internal `palmrun` variable `locat` is used to specify the location within the `palmrun` script where the error appeared. It can be used to restrict the execution of the error command, e.g.

```
EC:[[ \$locat = execution ]]  &&  cat  RUN_CONTROL
```

 In this example, the content of the run control file (if there is any) is output on the terminal or to the job protocol file, in case that PALM finished with a runtime error (e.g. segmentation fault). 








