---
title: IO Configuration File
---
#The File Connection File

---

PALM expects all input files to be in the temporary working directory from where it is executed. Also, it writes all output data to this directory. The file connnection file informs `palmrun` where it can find the user's input files, and where the output files shall be copied after the end of a run. The default file connnection file used by `palmrun` is stored in the installation folder under `.../packages/palm/model/share/config/.palm.iofiles`. In general, the file connection file does not need any modifications.

## Format of the file connection file

The file is in ASCII format (see [here](https://gitlab.palm-model.org/palm/model/-/blob/master/share/config/.palm.iofiles)) and it generally contains one line per input/output file, which is called a ***file connection statement***.

A file connection statement consists of a maximum of six columns. The first five columns are obligatory. The meaning of the columns are explained below based on the following example, which shows an excerpt of the default file:
```
#---------------------------------------------------------------------------------------------
# List of input-files
#---------------------------------------------------------------------------------------------
PARIN                    in:tr      d3#      $base_data/$run_identifier/INPUT          _p3d*
PARIN                    in:tr      d3r      $base_data/$run_identifier/INPUT          _p3dr*
TOPOGRAPHY_DATA          inopt:tr   d3#:d3r  $base_data/$run_identifier/INPUT          _topo*
DATA_1D_PR_NETCDF        inopt      prr      $base_data/$run_identifier/OUTPUT         _pr*    nc
#
#---------------------------------------------------------------------------------------------
# List of output-files
#---------------------------------------------------------------------------------------------
BINOUT*                  out:lnpe   restart  $fast_io_catalog/$run_identifier/RESTART  _d3d
DATA_3D_NETCDF*          out:tr     *        $base_data/$run_identifier/OUTPUT         _3d     nc
```

The **first column** gives the local filename in the temporary working directory. This is the name which appears in the corresponding Fortran OPEN statement. PALM generally expects all filenames in capital letters.

The **second column** gives the file attributes, to be separated by a colon (":") if more than one attribute is given. Two attributes are allowed. The **first attribute** defines, whether the file is an input (`in`) or an output file (`out`). `palmrun` immediately terminates, if it cannot find the input file, because it assumes that PALM cannot run without files with attribute `in`. Optional input files is shall not terminate a run if they are missing have the attribute `inopt`. The following options are allowed for the **second attribute**:

- `tr`<br>
In case of running jobs on remote hosts, transfer the file between the remote and local host via `scp`.
Input files will be transferred from the user's folder on the local host to the temporary working directory on the remote host, before the model is executed. This transfer is actually carried out by `palmrun` before the batch job to the remote host is submitted. Files will be stored on the remote host in folder `$fast_io_catalog/SOURCES_FOR_RUN_$run_identifier`. As soon as the batch job executes, files are copied from the `SOURCES_FOR_RUN`-folder to the temporary working directory.
Output files will be transferred from the temporary working directory on the remote host to the user's folder on the local host, after the model has finished execution. When running `palmrun` in interactive mode, `tr` is ignored and the file is copied on the local machine. The same holds for batch jobs on the local machine.

- `tra`<br>
Can only be used for ASCII output files. Same as `tr`, but the local file from the remote host will be appended to the file with highest cycle number, if a respective file already exists in the user's folder on the local host. In such a case, no new file cycle will be created.

- `pe`<br>
It is checked if there is a folder with the given name. The name is then interpreted as a folder name, and it is assumed that the respective folder contains multiple files, one for each MPI thread. Filenames in this folder are expected as `_######`, where `######` is the six-digit MPI-thread id. If the model is running on 4 cores (4 MPI threads), file names will be `_000000`, `_000001`, `_000002`, and `_000003`. `palmrun` will copy the complete folder. If there is no folder with that name, it is assumed to be a file.

- `lnpe`<br>
Same as `pe`, but `palmrun` tries to move the files using `ln` instead of `cp`, which reduces the time to copy the file drastically, but only if both folders (the local temporary folder and the user's folder) reside on the same file system. This attribute is given for restart files in the default file connection file, because they may have a huge size.

- `trpe`<br>
Same as `pe`, but the folders are transferred between the local and the remote host via `scp`.

The **third column** defines the so-called activation strings. Two or more strings need to be separated by colons (`:`). Activation strings are those given via the `palmrun` option `-a`. The respective file connection statement is only processed (i.e. the file is copied), if at least one of its activation strings is set via option `-a`. In case of the above excerpt of the default file connection file, `-a "d3#"` causes `palmrun` to provide the input files `PARIN` and `TOPOGRAPHY_DATA`. Option `-a "d3# restart"` additional copies/saves the output file `BINOUT`. Note that activation strings need to be separated by blanks when given via the option `-a`, and that they need to be embraced by `"`, even if only one string is given (otherwise, a `#` would cause the shell to interpret the remaining part of the command as a comment). For output files, a wildcard may be given in the third column. Respective local files will then always be copied, in case they exist, independently from the activation strings given with option `-a`. No warning will be issued if no file exists. Wildcards in the third column are not allowed for input files.

The path to the user's permanent file is given in the **fourth column**. It may contain environment variables that have been defined in the `palmrun`/`palmbuild` configuration file. `$run_identifier` will be replaced with the argument of `palmrun` option `-r`.

The **fifth column** gives the suffix of the user's permanent file.

The **sixth column** is optional and must only contain the string `nc`. This string has to be given if the user's permanent file is a netCDF file with suffix `.nc`.


## Determination of full file paths

The full path for input and output files is based on columns 4-6 and the argument of `palmrun` option `-r` 
```
<column 4>/<-r option><column 5>[.<column 6>]
```
As an example, the file connection statement
```
DATA_3D_NETCDF*    out:tr   *    $base_data/$run_identifier/OUTPUT   _3d   nc
```
plus the setting in the configuration file
```
%base_data         ~/palm/current_version/JOBS
```
plus the call `palmrun -r example_test ...` will create the output file
```
~/palm/current_version/JOBS/example_test/OUTPUT/example_test_3d.nc
```

## Using wildcards in filenames

In case of nested runs, masked 3d-output, wind turbine model output, or coupled atmosphere-ocean runs, multiple input/output files with additional suffixes are expected or created. As an example, a nested run with a parent and one child requires two parameter files for steering the simulation, where the second file is expected to have the suffix `_N02`. Each additional nest would require another parameter file with suffixes `_N03`, `_N04`, etc. In order to simplify the file connection file, only one file connection statement
```
PARIN     in:tr   d3#    $base_data/$run_identifier/INPUT     _p3d*
```
is required for this, where the wildcard in `_p3d*` causes `palmrun` to copy files with these suffixes too. Assuming an INPUT folder with files
```
example_test_p3d
example_test_p3d_N02
example_test_p3d_N03
```
`palmrun` will copy them to local files
```
PARIN
PARIN_N02
PARIN_N03
```
The same holds for output files. Here the wildcard must be given as part of the local filename (first column). The example given below assumes one nest (root domain plus one child domain), two output masks for the root domain, and three masks for the child. Assuming a file connection statement for an output file
```
DATA_MASK_NETCDF*   out:tr   *  $base_data/$run_identifier/OUTPUT   _masked   nc
```
and following local files that have been created by PALM
```
DATA_MASK_NETCDF_M01
DATA_MASK_NETCDF_M02
DATA_MASK_NETCDF_N02_M01
DATA_MASK_NETCDF_N02_M02
DATA_MASK_NETCDF_N02_M03
```
`palmrun` will copy these local temporary files to permanent files
```
....../OUTPUT/example_test_masked_M01.nc
....../OUTPUT/example_test_masked_M02.nc
....../OUTPUT/example_test_masked_N02_M01.nc
....../OUTPUT/example_test_masked_N02_M02.nc
....../OUTPUT/example_test_masked_N02_M03.nc
```
Instead of using wildcards in filenames, separate file connection statements may be used, one for each possible suffix.


## How to modify the file connection file

In specific cases the file connection file may need to be modified. For example, if a user-interface requires or creates additional input/output files, additional file connection statements need to be add to the file connection file. In such a case, the default file should be copied to the working directory, e.g.
```
cd ~/palm/current_version
cp release25.04/packages/palm/model/share/config/.palm.iofiles  .
```
and only the copied file shoulld be midified. `palmrun` always checks, if there is a `.palm.iofiles` in the working folder and will automatically use this file.

Beside adding new file connection statements, also existing statements may be modified, e.g. by defining new activation strings, or by changing the file paths given in the fourth column.









