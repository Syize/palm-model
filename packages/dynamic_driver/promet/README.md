# PROMET

PROMET is a tool that helps you create dynamic drivers for the PALM Model. It is currently able to work with ICON-D2/ICON-ART data as well as WRF/WRF-Chem data. (If you want to use COSMO data, please for now still use the tool INIFOR)

## Install

In order to install PROMET, please open the source directory and type the following two commands:

```bash
pip install -r requirements.txt
pip install .
```

## Usage

In order to run PROMET, you need a working PALM setup. PROMET reads the namelist file and the static driver file of the root domain of your PALM setup. Additionally, PROMET takes mesoscale output data from ICON-D2/ICON-ART or WRF/WRF-Chem as input from which the created dynamic driver will be populated. Please make sure the mesoscale data includes data about its grid. Finally, PROMET requires a dedicated config file in which you can specify what dynamic variable groups you would like to use and also details about the storage location of the mesocale data, the geographic projections of the data and the time period the created dynamic driver is supposed to cover.

PROMET has an extensive user-error detection logic, that you can use in order to be pushed towards a correct config file. The namelist file and the static driver file are specified using command line arguments (`--namelist` and `--static-driver`). Also the desired output path of the created dynamic driver can be specified using a the argument `--output-file`. In case you need an overview of all the available command line arguments you can type `promet --help`.

## Running the included testcase

A good starting point for using PROMET are the included testcases at `tests/integration`. There is one example for ICON and one for WRF. Please note that the actual mesoscale data files are not part of this repository. 

If you want to use the testcase for ICON, please make sure the directory `tests/integration/icon/input_data` contains your ICON output files in netcdf format that contain data on the respective original mesoscale model grid.

You can obtain this data for example via the pamore interface provided by the DWD. Here is an example command to request data via pamore:

```bash
pamore -G -F -d 2023070600 -hstart 0 -hstop 6 -hinc 1 -tflag best -ee P%HI,QC%HI,QV%HI,T%HI,U%HI,V%HI,W%H,HHL%H,W_SO%B,T_SO%B,CLAT%G,CLON%G,HSURF%G,PS%G -ires r19b07 -model ilam -ofmt netcdf
```

After storing/linking the data to the directory `tests/integration/icon/input_data` (make sure to unpack the date and if needed convert it to netcdf), you can execute the testcase as follows:

```bash
cd tests/integration/icon/
promet --config promet_test_01.yml --namelist promet_test_01_p3d --static-driver promet_test_01_static --output-file promet_test_01_dynamic --verbose --overwrite
```

## Capabilities Overview

PROMET uses mesoscale output data from ICON-D2/ICON-ART or WRF/WRF-Chem as input in order to create dynamic drivers for PALM. This is achieved by tranforming and interpolating the data. First, coordinate system transformations are carried out between the individual geodetic datum and projections of the input data and the desired PALM grid. Then the input data, which are assumed to be unstructured, are tessellated (triangulated) and a linear barycentric interpolation is performed on each triangle. This way, any structured, rotated and unstructured input data can always be interpolated onto the PALM grid in the same way. The interpolated data is then vertically compressed/stretched according to the difference between input orography and palm orography in order to achieve optimal vertical mapping and the most correct representation of vertical profiles close to the ground. Apart from that, a whole series of operations for data management and mapping of variables takes place internally to make sure the correct data from the input data set is baked into the newly created dynamic driver for PALM.


# License and Privacy Policy

This repository is free and open-source software licensed under the [AGPLv3](LICENSE). By downloading and executing this copy of PROMET, you agree to interact with the servers that host the software repository and downloadable files, the web front-end with its underlying api services as well as with the servers and apis that provide the content for the documentation, enable the software to look for new updates, send error reports to request extended error explanations based on internal error codes and request configuration- and setup-enhancement suggestions based on the setup used. If you do not register or otherwise transmit personal information, these servers collect the personal data that your browser or PROMET transmits by default. In this case the servers collect the following data, which is technically necessary to display websites to you, answer your api-requests and to ensure stability and security of the servers. The collected data contains IP address, date and time of the request, time zone difference to Greenwich Mean Time (GMT), content of the request (specific page URL, HTTP-header and body including data like internal error codes and configuration and setup parameters and their values transmitted in api-queries), access status/HTTP status code, amount of data transferred for each request, URL from which the request originates, browser/client software, operating system and its interface, language and version of the browser/client software and version of PROMET. This data is stored in log files on the servers, which are only accessible to the administrator of the servers. It is not stored together with other personal data. The legal basis for the temporary storage of data is Art. 6(1)(f) of the GDPR. The temporary processing of this data by the system is necessary to allow the servers to communicate with your terminal device and to answer your requests. For this purpose, your IP address must remain stored for the duration of the session. The storage in log files is done to ensure the functionality of the servers. In addition, the data is used to optimize the websites and api-services you used and to ensure the security of the underlying information technology systems. An evaluation of the data for marketing purposes does not take place. The data will be deleted as soon as they are no longer required to achieve the purpose for which they were collected.
