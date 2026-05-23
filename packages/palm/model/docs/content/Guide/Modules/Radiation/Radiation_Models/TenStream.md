# TenStream

Coupling TenStream solver to PALM model system

---

## Introduction

PALM-4U uses RRTMG as 1-D two-stream approach to provide the atmospheric heating and surface heating rates. The Radiative Transfer Model (RTM) is coupled to RRTMG to account for the three-dimensional radiation (3D) effects in the urban environment. While this coupled approach is computationally efficient, it is not capable to consider the effects of 3-D radiative transfer on atmospheric heating rates or dynamic heterogeneities, such as moving clouds or fog. Moreover, this coupled approach has limitations to estimate the atmospheric longwave radiation cooling of the near-surface air during nighttime ([Nunez and Oke, 1976](https://doi.org/10.1007/BF00229280); [Steeneveld et al., 2010](https://doi.org/10.1029/2009JD013074)).

To consider the 3-D radiation, an alternative approach using the TenStream solver. This scheme calculates more than ten streams for
each atmospheric grid cell. The TenStream considers the 3-D propagation of radiation in the atmosphere. The coupled model TenStream/PALM correctly considers interactions with atmospheric constituents in the urban canopy layer and above, such as water vapor, fog, and clouds.

## The TenStream solver

The TenStream radiative transfer model is a parallel approximate solver for the full 3-D radiative transfer equation. Unlike the traditional two-stream solver, e.g. RRTMG, it approximately solves the radiative transfer equation in 3-D and computes irradiance and heating rates from optical properties. Additional to the up- and downward fluxes, TenStream solver computes the radiative transfer coefficients for sideward streams.

Overview and detailed description of the TenStream solver are available at [Jakub and B. Mayer, 2015](https://doi.org/10.1007/BF00229280) and [Jakub and B. Mayer, 2016](https://doi.org/10.5194/gmd-9-1413-2016). However, the TenStream code is recently extended to consider the urban environment by considering the different urban surfaces as well as the vegetation.

## Installation

To install TenStream, prerequisites such as MPI, NetCDF libraries, and cmake are required. Additionally the Portable, Extensible Toolkit for Scientific Computation, PETSc ([Balay et al., 2014](https://doi.org/10.2172/1178109)), which offers a wide range of composable iterative solvers and matrix preconditioners, is required.

### Install with palm-cli

The most convenient method to install a PALM version with TenStream support is via the PALM model install script which will download and install PETSc and TenStream as part of the installation procedure and generate a PALM config file with appropriate options.
The palm {{ link_palm_repo_file('install', 'install') }} script, shall be invoked with the option `-T`, i.e. the Tenstream/PALM model can be installed with the following commands:

``` bash
export FC=mpif90
export CC=mpicc
export CXX=mpicxx
export install_prefix="<install-prefix>"
bash install -T -p ${install_prefix} -c $FC
export PATH=${install_prefix}/bin:${PATH}
```

please replace `<install-prefix>` with the desired installation directory.

!Note that the `TenStream submodule` will be cloned to `lib/tenstream/tenstream`, i.e. a script to download TenStream lookup tables is available at `lib/tenstream/tenstream/misc/download_LUT.sh`.


### Install manually and link with palmconfigure

The TenStream module can also be installed on its own as a standalone library and later linked against PALM.

#### Quick install guide of TenStream/standalone with gcc and openmpi

The source code of the TenStream solver is available for download at <https://github.com/tenstream>.
The following recipe may help you get started and the location of the TenStream and PETSc install will be referenced in the next step for the purpose of having an example installation.

``` bash
export FC=mpif90
export CC=mpicc
export CXX=mpicxx

export TENSTREAMROOT=$HOME/tenstream
git clone https://github.com/tenstream/tenstream.git $TENSTREAMROOT
cd $TENSTREAMROOT

# Build PETSc, LAPACK, HDF5 and NetCDF
misc/build_dependencies.sh
export PETSC_DIR=$(pwd)/external/petsc
export PETSC_ARCH=default

# build TenStream
mkdir build
cd build
cmake ..
make -j

# download TenStream lookup tables needed when running TenStream
misc/download_LUT.sh direct_3_10
misc/download_LUT.sh diffuse_10

# Run example from TenStream test suite
ctest -V -R rrtm_lw_sw
```

#### Link PALM against preinstalled TenStream library

To install the coupled Tenstream/PALM model, the wrapper code in PALM must be first activated by adjusting the palm configuration file
(`.palm.config file`) as follows:

-   Add the string `-D__tenstream` to the C pre-processor field (`%cpp_options`) in the palm configuration file. For example, this field may look like:

``` bash
%cpp_options        -cpp -D__gfortran -D__parallel -DMPI_REAL=MPI_DOUBLE_PRECISION -DMPI_2REAL=MPI_2DOUBLE_PRECISION -D__netcdf -D__netcdf4 -D__netcdf4_parallel -D__tenstream
```

-   Add the path of TenStream's library to the compiler options field (`%compiler_options`) as well as the linker options(`%linker_options`), as follow

``` bash
%tenstream_include   <path to the include folder of tenstream e.g. $TENSTREAMROOT/build/include>
%tenstream_lib       <path to the library folder of tenstream e.g. $TENSTREAMROOT/build/lib    >
%petsc_lib           <path to the library folder of petsc     e.g. $PETSC_DIR/$PETSC_ARCH/lib  >
%compiler_options    <...> -I$tenstream_include
%linker_options      <...> -L$tenstream_lib -ltenstream -L$PETSC_DIR/$PETSC_ARCH -lpetsc -lflapack -lfblas
```

## Job preparation

In order to run a case with TenStream as the radiation model, you need to provide the following:

-   **Look-Up-Tables (LUT)**. It is important to note here that the default tenstream solver is two-streams, i.e. one-dimension, (`2str`, see [TenStream solver options](#tenstream-solver-options) Section) and, if it is kept the same, you need **not** to provide LUT. However, if you want to include the main features of tenstream concerning the 3-D radiation interactions, you need to change the tenstream solver option and hence you will need to provide the LUT.

The LUT contains the voxel radiation-transport coefficients required for performing the radiative transfer processes. The solver expects the path to the LUT's location to be located at the variable `lut_basename`, which is by default the temporary folder created by PALM for each run. Calculating the LUT can be accomplished by running the following command
``` bash
mpirun ~/tenstream/build/createLUT 3_10
```
However, computing the LUT is computationally expensive and is only useful when custom details of the transfer coefficients are required.
You probably don't want to do this but rather download precomputed tables.

The default LUTs are available at <https://www.meteo.physik.uni-muenchen.de/~Fabian.Jakub/TenstreamLUT/>. The LUTs for 3 direct and 10 diffuse streams can be downloaded using these commands:
``` bash
<tenstream_source>/misc/download_LUT.sh direct_3_10
<tenstream_source>/misc/download_LUT.sh diffuse_10
```
It is recommended to store the LUT at a common location and provide its path via the `-lut_basename` commandline option or via an options file (see [TenStream solver options](#tenstream-solver-options) Section).
Note that if you installed with the PALM Model install script, the `TenStream submodule` will have been cloned to `lib/tenstream/tenstream`, i.e. a script to download TenStream lookup tables is available at `lib/tenstream/tenstream/misc/download_LUT.sh`.


-   **Background atmosphere file**. This is an ASCII file contains the background information of the atmosphere. It has the name `<run_identifier>_ts_back_atm` and contains the following columns:

``` bash
# AFGL atmospheric constituent profile. U.S. standard atmosphere 1976. ( AFGL-TR-86-0110)       
#     z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)
  30.000 1.197000e+01 2.265000e+02 3.827699e+17 2.509799e+12 8.004700e+16 1.809675e+12 1.263900e+14 2.359280e+09
  27.500 1.743000e+01 2.240000e+02 5.635873e+17 3.272892e+12 1.178760e+17 2.580300e+12 1.861200e+14 2.712840e+09
  24.000 2.972000e+01 2.206000e+02 9.757872e+17 4.518265e+12 2.040885e+17 4.198950e+12 3.222450e+14 2.988090e+09
  22.000 4.047000e+01 2.186000e+02 1.340895e+18 4.894274e+12 2.804780e+17 5.455230e+12 4.428600e+14 2.898720e+09
  20.000 5.529000e+01 2.167000e+02 1.847990e+18 4.768571e+12 3.864410e+17 7.211100e+12 6.101700e+14 2.570110e+09
  18.000 7.565000e+01 2.167000e+02 2.528494e+18 4.015110e+12 5.287700e+17 9.677249e+12 8.349000e+14 1.950630e+09
  16.000 1.035000e+02 2.167000e+02 3.459340e+18 3.012286e+12 7.235580e+17 1.367490e+13 1.142460e+15 1.104378e+09
  14.000 1.417000e+02 2.167000e+02 4.736121e+18 2.383717e+12 9.904510e+17 2.808805e+13 1.563870e+15 3.544772e+08
  12.000 1.940000e+02 2.167000e+02 6.484174e+18 2.008345e+12 1.356201e+18 1.236803e+14 2.141370e+15 2.044035e+08
  10.000 2.650000e+02 2.233000e+02 8.595457e+18 1.129443e+12 1.797818e+18 6.017959e+14 2.838660e+15 2.047276e+08
   9.000 3.080000e+02 2.297000e+02 9.711841e+18 8.910379e+11 2.031271e+18 1.538518e+15 3.207270e+15 2.254808e+08
   8.000 3.565000e+02 2.362000e+02 1.093179e+19 6.526804e+11 2.286460e+18 4.011698e+15 3.610200e+15 2.516200e+08
   7.000 4.111000e+02 2.427000e+02 1.226845e+19 6.151052e+11 2.566520e+18 7.024160e+15 4.052400e+15 2.824400e+08
   6.000 4.722000e+02 2.492000e+02 1.372429e+19 5.645776e+11 2.869570e+18 1.270574e+16 4.530900e+15 3.157900e+08
   5.000 5.405000e+02 2.557000e+02 1.531006e+19 5.772576e+11 3.201880e+18 2.140204e+16 5.055600e+15 3.523600e+08
   4.000 6.166000e+02 2.622000e+02 1.703267e+19 5.771448e+11 3.561360e+18 3.677232e+16 5.623200e+15 3.919200e+08
   3.000 7.012000e+02 2.687000e+02 1.890105e+19 6.274337e+11 3.952190e+18 6.017162e+16 6.240300e+15 4.349300e+08
   2.000 7.950000e+02 2.752000e+02 2.092331e+19 6.778279e+11 4.376460e+18 9.697315e+16 6.910200e+15 4.816201e+08
   1.000 8.988000e+02 2.817000e+02 2.310936e+19 6.779402e+11 4.834170e+18 1.404222e+17 7.632900e+15 5.319900e+08
   0.000 1.013000e+03 2.882000e+02 2.545818e+19 6.777680e+11 5.325320e+18 1.973426e+17 8.408400e+15 5.860400e+08
```

-   **TenStream options file**. Optionally, all the TenStream options (see [TenStream solver options](#tenstream-solver-options) Section) which are required to tune the solver for the current run should be added to a separate file. The file name should be `<run_identifier>_ts_options`.

-   **PALM Input/Output configuration file**. The input/output configuration file `.palm.iofiles` file should be adjusted to consider the input files needed to run TenStream as follow
``` bash
TS_BACKGROUND_ATM        inopt:tr   d3#:d3r  $base_data/$run_identifier/INPUT          _ts_back_atm
tenstream.options        inopt:tr   d3#:d3r  $base_data/$run_identifier/INPUT          _ts_options
```

-   **Radiation NAMELIST**. The NAMELIST for radiation in `<run_identifier>__p3d` file can optionally contain Tenstream related variables as follows:
``` fortran
 &radiation_parameters radiation_scheme = 'tenstream',
```
The debug messages of Tenstream will be reported automatically when the PALM global debug flag, i.e. `debug_output`, is set to `.TRUE.`.<br/>


## TenStream solver options

The TenStream radiation package allows to set a lot of solver options via runtime options which are given either in the local run directory in a file named `tenstream.options`
or in a global configure file located in `$HOME/.petscrc`.
Listing all available options here is out of the scope of this readme but the interested user is referred to the codebase.
``` bash
grep -r get_petsc_opt <tenstream_source_dir>
```
These options also allow you to steer the iterative matrix solvers within the PETSc library which may be necessary to achieve good parallel scaling if you run large simulations (the default solver is robust but does not scale well for high resolution simulations).

The following list comprises some options to get you started.

| Option name | Description |
| ---      | ---      |
| -ts_options_view | print available options and used values to STDOUT during runtime |
| -lut_basename | path prefix to where the lookup tables lie(has to be accessible on all compute nodes) |
| -specint | select one of the [Spectral Integration Schemes](#tenstream-spectral-integration-schemes) |
| -log_view | output performance timings at the end of the simulation |
| -solver 2str | use a 1D twostream solver to compute radiative transfer (default, fast but not 3D) |
| -solver 3_10 | 3 direct streams and 10 diffuse streams. This is the minimum set of streams to get reasonable 3D results |
| -pprts_1d_height 500 | treat all layers above 500m as 1D layers, keeping down the number of 3D layers makes the simulation cheaper but neglects horizontal transport |
| -pprts_view_geometry | output info about the geometry of the domain |
| -pprts_view_suninfo | output info on sun angles |
| -solver 3_10 // -solver 3_24 | select a different solver from commandline, to list available solvers, use -solver -1|
| -solar_dir_ksp_view | info on the solver setup for solar direct radiation |
| -solar_diff_ksp_view | info on the solver setup for solar diffuse radiation |
| -thermal_diff_ksp_view | info on the solver setup for thermal diffuse radiation |
| -thermal_diff_ksp_monitor | track residual of the solution process |
| -thermal_diff_ksp_converged_reason | print number of iterations needed |

Please notice the prefix of `-solar_dir_`, `-solar_diff_` and `-thermal_diff_` which can be used to select an iterative solver from the PETSc suite.
The possibilities are endless in that regard and finding a robust and scalable option is hard.
If you find the default solver to converge too slowly and want to experiment, I encourage you to get in touch with a TenStream developer and have a look at the PETSc web page and manual.

#### Advanced options

###### Multigrid
One example for a rather advanced PETSc solver would be to use Multigrid Preconditioning with GAMG which is slower for small simulations because it has quite some setup cost but if you find that the default solver takes many iterations to converge and you are running a high resolution simulation, this may help to improve performance:
E.g. use flexible GMRES with algebraic multigrid preconditioning. Use unsmoothed aggregation and drop coefficients below .2. On each coarse level use 5 iterations of SOR. Will coarsen problem until it solves it with direct LU on single processor...
```
 -solar_diff_ksp_type fgmres
 -solar_diff_pc_type gamg
 -solar_diff_pc_gamg_type agg 
 -solar_diff_pc_gamg_agg_nsmooths 0
 -solar_diff_pc_gamg_threshold 1e-3     # [.0-.2]
 -solar_diff_pc_gamg_sym_graph true
 -solar_diff_mg_levels_ksp_type richardson
 -solar_diff_mg_levels_ksp_max_it 5   # [1-8]
 -solar_diff_mg_levels_pc_type sor
```
You may play around a bit with the numbers in brackets to suite your network and/or problem.

###### Incomplete solves

Another option to speedup the solve is to allow for incomplete solves.
E.g. run only 2 iterations of SOR on the diffuse radiation.
This may introduce non-negligible errors in case that the scene is rapidly changing between two radiation calls.
However, if you call the radiation often enough we found that it does not introduce any biases but this depends on your simulation setup.
If you have any questions, don't hesitate to ask a TenStream dev.
Options could be something along:

```
 -solar_diff_explicit                   # use an explicit SOR solver, circumventing setup costs of the PETSc matrix
 -solar_diff_ksp_max_it 2               # run only 2 iterations of SOR
 -solar_diff_ksp_skip_residual          # dont evaluate residual, we run a fixed nr or iterations anyway
 -solar_diff_ksp_ignore_max_it  60      # ignore the max_it setting in the first 60 seconds of the simulation (i.e. do full solves for spinup)
 
 -thermal_diff_explicit                 # use an explicit SOR solver, circumventing setup costs of the PETSc matrix
 -thermal_diff_ksp_max_it 2             # run only 2 iterations of SOR
 -thermal_diff_ksp_skip_residual        # dont evaluate residual, we run a fixed nr or iterations anyway
 -thermal_diff_ksp_ignore_max_it  60    # ignore the max_it setting in the first 60 seconds of the simulation (i.e. do full solves for spinup)

 -accept_incomplete_solve               # dont throw an error if we havent fully converged
```


###### TenStream spectral integration schemes

For the purpose of atmospheric radiative transfer computations, we need to integrate over a wide range of the electromagnetic spectrum, i.e. the shortwave (ultra-violet to near infrared) and longwave (near infrared to cm wavelengths).
The highly non-linear absorption lines of participating trace-gases (often called gas optical properties) necessitates billions of monochromatic (one wavelength considered) radiation calls.
So called line-by-line computations are not feasible to do in a LES.
To reduce the number of radiation calls various methods exist to compute integrated radiative fluxes, the most prominent class being so called correlated-k methods of which most atmospheric models today use the **RRTMG** spectral integration scheme which accurately considers even minor trace gases (e.g. CFCs), relevant for climate modelling studies.
**RRTMG** for example only uses 260 monochromatic radiation calls to compute shortwave and longwave integrated fluxes.
Recent developments such as [**ECCKD**](https://doi.org/10.1029/2022MS003033) aim to reduce the number of monochromatic radiation calls even further e.g. by targeting only the most prominent trace-gases relevant to NWP and LES applications.
For most LES applications, we recommend to use the **ECCKD** spectral integration scheme which uses 32 monochromatic radiative transfer computations, yielding high accuracy at a speedup of about a factor of x4 over **RRTMG**.

The TenStream radiation package supports multiple spectral integration schemes which can be selected through TenStream options.

| Option name | Description |
| ---      | ---      |
| -specint ecckd  | [ECCKD (32 bands)](https://doi.org/10.1029/2022MS003033) |
| -specint rrtmg  | [RRTMG (260 bands)](http://rtweb.aer.com/rrtm_frame.html) |
| -specint repwvl | [representative wavelength (30 bands)](https://doi.org/10.3390/atmos14020332) |


###### Adaptive spectral integration

**Only relevant for `-specint RRTMG`)**
Each time the radiation scheme is called the solver integrates over the solar and thermal spectral range.
I.e. at each timestep the solver is called 260 times.
However, there are some spectral intervals that only vary slowly over time.
E.g. UV radiation that is already removed in high altitudes does not change near the surface over the course of a couple of timesteps.
One thing we can do is to estimate (fit an exponential to initial solver residuals) the error growth over time that happens in a spectral interval.
Then, if the estimate of change is small we use the values from the last timestep.

E.g. lets assume that we have a PALM radiation timestep of 30s then we could use options
```
-max_solution_time 300  # make sure that every band is computed every 5 minutes, irrespective of the estimated error
-max_solution_err 0.1   # skip the spectral band if the maxnorm error is estimated to be smaller
```

In case that you have a clearsky simulation, this may greatly reduce the number of solver calls.
But check that this is a viable option for your setup.
