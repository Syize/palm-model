---
title: GPU mode usage
---

# GPU mode usage
<br>

## Overview
The PALM core can be run on GPUs. The code porting is an ongoing effort which started many years ago. Porting is based on OpenACC directives (see [www.openacc.org](https://www.openacc.org/)). Parts of this work has been described in Knoop et al. (2018).

**The GPU mode of PALM is an experimental feature. Generally no support can be given in case of any problems that appear when using this mode. So far, only NVidia GPUs can be used.**


## Configuration settings

See below the recommended values of variable to be set in the configuration file (for example `.palm.config.gpu`):

```
%compiler_name       mpif90
%compiler_name_ser   nvfortran
%cpp_options         -cpp -DMPI_REAL=MPI_DOUBLE_PRECISION -DMPI_2REAL=MPI_2DOUBLE_PRECISION -D__parallel -D__cuda_fft -D__netcdf -D__netcdf4 -D__netcdf4_parallel
%make_options        -j 4
%compiler_options    -O3 -Mnofma -acc=verystrict -cuda -gpu=cc80 -Minfo=accel -I \\`nf-config --includedir\\`
%linker_options      -Wl,-rpath=\\$LD_RUN_PATH \\`nf-config --flibs\\` -O3 -Mnofma -acc=verystrict -gpu=cc80 -cuda -cudalib=cufft
%execute_command     mpirun -np {{mpi_tasks}} --map-by ppr:2:socket:pe=1 ./palm
```

These settings are for NVidia GPUs. The netCDF settings in the compiler- and linker options may require adjustments. Also the `execute_command` may need adjustments depending on the requirements of your system. Please note that we can not give any support for creating a configuration file that is working on your system.


## Namelist parameter settings

Some specific settings for namelist parameters are required / recommended.

- For [psolver](../../../Reference/LES_Model/Namelists/#initialization_parameters--psolver)=  *'poisfft'*, only the Temperton-FFT method and the FFTW are available. Set [fft_method](../../../Reference/LES_Model/Namelists/#initialization_parameters--fft_method) = *'temperton-algorithm'* or [fft_method](../../../Reference/LES_Model/Namelists/#initialization_parameters--fft_method) = *'system-specific'*. The latter setting invokes the CUDA FFTW library available on the GPU. It usually should show a much better performance than the Temperton-FFT. Be aware that CPP options `-D__cuda_fft` has to be set in the configuration file.
- To speed up the MPI communication when using multiple GPUs, set runtime parameter [use_contiguous_buffer](../../../Reference/LES_Model/Namelists/#runtime_parameters--use_contiguous_buffer) = *.TRUE.*. This usually improves the performance of ghost point exchange significantly.


## Running with palmrun

The GPU mode assumes that one GPU is attached to one MPI process, so the total numer of MPI processes must match the number of GPUs to be used on the system. In case of a cluster system equipped with 4 GPU cards per node, the `palmrun` command for running on 12 GPUs must be
```
palmrun -c gpu ...... -X12 -T4
```
meaning that PALM will be executed on 3 nodes.


## Performance issues

So far, no special focus has been given on performance optimization of the GPU mode. Our experience from test runs is that the price / performance ratio of state-of-the-art CPU and GPU is quite similar, but it should be noted that this heavily depends on the chosen setup. Carefully analyze the CPU time measurements (file in `MONITORING` folder with suffix `_cpu`) for any bottlenecks.


## Notes, shortcommings and open issues

- PALM can be run with single precision (32 bit floats) on GPU, too.
- There are plans to support usage of GPU devices via openMP5 offloading.

The GPU mode has the following restrictions:

- Beside the PALM core (dynamics and thermodynamics, except cloud physics) and the land surface model (LSM) no other modules have been ported, and there are no current plans for further porting.
- Both the direct Poisson-solver as well as the multigrid solver are available.
- Even when using just the PALM core, several settings may not work appropriately. Always compare GPU results with results from control runs carried out on CPUs (an easy way to run a setup purely on the CPU side is to set runtime parameter [enable_openacc](../../../Reference/LES_Model/Namelists/#runtime_parameters--enable_openacc) = *.FALSE.*.). Small differences in the run-control (`_rc`) output may always appear due to different round-off errors on the GPU. Because of the non-linear turbulence interactions, instantaneous flow fields may completely differ after some time, but averaged quantities should not be affected (provided that the averaging interval is long enough). 
- Some standard output quantities may not be available or the netCDF output files may contain wrong values. Please check any output carefully.


## References

**Knoop, H., T. Gronemeier, M. Sühring, P. Steinbach, M. Noack, F. Wende, T. Steinke, C. Knigge, S. Raasch, and K. Ketelsen. (2018):** Porting the MPI-parallelized LES model PALM to multi-GPU systems and many integrated core processors - an experience report,  Int. J. Computational Science and Engineering, 17(3), 297–309.


## Acknowledgements

The PALM developers acknowledge support from **natESM**, the national Earth System Modelling Strategy (funded by the German Federal Ministry of Education and Research, BMBF, grant no. 01LK2107A1), which provided Research Software Engineering support for this work.
