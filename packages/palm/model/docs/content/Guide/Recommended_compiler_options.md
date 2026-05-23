---
title: Recommended compiler options
---

# Recommended compiler options
<br>

The PALM code uses standard FORTRAN 2008 and should compile on any compiler conforming with this standard.


## General recommendation

Compilers currently used by the PALM group are GNU, Intel, Cray, and NVidia-Fortran. We recommend to use these compilers with options as given below (library and include-file options are omitted).

**Please always enable floating-point exceptions to abort runs.** Otherwise, runs will continue in case of floating point errors, results will contain `NaN`, and it will be impossible to determine the exact code location where the first error happened. If floating-point exceptions are not enabled, runs may abort with netCDF errors when `NaN` data is output. Such errors are only follow-up errors and do not give any information about the real cause of the problem.

**Please always switch on traceback options.** They are required to locate the code line where errors appear. Traceback options should also be given for optimized code, because they don't affect performance.


## GNU Fortran:

`gfortran -Ofast -g -ffpe-trap=invalid,zero,overflow -fbacktrace`

| Option | meaning |
|--------|---------|
| -Ofast | Disregard strict standards compliance. -Ofast enables all -O3 optimizations. It also enables optimizations that are not valid for all standard-compliant programs. It turns on -ffast-math and the Fortran-specific -fstack-arrays, unless -fmax-stack-var-size is specified, and -fno-protect-parens. |
| -g     | Produce debugging information in the operating system's native format. |
| -ffpe-trap= | Specifies a list of IEEE exceptions when a Floating Point Execption (FPE) should be raised. *invalid*: invalid floating point operation, such as SQRT(-1.0), *zero*: division by zero, *overflow*: overflow in a floating point operation. |
| -fbacktrace | Specify that, when a runtime error is encountered or a deadly signal is emitted (segmentation fault, illegal instruction, bus error or floating-point exception), the Fortran runtime library should output a backtrace of the error. |


## GNU Fortran for debugging:

`gfortran -O0 -g -ffpe-trap=invalid,zero,overflow -Wall -Wextra -pedantic -fcheck=all -fbacktrace`

| Option | meaning |
|--------|---------|
| -O0    | Reduce compilation time and make debugging produce the expected results. |
| -g     | Produce debugging information in the operating system's native format. |
| -ffpe-trap= | Specifies a list of IEEE exceptions when a Floating Point Execption (FPE) should be raised. *invalid*: invalid floating point operation, such as SQRT(-1.0), *zero*: division by zero, *overflow*: overflow in a floating point operation. |
|-Wall | Enables commonly used warning options pertaining to usage that we recommend avoiding and that we believe are easy to avoid. This currently includes -Waliasing, -Wampersand, -Wsurprising, -Wintrinsics-std, -Wno-tabs, -Wintrinsic-shadow and -Wline-truncation. |
| -Wextra | Enables some warning options for usages of language features which may be problematic. This currently includes -Wcompare-reals and -Wunused-parameter. |
| -pedantic | Issue warnings for uses of extensions to Fortran 95. |
| -fcheck=all | Enable all run-time tests of -fcheck. |
| -fbacktrace | Specify that, when a runtime error is encountered or a deadly signal is emitted (segmentation fault, illegal instruction, bus error or floating-point exception), the Fortran runtime library should output a backtrace of the error. |

**Importrant**: Older versions of the GNU compiler (4.8.X and older) do not support (or have problems with) some functions from the C-interface that are used in PALM.


## Intel (ifx and ifort) for optimized code:

`ifx -fpp -fpe0 -O3 -g -traceback -xHost -fp-model source -ftz -no-prec-div`

| Option | meaning |
|--------|---------|
| -fpp   | Runs the Fortran preprocessor on source files before compilation. |
| -fpe0  | Specifies the floating-point exception handling level. Floating-point invalid, divide-by-zero, and overflow exceptions are enabled. If any such exceptions occur, execution is aborted. This option causes denormalized floating-point results to be set to zero. |
| -O3 | Performs O2 optimizations and enables more aggressive loop transformations such as Fusion, Block-Unroll-and-Jam, and collapsing IF statements. |
| -g  | The compiler produces complete debug information. |
| -traceback | Tells the compiler to generate extra informa tion in the object file to provide source  file traceback information when a severe error occurs at runtime. |
| -xHost | This option tells the compiler to generate instructions for the highest instruction set available on the compilation host processor. |
| -fp-model source | Controls the semantics of floating-point calculations. Rounds intermediate results to source-defined precision. **This option is required to force jobs to deliver reproduceble results. If not set, two PALM runs with exactly identical setup may produce different results.** Code performance is slightly reduced when using this option. |
| -ftz | This option flushes denormal results to zero when the application is in the gradual underflow mode. It may improve performance if the denormal values are not critical to your application's behavior. |
| -no-prec-div | Enables optimizations that give slightly less precise results than full IEEE division. |


## Intel (ifx and ifort) for debugging:

`ifort -fpp -fpe0 -O0 -check -check nooutput_conversion -traceback -g`

| Option | meaning |
|--------|---------|
| -O0    | Disables all optimizations. |
| -check | All available runtime checks (e.g. array bounds) are switched on. |
| -check nooutput_conversion | No checking for the fit of data items within a designated format descriptor field. **If this option isn't set, compilation aborts with errors.** |
| -traceback | Tells the compiler to generate extra information in the object file to provide source file traceback information when a severe error occurs at run time. |
| -g | Tells the compiler to generate full debugging information in the object file. |

## Cray (version 8.5.6):

`ftn -eZ -K trap=fp -em -O3 -hnoomp -hnoacc -hfp3 -hdynamic -dynamic`

| Option | meaning |
|--------|---------|
| -eZ    | Perform source preprocessing and compilation on Fortran source files. |
| -K trap=fp | Enable traps on divz, inv, or ovf exceptions. |
| -em    | When this option is enabled, the compiler creates .mod files to hold module information for future compiles. |
| -O3    | Switch on most aggressive optimization level. |
| -hnoomp | Disable compiler recognition of OpenMP directives. |
| -hnoacc | Disables the compiler recognition of OpenACC accelerator directives. |
| -hfp3   | Controls the level of floating point optimizations. (highest level would be -hfp4) |
| -hdynamic -dynamik | Directs the compiler driver to link dynamic libraries at runtime. |

## Cray (version 8.5.6) for debugging:

`ftn -eZ -K trap=fp -eD -ei -em -O0 -hnoomp -hnoacc -hdynamic -dynamic`

| Option | meaning |
|--------|---------|
| -K trap=fp | Enable traps on divz, inv, or ovf exceptions. |
| -eD    | Enables all debugging options. |
| -ei    | Initializes all undefined local stack, static, and heap variables of type REAL or COMPLEX to an invalid value (signaling NaN). |
| -O0    | Disables all optimizations including floating point optimizations and OpenACC acceleration. |

**Important on Cray-systems:** Getting tracebacks requires setting of `export ATP_ENABLED=1` before execution of PALM. Use an input command in the configuration file for that:
`IC:export ATP_ENABLED=1`

## NVidia

Options are for compiling and running on NVidia GPU's using `openacc`-directives.

`nvfortran -cpp -O3 -acc=verystrict -cuda -gpu=cc80 -Minfo=accel`

| Option | meaning |
|--------|---------|
| -cpp   | Runs the Fortran preprocessor on source files before compilation. |
| -O3 | All -O1 and -O2 optimizations are performed. In addition, this level enables more aggressive code hoisting and scalar replacement optimizations that may or may not be profitable. |
| -acc=<br>verystrict | Enable OpenACC pragmas and directives to explicitly parallelize regions of code for execution by accelerator devices. Abort compilation with a fatal error when an accelerator directive is encountered which does not adhere to the OpenACC standard. |
| -cuda | Enable CUDA Fortran; add CUDA include paths, and link with the CUDA runtime libraries.  -cuda is required on the link line. If -cuda is used in compilation, it must also be used for linking. |
| -gpu=<br>cc80 | Generate code for a device with compute capability 80. Multiple compute capabilities can be specified, and one version will be generated for each. By default, the compiler will detect the compute capability for each installed GPU. |
| -Minfo=<br>accel | Emit information about accelerator region targeting. |
