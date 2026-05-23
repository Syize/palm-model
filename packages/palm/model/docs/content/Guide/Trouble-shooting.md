---
title: Trouble-shooting
---

# Trouble-shooting
<br>

This document helps to localize problems that may appear during the execution of PALM and to determine their reasons. The first section gives respective instructions. In case that PALM bugs have been identified, trouble-tickets can be submitted to the PALM developer group. The second section describes the rules for submitting a trouble-ticket.

## How to localize problems and determine reasons

- Check for missing or wrong options in the `palmrun` command.
- Check for errors in the [configuration file](../../../Reference/LES_Model/Scripts/config_file).
- Check for errors in the namelist parameter file. For example, always check if namelists are correctly closed with a `/`, every input is closed with a comma, namelists are correctly named, no comment lines exist in list-shaped input variales, etc.
- Carfully analyse the job-protocol for error messages or unexpected behavior.
- Switch-on namelist parameters [debug_output](../../../Reference/LES_Model/Namelists/#runtime_parameters--debug_output) and [debug_output_timestep](../../../Reference/LES_Model/Namelists/#runtime_parameters--debug_output_timestep) and check the debug output. This will be written to debug files located in the temporary working directory and are named `DEBUG_xxx` with one file for each computing core. This requires to set the `palmrun`-option `-B`, otherwise the temporary working directory will be deleted at the end of the run.
- In case of unexpected program terminations, re-compile and re-run the job with code with debug-, traceback-, and floating point error detection options as outlined [here](Recommended_compiler_options.md).
- Switch-off any user-code.


## How to submit a trouble-ticket

Trouble-tickets are answered on a voluntary basis by the PALM developers. They are only for reporting bugs that appear during the execution of PALM. Reports can be submitted at [https://www.palm-model.org/report-a-bug](https://www.palm-model.org/report-a-bug).

**Reports will only be accepted and processed if they satisfy the following requirements:**

- the most recent PALM release has to be used
- the software environment installed on the used computer / cluster-system must fulfil the [software requirements](../../../Get_Started#software-requirements)
- debug-, traceback-, and floating point error detection options have to be switched on: [see here](Recommended_compiler_options.md)
- the problem must be described precisely and complete error-reports (e.g. logfiles of batch jobs or the complete terminal output of the `palmrun` command) as well as the run-control file (suffix `_rc` in the `MONITORING` folder) must be attached to the ticket
- the complete setup (INPUT files like namelist files, static- and dynamic drivers, etc.) to reproduce the problem must be attached
- the complete `palmrun` command that has been used to start the simulation must be given
- reduce large setups so that they can run on a smaller number of cores but still can reproduce the problem
- for problems appearing in batch jobs, the jobs need to be generated via `palmrun` option `-b`
- a user-interface is not allowed to be switched on

**Please be aware that the trouble-tickets are not intended to provide tutorials, nor do they offer assistance in developing physically correct and functioning setups or provide any help with the verification of setups.** For these purposes, please try our [reddit channel](https://www.reddit.com/r/PALM_model_system/).

Since the online documentation of PALM is currently revised and extended, some information may still be missing. In such a case, respective questions can be asked via the trouble-ticket system.

