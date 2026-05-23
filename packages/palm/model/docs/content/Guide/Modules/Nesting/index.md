---
title: Overview
---
# Nesting

---

!!! warning
    This site is  Work in Progress.

## First activation of nesting in a restart run

In order to save computational resources, childs can be activated after the root domain has advanced for a specific time interval.

First prepare the namelist files for root and child domains for an initial run as well as a restart run.

Omit the [nesting_parameters](/Reference/LES_Model/Namelists/#nesting-parameters) namelist in the `_p3d` file of the root domain only, and choose [end_time](/Reference/LES_Model/Namelists/#runtime_parameters--end_time) as the time when the nesting shall be activated. The `_p3dr` file of the root domain must contain the [nesting_parameters](/Reference/LES_Model/Namelists/#nesting-parameters) namelist. Start the initial run via command 
```bash
palmrun .... -a "d3# restart"
```
Set [end_time](/Reference/LES_Model/Namelists/#runtime_parameters--end_time) in the `_p3dr` file of the root domain and the `_p3d` files of all childs to the desired value. Start a restart run via command
```bash
palmrun .... -a "d3r activate_nesting"
```
The run will take the `_p3dr` file of the root domain, but the `_p3d` files of all childs.

Further restart runs may be carried out via command

```bash
palmrun .... -a "d3r"
```

Keep in mind to always use activation string `restart` to store restart data for the next restart run.



## Notes, shortcommings and open issues

1. So far, using a spin-up for child domains is not allowed (see [spinup_time](/Reference/LES_Model/Namelists/#initialization_parameters--spinup_time)).
