---
title: Overview
---
# Nesting Overview

---

!!! warning
    This site is  Work in Progress.

PALM offers two different kinds of nesting.

- [self-nesting](self_nesting.md)<br>
  One or more instances of PALM (so called children or **childs**) with higher spatial resolution are nested inside the total domain (called the **root** domain). Recursive nesting (nests within nests) is allowed. A typical application for self-nesting is if high spatial resolution is only required at specific areas of interest, e.g. within the urban canopy layer or in the entrainment zone at the top of a mixed layer, and if computational resources do not allow to use that high resolution for the total tomain.<br><br>
- [mesoscale nesting](mesoscale_nesting.md)<br>
  Boundary layer processes contain a wide range of scales, ranging from the mesoscale, e.g. urban heat island, land-sea breeze, low-level jets, etc. down to the microscale, e.g. effects of single trees, building and roof shapes, local emissions, etc.. To consider both large model domains and a small grid size would often require huge computational resources. The idea of [mesoscale nesting](mesoscale_nesting.md) is to consider the effects of mesoscale processes in PALM via lateral (and top) boundary conditions provided by larger-scale models, and refine the grid within the domain of interest.

[Mesoscale nesting](mesoscale_nesting.md) and [self-nesting](self_nesting.md) can be used simultaneously. In a self-nesting setup, mesocale-nesting can be used for the root domain (but only for the root domain!).
