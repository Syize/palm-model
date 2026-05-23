# Terrain height

Surface height above sea level

---

The terrain height can be supplied as a raster file with the [`zt`](yaml.md#zt) parameter or as a column in a vector polygon file with the [`zt`](yaml.md#zt) parameter. If the terrain height is missing everywhere, a default value of 0 m is used. Partially missing values in a domain are not allowed.

![Terrain height illustration](Figures/terrain_height_Berlin.svg)  
*Terrain height raster with the terrain height in meters above sea level.*

palm_csd calculates the minimum terrain height of all domains, subtracts this value from all terrain heights and stores the result in the global `origin_z` attribute of the resulting static driver.

Setting [`interpolate_terrain`](yaml.md#interpolate_terrain) to `true` will enable bilinear interpolation of the terrain height of the parent domain to the current domain. This reduces the height differences due to different grid spacings between parent and child. Using this option requires setting [`use_palm_z_axis`](yaml.md#use_palm_z_axis) set to `True` to raster the input data on the z-grid of PALM for output. Note that PALM will convert continuous static driver data itself on its grid and apply additional filtering procedures. It is thus recommended to set this parameter to `False` in the non-nested case.

A typical set-up could look like this:

```yaml
input:
  files:
    zt: terrain_height.tif

domain_root:
  use_palm_z_axis: False
  interpolate_terrain: False

domain_N02:
  use_palm_z_axis: True
  interpolate_terrain: True
```
