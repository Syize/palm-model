---
title: Traffic Module Reference
---
# Traffic Module Reference

---

## General Approach

The Traffic module calculates additional wind and potential temperature tendences which are caused by the moving vehicles. The module does not simulate in detail resolution the airflow around the individual car, instead it uses the calculation based on the drag force and drag coefficient.

* The wind tendencies in the volume of the car are augmented by the drag forces caused by the difference between the car and air speed in the particular car sides (front, side, top).
* The potential tendency are augmented by the tendency caused by the heat produced by burning of the fuel in the car.
* Appliaction of this wind and temperature forcing can be done in two different ways: 
    * This wind and temperature forcing can be applied in time steps when car is running over the grid cell. Length of this time interval is given by the car speed and length. The breaks between forcing is given by the transportation intensity.
    * The application of the forcing can be done continuously in all timesteps with the proper normalization of the tendency. This approach decreases induced divergences in the pressure field which need to be eliminated in the pressure solve. 
    * Our tests showed that differences between time averaged modelled output values of wind speeds and pollutant concentrations in these two approaches are small.
* The module needs following inputs for eash grid cell where traffic effects apply:
    * car types and characteristics (car dimensions, aerodynamic drag coefficient, ...)
    * transportation line characteristics (direction vector, slope, width)
    * transportation intensity in each car category (cars per hour)
    * transportation flow speed
    * heat produced by the running car (in each category, grid cell, and time).
* Prescribed transportation intensities, speed, and produced heat are variable in space and time.


## Calculation of the wind velocity tendency induced by cars

The drag force $F_{d}$ of the car is given by the formula:

$$
\begin{equation}
F_{d} = 0.5 \cdot C_{d} \cdot S_{car} \cdot \rho_{air} \cdot (s_{car} - s_{wind})^2
\end{equation}
$$

where

* $C_{d}$ $\left[ 1\right]$ is the aerodynamic drag coefficient in given car dimension
* $S_{car}$ $\left[m^2\right]$ is the area of the vehicle in the given car dimension (front area, side cross-section, footprint)
* $\rho_{air}$ $\left[kg \, m^{-3}\right]$ is air density
* $s_{car}$ $\left[m \, s^{-1}\right]$ is vehicle velocity in the given car dimension
* $s_{wind}$ $\left[m \, s^{-1}\right]$ is wind velocity in the given car dimension

Note, that the calculation is done in the car coordinate system ($u$, $v$, $w$) = (front, side, top). 
The area $S_{car}$ means the area of the car projection in the given coordinate and the influence of the specific car shape to the drag force is included in the drag coeficient $F_{d}$. 
For simplicity reasons, we consider the car a cuboid with dimensions length, width, and heigth of the car in the following normalization calculation. 

The reactive force of the drag force applies to the air mass in the volume of the car in every car coordinate. The acceleration (=tendency) of the air can be calculated from the reactive drag force by:

$$
\begin{equation}
tend_{car} = \frac{F_{d}}{m_{air}} = \frac{F_{d}}{\rho_{air} \cdot S_{car} \cdot L_{car}} = \frac{0.5 \cdot C_{d} \cdot (s_{air} - s_{wind})^2}{L_{car}}
\end{equation}
$$

where 

* $tend_{car}$ $\left[ m \, s^{-2}\right]$ is the additional wind tendency acused by the moving car. 
* $L_{car}$ $\left[ m\right]$ is car length in a given car dimension (length, width, height).
* $m_{air}$ $\left[ kg\right]$ is the mass of the air contained in the volume of the car.

This calculation is again done in the car coordinate system and then rotated back to the model coordinate system.

The car forcing is applied for the time when the car occupies given area. Corresponding time calculation is done as 

$$
\begin{equation}
t_{car} = \frac{L_{car}}{s_{car}}
\end{equation}
$$

where

* $t_{car}$ $\left[s\right]$ is time of the one car tendency impulse.
* $L_{car}$ $\left[m\right]$ is the length of the car in the front car dimension. 
* $s_{car}$ $\left[m \, s^{-1}\right]$ is speed of the car in the front car dimension.

Note, that this coefficient is comon for calculation of the forcing in all car dimensions.  

Another approach is to normalize the tendency to all simulation time. This approach is less realistic, 
but it allows lower continuous tendencies what decreases the disbalances of the divergences of the pressure field
which have to be removed in the pressure solver. Performed tests suggestes small differences between these two 
approaches in resulting averaged velocity and air quality output fields.
In this second case, the normalization coefficient $time\_coef$ of the forcing is calculated as

$$
\begin{equation}
coef_{time} = \frac{L_{car}}{s_{car} \cdot interval_{car}}
\end{equation}
$$

where

* $interval_{car}$ $\left[ s\right]$ is the time interval between cars is calculated from the transportation intensity in the street.

This module parameterizes effects of the transportation to the wind flow in the street canyon, it is not intended to simulate in details flow around an individual vehicle. 
We consider cars moving randomly inside the transportation line and the forcing is thus normalized and applied in all transportation line width. 
The current implementation applies uniform distribution of the forcing induced by cars inside the transportation line. 
In reality, the cars are moving more often near the centre of the line and the corresponding distribution could be implemented instead, 
but the practical impact of this distribution to the results is supposed to be minimal. 

For grids on the borders of the transportation line, the fractional approach is applied. 
For computational efficiency reasons, the fraction coefficients are not calculated as fully 3D intersection between transportation line and grid cell,
but are split to horizontal and vertical intersection. 
This allows to precalculate the horizontal fraction coefficients in the preprocessing (e.g. in GIS tool) and supply the coefficients as the grid cell parameters. 
The vertical intersection is done inside the module based on grid and car heights. 
Considering typical street slopes, inaccuracies caused by this simplification are negligible.
Fraction can also be used for definition of the arbitrary distribution of the effects of the traffic inside the transportation line.
This can be derived e.g. from distribution of the car movements inside the line width if such information is available.

Similar approach is used for transformation of the wind speed from model coordinates to car coordinates and back transformation of the forcing from car coordinate system to model coordinate system.
Both transformations are done in 2D in horizontal level and the forcing in the vertical w car coordinate is applied to the z model direction. 
Acording testing simulations, the tendencies induced in side and top directions are two magnitudes lower than the front tendencies. 


## Calculation of the theta tendency produced by car


The change of the air temperature is related to the imposed heat by the formula. 

$$
\begin{equation}
\Delta T_{air} = \frac{heat_{car}}{m_{air} \cdot c_{p}} 
\end{equation}
$$

where

* $\Delta T_{air}$ $\left[K\right]$ is the change of the air temperature caused by the imposed car heat
* $heat_{car}$ $\left[J\right]$ is the heat produced by the car
* $m_{air}$ $\left[kg\right]$ is the mass of the affected air
* $c_{p}$ $\left[J \; kg^{-1} K^{-1}\right]$ is the air specific heat


The tendency of the air potential temperature imposed by car can be calculated:

$$
\begin{equation}
tend_{car} = \frac{power_{car}}{\rho_{air} \cdot c_{p} \cdot exner}
\end{equation}
$$

where

* $tend_{car}$ $\left[ K \; s^{-1} \right]$ is the tendency of the air potential temperature caused by car
* $power_{car}$ $\left[ W \right]$ is the heat produced by car per second
* $\rho_{air}$ $\left[kg \, m^{-3}\right]$ is the air density
* $c_p$ $\left[J \, kg^{-1} \, K^{-1}\right]$ is the heat capacity of air

This tendency applies to the air for the time when the car occupies given grid cell or in all time steps, the corresponding temporal normalization
as well as the normalization to the transportation line and fraction calculation for the grid is done similarly as in case of velocity tendency.

The heat produced by car per second is included in the traffic driver for each grid, time, and each car category.
This allows to apply values obtained by different ways (simulations by transportation model, observations,...).
A simple way to obtain these values is to derive them from the fuel consumed by car and from the combustion heat of the fuel:

$$
\begin{equation}
power_{car} = \frac{cons_{car} \cdot \rho_{fuel} \cdot 1000 \cdot comb\_heat_{fuel}}{100000 \cdot s_{car}}
\end{equation}
$$

where

* $cons_{car}$ $\left[l/100km\right]$ is consumption of the car of given category per 100 km in given conditions
* $s_{car}$ $\left[m \, s-1\right]$ is the car speed
* $\rho_{fuel}$ $\left[kg \, m^{-3}\right]$ is the density of the fuel
* $comb\_heat_{fuel}$ $\left[J \, kg^{-3}\right]$ is the combustion heat of the fuel





