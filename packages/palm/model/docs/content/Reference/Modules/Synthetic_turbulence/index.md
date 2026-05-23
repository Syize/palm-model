---
title: Synthetic Turbulence Generator Reference
---
# Synthetic Turbulence Generator Reference

---

## General Approach

In the following, the general approach of the employed filter method to generate synthetic turbulence is described.

To obtain turbulent flow components $u_{i,\mathrm{b}}$ on the lateral boundaries, spatially and temporally correlated disturbances $u^{\prime\prime}_i$
are imposed onto the preliminary velocity components $\overline{u}_{i,\mathrm{b}}$:

$$
\begin{equation}
$u_{i,\mathrm{b}} = \overline{u}_{i,\mathrm{b}} + a_{ij}\,u^{\prime\prime}_j\;,
\end{equation}
$$

with $i,j \in {1,2,3}$.
$a$ is the amplitude tensor that is calculated from the Reynolds stress tensor $r$. To consider cross-correlations between the velocity components, Lund et al. (1998) suggested a Cholesky decomposition to compute $a$ recursively by:

$$
\begin{eqnarray}
a_{11} &=& \sqrt{r_{11}} \\
a_{12} &=& 0 \\
a_{13} &=& 0 \\
a_{21} &=& r_{21} / a_{11} \\
a_{22} &=& \sqrt{r_{22} - a^2_{21}} \\
a_{23} &=& 0 \\
a_{31} &=& r_{31} / a_{11} \\
a_{32} &=& (r_{32} - a_{21}\,a_{31}) / a_{22} \\
a_{33} &=& \sqrt{r_{33} - a^2_{31} - a^2_{32}} \;.
\end{eqnarray}
$$

Depending on characteristic length scales $L$ and time scales $T$ of the flow, which are defined individually for each velocity component in each spatial direction, $u^{\prime\prime}_i$ computes as

$$
\begin{equation}
u^{\prime\prime}_{i}\left(t+\Delta t\right) = \Psi_{i}(t-\Delta t)\,\exp\left(-\frac{\pi\,\Delta t}{2 T_i}\right) + \Psi_{i}(t) \left[1-\exp\left(-\frac{\pi\,\Delta t}{T_i}\right)\right]^{0.5} \; ,
\end{equation}
$$

with $\Delta t$ being the actual LES time step and the two-dimensional spatially correlated disturbances

$$
\begin{equation}
$\Psi^{m,l}_i = \sum_{j=-N_i}^{N_i}\sum_{k=-N_i}^{N_i} b_j\,b_k\,\zeta^{m+j,l+k}_i \; .
\end{equation}
$$

The subscripts $m$ and $l$ indicate grid positions at the lateral boundary, $N_i = 2 L_i / \Delta x_i$, with $\Delta x_i$ being the grid spacing. $\zeta_i$ indicates a set of equally-distributed random velocities with zero mean and unit variance that are individually computed for each velocity component.
Finally, the spatial filter function computes as

$$
\begin{equation}
b_i = b^{*}_i \left( \sum_{k=-N_i}^{N_i} b^{*2}_k \right)^{-0.5}\,, %exp\left( -\frac{pi |k| \Delta x}{L}\right) \right)^{-0.5} \; ,
\end{equation}
$$

with

$$
\begin{equation}
b^{*}_i = \exp\left( -\frac{\pi |k| \Delta x_i}{L_i}\right) \; .
\end{equation}
$$

With this approach, the imposed $u^{\prime\prime}_i$ reflects the prescribed Reynolds stress as well as the spatial and temporal correlation according to $L_i$ and $T_i$, respectively.

From a mathematically point of view, the imposed fluctuations should have zero mean.
Due to a finite sample of random numbers and the finite number of discrete grid points, however, the fluctuations have mean values slightly different from zero in practice.
In order to overcome this, Kim et al. (2013) proposed a correction for the boundary normal flow component in order to maintain constant mass flux through the boundary.
In order to avoid that the perturbations imposed onto the non-normal components may have non-zero mean too, e.g. non-zero mean $w-$ and $v-$component at the western model boundary, we correct the imposed turbulent velocity components as

$$
\begin{equation}
u^{\prime\prime}_{i,\mathrm{corr}} = u^{\prime\prime}_i - \frac{1}{S}\int_{\partial S} u^{\prime\prime}_i dS \; ,
\end{equation}
$$

with $S$ being the surface area of the respective lateral boundary.

For non-stationary flows, such as in mesoscale nested setups, an inflow boundary can become an outflow boundary and vice versa.
Hence, the turbulence generator is applied at each lateral boundary simultaneously in this case, while at opposite boundaries (west and east, as well as north and south) we use the same $\Psi_i$ and thus the same set of random numbers (velocities).
By doing so, we save computational resources because the same set of $\Psi_i$ is already available on the west/east and north/south boundary according to our parallelization strategy. In case of idealized turbulent inflow conditions the synthetic turbulence generator is only applied at the left boundary. Further, perturbations are imposed at the end of each LES time step at the last Runge-Kutta substep right before the Poisson equation is solved to guarantee that the flow is free of divergence.

## Parametrization of Turbulence Profile

In order to create time- and height dependent synthetic turbulence, respective information about the Reynolds stresses, as well as turbulent length and time scales for the velocity components are required. For stationary flows these information can be deduced from observations or from cyclic precursor simulations (Xie and Castro, 2008) and can be stored in an ASCII-file.
However, for non-stationary flows with pronounced diurnal cycles and/or changing synoptic conditions running precursor simulations is practically not feasible.
Also, to take these information from the mesoscale-model output is not possible, since this detailed information is often neither available nor part of the operational output.
Hence, to allow for an adjustment of the synthetic inflow turbulence to changing atmospheric conditions, we parametrize the Reynolds stresses based on the time-dependent mesoscale inflow profiles. We follow the set of parametrizations presented by Rotach et al. (1996), which they employed in stochastic dispersion modelling. Please note, the following set of parametrizations refer to the stream- and spanwise components of the Reynolds stresses that are not necessarily parallel to the $x$- or $y$-axis, respectively. In order to emphasize this, we indicate stream- and spanwise components with a tilde in the following.

Rotach's et al. (1996) parameterizations are based on parametrizations by Brost et al. (1982) who derived them from observations in stratocumulus-topped marine boundary layers, which often differ in their vertical structure and turbulence production compared to boundary layers over land. However, since Rotach et al. (1996) have successfully validated the set of parametrizations against observations over land for a wide range of stability regimes, we are confident that the chosen set of parametrizations can be universally employed.
Based on the original formulation by Brost et al. (1982), the variance of the streamwise flow component $\widetilde{r}_{11}$ is parametrized following Rotach et al. (1996):

$$
\begin{equation}
\widetilde{r}_{11}(z) = u_\ast^2\,\left(\,0.35\,\left(\frac{-z_\mathrm{i}}{\kappa L_\mathrm{o}}\right)^{2/3} + \left( 5 - 4 \frac{z}{z_\mathrm{i}} \right) \,\right)\,, \quad \textrm{for} \quad z \le z_\mathrm{i} \; ,
\end{equation}
$$

who added a correction term to account for unstable near-surface stratification. Here, $u_\ast$ is the friction velocity, $\kappa = 0.4$ the von-Kármán constant, $L_\mathrm{o}$ the Obukhov length, $z_\mathrm{i}$ the boundary-layer depth and $z$ being the height above ground.
For neutral and stable situations the first term is simply ignored.
Similarly, we estimate the variance of the spanwise flow component by adding a correction term to the original formulation proposed by Brost et al. (1982):

$$
\begin{equation}
\widetilde{r}_{22}(z) = u_\ast^2\,\left(\,0.35\,\left(\frac{-z_\mathrm{i}}{\kappa L_\mathrm{o}}\right)^{2/3} + \left( 2 - \frac{z}{z_\mathrm{i}} \right) \,\right)\,, \quad \textrm{for} \quad z \le z_\mathrm{i} \; .
\end{equation}
$$

The profile of vertical velocity variance is taken from Gryning et al. (1987) as

$\widetilde{r}_{33}(z) = w_\ast^2\,\left(\,1.5\,\left(\frac{z}{z_\mathrm{i}}\right)^{2/3}\,\exp\left[-2\frac{z}{z_\mathrm{i}}\right] + \left( 1.7 - \frac{z}{z_\mathrm{i}} \right) \left(\frac{u_\ast}{w_\ast}\right)^2 \,\right)\,, \quad \textrm{for} \quad z \le z_\mathrm{i}\,,$

with $w_\ast$ being the convective velocity scale.
The vertical transport of horizontal streamwise and spanwise momentum is estimated by Brost et al. (1982) as

$$
\begin{equation}
\widetilde{r}_{31}(z) = -u_\ast^2\,\left( \frac{z}{z_\mathrm{i}} - 1 \right)\,, \quad \mathrm{for} \quad z \le z_\mathrm{i} \, ,
\end{equation}
$$

and

$$
\begin{equation}
\widetilde{r}_{32}(z) = -u_\ast^2\,\left(0.4\frac{z}{z_\mathrm{i}}\left(1-\frac{z}{z_\mathrm{i}}\right)\right)\,, \quad \textrm{for} \quad z \le z_\mathrm{i} \; ,
\end{equation}
$$

respectively. To our knowledge, no comparable formulation to estimate $\widetilde{r}_{21}$
exists in the literature. Hence, we decided to simply set

$$
\begin{equation}
\widetilde{r}_{21} = \sqrt{ \widetilde{r}_{31}^2 + \widetilde{r}_{32}^2} \; ,
\end{equation}
$$

assuming isotropy of horizontal and vertical transport of horizontal momentum.
To estimate the boundary-layer depth for a wide range of stability regimes, including buoyancy- and purely shear-driven boundary layers, we calculated $z_\mathrm{i}$ from a bulk Richardson number criterion according to Heinze et al. (2017) based on the bulk Richardson number

$$
\begin{equation}
Ri_\mathrm{b}(z) = \frac{g}{\theta_\mathrm{v,s}}\,\frac{\theta_\mathrm{v}(z) - \theta_\mathrm{v,s}}{u_\mathrm{h}(z)^2} \cdot z \; .
\end{equation}
$$

Starting at the surface, $z_\mathrm{i}$ is defined as the height where $Ri_\mathrm{b}$ first exceeds the critical bulk Richardson number $Ri_\mathrm{b,c} = 0.25$, which revealed to be a robust criterion to estimate the depth of the layer with significant turbulent transports caused by the presence of the surface (Heinze et al., 2017).
Here, $u_\mathrm{h}$ denotes the horizontal wind speed from mesoscale model input, $\theta_\mathrm{v}$ the virtual potential temperature, $\theta_\mathrm{v,s}$ the virtual potential surface temperature inferred from the second prognostic level above the surface, following Heinze et al. (2017), and $g$ is the acceleration of gravity.
In case of one-dimensional input from [dynamic input file](../../IO-Files/Drivers/dynamic.md) (lod = 1), $z_\mathrm{i}$ is determined based on the mean profiles prescribed at the lateral boundaries, while in case two-dimensional input from the from [dynamic input file](../../IO-Files/Drivers/dynamic.md) (lod = 2, $xz-$ and $yz-$ slices of boundary data), $z_\mathrm{i}$ is determined locally at each $(x,y)$-boundary grid point and averaged horizontally afterwards.

For more detailed scientific discussion about the synthetic turbulence generator please refer to [Kadasch et al. (2021)](https://doi.org/10.5194/gmd-14-5435-2021).

In general, the synthetic turbulence generator is designed to only impose perturbations onto the three velocity components. But we have extended the turbulence generator by an option to impose also perturbations onto the potential temperature. Tests have shown that the adjustment fetch can be significantly reduced that way.

## Namelist Parameters

For a list of all namelist parameters see [`&stg_par`](../../../../Reference/LES_Model/Namelists/#synthetic-turbulence-generator-parameters).

## References

- **Brost, R. A., Wyngaard, J. C., Lenschow, D. H.** 1982: Marine Stratocumulus Layers. Part II: Turbulence Budgets, J. Atmos. Sci., 39, 818–836, [https://doi.org/10.1175/1520-0469](https://doi.org/10.1175/1520-0469(1982)039<0818:MSLPIT>2.0.CO;2,)

- **Gryning, S., Holtslag, A., Irwin, J., Sivertsen, B.** 1987: Applied dispersion modelling based on meteorological scaling parameters, Atmos. Environ., 21, 79–89, [doi.org/10.1016/0004-6981(87)90273-3](https://doi.org/10.1016/0004-6981(87)90273-3)

- **Heinze, R., Moseley, C., Böske, L. N., Muppa, S. K., Maurer, V., Raasch, S., Stevens, B.** 2017: Evaluation of large-eddy simulations forced with mesoscale model output for a multi-week period during a measurement campaign, Atmos. Chem. Phys., 17, 7083–7109, [doi.org/10.5194/acp-17-7083-2017](https://doi.org/10.5194/acp-17-7083-2017)

- **Kadasch, E., Sühring, M., Gronemeier, T., Raasch, S.** 2021: Mesoscale nesting interface of the PALM model system 6.0. Geoscientific Model Development, 14: 5435–5465,. [10.5194/gmd-14-5435-2021](https://doi.org/10.5194/gmd-14-5435-2021)

- **Kim, Y., Castro, I. P., Xie, Z.-T.** 2013: Divergence-free turbulence inflow conditions for large-eddy simulations with incompressible flow solvers, Comput. Fluids, 84, 56–68, [doi.org/10.1016/j.compfluid.2013.06.001](https://doi.org/10.1016/j.compfluid.2013.06.001)

- **Lund, T. S., Wu, X., Squires, K. D.** 1998: Generation of Turbulent Inflow Data for Spatially-Developing Boundary Layer Simulations, J. Comput. Phys., 140, 233–258, [doi.org/10.1006/jcph.1998.5882](https://doi.org/10.1006/jcph.1998.5882)

- **Rotach, M. W., Gryning, S.-E., Tassone, C.** 1996: A two-dimensional Lagrangian stochastic dispersion model for daytime conditions, Q. J. Roy. Meteor. Soc., 122, 367–389, [doi.org/10.1002/qj.49712253004](https://doi.org/10.1002/qj.49712253004)

- **Xie, Z. and Castro, I.** 2008: Efficient Generation of Inflow Conditions for Large Eddy Simulation of Street-Scale Flows, Flow Turbul. Combust., 81, 449–470, [doi.org/10.1007/s10494-008-9151-5](https://doi.org/10.1007/s10494-008-9151-5)
