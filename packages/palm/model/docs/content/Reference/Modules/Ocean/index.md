---
title: Ocean Module Reference
---
# Ocean Module Reference

---

## General Features

PALM allows for studying the ocean mixed layer (OML) via the ocean mode, where the sea surface is defined at the top of the model with height `z` = *0.0*, and negative values of `z` indicate the depth. Hereafter, the term *surface* and the index *0* is used for variables at the sea surface (which is at the top of the ocean model). Fig. 1 shows the surface position in a vertical cross section of the PALM grid.

![ocean_grid](../../Figures/ocean_grid.png){width=50%} <br>
**Figure 1:** Vertical cross section of the PALM grid indicating the position of the ocean surface (*sea sfc*).

The ocean mode differs from the atmospheric mode by a few modifications, treated in the code by case distinctions, so that both versions share the same code base. In particular, seawater buoyancy and static stability in the ocean mode not only depend on $\theta$, but also on the salinity $\text{Sa}$. In order to account for the effect of salinity on density, a prognostic equation is added for $\text{Sa}$ (with unit PSU, *practical salinity unit*):

$\quad \begin{align*}
 \frac{\partial {\text{Sa}}}{\partial t} = -  \frac{\partial u_j
    {\text{Sa}}}{\partial x_j} - \frac{\partial}{\partial
    x_j}\left(\overline{u_j^{\prime\prime}{\text{Sa}}^{\prime\prime}}\right) +
  \Psi_{\text{Sa}}\;,
\end{align*}$

where $\Psi_{\text{Sa}}$ represents sources and sinks of salinity. Furthermore, the virtual potential temperature $\theta_\mathrm{v}$ is replaced by potential density $\rho_{\theta}$ in the buoyancy term of the momentum equation (see Sect. [wiki:doc/tec/gov governing equations])

$\quad \begin{align*}
 + g \frac{\theta_\mathrm{v} -
    \langle\theta_{\mathrm{v}}\rangle}{\langle\theta_{\mathrm{v}}\rangle}
  \delta_{i3}\,\rightarrow\,- g \frac{\rho_{\theta} -
    \langle\rho_{\theta}\rangle}{\langle\rho_{\theta}\rangle}
  \delta_{i3}\;,
\end{align*}$

in the stability-related term of the SGS-TKE equation (see Sect. [wiki:doc/tec/sgs turbulence closure])

$\quad \begin{align*}
+ \frac{g}{\theta_{\mathrm{v},0}}\overline{u_3^{\prime\prime}
    {\theta_{\mathrm{v}}}^{\prime\prime}}\,\rightarrow\,+
  \frac{g}{\rho_{\theta,0}}\overline{u_3^{\prime\prime}
    {\rho_{\theta}}^{\prime\prime}}\;,
\end{align*}$

as well as in the calculation of the mixing length (see Sect. [wiki:doc/tec/sgs turbulence closure])

$\quad \begin{align*}
  & \left(\frac{g}{\theta_{\mathrm{v},0}}
    \frac{\partial{\theta_{\mathrm{v}}}}{\partial
      z}\right)^{-\frac{1}{2}}\,\rightarrow\,\left(\frac{g}{\rho_{\theta,0}}
    \frac{\partial{\rho_{\theta}}}{\partial
      z}\right)^{-\frac{1}{2}}\;.
\end{align*}$

$\rho_{\theta}$ is calculated from the nonlinear equation of state of seawater after each time step using the algorithm proposed by Jackett et al. (2006). The algorithm is based on polynomials depending on $\text{Sa}$, $\theta$, and hydrostatic pressure $p_h$ (see Jackett et al. (2006), Table A2). In PALM, for  $p_h$ only the initial values enter this equation.

The ocean is driven by prescribed fluxes of momentum, heat and salinity at the top. The boundary conditions at the bottom of the model can be chosen as for atmospheric runs, including the possibility to use topography at the sea bottom.

The ocean mode can also account for the effect of surface waves (Langmuir circulation / Stokes drift and wave-breaking) as described in Noh et al. (2004).  Note that many older studies of the OML with PALM used the atmospheric mode, a subsequent inversion of the `z`-axis, and appropriate normalization of the results, instead of using the ocean mode (e.g., Noh et al., 2004, 2009).

The complete momentum equation including Stokes drift and wave breaking reads

$\quad \begin{align*}
\frac{\partial u_i}{\partial t} = &-  ( u_j + u_{s,j}) \frac{\partial u_i}{\partial x_j}
   - \epsilon_{ijk} f_j (u_k + u_{s,k}) + \epsilon_{i3j} f_3 u_{g,j} - \frac{\partial \pi^{*}}{\partial x_i} + \epsilon_{ijk} u_{s,j} \omega_k \\
   & - g \frac{\rho_{\theta} - \langle\rho_{\theta}\rangle}{\langle\rho_{\theta}\rangle} \delta_{i3} - \frac{\partial}{\partial x_j} \left( \overline{u_i^{\prime\prime} u_j^{\prime\prime}} - \frac{2}{3} e \delta_{ij} \right) + F_i \; ,
\end{align*}$

where $F_i$ is a random forcing representing the generation of small resolved-scale turbulence by wave breaking, $u_s$ the Stokes drift velocity

$\quad
 u_s = U_s \exp \left( \frac{4 \pi z}{\lambda_w} \right)$
 
with

$\quad
 U_s = \left( \frac{\pi z_w}{\lambda} \right)^2 \left( \frac{g \lambda_w}{2 \pi} \right)^{1/2}$

with wave height $z_w$ and wave length $\lambda_w$, and

$\quad
 \omega_i = \epsilon_{ijk} \frac{\partial u_k}{\partial x_j}$

the curl of the velocity.

The ocean mode features are also described in Maronga et al. (2019).


## Atmosphere-Ocean Coupling

This coupling has been developed in order to allow for studying interactions between
the atmospheric boundary layer and the ocean mixed layer. The coupling is realized via flux coupling by the exchange of respective fluxes at the sea surface (boundary conditions) between the atmosphere model and the ocean model. The atmospheric model uses a constant flux layer and transfers the
kinematic surface fluxes of heat and moisture as well as the momentum fluxes to the ocean model. Flux conservation between the ocean and the atmosphere requires to adjust the kinematic fluxes for the ocean by the density ratio of air and water $\rho_0 / \rho_{\mathrm{l},0}$:

$\quad \begin{align*}
  &                    \overline{w^{\prime\prime}u^{\prime\prime}}_0\big\vert_{\text{ocean}} = \frac{\rho_0}{\rho_{\mathrm{l},0}} \overline{w^{\prime\prime}u^{\prime\prime}}_0\;,\nonumber\\
  &
  \overline{w^{\prime\prime}v^{\prime\prime}}_0\big\vert_{\text{ocean}}
  = \frac{\rho_0}{\rho_{\mathrm{l},0}}
  \overline{w^{\prime\prime}v^{\prime\prime}}_0\;.
\end{align*}$

Since evaporation leads to a cooling of the surface water, the kinematic flux of heat in the ocean depends on both the atmospheric kinematic surface fluxes of heat and moisture and is calculated by

$\quad \begin{align*}
  &
  \overline{w^{\prime\prime}\theta^{\prime\prime}}_0\big\vert_{\text{ocean}}
  = \frac{\rho_0}{\rho_{\mathrm{l},0}} \frac{c_p}{c_{p, \mathrm{l}}}
  \left(\overline{w^{\prime\prime}\theta^{\prime\prime}}_0 +
    \frac{L_\mathrm{V}}{c_p}
    \overline{w^{\prime\prime}q^{\prime\prime}}_0 \right)\;.
\end{align*}$

Here, $c_{p, \mathrm{l}}$ = 4218 Jkg^-1^K^-1^ is the specific heat of water at constant pressure. Since salt does not evaporate, evaporation of water also leads to an increase in salinity at the ocean surface. This process is modeled after Steinhorn (1991) by a negative (downward) salinity flux at the ocean surface:

$\quad \begin{align*}
  \overline{w^{\prime\prime}S^{\prime\prime}}_0\big\vert_{\text{ocean}} = - \frac{\rho_0}{\rho_{\mathrm{l},0}} \frac{S}{1000\, \mathrm{PSU} - S} \overline{w^{\prime\prime}q^{\prime\prime}}_0\;.
\end{align*}$

Sea surface values of potential temperature and the horizontal velocity components are transferred back as surface boundary conditions to the atmosphere:

$\quad \begin{align*}
  & \theta_0 = \theta_0\big\vert_{\text{ocean}}\;,\,u_0 =
  u_0\big\vert_{\text{ocean}}\;,\,v_0 = v_0\big\vert_{\text{ocean}}.
\end{align*}$

If the horizontal grid spacings of the coupled models do not match, a two-way bi-linear interpolation is used for exchanging the 2d-surface data. A sketch of the data exchange between atmosphere and ocean is given in Fig. 2.

![flux_coupling](../../Figures/ocean_flux_coupling.png){width=80%} <br>
**Figure 2:** Data flow at the ocean-atmosphere interface. Blue arrows indicate quantities transferred from the atmosphere to the ocean, red arrows indicate quantities transferred from the ocean to the atmosphere. $z_p$ is the top of the constant flux layer in the atmosphere, which is the first computational grid point above the surface. $z_0$ is the height of the roughness layer. Respective roughness layer data are stored at vertical grid index `k` = *0*.


### Time Synchronization
The time steps for the atmosphere model and the ocean model are determined individually and don't need to be identical. The coupling is executed at a user-prescribed frequency (see [dt_coupling](../../../../Reference/LES_Model/Namelists/#runtime_parameters--dt_coupling)).


## Namelist Parameters

For a list of all namelist parameters see [`&ocean_parameters`](../../../../Reference/LES_Model/Namelists/#ocean-parameters).

## References

- **Esau, I.,** 2014: Indirect air-sea interactions simulated with a coupled turbulence-resolving model. Ocean Dynam. 64: 689–705. [doi:10.1007/s10236-014-0712-y](http://dx.doi.org/10.1007/s10236-014-0712-y)

- **Jackett, D.R., McDougall, T.J., Feistel, R., Wright, D.G., and Griffies, S.M.,** 2006: Algortihms for density, potential temperature, conservative temperature, and the freezing temperature of seawater. J. Atmos. Ocean. Tech. 23: 1709–1728.

- **Maronga, et. al.**, 2020: Overview of the PALM model system 6.0. Geosci. Model Dev., 13, 1335-1372.
[doi:10.5194/gmd-13-1335-2020](http://dx.doi.org/10.5194/gmd-13-1335-2020)

- **Noh, Y., Min, H.S., and Raasch S.,** 2004: Large eddy simulation of the ocean mixed layer: the effects of wave breaking and Langmuir circulation. J. Phys. Oceanogr. 34: 720–735.

- **Noh, Y., Goh, G., Raasch, S., and Gryschka, M.,** 2009: Formation of a diurnal thermocline in the ocean mixed layer simulated by LES. J. Phys. Oceanogr. 39: 1244–1257.

- **Steinhorn, I.,** 1991: Salt flux and evaporation. J. Phys. Oceanogr. 21: 1681–1683.

