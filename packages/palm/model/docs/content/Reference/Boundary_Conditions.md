---
title: Boundary conditions
---
# Boundary conditions

---

!!! warning
    This site is  Work in Progress.
    
    
## Basics


PALM offers a variety of boundary conditions. Dirichlet or Neumann boundary conditions can be chosen for $u$, $v$, $θ$, $q_v$, and $p^*$ at the bottom and top of the model. For the horizontal velocity components the choice of Neumann (Dirichlet) boundary conditions yields free-slip (no-slip) conditions. Neumann boundary conditions are also used for the SGS-TKE. Kinematic fluxes of heat and moisture can be prescribed at the surface instead of Dirichlet conditions for temperature and humidity. At the top of the model, Dirichlet boundary conditions can be used with given values of the geostrophic wind. By default, the lowest grid level `k = 0` for the scalar quantities and horizontal velocity components is not staggered vertically and defined at the surface `z = 0`. Vertical velocity is assumed to be zero at the surface and top boundaries, which implies using Neumann conditions for pressure. 

Following Monin-Obukhov similarity theory (MOST) a constant flux layer can be assumed as boundary condition between the surface and the first grid level where scalars and horizontal velocities are defined ($k = 1$, $z_{\mathrm{MO}} = 0.5 Δz$). It is then required to provide the roughness lengths for momentum $z_0$ and heat $z_\mathrm{0,h}$. Momentum and heat fluxes as well as the horizontal velocity components are calculated using the following framework. The formulation is theoretically only valid for horizontally-averaged quantities. In PALM we assume that MOST can be also applied locally and we therefore calculate local fluxes, velocities, and scaling parameters. 
 
 
Following MOST, the vertical profile of the horizontal wind velocity 

$$
\begin{eqnarray}
u_\mathrm{h} = (u^2 + v^2)^{\frac{1}{2}} 
\end{eqnarray}
$$

is given in the surface layer by

$$
\begin{eqnarray}
& \frac{\partial u_\mathrm{h}}{\partial z} = \frac{u_\ast}{\kappa z}\Phi_\mathrm{m}\left(\frac{z}{L}\right)\;,
\end{eqnarray}
$$

where $\kappa = 0.4$ is the Von Kármán constant and $\Phi_\mathrm{m}$ is the similarity function for momentum in the formulation of Businger-Dyer (see e.g. [Panofsky and Dutton 1984](#panofsky)):

$$
\begin{eqnarray}
& \Phi_\mathrm{m} =
 \begin{cases}
   1 + 5 \frac{z}{L} & \text{for~}  \frac{z}{L} \geq 0 \\
    \left(1 - 16 \frac{z}{L}\right)^{-\frac{1}{4}} & \text{for~}
    \frac{z}{L} < 0\;.
\end{cases}
\end{eqnarray}
$$

Here, $L$ is the Obukhov length, calculated as 

$$
\begin{eqnarray}
L = \frac{\theta_\mathrm{v}(z) u_\ast^2}{\kappa g
    \left[\theta_\ast + 0.61 \theta(z) q_\ast + 0.61
      q_\mathrm{v}(z) \theta_\ast\right]}\;.
\end{eqnarray}
$$

The scaling parameters $\theta_*$ and $q_*$ are defined by MOST as

$$
\begin{eqnarray}
& \theta_\ast = -
  \frac{\langle{w^{\prime\prime}\theta^{\prime\prime}}_0\rangle}{u_\ast},~q_\ast
  = -
  \frac{\langle{w^{\prime\prime}q_\mathrm{v}^{\prime\prime}}_0\rangle}{u_\ast}\;
\end{eqnarray}
$$

with the friction velocity $u_*$ defined as 

$$
\begin{eqnarray}
& u_\ast =
  \left[\left(\langle{u^{\prime\prime} w^{\prime\prime}}_0\rangle\right)^2
    + \left(\langle{v^{\prime\prime} w^{\prime\prime}}_0\rangle\right)^2
  \right]^{\frac{1}{4}}\;.
\end{eqnarray}
$$

In PALM, $u_*$ is calculated from $u_\mathrm{h}$ at $z_\mathrm{MO}$ by vertically integrating from $z_0$ to $z_\mathrm{MO}$. From the equations above it is possible to derive a formulation for the horizontal wind components, viz. 

$$
\begin{eqnarray}
\frac{\partial u}{\partial z} =
  \frac{-\langle{u^{\prime\prime} w^{\prime\prime}}_0\rangle}{u_\ast \kappa
    z}
  \Phi_\mathrm{m}\left(\frac{z}{L}\right)\,\text{and~}\,\frac{\partial
    v}{\partial z} = \frac{-\langle{v^{\prime\prime}
      w^{\prime\prime}}_0\rangle}{u_\ast \kappa z}
  \Phi_\mathrm{m}\left(\frac{z}{L}\right)\;.
\end{eqnarray}
$$

Vertically integrating the above equation from $z_0$ to $z_\mathrm{MO}$ then yields the surface momentum fluxes $\langle{u^{\prime\prime} w^{\prime\prime}}_0\rangle,\;\; \langle{v^{\prime\prime} w^{\prime\prime}}_0\rangle$.

The formulations above all require knowledge of the scaling parameters $\theta_*$ and $q_*$. These are deduced by vertically integrating 

$$
\begin{eqnarray}
\frac{\partial \theta}{\partial z} =
  \frac{\theta_\ast}{\kappa z}
  \Phi_\mathrm{h}\left(\frac{z}{L}\right)~{\text{and}}~\frac{\partial
    q_\mathrm{v}}{\partial z} = \frac{q_\ast}{\kappa z}
  \Phi_\mathrm{h}\left(\frac{z}{L}\right)
\end{eqnarray}
$$

from $z_\mathrm{0,h}$ to $z_\mathrm{MO}$. The similarity function $\Phi_\mathrm{h}$ is given by 

$$
\begin{eqnarray}
\Phi_\mathrm{h} =
  \begin{cases}
    1 + 5 \frac{z}{L} & \text{for}  \frac{z}{L} \geq 0 \\
    \left(1 - 16 \frac{z}{L}\right)^{-1/2} & \text{for} \frac{z}{L}
    < 0\;.
  \end{cases}
\end{eqnarray}
$$

The Obukhov length can be calculated by solving an implicit equation relating $L$ to the bulk Richardson number. This is achieved using a Newton iteration method. The surface fluxes are calculated based on the following sequence of actions:

1. calculate $u_h$ 
2. determine Obukhov length 
3. calculate $u_*$
4. calculate $θ_*$ and $q_*$
5. derive surface fluxes 

Depending on whether Neumann (prescribed fluxes) or Dirichlet boundary conditions are used for temperature and humidity, the bulk Richardson number is related to the Obukhov length via 

$$
\begin{eqnarray}
Ri_\mathrm{b,Di} = \dfrac{z}{L} \cdot \dfrac{[\phi_\mathrm{H}]}{[\phi_\mathrm{M}]^2} \;\;\;\;\textnormal{(Dirichlet conditions)} \;,
\end{eqnarray}

\begin{eqnarray}
Ri_\mathrm{b,Ne} = \dfrac{z}{L} \cdot \dfrac{1}{[\phi_\mathrm{M}]^3} \;\;\;\;\textnormal{(Neumann conditions)} \;,
\end{eqnarray}
$$

where

$$
\begin{eqnarray}
[\phi_\mathrm{H}] = \log\left(\dfrac{z_\mathrm{MO}}{z_\mathrm{0,h}}\right) - \Phi_\mathrm{H}\left(\dfrac{z_\mathrm{MO}}{L}\right) + \Phi_\mathrm{H}\left(\dfrac{z_\mathrm{0,h}}{L}\right)\;,
\end{eqnarray}
$$

$$
\begin{eqnarray}
[\phi_\mathrm{M}] = \log\left(\dfrac{z_\mathrm{MO}}{z_\mathrm{0}}\right) - \Phi_\mathrm{M}\left(\dfrac{z_\mathrm{MO}}{L}\right) + \Phi_\mathrm{M}\left(\dfrac{z_\mathrm{0}}{L}\right)\;,
\end{eqnarray}
$$

are the integrated universal profile stability functions of $\Phi_\mathrm{m}$ and $\Phi_\mathrm{h}$ (see [Paulson 1970](#paulson), [Holtslag and Bruin 1988](#holtslag)), and 

$$
\begin{eqnarray}
Ri_\mathrm{b,Di} = \dfrac{g z_\mathrm{MO} \left(\theta_\mathrm{v,1} - \theta_\mathrm{v,0}\right)}{u_\mathrm{h}^2 \theta_v}\;,
\end{eqnarray}
$$

$$
\begin{eqnarray}
Ri_\mathrm{b,Ne} = - \dfrac{g z_\mathrm{MO} \langle{w^{\prime\prime} \theta_\mathrm{v}^{\prime\prime}}_0\rangle}{\kappa^2 u_\mathrm{h}^3 \theta_v}\;.
\end{eqnarray}
$$

The above equations are solved for $L$ by Newton iteration, i.e. finding the root of the equation 

$$
\begin{eqnarray}
f = Ri_\mathrm{b} - \dfrac{z}{L} \cdot \dfrac{[\phi_\mathrm{H}]^x}{[\phi_\mathrm{M}]^y} \;,
\end{eqnarray}
$$

where $x$ and $y$ depend on the chosen boundary conditions (see above). The solution is given by iteration of

$$
\begin{eqnarray}
L^{t+1} = L^t - \dfrac{f(L^t)}{f'(L^t)}
\end{eqnarray}
$$

with

$$
\begin{eqnarray}
f'(L) = \dfrac{df}{dL}
\end{eqnarray}
$$

until $L$ meets a convergence criterion.

The flat bottom of the model can be replaced by a Cartesian topography (see Sect. [Topography](#Topography)). 

By default, lateral boundary conditions are set to be cyclic in both directions. Alternatively, it is possible to opt for non-cyclic conditions in one direction, i.e., a laminar or turbulent inflow boundary (see Sect. [Laminar and turbulent inflow boundary conditions](#Laminar)) and an open outflow boundary on the opposite site (see Sect. [Open outflow boundary conditions](#outflow)). The boundary conditions for the other direction have to remain cyclic. A complete overview about the non-cyclic lateral boundary conditions is given in Sect. [non-cyclic lateral boundary conditions ](#Noncyclic).

In order to prevent gravity waves from being reflected at the top boundary, a sponge layer (Rayleigh damping) can be applied to all prognostic variables in the upper part of the model domain ( [Klemp and Lilly, 1978](#klemp1978)). Such a sponge layer should be applied only within the free atmosphere, where no turbulence is present. 
 
The model cn be initialized by horizontally homogeneous vertical profiles of potential temperature, water vapor mixing ratio (or a passive scalar), and the horizontal wind velocities. The latter can be also provided from a 1-D precursor run (see Sect. [1-D model for precursor runs](../Modules/1D-Model)). Furthermore, initialization with prescribed arbitrary profiles and with data from larger scale models is possible, too. Uniformly distributed random perturbations with a user-defined amplitude can be imposed to the fields of the horizontal velocities components to initiate turbulence.


## Laminar and turbulent inflow boundary conditions <a name="Laminar"></a>

In case of laminar inflow, Dirichlet boundary conditions are used for all quantities, except for the SGS-TKE e and perturbation pressure $π^∗$ for which Neumann boundary conditions are used. Data from vertical profiles, as taken for the initialization of the simulation, are used for the Dirichlet boundary conditions. In order to allow for a fast onset of turbulence, random perturbations can be imposed on the velocity fields within a certain area behind the inflow boundary (inlet). These perturbations may be imposed throughout the entire simulation. For the purpose of preventing gravity waves from being reflected at the inlet, a relaxation area can be defined after [Davies (1976)](#davies1976). So far, it was found to be sufficient to implement this method for temperature only. It is realized by an additional term in the prognostic equation for $\theta$ (see third equation in Sect. [governing equation](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/gov)): 

$$
\begin{eqnarray}
\frac{\partial \theta}{\partial t} = \ldots - C_{\text{relax}}
 \left(\theta - \theta_{\text{inlet}}\,\right)\;.
\end{eqnarray}
$$

Here, $\theta_\mathrm{inlet}$ is the stationary inflow profile of $\theta$, and $C_\mathrm{relax}$ is a relaxation coefficient, depending on the distance $d$ from the inlet, viz. 

$$
\begin{eqnarray}
C_{\text{relax}}(d) =
\begin{cases}
   F_{\text{inlet}} \cdot \sin^2 \left(\frac{\pi}{2} \frac{D - d}{D}          \right) & \text{for~}  d < D \;,\\
   0   & \text{for~}  d \ge D \;,\\
 \end{cases}
\end{eqnarray}
$$
 
with $D$ being the length of the relaxation region and $F_\mathrm{inlet}$ being a damping factor. For more information about the inflow boundary see Sect. [inflow boundary](#Inflow). 


## Turbulence recycling

If non-cyclic horizontal boundary conditions are used, PALM offers the possibility of generating time-dependent turbulent inflow data by using a turbulence recycling method. The method follows the one described by [Lund et al. (1998)](#lund1998), with the modifications introduced by [Kataoka and Mizuno (2002)](#kataoka2002). Figure 1 gives an overview of the recycling method used in PALM. 


![Turbulence_recycling](Figures/Turbulence_recycling.png)
Figure 1: Schematic figure of the turbulence recycling method used for generation of turbulent inflow. The configuration represents exemplary conditions with a built-up analysis area (brown surface) and an open water recycling area (blue surface). The blue arrow indicates the flow direction. 


The turbulent signal $\varphi'(y, z, t)$ is taken from a recycling plane which is located at a fixed distance $x_\mathrm{recycle}$ from the inlet: 
 
$$
\begin{eqnarray}
\varphi^{\prime}(y, z, t) = \varphi(x_{\text{recycle}},y, z, t) -
  \langle \varphi\rangle_y(z, t) \;,
\end{eqnarray}
$$

where $<\varphi>y(z, t)$ is the line average of a prognostic variable $\varphi \in {u, v, w, \theta, e}$ along $y$ at $x = x_\mathrm{recycle}$. $\varphi'(y, z, t)$ is then added to the mean inflow profile $<\phi_\mathrm{inflow}>_y(z)$ at $x_\mathrm{inlet}$ after each time step: 

$$
\begin{eqnarray}
\varphi_{\text{inlet}}(y, z, t) = \langle
  \varphi_{\text{inlet}}\rangle_y(z) + \Phi(z) \varphi^{\prime}(y, z, t) \;,
\end{eqnarray}
$$

with the inflow damping function $\Phi(z)$, which has a value of $1$ below the initial boundary layer height, and which is linearly damped to $0$ above, in order to inhibit growth of the boundary layer depth. $<\varphi_\mathrm{inlet}>_y(z)$ is constant in time and either calculated from the results of the precursor run or prescribed by the user. The distance $x_\mathrm{recycle}$ has to be chosen much larger than the integral length scale of the respective turbulent flow. Otherwise, the same turbulent structures may be recycled repeatedly, which would affect the turbulence spectrum. It is thus recommended to use a precursor run for generating the initial turbulence field of the main run. The precursor run can have a comparatively small domain along the horizontal directions. In that case the domain of the main run is filled by cyclic repetition of the precursor run data. Note that the turbulence recycling is only implemented for $u, v, w, \theta, q, e$, and passive scalar $s$. 

Turbulence recycling is frequently used for simulations with urban topography. In such a case, topography / building elements should be placed sufficiently downstream of $x_\mathrm{recycle}$ to prevent effects on the turbulence at the inlet. Experience showed that for a neutral case with low wind speed and a large city domain, the distance $x_\mathrm{recycle}$ should be in the order of $1000\,\text{m}$ (see, e.g., [Gronemeier et al. (2017)](#gronemeier2017)).


## Open outflow boundary conditions <a name="outflow"></a>
At the outflow boundary (outlet), the velocity components $u_i$ meet radiation boundary conditions, viz. 

$$
\begin{eqnarray}
\frac{\partial u_i}{\partial t} + U_{u_i} \frac{\partial u_i}{\partial n} = 0 \;,
\end{eqnarray}
$$

as proposed by [Orlanski (1976)](#orlanski1976). Here $\partial/\partial n$ is the derivative normal to the outlet, and $U_{u_i}$ is a transport velocity, which includes wave propagation and advection. Rewriting the equation above yields the transport velocity 

$$
\begin{eqnarray}
U_{u_i} = -\left(\frac{\partial u_i}{\partial t}\right)\left(\frac{\partial u_i}{\partial n}\right)^{-1}
\end{eqnarray}
$$

that is calculated at interior grid points next to the outlet at the preceding time step for each velocity component. The transport velocity must be in the range $0 ≤ U_{u_i} ≤ \mathbf{Δ}/Δt$. In PALM, the upper threshold value of $U_{u_i} = \mathbf{Δ}/Δt$ is always used at the entire outlet. The above equations are discretized using an upstream method following [Miller and Thorpe (1981)](#miller1981). As the radiation boundary condition does not ensure conservation of mass, a mass flux correction can be applied at the outlet (see [mass flux correction](#Mass)). For more information about the outflow boundary see Sect. [outflow boundary](#Outflow). 


## Topography <a name="Topography"></a>
The Cartesian topography in PALM is generally based on the mask method ([Briscolini and Santangelo, 1989](#briscolini1989)) and allows for explicitly resolving solid obstacles such as buildings and orography. The implementation makes use of the following simplifications: 

1. the obstacle shape is approximated by (an appropriate number of) full grid cells to fit the grid, i.e., a grid cell is either $100\,\%$ fluid or $100\,\%$ obstacle, 
2.     the obstacles are fixed (not moving).

Overhanging structures as for example bridges, ceilings, or tunnels, are allowed, i.e. topography does not necessarily be surface-mounted. If no overhanging structures are present, the 3-D obstacle dimension reduces to a 2.5-D topography format, which is conform to the Digital Elevation Model (DEM) format. In case of overhanging structures, however, 3-D topography information is required to mask obstacles and their faces in PALM. 

The model grid separates into three cases (see [Fig. 2](#fig2)): 

1. grid points in free fluid without adjacent surfaces, where the standard PALM code is executed, 
2. grid points next to surface that require extra code (e.g., surface parametrization, boundary conditions), and 
3. grid points within obstacles, where prognostic equations are executed but results are multiplied by zero. 


![mask_method](Figures/mask_method.png)
Figure 2 <a name="fig2"></a>: Sketch of the topography implementation using the mask method (here for w). The yellow line represent the region where the prognostic equations for scalars and the w-component are masked, while the red line indicates the region where special surface-bounded code is masked. Additional surface-bounded code is executed at grid points in between the yellow and the red line (see [topography implementation](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/topography)). 


Additional code is executed in grid volumes affected by surfaces. The faces of the obstacles are always located in a way that the staggered-grid wall-normal velocity component ($u$, $v$, or $w$) is defined at the surface. (cf. Fig. 1 in Sect. [discretization](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/discret)). This allows to consider the impermeability boundary condition by setting the respective surface-normal velocity component to zero.

In case of the 5th-order advection scheme, the numerical stencil at grid points adjacent to obstacles would require data located within the obstacle. To avoid this, the order of the advection scheme is successively degraded at respective grid volumes adjacent to obstacles, i.e., from the 5th-order to 3rd-order at the second grid point above/beside an obstacle and from 3rd-order to 1st-order at grid points directly adjacent to an obstacle.

Surfaces in PALM can be aligned upward facing (e.g. bottom surface or rooftop), downward facing (e.g. undersurface of bridges), or vertical (facing north, east, south or west direction). At horizontal surfaces, PALM allows to either specify the surface values ($\theta$, $q_v$, $s$) or to prescribe their respective surface fluxes. The latter is the only option for vertically oriented surfaces. Simulations with topography require the application of MOST between each surface and the first computational grid point outside of the topography. For vertical and downward-facing surfaces, neutral stratification is assumed for MOST, even if, strictly speaking, MOST is derived for upward-facing surfaces only. This is simply attributed to the lack of knowledge in the literature about the best practice in this matter. 

The technical realization of the topography and treatment of surface-bounded grid cells will be outlined in Sect. [topography implementation](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/topography). 


## Non-cyclic lateral boundary conditions <a name="Noncyclic"></a>
[Figure 4](#fig4) shows the grid structure for non-cyclic boundary conditions at the left/right boundary **LB/RB** ([bc_lr](../Namelists/#initialization_parameters--bc_lr)) and [figure 5](#fig5) for non-cyclic boundary conditions at the north/south boundary **NB/SB** ([bc_ns](../Namelists/#initialization_parameters--bc_ns). The indices ($i$,$j$,$k$) represent the directions ($x$,$y$,$z$). The model domain extends from $-1:nx+1$ in the $x$-direction, from $-1:ny+1$ in the $y$-direction and from $0:nzt+1$ in the $z$-direction. For the advection scheme of Wicker and Skamarock, two more grid points are added at the lateral boundaries which are not needed for non-cyclic boundary conditions. The figures display the grid layer of the horizontal velocity components $u$ and $v$, and scalar $s$. The grid points of the vertical velocity $w$ are defined at the scalar position but shifted by one half grid spacing in vertical direction (not shown, detailed information about the grid structure in PALM can be found [here](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/discret)). The prognostic equations are solved at all inner grid points which are marked black. The grid points at the respective non-cyclic boundaries (blue) are treated as follows. **LB** is defined at $i = -1$ for $v$, $w$, $s$ and at $i = 0$ for $u$. **SB** is defined at $j = -1$ for $u$, $w$, $s$ and at $j = 0$ for $v$. **RB** is defined at $i = nx + 1$ and **NB** at $j = ny + 1$ for all quantities. **LB** and **SB** are treated this way so that the order and number of grid points for the streamwise velocity component and scalars is the same, independent of the flow direction. 

For technical reasons, the prognostic equations are first solved for $u$ at $i = 0$ ($v$ at $j = 0$), since these grid points technically belong to the inner grid, but afterwards, these results at $i = 0$ ($j = 0$) are replaced by the respective boundary condition in routine `boundary_conds`. In case of a Dirichlet condition, the values at $i = 0$ ($j = 0$) are taken from $i = -1$ ($j = -1$). In case of a radiation boundary condition, the solution of the Sommerfeld equation overwrites the prognostic values at $i = 0$ ($j = 0$). 


For non-cyclic lateral boundary conditions along only one of the horizontal directions, the multigrid solver has to be used to solve the Poisson equation (see [psolver](../Namelists/#initialization_parameters--psolver)). 


![grid_lr](Figures/grid_lr.png)
Figure 4 <a name="fig4"></a>: Grid structure at the lateral boundaries with non-cyclic lateral boundary conditions along the left-right direction. 


![grid_ns](Figures/grid_ns.png)
Figure 5 <a name="fig5"></a>: Grid structure of the lateral boundaries with non-cyclic lateral boundary conditions along the north-south direction. 


### Inflow boundary <a name="Inflow"></a>
At the inflow boundary, Dirichlet conditions are used for the three velocity components $\phi = \{u,v,w\}$ as well as for all scalar quantities $s$, and are implemented as follows (here e.g. for $s$ and a flow in positive $x$-direction):

$$
\begin{eqnarray}
s^{t + \Delta t}(k,j,-1) = s_\mathrm{init}(k) \; .
\end{eqnarray}
$$

$t$ denotes the time, $\Delta t$ the time step and $s_\mathrm{init}$ the initialization profile of the scalar quantities which is constant in time. The quantities at the inflow are either determined by the initial vertical profiles (see [initializing_actions](../Namelists/#initialization_parameters--initializing_actions)) or by data given via the dynamic driver. A Neumann condition is used for the subgrid-scale turbulent kinetic energy $e$ (here e.g. for a left-right flow): 

$$
\begin{eqnarray}
e^{t + \Delta t}(k,j,-1) = e^{t + \Delta t}(k,j,0) \; . 
\end{eqnarray}
$$

To prevent gravity waves from being reflected at the inflow, a relaxation term can be added to the prognostic equations for the potential temperature $theta$ [(Davies, 1976)](#davies1976): 

$$
\begin{eqnarray}
\theta^{t+1}(d) = ... - \Delta t \cdot K(d) \cdot \left( \theta^{t}(d) - \theta_\mathrm{init} \right) \; . 
\end{eqnarray}
$$

Here, $d$ is the distance normal to the wall and $\theta_\mathrm{init}$ the initial value of the potential temperature, which corresponds to the value at the inflow boundary. The damping or relaxation function $K$ depends only on the distance $d$ to the inflow. $K$ is calculated by 

$$
\begin{eqnarray}
K(d) =
\begin{cases}
d_\mathrm{f} \sin^2\left( \frac{\pi}{2} \frac{d_\mathrm{w} - d}{d_\mathrm{w}} \right) , \text{for } d < d_\mathrm{w} \\  \qquad\quad  0  \qquad \quad \;\;\; , \text{for } d \ge d_\mathrm{w}
\end{cases}
\end{eqnarray}
$$

where $d_\mathrm{f}$ is a damping factor to control the damping intensity, and $d_\mathrm{w}$ is the width of the relaxation region extending from the inflow. Quantities $d_\mathrm{f}$ and $d_\mathrm{w}$ can be set via parameters [pt_damping_factor](../Namelists/#initialization_parameters--pt_damping_factor) and [pt_damping_width](../Namelists/#initialization_parameters--pt_damping_width), respectively. Both parameters have to be set and must be adjusted case-by-case, because they depend on the numerical and physical conditions, so that application of universal default values is not possible.

As an example, simulations with gravity waves in case of cold air outbreaks show that they can grow in amplitude up to quite extreme values, if no damping is applied. In respective simulations, typical values for [pt_damping_factor](../Namelists/#initialization_parameters--pt_damping_factor) of *0.05* and for [pt_damping_width](../Namelists/initialization_parameters--pt_damping_width) of *25* km have been used in order to prevent the gravity waves from growing. 



### Outflow boundary <a name="Outflow"></a>
At the outflow, an open boundary condition is needed to ensure that disturbances of the mean flow can exit the model domain without effecting the flow upstream. For scalars, this is realized via Neumann boundary conditions at the outflow boundary. For the outflow plane normal velocity component a radiation boundary condition can be used (see [bc_lr](../Namelists/#initialization_parameters--bc_lr) and [bc_ns](../Namelists/#initialization_parameters--bc_ns)), for which the Sommerfeld radiation equation is solved at the outflow: 

$$
\begin{eqnarray}
\partial_t \psi  + c_{\psi} \partial_n \psi  = 0 \; .
\end{eqnarray}
$$

It considers flow disturbances propagating with the mean flow and by waves. Here $\psi$ is the transported quantity and $\partial_\mathrm{n}$ is the derivative normal to the outflow boundary. In general, following [Orlanski (1976)](#orlanski1976), the phase velocity $c_\psi$ should be calculated as 

$$
\begin{eqnarray}
c_{\psi} = - \dfrac{\partial_t \psi}{\partial_n \psi}  \; . 
\end{eqnarray}
$$

The discretized equation for the phase velocity (for outflow at the right boundary ($i = nx + 1$, Left-right flow, see [Fig. 4](#fig4)) reads

$$
\begin{eqnarray}
c_{\psi} = - c_\mathrm{max} \dfrac{\psi^t_{nx} - \psi^{t - \Delta t}_{nx} }{ \psi^{t-\Delta t}_{nx} - \psi^{t - \Delta t}_{nx-1} }   \; .
\end{eqnarray}
$$

The maximum phase velocity must not exceed

$$
\quad \begin{eqnarray}
c_\mathrm{max} = \dfrac{\Delta x}{\Delta t} \; . 
\end{eqnarray}
$$

Instead of calculating the phase velocity via the above equations, PALM always uses a constant phase velocity, which is assumed as the maximum velocity allowed by the CFL criterion, i.e. for a Courant number of one. Setting $c_\psi = c_\mathrm{max}$ leads to a simplified radiation boundary condition 
(here e.g. for a left-right flow along positive $x$-direction): 

$$\begin{eqnarray}
\psi^{t + \Delta t}(k,j,nx+1) = \psi^{t}(k,j,nx) \; ,
\end{eqnarray}
$$

with $\psi = \{u,v,w\}$. Although [Orlanski (1976)](#orlanski1976) suggested that this approach leads to reflection for waves smaller than $c_\mathrm{max}$ which may occur in complex geophysical flows, simulations of stable and convective boundary layers with background wind have shown no problems. Furthermore, this formulation of the radiation boundary condition saves computational time.


### Mass flux correction <a name="Mass"></a>
PALM offers the possibility of a mass flux correction at the outflow (e.g. [Tian, 2004](#tian2004)). Via setting parameter [conserve_volume_flow](../Namelists/#initialization_parameters--conserve_volume_flow) = *.T.*, the mass flux at the inflow and outflow is calculated by: 

$$
\begin{eqnarray}
\dot{m} = \sum_{k=1}^{nz-1} \Delta z(k) \sum_{l=0}^{nx_i} \psi(l,k) \Delta x_i \; , 
\end{eqnarray}
$$

where $\Delta x_i$ and $\psi$ is equal to $\Delta y (\Delta x$) and $u (v)$ in case of [bc_lr](../Namelists/#initialization_parameters--bc_lr) ([bc_ns](../Namelists/#initialization_parameters--bc_ns)). The correction factor for the outflow velocity, which eliminates the different mass flux at inflow and outflow, can be calculated by 

$$
\begin{eqnarray}
\psi_{corr} = \dfrac{\dot{m}_{inflow} - \dot{m}_{outflow}}{A} \; ,
\end{eqnarray}
$$

where $A$ is the area of the boundary plane as given by

$$
\begin{eqnarray}
A = \sum_{k=1}^{nz} \Delta z(k) \sum_0^{nx_i} \Delta x_i \; .
\end{eqnarray}
$$

The streamwise velocity at the outflow is corrected by adding $\psi_\mathrm{corr}$ at each grid point of the outflow between bottom and top ($k=1:nz$), in order to guarantee that the mass leaving the domain exactly balances the one entering it. 


## Synthetic Turbulence Generator
A synthetic turbulence generator is implemented to generate turbulence at the inflow plane (i.e. a turbulent signal is added to the velocity boundary values). The method is based on the work of [Xie and Castro (2008)](#xie2008) and [Kim et al. (2013)](#kim2013). Unscaled turbulent motions $u_{*j}$ are computed based on length scales along each direction and the amplitude tensor $a_{ij}$, which in turn is based on the Reynolds stress tensor. The calculated turbulence is then added to the horizontal mean inflow data of the velocity components $U_i$: 

$$
\begin{eqnarray}
u_i = U_i + a_{ij} u_{*j}.
\end{eqnarray}
$$

The amplitude tensor $a_{ij}$ depends on the Reynolds stress tensor $R_{ij}$ and is calculated using a Cholesky decomposition as suggested by [Lund et al. (1998)](#lund1998). The unscaled turbulent motions $u_{*j}$, which are calculated on the 2D inflow plane, depend on the prescribed time scales $t_{ij}$ and length scales $l_{ij}$: 

$$\begin{eqnarray}
u_{*i}(t+\Delta t) = u_{*i}(t) \exp\left(-\dfrac{C\Delta t}{T}\right) + \Psi_i(t,L)\left[1-\exp\left(-\dfrac{2C\Delta t}{T}\right)\right]^{0.5},
\end{eqnarray}
$$

where $\Psi$ denotes a part of the generated 2D signal which is correlated in space using the turbulent length scales $l_{ij}$ along the vertical and spanwise direction. Correlation along streamwise direction is assured via the time scale $t_{ij}$, which is estimated by $t_{ij}$ along streamwise direction using the $U_i$. 

After adding the turbulence to the mean inflow profiles, a mass flux correction suggested by [Kim et al. (2013)](#kim2013) is performed: 

$$
\begin{eqnarray}
u_{i,c} = \dfrac{U_{b,p}}{U_b} u_i \; ,
\end{eqnarray}
$$

where 

$$
\begin{eqnarray}
U_{b} = \dfrac{\int_S dS u_n}{S} \; ,
\end{eqnarray}
$$

and $u_{i,c}$ is the corrected wind velocity at the inflow boundary, $U_b$ and $U_{b,p}$ the instantaneous and prescribed bulk velocity at the inflow boundary, $S$ the surface area of the inflow boundary, and $u_n$ the component of $u_i$ normal to the inflow boundary.

The required length- and time scales, as well as the Reynolds strees tensor can be either prescribed (method 1), if known from previous simulations or measurements, or they can be parametrized (method 2). Please note, time and length scales as well as the components of the Reynolds stress tensor depend on the height and the horizontal location, particularly over heterogeneous surfaces and model domains with large horizontal extensions. For the sake of simplicity, only height dependent information of the Reynolds stress as well as length and time scales are considered. 

**Method 1**:<br>
If length- and time scales, as well as the Reynolds strees tensor is available, they can be provided via an ASCII file which contains all necessary information. This ASCII input file will be automatically read if an [stg_par](../Namelists/#synthetic-turbulence-generator-parameters) namelist is provided in the `_p3d` file.

Be sure that the required input file is added to the list of input files in the file connection file [.palm.iofiles](https://gitlab.palm-model.org/palm/model/-/blob/master/share/config/.palm.iofiles): 

```
STG_PROFILES    in:locopt    d3#:d3r    $base_data/$run_identifier/INPUT    _iprf
```

and named with the suffix `_iprf`. This file needs to be provided in the `INPUT` folder. Please have look at [STG_PROFILES](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles#STG_PROFILES) for a detailed format description of this input file. 

**Method 2**:<br>
In many cases detailed information about the Reynolds stress and turbulent length scales are not available, so that these information need to be parametrized. If no ASCII input file is provided in the input folder, this will be done automatically and the turbulence statistics at the inflow boundary will be estimated. Please note, the derived turbulence statistics will depend on the height above ground but not on the horizontal location. Parametrization of the Reynolds stress follows [Rotach et al. (1996)](#rotach1996). The diagonal components $R_{11}$, $R_{22}$, indicating the horizontal velocity variances, are estimated as follows: 

$$
\begin{eqnarray}
R_{1,1} = u_{*}^2\,( 0.35\,(-\frac{z_i}{\kappa\,L})^\frac{2}{3} + (5 - 4\,\frac{z}{z_i}) )\,, \: (i \in (1,2) \; ,
\end{eqnarray}
$$

with $u_*$ being the friction velocity, $\kappa$ the von-Kármán constant, $L$ the Obukhov length, and $z_i$ the mean boundary-layer depth. Please note, $u_*$, $L$ and $z_i$ are area-averaged values in this case. $z$ describes the height of the respective model grid level. $u_*$ is estimated from the mean horizontal wind speed at the first vertical grid point from the data provided at the lateral boundary using MOST. For the sake of simplicity, neutral conditions are assumed with $\Phi_m = 1$. $L$ is computed from the area-averaged surface temperature, surface sensible heat flux and roughness length in the model domain. $z_i$ is estimated from the bulk Richardson criterion, with $z_i$ being the height where the bulk Richardson first exceeds the critical Richardson number of $0.25$, according to [Heinze et al. (2017)](#heinze2017). In case of stable stratification ($L > 0$) or neutral stratification ($L = 0$), the first term is omitted in the computation of $R_{ii}$. 

Further, vertical velocity variances are parametrized as 

$$
\begin{eqnarray}
R_{3,3} = w_{m}^2\,( 1.5\,(\frac{z}{z_i})^\frac{2}{3}\,exp[-2\,\frac{z}{z_i}] + (1.7 - \frac{z}{z_i})\,(\frac{u_{*}}{w_m})^2)\,, \: (i \in (3) \; ,
\end{eqnarray}
$$

with

$$\begin{eqnarray}
w_m = (u_{*}^{3} + 0.6\,w_{*}^{3})^\frac{1}{3}
\end{eqnarray}
$$

being the momentum velocity scale, with the convective velocity scale $w_* $. In case of stable or neutral stratification, $w_m = u_* $. The remaining components $R_{31}$,$R_{32}$,$R_{21}$ are parametrized as 

$$
\begin{eqnarray}
R_{i,j} = u_{*}^2\,( 1 - exp[(3\,(\frac{z}{z_i} - 1)]) \; .
\end{eqnarray}
$$

In order to allow for a dynamic adjustment of the strength of the imposed turbulence, the Reynolds stress tensor is updated at regular intervals that can be adjusted individually. In other words, the turbulence generator considers for changing atmospheric stability in case of a multi-day simulation. For further information about this approach we refer to ​[Kadasch et al. (2021)](#kadasch2021).

At this point we emphasize that using the turbulence generator from [Xie and Castro (2008)](#xie2008) only generates turbulence which is correlated in space and time but not necessarily generate realistic turbulent structures. Large coherent structures like e.g. hexagonal patterns as typically observed in a convective boundary layer, cannot be generated by this method. Furthermore, please not that turbulence is only added to the three wind components. No perturbations are added to the subgrid-scale turbulent-kinetic energy and potential temperature.


## References
* **Briscolini M, and Santangelo P.** <a name="briscolini1989"></a> 1989. Development of the mask method for incompressible unsteady flows. J. Comput. Phys. 84: 57–75 [doi.org/10.1016/0021-9991(89)90181-2](https://doi.org/10.1016/0021-9991(89)90181-2).

* **Davies HC.** <a name="davies1976"></a> 1976. A lateral boundary formulation for multi-level prediction models. Q. J. Roy. Meteor. Soc. 102: 405–418 [doi.org/10.1002/qj.49710243210](https://doi.org/10.1002/qj.49710243210).

* **Gronemeier T, Raasch S, Ng E.** <a name="gronemeier2017"></a> 2017. Effects of Unstable Stratification on Ventilation in Hong Kong. Atmosphere. 8(9). 168 [doi.org/10.3390/atmos8090168](https://doi.org/10.3390/atmos8090168).

* **Heinze R, Moseley C, Böske L, Muppa S, Maurer V, Raasch S.** <a name="heinze2017"></a> 2017. Evaluation of large-eddy simulations forced with mesoscale model output for a multi-week period during a measurement campaign. Atmos. Chem. Phys. 17: 7083-7109 [doi.org/10.5194/acp-17-7083-2017](https://doi.org/10.5194/acp-17-7083-2017).

* **Holtslag AAM, Bruin HARD.** <a name="holtslag"></a> 1988. Applied modelling of the night-time surface energy balance over land. J. Appl. Meteorol. 27: 689–704 [doi.org/10.1175/1520-0450(1988)027<0689:AMOTNS>2.0.CO;2](https://doi.org/10.1175/1520-0450(1988)027<0689:AMOTNS>2.0.CO;2).

* **Kadasch E, Sühring M, Gronemeier T, Raasch S.**  <a name="kadasch2021"></a> 2021. Mesoscale nesting interface of the PALM model system 6.0. Geosci. Model Dev., 14: 5435–5465 [doi.org/10.5194/gmd-14-5435-2021](https://doi.org/10.5194/gmd-14-5435-2021).

* **Kataoka H, Mizuno M.**  <a name="kataoka2002"></a> 2002. Numerical flow computation around aerolastic 3d square
  cylinder using inflow turbulence. Wind Struct. 5: 379–392 [doi.org/10.12989/was.2002.5.2_3_4.379](https://doi.org/10.12989/was.2002.5.2_3_4.379).
 
* **Kim Y, Castro IP, Xie ZT.** <a name="kim2013"></a> 2013. Divergence-free turbulence inflow conditions for large-eddy simulations with incompressible flow solvers. Comput. Fluids. 84: 56–68 [doi.org/10.1016/j.compfluid.2013.06.001](https://doi.org/10.1016/j.compfluid.2013.06.001).

* **Klemp JB, Lilly DK.** <a name="klemp1978"></a> 1978. Numerical simulation of hydrostatic mountain waves. J. Atmos. Sci. 35: 78–107.

* **Lund TS, Wu X, Squires KD.** <a name="lund1998"></a> 1998. Generation of turbulent inflow data for spatially-developing boundary layer simulations. J. Comput. Phys. 140: 233–258 [doi.org/10.1006/jcph.1998.5882](https://doi.org/10.1006/jcph.1998.5882).

* **Miller MJ, Thorpe AJ.** <a name="miller1981"></a> 1981. Radiation conditions for the lateral boundaries of limited-area numerical models. Q. J. Roy. Meteor. Soc. 107: 615–628 [doi.org/10.1002/qj.49710745310](https://doi.org/10.1002/qj.49710745310).

* **Orlanski I.** <a name="orlanski1976"></a> 1976. A simple boundary condition for unbounded hyperbolic flows. J. Comput. Phys. 21: 251–269 [doi.org/10.1016/0021-9991(76)90023-1](https://doi.org/10.1016/0021-9991(76)90023-1).

* **Panofsky HA, Dutton JA.** <a name="panofsky"></a> 1984. Atmospheric Turbulence, Models and Methods for Engineering Applications. John Wiley & Sons. New York.

* **Paulson CA** <a name="paulson"></a> 1970. The mathematical representation of wind speed and temperature profiles in the unstable atmospheric surface layer. J. Appl. Meteorol. 9: 857–861 .

* **Rotach M, Gryning SE, Tassone C.** <a name="rotach1996"></a> 1996. A two-dimensional Lagrangian stochastic dispersion model for daytime condtions. Q.J.R. Meteorol. Soc. 122: 367–389 [doi.org/10.1002/qj.49712253004](https://doi.org/10.1002/qj.49712253004).

* **Tian W, Guo Z, Yu R.** <a name="tian2004"></a> 2004. Treatment of LBCs in 2D simulation of convection over hills. Adv. Atmos. Sci. 21: 573–586 [doi.org/10.1007/BF02915725](https://doi.org/10.1007/BF02915725).

* **Xie ZT, Castro IP.** <a name="xie2008"></a> 2008. Efficient generation of inflow conditions for large eddy simulation of street-scale flows. Flow, Turbul. Combust. 81: 449–470 [doi.org/10.1007/s10494-008-9151-5](https://doi.org/10.1007/s10494-008-9151-5).

