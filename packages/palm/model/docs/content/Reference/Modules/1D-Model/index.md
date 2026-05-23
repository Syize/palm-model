## 1-D model for providing initial vertical velocity profiles

The initial profiles of the horizontal wind components in PALM can be prescribed by the user by piecewise linear gradients or by directly using observational data. Alternatively, a 1-D model can be employed to calculate stationary boundary-layer wind profiles as a better first guess of initial conditions. The arrays of the 3-D variables are then initialized with the (stationary) solution of the 1-D model. These variables are $u_i$, with $i$ ∈ {1, 2}, $e$, $K_h$, $K_m$, and, with MOST applied between the surface and the first vertical grid level, also $L$, $u_∗$, as well as $\overline{u_i^{\prime\prime} u_3^{\prime\prime}}$ (with $i$ ∈ {1, 2}). The 1-D model is using a RANS turbulence parameterization. If the 3-D model is run in RANS mode, too, the first guess provided by the 1-D model is usually very good.

The 1-D model assumes the profiles of $\theta$ and $q_v$, which need to be prescribed, to be constant in time. The model solves the prognostic equations for $u_i$


$\quad \begin{align*}
  & \frac{\partial u_i}{\partial t} = -\varepsilon_{i3j}f_3 u_j +
  \varepsilon_{i3j}f_3 {u_{\mathrm{g},j}} - \frac{\partial
    \overline{u_i^{\prime\prime}u_3^{\prime\prime}}}{\partial x_3} \; .
\end{align*}$


The default turbulence parameterization is based on the turbulence kinetic energy $e$. The prognostic equation for $e$ reads

$\quad \begin{align*}
  & \frac{\partial e}{\partial t} = - \overline{u^{\prime\prime}w^{\prime\prime}} \frac{\partial u}{\partial z}
                                    - \overline{v^{\prime\prime}w^{\prime\prime}} \frac{\partial v}{\partial z}
                                    + \frac{g}{\theta} \overline{w^{\prime\prime}\theta^{\prime\prime}}
                                    - \frac{\partial \overline{w^{\prime\prime}e^{\prime\prime}}}{\partial z}
                                    - \epsilon\;.
\end{align*}$

The dissipation rate is parametrized by

$\quad \begin{align*}
  & \epsilon = 0.064 \frac{e^{\frac{3}{2}}}{l}
\end{align*}$

after Detering and Etling (1985). The mixing length is calculated after Blackadar (1997) as

$\quad \begin{align*}
  & l = \frac{\kappa z}{1 + \frac{\kappa
      z}{l_{\text{Bl}}}}\,\;\text{with}\,\;l_{\text{Bl}} = \frac{2.7\times10^{-4}}{f}\sqrt{u_{\mathrm{g}}^2 + v_{\mathrm{g}}^2}\;.
\end{align*}$

The turbulent fluxes are calculated using a gradient approach (first-order closure):

$\quad \begin{align*}
  & \overline{u_i^{\prime\prime}u_3^{\prime\prime}} = - K_\mathrm{m}
  \frac{\partial u_i}{\partial
    x_3}\;,\;\overline{w^{\prime\prime}\theta^{\prime\prime}} = -
  K_\mathrm{h} \frac{\partial \theta}{\partial z}
  \;,\,\overline{w^{\prime\prime}e^{\prime\prime}} = - K_\mathrm{m}
  \frac{\partial e}{\partial z}\;,
\end{align*}$

where $K_m$, and $K_h$ are calculated as

$\quad \begin{eqnarray}
 K_\mathrm{m} & = & c_\mathrm{m}\;\sqrt{e}\
  \left\{\begin{array}{ll}         l & \text{for} \,\,  \textit{Ri} \geq 0 \\
    l_{\text{Bl}} & \text{for} \,\, \textit{Ri} < 0
         \end{array}
  \right. \\
  K_\mathrm{h} & = & \frac{\Phi_\mathrm{h}}{\Phi_\mathrm{m}} K_\mathrm{m}
\end{eqnarray}$

with the similarity functions $\Phi_\mathrm{h}$, and $\Phi_\mathrm{m}$ (see Eqs. in Section [boundary conditions](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/bc)), and the gradient Richardson number

$\quad \begin{align*}
  & \textit{Ri} =
  \frac{\frac{g}{\theta_\mathrm{v}}\frac{\partial\theta}{\partial
      z}}{\left[\left(\frac{\partial u}{\partial z}\right)^2 +
      \left(\frac{\partial v}{\partial z}\right)^2 \right]} \cdot
  \begin{cases}
    1\,\;&\text{for~}\,\;\textit{Ri} \geq 0\;,\\
    (1 - 16 \cdot
    \textit{Ri})^{\frac{1}{4}}\,\;&\text{for~}\,\;\textit{Ri} < 0\;.
  \end{cases}
\end{align*}$

Note that the distinction of cases in the Eq. above is done using the value of $Ri$ from the previous time step.

As for the 3-D model, a turbulence model based on $e$ and a prognostic equation for dissipation $\epsilon$ is available, too.

Moreover, a Rayleigh damping can be switched on to speed up the damping of inertial oscillations. The 1-D model is discretized in space using finite differences. Discretization in time is achieved using the 3rd-order Runge--Kutta time-stepping scheme (Williamson, 1980). In order to avoid very small time steps forced by the diffusion time step criterion, the implicit Crank-Nicolson time step scheme can be used instead for treating all diffusion terms that appear in the prognostic equations. As a part of the Crank-Nicolson scheme, the algorithm of Stone (1973) is used to solve the tridiagonal systems of equations.

Dirichlet boundary conditions are used at the top and bottom boundaries of the model, except for $e$, for which Neumann conditions are set at the surface (see also Section [boundary conditions](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/bc)).


## References

- Blackadar, A.K. (1997): Turbulence and Diffusion in the Atmosphere. Springer. Berlin, Heidelberg, New York, 185 pp.

- Detering, H.W., Etling D. (1985): Application of the E-epsilon turbulence model to the atmospheric boundary layer. Boundary-Layer Meteorol., 33, 113–133.

- Stone, H.S. (1973): An Efficient Parallel Algorithm for the Solution of a Tridigonal Linear System of Equations. Journal of the Association for Computing Machinery, 20, 27-38.

- Williamson, J.H. (1980). Low-storage Runge–Kutta schemes. J. Comput. Phys., 35, 48–56.
