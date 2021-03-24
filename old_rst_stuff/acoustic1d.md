Acoustics in one dimension
==========================

The linearized acoustic equations can be derived from Navier-Stokes equations (which we solve in the example navierstokes).

Considering an ideal gas with negligible viscosity Navier-Stokes equations are reduced to the compressible Euler equations:

\[\begin{eqnarray}
&&\frac{D p}{Dt} +\gamma p \nabla \cdot \vec{v} = 0,\\
&&\frac{D \vec{v}}{Dt} + \frac{1}{\rho} \nabla p = 0, \\
&&\frac{D \rho}{Dt} + \rho \nabla \cdot \vec{v} = 0. 
\end{eqnarray}\]

Here we consider the one dimensional case and small perturbations around a known flow state \((1/\gamma, a(x), \rho_0(x))\). That is we replace

\[\begin{equation}
p \rightarrow 1/\gamma + p(x,t), \ \ 
v \rightarrow a(x) + v(x,t), \ \
\rho \rightarrow \rho_0(x) + \rho(x,t),  
\end{equation}\]

and neglect small quadratic terms. The resulting linearized equations are

\[\begin{eqnarray}
p_t + a(x) p_x + v_x &=& - a'(x), \\ 
v_t + a(x) v_x + \frac{p_x}{\rho_0(x)} &=& - a(x) a'(x), \\
\rho_t + a(x) \rho_x + \rho_0(x) v_x  &=& - \rho_0'(x) a(x) - a'(x)
\rho_0(x).
\end{eqnarray}\]

This module considers three different cases

> 1.  Zero mean flow, \(a(x) = 0\) and constant mean density, \(\rho_0(x) = 1\).
> 2.  Constant mean flow, \(a(x) = a\) and constant density, \(\rho_0(x) = 1\).
> 3.  Variable mean flow with variable density.

Case 1, zero mean flow and constant mean density
------------------------------------------------

For this case the equations for the pressure and velocity are reduced to

\[\begin{eqnarray}
p_t  &=& - v_x, \\ 
v_t  &=& - p_x.
\end{eqnarray}\]

The evolution of the above system using wall boundary conditions \(v = 0\) on both sides (which by using the PDE also means \(p_x = 0\)) and with the exact solution

\[\begin{eqnarray}
&&p(x,t) = \cos(k \pi t) \cos(k \pi x), \\
&&v(x,t) = \sin(k \pi t) \sin(k \pi x), 
\end{eqnarray}\]

is implemented in the file `/chides/acoustic/1D/acoustics1d_wall.f90`

> 1.  **Add description of PDE.**
> 2.  **Add description of enforcing BC.**

Case 2, Constant mean flow, \(a(x) = a\) and constant density \(\rho_0(x) = 1\)
-------------------------------------------------------------------------------

For this case the equations for the pressure and velocity are reduced to

\[\begin{eqnarray}
p_t + a p_x + v_x &=& 0, \\ 
v_t + a v_x + p_x &=& 0, \\
\end{eqnarray}\]

For this problem we consider periodic boundary conditions, as implemented in `/chides/acoustic/1D/acoustics1d_per.f90` and characteristics based non-reflecting boundary conditions as implemented in `/chides/acoustic/1D/acoustics1d_cbc.f90`.

> 1.  **Add description of CBC, equations.**
> 2.  **Add description of CBC, how the extrapolation is done.**

Case 3, Varying mean flow and density with source terms
-------------------------------------------------------

For this case the equations for the pressure and velocity are reduced to

\[\begin{eqnarray}
p_t + a(x) p_x + v_x &=& - a'(x) + \tilde{f}_p(x,t), \\ 
v_t + a(x) v_x + \frac{p_x}{\rho_0(x)} &=& - a(x) a'(x) + \tilde{f}_v(x,t), \\
\end{eqnarray}\]

which we reformulate as

\[\begin{eqnarray}
p_t + a(x) p_x + v_x &=& f_p(x,t), \\ 
v_t + a(x) v_x + s(x) p_x &=& f_v(x,t), 
\end{eqnarray}\]

For this problem we consider periodic boundary conditions, as implemented in `/chides/acoustic/1D/acoustics1d_var.f90`.

> 1.  **Add description of how** `cofs(a,s,f_p,f_v,mp,dx,x,t)` **is used.**
> 2.  **Explain polynomial multiplication.**
> 3.  **Explain** `point_to_Taylor`.
> 4.  **Explain the RK4 procedure.**

