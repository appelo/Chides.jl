The compressible Navier-Stokes equations
========================================

In two dimensions the compressible Navier-Stokes equations can be written as

\[\begin{equation}
   \vec{Q}_t + \frac{\partial}{\partial x}(\vec{f}(\vec{Q})+\vec{F}(\vec{Q})) +
   \frac{\partial}{\partial y} (\vec{g}(\vec{Q})+\vec{G}(\vec{Q})) = 0.
\end{equation}\]

Where the conservative state variables are \(\vec{Q} = [\rho,
\rho u, \rho v, E]^T\) and the hyperbolic fluxes are

\[\begin{equation}
  \vec{f}(\vec{Q}) = \left( \begin{array}{c}
   \rho u\\
   \rho u^2 + p \\
   \rho u v \\
   u (E+p)      
  \end{array} \right), \ \ \  \vec{g}(\vec{Q}) = \left( \begin{array}{c}
   \rho v\\
   \rho u v  \\
   \rho v^2 + p \\
   v (E+p)      
  \end{array} \right)
\end{equation}\]

Here we assume an equation of state \(p = (\gamma-1) \rho e\) and note that we can relate \(E\) and \(p\) through \(\rho e = E - \rho (u^2+v^2)/2\).

The viscous fluxes are:

\[\begin{equation}
  \vec{F}(\vec{Q}) = - \left( \begin{array}{c}
   0\\
  \frac{4 \mu}{3} u_x+(\frac{\mu}{3}+\mu_b) v_y \\
   \mu (u_y+v_x) \\
   u F_2 + v F_3 + \frac{\mu \gamma}{Pr} T_x     
  \end{array} \right), \ \ \  \vec{G}(\vec{Q}) = - \left( \begin{array}{c}
   0\\
   \mu (u_y+v_x) \\
   \frac{4 \mu}{3} v_y + (\frac{mu}{3} + \mu_b) u_x\\
   u G_2 + v G_3 + \frac{\mu \gamma}{Pr} T_y     
  \end{array} \right).
\end{equation}\]

Here \(\rho\) is the density, \((u,v)\) is the velocity vector, \(E\) is the total energy, \(e\) is the specific internal energy, \(\mu\) is the dynamic shear viscosity, \(\mu_b\) is the bulk viscosity, \(\gamma\) is the ratio of specific heats, \(Pr\) is the Prandtl number, \(T =
(\gamma-1)e/R\) is the temperature and \(R\) is the gas constant.

Isentropic vortex
-----------------

Couette flow
------------
