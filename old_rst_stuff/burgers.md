Viscous Burgers' equation in one dimension
==========================================

Consider the non-linear conservation law:

\[\begin{equation}
u_t+f(u)_x=0,
\end{equation}\]

with the viscous Burgers' equation as a special case with \(f(x)
= \frac{u^2}{2}-\nu u_x\). As before we will take derivatives in space and time of the equation and insist that the resulting relations hold for the approximation. Typically the flux function will be non-linear and we must find an polynomial approximation of it using the polynomial manipulation routines that are part of CHIDES. As the non-linearity depends on time it becomes less favorable to use a Taylor series expansion in time and instead we opt to use a one-step method to evolve the solution and its derivatives. Among other things this example illustrate:

> 1.  How to use one-step time-integrators, for example a Runge-Kutta method together with Hermite interpolation.
> 2.  How to improve the accuracy and time-step restrictions of the method by taking multiple sub-timesteps per half timestep.
> 3.  How to discretize second order derivatives.
> 4.  Semi-conservative and non-conservative discretizations of conservation laws.

The source files reside in `chides/burgers/1D`

To find an approximation to the solution of the above conservation law we assume we have interpolated two degree \(m\) polynomials at two adjacent gridpoints on the primal grid. That is, at a dual gridpoint \(x_{i+\frac{1}{2}}\) is to we start with the degree \(2m+1\) polynomial

\[\begin{equation}
p_{i-\frac{1}{2}}(x,t)= \sum_{l = 0}^{2m+1} c_l(t) \left( \frac{x - x_{i-\frac{1}{2} }}{h_x}\right)^l,\ \  x \in [x_{i-1},x_{i}], \ \ i = 1,\ldots,n_x.
\end{equation}\]

Note that we now assume that the coefficients \(c_l(t)\) depend on time so that if we compute one time derivative and \(l\) derivatives in space and evaluate at \(x = x_{i-\frac{1}{2}}\) we get terms of the kind \(c^\prime_l(t)\) that we can be identified to be scaled approximations of \(\frac{\partial^l u_t(x,t)}{\partial x^l}\) at the dual grid point.

Evaluating the non-linearity
----------------------------

As we obtained the time derivatives \(c^\prime_l(t)\) by differentiating the approximation it is natural evolve them according to the system of ODE that is obtained by differentiating the conservation law in space and insisting that the approximation satisfies those equations at \(x = x_{i-\frac{1}{2}}\). At first this may appear as a cumbersome proposition as we must compute multiple derivatives of \(f(u)_x\). However, if we consider the flux to be approximated by a polynomial:

\[\begin{equation}
f(u) \approx \sum_{l = 0}^{2m+1} a_l(t) \left( \frac{x - x_{i-\frac{1}{2} }}{h_x}\right)^l,\ \  x \in [x_{i-1},x_{i}], \ \ i = 1,\ldots,n_x,
\end{equation}\]

we can immediately see that the system of ODE governing \(c_l(t)\) is:

\[\begin{eqnarray}
   c^\prime_l(t) & = & -\frac{(l+1)}{h_x} a_{l+1}(t),\ \ l = 0,\ldots,2m, \\
   c^\prime_{2m+1}(t) & = & 0.
\end{eqnarray}\]

For the particular example we consider here we have that \(f(u) =
\frac{u^2}{2}-\nu u_x\) which is computed in the routine `compute_flux(u,ut,mp,nvar,hx)`

``` {.sourceCode .fortran}
SUBROUTINE compute_flux(u,flx,mp,nvar,hx)
 !
 ! compute (f(u)-nu*u_x)
 !
 implicit none
 INTEGER, INTENT(IN) :: mp,nvar
 DOUBLE PRECISION, DIMENSION(0:mp,nvar), INTENT(IN) :: u
 DOUBLE PRECISION, DIMENSION(0:mp,nvar), INTENT(OUT) :: flx
 DOUBLE PRECISION, INTENT(IN) :: hx
 INTEGER :: idx
 DOUBLE PRECISION, DIMENSION(0:mp) :: q1,p1
 ! 
 ! Compute the hyperbolic part of the flux
 CALL tpolymul(u(:,1),u(:,1),p1,mp)
 flx(:,1) = 0.5d0*p1
 ! Then add the viscous part
 DO idx=0,mp-1
    q1(idx) = dble(idx+1)*u(idx+1,1)
 END DO
 q1(mp) = 0.d0 
 flx(:,1) = flx(:,1) - (nu/hx)*q1
 ! 
END SUBROUTINE compute_flux
```

In the subroutine we see the use of `tpolymul(u(:,1),u(:,1),p1,mp)` which is used to compute the square of the polynomial approximation.

After having found the flux it is straightforward to find (\(\Delta t\) times) the coefficients \(a_l\) by differentiating the flux as is done in the subroutine `pde(fp,flx,mp,nvar,cflx)`

``` {.sourceCode .fortran}
fp = 0.d0
DO idx=0,mp-1
   fp(idx,1)  = -cflx*DBLE(idx+1)*flx(idx+1,1)
END DO
```

The output of this subroutine can then be used to evolve the interpolant \(p(x,t)\) over half a timestep.

The driver `herm1d_non_cons`
----------------------------

An implementation based on the above description can be found in `/chides/burgers/1d/herm1d_non_cons.f90`. As before the basic steps of the algorithm are:

> 1.  Setup the initial data
>
> > ``` {.sourceCode .fortran}
> > CALL indat(u,x,hx,nx,m,nvar) 
> > ```
>
> 1.  Form the Hermite interpolant.
>
>     ``` {.sourceCode .fortran}
>     DO i=1,nx
>     ```
>
>     > ! Interpolate solution and fluxes to the cell center call Taylor\_to\_Hermite1d(uhermite,u(:,:,i-1:i),Aintx,m,nvar) ...
>
> 2.  And evolution of `uhermite` using a Runge-Kutta method
>
>     ``` {.sourceCode .fortran}
>     call Taylor_to_Hermite1d(uhermite,u(:,:,i-1:i),Aintx,m,nvar)
>       do irk_sub = 1,rksub
>          urk = uhermite
>          do irk=1,rkstages ! do a local RK loop
>             urk = uhermite
>             ! if irk > 1 we add on previous data
>             do ii = 1,irk-1
>                urk = urk + ButA(irk,ii)*k_stages(:,:,ii)
>             end do
>             call compute_flux(urk,flxcofs,2*m+1,nvar,hx)
>             call pde(ut,flxcofs,2*m+1,nvar,cflx)
>             k_stages(:,:,irk) = ut
>          end do
>       ! At the end we add up the stages
>       do ii = 1,rkstages
>          uhermite = uhermite + ButB(ii)*k_stages(:,:,ii)
>       end do
>     end do
>     uh(:,:,i) = uhermite(0:m,1:nvar)
>     ```
>
>     This section of code evolves the interpolant using a the Runge-Kutta method defined by the Butcher arrays `ButA` and `ButB` that have been initialized by the `call set_butcher_tables`. The time stepping performs `rksub` local timesteps with the local timestep `dts = 0.5d0*dt/DBLE(rksub)` to update `uh`.
>
Local timestepping and CFL constraints
--------------------------------------

The purpose of this local timestepping is to maximize the computation-to-communication ratio. For linear problems with exact Taylor timestepping we have already seen that we can run the solvers at the full physically allowable CFL number. For non-linear problems with approximate timestepping we can use large timesteps if the local evolution is accurate enough.

For this example we can choose between three different RK solvers in the module `rkmodule.f90`, the classic fourth order Runge-Kutta, the eight order DORPRI8 and the third order SSP(3,3) by Shu and Osher. To illustrate the effect of local timestepping we consider the problem defined by:

> ``` {.sourceCode .fortran}
> INTEGER, PARAMETER :: m = 9 !spatial order 2m+1
> INTEGER, PARAMETER :: nx = 100
> DOUBLE PRECISION, PARAMETER :: xl = -1d0
> DOUBLE PRECISION, PARAMETER :: xr =  1d0
> DOUBLE PRECISION, PARAMETER :: ttot = 0.2d0
> DOUBLE PRECISION, PARAMETER :: nu = 0.0001d0
> ```

for different `CFL` and `rksub`. With DORPRI8 we can solve the problem with one single local timestep at a `CFL=0.65`, if we use RK4 we must take `CFL=0.30` and if we use SSP(3,3) we must take `CFL=0.20`. However if we perform 6 local timesteps for SSP(3,3) or 4 local timesteps with RK4 or 2 local timesteps with DORPRI8 we can use `CFL=0.80`.

As this problem also has a viscous term we have choose the timestep according to

> ``` {.sourceCode .fortran}
> dt = min(CFL*hx/maxval(abs(u(0,1,:))),CFL_V*hx**2/(nu*m))
> ```

where `CFL_V` typically be chosen as 0.1.

The driver `herm1d_semi_cons`
-----------------------------

The driver `herm1d_semi_cons.f90` illustrates an alternative approach to discretize conservation laws. The alterative discretization is motivated by the control volume statement of the conservation law obtained by integrating over a space-time cell

\[\begin{equation}   
\int_{0}^{\Delta t } \int_{x_i}^{x_{i+1}} u_t + f(u)_x \,  dx dt = 0,
\end{equation}\]

which can be simplified into:

\[\begin{equation}
\int_{x_i}^{x_{i+1}} u(\Delta t,x) - u(0,x) dx = - \int_{0}^{\Delta t } \left[ f(u(t,x_{i+1})) -f(u(t,x_{i})) \right] \,dt. 
\end{equation}\]

Now for a problem posed on a periodic domain (or with zero flux boundary conditions) we have that the global mass is conserved

\[\begin{equation}
\sum_i \int_{x_i}^{x_{i+1}} u(\Delta t,x) dx = \sum_i \int_{x_i}^{x_{i+1}} u(0,x) dx.
\end{equation}\]

For a numerical method to be able to have a similar conservation property we must require that the flux from the left and right of an interface coincide. Unless the flux is a linear function this is not true in general for the above formulation. In the alternative formulation in `herm1d_semi_cons.f90` we compute the fluxes at the faces of each cell and interpolate them to the cell center. This will guarantee that the fluxes are continuous at interfaces. The interpolated fluxes are then differentiated to compute the first slope in the Runge-Kutta method.

``` {.sourceCode .fortran}
DO i=0,nx
     call compute_flux(u(:,:,i),flx(:,:,i),m,nvar,hx)
  END DO
  ! Add the viscous flux
  DO i=1,nx
     ! Semi-conservative
     ! Interpolate solution and fluxes to the cell center
     call Taylor_to_Hermite1d(uhermite,u(:,:,i-1:i),Aintx,m,nvar)
     call Taylor_to_Hermite1d(flxcofs,flx(:,:,i-1:i),Aintx,m,nvar)
     do irk_sub = 1,rksub
        if(irk_sub.gt.1) then
           urk = uhermite
           call compute_flux(urk,flxcofs,2*m+1,nvar,hx)
        end if
        ! Get the first slope
        call pde(ut,flxcofs,2*m+1,nvar,cflx)
        k_stages(:,:,1) = ut
        DO irk=2,rkstages ! do a local RK loop
           urk = uhermite
           ! if irk > 1 we add on previous data
           do ii = 1,irk-1
              urk = urk + ButA(irk,ii)*k_stages(:,:,ii)
           end do
           call compute_flux(urk,flxcofs,2*m+1,nvar,hx)
           call pde(ut,flxcofs,2*m+1,nvar,cflx)
           k_stages(:,:,irk) = ut
        end do
        ! At the end we add up the stages
        do ii = 1,rkstages
           uhermite = uhermite + ButB(ii)*k_stages(:,:,ii)
        end do
     end do
     uh(:,:,i) = uhermite(0:m,1:nvar)
  END DO
```

Although this formulation requires one additional interpolation it is more robust than the non-conservative formulation. For the latter multi dimensional examples with Euler's equations and compressible Navier-Stokes equations, this is the formulation we will use.
