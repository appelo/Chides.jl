Quickstart
==========

This section describes the plain vanilla Hermite-Taylor method for solving the scalar transport equation with periodic boundary conditions

\[\begin{eqnarray}
  && u_t = a u_x, \ \ x \in [x_L,x_R], \, t>0, \\
  && u(x,0) = f(x), \ \ u(x_L,t) = u(x_R,t).
\end{eqnarray}\]

For this particular example we choose \(x_L = 0, \, x_R = 1, \, f(x) = \sin(2 \pi x)\).

The implementation of this basic example has the same structure as the more involved examples discussed later. The main program is contained in the file [chides/transport/1d/herm1d.f90](https://bitbucket.org/appelo/chides). The main program will call subroutines that are specific to the above problem and contained in the module `pdemodule` implemented in the file [chides/transport/1d/transport1d.f90](https://bitbucket.org/appelo/chides). All examples will have these two elements, a main program (`hermXd.f90`, X being the number of dimensions) and a PDE module (`pdename.f90`). In addition, a typical implementation will call various `service` routines for polynomial manipulation, input-output, etc. For this first example we try to keep the complexity low and only use two service routines [chides/service/Hermite.f90](https://bitbucket.org/appelo/chides) related to the central task of forming a Hermite interpolating polynomial, and [chides/service/get\_Taylor1d.f90](https://bitbucket.org/appelo/chides) to turn point values into Taylor polynomials.

Grids
-----

The equations will be advanced one half time step \(\frac{\Delta t}{2}\) from a primal grid

\[\begin{equation}
   x_i = x_L + i h_x, \ \ i = 0, \ldots, n_x, \ \ h_x = (x_R-x_L)/n_x, 
\end{equation}\]

to a dual grid

\[\begin{equation}
   x_i = x_L + (i+1/2) h_x, \ \ i = 1, \ldots, n_x, \ \ h_x = (x_R-x_L)/n_x. 
\end{equation}\]

In `herm1d.f90` the generation of the grids is done in the code block:

``` {.sourceCode .fortran}
hx=(xr-xl)/DBLE(nx)
! Primal grid
DO i=0,nx
   x(i)=xl+hx*DBLE(i)
END DO
! Dual grid
DO i=1,nx
   xh(i)=xl+hx*(DBLE(i)-0.5d0)
END DO
```

Approximating the solution by a piecewise polynomial
----------------------------------------------------

Hermite methods are based on the representation of the solution \(u(x,t_n)\) at some fixed time level \(t_n\) in terms of degree \(m+1\) nodal based scaled polynomials centered around each grid point. For example on the primal grid we have that

\[\begin{equation}
u(x,t_n) \approx p_i(x) = \sum_{l = 0}^{m} c_l \left( \frac{x - x_i }{h_x}\right)^l, x \in [x_{i-\frac{1}{2}},x_{i+\frac{1}{2}}], \ \ i = 0,\ldots,n_x.
\end{equation}\]

Note that at the grid point \(x_i\) the coefficients of the polynomial are related to the derivatives of the polynomial in the following way

\[\begin{equation}
c_l = \frac{h_x^l}{l!} \left. \frac{d^l p}{dx^l} \right|_{x = x_i}.
\end{equation}\]

Thus, the degrees of freedom in the method can (interchangeably) be thought of as derivatives of the solution or coefficients of a local polynomial.

In CHIDES the coefficients \(c_l\) on each node on the primal and dual grid are stored in two arrays, again see `herm1d.f90`:

``` {.sourceCode .fortran}
DOUBLE PRECISION, DIMENSION(0:m,nvar,0:nx) :: u
DOUBLE PRECISION, DIMENSION(0:m,nvar,1:nx) :: uh
```

The leading dimension refers to the number of coefficients (or derivatives) at each node, nvar is the number of variables, i.e. 1 for a scalar problem and \>1 for a system of equations. The third dimension refers to the index along the grids. The array u holds the solution (represented by expansion coefficients / collocated derivatives) on the primal grid where the index 0 corresponds to the left boundary and the index nx to thr right boundary. The array uh holds the solution on the dual grid.

Initial data
------------

At the beginning of the computation the initial data is used to initialize the coefficients \(c_l\) by the call `CALL indat(u,x,hx,nx,m,nvar)`, with `indat` being specified in `transport1d.f90`. Here, as the initial data is a sin function computing derivatives is straightforward, but for a more complex function it may be difficult to compute the derivatives analytically. Thus we use the general purpose service routine `get_Taylor1d.f90` which computes the coefficients using polynomial interpolation:

``` {.sourceCode .fortran}
u = 0.d0
DO j=0,nx
   ! 
   CALL get_Taylor1d(u(:,1,j),x(j),m,2*m+1,hx,infun) 
   !
END DO
```

Here the inputs to get\_Taylor1d are the local coordinate `x(j)`, the degree of the Taylor polynomial `m`, the degree of the interpolant `2*m+1` chosen to guarantee convergence at the design order, the cell width `hx`, and finally the function evaluating the initial data `infun`. In this case

``` {.sourceCode .fortran}
!
  infun=SIN(freq*pi*x)
!
```

The output is the array of scaled Taylor coefficients at the jth node. The interpolation `get_Taylor1d.f90` is carried out on a Chebyshev grid of width `hx` centered at `x(j)`. .. See the description in Chapter 4 for details.

Forming the Hermite interpolating polynomial
--------------------------------------------

The first step in evolving the approximate solution to the next half timestep at a dual gridpoint \(x_{i+\frac{1}{2}}\) is to form a degree \(2m+1\) polynomial

\[\begin{equation}
p_{i-\frac{1}{2}}(x) = \sum_{l = 0}^{2m+1} c_l \left( \frac{x - x_{i-\frac{1}{2} }}{h_x}\right)^l,\ \  x \in [x_{i-1},x_{i}], \ \ i = 1,\ldots,n_x,
\end{equation}\]

from the polynomials representing the initial data at the two nearby primal gridpoints. That is, the coefficients of the new polynomial are determined by insisting that it interpolates the nearby polynomials. Functionally, this amounts to applying the linear map corresponding to Hermite interpolation, say \(A\), to the nearby polynomials to find \(p_{i-\frac{1}{2}} = A(p_{i-1},p_{i})\). In the code this just corresponds to assembling the nearby coefficients into a \((2m+2) \times 1\) vector which is multiplied by a matrix `Aintx` corresponding to \(A\). For example in the section of the code where the initial data is output to file we have:

``` {.sourceCode .fortran}
DOUBLE PRECISION, DIMENSION(0:2*m+1,0:2*m+1) :: Aintx
DOUBLE PRECISION, DIMENSION(0:2*m+1) :: udata, pplot

. 
.

   udata(0:m)       = u(0:m,1,i)
   udata(m+1:2*m+1) = u(0:m,1,i+1)
   pplot = MATMUL(Aintx,udata)
```

The matrix `Aintx` is computed in the program `herm1d.f90`

``` {.sourceCode .fortran}
CALL Hermite(m,node_left,node_right,node_interp,Aintx,0)
```

using the subroutine `Hermite` from in [chides/service/Hermite.f90](https://bitbucket.org/appelo/chides).

Evolving the Hermite interpolating polynomial
---------------------------------------------

To evolve the above polynomial we consider the space-time approximation

\[\begin{equation}
p^n_{i-\frac{1}{2}}(x,t) = \sum_{l=0}^{2m+1} \sum_{s=0}^{q} d_{ls} \left( \frac{x - x_{i+\frac{1}{2}} }{h_x}\right)^l \left( \frac{t - t_n }{\Delta t}\right)^s,
\end{equation}\]

which when evaluated at time \(t = t_n\) reduces to

\[\begin{equation}
p_{i-\frac{1}{2}}(x,t_n) = \sum_{l = 0}^{2m+1} d_{l0} \left( \frac{x - x_{i-\frac{1}{2} }}{h_x}\right)^l,\ \  x \in [x_{i-1},x_{i}], \ \ i = 1,\ldots,n_x,
\end{equation}\]

and we thus have \(d_{l0} = c_l\).

To find the remaining coefficients \(d_{ls}, s > 0\) we will use the PDE. In particular we will insist that the approximation \(p_{i-\frac{1}{2}}(x,t)\) satisfies the PDE at \((x^\ast,t_n)\). For notational convenience consider the approximation at one dual gridpoint, say \(x^\ast\) and denote it by \(p(x^\ast,t_n)\). Then we have that

\[\begin{eqnarray}
  \frac{\partial p(x^\ast,t_n)}{\partial x} = \frac{d_{10}}{h_x},\\
  \frac{\partial p(x^\ast,t_n)}{\partial t} = \frac{d_{01}}{\Delta t},
\end{eqnarray}\]

thus from the PDE we find that we may express the first time-derivative in the Taylor series

\[\begin{equation}
  \frac{\partial p(x^\ast,t_n)}{\partial t} = a \frac{\partial p(x^\ast,t_n)}{\partial x} \ \ \Rightarrow \ \ d_{01} = a \frac{{\Delta t}}{h_x} d_{10}.
\end{equation}\]

By repeated differentiation in space and time of the PDE

\[\begin{equation}
  \frac{\partial^{s+l} u}{\partial t^s \partial^l} = a \frac{\partial^{s-1+l+1} u}{\partial t^{s-1} \partial x^{l+1}},
\end{equation}\]

which upon replacing \(u\) with \(p\) at \((x^\ast,t_n)\) yields the recursion relation:

\[\begin{equation}
  d_{sl} = a \frac{{\Delta t}}{h_x} \frac{l+1}{s} d_{s-1,l+1},, \ \ s = 1,\ldots, \ \ l = 0,\ldots,
\end{equation}\]

specifying all time derivatives in the Taylor series in terms of spatial derivatives.

In the code this update is done in the main time-stepping loop

``` {.sourceCode .fortran}
DO it=1,nsteps
   ! first half step
   t = DBLE(it-1)*dt
   ! At the end of the first half step all entries of
   ! uh will have been updated
   uh = 0.d0
   ! Loop over all dual gridpoints
   DO i=1,nx

      ! tcofs contains the coefficients c_{l,s}
      ! l = 0,2*m+1, s = 0,q
      tcofs=0.d0

      ! c_{l,0}, l = 0,2*m+1 are found by interpolating data
      ! from adjacent primal nodes
      udata(0:m)       = u(:,1,i-1)
      udata(m+1:2*m+1) = u(:,1,i)
      tcofs(0:2*m+1,0,1)=MATMUL(Aintx,udata)

      ! The rest of the coefficients are found via the
      ! recursion:
      ! c_{l,s} = dt/h*(l+1)/s * c_{l+1,s-1}
      ! This is done in the subroutine pde
      CALL pde(tcofs,2*m+1,q,nvar,cflx,etax,xh(i),t)
```

Filling in all the ''time coefficients'' in `tcofs` is done inside the subroutine `pde` from [chides/transport/1d/herm1d.f90](https://bitbucket.org/appelo/chides):

``` {.sourceCode .fortran}
SUBROUTINE pde(tcofs,mp,q,nvar,cflx,etax,x,t)
  !
  ! This pde implements the recursion relation
  ! associated with the linear transport equation 
  ! u_t = a*u_x
  !
  implicit none
  INTEGER, INTENT(IN) :: mp,q,nvar
  DOUBLE PRECISION, DIMENSION(0:mp,0:q,nvar), INTENT(INOUT) :: tcofs
  DOUBLE PRECISION, INTENT(IN) :: cflx,etax,x,t
  INTEGER :: idx,idt
  DOUBLE PRECISION :: a_dt_h
  !
  tcofs(:,1,:)=0.d0
  a_dt_h = a*cflx
  DO idt=1,q
     DO idx=0,mp-idt
        tcofs(idx,idt,1)=a_dt_h*DBLE(idx+1)/DBLE(idt)*tcofs(idx+1,idt-1,1)
     END DO
  END DO
END SUBROUTINE pde
```

With all the time-coefficients available we can evaluate the approximation \(p(x,t)\) at a later time \((x^\ast,t_n+\frac{\Delta t}{2})\) for the first \(m\) spatial derivatives yielding the new approximation of \(u(x,t_{n+\frac{1}{2}})\):

\[\begin{equation}
     u(x,t_{n+\frac{1}{2}}) \approx  p^n_{i-\frac{1}{2}}(x) = \sum_{l=0}^{m} \sum_{s=0}^{q} {\tilde{c}}_{l} \left( \frac{x - x_{i+\frac{1}{2}} }{h_x}\right)^l,  \ \  x \in [x_{i-1},x_{i}], \ \ i = 1,\ldots,n_x,
\end{equation}\]

where the new coefficients \(\tilde{c}_l\) are

\[\begin{equation}
  \tilde{c}_l = \sum_{s = 0}^{q}  d_{ls} \left( \frac{1}{2} \right)^s, \ \  l =0,\ldots,m. 
\end{equation}\]

In the code this looks like:

``` {.sourceCode .fortran}
! When all coefficients are known
! the Taylor series can be evaluated
! and uh(i) updated
DO idx=0,m
   uh(idx,1,i) = tcofs(idx,0,1)
   dfact = 1.d0
   DO l = 1,q
      dfact = dfact*0.5d0
      uh(idx,1,i) = uh(idx,1,i) &
           + dfact*tcofs(idx,l,1)
   END DO
END DO
```

Second half timestep and boundary conditions
--------------------------------------------

The second half timestep follows the same procedure as described above:

> 1.  Form an interpolating polynomial
>
>     ``` {.sourceCode .fortran}
>     ! second half step
>     ! Same thing as above but we now update
>     ! the interor of u
>     t = t + 0.5d0*dt
>     DO i = 1,nx-1
>        tcofs = 0.d0
>        udata(0:m)       = uh(:,1,i)
>        udata(m+1:2*m+1) = uh(:,1,i+1)
>        tcofs(0:2*m+1,0,1) = MATMUL(Aintx,udata)
>     ```
>
> 2.  Then use the PDE and update the approximation at the full timestep
>
>     ``` {.sourceCode .fortran}
>     CALL pde(tcofs,2*m+1,q,nvar,cflx,etax,x(i),t)
>
>      ! Add up the Taylor series
>      DO idx=0,m
>         u(idx,1,i)=tcofs(idx,0,1)
>         dfact=1.d0
>         DO l=1,q
>            dfact=dfact*.5d0
>            u(idx,1,i)=u(idx,1,i) &
>                 +dfact*tcofs(idx,l,1)
>         END DO
>      END DO
>     END DO
>     ```
>
For this example the boundary conditions are taken to be periodic. This means that the interior update procedure can be applied on the boundary as well, with the slight difference that when the interpolating polynomial is formed, data ''from outside the domain'' is replaced by the periodic copy:

``` {.sourceCode .fortran}
! now the left and right bcs
! left
tcofs=0.d0
udata(0:m)=uh(:,1,nx)
udata(m+1:2*m+1)=uh(:,1,1)
tcofs(0:2*m+1,0,1)=MATMUL(Aintx,udata)
CALL pde(tcofs,2*m+1,q,nvar,cflx,etax,x(0),t)
! Add up the Taylor series
DO idx=0,m
   u(idx,1,0)=tcofs(idx,0,1)
   dfact=1.d0
   DO l=1,q
      dfact=dfact*.5d0
      u(idx,1,0)=u(idx,1,0) &
           +dfact*tcofs(idx,l,1)
   END DO
END DO
! Right
u(:,:,nx)=u(:,:,0)
```

Compiling and Running the program
---------------------------------

In the directory `chides/transport/1D` there is a simple `makefile` that can be used to compile the program by simply typing `make` in a terminal window. This assumies you are standing in the directory `chides/transport/1D` and that your system has `gfortran` installed. If you are using another compiler simply change the makefile.

If you are successful in building the solver you should have an executable named `herm1d.x` that you can run.

Displaying the results
----------------------

After running the program you should see something like:

``` {.sourceCode .none}
bash-3.2$ ./herm1d.x 
     Space step,     Max-error
     5.00000E-02     6.60343E-10
```

The program should also have generated a number of snapshots of the solution, `u0000000.txt`, `u0000001.txt`, etc. These can be displayed in Matlab using the script `plot_sol.m`.

Changing the order and number of gridpoints
-------------------------------------------

All the parameters associated with a simulation are contained in the top block of the `pdemodule` (residing in `transport1d.f90`):

``` {.sourceCode .fortran}
MODULE pdemodule
 !
 !
 !  Problem setup for solving 1d transport equation:
 !  u_t = a u_x
 !  With periodic boundary conditions
 !
 ! 
 IMPLICIT NONE
 INTEGER, PARAMETER :: m = 3        ! spatial order 2m+1
 INTEGER, PARAMETER :: q = 2*m+2    ! terms in the Taylor series
 INTEGER, PARAMETER :: nx = 20      ! Discretization uses nx+1 points from xl to xr 
 INTEGER, PARAMETER :: nxp = 500    ! Number of points in the output 

 DOUBLE PRECISION, PARAMETER :: CFL = 0.95d0 ! Ratio dt/h
 INTEGER, PARAMETER :: nout = 1     ! output the solution every nout times
 INTEGER, PARAMETER :: nvar = 1     ! 

 DOUBLE PRECISION, PARAMETER :: a = 1.0d0 ! Wave speed 

 DOUBLE PRECISION, PARAMETER :: ttot = 10.0d0 ! Final time
 DOUBLE PRECISION, PARAMETER :: xl = 0.0d0   ! Domain boundaries
 DOUBLE PRECISION, PARAMETER :: xr = 1.0d0
 ! 
 ! Other constants
 DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932385d0
 DOUBLE PRECISION, PARAMETER :: freq = 2.d0
 ! 
CONTAINS
```

This is where you would change the order and the number of gridpoints etc.
