The acoustic equations in two dimensions
========================================

This example illustrates how to implement a two dimensional version of the Hermite-Taylor method (see quickstart for an example in 1D) for solving a linear first order hyperbolic system. The source files reside in `chides/acoustic/2D`

We will consider the equations of acoustics which can be derived from the compressible Navier-Stokes equations. Assuming an isentropic, inviscid, quiscent flow and small perturbations the conservation of mass and momentum are reduced to

\[\begin{eqnarray}
  \frac{\partial \rho}{\partial t} + \rho_0 \nabla \cdot \vec{u} = 0,\\
  \frac{\partial \vec{u}}{\partial t} + \frac{c_0^2}{\rho_0} \nabla \rho = 0.
\end{eqnarray}\]

Here we will take \(\rho_0 = c_0 = 1\).

To make the equations look more like the implementation we write the above system in terms of the density and the velocities \(u\) and \(v\)

\[\begin{eqnarray}
  &&\frac{\partial \rho}{\partial t} + \frac{\partial u}{\partial x}+ \frac{\partial v}{\partial y} = 0,\\
  &&\frac{\partial u}{\partial t} + \frac{\partial \rho}{\partial x}
  = 0, \\
  &&\frac{\partial v}{\partial t} + \frac{\partial \rho}{\partial y} = 0.
\end{eqnarray}\]

A first example in two dimensions with periodic boundary conditions
-------------------------------------------------------------------

As a first example of how CHIDES can be used to construct multi-dimensional solvers we consider the above equations in a periodic box \((x,y) \in [-\pi,\pi]^2\). Then an exact solution is

\[\begin{eqnarray}
  p(x,y,t) &=& \sin(x) \sin(y) \sin(\sqrt{2} t),  \\
  u(x,y,t) &=& \frac{1}{\sqrt{2}} \cos(x) \sin(y) \cos(\sqrt{2} t),\\
  v(x,y,t) &=& \frac{1}{\sqrt{2}} \sin(x) \cos(y) \cos(\sqrt{2} t).
\end{eqnarray}\]

The code structure is similar to the one-dimensional case. Here the main program is contained in the file [chides/acoustic/2D/herm2d\_per.f90](https://bitbucket.org/appelo/chides/src) and the module describing the PDE, initial data and discretization related parameters resides in the file [chides/acoustic/2D/acoustic\_per.f90](https://bitbucket.org/appelo/chides/src).

### Grids

As in the one-dimensional case the equations will be advanced one half time step \(\frac{\Delta t}{2}\) from a primal grid

\[\begin{eqnarray}
   && x_i = x_L + i h_x, \ \ i = 0, \ldots, n_x, \ \ h_x = (x_R-x_L)/n_x,\\ 
   && x_j = y_B + j h_y, \ \ j = 0, \ldots, n_y, \ \ h_y = (y_T-y_B)/n_y, 
\end{eqnarray}\]

to a dual grid

\[\begin{eqnarray}
   x_i = x_L + (i+1/2) h_x, \ \ i = 1, \ldots, n_x, \ \ h_x = (x_R-x_L)/n_x,\\
   y_j = y_B + (j+1/2) h_y, \ \ j = 1, \ldots, n_y, \ \ h_y = (y_T-y_B)/n_y. 
\end{eqnarray}\]

In `herm2d_per.f90` this corresponds to the code block:

``` {.sourceCode .fortran}
hx=(xr-xl)/DBLE(nx)
DO i=0,nx
   x(i)=xl+hx*DBLE(i)
   IF (i > 0) THEN
      xh(i)=xl+hx*(DBLE(i)-.5d0)
   END IF
END DO
hy=(yt-yb)/DBLE(ny)
DO i=0,ny
   y(i)=yb+hy*DBLE(i)
   IF (i > 0) THEN
      yh(i)=yb+hy*(DBLE(i)-.5d0)
   END IF
END DO
```

### Approximating the solution by a piecewise polynomial

In two dimensions we represent the solution of one of the fields (the density and the velocities in the x and y dierction) by a polynomial. For example for the velocity in the x direction we have \(u(x,y,t_n)\) at a fixed time level \(t_n\) in terms of degree \(m+1\) nodal based scaled tensor product polynomials centered around each grid point. For example on the primal grid we have that

\[\begin{equation}
u(x,y,t_n) \approx p_{ij}(x,y) = \sum_{k = 0}^{m} \sum_{l = 0}^{m}
c_{kl} \left( \frac{x - x_i }{h_x}\right)^k
\left( \frac{y - y_j }{h_y}\right)^l.
\end{equation}\]

Again, at a grid point \((x_i,y_j)\) the coefficients of the polynomial can be expressed in terms of derivatives

\[\begin{equation}
c_{kl} = \frac{h_x^k}{k!} \frac{h_y^l}{l!}\left. \frac{d^{k+l} p}{dx^k dy^l} \right|_{x = x_i,y = y_j}.
\end{equation}\]

In CHIDES the coefficients \(c_{kl}\) on each node on the primal and dual grid are stored in two arrays, see `herm2d_per.f90`:

``` {.sourceCode .fortran}
DOUBLE PRECISION, DIMENSION(0:m,0:m,nvar,0:nx,0:ny) :: u
DOUBLE PRECISION, DIMENSION(0:m,0:m,nvar,0:nx+1,0:ny+1) :: uh
```

The two leading dimension refers to the number of coefficients (or derivatives) at each node, nvar is the number of variables, e.g. 3 for the acoustic system. The last two dimensions refers to the indicies along the grids.

### Initial data

At the beginning of the computation the initial data is used to initialize the coefficients \(c_{kl}\) by the call `call indat(u,x,y,t,hx,hy,nx,ny,m,nvar)`, with `indat` being specified in `acoustic2d_per.f90`. Here, as the initial data are composed of trigonometric functions we use the same approach as described in quickstart.

### Forming the Hermite interpolating polynomial

The first step in evolving the approximate solution to the next half timestep at a dual gridpoint \((x_{i+\frac{1}{2}},y_{j+\frac{1}{2}})\) is to form a degree \(2m+1\) tensor product polynomial

\[\begin{equation}
p_{i-\frac{1}{2},j-\frac{1}{2}}(x) = \sum_{k = 0}^{2m+1} \sum_{l = 0}^{2m+1}
c_{kl} \left( \frac{x - x_{i-\frac{1}{2} } }{h_x}\right)^k
\left( \frac{y - y_{j-\frac{1}{2} }}{h_y}\right)^l.
\end{equation}\]

from the polynomials representing the initial data at the four nearby primal gridpoints. That is, the coefficients of the new polynomial are determined by insisting that it interpolates the nearby polynomials. For each variable we form the Hermite interpolant at the dual gridpoint `(xh(i),yh(j))`, i.e. at the center of the cell with corners `(x(i-1),y(j))`, `(x(i),y(j))`, `(x(i-1),y(j-1))`, `(x(i),y(j-1))`. As the approximation is a tensor product polynomial we first interpolate on the top and bottom of the cell to the two face midpoints. Then we interpolate from the to the cell center. In the code this is done for one cell at a time by a simple call to `Taylor_to_Hermite2d`

``` {.sourceCode .fortran}
DO j=1,ny
  DO i=1,nx
    tcofs = 0.d0
    call Taylor_to_Hermite2d(uhermite,u(:,:,:,i-1:i,j-1:j),cintx,cinty,m,nvar)
    ! We store the Hermite interpolant in the array tcofs. The third index is for the
    ! time derivatives of the variables. The Hermite interpolant is the 0th time-
    ! derivative.
    do ivar = 1,nvar
    tcofs(0:2*m+1,0:2*m+1,0,ivar) = uhermite(:,:,ivar) 
    end do
....
```

### Evolving the Hermite interpolating polynomial

To evolve in time we first expand the solution in a Taylor series in time:

\[\begin{equation}
p^n_{i-\frac{1}{2},j-\frac{1}{2}}(x) = \sum_{k = 0}^{2m+1} \sum_{l =
0}^{2m+1} \sum_{s=0}^{q}
c_{kls} \left( \frac{x - x_{i-\frac{1}{2} } }{h_x}\right)^k
\left( \frac{y - y_{j-\frac{1}{2} }}{h_y}\right)^l \left( \frac{t - t_n }{\Delta t}\right)^s.
\end{equation}\]

The evolution in time requires that we know all coefficients \(c_{kls}\). In the code these are computed by a call to the routine `pde`. This routine fills in the time derivatives \(c_{kls}, \, s =1,\ldots,q\) and returns them in the array `tcofs(:,:,1:q,:)`.

``` {.sourceCode .fortran}
....
    ! The rest of the time-derivatives are filled in using the PDE.
    ! Here the PDE is autonomous so we do not need to pass in coordinates.  
    CALL pde(tcofs,2*m+1,q,nvar,cflx,cfly)
```

Once the time derivatives are known we compute the Taylor series at the next half timestep and store it in `uh`.

``` {.sourceCode .fortran}
....
    ! With all the time derivatives we finally add up the Taylor series
    ! for our degrees of freedom, i.e. derivatives up to degree m.
    uh(0:m,0:m,1:nvar,i,j) = tcofs(0:m,0:m,0,1:nvar)
    dfact=1.d0
    DO l=1,q
        dfact=dfact*.5d0
        uh(0:m,0:m,1:nvar,i,j) = uh(0:m,0:m,1:nvar,i,j) &
             +dfact*tcofs(0:m,0:m,l,1:nvar)
      END DO
  END DO
END DO
```

The `pde` routine residing in `acoustic2d_per.f90` has the same structure as in the one dimensional case

``` {.sourceCode .fortran}
SUBROUTINE pde(tcofs,mp,q,nvar,cflx,cfly)
 !
 ! time recursion coefficients as determined by the pde
 !
 IMPLICIT NONE

 INTEGER, INTENT(IN) :: nvar,mp,q
 DOUBLE PRECISION, DIMENSION(0:mp,0:mp,0:q,nvar), INTENT(INOUT) :: tcofs
 DOUBLE PRECISION, INTENT(IN) :: cflx,cfly
 INTEGER :: idx,idy,idt
 DOUBLE PRECISION, DIMENSION(0:mp,0:mp) :: ux,px,vy,py
 !
 DO idt=1,q
    ! Compute the "x-derivatives" 
    DO idx=0,mp-1
       px(idx,:) = cflx*DBLE(idx+1)*tcofs(idx+1,:,idt-1,1)
       ux(idx,:) = cflx*DBLE(idx+1)*tcofs(idx+1,:,idt-1,2)
    END DO
    px(mp,:) = 0.d0
    ux(mp,:) = 0.d0
    ! Compute the "y-derivatives" 
    DO idy=0,mp-1
       py(:,idy) = cfly*DBLE(idy+1)*tcofs(:,idy+1,idt-1,1)
       vy(:,idy) = cfly*DBLE(idy+1)*tcofs(:,idy+1,idt-1,3)
    END DO
    py(:,mp) = 0.d0
    vy(:,mp) = 0.d0
    ! 
    tcofs(:,:,idt,1) = -ux - vy
    tcofs(:,:,idt,2) = -px
    tcofs(:,:,idt,3) = -py
    tcofs(:,:,idt,:) = tcofs(:,:,idt,:)/DBLE(idt)
 END DO
```

> END SUBROUTINE pde

The main difference here is that the equations are coupled and the recursion relations needs to be evaluated one time derivative at a time. The derivation of the recursion is analogous to the one dimensional case (described in detail in quickstart).

### Second half timestep and boundary conditions

The second half timestep first enforces the periodic boundary conditions and then evolves the solution in the same way as above:

> 1.  Enforce periodic boundary conditions
>
>     ``` {.sourceCode .fortran}
>     ```
>
> > ! second half step t = t + 0.5d0\*dt
> >
> > ! We first enforce periodic boundary conditions by copying ! the periodic values to the "ghost points" uh(:,:,:,0,:) = uh(:,:,:,nx,:) uh(:,:,:,nx+1,:) = uh(:,:,:,1,:) uh(:,:,:,:,0) = uh(:,:,:,:,ny) uh(:,:,:,:,ny+1) = uh(:,:,:,:,1)
>
> 2. Form an interpolating polynomial
>
> > ``` {.sourceCode .fortran}
> > ! Then we proceed as in the first half step.
> > ```
>
> > DO j=0,ny  
> > DO i=0,nx  
> > tcofs=0.d0 ! Form the interpolant at a primal gridpoint ! (x(i),y(j)) call Taylor\_to\_Hermite2d(uhermite,uh(0:m,0:m,1:nvar,i:i+1,j:j+1),cintx,cinty,m,nvar) do ivar = 1,nvar tcofs(0:2\*m+1,0:2\*m+1,0,ivar) = uhermite(:,:,ivar) end do
> >
> > ...
>
> 1.  Then use the PDE and update the approximation at the full timestep
>
>     ``` {.sourceCode .fortran}
>     ```
>
> > ...
> >
> > > ! Compute time derivatives CALL pde(tcofs,2\*m+1,q,nvar,cflx,cfly) ! Evaluate the Taylor series and update the solution ! at a primal gridpoint. u(0:m,0:m,1:nvar,i,j)=tcofs(0:m,0:m,0,1:nvar) dfact=1.d0 DO l=1,q dfact=dfact\*.5d0 u(0:m,0:m,1:nvar,i,j)=u(0:m,0:m,1:nvar,i,j) & +dfact\*tcofs(0:m,0:m,l,1:nvar) END DO
> >
> > > END DO
> >
> > > END DO

### Compiling and Running the program

In the directory `chides/acoustic/2D` there is a simple `makefile` that can be used to compile the program by simply typing `make` in a terminal window. This assumies you are standing in the directory `chides/acoustic/2D` and that your system has `gfortran` installed. If you are using another compiler simply change the makefile.

If you are successful in building the solver you should have an executable named `periodic.x` that you can run.

### Displaying the results

After running the program you should see something like:

``` {.sourceCode .none}
bash-3.2$ ./periodic.x 
Computation finished:
-------------------------------------------------------------------------
m          dt          hx      maxerr     l2err-p     l2err-u     l2err-v
-------------------------------------------------------------------------
05  0.1667E+01  0.2094E+01  0.4164E-07  0.1215E-06  0.1471E-07  0.1471E-07
-------------------------------------------------------------------------
```

The program should also have generated a number of snapshots of the solution, `u0000000.txt`, `u0000001.txt`, etc. These can be displayed in Matlab using the script `plot_sol2D.m`.

### Changing the order and number of gridpoints

All the parameters associated with a simulation are contained in the top block of the `pdemodule` (residing in `acoustid2d_per.f90`), this is where you would change the order and the number of gridpoints etc.

**NOTE** that as the makefile structure is as simple as possible you should do `make clean; make` when you change the parameters. This ensures that all dependencies between modules are up to date.

``` {.sourceCode .fortran}
MODULE pdemodule
!
! This module solves the pde
! rho_t + u_x + v_y = 0,
! u_t + p_x = 0,
! v_t + v_y = 0,
!
! On a domain (x,y) \in [-\pi,pi]^2,
! with periodic boundary conditions.
!
! A solution is, p(x,y,t) = \sin(x) \sin(y) \sin(\omega t), with \omega = \sqrt(2),
!                u(x,y,t) = 1/\omega \cos(x) \sin(y) \cos(\omega t),
!                v(x,y,t) = 1/\omega \sin(x) \cos(y) \cos(\omega t).
!
!
IMPLICIT NONE
INTEGER, PARAMETER :: m = 5        ! m is the number of derivatives on each node
                                   ! the spatial order of the method is 2m+1
INTEGER, PARAMETER :: q = 4*m+2    ! terms in the Taylor series
!
DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932385d0
DOUBLE PRECISION, PARAMETER :: xl = -pi
DOUBLE PRECISION, PARAMETER :: xr =  pi
DOUBLE PRECISION, PARAMETER :: yb = -pi
DOUBLE PRECISION, PARAMETER :: yt =  pi
DOUBLE PRECISION, PARAMETER :: ttot = 10.d0

INTEGER, PARAMETER :: nx = 3      ! Discretization uses nx+1 points from xl to xr
INTEGER, PARAMETER :: ny = nx
DOUBLE PRECISION, PARAMETER :: CFL = 0.95d0 ! Ratio dt/h
! 
INTEGER, PARAMETER :: nref = 50   ! refinement factor for i-o and error computation.
!
INTEGER, PARAMETER :: nout = 1    ! output the solution every nout times
INTEGER, PARAMETER :: nvar = 3    ! Number of variables in solution
```

An example with boundary conditions
-----------------------------------

We next consider a case were we impose boundary conditions. This example is described in the files `acoustic2d_mirror.f90` and `herm2d_mirror.f90`. The domain is still a rectangle \((x,y) \in [x_L,x_R] \times [y_B,y_T]\) but now we impose zero Dirichlet boundary conditions for the density. Differentiating the boundary condition in time and in the tangential directions together with the governing equations yields the boundary conditions:

\[\begin{eqnarray}
  && \rho = 0, \ \ x={x_L,x_R}, y={y_L,y_T},\\
  && u = 0, \ \ v_y = 0, \ \ y={y_L,y_T}, \\
  && v = 0, \ \ u_x = 0, \ \ x={x_L,x_R}.
\end{eqnarray}\]

The boundary conditions are imposed at the half level by filling in the ghost point values by mirroring the polynomials so that they are either odd around the boundary and thus agree with the homogenous Dirichlet boundary conditions or so that they are even and satisfy homogenous Neumann conditions. The update of the boundary conditions is done in `herm2d_mirror.f90` by calling the subroutine below:

``` {.sourceCode .fortran}
SUBROUTINE update_bc(uh,nx,ny,m,nvar)
  implicit none
  INTEGER, INTENT(IN) :: nx,ny,m,nvar
  DOUBLE PRECISION, DIMENSION(0:m,0:m,nvar,0:nx+1,0:ny+1), INTENT(INOUT) :: uh
  DOUBLE PRECISION, DIMENSION(0:m) :: DIR,NEU
  INTEGER :: i,j,idx,idy
  do i = 0,m
     NEU(i) = (-1.d0)**i
     DIR(i) = (-1.d0)**(i+1)
  end do
  NEU(0) = 1.d0

  ! The pde is 
  ! (1) rho_t = -u_x - v_y,
  ! (2) u_t = -rho_x,
  ! (3) v_t = -rho_y.

  ! We implement rho = 0 on all boundaries.
  ! By (2) this implies u = 0 on top and bottom.
  ! By (1)+(2) we have v_y = 0 on top and bottom.
  ! Also:
  ! By (3) this implies v = 0 on left and right.
  ! By (1)+(3) we have u_x = 0 on left and right.

  ! Left and right side:
  do j = 1,ny
     do idy = 0,m
        uh(0:m,idy,1,0,j)    = DIR*uh(0:m,idy,1,1,j)
        uh(0:m,idy,1,nx+1,j) = DIR*uh(0:m,idy,1,nx,j)
        uh(0:m,idy,2,0,j)    = NEU*uh(0:m,idy,2,1,j)
        uh(0:m,idy,2,nx+1,j) = NEU*uh(0:m,idy,2,nx,j)
        uh(0:m,idy,3,0,j)    = DIR*uh(0:m,idy,3,1,j)
        uh(0:m,idy,3,nx+1,j) = DIR*uh(0:m,idy,3,nx,j)
     end do
  end do
  ! Top and bottom
  do i = 0,nx+1
     do idx = 0,m
        uh(idx,0:m,1,i,0)    = DIR*uh(idx,0:m,1,i,1)
        uh(idx,0:m,1,i,ny+1) = DIR*uh(idx,0:m,1,i,ny)
        uh(idx,0:m,2,i,0)    = DIR*uh(idx,0:m,2,i,1)
        uh(idx,0:m,2,i,ny+1) = DIR*uh(idx,0:m,2,i,ny)
        uh(idx,0:m,3,i,0)    = NEU*uh(idx,0:m,3,i,1)
        uh(idx,0:m,3,i,ny+1) = NEU*uh(idx,0:m,3,i,ny)
     end do
  end do
```

An example in a duct with forcing and Taylor series time stepping
-----------------------------------------------------------------

This example is contained in the files `herm2d_duct.f90` and `acoustic2d_duct.f90`.

``` {.sourceCode .fortran}
! This module solves the pde
! p_t = u_x + v_y + f,
! u_t = p_x,
! v_t = p_y,
!
! On a domain (x,y) \in [-10,10] X [0,1] 
!
! With wall boundary conditions on y=0,1 and periodic
! boundary conditions in x 
!
```
