Service Routines
================

Hermite interpolation
---------------------

As was already mentioned in the quickstart example, Hermite methods are based on the representation of the solution \(u(x,t)\) in terms of degree \(m+1\) nodal based scaled polynomials centered around each grid point. For example on the primal grid we have that

\[\begin{equation}
u(x,t) \approx p_i(x) = \sum_{l = 0}^{m} c_l \left( \frac{x - x_i }{h_x}\right)^l, x \in [x_{i-\frac{1}{2}},x_{i+\frac{1}{2}}], \ \ i = 0,\ldots,n_x.
\end{equation}\]

### Forming the Hermite interpolating polynomial in one dimension

Forming the Hermite interpolant at \(x_{i+\frac{1}{2}}\) from the polynomials representing the initial data at the two nearby primal gridpoints can be done by first assembling the data in a vector and then multiplying that vector with a matrix. That is, the coefficients of the new polynomial are determined by insisting that it interpolates the nearby polynomials. Functionally, this amounts to applying a the linear map corresponding to Hermite interpolation, say \(A\), to the nearby polynomials to find \(p_{i-\frac{1}{2}} =
A(p_{i-1},p_{i})\).

The matrix \(A\) is obtained by a call

``` {.sourceCode .fortran}
CALL Hermite(m,node_left,node_right,node_interp,A,0)
```

using the subroutine `Hermite` from [chides/service/Hermite.f90](https://bitbucket.org/appelo/chides). As is perhaps obvious from the call to `Hermite` it takes four arguments, the location of the nodes, the location where the Hermite interpolant is centered and an integer which is 0 if the center point is distinct from the nodes, -1 if it coincides with the left node and 1 if it coincides with the right.

Aside from the routine `Hermite` there are also three routines `Taylor_to_HermiteXd`, with `X = 1,2,3`. These return the Hermite interpolant at the cell center given nodal data at the vertexes. See also the acoustics examples in 1D \<acoustic1d\> and 2D \<acoustic2d\> for examples of how they are used.

The Polyfuns libraries
----------------------

### Modal multiplication - direct convolution

### Modal multiplication - FFT

### Pointwise multiplication
