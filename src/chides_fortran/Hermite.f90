SUBROUTINE Hermite(m,xl,xr,xc,tmat,icase)

 ! This subroutine computes the coefficient matrix tmat which
 ! transfers derivative data of order through m at xl and xr
 ! to the function and derivative data at xc
 ! That is, if p is a polynomial of degree 2m+1
 !
 !  h^k D^k p (xc)/k! = sum_(j=0)^m tmat(k,j) h^j D^j p(xl)/j!
 !                                + tmat(k,j+m+1) h^j D^j p(xr)/j!
 !
 !  icase < 0 => xc=xl (left boundary case)
 !  icase = 0 => xl < xc < xr
 !  icase > 0 => xc=xr (right boundary case)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: m,icase
 DOUBLE PRECISION, INTENT(IN) :: xl,xr,xc
 DOUBLE PRECISION, DIMENSION(0:2*m+1,0:2*m+1), INTENT(OUT) :: tmat
 DOUBLE PRECISION, DIMENSION(0:m+1,0:m+1) :: bcofs
 DOUBLE PRECISION :: h,z,zc,adl,adr,sign,c1l,c1r,c2l,c2r
 INTEGER :: i,j,k

 ! Compute in normalized coordinates

 h=xr-xl
 z=(xc-xl)/h
 zc=z-1.d0

 IF (icase > 0) THEN
  z=1.d0
  zc=0.d0
 ELSE IF (icase < 0) THEN
  z=0.d0
  zc=-1.d0
 END IF

 call binomial(bcofs,m+1)

 !
 ! We begin with the Hermite-Lagrange interpolants:
 !
 !   Q_j (z) = z^(m+1) sum_{k=j}^m h_{kj} (z-1)^k,
 !
 !   j=0, ... , m
 !
 !   satisfying Q_j = (z-1)^j + O((z-1)^(m+1)) at z=1
 !
 !   After some algebra once can show:
 !
 !   h_jj = 1,  h_kj = -sum_{p=j}^{k-1} b_(k-p)^(m+1) h_pj ,
 !              for k>j
 !
 !   here b_(k-p)^(m+1) is the binomial coefficient (m+1)!/((k-p)!(m+1-k+p)!)
 !
 ! To construct the matrix we
 ! now evaluate the interpolant and its derivatives at z
 !
 ! Data to the left is handled by a simple change of variables z -> 1-z
 !
 ! Start with the last column - note that we directly use the recursive
 ! definition of the h's to fold previously computed columns into old
 ! ones. Note that the polynomial is now centered about the midpoint
 !

 DO i=0,2*m+1
  adl=0.d0
  adr=0.d0
  DO j=MAX(0,i-m),MIN(i,m+1)
   if ((m-i+j).eq. 0 ) then
    c2l = 1.d0
    c2r = 1.d0
   else
    c2l = z**(m-i+j)
    c2r = zc**(m-i+j)
   end if
   if ((m+1-j).eq.0) then
    c1l = 1.d0
    c1r = 1.d0
   else
    c1l = (zc**(m+1-j))
    c1r = (z**(m+1-j))
   end if
   adr=adr+bcofs(m+1,j)*bcofs(m,i-j)*c1r*c2r
   adl=adl+bcofs(m+1,j)*bcofs(m,i-j)*c1l*c2l
  END DO
  tmat(i,2*m+1)=adr
  tmat(i,m)=((-1.d0)**(m+1))*adl
 END DO

 ! Now loop over the other columns backwards

 DO k=m-1,0,-1
  DO i=0,2*m+1
   adl=0.d0
   adr=0.d0
   DO j=MAX(0,i-k),MIN(i,m+1)
    if ((k-i+j).eq. 0 ) then
     c2l = 1.d0
     c2r = 1.d0
    else
     c2l = z**(k-i+j)
     c2r = zc**(k-i+j)
    end if
    if ((m+1-j).eq.0) then
     c1l = 1.d0
     c1r = 1.d0
    else
     c1l = (zc**(m+1-j))
     c1r = (z**(m+1-j))
    end if
    adr=adr+bcofs(m+1,j)*bcofs(k,i-j)*c1r*c2r
    adl=adl+bcofs(m+1,j)*bcofs(k,i-j)*c1l*c2l
   END DO
   tmat(i,k+m+1)=adr
   tmat(i,k)=((-1.d0)**(m+1))*adl
   sign=1.d0
   DO j=k+1,m
    sign=-sign
    tmat(i,k)=tmat(i,k)-sign*bcofs(m+1,j-k)*tmat(i,j)
    tmat(i,k+m+1)=tmat(i,k+m+1)-bcofs(m+1,j-k)*tmat(i,j+m+1)
   END DO
  END DO
 END DO

END SUBROUTINE Hermite

SUBROUTINE binomial(coeffs,m)
 ! Computes the binomial coefficients of order up through m
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: m
 DOUBLE PRECISION, DIMENSION(0:m,0:m), INTENT(OUT) :: coeffs
 INTEGER :: i,j

 coeffs(0,0)=1.d0

 DO i=1,m
  coeffs(i,0)=1.d0
  DO j=1,i
   coeffs(i,j)=coeffs(i,j-1)*DBLE(i-j+1)/DBLE(j)
  END DO
 END DO

END SUBROUTINE binomial

SUBROUTINE point_to_Taylor(q,y,p)
 !
 ! Takes data on q+1 Chebyshev nodes and
 ! uses the Bjorck-Peyreyra algorithm to express
 ! the interpolant as a degree q Taylor series
 ! (see Dahlquist-Bjorck volume 1 page 375) 
 !
 INTEGER, INTENT(IN) :: q
 DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: y
 DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: p 
!
 DOUBLE PRECISION, DIMENSION(0:q) :: z
 DOUBLE PRECISION :: pi=3.1415926535897932385d0
 INTEGER :: j,k
!
 DO j=0,q
   z(j)=-.5d0*COS(pi*DBLE(j)/DBLE(q))
 END DO
!
! First divided differences
!
 p=y
 DO k=1,q
   DO j=q,k,-1
     p(j)=(p(j)-p(j-1))/(z(j)-z(j-k))
   END DO
 END DO
!
! Now Horner
!
  DO k=q-1,0,-1
    DO j=k,q-1
      p(j)=p(j)-z(k)*p(j+1)
    END DO
  END DO
!
END SUBROUTINE point_to_Taylor 

