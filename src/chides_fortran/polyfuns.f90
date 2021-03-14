SUBROUTINE tpolydiv(a,b,c,q)
  !
  !  truncated division of degree q polynomials
  !  with coefficients a,b. On return c will contain the
  !  coefficients of a/b truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  INTEGER :: j,k
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: a,b
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: c
  !
  c=0.d0
  c(0)=a(0)/b(0)
  DO j=1,q
     c(j)=a(j)
     DO k=1,j
        c(j)=c(j)-c(j-k)*b(k)
     END DO
     c(j)=c(j)/b(0)
  END DO

END SUBROUTINE tpolydiv

SUBROUTINE tpolypow(a,c,r,q)
  !
  !  truncated power function of degree q polynomial
  !  with coefficients a. On return c will contain the
  !  coefficients of a**r truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  INTEGER :: j,k
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: c
  DOUBLE PRECISION, INTENT(IN) :: r
  !
  c=0.d0
  c(0)=a(0)**r
  DO j=1,q
     DO k=1,j
        c(j)=c(j)+a(k)*c(j-k)*(DBLE(k)*(1.d0+r)-DBLE(j))
     END DO
     c(j)=c(j)/(DBLE(j)*a(0))
  END DO

END SUBROUTINE tpolypow

SUBROUTINE tpolypow2(a,c,r,q)
  !
  !  truncated power function of a tensor product degree
  !  (q,q) polynomial
  !  with coefficients a. On return c will contain the
  !  coefficients of a**r truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(OUT) :: c
  DOUBLE PRECISION, INTENT(IN) :: r
  INTEGER :: j1,j2,k1,k2
  !
  c=0.d0
  c(0,0)=a(0,0)**r
  DO k1=1,q
     DO k2=1,k1
        c(0,k1)=c(0,k1)+a(0,k2)*c(0,k1-k2)*(DBLE(k2)*(1.d0+r)-DBLE(k1))
     END DO
     c(0,k1)=c(0,k1)/(DBLE(k1)*a(0,0))
  END DO
  DO j1=1,q
     DO j2=1,j1
        c(j1,0)=c(j1,0)+a(j2,0)*c(j1-j2,0)*(DBLE(j2)*(1.d0+r)-DBLE(j1))
     END DO
     c(j1,0)=c(j1,0)/(DBLE(j1)*a(0,0))
     DO k1=1,q
        DO j2=1,j1
           DO k2=0,k1
              c(j1,k1)=c(j1,k1)+   &
                   a(j2,k2)*c(j1-j2,k1-k2)*(DBLE(j2)*(1.d0+r)-DBLE(j1))
           END DO
        END DO
        DO k2=1,k1
           c(j1,k1)=c(j1,k1)-DBLE(j1)*a(0,k2)*c(j1,k1-k2)
        END DO
        c(j1,k1)=c(j1,k1)/(DBLE(j1)*a(0,0))
     END DO
  END DO

END SUBROUTINE tpolypow2

SUBROUTINE tpolyexp(a,c,q)
  !
  !  truncated exponentiation of degree q polynomial
  !  with coefficients a. On return c will contain the
  !  coefficients of exp(a) truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  INTEGER :: j,k
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: c
  !
  c=0.d0
  c(0)=EXP(a(0))
  DO j=1,q
     DO k=1,j
        c(j)=c(j)+DBLE(k)*a(k)*c(j-k)
     END DO
     c(j)=c(j)/DBLE(j)
  END DO

END SUBROUTINE tpolyexp

SUBROUTINE tpolyexp2(a,c,q)
  !
  !  truncated exponential function of a tensor product degree
  !  (q,q) polynomial
  !  with coefficients a. On return c will contain the
  !  coefficients of exp(a) truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(OUT) :: c
  INTEGER :: j1,j2,k1,k2
  !
  c=0.d0
  c(0,0)=EXP(a(0,0))
  DO k1=1,q
     DO k2=1,k1
        c(0,k1)=c(0,k1)+a(0,k2)*c(0,k1-k2)*DBLE(k2)
     END DO
     c(0,k1)=c(0,k1)/DBLE(k1)
  END DO
  DO j1=1,q
     DO j2=1,j1
        c(j1,0)=c(j1,0)+a(j2,0)*c(j1-j2,0)*DBLE(j2)
     END DO
     c(j1,0)=c(j1,0)/DBLE(j1)
     DO k1=1,q
        DO j2=1,j1
           DO k2=0,k1
              c(j1,k1)=c(j1,k1)+a(j2,k2)*c(j1-j2,k1-k2)*DBLE(j2)
           END DO
        END DO
        c(j1,k1)=c(j1,k1)/DBLE(j1)
     END DO
  END DO

END SUBROUTINE tpolyexp2

SUBROUTINE tpolyexp3(a,c,q)
  !
  !  truncated exponential function of a tensor product degree
  !  (q,q,q) polynomial
  !  with coefficients a. On return c will contain the
  !  coefficients of exp(a) truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q,0:q), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:q,0:q,0:q), INTENT(OUT) :: c
  DOUBLE PRECISION, DIMENSION(0:q,0:q) :: s2
  INTEGER :: j1,j2,k1,k2
  !
  c=0.d0
  c(0,0,0)=EXP(a(0,0,0))
  DO k1=1,q
     DO k2=1,k1
        c(0,0,k1)=c(0,0,k1)+a(0,0,k2)*c(0,0,k1-k2)*DBLE(k2)
     END DO
     c(0,0,k1)=c(0,0,k1)/DBLE(k1)
  END DO
  DO j1=1,q
     DO j2=1,j1
        c(0,j1,0)=c(0,j1,0)+a(0,j2,0)*c(0,j1-j2,0)*DBLE(j2)
     END DO
     c(0,j1,0)=c(0,j1,0)/DBLE(j1)
     DO k1=1,q
        DO j2=1,j1
           DO k2=0,k1
              c(0,j1,k1)=c(0,j1,k1)+a(0,j2,k2)*c(0,j1-j2,k1-k2)*DBLE(j2)
           END DO
        END DO
        c(0,j1,k1)=c(0,j1,k1)/DBLE(j1)
     END DO
  END DO
  !
  DO j1=1,q
     DO j2=1,j1
        CALL tpolymul2(a(j2,:,:),c(j1-j2,:,:),s2,q)
        c(j1,:,:)=c(j1,:,:)+DBLE(j2)*s2
     END DO
     c(j1,:,:)=c(j1,:,:)/DBLE(j1)
  END DO

END SUBROUTINE tpolyexp3

SUBROUTINE tpolyinv3(a,b,q)
  !
  !  multiplicative inverse of a
  !  degree q polynomials in 3 variables
  !  with coefficients a. On return b will contain the
  !  coefficients of "1/a" truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q,0:q), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:q,0:q,0:q), INTENT(OUT) :: b
  INTEGER :: j1,j2,j3,k1,k2,k3
  DOUBLE PRECISION :: ab
  !
  b = 0.d0
  b(0,0,0) = 1.d0/a(0,0,0)
  DO j3=0,q
     DO j2=0,q
        DO j1=0,q
           ab = 0.d0
           b(0,0,0) = 1.d0/a(0,0,0)
           DO k3=0,j3
              DO k2=0,j2
                 DO k1=0,j1
                    ab = ab + a(k1,k2,k3)*b(j1-k1,j2-k2,j3-k3)
                 END DO
              END DO
           END DO
           b(j1,j2,j3) = -ab/a(0,0,0)
           ! if( (j1.ne.0).and.(j2.ne.0).and.(j3.ne.0) )
        END DO
     END DO
  END DO

END SUBROUTINE tpolyinv3

SUBROUTINE extrapolate(p,z,q)
  !
  ! recenters a polynomial from 0 to z
  !
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(INOUT) :: p
  DOUBLE PRECISION, INTENT(IN) :: z
  DOUBLE PRECISION, DIMENSION(0:q) :: ploc
  INTEGER :: j,k
  !
  ploc=0.d0
  ploc(0)=p(q)
  DO j=q-1,0,-1
     DO k=q-j,1,-1
        ploc(k)=z*ploc(k)+ploc(k-1)
     END DO
     ploc(0)=z*ploc(0)+p(j)
  END DO
  !
  p=ploc

END SUBROUTINE extrapolate
