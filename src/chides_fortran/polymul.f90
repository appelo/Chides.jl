SUBROUTINE tpolymul(a,b,c,q)
  !
  !  truncated multiplication of degree q polynomials
  !  with coefficients a,b. On return c will contain the
  !  coefficients of the product truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  INTEGER :: j,k
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: a,b
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: c
  !
  c=0.d0
  DO j=0,q
     DO k=0,j
        c(j)=c(j)+a(k)*b(j-k)
     END DO
  END DO

END SUBROUTINE tpolymul

SUBROUTINE tpolymul2(a,b,c,q)
  !
  !  truncated multiplication of tensor product
  !  degree q polynomials in 2 variables
  !  with coefficients a,b. On return c will contain the
  !  coefficients of the product truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(IN) :: a,b
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(OUT) :: c
  INTEGER :: j1,j2,k1,k2
  !
  c=0.d0
  DO j2=0,q
     DO j1=0,q
        DO k2=0,j2
           DO k1=0,j1
              c(j1,j2)=c(j1,j2)+a(k1,k2)*b(j1-k1,j2-k2)
           END DO
        END DO
     END DO
  END DO
END SUBROUTINE tpolymul2

SUBROUTINE tpolymul3(a,b,c,q)
  !
  !  truncated multiplication of tensor product
  !  degree q polynomials in 3 variables
  !  with coefficients a,b. On return c will contain the
  !  coefficients of the product truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q,0:q), INTENT(IN) :: a,b
  DOUBLE PRECISION, DIMENSION(0:q,0:q,0:q), INTENT(OUT) :: c
  INTEGER :: j1,j2,j3,k1,k2,k3
  !
  c=0.d0
  DO j3=0,q
     DO j2=0,q
        DO j1=0,q
           DO k3=0,j3
              DO k2=0,j2
                 DO k1=0,j1
                    c(j1,j2,j3)=c(j1,j2,j3)+a(k1,k2,k3)*b(j1-k1,j2-k2,j3-k3)
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
END SUBROUTINE tpolymul3

