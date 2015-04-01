      SUBROUTINE D02TKX(RHO,COEF,KCOL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C**********************************************************************
C
C   Purpose:
C      solve Vandermonde system v * x = e
C      with  v(i,j) = rho(j)**(i-1)/(i-1)!,   i,j=1,2,...,k
C      and e = [0,...,0,1]**T.
C
C   Arguments:
C      RHO    - contains the coefficients of the Vandermonde matrix
C      COEF   - returns the solution vector
C      KCOL   - the dimension of the system (the number of collocation
C               points
C
C    Author:
C       R.W. Brankin, NAG Ltd., August 1994
C       (straight from COLNEW routine VMONDE)
C
C    Comment:
C       Vandermone matrices become increasingly ill-conditioned the
C       larger the dimension. However, the max size considered here
C       is 7. Coefficients are positive and stored in increasing
C       magnitude which provides for accurate computation (despite
C       ill-condition, see Higham NJ, Error analysis of the Bjork-Peryra
C       algorithms for solving van der Monde systems, Numer Math 50
C       (1987), 613-632)
C
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER           KCOL
C     .. Array Arguments ..
      DOUBLE PRECISION  COEF(KCOL), RHO(KCOL)
C     .. Local Scalars ..
      INTEGER           I, IFAC, J, KM1, KMI
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IF (KCOL.GT.1) THEN
         KM1 = KCOL - 1
         DO 40 I = 1, KM1
            KMI = KCOL - I
            DO 20 J = 1, KMI
               COEF(J) = (COEF(J+1)-COEF(J))/(RHO(J+I)-RHO(J))
   20       CONTINUE
   40    CONTINUE
C
         IFAC = 1
         DO 80 I = 1, KM1
            KMI = KCOL + 1 - I
            DO 60 J = 2, KMI
               COEF(J) = COEF(J) - RHO(J+I-1)*COEF(J-1)
   60       CONTINUE
            COEF(KMI) = DBLE(IFAC)*COEF(KMI)
            IFAC = IFAC*I
   80    CONTINUE
         COEF(1) = DBLE(IFAC)*COEF(1)
      END IF
      RETURN
      END
