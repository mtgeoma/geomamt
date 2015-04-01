      SUBROUTINE Y90RQY(N,Q,LDQ)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =====================================================
C         *  Y90RQF :  Utility for Toeplitz Matrix Generator  *
C         =====================================================
C
C
C     -- Written on 10-October-1990.
C     Sven Hammarling, Nag Ltd.
C
C
C     Purpose
C     =======
C
C     Y90RQY  forms  an  initial orthogonal matrix  Q  for the routine
C     Y90RQF.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D+0)
C     .. Scalar Arguments ..
      INTEGER           LDQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST1, CONST2, CONSTJ
      INTEGER           I, J
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         SIN, SQRT
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Execution
C
C-----------------------------------------------------------------------
      CONST1 = X01AAF(CONST1)/(N+1)
      CONST2 = SQRT(TWO/(N+1))
      DO 40 J = 1, N
         CONSTJ = J*CONST1
         DO 20 I = 1, N
            Q(I,J) = CONST2*SIN(I*CONSTJ)
   20    CONTINUE
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RQY.
C
C-----------------------------------------------------------------------
      RETURN
      END
