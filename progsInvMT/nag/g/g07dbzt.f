      SUBROUTINE G07DBZ(C,BETA)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0,TWO=2.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA, C
C     .. Local Scalars ..
      DOUBLE PRECISION  C2, DC, EXMIN, PC
      INTEGER           IFAIL
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X01AAF, X02AMF
      EXTERNAL          S15ABF, X01AAF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SQRT
C     .. Executable Statements ..
C     Compute BETA as a function of C
C
      C2 = -C*C/TWO
      IFAIL = 0
      PC = S15ABF(C,IFAIL)
      EXMIN = LOG(X02AMF())
      IF (C2.GE.EXMIN) THEN
         DC = EXP(C2)/SQRT(2.0D0*X01AAF(0.0D0))
      ELSE
         DC = ZERO
      END IF
      BETA = -C*DC + PC - HALF - C2*(ONE-PC)*TWO
      RETURN
      END
