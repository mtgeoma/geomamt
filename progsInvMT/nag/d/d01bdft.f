      SUBROUTINE D01BDF(F,A,B,EPSABS,EPSREL,RESULT,ABSERR)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QNG.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, B, EPSABS, EPSREL, RESULT
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      INTEGER           IER, NEVAL
C     .. External Subroutines ..
      EXTERNAL          D01BDV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      CALL D01BDV(F,A,B,ABS(EPSABS),ABS(EPSREL),RESULT,ABSERR,NEVAL,IER)
C
      RETURN
      END
