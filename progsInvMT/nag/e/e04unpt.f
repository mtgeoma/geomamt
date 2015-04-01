      SUBROUTINE E04UNP(MODE,N,M,Y,F,FJAC,LDFJ,OBJF,GRAD,YF)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     ==================================================================
C     E04UNP  loads the objective and gradient of the function
C                      (1/2)(y - f(x))'(y - f(x)).
C
C     Original version written 12-Jul-94.
C     This version of  E04UNP  dated 12-Jul-94.
C     ==================================================================
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJF
      INTEGER           LDFJ, M, MODE, N
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), FJAC(LDFJ,*), GRAD(N), Y(M), YF(M)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV
C     .. Executable Statements ..
C
      CALL DCOPY(M,Y,1,YF,1)
      CALL DAXPY(M,(-ONE),F,1,YF,1)
      IF (MODE.NE.1) THEN
         OBJF = HALF*(DNRM2(M,YF,1))**2
      END IF
C
      IF (MODE.NE.0) THEN
         CALL DGEMV('Transpose',M,N,(-ONE),FJAC,LDFJ,YF,1,ZERO,GRAD,1)
      END IF
C
C     End of E04UNP. (NLOBJF)
C
      RETURN
      END
