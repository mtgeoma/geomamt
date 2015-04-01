      SUBROUTINE D02EJU(XOUT,OUTPUT,CHK,N,YSAV,ACOR,X,NQU,HU,H,IFLAG,
     *                  XLAST,XNEW,DIR)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DIR, H, HU, X, XLAST, XNEW, XOUT
      INTEGER           IFLAG, N, NQU
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(N), CHK(N), YSAV(N,6)
C     .. Subroutine Arguments ..
      EXTERNAL          OUTPUT
C     .. Local Scalars ..
      DOUBLE PRECISION  XSAVE
      LOGICAL           INTERP
C     .. External Subroutines ..
      EXTERNAL          D02XKF
C     .. Executable Statements ..
   20 INTERP = DIR*XLAST .LE. DIR*XOUT .AND. DIR*XOUT .LE. DIR*XNEW
      IF ( .NOT. INTERP) RETURN
      IFLAG = 1
      CALL D02XKF(XOUT,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
      IF (IFLAG.NE.0) RETURN
      XSAVE = XOUT
      CALL OUTPUT(XOUT,CHK)
      IF (DIR*XOUT.LT.DIR*XSAVE) THEN
         IFLAG = 5
         RETURN
      ELSE
         GO TO 20
      END IF
      END
