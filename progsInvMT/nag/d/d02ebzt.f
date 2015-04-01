      SUBROUTINE D02EBZ(OUTPUT,N,YSAV,ACOR,CHK,ISAVE,HU,H,X,XOUT,IMON,
     *                  DIR,RNQU)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DIR, H, HU, RNQU, X, XOUT
      INTEGER           IMON, ISAVE, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(N), CHK(N), YSAV(N,6)
C     .. Subroutine Arguments ..
      EXTERNAL          OUTPUT
C     .. Local Scalars ..
      DOUBLE PRECISION  XSAVE
      INTEGER           IFLAG, NQU
      LOGICAL           INTERP
C     .. External Subroutines ..
      EXTERNAL          D02XKF
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Executable Statements ..
      IF (IMON.NE.1) RETURN
   20 INTERP = DIR*XOUT .LE. DIR*X
      IF ( .NOT. INTERP) RETURN
      NQU = INT(RNQU)
      IFLAG = 1
      CALL D02XKF(XOUT,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
      IF (IFLAG.NE.0) THEN
         IMON = -2
         ISAVE = 7
         RETURN
      END IF
      XSAVE = XOUT
      CALL OUTPUT(XOUT,CHK)
      IF (DIR*XOUT.LT.DIR*XSAVE) THEN
         IMON = -2
         ISAVE = 5
      ELSE
         GO TO 20
      END IF
      RETURN
      END
