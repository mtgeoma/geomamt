      DOUBLE PRECISION FUNCTION E04DGW(OBJF,FGUESS,GTP,SMAX)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     E04DGW calculates an initial step length. Refer
C     to the formula on page 10 of the CG paper.
C     FGUESS here is F(est) there, SMAX is a maximum step
C     length.
C
C     -- Written on 4-June-1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, RDUMMY
      PARAMETER                        (ONE=1.0D+0,RDUMMY=-11111.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 FGUESS, GTP, OBJF, SMAX
C     .. Arrays in Common ..
      DOUBLE PRECISION                 WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION                 D, EPSMCH
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Common blocks ..
      COMMON                           /AX02ZA/WMACH
C     .. Save statement ..
      SAVE                             /AX02ZA/
C     .. Executable Statements ..
      EPSMCH = WMACH(3)
      IF (FGUESS.EQ.RDUMMY) THEN
         D = 0.0D0
      ELSE
         D = ABS(OBJF-FGUESS)
      END IF
      E04DGW = ONE
      IF ((2*D.LE.(-GTP)) .AND. (D.GE.EPSMCH)) E04DGW = -2*D/GTP
      IF (E04DGW.GE.SMAX) E04DGW = SMAX
      RETURN
C
C     End of E04DGW (STEP1).
C
      END
