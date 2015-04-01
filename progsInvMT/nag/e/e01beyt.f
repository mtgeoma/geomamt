      DOUBLE PRECISION FUNCTION E01BEY(ARG1,ARG2)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C         PCHST:  Sign-testing Routine.
C
C     Returns:
C        -1. if ARG1 And ARG2 are of opposite sign.
C         0. if either argument is zero.
C        +1. if ARG1 and ARG2 are of the same sign.
C
C     The object is to do this without multiplying ARG1*ARG2, to avoid
C     possible over/underflow problems.
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     ------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, ONE
      PARAMETER                        (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 ARG1, ARG2
C     .. Intrinsic Functions ..
      INTRINSIC                        SIGN
C     .. Executable Statements ..
C     Perform the test.
      E01BEY = SIGN(ONE,ARG1)*SIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO)) E01BEY = ZERO
      RETURN
C
      END
