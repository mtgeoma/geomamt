      DOUBLE PRECISION FUNCTION D01AHU(A)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-820 (DEC 1989).
C     USED TO AVOID ZERO DIVISIONS.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A
C     .. Scalars in Common ..
      DOUBLE PRECISION                 AFLOW, EPMACH, UFLOW
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SIGN
C     .. Common blocks ..
      COMMON                           /DD01AH/EPMACH, UFLOW, AFLOW
C     .. Executable Statements ..
      D01AHU = 10.0D0*UFLOW
      IF (A.EQ.0.0D0) RETURN
      D01AHU = SIGN(MAX(ABS(A),10.0D0*UFLOW),A)
      RETURN
      END
