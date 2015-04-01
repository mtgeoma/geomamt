      DOUBLE PRECISION FUNCTION G05DCF(A,B)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER, LOGISTICALLY DISTRIBUTED WITH
C     MEAN A AND STANDARD DEVIATION B*PI/SQRT(3).
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
C     .. Local Scalars ..
      DOUBLE PRECISION                 ONE, X
C     .. External Functions ..
      DOUBLE PRECISION                 G05CAF
      EXTERNAL                         G05CAF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Data statements ..
      DATA                             ONE/1.0D0/
C     .. Executable Statements ..
      X = G05CAF(X)
      G05DCF = A + B*LOG(X/(ONE-X))
      RETURN
      END
