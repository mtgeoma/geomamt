      DOUBLE PRECISION FUNCTION G05DBF(A)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER, NEGATIVE EXPONENTIALLY
C     DISTRIBUTED WITH MEAN A.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A
C     .. Local Scalars ..
      DOUBLE PRECISION                 X
C     .. External Functions ..
      DOUBLE PRECISION                 G05CAF
      EXTERNAL                         G05CAF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, LOG
C     .. Executable Statements ..
      G05DBF = -ABS(A)*LOG(G05CAF(X))
      RETURN
      END
