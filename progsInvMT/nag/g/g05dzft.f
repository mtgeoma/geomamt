      LOGICAL FUNCTION G05DZF(P)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A LOGICAL RESULT, TRUE WITH PROBABILITY P.
C     .. Scalar Arguments ..
      DOUBLE PRECISION        P
C     .. Local Scalars ..
      DOUBLE PRECISION        X
C     .. External Functions ..
      DOUBLE PRECISION        G05CAF
      EXTERNAL                G05CAF
C     .. Executable Statements ..
      G05DZF = P .GT. G05CAF(X)
      RETURN
      END
