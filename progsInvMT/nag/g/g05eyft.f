      INTEGER FUNCTION G05EYF(R,NR)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS AN INTEGER FROM THE REFERENCE VECTOR IN R.
C     .. Scalar Arguments ..
      INTEGER                 NR
C     .. Array Arguments ..
      DOUBLE PRECISION        R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION        X
      INTEGER                 N
C     .. External Functions ..
      DOUBLE PRECISION        G05CAF
      EXTERNAL                G05CAF
C     .. Intrinsic Functions ..
      INTRINSIC               INT
C     .. Executable Statements ..
      X = G05CAF(X)
      N = INT(X*R(1))
      N = INT(R(N+3))
      IF (X.LE.R(N)) GO TO 40
   20 N = N + 1
      IF (X.GT.R(N)) GO TO 20
   40 G05EYF = N + INT(R(2))
      RETURN
      END
