      DOUBLE PRECISION FUNCTION G01BDZ(A,B)
C     MARK 7 RELEASE, NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CALLED BY NAG ROUTINE G01BDF
C
C     EVALUATES NATURAL LOGARITHM OF COMPLETE
C     BETA FUNCTION B(A,B), USING S14AAF IF
C     A,B,A+B ARE SMALL, AND APPROXIMATION
C     6.1.41 ON P.257 OF ABRAMOWITZ AND STEGUN
C     IF ARGUMENTS ARE LARGE.
C     ON ENTRY,A AND B ARE POSITIVE.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
C     .. Local Scalars ..
      DOUBLE PRECISION                 G, GAL, GBL, TWOPIL, X, Z
      INTEGER                          IFA, N
C     .. External Functions ..
      DOUBLE PRECISION                 S14AAF, X01AAF
      EXTERNAL                         S14AAF, X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      TWOPIL = LOG(2.0D0*X01AAF(X))
      IFA = 1
      G = S14AAF(A,IFA)
      IF (IFA.EQ.0) GO TO 20
      N = 1
      X = A
      GO TO 100
   20 G = LOG(G)
   40 GAL = G
      IFA = 1
      G = S14AAF(B,IFA)
      IF (IFA.EQ.0) GO TO 60
      N = 2
      X = B
      GO TO 100
   60 G = LOG(G)
   80 GBL = G
      IFA = 1
      G = S14AAF(A+B,IFA)
      IF (IFA.EQ.0) GO TO 120
      N = 3
      X = A + B
  100 Z = 1.0D0/(X*X)
      G = (X-0.5D0)*LOG(X) - X + 0.5D0*TWOPIL + (((-3.0D0*Z+4.0D0)
     *    *Z-14.0D0)*Z+420.0D0)/(5040.0D0*X)
      GO TO (40,80,140) N
  120 G = LOG(G)
  140 G01BDZ = GAL + GBL - G
      RETURN
      END
