      DOUBLE PRECISION FUNCTION G05DEF(A,B)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER, LOGNORMALLY DISTRIBUTED WITH
C     MEAN EXP(A+B*B/2) AND VARIANCE EXP(2*A+B*B)*(EXP(B*B)-1).
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
C     .. External Functions ..
      DOUBLE PRECISION                 G05DDF
      EXTERNAL                         G05DDF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP
C     .. Executable Statements ..
      G05DEF = EXP(G05DDF(A,B))
      RETURN
      END
