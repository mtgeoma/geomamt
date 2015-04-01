      DOUBLE PRECISION FUNCTION F04JGW(N,X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (V2NORM)
C
C     REAL FUNCTION F04JGW RETURNS THE VALUE OF THE EUCLIDEAN
C     LENGTH OF THE N ELEMENT VECTOR X. F04JGW IS DEFINED AS
C
C     F04JGW = SQRT(X(1)**2+X(2)**2+...+X(N)**2).
C
C     .. Scalar Arguments ..
      INTEGER                          N
C     .. Array Arguments ..
      DOUBLE PRECISION                 X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 BIG, SMALL, TINY
C     .. External Functions ..
      DOUBLE PRECISION                 F04JGV, X02AMF
      EXTERNAL                         F04JGV, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      BIG = 1.0D0/SMALL
      TINY = SQRT(SMALL)
C
      F04JGW = F04JGV(N,X,TINY,BIG)
C
      RETURN
      END
