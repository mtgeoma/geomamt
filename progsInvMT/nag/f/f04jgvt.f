      DOUBLE PRECISION FUNCTION F04JGV(N,X,TINY,BIG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (VENORM)
C
C     REAL FUNCTION F04JGV RETURNS THE VALUE OF THE EUCLIDEAN
C     LENGTH OF THE N ELEMENT VECTOR X. F04JGV IS DEFINED AS
C
C     F04JGV = SQRT(X(1)**2+X(2)**2+...+X(N)**2).
C
C     TINY AND BIG MUST BE GIVEN BY
C
C     TINY = SQRT(X02AMF)   AND   BIG = 1.0/X02AMF
C
C     WHERE X02AMF IS THE SMALL REAL VALUE RETURNED FROM
C     ROUTINE X02AMF.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 BIG, TINY
      INTEGER                          N
C     .. Array Arguments ..
      DOUBLE PRECISION                 X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUMSQ
C     .. External Functions ..
      DOUBLE PRECISION                 F04JGU
      LOGICAL                          X02DAF
      EXTERNAL                         F04JGU, X02DAF
C     .. External Subroutines ..
      EXTERNAL                         F04JGT
C     .. Executable Statements ..
      SCALE = 0.0D0
      SUMSQ = 1.0D0
C
      CALL F04JGT(N,X,SCALE,SUMSQ,TINY,X02DAF(0.0D0))
C
      F04JGV = F04JGU(SCALE,SUMSQ,BIG)
C
      RETURN
      END
