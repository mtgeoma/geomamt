      DOUBLE PRECISION FUNCTION F04JGS(N,A,NRA)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (UPPNRM)
C
C     F04JGS RETURNS THE EUCLIDEAN NORM OF THE N*N UPPER TRIANGULAR
C     MATRIX A.
C
C     NRA MUST BE THE ACTUAL ROW DIMENSION OF A AS DECLARED IN THE
C     CALLING PROGRAM AND MUST BE AT LEAST N.
C
C     ONLY THE UPPER TRIANGULAR PART OF A IS REFERENCED.
C
C     .. Scalar Arguments ..
      INTEGER                          N, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(NRA,N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 BIG, SCALE, SMALL, SUMSQ, TINY
      INTEGER                          J
      LOGICAL                          UNDFLW
C     .. External Functions ..
      DOUBLE PRECISION                 F04JGU, X02AMF
      LOGICAL                          X02DAF
      EXTERNAL                         F04JGU, X02AMF, X02DAF
C     .. External Subroutines ..
      EXTERNAL                         F04JGT
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      TINY = SQRT(SMALL)
      BIG = 1.0D0/SMALL
      UNDFLW = X02DAF(0.0D0)
C
      SCALE = 0.0D0
      SUMSQ = 1.0D0
C
      DO 20 J = 1, N
C
         CALL F04JGT(J,A(1,J),SCALE,SUMSQ,TINY,UNDFLW)
C
   20 CONTINUE
C
      F04JGS = F04JGU(SCALE,SUMSQ,BIG)
C
      RETURN
      END
