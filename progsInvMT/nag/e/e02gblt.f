      SUBROUTINE E02GBL(NUM,SCAL,V,W)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     ADD  SCAL*V  TO  W.
C     STORE THE RESULT IN  W.
C     ***************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SCAL
      INTEGER           NUM
C     .. Array Arguments ..
      DOUBLE PRECISION  V(NUM), W(NUM)
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
      DO 20 I = 1, NUM
         W(I) = W(I) + SCAL*V(I)
   20 CONTINUE
      RETURN
      END
