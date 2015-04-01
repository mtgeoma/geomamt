      SUBROUTINE G02GCW(N,FV,T,VAR,WT,NO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C
C     CALCULATES VARIANCE FUNCTION FOR POISSON GLM
C
C     .. Scalar Arguments ..
      INTEGER           N, NO
C     .. Array Arguments ..
      DOUBLE PRECISION  FV(N), T(*), VAR(N), WT(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IF (N.NE.NO) THEN
         DO 20 I = 1, N
            IF (WT(I).EQ.0.0D0) THEN
               VAR(I) = 0.0D0
            ELSE
               VAR(I) = 1.0D0/SQRT(FV(I))
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = 1, N
            VAR(I) = 1.0D0/SQRT(FV(I))
   40    CONTINUE
      END IF
      RETURN
      END
