      SUBROUTINE G02GBW(N,FV,T,VAR,WT,NO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C
C     CALCULATES THE VARIANCE FUNCTION FOR BINOMIAL GLM
C
C     .. Scalar Arguments ..
      INTEGER           N, NO
C     .. Array Arguments ..
      DOUBLE PRECISION  FV(N), T(N), VAR(N), WT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  D
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IF (N.NE.NO) THEN
         DO 20 I = 1, N
            IF (WT(I).EQ.0.0D0 .OR. T(I).EQ.0.0D0) THEN
               VAR(I) = 0.0D0
            ELSE
               D = FV(I)*(T(I)-FV(I))
               VAR(I) = SQRT(T(I)/D)
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = 1, N
            D = FV(I)*(T(I)-FV(I))
            VAR(I) = SQRT(T(I)/D)
   40    CONTINUE
      END IF
      RETURN
      END
