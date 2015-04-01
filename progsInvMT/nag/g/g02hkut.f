      SUBROUTINE G02HKU(T,B,U,UD,W,WD)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, U, UD, W, WD
C     .. Array Arguments ..
      DOUBLE PRECISION  B(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, B2, C, T2
C     .. Executable Statements ..
C
C  UCV FUNCTION
C
      A2 = B(1)
      B2 = B(2)
      U = 1.0D0
      UD = 0.0D0
      IF (T.NE.0) THEN
         T2 = T*T
         IF (T2.GT.B2) THEN
            U = B2/T2
            UD = -2.0D0*U/T
         ELSE IF (T2.LT.A2) THEN
            U = A2/T2
            UD = -2.0D0*U/T
         END IF
      END IF
C
C     W FUNCTION
C
      C = B(3)
      IF (T.GT.C) THEN
         W = C/T
         WD = -W/T
      ELSE
         W = 1.0D0
         WD = 0.0D0
      END IF
      RETURN
      END
