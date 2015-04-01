      SUBROUTINE E02GAZ(V1,V2,MLT,M1,IOUT)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  MLT
      INTEGER           IOUT, M1
C     .. Array Arguments ..
      DOUBLE PRECISION  V1(M1), V2(M1)
C     .. Local Scalars ..
      INTEGER           I, IOUTM1, IOUTP1
C     .. Executable Statements ..
      IOUTP1 = IOUT + 1
      IF (IOUT.EQ.1) GO TO 40
      IOUTM1 = IOUT - 1
      DO 20 I = 1, IOUTM1
         V1(I) = V1(I) + V2(I)*MLT
   20 CONTINUE
   40 DO 60 I = IOUTP1, M1
         V1(I) = V1(I) + V2(I)*MLT
   60 CONTINUE
      RETURN
      END
