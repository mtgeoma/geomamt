      SUBROUTINE G03FAZ(N,X,ROW,TOTAL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOTAL
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  ROW(N), X(N*(N+1)/2)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, IJ, J
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      TOTAL = 0.0D0
      DO 20 I = 1, N
         ROW(I) = 0.0D0
   20 CONTINUE
      IJ = 0
      DO 60 J = 1, N
         DO 40 I = 1, J - 1
            IJ = IJ + 1
            ROW(I) = ROW(I) + X(IJ)
            ROW(J) = ROW(J) + X(IJ)
   40    CONTINUE
         IJ = IJ + 1
         ROW(J) = ROW(J) + X(IJ)
   60 CONTINUE
      DO 80 I = 1, N
         TOTAL = TOTAL + ROW(I)
         ROW(I) = ROW(I)/DBLE(N)
   80 CONTINUE
      TOTAL = TOTAL/DBLE(N*N)
      IJ = 0
      DO 120 J = 1, N
         TEMP = ROW(J) - TOTAL
         DO 100 I = 1, J
            IJ = IJ + 1
            X(IJ) = X(IJ) - ROW(I) - TEMP
  100    CONTINUE
  120 CONTINUE
      RETURN
      END
