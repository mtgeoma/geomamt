      SUBROUTINE G13CAY(XG,NX,MTX)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13CAY MEAN OR TREND CORRECTS THE DATA
C
C     XG    - DATA ARRAY
C     NX    - NUMBER OF DATA POINTS
C     MTX   - MEAN TREND CORRECTION INDICATOR
C     0 NO CORRECTION
C     1 MEAN CORRECTION
C     2 TREND CORRECTION
C
C
C     .. Scalar Arguments ..
      INTEGER           MTX, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  XG(NX)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, S1, S2, S3
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IF (MTX-1) 140, 80, 20
   20 S1 = 0.0D0
      S2 = 0.0D0
      DO 40 I = 1, NX
         S1 = S1 + XG(I)
         S2 = S2 + DBLE(I)*XG(I)
   40 CONTINUE
      S3 = DBLE(NX)*DBLE(NX+1)/2.0D0
      A = S3*DBLE(2*NX+1)/3.0D0
      B = DBLE(NX)*A - S3*S3
      A = (A*S1-S3*S2)/B
      B = (DBLE(NX)*S2-S3*S1)/B
      DO 60 I = 1, NX
         XG(I) = XG(I) - A - B*DBLE(I)
   60 CONTINUE
      GO TO 140
   80 S1 = 0.0D0
      DO 100 I = 1, NX
         S1 = S1 + XG(I)
  100 CONTINUE
      A = S1/DBLE(NX)
      DO 120 I = 1, NX
         XG(I) = XG(I) - A
  120 CONTINUE
  140 CONTINUE
      RETURN
      END
