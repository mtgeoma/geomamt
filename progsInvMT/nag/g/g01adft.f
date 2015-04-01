      SUBROUTINE G01ADF(K,X,IFREQ,XMEAN,S2,S3,S4,N,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 5C REVISED IER-82
C     MARK 9 REVISED. IER-331 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. IER-619 (APR 1988).
C
C     G01ADF - THE ROUTINE COMPUTES MEAN, STD DEVN, AND COEFFTS
C     OF SKEWNESS AND KURTOSIS FOR GROUPED DATA
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01ADF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S2, S3, S4, XMEAN
      INTEGER           IFAIL, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(K)
      INTEGER           IFREQ(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  ATEMP, BTEMP, FREQ, FTEMP, WTEMP, XT2, XTEMP
      INTEGER           I, ICOUNT, IERR, IJ, J, KK, KL1, KM1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      KM1 = K - 1
      IF (KM1-1) 20, 40, 40
   20 IERR = 1
      GO TO 360
C     CHECK CLASS INTERVALS IN ASCENDING ORDER
   40 XTEMP = X(1)
      DO 100 I = 2, K
         XT2 = X(I)
         IF (XTEMP-XT2) 80, 80, 60
   60    IERR = 2
         GO TO 360
   80    XTEMP = XT2
  100 CONTINUE
C     DETERMINE FIRST NON-ZERO FREQ, AND OTHERS GE 0
      DO 240 I = 1, KM1
         IJ = IFREQ(I)
         IF (IJ) 340, 240, 140
  140    ICOUNT = 1
         IF (I-KM1) 160, 220, 220
  160    KK = I + 1
         DO 200 J = KK, KM1
            IF (IFREQ(J)) 340, 200, 180
  180       ICOUNT = ICOUNT + 1
  200    CONTINUE
  220    KL1 = I
         GO TO 260
  240 CONTINUE
      GO TO 340
C     INITIALISE FOR MAIN CALCULATION LOOP
  260 FREQ = 0.0D0
      XMEAN = 0.0D0
      S2 = 0.0D0
      S3 = 0.0D0
      S4 = 0.0D0
C     LOOP FOR MEANS ETC
      DO 280 I = KL1, KM1
         FTEMP = IFREQ(I)
         IF (FTEMP.LE.0.0D0) GO TO 280
         FREQ = FTEMP + FREQ
         XTEMP = (X(I)+X(I+1))*0.5D0 - XMEAN
         ATEMP = FTEMP*XTEMP
         BTEMP = ATEMP/FREQ
         ATEMP = (XTEMP-BTEMP)*ATEMP
         WTEMP = 3.0D0*BTEMP*S2
         S4 = BTEMP*(2.0D0*WTEMP-4.0D0*S3) +
     *        ATEMP*(-3.0D0*ATEMP/FREQ+XTEMP*XTEMP) + S4
         S3 = -WTEMP + ATEMP*(-2.0D0*BTEMP+XTEMP) + S3
         S2 = ATEMP + S2
         XMEAN = BTEMP + XMEAN
  280 CONTINUE
      N = 0.5D0 + FREQ
      IF (FREQ.LT.1.2D0) GO TO 320
      IFAIL = 0
      IF (ICOUNT.LE.1) GO TO 300
      ATEMP = FREQ - 1.0D0
      WTEMP = S2/ATEMP
      S2 = SQRT(WTEMP)
      S3 = S3/ATEMP/(S2*WTEMP)
      S4 = S4/ATEMP/(WTEMP*WTEMP) - 3.0D0
      RETURN
  300 S2 = 0.0D0
      S3 = 0.0D0
      S4 = 0.0D0
      RETURN
  320 IERR = 4
      GO TO 360
  340 IERR = 3
  360 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
