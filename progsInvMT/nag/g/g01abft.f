      SUBROUTINE G01ABF(N,X1,X2,IWT,WT,RES,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G01ABF - THE ROUTINE COMPUTES MEANS, STD DEVNS, MINIMA,
C     MAXIMA, CORRECTED SUMS OF SQUARES AND PRODUCTS,
C     AND THE PRODUCT MOMENT CORRELATION COEFFICIENT
C     FOR TWO VARIABLES IN THE ARRAYS X1 AND X2 -
C     OPTIONALLY USING WEIGHTS IN THE ARRAY WT
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01ABF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IWT, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RES(13), WT(N), X1(N), X2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AT1, AT2, SS1, SS2, SSP, T1, T2, W, WS, WSUM,
     *                  WTEMP, XN1, XN2, XTEMP1, XTEMP2, XX1, XX2
      INTEGER           I, IERR, K, LOW
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, SQRT
C     .. Executable Statements ..
      IF (N) 20, 20, 80
   20 IERR = 1
   40 IWT = 0
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
   80 IF (IWT-1) 100, 140, 100
C     INITIALISE WEIGHTS IF NOT GIVEN
  100 DO 120 I = 1, N
         WT(I) = 1.0D0
  120 CONTINUE
      LOW = 1
      WTEMP = 1.0D0
      IWT = N
      GO TO 300
C     FIND FIRST POSITIVE WEIGHT -- CHECK OTHERS NON-NEG
  140 IWT = 0
      DO 260 I = 1, N
         WTEMP = WT(I)
         IF (WTEMP) 280, 260, 160
  160    IF (I-N) 180, 240, 240
  180    LOW = I + 1
         DO 220 K = LOW, N
            IF (WT(K)) 280, 220, 200
  200       IWT = IWT + 1
  220    CONTINUE
  240    LOW = I
         IWT = IWT + 1
         GO TO 300
  260 CONTINUE
  280 IERR = 3
      GO TO 40
C     INITIALISE MEANS ETC
  300 WSUM = WTEMP
      T1 = X1(LOW)
      T2 = X2(LOW)
      SS1 = 0.0D0
      SSP = 0.0D0
      SS2 = 0.0D0
      XN1 = T1
      XX1 = T1
      XN2 = T2
      XX2 = T2
      IF (IWT-1) 280, 380, 320
  320 LOW = LOW + 1
      W = 0.0D0
C     LOOP FOR CALCULATIONS
      DO 360 I = LOW, N
         WS = WSUM
         WTEMP = WT(I)
         IF (WTEMP) 280, 360, 340
  340    W = 2.0D0*WS*WTEMP + W
         XTEMP1 = X1(I)
         XTEMP2 = X2(I)
         IF (XTEMP1.LT.XN1) XN1 = XTEMP1
         IF (XTEMP2.LT.XN2) XN2 = XTEMP2
         IF (XTEMP1.GT.XX1) XX1 = XTEMP1
         IF (XTEMP2.GT.XX2) XX2 = XTEMP2
         WSUM = WTEMP + WS
         XTEMP1 = XTEMP1 - T1
         XTEMP2 = XTEMP2 - T2
         AT1 = WTEMP*XTEMP1/WSUM
         AT2 = WTEMP*XTEMP2/WSUM
         XTEMP2 = XTEMP2*WS
         SS1 = AT1*WS*XTEMP1 + SS1
         SSP = AT1*XTEMP2 + SSP
         SS2 = AT2*XTEMP2 + SS2
         T1 = AT1 + T1
         T2 = AT2 + T2
  360 CONTINUE
  380 RES(1) = T1
      RES(2) = T2
      SS1 = MAX(SS1,0.0D0)
      RES(5) = SS1
      RES(6) = SSP
      SS2 = MAX(0.0D0,SS2)
      RES(7) = SS2
      RES(9) = XN1
      RES(10) = XX1
      RES(11) = XN2
      RES(12) = XX2
      RES(13) = WSUM
      IF (IWT-1) 280, 400, 420
  400 IERR = 2
      GO TO 60
  420 WTEMP = WSUM/W
      RES(3) = SQRT(SS1*WTEMP)
      RES(4) = SQRT(WTEMP*SS2)
      W = 0.0D0
      IF (SS1.GT.0.0D0 .AND. SS2.GT.0.0D0) W = SSP/SQRT(SS1*SS2)
      RES(8) = W
      IFAIL = 0
      RETURN
      END
