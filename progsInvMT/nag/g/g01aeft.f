      SUBROUTINE G01AEF(N,K2,X,ICLASS,CINT,IFREQ,XMIN,XMAX,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED IER-49/42
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-510 (AUG 1986).
C
C     G01AEF - THE ROUTINE CONSTRUCTS A FREQUENCY DISTRIBUTION
C     FOR VALUES IN AN ARRAY X WITH EITHER ROUTINE
C     CALCULATED OR USER SUPPLIED CLASS INTERVALS
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01AEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           ICLASS, IFAIL, K2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CINT(K2), X(N)
      INTEGER           IFREQ(K2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, CR, STEP, XTEMP
      INTEGER           I, IERR, J, JJ, K0, LDN, LUP, NDN, NOB, NUP
      LOGICAL           LK0
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      IF (K2-2) 480, 20, 20
   20 IF (N-1) 500, 40, 40
C     ZERO  FREQUENCIES
   40 DO 60 I = 1, K2
         IFREQ(I) = 0
   60 CONTINUE
      NUP = 0
      NDN = 0
      K0 = K2 - 2
      LK0 = K0 .EQ. 0
C     DETERMINE  MIN  AND  MAX
      XMIN = X(1)
      XMAX = XMIN
      IF (N.EQ.1) GO TO 100
      DO 80 I = 2, N
         XTEMP = X(I)
         IF (XTEMP.LT.XMIN) XMIN = XTEMP
         IF (XTEMP.GT.XMAX) XMAX = XTEMP
   80 CONTINUE
  100 IF (ICLASS.EQ.1) GO TO 200
      IF (LK0) GO TO 420
C     CALCULATE CLASS INTERVALS WITH EQUAL SPACING
      CR = XMAX - XMIN
      ALPHA = 0.001D0
      IF (CR.GT.0.0D0) GO TO 160
      IF (XMIN) 120, 140, 120
  120 CR = ABS(XMIN)*ALPHA
      GO TO 160
  140 CR = ALPHA
  160 STEP = CR*(1.0D0+ALPHA)/DBLE(K0)
      CR = XMIN - 0.5D0*ALPHA*CR
      CINT(1) = CR
      DO 180 I = 1, K0
         XTEMP = CR
         CR = XTEMP + STEP
         CINT(I+1) = CR
  180 CONTINUE
      GO TO 260
C     CLASS INTERVALS SUPPLIED--CHECK IN ASCENDING ORDER
  200 CR = CINT(1)
      IF (LK0) GO TO 440
      JJ = K0 + 1
      DO 240 I = 2, JJ
         XTEMP = CINT(I)
         IF (CR.LT.XTEMP) GO TO 220
         IERR = 3
         GO TO 520
  220    CR = XTEMP
  240 CONTINUE
C     DETERMINE CLASS TO WHICH EACH CASE BELONGS
  260 NOB = K2/2
      CR = CINT(NOB)
      LDN = NOB + 1
      LUP = K2 - 1
      K0 = NOB - 1
      LK0 = K0 .EQ. 0
      DO 380 I = 1, N
         XTEMP = X(I)
         IF (XTEMP-CR) 320, 280, 280
C        VARIATE VALUE ABOVE MIDDLE  INTERVAL
  280    DO 300 J = LDN, LUP
            IF (XTEMP.GE.CINT(J)) GO TO 300
            IFREQ(J) = IFREQ(J) + 1
            GO TO 380
  300    CONTINUE
         NUP = NUP + 1
         GO TO 380
C        VARIATE VALUE BELOW MIDDLE  INTERVAL
  320    IF (LK0) GO TO 360
         DO 340 J = 1, K0
            JJ = NOB - J
            IF (XTEMP.LT.CINT(JJ)) GO TO 340
            JJ = JJ + 1
            IFREQ(JJ) = IFREQ(JJ) + 1
            GO TO 380
  340    CONTINUE
  360    NDN = NDN + 1
  380 CONTINUE
  400 IFREQ(1) = NDN
      IFREQ(K2) = NUP
      IFAIL = 0
      RETURN
  420 CR = (XMAX+XMIN)*0.5D0
      CINT(1) = CR
  440 DO 460 I = 1, N
         IF (X(I).LT.CR) NDN = NDN + 1
  460 CONTINUE
      NUP = N - NDN
      GO TO 400
  480 IERR = 1
      GO TO 520
  500 IERR = 2
  520 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
