      SUBROUTINE G01FMZ(P,V,IR,QL,QU,PL,PU,IERR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     COMPUTES AN APPPROXIMATE INTERVAL FOR THE QUANTILE.
C     FROM THE STUDENTIZED RANGE DISTRIBUTION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, PL, PU, QL, QU, V
      INTEGER           IERR, IR
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, C3, C4, C5, FOUR, HALF, ONE, Q, R, STEP,
     *                  T, VMAX
      INTEGER           IFAIL, K, KMAX
      LOGICAL           LESS
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF, G01EMF
      EXTERNAL          G01CEF, G01EMF
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, MAX, DBLE
C     .. Data statements ..
      DATA              VMAX, HALF, ONE, FOUR, C1, C2, C3, C4,
     *                  C5/120.0D0, 0.5D0, 1.0D0, 4.0D0, 0.8843D0,
     *                  0.2368D0, 1.214D0, 1.208D0, 1.4142D0/
C     .. Executable Statements ..
C
      LESS = .FALSE.
      R = DBLE(IR)
C
      IFAIL = 0
      T = G01CEF(HALF+HALF*P,IFAIL)
      IF (V.LT.VMAX) T = T + (T*T*T+T)/(V*FOUR)
      Q = C1 - C2*T
      IF (V.LT.VMAX) Q = Q - C3/V + C4*T/V
      QL = T*(Q*LOG(R-ONE)+C5)
      IFAIL = 1
      PL = G01EMF(QL,V,IR,IFAIL)
      IF (V.LE.5.0D0) THEN
         STEP = MAX(1.0D0,QL)
      ELSE IF (P.LE.0.4D0) THEN
         IF (IR.LE.5) THEN
            STEP = 1.0D0
         ELSE IF (IR.LE.13) THEN
            STEP = 2.0D0
         ELSE IF (IR.LE.40) THEN
            STEP = 3.0D0
         ELSE IF (IR.LE.100) THEN
            STEP = 4.0D0
         ELSE
            STEP = 5.0D0
         END IF
      ELSE IF (P.LE.0.6D0) THEN
         IF (IR.LE.6) THEN
            STEP = 1.0D0
         ELSE
            STEP = 2.0D0
         END IF
      ELSE IF (IR.LE.5) THEN
         STEP = 0.5D0
      ELSE
         STEP = 1.0D0
      END IF
      K = 0
      KMAX = 100
   20 CONTINUE
      K = K + 1
      IF (K.GT.KMAX) THEN
         IERR = 2
         GO TO 40
      END IF
      IF (PL.LT.P) THEN
         QU = QL + STEP
         LESS = .TRUE.
      ELSE
         QU = QL - STEP
      END IF
      IFAIL = 1
      PU = G01EMF(QU,V,IR,IFAIL)
      IF (PU.LT.P .AND. (LESS)) THEN
         QL = QU
         PL = PU
         GO TO 20
      ELSE IF (PU.GE.P .AND. .NOT. (LESS)) THEN
         QL = QU
         PL = PU
         GO TO 20
      END IF
   40 RETURN
      END
