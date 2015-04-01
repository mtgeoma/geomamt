      SUBROUTINE D01GYF(N,NPTS,VK,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     CALCULATION OF THE KOROBOV OPTIMAL COEFFICIENTS
C     FOR A PRIME NUMBER OF POINTS.
C
C     SEE KOROBOV,NUMBER THEORETIC METHODS IN APPROXIMATE ANALYSIS,
C     FIZMATGIZ,MOSCOW,1963.
C
C     INPUT ARGUMENTS
C
C     N        -INTEGER NUMBER OF DIMENSIONS.
C
C     NPTS     -PRIME INTEGER SPECIFYING THE NUMBER OF POINTS
C               TO BE USED.
C
C     IFAIL    -INTEGER VARIABLE SPECIFYING THE ERROR
C               REPORTING OPTION. IFAIL IS SET TO 0 FOR HARD
C               FAIL AND TO 1 FOR SOFT FAIL REPORT
C
C     OUTPUT ARGUMENTS
C
C     VK       -REAL ARRAY WHOSE FIRST N ELEMENTS CONTAIN
C               THE OPTIMAL COEFFICIENTS.
C
C     IFAIL    -THIS REPORTS THE ERROR CONDITIONS AS FOLLOWS
C               IFAIL=0   NO ERROR REPORTED
C               IFAIL=1   N .LT. 1
C               IFAIL=2   NPTS .LT. 5
C               IFAIL=3   NPTS IS NOT PRIME
C               IFAIL=4   INSUFFICIENT PRECISION
C
C     VALIDITY CHECK
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01GYF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  VK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  HMIN, PROD, SUM, T, TLIM, XJ, XN, XP
      INTEGER           I, ICHECK, IP, J, JMIN, K, L, L1, L2, LR, LS, M,
     *                  M1, M2, NDIM
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02ALF
      INTEGER           D01GYZ, P01ABF
      EXTERNAL          X02AJF, X02ALF, D01GYZ, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         AINT, MIN, MOD, DBLE, INT
C     .. Executable Statements ..
      ICHECK = 1
      IF (N.LT.1) GO TO 120
      ICHECK = 2
      IF (NPTS.LT.5) GO TO 120
      ICHECK = 3
      IF (D01GYZ(NPTS).NE.0) GO TO 120
      ICHECK = 0
      NDIM = N
      VK(1) = 1.0D0
      IF (NDIM.EQ.1) GO TO 120
      IP = NPTS/2
      XN = DBLE(NPTS)
      XP = 1.0D0/XN
      HMIN = X02ALF()
      TLIM = 1.0D0/X02AJF()
      JMIN = 1
C     SCAN FOR MINIMUM
      DO 80 J = 2, IP
C        SOLVE J*X=1(MOD NPTS)
         L1 = NPTS
         L2 = J
         M1 = 0
         M2 = 1
   20    LS = L1/L2
         LR = L1 - LS*L2
         M = LS*M2 + M1
         M1 = M2
         M2 = M
         L1 = L2
         L2 = LR
         IF (L2.NE.1) GO TO 20
         IF (J.GT.MIN(M,NPTS-M)) GO TO 80
C        CALCULATE H(J) FUNCTION
         SUM = 0.0D0
         XJ = DBLE(J)
         DO 60 K = 1, IP
            PROD = 1.0D0
            L1 = K
            DO 40 L = 1, NDIM
               PROD = PROD*DBLE(NPTS-L1-L1)*XP
               T = DBLE(L1)*XJ
               IF (T.GT.TLIM) GO TO 140
               L1 = INT(MOD(T,XN)+0.5D0)
   40       CONTINUE
            SUM = SUM + PROD*PROD
            IF (SUM.GE.HMIN) GO TO 80
   60    CONTINUE
         HMIN = SUM
         JMIN = J
   80 CONTINUE
C     GENERATE COEFFICIENT VECTOR
      VK(2) = DBLE(JMIN)
      IF (NDIM.EQ.2) GO TO 120
      XJ = DBLE(JMIN)
      DO 100 I = 3, NDIM
         VK(I) = AINT(MOD(VK(I-1)*XJ,XN)+0.5D0)
  100 CONTINUE
  120 IFAIL = P01ABF(IFAIL,ICHECK,SRNAME,0,P01REC)
      RETURN
  140 ICHECK = 4
      GO TO 120
      END
