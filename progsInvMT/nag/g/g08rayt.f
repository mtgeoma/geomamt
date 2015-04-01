      SUBROUTINE G08RAY(ZIN,ETA,VAPVEC,Y,N,N1,TOL,INDW)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CORRECTS EXPECTED VALUES AND VARIANCE-COVARIANCE MATRIX OF
C     PARTICULAR FUNCTION OF ORDER STATISTICS WHEN SOME
C     OBSERVATIONS ARE TIED.
C
C     PETTITT, A.N.P. INFERENCE FOR THE LINEAR MODEL USING A LIKELIHOOD
C                     BASED ON RANKS.
C                     JRSS, B, 44, PP 234-243.
C
C     ARGUMENTS :
C               ZIN - ON ENTRY, CONTAINS EXPECTED VALUES OF FUNCTION
C                     OF ORDER STATISTICS AND, ON EXIT, THE
C                     ADJUSTED VALUES.
C               ETA - ON ENTRY, CONTAINS EXPECTED VALUES OF
C                     DERIVATIVE OF FUNCTION OF ORDER STATISTICS AND,
C                     ON EXIT, THE ADJUSTED VALUES.
C            VAPVEC - ON ENTRY, CONTAINS THE VARIANCE-COVARIANCE
C                     MATRIX OF THE FUNCTION OF THE ORDER STATISTICS
C                     AND, ON EXIT, THE ADJUSTED MATRIX.
C                 Y - VECTOR OF OBSERVATIONS ARRANGED IN ASCENDING
C                     ORDER.
C                 N - SAMPLE SIZE.
C                N1 - DIMENSION OF VECTOR REQUIRED TO STORE UPPER
C                     TRIANGLE OF SQUARE SYMMETRIC MATRIX OF DIMENSION
C                     N. N1 .GE. N(N+1)/2.
C               TOL - TOLERANCE CRITERION FOR TIES.
C              INDW - INTEGER ARRAY OF DIMENSION N USED AS WORKSPACE.
C
C
C     FIND LOCATION AND LENGTH OF SETS OF TIES
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), VAPVEC(N1), Y(N), ZIN(N)
      INTEGER           INDW(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, CRCT, EMOD, SUM1, SUM2, SUM3, SUM4, VMOD,
     *                  ZMOD
      INTEGER           I, IC, ICO, II, IS, J, JC, JJ, JS, M
C     .. External Functions ..
      INTEGER           G01DCU
      EXTERNAL          G01DCU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      DO 20 I = 1, N
         INDW(I) = 0
   20 CONTINUE
      ICO = 0
      I = 0
   40 I = I + 1
      IF (I.EQ.N) GO TO 100
      A = Y(I)
      DO 60 J = I + 1, N
         IF (ABS(A-Y(J)).GT.TOL) GO TO 80
   60 CONTINUE
      J = N + 1
   80 IF ((J-I).EQ.1) GO TO 40
      ICO = ICO + 1
      INDW(2*(ICO-1)+1) = I
      INDW(2*(ICO-1)+2) = J - I
      I = J - 1
      IF (I.LE.N-1) GO TO 40
  100 CONTINUE
C
C     IF THERE ARE NO TIES OR ALL OBSERVATIONS ARE TIED FINISH.
C
      IF (ICO.EQ.0 .OR. INDW(2).EQ.N) RETURN
C
      DO 400 II = 1, ICO
C
C        CHANGE MEANS AND VARIANCES OF EACH SET OF TIED OBSERVATIONS.
C
         IS = INDW(2*(II-1)+1)
         IC = INDW(2*(II-1)+2)
         SUM1 = 0.0D0
         SUM2 = 0.0D0
         SUM3 = 0.0D0
         SUM4 = 0.0D0
         DO 120 I = 1, IC
            SUM1 = SUM1 + ZIN(I+IS-1)
            SUM2 = SUM2 + ETA(I+IS-1)
            M = G01DCU(I+IS-1,I+IS-1)
            SUM3 = SUM3 + VAPVEC(M)
            SUM4 = SUM4 + ZIN(I+IS-1)**2
  120    CONTINUE
         CRCT = SUM4/IC - (SUM1/IC)*(SUM1/IC)
         ZMOD = SUM1/IC
         EMOD = SUM2/IC
         VMOD = SUM3/IC + CRCT
         DO 140 I = 1, IC
            ZIN(I+IS-1) = ZMOD
            ETA(I+IS-1) = EMOD
            M = G01DCU(I+IS-1,I+IS-1)
            VAPVEC(M) = VMOD
  140    CONTINUE
C
C        CHANGE COVARIANCES OF TIED OBSERVATIONS IN THE SAME SET
C
         SUM3 = 0.0D0
         DO 180 I = 1, IC
            DO 160 J = I + 1, IC
               M = G01DCU(I+IS-1,J+IS-1)
               SUM3 = SUM3 + VAPVEC(M)
  160       CONTINUE
  180    CONTINUE
         VMOD = 2*SUM3/(IC*(IC-1)) - CRCT/(IC-1)
         DO 220 I = 1, IC
            DO 200 J = I + 1, IC
               M = G01DCU(I+IS-1,J+IS-1)
               VAPVEC(M) = VMOD
  200       CONTINUE
  220    CONTINUE
C
C        CHANGE COVARIANCES OF DIFFERENT SETS OF TIED OBSERVATIONS.
C
         DO 320 JJ = II + 1, ICO
            JS = INDW(2*(JJ-1)+1)
            JC = INDW(2*(JJ-1)+2)
            SUM3 = 0.0D0
            DO 260 I = 1, IC
               DO 240 J = 1, JC
                  M = G01DCU(I+IS-1,J+JS-1)
                  SUM3 = SUM3 + VAPVEC(M)
  240          CONTINUE
  260       CONTINUE
            VMOD = SUM3/(IC*JC)
            DO 300 I = 1, IC
               DO 280 J = 1, JC
                  M = G01DCU(I+IS-1,J+JS-1)
                  VAPVEC(M) = VMOD
  280          CONTINUE
  300       CONTINUE
  320    CONTINUE
C
C        CHANGE COVARIANCES OF TIED ELEMENTS WITH UNTIED.
C
         J = 0
         JJ = 1
  340    J = J + 1
         IF (J.GT.N) GO TO 400
         IF (J.EQ.INDW(2*(JJ-1)+1)) THEN
            J = J + INDW(2*(JJ-1)+2) - 1
            IF (JJ.LT.ICO) JJ = JJ + 1
            GO TO 340
         END IF
C
C        HAVE FOUND AN UNTIED OBSERVATION
C
         SUM3 = 0.0D0
         DO 360 I = 1, IC
            M = G01DCU(I+IS-1,J)
            SUM3 = SUM3 + VAPVEC(M)
  360    CONTINUE
         VMOD = SUM3/IC
         DO 380 I = 1, IC
            M = G01DCU(I+IS-1,J)
            VAPVEC(M) = VMOD
  380    CONTINUE
         GO TO 340
  400 CONTINUE
      RETURN
      END
