      SUBROUTINE D01GCF(N,F,REGION,NPTS,VK,NRAND,ITRANS,RES,ERR,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     MULTIPLE INTEGRATION
C
C     THIS SUBROUTINE CALCULATES AN APPROXIMATION TO
C            D1      DN
C            I ......I F(X1,...,XN) DX1 .. DXN
C            C1      CN
C
C     USING THE KOROBOV-CONROY NUMBER-THEORETIC METHOD.
C     (N.M.KOROBOV,NUMBER THEORETIC METHODS IN APPROXIMATE ANALYSIS,
C     FIZMATGIZ,MOSCOW,1963, H.CONROY,J.CHEM.PHYS.47,(1967),
C     5307-5813)
C
C     INPUT ARGUMENTS
C
C     N        -INTEGER SPECIFYING THE NUMBER OF DIMENSIONS.
C               1 .LE. N .LE. 20.
C
C     F        -EXTERNALLY DECLARED REAL USER FUNCTION INTEGRAND
C               HAVING ARGUMENTS (N,X) WHERE X IS A REAL ARRAY
C               OF DIMENSION N WHOSE ELEMENTS ARE THE VARIABLE
C               VALUES.
C
C     REGION   -EXTERNALLY DECLARED USER SUBROUTINE WITH ARGUMENTS
C               (N,X,J,C,D) WHICH CALCULATES THE LOWER LIMIT C
C               AND THE UPPER LIMIT D CORRESPONDING TO THE ARRAY
C               VARIABLE X(J). C AND D MAY DEPEND ON X(1)..X(J-1).
C
C     NPTS     -INTEGER VARIABLE WHICH SPECIFIES THE KOROBOV RULE
C               TO BE USED. THERE ARE TWO ALTERNATIVES DEPENDING
C               ON THE VALUE OF NPTS
C               (A) 1 .LE. NPTS .LE. 6
C                   IN THIS CASE ONE OF THE SIX PRESET RULES IS
C                   CHOSEN USING 2129, 5003, 10007, 20011, 40009 OR
C                   80021 POINTS ACCORDING TO THE RESPECTIVE VALUE
C                   OF NPTS.
C               (B) NPTS .GT. 6
C                   NPTS IS THE NUMBER OF ACTUAL POINTS TO BE USED
C                   WITH CORRESPONDING OPTIMAL COEFFICIENTS SUPPLIED
C                   BY THE USER IN ARRAY VK.
C
C     VK       -REAL ARRAY CONTAINING N OPTIMAL COEFFICIENTS
C               FOR THE CASE NPTS .GT. 6.
C
C     NRAND    -INTEGER NUMBER SPECIFYING THE NUMBER OF RANDOM
C               SAMPLES TO BE GENERATED IN THE ERROR ESTIMATION.
C               (GENERALLY A SMALL VALUE ,SAY 3 TO 5, IS SUFFICIENT)
C               THE TOTAL NUMBER OF INTEGRAND EVALUATIONS WILL BE
C               NRAND*NPTS.
C
C     ITRANS   -THIS SHOULD BE SET TO 0 OR 1. THE NORMAL SETTING
C               IS 0. WHEN SET TO 1 THE PERIODIZING TRANSFORMATION
C               IS SUPPRESSED (TO COVER CASES WHERE THE INTEGRAND IS
C               ALREADY PERIODIC OR WHERE THE USER DESIRES TO
C               SPECIFY A PARTICULAR TRANSFORMATION IN THE
C               DEFINITION OF F).
C
C     IFAIL    -INTEGER VARIABLE SPECIFYING THE ERROR REPORTING
C               OPTION. IFAIL IS SET TO 0 FOR HARD FAIL AND TO 1 FOR
C               SOFT FAIL REPORT.
C
C
C     OUTPUT ARGUMENTS
C
C     VK       -WHEN THE SUBROUTINE IS CALLED WITH 1 .LE. NPTS
C               .LE. 6 THEN N ELEMENTS OF ARRAY VK WILL BE SET
C               TO THE OPTIMAL COEFFICIENTS USED BY THE PRESET RULE.
C
C     RES      -APPROXIMATION TO THE VALUE OF THE INTEGRAL.
C
C     ERR      -STANDARD ERROR AS COMPUTED FROM NRAND SAMPLE VALUES.
C               IF NRAND IS SET TO 1 THEN ERR WILL BE RETURNED WITH
C               VALUE 0.
C
C     IFAIL    -THIS REPORTS THE ERROR CONDITIONS AS FOLLOWS
C               IFAIL=0  NO ERROR DETECTED.
C                    =1  N FAILS TO SATISFY 1 .LE. N .LE. 20
C                    =2  NPTS .LT. 1.
C                    =3  NRAND .LT. 1.
C
C     ************* IMPLEMENTATION-DEPENDENT CODING ****************
C     THIS ROUTINE REQUIRES AT ONE POINT TO PERFORM EXACT INTEGER
C     ARITHMETIC WITH AT LEAST 32 BITS AVAILABLE FOR INTERMEDIATE
C     RESULTS. THE INTEGERS ARE STORED AS FLOATING-POINT VARIABLES
C     (IN ORDER TO AVOID INTEGER OVERFLOW ON SOME MACHINES). IF THE
C     MACHINE DOES NOT PROVIDE AT LEAST 32 BITS OF BASIC FLOATING-
C     POINT PRECISION, THEN ONE STATEMENT IN THE CODE MUST BE
C     CONVERTED TO ADDITIONAL PRECISION (SEE COMMENTS IN CODE).
C     **************************************************************
C
C     REGION
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01GCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERR, RES
      INTEGER           IFAIL, ITRANS, N, NPTS, NRAND
C     .. Array Arguments ..
      DOUBLE PRECISION  VK(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Subroutine Arguments ..
      EXTERNAL          REGION
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, DPK, GRAD, PK, PTS, SUM, SUM2, WT, XA, XX
      INTEGER           I, ICHECK, J, K, M, MAXDIM, MAXRUL, NDIM, NP
C     .. Local Arrays ..
      DOUBLE PRECISION  ALPHA(20), KPTS(6), KR(6,20), X(20)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, AINT, MOD, DBLE, SQRT
C     .. Data statements ..
C
C     KOROBOV COEFFICIENTS FOR 1 TO 20 DIMENSIONS.
C
C     (ACKNOWLEDGEMENT-CALCULATED  ON THE CRAY I COMPUTER AT SANDIA
C     NATIONAL LABORATORY, LIVERMORE, CALIFORNIA, USA, THROUGH
C     THE KIND ASSISTANCE OF DR. R.E. HUDDLESTON)
C
C     KR(I,J) CONTAINS KOROBOV COEFFICIENT FOR RULE I DIMENSION J.
C     ASSOCIATED NO. OF POINTS IS GIVEN BY KPTS(I)
C
C
C     KPTS(I) = NO. OF POINTS ASSOCIATED WITH RULE I.
C
      DATA              MAXRUL, MAXDIM/6, 20/
      DATA              KR(1,1), KR(2,1), KR(3,1), KR(4,1), KR(5,1),
     *                  KR(6,1)/1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0/
      DATA              KR(1,2), KR(2,2), KR(3,2), KR(4,2), KR(5,2),
     *                  KR(6,2)/780.D0, 1850.D0, 3822.D0, 6103.D0,
     *                  15152.D0, 30954.D0/
      DATA              KR(1,3), KR(2,3), KR(3,3), KR(4,3), KR(5,3),
     *                  KR(6,3)/359.D0, 1476.D0, 544.D0, 4104.D0,
     *                  16592.D0, 19394.D0/
      DATA              KR(1,4), KR(2,4), KR(3,4), KR(4,4), KR(5,4),
     *                  KR(6,4)/766.D0, 792.D0, 1206.D0, 6016.D0,
     *                  9023.D0, 15710.D0/
      DATA              KR(1,5), KR(2,5), KR(3,5), KR(4,5), KR(5,5),
     *                  KR(6,5)/618.D0, 840.D0, 198.D0, 6019.D0,
     *                  12216.D0, 2302.D0/
      DATA              KR(1,6), KR(2,6), KR(3,6), KR(4,6), KR(5,6),
     *                  KR(6,6)/41.D0, 2037.D0, 2240.D0, 4167.D0,
     *                  4902.D0, 9227.D0/
      DATA              KR(1,7), KR(2,7), KR(3,7), KR(4,7), KR(5,7),
     *                  KR(6,7)/596.D0, 229.D0, 2304.D0, 3851.D0,
     *                  12506.D0, 3420.D0/
      DATA              KR(1,8), KR(2,8), KR(3,8), KR(4,8), KR(5,8),
     *                  KR(6,8)/86.D0, 1578.D0, 436.D0, 4138.D0,
     *                  7824.D0, 3824.D0/
      DATA              KR(1,9), KR(2,9), KR(3,9), KR(4,9), KR(5,9),
     *                  KR(6,9)/636.D0, 526.D0, 470.D0, 259.D0, 6093.D0,
     *                  22300.D0/
      DATA              KR(1,10), KR(2,10), KR(3,10), KR(4,10),
     *                  KR(5,10), KR(6,10)/287.D0, 431.D0, 1554.D0,
     *                  1117.D0, 12088.D0, 5130.D0/
      DATA              KR(1,11), KR(2,11), KR(3,11), KR(4,11),
     *                  KR(5,11), KR(6,11)/707.D0, 1485.D0, 480.D0,
     *                  1188.D0, 2399.D0, 11222.D0/
      DATA              KR(1,12), KR(2,12), KR(3,12), KR(4,12),
     *                  KR(5,12), KR(6,12)/707.D0, 1450.D0, 1004.D0,
     *                  173.D0, 8764.D0, 17698.D0/
      DATA              KR(1,13), KR(2,13), KR(3,13), KR(4,13),
     *                  KR(5,13), KR(6,13)/96.D0, 1001.D0, 684.D0,
     *                  2919.D0, 5491.D0, 7057.D0/
      DATA              KR(1,14), KR(2,14), KR(3,14), KR(4,14),
     *                  KR(5,14), KR(6,14)/49.D0, 1001.D0, 684.D0,
     *                  235.D0, 9274.D0, 28739.D0/
      DATA              KR(1,15), KR(2,15), KR(3,15), KR(4,15),
     *                  KR(5,15), KR(6,15)/373.D0, 1001.D0, 1447.D0,
     *                  3043.D0, 3054.D0, 33207.D0/
      DATA              KR(1,16), KR(2,16), KR(3,16), KR(4,16),
     *                  KR(5,16), KR(6,16)/613.D0, 2.D0, 857.D0,
     *                  1249.D0, 2648.D0, 27717.D0/
      DATA              KR(1,17), KR(2,17), KR(3,17), KR(4,17),
     *                  KR(5,17), KR(6,17)/373.D0, 2.D0, 2.D0, 1249.D0,
     *                  2648.D0, 33207.D0/
      DATA              KR(1,18), KR(2,18), KR(3,18), KR(4,18),
     *                  KR(5,18), KR(6,18)/2.D0, 2.D0, 2.D0, 2.D0,
     *                  2648.D0, 1420.D0/
      DATA              KR(1,19), KR(2,19), KR(3,19), KR(4,19),
     *                  KR(5,19), KR(6,19)/2.D0, 2.D0, 2.D0, 2.D0, 2.D0,
     *                  1420.D0/
      DATA              KR(1,20), KR(2,20), KR(3,20), KR(4,20),
     *                  KR(5,20), KR(6,20)/2.D0, 2.D0, 2.D0, 2.D0, 2.D0,
     *                  2.D0/
      DATA              KPTS(1), KPTS(2), KPTS(3), KPTS(4), KPTS(5),
     *                  KPTS(6)/2129.D0, 5003.D0, 10007.D0, 20011.D0,
     *                  40009.D0, 80021.D0/
C     .. Executable Statements ..
      NDIM = N
      PTS = DBLE(NPTS)
C     VALIDITY CHECK
      ICHECK = 1
      IF (NDIM.LT.1 .OR. NDIM.GT.MAXDIM) GO TO 160
      ICHECK = 2
      IF (NPTS.LT.1) GO TO 160
      ICHECK = 3
      IF (NRAND.LT.1) GO TO 160
      ICHECK = 0
      IF (NPTS.GT.MAXRUL) GO TO 40
C     SELECT KOROBOV VECTOR
      NP = NPTS
      PTS = KPTS(NPTS)
      VK(1) = 1.0D0
      IF (NDIM.EQ.1) GO TO 40
      VK(2) = KR(NP,NDIM)
      IF (NDIM.EQ.2) GO TO 40
      XA = VK(2)
      DO 20 I = 3, NDIM
         VK(I) = AINT(MOD(VK(I-1)*XA,PTS)+0.5D0)
C        ************* IMPLEMENTATION-DEPENDENT CODING ****************
C        IF THE MACHINE DOES NOT PROVIDE AT LEAST 32 BITS OF BASIC
C        PRECISION, THEN INTERMEDIATE WORKING IN THE PRECEDING
C        STATEMENT MUST BE CHANGED TO ADDITIONAL PRECISION, E.G. IN A
C        SINGLE PRECISION IMPLEMENTATION CHANGED TO
C        VK(I) = AINT(SNGL(DMOD(DBLE(VK(I-1))*DBLE(XA),DBLE(PTS)))+
C        *    0.5)
C        THE FUNCTION DMOD MUST BE SUPPLIED IF NECESSARY.
C        **************************************************************
   20 CONTINUE
C     BEGIN INTEGRATION
   40 SUM = 0.0D0
      SUM2 = 0.0D0
      DO 140 M = 1, NRAND
         RES = 0.0D0
C        CALCULATE RANDOM SHIFT
         DO 60 K = 1, NDIM
            ALPHA(K) = G05CAF(ALPHA(K))
   60    CONTINUE
C        CALCULATE TRANSFORMED INTEGRAND
         PK = 1.0D0
         DPK = 1.0D0/PTS
   80    WT = DPK
         DO 120 J = 1, NDIM
            XX = MOD(ALPHA(J)+VK(J)*PK*DPK,1.0D0)
            CALL REGION(NDIM,X,J,C,D)
            GRAD = D - C
            WT = WT*GRAD
            IF (ITRANS.NE.0) GO TO 100
C           PERIODIZING TRANSFORMATION
            WT = WT*6.0D0*XX*(1.0D0-XX)
            XX = XX*XX*(3.0D0-2.0D0*XX)
  100       X(J) = C + GRAD*XX
  120    CONTINUE
         RES = F(NDIM,X)*WT + RES
         PK = AINT(PK+1.5D0)
         IF (PK.LE.PTS) GO TO 80
         SUM = SUM + RES
         SUM2 = SUM2 + RES*RES
  140 CONTINUE
C     MEAN
      RES = SUM/DBLE(NRAND)
      ERR = 0.0D0
      IF (NRAND.EQ.1) GO TO 160
C     STANDARD ERROR
      ERR = SQRT(ABS((SUM2-DBLE(NRAND)*RES*RES)/DBLE(NRAND*(NRAND-1))))
  160 IFAIL = P01ABF(IFAIL,ICHECK,SRNAME,0,P01REC)
      RETURN
      END
