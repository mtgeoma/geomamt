      SUBROUTINE D01BCF(ITYPE,AA,BB,CC,DD,NPNTS,WEIGHT,ABSCIS,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9C REVISED. IER-370 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14A REVISED. IER-677 (DEC 1989).
C     MARK 14B REVISED. IER-840 (MAR 1990).
C     SUBROUTINE FOR THE DETERMINATION OF GAUSSIAN QUADRATURE RULES
C     **************************************************************
C
C     INPUT PARAMETERS
C
C     ITYPE  INTEGER WHICH SPECIFIES THE RULE TYPE CHOSEN
C               WEIGHT W(X)          INTERVAL     RESTRICTIONS
C     0            1                    A,B          B.GT.A
C     1    (B-X)**C*(X-A)**D            A,B     B.GT.A,C,D.GT.-1
C     2   ABS(X-0.5*(A+B))**C           A,B     C.GT.-1,B.GT.A
C     3  ABS(X-A)**C*EXP(-B*X)          A,INF   C.GT.-1,B.GT.0
C     3  ABS(X-A)**C*EXP(-B*X)       -INF,A     C.GT.-1,B.LT.0
C     4 ABS(X-A)**C*EXP(-B*(X-A)**2) -INF,INF   C.GT.-1,B.GT.0
C     5  ABS(X-A)**C/ABS(X+B)**D      A,INF A.GT.-B,C.GT.-1,D.GT.C+1
C     5  ABS(X-A)**C/ABS(X+B)**D     -INF,A A.LT.-B,C.GT.-1,D.GT.C+1
C     ABS(ITYPE) MUST BE LESS THAN 6. IF ITYPE IS GIVEN LESS THAN
C     ZERO THEN THE ADJUSTED WEIGHTS ARE CALCULATED. IF NPNTS IS
C     ODD AND ITYPE EQUALS -2 OR -4 AND C IS NOT ZERO, THERE MAY BE
C     PROBLEMS.
C
C     AA     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     BB     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     CC     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     DD     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     IFAIL  NAG FAILURE PARAMETER.  SEE NAG DOCUMENTATION.
C
C     NPNTS  INTEGER THAT DETERMINES DIMENSION OF WEIGHT AND ABSCIS
C
C     OUTPUT PARAMETERS
C
C     WEIGHT  REAL ARRAY OF DIMENSION NPNTS WHICH CONTAINS
C     RULE WEIGHTS
C
C     ABSCIS  REAL ARRAY OF DIMENSION NPNTS WHICH CONTAINS
C     RULE ABSCISSAE
C
C     IFAIL INTEGER NAG FAILURE PARAMETER
C       IFAIL=0 FOR NORMAL EXIT
C       IFAIL=1 FOR FAILURE IN NAG ROUTINE F02AVF
C       IFAIL=2 FOR PARAMETER NPNTS OR ITYPE OUT OF RANGE
C       IFAIL=3 FOR PARAMETER AA OR BB OR CC OR DD OUT OF
C               ALLOWED RANGE
C       IFAIL=4 FOR OVERFLOW IN CALCULATION OF WEIGHTS
C       IFAIL=5 FOR UNDERFLOW IN CALCULATION OF WEIGHTS
C       IFAIL=6 FOR ITYPE=-2 OR -4, NPNTS ODD, C NOT ZERO
C
C     **************************************************************
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AA, BB, CC, DD
      INTEGER           IFAIL, ITYPE, NPNTS
C     .. Array Arguments ..
      DOUBLE PRECISION  ABSCIS(NPNTS), WEIGHT(NPNTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ABSPNC, B, BN, C, CN, CNO, D, FACN, FN, FOUR,
     *                  GAMMA, GAMMAB, GAMMB, HALF, ONE, PNA, PNB, PNC,
     *                  PONORM, PSQRD, REALMX, SMALL, SQRTCN, STORE,
     *                  TWNAPB, TWO, WTSUM, Y, ZERO
      INTEGER           IERROR, ISUB, J, MITYPE, N, NBUG, NFAC, NHALF
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S14AAF, X02AJF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          S14AAF, X02AJF, X02ALF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02AVF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, LOG, EXP, DBLE, SQRT, INT
C     .. Data statements ..
      DATA              ZERO, ONE, TWO, FOUR/0.0D0, 1.0D0, 2.0D0, 4.0D0/
      DATA              HALF/0.5D0/
C     .. Executable Statements ..
C
C        INITIALISATION AND PARAMETER CHECKING
C
      SMALL = X02AJF()
      IF (NPNTS.LE.0) GO TO 780
      DO 20 J = 1, NPNTS
         ABSCIS(J) = ZERO
         WEIGHT(J) = ZERO
   20 CONTINUE
      MITYPE = ABS(ITYPE) + 1
      IF (MITYPE.GT.6) GO TO 780
      A = AA
      B = BB
      C = CC
      D = DD
      GO TO (40,60,100,120,140,160) MITYPE
   40 C = ZERO
      D = ZERO
   60 IF (C.LE.-ONE .OR. D.LE.-ONE) GO TO 800
      IF (B.LE.A) GO TO 800
      PONORM = (HALF*(B-A))**(C+D+ONE)
      IF (ITYPE.LT.0) PONORM = PONORM/(HALF*(B-A))**(C+D)
   80 IERROR = 1
      GAMMA = S14AAF(C+ONE,IERROR)
      IF (IERROR.GT.0) GO TO 800
      IERROR = 1
      GAMMB = S14AAF(D+ONE,IERROR)
      IF (IERROR.GT.0) GO TO 800
      IERROR = 1
      GAMMAB = S14AAF(C+D+TWO,IERROR)
      IF (IERROR.GT.0) GO TO 800
      PONORM = PONORM*TWO**(C+D+ONE)*GAMMA*GAMMB/GAMMAB
      ABSCIS(1) = (D-C)/(C+D+TWO)
      GO TO 180
  100 IF (C.LE.-ONE .OR. B.LE.A) GO TO 800
      PONORM = TWO*(HALF*(B-A))**(C+ONE)/(C+ONE)
      IF (ITYPE.LT.0) PONORM = PONORM/(HALF*(B-A))**C
      GO TO 180
  120 IF (C.LE.-ONE .OR. B.EQ.ZERO) GO TO 800
      IERROR = 1
      PONORM = S14AAF(C+ONE,IERROR)*EXP(-B*A)/ABS(B)**(C+ONE)
      IF (ITYPE.LT.0) PONORM = PONORM/EXP(-B*A)*ABS(B)**C
      IF (IERROR.GT.0) GO TO 800
      ABSCIS(1) = C + ONE
      GO TO 180
  140 IF (C.LE.-ONE .OR. B.LE.ZERO) GO TO 800
      IERROR = 1
      PONORM = S14AAF((C+ONE)/TWO,IERROR)/B**((C+ONE)/TWO)
      IF (ITYPE.LT.0) PONORM = PONORM*B**(C/TWO)
      IF (IERROR.GT.0) GO TO 800
      GO TO 180
  160 IF (A+B.EQ.ZERO) GO TO 800
      IF (C.LE.-ONE .OR. D.LE.C+ONE) GO TO 800
      D = D - C - TWO
      PONORM = ONE/(TWO**(C+D+ONE))/(ABS(A+B)**(D+ONE))
      IF (ITYPE.LT.0) PONORM = PONORM*(TWO**(C+D+TWO))*(ABS(A+B)
     *                         **(D+TWO))
      GO TO 80
C
C       COMPUTE DIAGONAL AND OFF-DIAGONAL OF SYMMETRIC TRI-DIAGONAL
C         MATRIX WHICH HAS ABSCISSAE AS EIGENVALUES
C
  180 IF (NPNTS.EQ.1) GO TO 320
      DO 300 N = 2, NPNTS
         FN = N - 1
         GO TO (200,200,220,240,260,200) MITYPE
  200    TWNAPB = FN + FN + C + D
         ABSCIS(N) = (D+C)*(D-C)/(TWNAPB*(TWNAPB+TWO))
         CN = FOUR*(FN+C)*(FN+D)*FN/(TWNAPB**2*(TWNAPB+ONE))
         IF (N.GT.2) CN = CN*((C+D+FN)/(TWNAPB-ONE))
         GO TO 280
  220    ABSCIS(N) = ZERO
         CN = (FN+C*MOD(FN,TWO))**2/((FN+FN+C)**2-ONE)
         GO TO 280
  240    ABSCIS(N) = C + FN + FN + ONE
         CN = FN*(C+FN)
         GO TO 280
  260    ABSCIS(N) = ZERO
         CN = (FN+C*MOD(FN,TWO))/TWO
  280    WEIGHT(N) = SQRT(CN)
  300 CONTINUE
C
C        USE NAG ROUTINE TO FIND EIGENVALUES WHICH ARE ABSCISSAE
C
  320 IERROR = 1
      CALL F02AVF(NPNTS,X02AJF(),ABSCIS,WEIGHT,IERROR)
      IF (IERROR.GT.0) GO TO 760
C
C        LOOP TO DETERMINE WEIGHTS
C             EVALUATE EACH ORTHONORMAL POLYNOMIAL OF DEGREE
C         LESS THAN NPNTS AT ABSCIS(J) AND SUM SQUARES OF
C         RESULTS TO DETERMINE WEIGHT(J)
      IERROR = 0
      REALMX = X02ALF()
      DO 700 J = 1, NPNTS
C
C        INITIALISE INNER LOOP AND SCALE WEIGHT(J) AND ABSCIS(J)
C        DIVIDE EXPONENTIAL TERMS INTO FACTORS THAT DON'T UNDERFLOW
C
         WEIGHT(J) = ZERO
         Y = ABSCIS(J)
         PNA = ZERO
         CNO = ZERO
         NFAC = 0
         PNB = ONE/SQRT(PONORM)
         GO TO (340,340,360,400,420,440) MITYPE
  340    ABSCIS(J) = Y*(HALF*(B-A)) + (HALF*(A+B))
         IF (ITYPE.GT.0) GO TO 460
         PNB = PNB*(ONE-Y)**(C*HALF)*(ONE+Y)**(D*HALF)
         GO TO 460
  360    ABSCIS(J) = Y*(HALF*(B-A)) + (HALF*(A+B))
         IF (ITYPE.GT.0 .OR. C.EQ.ZERO) GO TO 460
         IF (Y.EQ.ZERO .AND. C.GT.ZERO) GO TO 660
         IF (C.GT.ZERO) GO TO 380
         IF (PONORM.GE.ONE) GO TO 380
         IF (ABS(Y).LE.(ONE/(REALMX*PONORM))**(-ONE/C)) GO TO 680
  380    PNB = PNB*ABS(Y)**(C*HALF)
         GO TO 460
  400    ABSCIS(J) = Y/B + A
         IF (ITYPE.GT.0) GO TO 460
         PNB = PNB*Y**(C*HALF)
         NFAC = INT(Y/LOG(HALF*REALMX)) + 1
         FACN = EXP(-HALF*Y/DBLE(NFAC))
         GO TO 460
  420    ABSCIS(J) = Y/SQRT(B) + A
         IF (ITYPE.GT.0) GO TO 460
         NFAC = INT(Y*Y/LOG(HALF*REALMX)) + 1
         FACN = EXP(-HALF*Y*Y/DBLE(NFAC))
         IF (C.EQ.ZERO) GO TO 460
         IF (Y.EQ.ZERO .AND. C.GT.ZERO) GO TO 660
         IF (Y.EQ.ZERO .AND. C.LT.ZERO) GO TO 680
         PNB = PNB*ABS(Y)**(C*HALF)
         GO TO 460
  440    ABSCIS(J) = TWO*(A+B)/(Y+ONE) - B
         IF (ITYPE.GT.0) GO TO 460
         PNB = PNB*(ONE-Y)**(C*HALF)*(ONE+Y)**(HALF*(D+TWO))
  460    WTSUM = PNB*PNB
         IF (NPNTS.EQ.1) GO TO 640
C
C          LOOP TO EVALUATE ORTHONORMAL POLYNOMIALS USING THREE
C           TERM RECURRENCE RELATION.
C
         DO 620 N = 2, NPNTS
            FN = N - 1
            GO TO (480,480,500,520,540,480) MITYPE
  480       TWNAPB = FN + FN + C + D
            BN = (D-C)/TWNAPB
            IF (N.GT.2) BN = BN*(C+D)/(TWNAPB-TWO)
            CN = FOUR*FN*(C+FN)*(D+FN)/(TWNAPB**2*(TWNAPB+ONE))
            IF (N.GT.2) CN = CN*((C+D+FN)/(TWNAPB-ONE))
            GO TO 560
  500       BN = ZERO
            CN = (FN+C*MOD(FN,TWO))**2/((FN+FN+C)**2-ONE)
            GO TO 560
  520       BN = C + FN + FN - ONE
            CN = FN*(FN+C)
            GO TO 560
  540       BN = ZERO
            CN = (FN+C*MOD(FN,TWO))/TWO
  560       SQRTCN = SQRT(CN)
            PNC = ((Y-BN)*PNB-CNO*PNA)/SQRTCN
            CNO = SQRTCN
            ABSPNC = ABS(PNC)
            IF (ABSPNC.LE.ONE) GO TO 580
            IF (ABSPNC.LE.REALMX/ABSPNC) GO TO 580
            IF (ITYPE.GT.0) GO TO 680
            IF (NFAC.LE.0) GO TO 680
            PNB = PNB*FACN
            PNC = PNC*FACN
            WTSUM = WTSUM*FACN*FACN
            NFAC = NFAC - 1
  580       PSQRD = PNC*PNC
            IF (WTSUM.LE.REALMX-PSQRD) GO TO 600
            IF (ITYPE.GT.0) GO TO 680
            IF (NFAC.LE.0) GO TO 680
            PNB = PNB*FACN
            PNC = PNC*FACN
            WTSUM = WTSUM*FACN*FACN
            PSQRD = PSQRD*FACN*FACN
            NFAC = NFAC - 1
  600       WTSUM = WTSUM + PSQRD
            PNA = PNB
            PNB = PNC
  620    CONTINUE
C
C          END LOOP FOR POLYNOMIAL EVALUATION
C
C Richard Brankin - NAG, Oxford - 26th July 1989
C replaced the following line ....
C
C  640    IF (NFAC.GT.0) WTSUM = WTSUM*FACN**(2*NFAC)
C
C so as not to get needless underflow to zero when powering up FACN
C for 0.0 < FACN << 1.0. The error was brought to light in a VAX
C double precision implementation when a user tried to compute modified
C Laguerre weights (ITYPE = -3) for more than 25 abscissae (N > 25).
C As a result, before the assignment in the above line
C   WTSUM = O(1.0e+38), FACN = O(1.0e-10), NFAC = 2
C WTSUM was assigned a value of 0.0 since O(1.0e-10)**4 underflows
C although WTSUM should have been assigned O(1.0e+2). This correction
C also applies for other values of ITYPE.
C
  640    IF (NFAC.GT.0) THEN
            DO 650 NBUG = 1, 2*NFAC
               WTSUM = WTSUM*FACN
  650       CONTINUE
         END IF
C
C End of correction
C
         IF (WTSUM.EQ.ZERO) GO TO 660
         WEIGHT(J) = ONE/WTSUM
         GO TO 700
  660    IERROR = 4
         WEIGHT(J) = REALMX
         GO TO 700
  680    IERROR = 5
  700 CONTINUE
C
C        END LOOP FOR WEIGHTS
C
C        REVERSE RATIONAL OR LAGUERRE POINTS
C
      IF ((MITYPE.NE.6 .OR. A+B.LT.ZERO)
     *    .AND. (MITYPE.NE.4 .OR. B.GT.ZERO)) GO TO 740
      NHALF = NPNTS/2
      IF (NHALF.LE.1) GO TO 740
      DO 720 J = 1, NHALF
         ISUB = NPNTS + 1 - J
         STORE = ABSCIS(J)
         ABSCIS(J) = ABSCIS(ISUB)
         ABSCIS(ISUB) = STORE
         STORE = WEIGHT(J)
         WEIGHT(J) = WEIGHT(ISUB)
         WEIGHT(ISUB) = STORE
  720 CONTINUE
C
C        ASSIGNMENT OF IFAIL PARAMETER
C
  740 IF ((ITYPE.EQ.-2 .OR. ITYPE.EQ.-4) .AND. MOD(NPNTS,2)
     *    .EQ.1 .AND. C.NE.ZERO) IERROR = 6
      GO TO 820
  760 IERROR = 1
      GO TO 820
  780 IERROR = 2
      GO TO 820
  800 IERROR = 3
  820 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
