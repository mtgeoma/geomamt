      SUBROUTINE E01DAY(NORDER,NKNOTS,M,X,W,IW1,LW,IFNZWT,MDNZWT,LAMBDA,
     *                  LLMBDA)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE B1K       PROVIDE TRIAL KNOT SET, GIVEN
C     ==============       ABSCISSAE AND WEIGHTS
C                          - BASE VERSION
C
C     CREATED 13 10 80.  UPDATED 16 03 83.  RELEASE 00/07
C
C     AUTHORS ... MAURICE G. COX, PAULINE E. M. CURTIS
C                 AND MICHAEL A. SINGER.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1980-1982
C
C     **********************************************************
C
C     B1K.  GIVEN A SET OF ABSCISSA (INDEPENDENT VARIABLE)
C     VALUES AND CORRESPONDING WEIGHTS, DETERMINES A SET OF
C     INTERIOR KNOTS SUCH THAT THERE IS AN APPROXIMATELY
C     EQUAL NUMBER OF ABSCISSAE WITH NONZERO WEIGHTS IN
C     EACH KNOT INTERVAL (EXCEPT IN THE FIRST AND LAST
C     INTERVALS WHERE THERE ARE APPROXIMATELY  NORDER/2
C     TIMES AS MANY VALUES).
C
C     INPUT PARAMETERS
C        NORDER   ORDER (DEGREE + 1) OF SPLINE  S
C        NKNOTS   NUMBER OF INTERIOR KNOTS
C        M        NUMBER OF INDEPENDENT VARIABLE VALUES
C        X        INDEPENDENT VARIABLE VALUES
C        W        WEIGHTS (INVERSELY PROPORTIONAL TO
C                    EXPECTED ERRORS IN F-VALUES)
C        IW1      INDEX INCREMENT OF  W
C        LW       DIMENSION OF  W
C        IFNZWT   INDEX OF FIRST DATA POINT WITH NONZERO
C                    WEIGHT
C        MDNZWT   NUMBER OF DATA POINTS WITH DISTINCT
C                    ABSCISSA VALUES AND NONZERO WEIGHT
C
C     OUTPUT (AND ASSOCIATED) PARAMETERS
C        LAMBDA   INTERIOR KNOTS
C        LLMBDA   DIMENSION OF  LAMBDA.
C                    .GE. MAX(NKNOTS, 1).
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IFNZWT, IW1, LLMBDA, LW, M, MDNZWT, NKNOTS,
     *                  NORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  LAMBDA(LLMBDA), W(LW), X(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, C1, DIST, HALF, ONE, PPJ, RI1, RM, RQ,
     *                  SI, XI, XINEXT, ZERO
      INTEGER           I, I1, I1L, IP1, ISTRT, IW, J, K, P
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MOD
C     .. Data statements ..
C
      DATA              ZERO, HALF, ONE/0.0D+0, 0.5D+0, 1.0D+0/
C     .. Executable Statements ..
C
      IF (NKNOTS.LE.0) GO TO 120
      RM = MDNZWT
      RQ = NKNOTS + NORDER
      P = NORDER/2
      IF (MOD(NORDER,2).NE.0) GO TO 20
C
C        NORDER  EVEN
C
      DIST = (RM-ONE)/(RQ-ONE)
      C1 = ONE - DIST
      GO TO 40
C
C        NORDER  ODD
C
   20 DIST = RM/RQ
      C1 = HALF
C
C     NOW FORM THE KNOTS
C
   40 I1L = 0
      I = IFNZWT
      XINEXT = X(I)
      IW = 1 + (IFNZWT-1)*IW1
      DO 100 J = 1, NKNOTS
         PPJ = P + J
         SI = C1 + PPJ*DIST
         I1 = INT(SI)
         RI1 = I1
         ALPHA = SI - RI1
         ISTRT = I1L + 1
         I1L = I1
         DO 80 K = ISTRT, I1
            XI = XINEXT
            IP1 = I + 1
C
C           DETERMINE  XINEXT,  THE NEXT X-VALUE WITH
C           NONZERO WEIGHT THAT IS DISTINCT FROM THE
C           PREVIOUS SUCH X-VALUE
C
            DO 60 I = IP1, M
               XINEXT = X(I)
               IW = IW + IW1
               IF (W(IW).NE.ZERO .AND. XINEXT.NE.XI) GO TO 80
   60       CONTINUE
   80    CONTINUE
         LAMBDA(J) = (ONE-ALPHA)*XI + ALPHA*XINEXT
  100 CONTINUE
  120 CONTINUE
      RETURN
C
C     END E01DAY
C
      END
