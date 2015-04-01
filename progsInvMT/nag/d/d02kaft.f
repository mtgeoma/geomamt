      SUBROUTINE D02KAF(XL,XR,COEFFN,BCOND,K,TOL,ELAM,DELAM,MONIT,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8A REVISED. IER-257 (AUG 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02KAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DELAM, ELAM, TOL, XL, XR
      INTEGER           IFAIL, K
C     .. Array Arguments ..
      DOUBLE PRECISION  BCOND(3,2)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFFN, MONIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, HMX, LAMDA, MINSC, ONE, PI, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      INTEGER           I, IND, IND1, MAXFUN, MAXIT
C     .. Local Arrays ..
      DOUBLE PRECISION  HMAX(2,5), XPOINT(5)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02KAZ, D02KDF
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, HMX, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      IF (K.GE.0 .AND. TOL.GT.0.D0) GO TO 20
      IND = 1
      GO TO 300
   20 XPOINT(1) = XL
      XPOINT(2) = XL
      XPOINT(3) = .5D0*(XL+XR)
      XPOINT(4) = XR
      XPOINT(5) = XR
      HMAX(1,1) = 0.D0
      HMAX(1,2) = 0.D0
      IND1 = 1
      MAXFUN = 0
      MAXIT = 15
      DO 40 I = 1, 2
         YL(I) = BCOND(I,1)
         YR(I) = BCOND(I,2)
   40 CONTINUE
C
      CALL D02KDF(XPOINT,5,COEFFN,D02KAZ,K,TOL,ELAM,DELAM,HMAX,MAXIT,
     *            MAXFUN,MONIT,IND1)
C
      IND = 0
      IF (IND1.EQ.0) GO TO 280
      IF (IND1.GT.0 .AND. IND1.LT.13)
     *    GO TO (60,80,100,120,140,60,160,180,200,
     *    220,240,260) IND1
   60 IND = 12
      GO TO 280
   80 IND = 2
      GO TO 280
  100 IND = 3
      GO TO 280
  120 IND = 4
      GO TO 280
  140 IND = 5
      GO TO 280
  160 IND = 6
      GO TO 280
  180 IND = 7
      GO TO 280
  200 IND = 8
      GO TO 280
  220 IND = 9
      GO TO 280
  240 IND = 10
      GO TO 280
  260 IND = 11
  280 BCOND(3,1) = HMAX(1,4)
      BCOND(3,2) = HMAX(1,5)
  300 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
