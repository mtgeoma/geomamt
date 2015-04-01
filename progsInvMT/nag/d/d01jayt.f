      SUBROUTINE D01JAY(R,RR,RO,WEIGHT,N,NNMIN,RADIUS,ITRANS,CONST,RMID)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     D01JAY COMPUTES THE TRANSFORMATION- RO AND WEIGHT,
C     CORRESPONDING TO R.
C     THE VALUE OF ITRANS DETERMINES WHICH TRANSFORMATION IS USED.
C     FOR FULL DESCRIPTION OF THE PARAMETERS, SEE D01JAF.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONST, R, RADIUS, RMID, RO, RR, WEIGHT
      INTEGER           ITRANS, N, NNMIN
C     .. Scalars in Common ..
      DOUBLE PRECISION  PI
C     .. Local Scalars ..
      DOUBLE PRECISION  ARG, CHARG, CHISH, EARG, EARGI, ESH, ESHI, PI2,
     *                  PI4, SHARG, Z
C     .. Intrinsic Functions ..
      INTRINSIC         EXP
C     .. Common blocks ..
      COMMON            /CD01JA/PI
C     .. Executable Statements ..
      PI2 = PI/2.0D0
      PI4 = PI2/2.0D0
      IF (R.EQ.0.D0) GO TO 160
      IF (ITRANS.GE.5) GO TO 80
      IF (ITRANS.EQ.1 .OR. ITRANS.EQ.3) GO TO 40
      ARG = R
   20 EARG = EXP(ARG)
      EARGI = 1.D0/EARG
      CHARG = (EARG+EARGI)*0.5D0
      SHARG = PI4*(EARG-EARGI)
      ESH = EXP(SHARG)
      ESHI = 1.D0/ESH
      CHISH = 2.D0/(ESH+ESHI)
      IF (ITRANS.GE.5) GO TO 120
      IF (ITRANS.GE.3) GO TO 60
      Z = (ESH-ESHI)/(ESH+ESHI)
      RO = Z*RADIUS
      WEIGHT = (Z/R)**NNMIN*PI2*CHARG*CHISH*CHISH
      IF (ITRANS.EQ.2) RETURN
      WEIGHT = WEIGHT*CONST*(1.D0+RR)/(1.D0-RR)**2
      RETURN
   40 ARG = CONST*R/(1.D0-RR)
      GO TO 20
   60 Z = ESHI*CHISH
      RO = Z*RADIUS
      WEIGHT = ((1.D0-Z)/R)**NNMIN*PI2*CHARG*CHISH*CHISH
      IF (ITRANS.EQ.4) RETURN
      WEIGHT = WEIGHT*CONST*(1.D0+RR)/(1.D0-RR)**2
      RETURN
   80 IF (ITRANS.EQ.6) GO TO 100
      ARG = CONST*(R-RMID)/R/(1.D0-R)
      GO TO 20
  100 ARG = CONST*(R-RMID)*(R+1.D0)/R
      GO TO 20
  120 Z = ESH*CHISH*0.5D0
      RO = Z*RADIUS
      WEIGHT = (Z/R)**NNMIN*CHARG*CHISH*CHISH
      IF (ITRANS.EQ.6) GO TO 140
      WEIGHT = WEIGHT*PI4*CONST*(RR-2.D0*R*RMID+RMID)/(RR*(1.D0-R)**2)
      RETURN
  140 WEIGHT = WEIGHT*PI4*CONST*((2.D0*R-RMID+1.D0)/R-(R-RMID)*(R+1.D0)
     *         /RR)
      RETURN
  160 GO TO (180,200,180,200,220,220) ITRANS
  180 WEIGHT = (CONST*PI2)**N
      RO = 0.D0
      IF (ITRANS.EQ.3) RO = RADIUS
      RETURN
  200 WEIGHT = PI2**N
      RO = 0.D0
      IF (ITRANS.EQ.4) RO = RADIUS
      RETURN
  220 WEIGHT = 0.D0
      RO = 0.D0
      RETURN
      END
