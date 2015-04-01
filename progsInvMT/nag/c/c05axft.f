      SUBROUTINE C05AXF(X,FX,TOL,IR,SCALE,C,IND,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-301 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     USES METHOD OF OF SWIFT AND LINDFIELD,
C     C.J.VOL 21. MINK=C(1),MAXK=C(2).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05AXF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FX, SCALE, TOL, X
      INTEGER           IFAIL, IND, IR
C     .. Array Arguments ..
      DOUBLE PRECISION  C(26)
C     .. Local Scalars ..
      DOUBLE PRECISION  AB, OLDTH, REL, SE, SF
      INTEGER           I
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      IF (TOL.LE.0.0D0) GO TO 440
      IF (X+SCALE.EQ.X) GO TO 480
      IF (IND.LT.-1 .OR. IND.GT.4) GO TO 460
      IF (IR.LT.0 .OR. IR.GT.2) GO TO 440
      AB = 1.0D0
      REL = 1.0D0
      IF (IR.EQ.1) REL = 0.0D0
      IF (IR.EQ.2) AB = 0.0D0
      I = IND + 2
      GO TO (40,460,20,40,80,260) I
   20 IND = 2
      RETURN
   40 IF (FX.EQ.0.0D0) GO TO 580
      C(1) = 3.0D0
      C(2) = 8.0D0
      C(4) = 1.0D0
      C(6) = 1.0D0
      C(7) = 0.5D0
      C(8) = FX
      C(22) = X
      C(25) = X02AJF()
      C(19) = SQRT(C(25))
      C(5) = 1.0D0 - 0.5D0*C(19)
      C(23) = C(19)*SCALE
      C(20) = 0.D0
   60 X = C(22) + C(23)
      IF (X.EQ.C(22)) GO TO 100
      IND = 3
      RETURN
   80 IF (FX.EQ.0.0D0) GO TO 580
      C(24) = FX - C(8)
      X = C(22)
      IF (ABS(C(24)).GT.100.0D0*C(25)*MAX(ABS(FX),ABS(C(8))))
     *    GO TO 120
  100 C(23) = 10.0D0*C(23)
      IF (ABS(C(23)).GT.ABS(SCALE)) GO TO 480
      GO TO 60
  120 C(12) = C(23)/C(24)
C
C     CONTINUATION OUTER LOOP
C
  140 C(14) = SQRT(TOL)
      IF (C(5).EQ.0.0D0) C(14) = TOL
      C(3) = 0.0D0
      C(9) = X
      C(13) = C(12)
      C(15) = C(6) - C(5)
      IF (C(15).LT.100.0D0*C(25)) GO TO 500
      C(21) = C(8)*C(5)
      C(10) = C(15)*C(8)
      C(11) = C(10)
      C(26) = 0.0D0
C
C     CONTINUATION INNER LOOP
C
  160 C(17) = C(12)*C(10)
      SF = C(26)
      C(26) = 0.0D0
      IF (ABS(C(17)).LT.C(14)*MAX(REL*ABS(X),AB)) C(26) = SF + 1.0D0
      IF (C(26).GT.0.0D0 .AND. ABS(C(10)).LT.ABS(C(11)) .AND. C(3)
     *    .GT.1.5D0) GO TO 280
      IF (ABS(C(3)-C(2)).GT.0.5D0) GO TO 240
      IF (C(26).GT.C(2)-2.5D0) GO TO 280
      SE = (C(10)+C(21))/C(8)
      IF (SE.GT.C(5) .AND. SE.LT.C(6)) GO TO 200
      C(12) = C(13)
      C(7) = 0.5D0*C(7)
      IF (C(6)-C(5).GT.100.0D0*C(25)) GO TO 180
      IF (C(4).LT.1.5D0) GO TO 500
      GO TO 400
  180 C(5) = 0.5D0*(C(5)+C(6))
      X = C(9)
      GO TO 140
  200 C(21) = SE*C(8)
      C(15) = C(6) - SE
      IF (C(5).NE.0.0D0) GO TO 220
      IF (C(15).LT.100.0D0*C(25)) GO TO 500
  220 C(5) = SE
      GO TO 280
  240 X = X - C(17)
      IND = 4
      RETURN
  260 IF (FX.EQ.0.0D0) GO TO 580
      SE = FX - C(21)
      SF = SE - C(10)
      IF (ABS(SF).LE.2.D0*C(25)*MAX(ABS(SE),ABS(C(10)),ABS(C(21))))
     *    GO TO 280
      C(10) = SE
      C(3) = C(3) + 1.0D0
      IF (SF.EQ.0.0D0) GO TO 160
      C(12) = -C(17)/SF
      GO TO 160
C
C     SUCCESSFUL INNER LOOP
C
  280 IF (C(5).EQ.0.0D0) GO TO 580
      C(17) = C(9) - X
      OLDTH = C(6)
      C(6) = C(5)
      IF (ABS(C(4)-1.0D0).GT.0.5D0) GO TO 300
      C(5) = 1.0D0 - C(19)
      GO TO 420
  300 SF = C(17)/C(15)
      SE = (C(18)/C(16)-SF)/(C(15)+C(16))
      IF (ABS(SF)*C(19).GT.ABS(C(15)) .OR. SE.EQ.0.D0 .OR. ABS(SE)*C(19)
     *    .GT.ABS(C(15)+C(16))) GO TO 320
      IF (C(3).LE.C(1)+0.5D0 .AND. C(4).GT.2.5D0) C(7) = 2.0D0*C(7)
      IF (ABS(SE)*C(5).GT.C(7)*ABS(SF-SE*C(15))) GO TO 360
      C(5) = 0.D0
      GO TO 420
  320 IF (C(5).LE.0.9D0) GO TO 340
      C(5) = 0.9D0
      GO TO 380
  340 C(5) = C(5)*0.5D0
      GO TO 380
  360 C(5) = C(5) - C(7)*ABS(SF/SE-C(15))
      IF (C(5).LT.0.0D0) C(5) = 0.0D0
      IF (C(5).EQ.0.0D0) GO TO 420
  380 CONTINUE
      IF (OLDTH-C(5).GE.100.0D0*C(25) .AND. C(5).GE.C(19)) GO TO 420
  400 CONTINUE
      IF (C(4).GT.2.5D0) GO TO 500
      C(19) = 0.1D0*C(19)
      IF (C(19).LT.100.0D0*C(25)) GO TO 540
      C(6) = 1.0D0
      C(5) = 1.0D0 - 0.5D0*C(19)
      C(20) = 0.0D0
      X = C(9)
      C(4) = 1.0D0
      C(7) = 0.5D0
      GO TO 140
  420 C(16) = C(15)
      C(18) = C(17)
      C(4) = C(4) + 1.0D0
      C(20) = 1.0D0
      IF (C(5).EQ.0.0D0) C(20) = 2.0D0
      GO TO 140
C     INPUT ERROR
  440 I = 1
      GO TO 600
C     WRONG IND ON ENTRY
  460 I = 2
      GO TO 600
C     WRONG SCALE
  480 I = 3
      GO TO 600
C     SINGULARITY ON CONTINUATION PATH
  500 IF (C(20).NE.1.0D0) GO TO 520
      I = 4
      GO TO 600
C     CANNOT GET STARTED
  520 IF (C(20).NE.0.0D0) GO TO 560
  540 I = 5
      GO TO 600
C     CANNOT FINISH
  560 I = 6
      GO TO 600
C     FINISHED
  580 I = 0
  600 IND = 0
      IFAIL = P01ABF(IFAIL,I,SRNAME,0,P01REC)
      RETURN
      END
