      SUBROUTINE G01AGY(X,Y,XMIN,XSTEP,NSTEPX,YMIN,YSTEP,NSTEPY,NOBS,
     *                  MAXA,MAXB,KEY,MAXYA,MAXYB)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     G01AGY DOES THE WORK OF PLOTTING A GRAPH
C     REQUESTED BY A CALL TO G01AGF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMIN, XSTEP, YMIN, YSTEP
      INTEGER           MAXA, MAXB, MAXYA, MAXYB, NOBS, NSTEPX, NSTEPY
C     .. Array Arguments ..
      DOUBLE PRECISION  X(NOBS), Y(NOBS)
      INTEGER           KEY(NOBS)
C     .. Local Scalars ..
      DOUBLE PRECISION  XDIFF, YAXIS, YDIFF, ZMP
      INTEGER           I, IAY, IFOR, II, IR, ISP1, ISP2, ISPACE, IXAX,
     *                  IXZERO, IYZERO, J, J1, K, MAXF, MAXYF, NCOUNT,
     *                  NDEV, NEND, NSTPX2
      LOGICAL           LXE, LYE
      CHARACTER*135     REC
C     .. Local Arrays ..
      DOUBLE PRECISION  XAXIS(28)
      INTEGER           IA(133)
      CHARACTER         IAC(133), ICODE(39), IFMA(16), IFMB(11),
     *                  IFMO(17), IXLINE(133), NUMB(10)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, AINT, INT, MIN, MOD, DBLE, SIGN
C     .. Data statements ..
      DATA              IFMO(1), IFMO(2), IFMO(3), IFMO(4), IFMO(5),
     *                  IFMO(6), IFMO(7), IFMO(8), IFMO(9), IFMO(10),
     *                  IFMO(11), IFMO(12), IFMO(13), IFMO(14),
     *                  IFMO(15), IFMO(16), IFMO(17)/'(', ' ', ' ', 'X',
     *                  ',', '2', '8', '(', 'F', ' ', '.', ' ', ',',
     *                  ' ', 'X', ')', ')'/
      DATA              IFMA(1), IFMA(2), IFMA(3), IFMA(4), IFMA(5),
     *                  IFMA(6), IFMA(7), IFMA(8), IFMA(9), IFMA(10),
     *                  IFMA(11), IFMA(12), IFMA(13), IFMA(14),
     *                  IFMA(15), IFMA(16)/'(', '1', 'X', ',', 'F', '1',
     *                  '1', '.', ' ', ',', ' ', ' ', ' ', 'A', '1',
     *                  ')'/
      DATA              IFMB(1), IFMB(2), IFMB(3), IFMB(4), IFMB(5),
     *                  IFMB(6), IFMB(7), IFMB(8), IFMB(9), IFMB(10),
     *                  IFMB(11)/'(', '1', '2', 'X', ',', ' ', ' ', ' ',
     *                  'A', '1', ')'/
      DATA              ICODE(1), ICODE(2), ICODE(3), ICODE(4),
     *                  ICODE(5), ICODE(6), ICODE(7), ICODE(8),
     *                  ICODE(9), ICODE(10), ICODE(11), ICODE(12),
     *                  ICODE(13), ICODE(14), ICODE(15), ICODE(16),
     *                  ICODE(17), ICODE(18), ICODE(19), ICODE(20),
     *                  ICODE(21), ICODE(22), ICODE(23), ICODE(24),
     *                  ICODE(25), ICODE(26), ICODE(27), ICODE(28),
     *                  ICODE(29), ICODE(30), ICODE(31), ICODE(32),
     *                  ICODE(33), ICODE(34), ICODE(35), ICODE(36),
     *                  ICODE(37), ICODE(38), ICODE(39)/'.', '+', ' ',
     *                  '1', '2', '3', '4', '5', '6', '7', '8', '9',
     *                  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
     *                  'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
     *                  'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '*'/
      DATA              NUMB(1), NUMB(2), NUMB(3), NUMB(4), NUMB(5),
     *                  NUMB(6), NUMB(7), NUMB(8), NUMB(9),
     *                  NUMB(10)/'0', '1', '2', '3', '4', '5', '6', '7',
     *                  '8', '9'/
C     .. Executable Statements ..
      CALL X04ABF(0,NDEV)
      ZMP = 1.0D0 + 2.0D0*X02AJF()
      XDIFF = 0.0D0
      YDIFF = 0.0D0
      IXZERO = 0
      IYZERO = 0
C     ADJUST SPAN OF X-AXIS TO START AND END ON GRADATIONS
      IF (XMIN.EQ.0.0D0) GO TO 40
      IF (ABS((XMIN/XSTEP)-AINT(ZMP*(XMIN/XSTEP)+SIGN(0.5D0,XMIN)))
     *    .GE.0.00005D0) XDIFF = MOD(XMIN,XSTEP)
      IF (XDIFF.EQ.0.0D0) GO TO 20
      IF (XDIFF.LT.0.0D0) XDIFF = XSTEP + XDIFF
      XMIN = XMIN - XDIFF
      NSTEPX = NSTEPX + 1
   20 IF (XMIN.EQ.0.0D0) GO TO 40
      IXZERO = -INT(ZMP*(XMIN/XSTEP)+SIGN(0.5D0,XMIN))
C     ADJUST SPAN OF Y-AXIS TO START AND END ON GRADATIONS
   40 IF (YMIN.EQ.0.0D0) GO TO 80
      IF (ABS((YMIN/YSTEP)-AINT(ZMP*(YMIN/YSTEP)+SIGN(0.5D0,YMIN)))
     *    .GE.0.00005D0) YDIFF = MOD(YMIN,YSTEP)
      IF (YDIFF.EQ.0.0D0) GO TO 60
      IF (YDIFF.LT.0.0D0) YDIFF = YSTEP + YDIFF
      YMIN = YMIN - YDIFF
      NSTEPY = NSTEPY + 1
   60 IF (YMIN.EQ.0.0D0) GO TO 80
      IYZERO = -INT(ZMP*(YMIN/YSTEP)+SIGN(0.5D0,YMIN))
C     DRAW BOX AROUND PLOT AND AXES
   80 NSTPX2 = NSTEPX + 2
      DO 100 I = 1, NSTPX2
         IXLINE(I) = ICODE(1)
         IF (MOD(ABS(IXZERO-I+1),5).NE.0) GO TO 100
         IXLINE(I) = ICODE(2)
         NEND = (I-1)/5 + 1
         XAXIS(NEND) = XMIN + XSTEP*DBLE(I-1)
  100 CONTINUE
C     SET UP FORMATS FOR ANNOTATION OF Y AXIS AND GRAPH
      IFOR = MOD(IXZERO,5)
      MAXYF = MAXYA + MAXYB + 2
      IF (MAXYF.GT.11) GO TO 120
      IFMA(9) = NUMB(MAXYB+1)
      LYE = .FALSE.
      GO TO 140
  120 IFMA(5) = ICODE(17)
      IFMA(9) = NUMB(4)
      LYE = .TRUE.
  140 J1 = NSTPX2/100
      J = J1 + 1
      IFMA(11) = NUMB(J)
      IFMB(6) = NUMB(J)
      J = NSTPX2/10 + 1 - J1*10
      IFMA(12) = NUMB(J)
      IFMB(7) = NUMB(J)
      J = MOD(NSTPX2,10) + 1
      IFMA(13) = NUMB(J)
      IFMB(8) = NUMB(J)
      IF (MOD(ABS(IYZERO-(NSTEPY+1)),5).EQ.0) GO TO 160
      WRITE (REC,FMT=IFMB) (IXLINE(I),I=1,NSTPX2)
      CALL X04BAF(NDEV,REC)
      GO TO 180
  160 YAXIS = YMIN + YSTEP*DBLE(NSTEPY+1)
      WRITE (REC,FMT=IFMA) YAXIS, (IXLINE(I),I=1,NSTPX2)
      CALL X04BAF(NDEV,REC)
  180 NCOUNT = 1
      DO 380 J = 1, NSTEPY
         K = NSTEPY - J + 1
         DO 200 I = 1, NSTEPX
            IA(I) = 3
  200    CONTINUE
  220    IF (NCOUNT.GT.NOBS) GO TO 280
         IAY = INT(ZMP*((Y(NCOUNT)-YMIN)/YSTEP)+0.5D0)
         IF (IAY.EQ.0) IAY = 1
         IF (IAY.LT.K) GO TO 280
         II = KEY(NCOUNT)
         IXAX = INT(ZMP*((X(II)-XMIN)/XSTEP)+0.5D0)
         IF (IXAX.EQ.0) IXAX = 1
         IF (IA(IXAX).LE.3) GO TO 240
         IA(IXAX) = IA(IXAX) + 1
         GO TO 260
  240    IA(IXAX) = 4
  260    NCOUNT = NCOUNT + 1
         GO TO 220
  280    DO 340 I = 1, NSTEPX
            IF (IA(I).NE.3) GO TO 320
            IF (K.NE.IYZERO) GO TO 300
            IA(I) = 1
            IF (MOD(ABS(IXZERO-I),5).EQ.0) IA(I) = 2
            GO TO 320
  300       IF (I.NE.IXZERO) GO TO 320
            IA(I) = 1
            IF (MOD(ABS(IYZERO-K),5).EQ.0) IA(I) = 2
  320       II = MIN(IA(I),39)
            IAC(I) = ICODE(II)
  340    CONTINUE
         IF (MOD(ABS(IYZERO-K),5).EQ.0) GO TO 360
         WRITE (REC,FMT=IFMB) ICODE(1), (IAC(I),I=1,NSTEPX), ICODE(1)
         CALL X04BAF(NDEV,REC)
         GO TO 380
  360    YAXIS = YMIN + YSTEP*DBLE(K)
         WRITE (REC,FMT=IFMA) YAXIS, ICODE(2), (IAC(I),I=1,NSTEPX),
     *     ICODE(2)
         CALL X04BAF(NDEV,REC)
  380 CONTINUE
      IF (MOD(ABS(IYZERO),5).EQ.0) GO TO 400
      WRITE (REC,FMT=IFMB) (IXLINE(I),I=1,NSTPX2)
      CALL X04BAF(NDEV,REC)
      GO TO 420
  400 YAXIS = YMIN
      WRITE (REC,FMT=IFMA) YAXIS, (IXLINE(I),I=1,NSTPX2)
      CALL X04BAF(NDEV,REC)
C     SET UP FORMATS FOR ANNOTATION OF X AXIS
  420 IFOR = MOD(IXZERO,5)
      IF (IFOR.LT.0) IFOR = IFOR + 5
      ISP1 = 11 + IFOR
      IF (MAXA+MAXB+2.LE.9) ISP1 = ISP1 - MAXA
      IR = NSTEPX/10 + 1
      I = ISP1/10 + 1
      IFMO(2) = NUMB(I)
      I = MOD(ISP1,10) + 1
      IFMO(3) = NUMB(I)
      MAXF = MAXA + MAXB + 2
      IF (MAXF.GT.9) GO TO 440
      ISPACE = 10 - MAXF
      LXE = .FALSE.
      IFMO(10) = NUMB(MAXF+1)
      IFMO(12) = NUMB(MAXB+1)
      IFMO(14) = NUMB(ISPACE+1)
      GO TO 460
  440 IFMO(9) = ICODE(17)
      IFMO(10) = NUMB(10)
      IFMO(12) = NUMB(3)
      IFMO(14) = NUMB(2)
      LXE = .TRUE.
  460 WRITE (REC,FMT=IFMO) (XAXIS(II),II=1,NEND,2)
      CALL X04BAF(NDEV,REC)
      ISP2 = ISP1 + 5
      I = ISP2/10 + 1
      IFMO(2) = NUMB(I)
      I = MOD(ISP2,10) + 1
      IFMO(3) = NUMB(I)
      WRITE (REC,FMT=IFMO) (XAXIS(II),II=2,NEND,2)
      CALL X04BAF(NDEV,REC)
      IF (LXE) IFMO(9) = ICODE(18)
      IF (LYE) IFMA(5) = ICODE(18)
      RETURN
      END
