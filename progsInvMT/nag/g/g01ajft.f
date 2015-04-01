      SUBROUTINE G01AJF(X,N,NSTEPX,NSTEPY,ITYPE,ISPACE,XMIN,XMAX,XSTEP,
     *                  N1,MULTY,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G01AJF PRINTS A HISTOGRAM OF A VECTOR OF REAL DATA,
C     WITH A CHOSEN NUMBER OF CHARACTER POSITIONS IN EACH DIRECTION.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01AJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN, XSTEP
      INTEGER           IFAIL, ISPACE, ITYPE, MULTY, N, N1, NSTEPX,
     *                  NSTEPY
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XI, XMN, XMX
      INTEGER           I, IFA, II, IS1, IS2, IS3, J, MAXF, MINF, NDEV,
     *                  NN, NSPX2, NSX2
      CHARACTER*1       IBLANK, IDASH, IDOT, ISTAR
      CHARACTER*135     REC
C     .. Local Arrays ..
      INTEGER           IFREQ(99)
      CHARACTER*1       IF1(10), IF2(14), IF3(20), IF4(22), IH1(10),
     *                  IH2(10), K(10), LINE(101), P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DBLE
C     .. Data statements ..
      DATA              IBLANK, IDOT, ISTAR, IDASH/' ', '.', '*', '-'/
      DATA              K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8),
     *                  K(9), K(10)/'0', '1', '2', '3', '4', '5', '6',
     *                  '7', '8', '9'/
      DATA              IF1(1), IF1(2), IF1(3), IF1(4), IF1(5), IF1(6),
     *                  IF1(7), IF1(8), IF1(9), IF1(10)/'(', ' ', ' ',
     *                  'X', ',', '1', '0', 'A', '1', ')'/
      DATA              IH1(1), IH1(2), IH1(3), IH1(4), IH1(5), IH1(6),
     *                  IH1(7), IH1(8), IH1(9), IH1(10)/'F', 'R', 'E',
     *                  'Q', 'U', 'E', 'N', 'C', 'Y', ' '/
      DATA              IH2(1), IH2(2), IH2(3), IH2(4), IH2(5), IH2(6),
     *                  IH2(7), IH2(8), IH2(9), IH2(10)/'C', 'U', 'M',
     *                  '.', ' ', 'F', 'R', 'E', 'Q', '.'/
      DATA              IF2(1), IF2(2), IF2(3), IF2(4), IF2(5), IF2(6),
     *                  IF2(7), IF2(8), IF2(9), IF2(10), IF2(11),
     *                  IF2(12), IF2(13), IF2(14)/'(', ' ', ' ', 'X',
     *                  ',', '7', 'X', ',', ' ', ' ', ' ', 'A', '1',
     *                  ')'/
      DATA              IF3(1), IF3(2), IF3(3), IF3(4), IF3(5), IF3(6),
     *                  IF3(7), IF3(8), IF3(9), IF3(10), IF3(11),
     *                  IF3(12), IF3(13), IF3(14), IF3(15), IF3(16),
     *                  IF3(17), IF3(18), IF3(19), IF3(20)/'(', ' ',
     *                  ' ', 'X', ',', 'I', '5', ',', '2', 'X', ',',
     *                  ' ', ' ', ' ', 'A', '1', ',', 'I', '5', ')'/
      DATA              IF4(1), IF4(2), IF4(3), IF4(4), IF4(5), IF4(6),
     *                  IF4(7), IF4(8), IF4(9), IF4(10), IF4(11),
     *                  IF4(12), IF4(13), IF4(14), IF4(15), IF4(16),
     *                  IF4(17), IF4(18), IF4(19), IF4(20), IF4(21),
     *                  IF4(22)/'(', ' ', ' ', 'X', ',', '3', 'X', ',',
     *                  'F', '8', '.', '2', ',', ' ', ' ', 'X', ',',
     *                  'F', '8', '.', '2', ')'/
C     .. Executable Statements ..
C
C     CHECK PARAMETERS OF PLOT SIZE AND POSITIONING
C
      IFA = 1
      IF (N.LT.1) GO TO 380
      IF (NSTEPX.LT.10) NSTEPX = 10
      IF (NSTEPX.GT.99) NSTEPX = 99
      IF (NSTEPY.LT.10) NSTEPY = 10
      IF (NSTEPY.GT.99) NSTEPY = 99
      IF (ISPACE.LT.0) ISPACE = 0
      IF (ISPACE+NSTEPX+14.GT.120) ISPACE = 0
C
C     CLEAR VECTOR OF FREQUENCY COUNTS
C
      DO 20 I = 1, NSTEPX
         IFREQ(I) = 0
   20 CONTINUE
      IF (XMIN.GE.XMAX) GO TO 60
C
C     MIN AND MAX SUPPLIED - COMPUTE X STEP LENGTH
C
      XSTEP = (XMAX-XMIN)/DBLE(NSTEPX)
C
C     FORM FREQUENCY COUNT FOR EACH CLASS INTERVAL
C
      N1 = 0
      DO 40 I = 1, N
         XI = X(I)
         IF (XI.LT.XMIN .OR. XI.GT.XMAX) GO TO 40
         J = (XI-XMIN)/XSTEP
         IF (XI.LT.XMAX) J = J + 1
         IFREQ(J) = IFREQ(J) + 1
         N1 = N1 + 1
   40 CONTINUE
      XMX = XMAX
      XMN = XMIN
      GO TO 140
C
C     MIN AND MAX NOT SUPPLIED - FIND THEM
C
   60 XMIN = X(1)
      XMAX = XMIN
      IF (N.EQ.1) GO TO 120
      DO 80 I = 2, N
         XMIN = MIN(XMIN,X(I))
         XMAX = MAX(XMAX,X(I))
   80 CONTINUE
      IF (XMIN.EQ.XMAX) GO TO 120
C
C     DATA HAS NON-ZERO RANGE - FORM FREQUENCY COUNTS
C
      XSTEP = (XMAX-XMIN)/DBLE(NSTEPX)
      N1 = N
      DO 100 I = 1, N
         XI = X(I)
         J = (XI-XMIN)/XSTEP
         IF (XI.LT.XMAX) J = J + 1
         IFREQ(J) = IFREQ(J) + 1
  100 CONTINUE
      XMN = XMIN
      XMX = XMAX
      GO TO 140
C
C     DATA HAS ZERO RANGE - PUT ALL OBSERVATIONS IN CENTRE
C
  120 NSX2 = NSTEPX/2
      IFREQ(NSX2) = N
      XSTEP = 0.01D0
      XMN = XMIN - 0.01D0*DBLE(NSX2-1)
      XMX = XMAX + 0.01D0*DBLE(NSTEPX-NSX2)
C
C     FREQUENCY COUNTS ARE NOW FORMED
C     FORM CUMULATIVE FREQUENCIES IF REQUESTED
C
  140 IF (ITYPE.EQ.0) GO TO 180
      DO 160 I = 2, NSTEPX
         IFREQ(I) = IFREQ(I) + IFREQ(I-1)
  160 CONTINUE
C
C     FIND MAXIMUM FREQUENCY
C
      MAXF = IFREQ(NSTEPX)
      IF (ITYPE.NE.0) GO TO 220
  180 MAXF = IFREQ(1)
      DO 200 I = 2, NSTEPX
         MAXF = MAX(MAXF,IFREQ(I))
  200 CONTINUE
C
C     COMPUTE COUNTS PER Y CHARACTER POSITION
C
  220 MULTY = 1
  240 IF (MAXF.LE.NSTEPY*MULTY) GO TO 260
      MULTY = MULTY + 1
      GO TO 240
C
C     MODIFY HEADING FORMAT AND PRINT HEADING
C
  260 IS1 = (ISPACE+1)/10 + 1
      IS2 = ISPACE + 2 - 10*(IS1-1)
      IF1(2) = K(IS1)
      IF1(3) = K(IS2)
C
C     OBTAIN LOGICAL UNIT NUMBER FOR OUTPUT
C
      CALL X04ABF(0,NDEV)
      IF (ITYPE.EQ.0) WRITE (REC,FMT=IF1) (IH1(II),II=1,10)
      IF (ITYPE.NE.0) WRITE (REC,FMT=IF1) (IH2(II),II=1,10)
      CALL X04BAF(NDEV,REC)
C
C     MODIFY ROW FORMAT TO ALLOW FOR SPACING AND PLOT WIDTH
C
      IF3(2) = K(IS1)
      IF3(3) = K(IS2)
      IF4(2) = K(IS1)
      IF4(3) = K(IS2)
      NSPX2 = NSTEPX + 2
      NN = NSPX2
      IS1 = NN/100
      NN = NN - 100*IS1
      IS2 = NN/10
      IS3 = NN - 10*IS2
      IF3(12) = K(IS1+1)
      IF3(13) = K(IS2+1)
      IF3(14) = K(IS3+1)
C
C     CLEAR ROW FORMAT
C
      DO 280 I = 1, NSPX2
         LINE(I) = IBLANK
  280 CONTINUE
      LINE(1) = IDOT
      LINE(NSPX2) = IDOT
C
C     COMPUTE AND PRINT ROWS OF HISTOGRAM
C
      MINF = NSTEPY*MULTY
      DO 320 I = 1, NSTEPY
         DO 300 J = 1, NSTEPX
            IF (IFREQ(J).LT.MINF) GO TO 300
            IFREQ(J) = 0
            LINE(J+1) = ISTAR
  300    CONTINUE
         WRITE (REC,FMT=IF3) MINF, (LINE(II),II=1,NSPX2), MINF
         CALL X04BAF(NDEV,REC)
         MINF = MINF - MULTY
  320 CONTINUE
C
C     COMPLETE THE HISTOGRAM WITH LINE OF HYPHENS
C
      IF2(2) = IF3(2)
      IF2(3) = IF3(3)
      IF2(9) = IF3(12)
      IF2(10) = IF3(13)
      IF2(11) = IF3(14)
      DO 340 I = 1, NSPX2
         LINE(I) = IDASH
  340 CONTINUE
      WRITE (REC,FMT=IF2) (LINE(II),II=1,NSPX2)
      CALL X04BAF(NDEV,REC)
C
C     FOLLOWED BY XMIN AND XMAX
C
      IF (MAX(ABS(XMN),ABS(XMX)).GT.9999.99D0) GO TO 360
      NN = NSTEPX - 9
      IS1 = NN/10
      IS2 = NN - 10*IS1
      IF4(14) = K(IS1+1)
      IF4(15) = K(IS2+1)
      WRITE (REC,FMT=IF4) XMN, XMX
      CALL X04BAF(NDEV,REC)
  360 IFAIL = 0
      GO TO 400
  380 IFAIL = P01ABF(IFAIL,IFA,SRNAME,0,P01REC)
  400 RETURN
      END
