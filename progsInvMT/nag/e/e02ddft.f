      SUBROUTINE E02DDF(START,M,X,Y,F,W,S,NXEST,NYEST,NX,LAMDA,NY,MU,C,
     *                  FP,RANK,WRK,LWRK,IWRK,LIWRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      INTEGER           KX, KY, MAXIT
      DOUBLE PRECISION  TOL, ZERO
      PARAMETER         (KX=3,KY=3,MAXIT=20,TOL=0.1D-2,ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02DDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FP, S
      INTEGER           IFAIL, LIWRK, LWRK, M, NX, NXEST, NY, NYEST,
     *                  RANK
      CHARACTER*1       START
C     .. Array Arguments ..
      DOUBLE PRECISION  C((NXEST-4)*(NYEST-4)), F(M), LAMDA(NXEST),
     *                  MU(NYEST), W(M), WRK(LWRK), X(M), Y(M)
      INTEGER           IWRK(LIWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, XB, XE, YB, YE
      INTEGER           I, IB1, IB3, IER, IOPT, IW, JB1, KI, KM1, KM2,
     *                  KMAX, KN, KWEST, KX1, KY1, LA, LBX, LBY, LCO,
     *                  LF, LFF, LFP, LH, LKNOTX, LKNOTY, LQ, LSX, LSY,
     *                  LWEST, LWRK2, MM, NC, NEK, NEST, NMINX, NMINY,
     *                  NMX, NMY, NREC, NREG, NRINT, NXK, NYK
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02DDZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C     KX and KY are the degrees of the bivariate spline. They must be
C     in the range 0 - 5.
      IER = 1
      KX1 = KX + 1
      KY1 = KY + 1
      KMAX = MAX(KX,KY)
      KM1 = KMAX + 1
      KM2 = KM1 + 1
      NMINX = 2*KX1
      NMINY = 2*KY1
      NEST = MAX(NXEST,NYEST)
      NXK = NXEST - KX1
      NYK = NYEST - KY1
      NC = NXK*NYK
      NMX = NXEST - NMINX + 1
      NMY = NYEST - NMINY + 1
      NRINT = NMX + NMY
      NREG = NMX*NMY
      IB1 = KX*NYK + KY1
      JB1 = KY*NXK + KX1
      IB3 = KX1*NYK + 1
      IF (IB1.GT.JB1) THEN
         IB1 = JB1
         IB3 = KY1*NXK + 1
      END IF
      IW = MAX(NXK,NYK)
      LWEST = (7*NXK*NYK+25*IW)*(IW+1) + 2*(NXK+NYK+4*M) + 23*IW + 56
      KWEST = M + NREG + NREG
      IF (START.EQ.'C' .OR. START.EQ.'c') THEN
         IOPT = 0
      ELSE IF (START.EQ.'W' .OR. START.EQ.'w') THEN
         IOPT = 1
      ELSE
         IOPT = -1
      END IF
      NREC = 1
      IF (IOPT.EQ.-1) THEN
         WRITE (REC,FMT=99995) START
      ELSE IF (NXEST.LT.NMINX) THEN
         WRITE (REC,FMT=99993) 'NXEST', NMINX, 'NXEST', NXEST
      ELSE IF (NYEST.LT.NMINY) THEN
         WRITE (REC,FMT=99993) 'NYEST', NMINY, 'NYEST', NYEST
      ELSE IF (LWRK.LT.LWEST) THEN
         WRITE (REC,FMT=99994) 'LWRK', LWEST, 'LWRK', LWRK
      ELSE IF (LIWRK.LT.KWEST) THEN
         WRITE (REC,FMT=99994) 'LIWRK', KWEST, 'LIWRK', LIWRK
      ELSE IF (S.LE.ZERO) THEN
         WRITE (REC,FMT=99996) S
      ELSE
         IF (IOPT.EQ.1) THEN
C           Copy the current knots into workspace.
            DO 20 I = 1, NX
               WRK(I+1) = LAMDA(I)
   20       CONTINUE
            DO 40 I = 1, NY
               WRK(NEST+I+1) = MU(I)
   40       CONTINUE
C           We fetch the boundaries of the approximation domain.
            XB = LAMDA(KX1)
            NXK = NX - KX
            XE = LAMDA(NXK)
            YB = MU(KY1)
            NYK = NY - KY
            YE = MU(NYK)
         ELSE
C           We determine the boundaries of the approximation domain.
            XB = X(1)
            XE = XB
            YB = Y(1)
            YE = YB
            MM = 0
            DO 60 I = 1, M
               IF (W(I).NE.ZERO) MM = MM + 1
               IF (X(I).LT.XB) XB = X(I)
               IF (X(I).GT.XE) XE = X(I)
               IF (Y(I).LT.YB) YB = Y(I)
               IF (Y(I).GT.YE) YE = Y(I)
   60       CONTINUE
            IF (MM.LT.16) THEN
               NREC = 2
               WRITE (REC,FMT=99997) MM
               GO TO 120
            ELSE IF (XB.EQ.XE) THEN
               IER = 2
               WRITE (REC,FMT=99999)
               GO TO 120
            ELSE IF (YB.EQ.YE) THEN
               IER = 2
               WRITE (REC,FMT=99998)
               GO TO 120
            END IF
         END IF
         IER = 0
         NREC = 0
         EPS = X02AJF()
C        We partition the working space and determine the spline
C        approximation
         KN = 1
         KI = KN + NREG
         LKNOTX = 2
         LKNOTY = LKNOTX + NEST
         LQ = LKNOTY + NEST
         LA = LQ + (NC+IB3)*IB3
         LF = LA + (NC+IB1)*IB1
         LFF = LF + NC
         LFP = LFF + NC
         LCO = LFP + NRINT
         LH = LCO + NRINT
         LBX = LH + IB3
         NEK = NEST*KM2
         LBY = LBX + NEK
         LSX = LBY + NEK
         LSY = LSX + M*KM1
         LWRK2 = LWRK - LWEST + 1
         CALL E02DDZ(IOPT,M,X,Y,F,W,XB,XE,YB,YE,KX,KY,S,NXEST,NYEST,EPS,
     *               TOL,MAXIT,NEST,KM1,KM2,IB1,IB3,NC,NRINT,NREG,NX,
     *               WRK(LKNOTX),NY,WRK(LKNOTY),C,FP,WRK(1),RANK,
     *               WRK(LFP),WRK(LCO),WRK(LF),WRK(LFF),WRK(LA),WRK(LQ),
     *               WRK(LBX),WRK(LBY),WRK(LSX),WRK(LSY),WRK(LH),
     *               IWRK(KI),IWRK(KN),WRK(LWEST),LWRK2,IER)
C        Copy knots from workspace
         DO 80 I = 1, NX
            LAMDA(I) = WRK(I+1)
   80    CONTINUE
         DO 100 I = 1, NY
            MU(I) = WRK(NEST+I+1)
  100    CONTINUE
         IF (IER.EQ.3) THEN
            NREC = 3
            WRITE (REC,FMT=99992) NXEST, NYEST, S
         ELSE IF (IER.EQ.4) THEN
            NREC = 3
            WRITE (REC,FMT=99991) M, S
         ELSE IF (IER.EQ.5) THEN
            NREC = 3
            WRITE (REC,FMT=99990) S
         ELSE IF (IER.EQ.6) THEN
            NREC = 2
            WRITE (REC,FMT=99989) S
         ELSE IF (IER.GE.7) THEN
            IWRK(1) = IER + LWEST - 1
            IWRK(2) = LWEST + 4*NXK*NYK*IW + 2*NXK*NYK + 4*IW
            IER = 7
            NREC = 3
            WRITE (REC,FMT=99988) LWRK, IWRK(1),
     *        LWEST + 4*NXK*NYK*IW + 2*NXK*NYK + 4*IW
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, all the values in array X are equal.')
99998 FORMAT (1X,'** On entry, all the values in array Y are equal.')
99997 FORMAT (1X,'** On entry, MM .lt. 16, where MM is the number of p',
     *       'oints with',/4X,'non-zero weight: MM =',I16,'.')
99996 FORMAT (1X,'** On entry, S .le. 0.0: S =',1P,D13.5,' .')
99995 FORMAT (1X,'** On entry, START .ne. ''C'' or ''W'': START = ''',A,
     *       '''.')
99994 FORMAT (1X,'** On entry, ',A,' .lt.',I16,': ',A,' =',I16,'.')
99993 FORMAT (1X,'** On entry, ',A,' .lt.',I2,': ',A,' =',I16,'.')
99992 FORMAT (1X,'** The number of knots needed in one direction is gr',
     *       'eater than NXEST or NYEST:',/4X,'NXEST, NYEST =',I8,',',
     *       I8,'.',/4X,'Possibly S is too small: S =',1P,D13.5,' .')
99991 FORMAT (1X,'** No more knots added; the number of B-spline coeff',
     *       'icients already exceeds M.',/4X,'Either M or S is probab',
     *       'ly too small:',/4X,'M =',I10,', S =',1P,D13.5,' .')
99990 FORMAT (1X,'** No more knots added; the additional knot would co',
     *       'incide with an old one.',/4X,'Possibly an inaccurate dat',
     *       'a point has too large a weight, or S is',/4X,'too small:',
     *       ' S =',1P,D13.5,' .')
99989 FORMAT (1X,'** The iterative process has failed to converge.',/4X,
     *       'Possibly S is too small: S =',1P,D13.5,' .')
99988 FORMAT (1X,'** LWRK is too small to compute the minimal least-sq',
     *       'uares solution:',/4X,'LWRK =',I10,'. Minimum requested v',
     *       'alue for LWRK is',I10,';',/4X,'Safe requested value for ',
     *       'LWRK is',I10,'.')
      END
